
#include "pulseqlib_internal.h"

/* Helper: select the canonical TR window for a given canonical_tr_idx */
static void pulseqlib__select_canonical_tr_window_idx(
    const struct pulseqlib_sequence_descriptor* desc,
    int canonical_tr_idx,
    int* start_block,
    int* block_count,
    int* amplitude_mode,
    int* num_instances,
    float* tr_duration_us)
{
    const struct pulseqlib_tr_descriptor* trd = &desc->tr_descriptor;
    int has_nd_prep = (trd->num_prep_blocks > 0 && !trd->degenerate_prep);
    int has_nd_cool = (trd->num_cooldown_blocks > 0 && !trd->degenerate_cooldown);

    if (has_nd_prep || has_nd_cool) {
        /* Non-degenerate: canonical_tr_idx selects pass.
         * Pass-expanded waveform already bakes in num_averages via
         * _build_pass_expanded_block_order; num_instances counts how
         * many times this pass-expanded waveform repeats. */
        int pass_len = desc->pass_len;
        *start_block = canonical_tr_idx * pass_len;
        *block_count = pass_len;
        *amplitude_mode = PULSEQLIB_AMP_MAX_POS;
        *num_instances = (desc->num_passes > 1) ? desc->num_passes : 1;
        *tr_duration_us = 0.0f;
        {
            int n;
            for (n = 0; n < pass_len; ++n) {
                int idx;
                const struct pulseqlib_block_table_element* bte;
                const struct pulseqlib_block_definition* bdef;
                idx = *start_block + n;
                bte = &desc->block_table[idx];
                bdef = &desc->block_definitions[bte->id];
                *tr_duration_us += (bte->duration_us >= 0)
                    ? (float)bte->duration_us
                    : (float)bdef->duration_us;
            }
        }
        return;
    }

    /* Degenerate: canonical_tr_idx selects imaging TR.
     * No pass expansion occurs, so averages must be accounted for here. */
    {
        int num_avgs = (desc->num_averages > 1) ? desc->num_averages : 1;
        *start_block = trd->num_prep_blocks + trd->imaging_tr_start + canonical_tr_idx * trd->tr_size;
        *block_count = trd->tr_size;
        *amplitude_mode = PULSEQLIB_AMP_MAX_POS;
        *num_instances = trd->num_trs * num_avgs;
        *tr_duration_us = trd->tr_duration_us;
    }
}

static int pulseqlib__build_pass_expanded_block_order(
    const struct pulseqlib_sequence_descriptor* desc,
    int pass_base,
    int** out_block_order,
    int* out_block_count,
    float* out_duration_us)
{
    const struct pulseqlib_tr_descriptor* trd;
    int prep_blk, img_len, cool_blk, num_avgs, exp_count;
    int* block_order;
    int avg_i, pos_i, n;

    if (!desc || !out_block_order || !out_block_count) {
        return PULSEQLIB_ERR_INVALID_ARGUMENT;
    }

    trd = &desc->tr_descriptor;
    prep_blk = trd->num_prep_blocks;
    img_len = trd->num_trs * trd->tr_size;
    cool_blk = trd->num_cooldown_blocks;
    num_avgs = (desc->num_averages > 0) ? desc->num_averages : 1;
    exp_count = prep_blk + num_avgs * img_len + cool_blk;

    if (exp_count <= 0) {
        return PULSEQLIB_ERR_TR_NO_BLOCKS;
    }

    block_order = (int*)PULSEQLIB_ALLOC((size_t)exp_count * sizeof(int));
    if (!block_order) {
        return PULSEQLIB_ERR_ALLOC_FAILED;
    }

    n = 0;
    for (pos_i = 0; pos_i < prep_blk; ++pos_i)
        block_order[n++] = pass_base + pos_i;
    for (avg_i = 0; avg_i < num_avgs; ++avg_i)
        for (pos_i = 0; pos_i < img_len; ++pos_i)
            block_order[n++] = pass_base + prep_blk + pos_i;
    for (pos_i = 0; pos_i < cool_blk; ++pos_i)
        block_order[n++] = pass_base + prep_blk + img_len + pos_i;

    if (out_duration_us) {
        *out_duration_us = 0.0f;
        for (pos_i = 0; pos_i < exp_count; ++pos_i) {
            int blk_idx;
            const struct pulseqlib_block_table_element* bte;
            const struct pulseqlib_block_definition* bdef;

            blk_idx = block_order[pos_i];
            if (blk_idx < 0 || blk_idx >= desc->num_blocks) {
                PULSEQLIB_FREE(block_order);
                return PULSEQLIB_ERR_INVALID_ARGUMENT;
            }
            bte = &desc->block_table[blk_idx];
            bdef = &desc->block_definitions[bte->id];
            *out_duration_us += (bte->duration_us >= 0)
                ? (float)bte->duration_us
                : (float)bdef->duration_us;
        }
    }

    *out_block_order = block_order;
    *out_block_count = exp_count;
    return PULSEQLIB_SUCCESS;
}

/* pulseqlib_safety.c -- safety checks, mechanical resonance analysis, PNS, segment timing
 *
 * Public functions:
 *   pulseqlib_check_safety
 *   pulseqlib_calc_mech_resonances   / _free
 *   pulseqlib_calc_pns               / _free
 *
 * Internal:
 *   pulseqlib__calc_segment_timing
 *   check_max_grad / check_grad_continuity / check_max_slew
 */

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "pulseqlib_internal.h"
#include "pulseqlib_methods.h"

#include "external_kiss_fft.h"
#include "external_kiss_fftr.h"

/* ================================================================== */
/*  File-scope constants                                               */
/* ================================================================== */
#define PNS_KERNEL_DURATION_FACTOR 20.0f


/* ================================================================== */
/*  K-space trajectory from uniform gradient waveforms                */
/* ================================================================== */

/*
 * Computes k-space trajectory by cumulative trapezoidal integration of
 * gradient waveforms (already on uniform raster).
 *
 * Output arrays (kx, ky, kz, krss) must be caller-allocated with at
 * least waveforms->gx.num_samples elements.  dt_us is returned for
 * convenience.
 */
static int compute_kspace_trajectory(
    const pulseqlib__uniform_grad_waveforms* waveforms,
    float* kx, float* ky, float* kz, float* krss,
    float* dt_us,
    const int* refocus_samples, int num_refocus,
    const int* excite_samples,  int num_excite)
{
    int i, n, r, e;
    float dt_s;
    float cum_x, cum_y, cum_z;
    float v;

    if (!waveforms || !waveforms->gx || waveforms->num_samples < 2)
        return PULSEQLIB_ERR_INVALID_ARGUMENT;

    n = waveforms->num_samples;
    *dt_us = waveforms->raster_us;
    dt_s   = (*dt_us) * 1e-6f;

    /* cumulative trapezoidal integration */
    cum_x = 0.0f; cum_y = 0.0f; cum_z = 0.0f;
    kx[0] = 0.0f; ky[0] = 0.0f; kz[0] = 0.0f;

    r = 0;  /* index into refocus_samples */
    e = 0;  /* index into excite_samples  */

    for (i = 1; i < n; ++i) {
        cum_x += 0.5f * (waveforms->gx[i - 1] + waveforms->gx[i]) * dt_s;
        cum_y += 0.5f * (waveforms->gy[i - 1] + waveforms->gy[i]) * dt_s;
        cum_z += 0.5f * (waveforms->gz[i - 1] + waveforms->gz[i]) * dt_s;

        /* reset k=0 at excitation RF isocenter (90 deg pulse) */
        if (excite_samples && e < num_excite &&
            i == excite_samples[e]) {
            cum_x = 0.0f;
            cum_y = 0.0f;
            cum_z = 0.0f;
            e++;
        }

        /* negate k at refocusing RF isocenter (180 deg pulse) */
        if (refocus_samples && r < num_refocus &&
            i == refocus_samples[r]) {
            cum_x = -cum_x;
            cum_y = -cum_y;
            cum_z = -cum_z;
            r++;
        }

        kx[i] = cum_x;
        ky[i] = cum_y;
        kz[i] = cum_z;
    }

    /* RSS magnitude */
    for (i = 0; i < n; ++i) {
        v = kx[i] * kx[i] + ky[i] * ky[i] + kz[i] * kz[i];
        krss[i] = (v > 0.0f) ? (float)sqrt((double)v) : 0.0f;
    }

    return PULSEQLIB_SUCCESS;
}

/* ================================================================== */
/*  Find k-space zero crossings (k=0 passages)                       */
/* ================================================================== */

/*
 * A zero crossing is a local minimum of krss that is <= threshold.
 * Floor convention: for symmetric plateaus the leftmost sample is kept.
 *
 * Two-pass protocol:
 *   find_kspace_zero_crossings(krss, n, thr, NULL, &cnt, 1);
 *   indices = PULSEQLIB_ALLOC(cnt * sizeof(int));
 *   find_kspace_zero_crossings(krss, n, thr, indices, &cnt, 0);
 */
static int find_kspace_zero_crossings(
    const float* krss, int n, float threshold,
    int* zero_indices, int* out_count, int count_only)
{
    int i, cnt;
    float prev, curr, next;

    if (!krss || n < 2 || !out_count)
        return PULSEQLIB_ERR_INVALID_ARGUMENT;

    cnt = 0;

    for (i = 0; i < n; ++i) {
        curr = krss[i];
        if (curr > threshold) continue;

        prev = (i > 0)     ? krss[i - 1] : curr + 1.0f;
        next = (i < n - 1) ? krss[i + 1] : curr + 1.0f;

        /* curr <= both neighbors => local minimum */
        if (curr <= prev && curr <= next) {
            /* plateau dedup: skip if left neighbor equals curr and
             * was itself a local min (already recorded) */
            if (i > 0 && krss[i - 1] == curr) {
                float pprev = (i > 1) ? krss[i - 2] : curr + 1.0f;
                if (krss[i - 1] <= pprev)
                    continue;
            }
            if (!count_only && zero_indices)
                zero_indices[cnt] = i;
            cnt++;
        }
    }

    *out_count = cnt;
    return PULSEQLIB_SUCCESS;
}

/* ================================================================== */
/*  Compute segment timing anchors                                    */
/* ================================================================== */

int pulseqlib__calc_segment_timing(
    pulseqlib_sequence_descriptor* desc,
    pulseqlib_diagnostic* diag)
{
    pulseqlib_diagnostic local_diag;
    int seg_idx, blk, block_idx, result;
    const pulseqlib_tr_segment* seg;
    const pulseqlib_block_table_element* bte;
    const pulseqlib_block_definition* bdef;
    int rf_count, adc_count;
    float t_accum, block_dur_us;
    int rf_raw;
    const pulseqlib_rf_definition* rdef;
    const pulseqlib_adc_definition* adef;
    const pulseqlib_rf_table_element* rte;
    pulseqlib_segment_rf_anchor* rf_arr;
    pulseqlib_segment_adc_anchor* adc_arr;
    int rf_def_id, adc_def_id;
    float adc_dur_us;

    /* k-space trajectory variables */
    pulseqlib__uniform_grad_waveforms min_waveforms;
    float *kx, *ky, *kz, *krss;
    int *kzero_indices;
    int num_kzero, n_samples;
    float dt_us, k_threshold;

    /* refocusing RF detection variables */
    int *refocus_samples;
    int  num_refocus;

    /* excitation RF detection variables */
    int *excite_samples;
    int  num_excite;

    /* ADC-to-kzero mapping variables */
    int a, s, kz_sample;
    float seg_time_offset;
    float kzero_in_adc;
    int pos_in_tr, num_prep;
    int tr_size;
    int adc_raster_start, adc_raster_end, min_raster_idx;
    float min_krss_val;

    int has_kspace;

    if (!diag) { pulseqlib_diagnostic_init(&local_diag); diag = &local_diag; }

    memset(&min_waveforms, 0, sizeof(min_waveforms));
    kx = NULL; ky = NULL; kz = NULL; krss = NULL;
    kzero_indices = NULL;
    refocus_samples = NULL;
    num_refocus = 0;
    excite_samples = NULL;
    num_excite = 0;
    num_kzero = 0;
    n_samples = 0;
    dt_us = 0.0f;
    has_kspace = 0;

    if (!desc || desc->num_unique_segments <= 0)
        return PULSEQLIB_SUCCESS;

    num_prep = desc->tr_descriptor.num_prep_blocks;
    tr_size  = desc->tr_descriptor.tr_size;

    /* ---- Step A: build min-amplitude k-space trajectory ---- */
    if (tr_size > 0) {
        result = pulseqlib__get_gradient_waveforms_range(desc, &min_waveforms, diag,
            num_prep, tr_size, PULSEQLIB_AMP_ZERO_VAR, NULL, 0, NULL);

        if (!PULSEQLIB_FAILED(result) && min_waveforms.num_samples >= 2) {
            n_samples = min_waveforms.num_samples;

            /* ---- Step A.1: find refocusing RF isocenters in TR ---- */
            /*
             * Walk all blocks in the main TR range and identify
             * refocusing pulses.  If rf_use is tagged in the file,
             * use it directly; otherwise auto-detect from flip
             * angle (|flip| within 10% of 180 deg).
             */
            {
                float rf_t_accum = 0.0f;
                int   rf_cap = 0, rb;

                /* first pass: count */
                for (rb = 0; rb < tr_size; ++rb) {
                    int bi = num_prep + rb;
                    if (bi < 0 || bi >= desc->num_blocks) continue;
                    bte = &desc->block_table[bi];
                    rf_raw = bte->rf_id;
                    if (rf_raw >= 0 && rf_raw < desc->rf_table_size) {
                        rte = &desc->rf_table[rf_raw];
                        rf_def_id = rte->id;
                        if (rf_def_id >= 0 && rf_def_id < desc->num_unique_rfs) {
                            int use = rte->rf_use;
                            if (use == PULSEQLIB_RF_USE_UNKNOWN) {
                                /* auto-detect: actual_flip = base_flip *
                                 * |amplitude| / base_amplitude */
                                rdef = &desc->rf_definitions[rf_def_id];
                                if (rdef->stats.base_amplitude_hz > 0.0f) {
                                    float ratio = (float)fabs((double)rte->amplitude) /
                                                  rdef->stats.base_amplitude_hz;
                                    float actual_flip = rdef->stats.flip_angle_deg * ratio;
                                    if (actual_flip > 162.0f && actual_flip < 198.0f)
                                        use = PULSEQLIB_RF_USE_REFOCUSING;
                                }
                            }
                            if (use == PULSEQLIB_RF_USE_REFOCUSING)
                                rf_cap++;
                        }
                    }
                }

                /* second pass: collect isocenter sample indices */
                if (rf_cap > 0) {
                    refocus_samples = (int*)PULSEQLIB_ALLOC(
                        (size_t)rf_cap * sizeof(int));
                }
                if (refocus_samples) {
                    rf_t_accum = 0.0f;
                    num_refocus = 0;
                    for (rb = 0; rb < tr_size; ++rb) {
                        int bi = num_prep + rb;
                        if (bi < 0 || bi >= desc->num_blocks) continue;
                        bte  = &desc->block_table[bi];
                        bdef = &desc->block_definitions[bte->id];
                        block_dur_us = (bte->duration_us >= 0)
                            ? (float)bte->duration_us
                            : (float)bdef->duration_us;

                        rf_raw = bte->rf_id;
                        if (rf_raw >= 0 && rf_raw < desc->rf_table_size) {
                            rte = &desc->rf_table[rf_raw];
                            rf_def_id = rte->id;
                            if (rf_def_id >= 0 && rf_def_id < desc->num_unique_rfs) {
                                int use = rte->rf_use;
                                if (use == PULSEQLIB_RF_USE_UNKNOWN) {
                                    rdef = &desc->rf_definitions[rf_def_id];
                                    if (rdef->stats.base_amplitude_hz > 0.0f) {
                                        float ratio = (float)fabs((double)rte->amplitude) /
                                                      rdef->stats.base_amplitude_hz;
                                        float actual_flip = rdef->stats.flip_angle_deg * ratio;
                                        if (actual_flip > 162.0f && actual_flip < 198.0f)
                                            use = PULSEQLIB_RF_USE_REFOCUSING;
                                    }
                                }
                                if (use == PULSEQLIB_RF_USE_REFOCUSING) {
                                    float iso_us;
                                    int   iso_sample;
                                    rdef = &desc->rf_definitions[rf_def_id];
                                    iso_us = rf_t_accum + (float)rdef->delay +
                                             (float)rdef->stats.isodelay_us;
                                    iso_sample = (int)(iso_us / min_waveforms.raster_us + 0.5f);
                                    if (iso_sample < 0) iso_sample = 0;
                                    if (iso_sample >= n_samples) iso_sample = n_samples - 1;
                                    refocus_samples[num_refocus++] = iso_sample;
                                }
                            }
                        }
                        rf_t_accum += block_dur_us;
                    }
                }
            }

            /* ---- Step A.1b: find excitation RF isocenters in TR ---- */
            {
                float rf_t_accum = 0.0f;
                int   ex_cap = 0, rb;

                /* first pass: count excitation pulses */
                for (rb = 0; rb < tr_size; ++rb) {
                    int bi = num_prep + rb;
                    if (bi < 0 || bi >= desc->num_blocks) continue;
                    bte = &desc->block_table[bi];
                    rf_raw = bte->rf_id;
                    if (rf_raw >= 0 && rf_raw < desc->rf_table_size) {
                        rte = &desc->rf_table[rf_raw];
                        rf_def_id = rte->id;
                        if (rf_def_id >= 0 && rf_def_id < desc->num_unique_rfs) {
                            if (rte->rf_use == PULSEQLIB_RF_USE_EXCITATION)
                                ex_cap++;
                        }
                    }
                }

                /* second pass: collect isocenter sample indices */
                if (ex_cap > 0) {
                    excite_samples = (int*)PULSEQLIB_ALLOC(
                        (size_t)ex_cap * sizeof(int));
                }
                if (excite_samples) {
                    rf_t_accum = 0.0f;
                    num_excite = 0;
                    for (rb = 0; rb < tr_size; ++rb) {
                        int bi = num_prep + rb;
                        if (bi < 0 || bi >= desc->num_blocks) continue;
                        bte  = &desc->block_table[bi];
                        bdef = &desc->block_definitions[bte->id];
                        block_dur_us = (bte->duration_us >= 0)
                            ? (float)bte->duration_us
                            : (float)bdef->duration_us;

                        rf_raw = bte->rf_id;
                        if (rf_raw >= 0 && rf_raw < desc->rf_table_size) {
                            rte = &desc->rf_table[rf_raw];
                            rf_def_id = rte->id;
                            if (rf_def_id >= 0 && rf_def_id < desc->num_unique_rfs &&
                                rte->rf_use == PULSEQLIB_RF_USE_EXCITATION) {
                                float iso_us;
                                int   iso_sample;
                                rdef = &desc->rf_definitions[rf_def_id];
                                iso_us = rf_t_accum + (float)rdef->delay +
                                         (float)rdef->stats.isodelay_us;
                                iso_sample = (int)(iso_us / min_waveforms.raster_us + 0.5f);
                                if (iso_sample < 0) iso_sample = 0;
                                if (iso_sample >= n_samples) iso_sample = n_samples - 1;
                                excite_samples[num_excite++] = iso_sample;
                            }
                        }
                        rf_t_accum += block_dur_us;
                    }
                }
            }

            /* ---- Step A.2: compute k-space trajectory ---- */
            kx   = (float*)PULSEQLIB_ALLOC((size_t)n_samples * sizeof(float));
            ky   = (float*)PULSEQLIB_ALLOC((size_t)n_samples * sizeof(float));
            kz   = (float*)PULSEQLIB_ALLOC((size_t)n_samples * sizeof(float));
            krss = (float*)PULSEQLIB_ALLOC((size_t)n_samples * sizeof(float));
            if (kx && ky && kz && krss) {
                result = compute_kspace_trajectory(&min_waveforms,
                    kx, ky, kz, krss, &dt_us,
                    refocus_samples, num_refocus,
                    excite_samples, num_excite);
                if (!PULSEQLIB_FAILED(result)) {
                    /* threshold = 1% of max |k| */
                    k_threshold = 0.0f;
                    for (a = 0; a < n_samples; ++a)
                        if (krss[a] > k_threshold) k_threshold = krss[a];
                    k_threshold *= 0.01f;
                    if (k_threshold < 1e-10f) k_threshold = 1e-10f;

                    find_kspace_zero_crossings(krss, n_samples, k_threshold,
                                               NULL, &num_kzero, 1);
                    if (num_kzero > 0) {
                        kzero_indices = (int*)PULSEQLIB_ALLOC(
                            (size_t)num_kzero * sizeof(int));
                        if (kzero_indices) {
                            find_kspace_zero_crossings(krss, n_samples,
                                k_threshold, kzero_indices, &num_kzero, 0);
                            has_kspace = 1;
                        }
                    } else {
                        has_kspace = 1;  /* valid trajectory, just no crossings */
                    }
                }
            }
        }
    }

    /* ---- Step B: for each segment, collect RF and ADC anchors ---- */
    for (seg_idx = 0; seg_idx < desc->num_unique_segments; ++seg_idx) {
        seg = &desc->segment_definitions[seg_idx];

        /* count RF and ADC events (use block_definitions, not
           block_table, because start_block may point to a dummy TR
           instance where ADC/RF events are inactive)                */
        rf_count  = 0;
        adc_count = 0;
        for (blk = 0; blk < seg->num_blocks; ++blk) {
            int scan_idx = seg->max_energy_start_block + blk;
            if (scan_idx >= 0 && scan_idx < desc->scan_table_len)
                block_idx = desc->scan_table_block_idx[scan_idx];
            else
                block_idx = seg->start_block + blk;

            if (block_idx < 0 || block_idx >= desc->num_blocks) continue;
            bte  = &desc->block_table[block_idx];
            bdef = &desc->block_definitions[bte->id];
            if (bdef->rf_id >= 0)  rf_count++;
            if (bdef->adc_id >= 0) adc_count++;
        }

        /* allocate anchor arrays */
        rf_arr  = NULL;
        adc_arr = NULL;
        if (rf_count > 0) {
            rf_arr = (pulseqlib_segment_rf_anchor*)PULSEQLIB_ALLOC(
                (size_t)rf_count * sizeof(pulseqlib_segment_rf_anchor));
            if (!rf_arr) goto timing_fail;
        }
        if (adc_count > 0) {
            adc_arr = (pulseqlib_segment_adc_anchor*)PULSEQLIB_ALLOC(
                (size_t)adc_count * sizeof(pulseqlib_segment_adc_anchor));
            if (!adc_arr) {
                if (rf_arr) PULSEQLIB_FREE(rf_arr);
                goto timing_fail;
            }
        }

        /* compute segment start time within TR (for kzero mapping) */
        seg_time_offset = 0.0f;
        if (has_kspace && seg->start_block >= num_prep &&
            seg->start_block < num_prep + tr_size) {
            pos_in_tr = seg->start_block - num_prep;
            for (s = 0; s < pos_in_tr; ++s) {
                bte = &desc->block_table[num_prep + s];
                bdef = &desc->block_definitions[bte->id];
                seg_time_offset += (bte->duration_us >= 0)
                    ? (float)bte->duration_us
                    : (float)bdef->duration_us;
            }
        }

        /* Kzero estimation uses the full-TR canonical k-space trajectory
         * (PULSEQLIB_AMP_ZERO_VAR, computed above in min_waveforms/krss).
         * Variable-amplitude gradients (phase encodes, spiral readouts) are
         * zeroed; constant gradients (readout prephaser, slice select) are
         * kept — so the k-space trajectory correctly identifies the canonical
         * k-space zero crossing (N/2 for Cartesian, sample 0 for spiral).
         * seg_time_offset maps the segment's ADC window into the full-TR
         * raster.  No per-segment PULSEQLIB_AMP_ACTUAL extraction is used;
         * that would give instance-specific k and not the canonical result.  */

        /* fill anchors */
        rf_count  = 0;
        adc_count = 0;
        t_accum   = 0.0f;

        for (blk = 0; blk < seg->num_blocks; ++blk) {
            int scan_idx = seg->max_energy_start_block + blk;
            if (scan_idx >= 0 && scan_idx < desc->scan_table_len)
                block_idx = desc->scan_table_block_idx[scan_idx];
            else
                block_idx = seg->start_block + blk;

            if (block_idx < 0 || block_idx >= desc->num_blocks) continue;
            bte  = &desc->block_table[block_idx];
            bdef = &desc->block_definitions[bte->id];
            block_dur_us = (bte->duration_us >= 0)
                ? (float)bte->duration_us
                : (float)bdef->duration_us;

            /* RF anchor */
            rf_raw = bte->rf_id;
            if (rf_raw >= 0 && rf_raw < desc->rf_table_size) {
                rte = &desc->rf_table[rf_raw];
                rf_def_id = rte->id;
                if (rf_def_id >= 0 && rf_def_id < desc->num_unique_rfs) {
                    int use;
                    rdef = &desc->rf_definitions[rf_def_id];
                    rf_arr[rf_count].block_offset = blk;
                    rf_arr[rf_count].start_us     = t_accum + (float)rdef->delay;
                    rf_arr[rf_count].end_us        = t_accum + (float)rdef->delay +
                                                     rdef->stats.duration_us;
                    rf_arr[rf_count].isocenter_us   = t_accum + (float)rdef->delay +
                                                      (float)rdef->stats.isodelay_us;
                    rf_arr[rf_count].base_amplitude_hz = rte->amplitude;

                    /* rf_use: from file tag, or auto-detect from flip angle */
                    use = rte->rf_use;
                    if (use == PULSEQLIB_RF_USE_UNKNOWN &&
                        rdef->stats.base_amplitude_hz > 0.0f) {
                        float ratio = (float)fabs((double)rte->amplitude) /
                                      rdef->stats.base_amplitude_hz;
                        float actual_flip = rdef->stats.flip_angle_deg * ratio;
                        if (actual_flip > 162.0f && actual_flip < 198.0f)
                            use = PULSEQLIB_RF_USE_REFOCUSING;
                        else
                            use = PULSEQLIB_RF_USE_EXCITATION;
                    }
                    rf_arr[rf_count].rf_use = use;
                    rf_count++;
                }
            }

            /* ADC anchor */
            adc_def_id = bdef->adc_id;
            if (adc_def_id >= 0 && adc_def_id < desc->num_unique_adcs) {
                    adef = &desc->adc_definitions[adc_def_id];
                    adc_dur_us = (float)adef->num_samples *
                                 (float)adef->dwell_time * 1e-3f;

                    adc_arr[adc_count].block_offset = blk;
                    adc_arr[adc_count].start_us = t_accum + (float)adef->delay;
                    adc_arr[adc_count].end_us   = t_accum + (float)adef->delay +
                                                  adc_dur_us;

                    /* default: N/2 */
                    adc_arr[adc_count].kzero_index = adef->num_samples / 2;
                    adc_arr[adc_count].kzero_us    =
                        adc_arr[adc_count].start_us +
                        (float)(adef->num_samples / 2) *
                        (float)adef->dwell_time * 1e-3f;

                    /* full-TR canonical kRSS (ZERO_VAR): find min within ADC.
                     * Do not gate on seg->start_block being inside the first
                     * main TR: merged/full-pass segments (e.g. average-
                     * expanded multipass scans) still need kzero anchors for
                     * their ADC blocks. */
                    if (has_kspace && krss && n_samples > 0) {
                        adc_raster_start = (int)((seg_time_offset +
                            adc_arr[adc_count].start_us) / dt_us);
                        adc_raster_end   = (int)((seg_time_offset +
                            adc_arr[adc_count].end_us) / dt_us);
                        if (adc_raster_start < 0) adc_raster_start = 0;
                        if (adc_raster_end >= n_samples)
                            adc_raster_end = n_samples - 1;

                        if (adc_raster_start <= adc_raster_end) {
                            min_raster_idx = adc_raster_start;
                            min_krss_val   = krss[adc_raster_start];
                            for (a = adc_raster_start + 1;
                                 a <= adc_raster_end; ++a) {
                                if (krss[a] < min_krss_val) {
                                    min_krss_val   = krss[a];
                                    min_raster_idx = a;
                                }
                            }

                            kzero_in_adc = (float)min_raster_idx * dt_us -
                                (seg_time_offset + adc_arr[adc_count].start_us);
                            kz_sample = (int)(kzero_in_adc /
                                ((float)adef->dwell_time * 1e-3f));
                            if (kz_sample < 0) kz_sample = 0;
                            if (kz_sample >= adef->num_samples)
                                kz_sample = adef->num_samples - 1;
                            adc_arr[adc_count].kzero_index = kz_sample;
                            adc_arr[adc_count].kzero_us    =
                                (float)min_raster_idx * dt_us - seg_time_offset;
                        }
                    }

                    adc_count++;
            }

            t_accum += block_dur_us;
        }

        /* store timing */
        ((pulseqlib_tr_segment*)seg)->timing.num_rf_anchors  = rf_count;
        ((pulseqlib_tr_segment*)seg)->timing.rf_anchors      = rf_arr;
        ((pulseqlib_tr_segment*)seg)->timing.num_adc_anchors = adc_count;
        ((pulseqlib_tr_segment*)seg)->timing.adc_anchors     = adc_arr;
        ((pulseqlib_tr_segment*)seg)->timing.num_kzero_crossings = num_kzero;
        ((pulseqlib_tr_segment*)seg)->timing.kzero_crossing_indices = NULL;

        if (num_kzero > 0 && kzero_indices) {
            int* copy = (int*)PULSEQLIB_ALLOC((size_t)num_kzero * sizeof(int));
            if (copy) {
                int ci;
                for (ci = 0; ci < num_kzero; ++ci) copy[ci] = kzero_indices[ci];
                ((pulseqlib_tr_segment*)seg)->timing.kzero_crossing_indices = copy;
            }
        }
    }

    /* cleanup */
    if (kx)   PULSEQLIB_FREE(kx);
    if (ky)   PULSEQLIB_FREE(ky);
    if (kz)   PULSEQLIB_FREE(kz);
    if (krss) PULSEQLIB_FREE(krss);
    if (kzero_indices) PULSEQLIB_FREE(kzero_indices);
    if (refocus_samples) PULSEQLIB_FREE(refocus_samples);
    if (excite_samples)  PULSEQLIB_FREE(excite_samples);
    pulseqlib__uniform_grad_waveforms_free(&min_waveforms);

    return PULSEQLIB_SUCCESS;

timing_fail:
    if (kx)   PULSEQLIB_FREE(kx);
    if (ky)   PULSEQLIB_FREE(ky);
    if (kz)   PULSEQLIB_FREE(kz);
    if (krss) PULSEQLIB_FREE(krss);
    if (kzero_indices) PULSEQLIB_FREE(kzero_indices);
    if (refocus_samples) PULSEQLIB_FREE(refocus_samples);
    if (excite_samples)  PULSEQLIB_FREE(excite_samples);
    pulseqlib__uniform_grad_waveforms_free(&min_waveforms);
    return PULSEQLIB_ERR_ALLOC_FAILED;
}

/* ================================================================== */
/*  Acoustic spectra free                                             */
/* ================================================================== */

void pulseqlib_mech_resonances_spectra_free(pulseqlib_mech_resonances_spectra* s)
{
    if (!s) return;

    if (s->spectrogram_gx)    PULSEQLIB_FREE(s->spectrogram_gx);
    if (s->spectrogram_gy)    PULSEQLIB_FREE(s->spectrogram_gy);
    if (s->spectrogram_gz)    PULSEQLIB_FREE(s->spectrogram_gz);
    if (s->peaks_gx)          PULSEQLIB_FREE(s->peaks_gx);
    if (s->peaks_gy)          PULSEQLIB_FREE(s->peaks_gy);
    if (s->peaks_gz)          PULSEQLIB_FREE(s->peaks_gz);
    if (s->spectrum_full_gx)  PULSEQLIB_FREE(s->spectrum_full_gx);
    if (s->spectrum_full_gy)  PULSEQLIB_FREE(s->spectrum_full_gy);
    if (s->spectrum_full_gz)  PULSEQLIB_FREE(s->spectrum_full_gz);
    if (s->peaks_full_gx)    PULSEQLIB_FREE(s->peaks_full_gx);
    if (s->peaks_full_gy)    PULSEQLIB_FREE(s->peaks_full_gy);
    if (s->peaks_full_gz)    PULSEQLIB_FREE(s->peaks_full_gz);
    if (s->spectrum_seq_gx)   PULSEQLIB_FREE(s->spectrum_seq_gx);
    if (s->spectrum_seq_gy)   PULSEQLIB_FREE(s->spectrum_seq_gy);
    if (s->spectrum_seq_gz)   PULSEQLIB_FREE(s->spectrum_seq_gz);
    if (s->peaks_seq_gx)     PULSEQLIB_FREE(s->peaks_seq_gx);
    if (s->peaks_seq_gy)     PULSEQLIB_FREE(s->peaks_seq_gy);
    if (s->peaks_seq_gz)     PULSEQLIB_FREE(s->peaks_seq_gz);

    if (s->analytical_peak_freqs)     PULSEQLIB_FREE(s->analytical_peak_freqs);
    if (s->analytical_peak_amp_gx)    PULSEQLIB_FREE(s->analytical_peak_amp_gx);
    if (s->analytical_peak_amp_gy)    PULSEQLIB_FREE(s->analytical_peak_amp_gy);
    if (s->analytical_peak_amp_gz)    PULSEQLIB_FREE(s->analytical_peak_amp_gz);
    if (s->analytical_peak_phase_gx)  PULSEQLIB_FREE(s->analytical_peak_phase_gx);
    if (s->analytical_peak_phase_gy)  PULSEQLIB_FREE(s->analytical_peak_phase_gy);
    if (s->analytical_peak_phase_gz)  PULSEQLIB_FREE(s->analytical_peak_phase_gz);
    if (s->analytical_peak_widths_hz) PULSEQLIB_FREE(s->analytical_peak_widths_hz);
    if (s->candidate_freqs)           PULSEQLIB_FREE(s->candidate_freqs);
    if (s->candidate_amps_gx)         PULSEQLIB_FREE(s->candidate_amps_gx);
    if (s->candidate_amps_gy)         PULSEQLIB_FREE(s->candidate_amps_gy);
    if (s->candidate_amps_gz)         PULSEQLIB_FREE(s->candidate_amps_gz);
    if (s->candidate_grad_amps)       PULSEQLIB_FREE(s->candidate_grad_amps);
    if (s->candidate_grad_amps_gx)    PULSEQLIB_FREE(s->candidate_grad_amps_gx);
    if (s->candidate_grad_amps_gy)    PULSEQLIB_FREE(s->candidate_grad_amps_gy);
    if (s->candidate_grad_amps_gz)    PULSEQLIB_FREE(s->candidate_grad_amps_gz);
    if (s->candidate_violations)      PULSEQLIB_FREE(s->candidate_violations);
    if (s->component_freqs_hz)        PULSEQLIB_FREE(s->component_freqs_hz);
    if (s->component_amps)            PULSEQLIB_FREE(s->component_amps);
    if (s->component_phases_rad)      PULSEQLIB_FREE(s->component_phases_rad);
    if (s->component_widths_hz)       PULSEQLIB_FREE(s->component_widths_hz);
    if (s->component_axes)            PULSEQLIB_FREE(s->component_axes);
    if (s->component_def_ids)         PULSEQLIB_FREE(s->component_def_ids);
    if (s->component_contrib_ids)     PULSEQLIB_FREE(s->component_contrib_ids);
    if (s->component_run_ids)         PULSEQLIB_FREE(s->component_run_ids);
    if (s->surviving_freqs_hz)        PULSEQLIB_FREE(s->surviving_freqs_hz);

    memset(s, 0, sizeof(*s));
}

/* ================================================================== */
/*  Acoustic support structure                                        */
/* ================================================================== */

typedef struct {
    int   nwin;
    int   nfft;
    int   nfreq;
    int   output_freq_bins;
    int   num_windows;
    int   hop_size;
    float grad_raster_us;
    float freq_resolution;
    float max_frequency_hz;
    float* cos_window;
    float* work_buffer;
    kiss_fftr_cfg fft_cfg;
    kiss_fft_cpx* fft_out;
} acoustic_support;

typedef struct {
    int    num_samples;
    int    num_samples_original;
    float* samples;
    int    owns_memory;
} acoustic_waveform;

static int acoustic_support_init(
    acoustic_support* sup,
    int num_samples, int target_window_size,
    float target_spectral_resolution_hz,
    float grad_raster_us, float max_frequency_hz)
{
    int nwin, nfft, nfreq, output_freq_bins;
    int hop_size, num_windows, padded_len, min_nfft, max_idx, i;
    float freq_res;
    float* cos_win  = NULL;
    float* work     = NULL;
    kiss_fft_cpx* fft_out = NULL;
    kiss_fftr_cfg cfg = NULL;

    if (!sup || num_samples <= 0 || target_window_size <= 0 ||
        target_spectral_resolution_hz <= 0.0f)
        return PULSEQLIB_ERR_INVALID_ARGUMENT;

    memset(sup, 0, sizeof(*sup));

    nwin = (num_samples >= target_window_size) ? target_window_size : num_samples;

    min_nfft = (int)ceil((double)(1.0e6f / (grad_raster_us * target_spectral_resolution_hz)));
    nfft = (min_nfft < nwin) ? nwin : (int)pulseqlib__next_pow2((size_t)min_nfft);
    if (nfft < nwin) nfft = (int)pulseqlib__next_pow2((size_t)nwin);

    nfreq = nfft / 2 + 1;
    freq_res = (float)(1.0e6 / (grad_raster_us * (double)nfft));

    if (max_frequency_hz < 0.0f) {
        output_freq_bins = nfreq;
    } else {
        max_idx = (int)(max_frequency_hz / freq_res + 0.5f);
        if (max_idx >= nfreq) output_freq_bins = nfreq;
        else if (max_idx < 1) output_freq_bins = 1;
        else output_freq_bins = max_idx + 1;
    }

    hop_size = nwin / 2;
    if (hop_size < 1) hop_size = 1;

    if (num_samples <= nwin) {
        num_windows = 1;
    } else {
        padded_len  = ((num_samples + nwin - 1) / nwin) * nwin;
        num_windows = (padded_len - nwin) / hop_size + 1;
    }

    cos_win = (float*)PULSEQLIB_ALLOC(nwin * sizeof(float));
    work    = (float*)PULSEQLIB_ALLOC(nfft * sizeof(float));
    fft_out = (kiss_fft_cpx*)PULSEQLIB_ALLOC(nfreq * sizeof(kiss_fft_cpx));
    if (!cos_win || !work || !fft_out) goto fail;

    cfg = kiss_fftr_alloc(nfft, 0, NULL, NULL);
    if (!cfg) goto fail;

    for (i = 0; i < nwin; ++i)
        cos_win[i] = 0.5f * (1.0f - (float)cos(2.0 * M_PI * (double)(i + 1) / (double)nwin));

    sup->nwin             = nwin;
    sup->nfft             = nfft;
    sup->nfreq            = nfreq;
    sup->output_freq_bins = output_freq_bins;
    sup->num_windows      = num_windows;
    sup->hop_size         = hop_size;
    sup->grad_raster_us   = grad_raster_us;
    sup->freq_resolution  = freq_res;
    sup->max_frequency_hz = max_frequency_hz;
    sup->cos_window       = cos_win;
    sup->work_buffer      = work;
    sup->fft_cfg          = cfg;
    sup->fft_out          = fft_out;
    return PULSEQLIB_SUCCESS;

fail:
    if (cos_win) PULSEQLIB_FREE(cos_win);
    if (work)    PULSEQLIB_FREE(work);
    if (fft_out) PULSEQLIB_FREE(fft_out);
    if (cfg)     kiss_fftr_free(cfg);
    return PULSEQLIB_ERR_ALLOC_FAILED;
}

static void acoustic_support_free(acoustic_support* sup)
{
    if (!sup) return;
    if (sup->cos_window)   PULSEQLIB_FREE(sup->cos_window);
    if (sup->work_buffer)  PULSEQLIB_FREE(sup->work_buffer);
    if (sup->fft_out)      PULSEQLIB_FREE(sup->fft_out);
    if (sup->fft_cfg)      kiss_fftr_free(sup->fft_cfg);
    memset(sup, 0, sizeof(*sup));
}

static int acoustic_waveform_init(
    acoustic_waveform* aw,
    const acoustic_support* sup,
    const float* waveform, int num_samples, int padded_len)
{
    float* buf;
    int i;

    buf = NULL;

    if (!aw || !sup || !waveform || num_samples <= 0)
        return PULSEQLIB_ERR_NULL_POINTER;

    memset(aw, 0, sizeof(*aw));
    aw->num_samples_original = num_samples;

    buf = (float*)PULSEQLIB_ALLOC(padded_len * sizeof(float));
    if (!buf) return PULSEQLIB_ERR_ALLOC_FAILED;

    for (i = 0; i < num_samples; ++i) buf[i] = waveform[i];
    for (i = num_samples; i < padded_len; ++i) buf[i] = 0.0f;

    aw->num_samples = padded_len;
    aw->samples     = buf;
    aw->owns_memory = 1;
    return PULSEQLIB_SUCCESS;
}

static void acoustic_waveform_free(acoustic_waveform* aw)
{
    if (!aw) return;
    if (aw->owns_memory && aw->samples) PULSEQLIB_FREE(aw->samples);
    memset(aw, 0, sizeof(*aw));
}

/* ================================================================== */
/*  Single window spectrum                                            */
/* ================================================================== */

static int compute_window_spectrum(
    acoustic_support* sup, float* spectrum,
    const acoustic_waveform* aw, int window_index)
{
    int i, start_idx;
    int nwin, nfft, nfreq, out_bins;
    float* work;
    float* cos_win;
    const float* samples;
    float mean, fft_norm;
    kiss_fft_cpx* fft_out;

    if (!spectrum || !sup || !aw || !aw->samples)
        return PULSEQLIB_ERR_NULL_POINTER;
    if (window_index < 0 || window_index >= sup->num_windows)
        return PULSEQLIB_ERR_INVALID_ARGUMENT;

    nwin     = sup->nwin;
    nfft     = sup->nfft;
    nfreq    = sup->nfreq;
    out_bins = sup->output_freq_bins;
    work     = sup->work_buffer;
    cos_win  = sup->cos_window;
    fft_out  = sup->fft_out;
    samples  = aw->samples;

    start_idx = window_index * sup->hop_size;
    if (start_idx + nwin > aw->num_samples)
        return PULSEQLIB_ERR_INVALID_ARGUMENT;

    for (i = 0; i < nwin; ++i)  work[i] = samples[start_idx + i];
    for (i = nwin; i < nfft; ++i) work[i] = 0.0f;

    mean = 0.0f;
    for (i = 0; i < nwin; ++i) mean += work[i];
    mean /= (float)nwin;
    for (i = 0; i < nwin; ++i) work[i] -= mean;

    for (i = 0; i < nwin; ++i) work[i] *= cos_win[i];

    kiss_fftr(sup->fft_cfg, work, fft_out);

    fft_norm = 1.0f / (float)nfft;
    for (i = 0; i < nfreq; ++i) { fft_out[i].r *= fft_norm; fft_out[i].i *= fft_norm; }

    for (i = 0; i < out_bins; ++i)
        spectrum[i] = (float)sqrt((double)(fft_out[i].r * fft_out[i].r +
                                           fft_out[i].i * fft_out[i].i));
    return PULSEQLIB_SUCCESS;
}

/* ================================================================== */
/*  Resonance peak detection                                          */
/* ================================================================== */

static int compare_float(const void* a, const void* b)
{
    float fa = *(const float*)a;
    float fb = *(const float*)b;
    return (fa > fb) - (fa < fb);
}

/* FIX: output first */
static void detect_resonances(
    int* peaks,
    const float* mag,
    int n,
    float peak_log10_threshold,
    float peak_norm_scale,
    float peak_eps,
    float peak_prominence)
{
    int i;
    float max_val, median_log, norm;
    float* log_buf;
    float left_min, right_min, prom;

    if (n <= 0 || !mag || !peaks) return;
    for (i = 0; i < n; ++i) peaks[i] = 0;
    if (n < 3) return;

    max_val = 0.0f;
    for (i = 0; i < n; ++i) if (mag[i] > max_val) max_val = mag[i];
    if (max_val <= 0.0f) return;

    /* compute log-scaled values and median baseline */
    log_buf = (float*)PULSEQLIB_ALLOC((size_t)n * sizeof(float));
    if (!log_buf) return;

    for (i = 0; i < n; ++i) {
        norm = (mag[i] / max_val + peak_eps) * peak_norm_scale;
        log_buf[i] = (float)log10((double)norm);
    }
    {
        float* sorted = (float*)PULSEQLIB_ALLOC((size_t)n * sizeof(float));
        if (!sorted) { PULSEQLIB_FREE(log_buf); return; }
        memcpy(sorted, log_buf, (size_t)n * sizeof(float));
        qsort(sorted, (size_t)n, sizeof(float), compare_float);
        if (n % 2 == 1)
            median_log = sorted[n / 2];
        else
            median_log = 0.5f * (sorted[n / 2 - 1] + sorted[n / 2]);
        PULSEQLIB_FREE(sorted);
    }

    for (i = 1; i < n - 1; ++i) {
        if (mag[i] > mag[i-1] && mag[i] > mag[i+1]) {
            if (log_buf[i] - median_log > peak_log10_threshold) {
                if (peak_prominence > 0.0f) {
                    /* left trough: lowest log value between previous higher peak and i */
                    left_min = log_buf[i];
                    {
                        int j;
                        for (j = i - 1; j >= 0; --j) {
                            if (log_buf[j] < left_min) left_min = log_buf[j];
                            if (log_buf[j] > log_buf[i]) break;
                        }
                    }
                    /* right trough: lowest log value between i and next higher peak */
                    right_min = log_buf[i];
                    {
                        int j;
                        for (j = i + 1; j < n; ++j) {
                            if (log_buf[j] < right_min) right_min = log_buf[j];
                            if (log_buf[j] > log_buf[i]) break;
                        }
                    }
                    prom = log_buf[i] - (left_min > right_min ? left_min : right_min);
                    if (prom < peak_prominence) continue;
                }
                peaks[i] = 1;
            }
        }
    }

    PULSEQLIB_FREE(log_buf);
}

/* ================================================================== */
/*  Structural acoustic analysis — data types                         */
/* ================================================================== */

/** Max vertices in piecewise-linear waveform envelope. */
#define SA_MAX_PWL_VERTICES 16

/** Minimum sub-period repetitions for decomposition. */
#define SA_MIN_SUB_PERIOD_REPS 3

/** Autocorrelation peak threshold (normalised 0..1). */
#define SA_AUTOCORR_THRESHOLD 0.5f

/** Coherence ratio threshold for Tier 1 candidate detection. */
#define SA_COHERENCE_RATIO_THRESHOLD 2.0f

/** Absolute amplitude floor for Tier 1 (Hz/m). */
#define SA_AMPLITUDE_FLOOR 1.0f

/** FFT power gate: candidate requires FFT RSS magnitude above this fraction
 *  of the peak RSS magnitude.  Prevents flagging TR harmonics where the
 *  actual gradient waveform has negligible spectral energy. */
#define SA_FFT_GATE_FRACTION 0.10f

/** FFT peak promotion: a dense-FFT local maximum must exceed this fraction
 *  of the peak RSS magnitude to be promoted into the evaluation set when
 *  it falls between TR harmonics.  Higher than the gate to avoid adding
 *  many minor spectral features as candidates. */
#define SA_FFT_PROMOTION_FRACTION 0.25f

/** Absolute amplitude safety threshold for single-event Tier 2 (Hz/m).
 * Only applies when an axis has exactly one event (coherence ratio is
 * undefined for N=1).  Value is in the same units as the coherent
 * magnitude — roughly amp_Hz_per_m × h(f) × N_events. */
#define SA_AMPLITUDE_SAFETY 1.0e6f

/**
 * A gradient event within the canonical TR.
 * Each event is a single waveform (or sub-event from decomposed arbitrary)
 * played at a specific time with a specific amplitude.
 */
typedef struct {
    int   def_id;           /**< gradient definition id                */
    int   def_index;        /**< index into grad_definitions[]         */
    float start_time_us;    /**< event start time within the TR (us)   */
    float amplitude;        /**< worst-case positional amplitude (Hz/m), signed */
    int   pwl_num_vertices; /**< >0 -> use piecewise-linear W(f)       */
    float pwl_times_us[SA_MAX_PWL_VERTICES]; /**< vertex times (us from event start) */
    float pwl_values[SA_MAX_PWL_VERTICES];   /**< vertex amplitudes (normalised)     */
} sa_event;

/** Per-axis event list. */
typedef struct {
    int num_events;
    sa_event* events;  /**< allocated [num_events] */
} sa_axis_events;

/** Structural analysis: event lists for all three axes. */
typedef struct {
    sa_axis_events axes[3]; /**< 0=gx, 1=gy, 2=gz */
} sa_structural_events;

static void sa_free_structural_events(sa_structural_events* se)
{
    int ax;
    if (!se) return;
    for (ax = 0; ax < 3; ++ax) {
        if (se->axes[ax].events)
            PULSEQLIB_FREE(se->axes[ax].events);
        se->axes[ax].events = NULL;
        se->axes[ax].num_events = 0;
    }
}

/* ================================================================== */
/*  Structural acoustic analysis — occurrence extraction              */
/* ================================================================== */

/**
 * Extract grad-def occurrences for one axis within the canonical TR.
 * Collects (def_id, def_index, start_time_us, amplitude) for each block
 * that has a gradient event on the given axis.
 * @param axis  0=gx, 1=gy, 2=gz
 * Returns number of occurrences, or -1 on allocation failure.
 * Caller must free *out_events.
 */
typedef struct {
    int   def_id;
    int   def_index;
    float start_time_us;
    float amplitude;
} sa_raw_occurrence;

static int sa_extract_raw_occurrences(
    sa_raw_occurrence** out_occ,
    const struct pulseqlib_sequence_descriptor* desc,
    int start_block, int block_count, int axis)
{
    int i, n, cap;
    float time_us;
    sa_raw_occurrence* occ;

    cap = (block_count > 0) ? block_count : 16;
    occ = (sa_raw_occurrence*)PULSEQLIB_ALLOC((size_t)cap * sizeof(sa_raw_occurrence));
    if (!occ) return -1;
    n = 0;
    time_us = 0.0f;

    for (i = 0; i < block_count; ++i) {
        const struct pulseqlib_block_table_element* bte =
            &desc->block_table[start_block + i];
        const struct pulseqlib_block_definition* bdef =
            &desc->block_definitions[bte->id];
        int grad_table_idx;
        float blk_dur;

        switch (axis) {
            case 0: grad_table_idx = bte->gx_id; break;
            case 1: grad_table_idx = bte->gy_id; break;
            case 2: grad_table_idx = bte->gz_id; break;
            default: grad_table_idx = -1; break;
        }

        if (grad_table_idx >= 0) {
            const struct pulseqlib_grad_table_element* gte =
                &desc->grad_table[grad_table_idx];
            if (n >= cap) {
                sa_raw_occurrence* tmp;
                cap *= 2;
                tmp = (sa_raw_occurrence*)PULSEQLIB_ALLOC((size_t)cap * sizeof(sa_raw_occurrence));
                if (!tmp) { PULSEQLIB_FREE(occ); return -1; }
                memcpy(tmp, occ, (size_t)n * sizeof(sa_raw_occurrence));
                PULSEQLIB_FREE(occ);
                occ = tmp;
            }
            occ[n].def_id        = gte->id;
            occ[n].def_index     = gte->id;
            occ[n].start_time_us = time_us;
            occ[n].amplitude     = gte->amplitude;
            ++n;
        }

        blk_dur = (bte->duration_us >= 0) ? (float)bte->duration_us
                                           : (float)bdef->duration_us;
        time_us += blk_dur;
    }

    *out_occ = occ;
    return n;
}

/* ================================================================== */
/*  Structural acoustic analysis — sub-period detection               */
/* ================================================================== */

/**
 * Detect dominant sub-period in an arbitrary waveform via autocorrelation
 * of the rectified signal.  Returns detected period in samples, or 0 if
 * no periodicity found.
 *
 * @param wave      Decompressed waveform samples.
 * @param n         Number of samples.
 * @param out_period  Detected period in samples (0 if none).
 * @param out_reps    Number of complete repetitions.
 * @return 1 if sub-period found, 0 otherwise.
 */
static int sa_detect_sub_period(
    const float* wave, int n,
    int* out_period, int* out_reps)
{
    int tau, max_lag, best_tau, reps, k;
    float mean_val, sum_r0, sum_rt, r_norm, best_r;
    float* buf;

    *out_period = 0;
    *out_reps   = 0;

    if (n < 6) return 0;  /* Too short for periodicity */

    buf = (float*)PULSEQLIB_ALLOC((size_t)n * sizeof(float));
    if (!buf) return 0;

    /* Remove mean (no rectification — preserve sign for bipolar detection) */
    mean_val = 0.0f;
    for (k = 0; k < n; ++k)
        mean_val += wave[k];
    mean_val /= (float)n;
    for (k = 0; k < n; ++k)
        buf[k] = wave[k] - mean_val;

    /* Compute R[0] for normalisation */
    sum_r0 = 0.0f;
    for (k = 0; k < n; ++k)
        sum_r0 += buf[k] * buf[k];
    if (sum_r0 < 1.0e-20f) {
        PULSEQLIB_FREE(buf);
        return 0;
    }

    /* Search for first significant peak in normalised autocorrelation.
     * Skip lag 0 (trivially 1.0).  Start from lag 2 to avoid edge effects. */
    max_lag = n / 2;
    best_tau = 0;
    best_r = SA_AUTOCORR_THRESHOLD;  /* Minimum acceptable peak height */

    /* First, find the first valley (autocorrelation must dip below threshold
     * before we accept a peak, to avoid the initial decay region). */
    {
        int found_valley = 0;
        float prev_r = 1.0f;
        for (tau = 2; tau < max_lag; ++tau) {
            sum_rt = 0.0f;
            for (k = 0; k < n - tau; ++k)
                sum_rt += buf[k] * buf[k + tau];
            r_norm = sum_rt / sum_r0;
            if (r_norm < SA_AUTOCORR_THRESHOLD * 0.5f) {
                found_valley = 1;
            }
            if (found_valley && r_norm > best_r &&
                r_norm > prev_r) {
                /* Check it's a local max: look one step ahead */
                float r_next = 0.0f;
                int tau_next = tau + 1;
                if (tau_next < max_lag) {
                    int kk;
                    for (kk = 0; kk < n - tau_next; ++kk)
                        r_next += buf[kk] * buf[kk + tau_next];
                    r_next /= sum_r0;
                }
                if (r_norm >= r_next) {
                    best_tau = tau;
                    best_r = r_norm;
                    break;  /* First significant peak after valley */
                }
            }
            prev_r = r_norm;
        }
    }

    PULSEQLIB_FREE(buf);

    if (best_tau < 2) return 0;

    /* Validate: check that R[2*tau] is also a peak (periodic confirmation) */
    {
        float r2;
        int tau2 = 2 * best_tau;
        if (tau2 < max_lag) {
            float* buf2 = (float*)PULSEQLIB_ALLOC((size_t)n * sizeof(float));
            if (buf2) {
                for (k = 0; k < n; ++k)
                    buf2[k] = wave[k] - mean_val;
                sum_rt = 0.0f;
                for (k = 0; k < n - tau2; ++k)
                    sum_rt += buf2[k] * buf2[k + tau2];
                r2 = sum_rt / sum_r0;
                PULSEQLIB_FREE(buf2);
                if (r2 < SA_AUTOCORR_THRESHOLD * 0.5f)
                    return 0;  /* Second harmonic too weak — not periodic */
            }
        }
    }

    reps = n / best_tau;
    if (reps < SA_MIN_SUB_PERIOD_REPS) return 0;

    *out_period = best_tau;
    *out_reps   = reps;
    return 1;
}

/* ================================================================== */
/*  Structural acoustic analysis — build events per axis              */
/* ================================================================== */

static int sa_compare_int(const void* a, const void* b)
{
    return (*(const int*)a) - (*(const int*)b);
}

/**
 * Build event list for one axis from raw occurrences.
 *
 * For each unique def_id:
 *   - Trapezoid: one event per occurrence, PWL from rise/flat/fall.
 *   - Arbitrary with few samples: one event per occurrence, PWL approximation.
 *   - Arbitrary with many samples + detected sub-period: M sub-events per
 *     occurrence, each with PWL of one period.
 *   - Arbitrary with many samples + no sub-period: one event per occurrence,
 *     PWL approximation of full waveform (keeping first/last + zero crossings).
 */
static int sa_build_axis_events(
    sa_axis_events* ae,
    const sa_raw_occurrence* occ, int num_occ,
    const struct pulseqlib_sequence_descriptor* desc)
{
    int i, j, n_events, cap, did, idx, s_idx, nv;
    int* unique_ids;
    int num_unique;
    float raster;

    ae->num_events = 0;
    ae->events = NULL;
    if (num_occ <= 0) return PULSEQLIB_SUCCESS;

    raster = desc->grad_raster_us;

    /* Collect unique def_ids */
    unique_ids = (int*)PULSEQLIB_ALLOC((size_t)num_occ * sizeof(int));
    if (!unique_ids) return PULSEQLIB_ERR_ALLOC_FAILED;
    for (i = 0; i < num_occ; ++i) unique_ids[i] = occ[i].def_id;
    qsort(unique_ids, (size_t)num_occ, sizeof(int), sa_compare_int);
    num_unique = 1;
    for (i = 1; i < num_occ; ++i)
        if (unique_ids[i] != unique_ids[i - 1])
            unique_ids[num_unique++] = unique_ids[i];

    /* Initial allocation — may grow for decomposed arbitraries */
    cap = num_occ * 2;
    ae->events = (sa_event*)PULSEQLIB_ALLOC((size_t)cap * sizeof(sa_event));
    if (!ae->events) { PULSEQLIB_FREE(unique_ids); return PULSEQLIB_ERR_ALLOC_FAILED; }
    n_events = 0;

    for (idx = 0; idx < num_unique; ++idx) {
        const struct pulseqlib_grad_definition* gdef;
        /* Shared PWL template and sub-event decomposition info */
        int   pwl_nv;
        float pwl_t[SA_MAX_PWL_VERTICES];
        float pwl_v[SA_MAX_PWL_VERTICES];
        int   decompose;
        int   sub_reps;
        int   sub_pwl_nv;
        float sub_pwl_t[SA_MAX_PWL_VERTICES];
        float sub_pwl_v[SA_MAX_PWL_VERTICES];
        float sub_period_us;

        did = unique_ids[idx];
        gdef = NULL;
        pwl_nv = 0;
        decompose = 0;
        sub_reps = 0;
        sub_pwl_nv = 0;
        sub_period_us = 0.0f;

        /* Find first occurrence of this def to get gdef */
        for (j = 0; j < num_occ; ++j) {
            if (occ[j].def_id == did) {
                gdef = &desc->grad_definitions[occ[j].def_index];
                break;
            }
        }
        if (!gdef) continue;

        /* --- Compute shared PWL parameters for this def --- */

        if (gdef->type == 0) {
            /* Trapezoid: 4-vertex PWL */
            float tr = (float)gdef->rise_time_or_unused;
            float tf = (float)gdef->flat_time_or_unused;
            float td = (float)gdef->fall_time_or_num_uncompressed_samples;
            pwl_nv = 4;
            pwl_t[0] = 0.0f;       pwl_v[0] = 0.0f;
            pwl_t[1] = tr;          pwl_v[1] = 1.0f;
            pwl_t[2] = tr + tf;     pwl_v[2] = 1.0f;
            pwl_t[3] = tr + tf + td; pwl_v[3] = 0.0f;
        }
        else if (gdef->type == 1) {
            /* Arbitrary waveform */
            int shape_id = gdef->shot_shape_ids[0];
            int time_shape_id = gdef->unused_or_time_shape_id;

            if (shape_id > 0 && shape_id <= desc->num_shapes) {
                pulseqlib_shape_arbitrary decomp_wave;
                decomp_wave.samples = NULL;
                decomp_wave.num_samples = 0;
                decomp_wave.num_uncompressed_samples = 0;

                if (pulseqlib__decompress_shape(&decomp_wave,
                        &desc->shapes[shape_id - 1], 1.0f)
                    && decomp_wave.num_uncompressed_samples > 0)
                {
                    int num_samp = decomp_wave.num_uncompressed_samples;
                    float* wave_samp = decomp_wave.samples;

                    /* Decompress time shape if present */
                    pulseqlib_shape_arbitrary decomp_time;
                    float* time_us = NULL;
                    int has_time = 0;
                    decomp_time.samples = NULL;
                    decomp_time.num_samples = 0;
                    decomp_time.num_uncompressed_samples = 0;

                    if (time_shape_id > 0 && time_shape_id <= desc->num_shapes) {
                        if (pulseqlib__decompress_shape(&decomp_time,
                                &desc->shapes[time_shape_id - 1], raster)
                            && decomp_time.num_uncompressed_samples > 0)
                        {
                            has_time = 1;
                            time_us = decomp_time.samples;
                        }
                    }

                    if (num_samp < PULSEQLIB_MIN_ARBITRARY_SAMPLES) {
                        /* Few-sample arb: PWL from vertices */
                        nv = num_samp;
                        if (nv > SA_MAX_PWL_VERTICES) nv = SA_MAX_PWL_VERTICES;
                        pwl_nv = nv;
                        for (s_idx = 0; s_idx < nv; ++s_idx) {
                            if (has_time && time_us)
                                pwl_t[s_idx] = time_us[s_idx];
                            else
                                pwl_t[s_idx] = 0.5f * raster + (float)s_idx * raster;
                            pwl_v[s_idx] = wave_samp[s_idx];
                        }
                    }
                    else {
                        /* Many-sample arb: try sub-period detection */
                        int sp, sr;
                        if (sa_detect_sub_period(wave_samp, num_samp, &sp, &sr)) {
                            /* Sub-period found: build PWL for one period */
                            decompose = 1;
                            sub_reps = sr;
                            sub_period_us = (has_time && time_us && sp < num_samp)
                                ? (time_us[sp] - time_us[0])
                                : (float)sp * raster;

                            /* Build PWL for one sub-period */
                            nv = sp;
                            if (nv > SA_MAX_PWL_VERTICES) nv = SA_MAX_PWL_VERTICES;
                            sub_pwl_nv = nv;
                            if (nv == sp) {
                                /* Use all samples */
                                for (s_idx = 0; s_idx < nv; ++s_idx) {
                                    if (has_time && time_us)
                                        sub_pwl_t[s_idx] = time_us[s_idx] - time_us[0];
                                    else
                                        sub_pwl_t[s_idx] = (float)s_idx * raster;
                                    sub_pwl_v[s_idx] = wave_samp[s_idx];
                                }
                            } else {
                                /* Subsample uniformly */
                                for (s_idx = 0; s_idx < nv; ++s_idx) {
                                    int si = (s_idx * (sp - 1)) / (nv - 1);
                                    if (has_time && time_us)
                                        sub_pwl_t[s_idx] = time_us[si] - time_us[0];
                                    else
                                        sub_pwl_t[s_idx] = (float)si * raster;
                                    sub_pwl_v[s_idx] = wave_samp[si];
                                }
                            }
                        }
                        else {
                            /* No sub-period: build PWL approximation of full waveform.
                             * Use SA_MAX_PWL_VERTICES uniformly sampled points. */
                            nv = SA_MAX_PWL_VERTICES;
                            if (nv > num_samp) nv = num_samp;
                            pwl_nv = nv;
                            for (s_idx = 0; s_idx < nv; ++s_idx) {
                                int si = (nv > 1) ? (s_idx * (num_samp - 1)) / (nv - 1) : 0;
                                if (has_time && time_us)
                                    pwl_t[s_idx] = time_us[si];
                                else
                                    pwl_t[s_idx] = 0.5f * raster + (float)si * raster;
                                pwl_v[s_idx] = wave_samp[si];
                            }
                        }
                    }

                    if (decomp_time.samples) PULSEQLIB_FREE(decomp_time.samples);
                }
                if (decomp_wave.samples) PULSEQLIB_FREE(decomp_wave.samples);
            }
        }

        /* --- Emit events for each occurrence of this def --- */
        for (j = 0; j < num_occ; ++j) {
            if (occ[j].def_id != did) continue;

            if (decompose && sub_reps > 0) {
                /* Decompose into sub-events */
                int m;
                for (m = 0; m < sub_reps; ++m) {
                    if (n_events >= cap) {
                        sa_event* tmp;
                        cap *= 2;
                        tmp = (sa_event*)PULSEQLIB_ALLOC((size_t)cap * sizeof(sa_event));
                        if (!tmp) { PULSEQLIB_FREE(unique_ids); PULSEQLIB_FREE(ae->events); ae->events = NULL; return PULSEQLIB_ERR_ALLOC_FAILED; }
                        memcpy(tmp, ae->events, (size_t)n_events * sizeof(sa_event));
                        PULSEQLIB_FREE(ae->events);
                        ae->events = tmp;
                    }
                    ae->events[n_events].def_id = did;
                    ae->events[n_events].def_index = occ[j].def_index;
                    ae->events[n_events].start_time_us = occ[j].start_time_us
                        + (float)m * sub_period_us
                        + (float)gdef->delay;
                    ae->events[n_events].amplitude = occ[j].amplitude;
                    ae->events[n_events].pwl_num_vertices = sub_pwl_nv;
                    memcpy(ae->events[n_events].pwl_times_us, sub_pwl_t,
                           (size_t)sub_pwl_nv * sizeof(float));
                    memcpy(ae->events[n_events].pwl_values, sub_pwl_v,
                           (size_t)sub_pwl_nv * sizeof(float));
                    n_events++;
                }
            }
            else {
                /* Single event */
                if (n_events >= cap) {
                    sa_event* tmp;
                    cap *= 2;
                    tmp = (sa_event*)PULSEQLIB_ALLOC((size_t)cap * sizeof(sa_event));
                    if (!tmp) { PULSEQLIB_FREE(unique_ids); PULSEQLIB_FREE(ae->events); ae->events = NULL; return PULSEQLIB_ERR_ALLOC_FAILED; }
                    memcpy(tmp, ae->events, (size_t)n_events * sizeof(sa_event));
                    PULSEQLIB_FREE(ae->events);
                    ae->events = tmp;
                }
                ae->events[n_events].def_id = did;
                ae->events[n_events].def_index = occ[j].def_index;
                ae->events[n_events].start_time_us = occ[j].start_time_us
                    + (float)gdef->delay;
                ae->events[n_events].amplitude = occ[j].amplitude;
                ae->events[n_events].pwl_num_vertices = pwl_nv;
                if (pwl_nv > 0) {
                    memcpy(ae->events[n_events].pwl_times_us, pwl_t,
                           (size_t)pwl_nv * sizeof(float));
                    memcpy(ae->events[n_events].pwl_values, pwl_v,
                           (size_t)pwl_nv * sizeof(float));
                }
                n_events++;
            }
        }
    }

    ae->num_events = n_events;
    PULSEQLIB_FREE(unique_ids);
    return PULSEQLIB_SUCCESS;
}

/**
 * Build event lists for all three axes.
 */
static int sa_build_structural_events(
    sa_structural_events* se,
    const struct pulseqlib_sequence_descriptor* desc,
    int start_block, int block_count)
{
    int ax, result;
    sa_raw_occurrence* occ;
    int num_occ;

    memset(se, 0, sizeof(*se));

    for (ax = 0; ax < 3; ++ax) {
        occ = NULL;
        num_occ = sa_extract_raw_occurrences(&occ, desc, start_block, block_count, ax);
        if (num_occ < 0) {
            sa_free_structural_events(se);
            return PULSEQLIB_ERR_ALLOC_FAILED;
        }
        result = sa_build_axis_events(&se->axes[ax], occ, num_occ, desc);
        PULSEQLIB_FREE(occ);
        if (PULSEQLIB_FAILED(result)) {
            sa_free_structural_events(se);
            return result;
        }
    }
    return PULSEQLIB_SUCCESS;
}

/* ================================================================== */
/*  Structural acoustic analysis — analytical spectrum evaluation     */
/* ================================================================== */

/**
 * Evaluate |G(f)| / |G(0)| for a piecewise-linear waveform defined by
 * n_vtx vertices (t_us[i], v[i]).  Times are in microseconds.
 *
 * For each linear segment [t_k, t_{k+1}] with values [v_k, v_{k+1}]:
 *   G_k(f) = e^{-j*omega*t_k} * [a*I_0 + b*I_1]
 * where a = v_k, b = (v_{k+1}-v_k)/T_k, omega = 2*pi*f,
 *   I_0 = integral_0^T e^{-j*omega*tau} dtau
 *   I_1 = integral_0^T tau * e^{-j*omega*tau} dtau
 *
 * Returns |G(f)| / |G(0)|.
 */
static float sa_eval_pwl_response(
    float f_hz, const float* t_us, const float* v, int n_vtx)
{
    double omega, g0, g_re, g_im, mag;
    int k;

    if (n_vtx < 2) return 1.0f;

    g0 = 0.0;
    for (k = 0; k < n_vtx - 1; ++k) {
        double dt = (double)(t_us[k + 1] - t_us[k]);
        g0 += 0.5 * ((double)v[k] + (double)v[k + 1]) * dt;
    }
    if (g0 < 1.0e-30 && g0 > -1.0e-30) return 1.0f;
    g0 = fabs(g0);

    omega = 2.0 * M_PI * (double)f_hz * 1.0e-6;
    if (omega < 1.0e-12 && omega > -1.0e-12) return 1.0f;

    g_re = 0.0;
    g_im = 0.0;

    for (k = 0; k < n_vtx - 1; ++k) {
        double tk  = (double)t_us[k];
        double Tk  = (double)(t_us[k + 1] - t_us[k]);
        double ak  = (double)v[k];
        double bk  = ((double)v[k + 1] - ak) / Tk;
        double wTk = omega * Tk;
        double cos_wTk = cos(wTk);
        double sin_wTk = sin(wTk);
        double cos_ph  = cos(omega * tk);
        double sin_ph  = sin(omega * tk);
        double c0_re, c0_im, c1_re, c1_im;
        double x_re, x_im, seg_re, seg_im;

        if (Tk < 1.0e-12) continue;

        c0_re = sin_wTk / omega;
        c0_im = (cos_wTk - 1.0) / omega;

        {
            double I0mT_re = c0_re - Tk * cos_wTk;
            double I0mT_im = c0_im + Tk * sin_wTk;
            c1_re =  I0mT_im / omega;
            c1_im = -I0mT_re / omega;
        }

        x_re = ak * c0_re + bk * c1_re;
        x_im = ak * c0_im + bk * c1_im;

        seg_re = x_re * cos_ph + x_im * sin_ph;
        seg_im = x_im * cos_ph - x_re * sin_ph;

        g_re += seg_re;
        g_im += seg_im;
    }

    mag = sqrt(g_re * g_re + g_im * g_im);
    return (float)(mag / g0);
}

/**
 * Evaluate the waveform response W_k(f) for a single event.
 * Returns real-valued |W(f)| / |W(0)| (normalised).
 */
static float sa_eval_event_response(const sa_event* ev, float f_hz)
{
    if (ev->pwl_num_vertices >= 2)
        return sa_eval_pwl_response(f_hz, ev->pwl_times_us, ev->pwl_values,
                                    ev->pwl_num_vertices);
    return 1.0f;  /* Broadband fallback */
}

/**
 * Evaluate the complex spectral contribution of one event at frequency f:
 *   a_k(f) = A_k * W_k(f) * e^{-j 2 pi f t_k}
 */
static void sa_eval_event_spectrum(
    float f_hz, const sa_event* ev,
    float* out_re, float* out_im)
{
    float amp, h;
    double phase, cos_ph, sin_ph;

    amp = ev->amplitude;
    h = sa_eval_event_response(ev, f_hz);

    phase = -2.0 * M_PI * (double)f_hz * (double)ev->start_time_us * 1.0e-6;
    cos_ph = cos(phase);
    sin_ph = sin(phase);

    *out_re = amp * h * (float)cos_ph;
    *out_im = amp * h * (float)sin_ph;
}

/**
 * Evaluate the inner TR spectrum magnitude for one axis at frequency f:
 *   |H(f)| = |sum_k a_k(f)|
 * Also returns the incoherent sum sqrt(sum |a_k|^2) for coherence ratio.
 */
static void sa_eval_inner_spectrum(
    float f_hz,
    const sa_axis_events* ae,
    float* out_coherent_mag,
    float* out_incoherent_mag)
{
    float sum_re, sum_im;
    float inc_sum;
    int k;

    sum_re = 0.0f;
    sum_im = 0.0f;
    inc_sum = 0.0f;

    for (k = 0; k < ae->num_events; ++k) {
        float d_re, d_im, d_mag2;
        sa_eval_event_spectrum(f_hz, &ae->events[k], &d_re, &d_im);
        sum_re += d_re;
        sum_im += d_im;
        d_mag2 = d_re * d_re + d_im * d_im;
        inc_sum += d_mag2;
    }

    *out_coherent_mag = (float)sqrt((double)(sum_re * sum_re + sum_im * sum_im));
    *out_incoherent_mag = (float)sqrt((double)inc_sum);
}

/**
 * Compute spectral-weighted effective gradient amplitude G_eff at frequency f:
 *   G_eff(f) = sum_k |A_k| * |a_k(f)| / sum_k |a_k(f)|
 */
static float sa_eval_geff(
    float f_hz,
    const sa_axis_events* ae)
{
    float w_sum, w_amp_sum;
    int k;

    w_sum = 0.0f;
    w_amp_sum = 0.0f;

    for (k = 0; k < ae->num_events; ++k) {
        float d_re, d_im, d_mag;
        float abs_amp;
        sa_eval_event_spectrum(f_hz, &ae->events[k], &d_re, &d_im);
        d_mag = (float)sqrt((double)(d_re * d_re + d_im * d_im));
        abs_amp = ae->events[k].amplitude;
        if (abs_amp < 0.0f) abs_amp = -abs_amp;
        w_sum += d_mag;
        w_amp_sum += d_mag * abs_amp;
    }

    if (w_sum <= 1.0e-20f) return 0.0f;
    return w_amp_sum / w_sum;
}

/* ================================================================== */
/*  Structural acoustic analysis — top-level violation check          */
/* ================================================================== */

/**
 * Run full structural acoustic analysis with event-level analytical framework.
 *
 * For each axis:
 *   1. Extract block-level occurrences, decompose arbitrary waveforms via
 *      autocorrelation if sub-periodicity detected.
 *   2. Build candidate frequency grid (TR harmonics; K is clamped to
 *      min 2 so the Fourier comb works even for single-TR sequences).
 *   3. Evaluate inner TR spectrum H(f) at each candidate via coherent sum
 *      of per-event contributions.
 *   4. Two-tier candidate selection with FFT power gate:
 *      - Tier 1: coherence ratio > threshold AND |H| > floor AND FFT gate
 *      - Tier 2: |H| > safety threshold AND FFT gate
 *   5. For each candidate: compute G_eff(f) and check forbidden bands.
 */
static int sa_check_structural_violations(
    pulseqlib_mech_resonances_spectra* spectra,
    const struct pulseqlib_sequence_descriptor* desc,
    int start_block, int block_count,
    int num_instances, float tr_duration_us,
    int num_forbidden_bands,
    const pulseqlib_forbidden_band* forbidden_bands,
    float peak_log10_threshold, float peak_norm_scale,
    float peak_eps, float peak_prominence)
{
    sa_structural_events se;
    int result, ax, i, b, ci;
    int num_freqs, num_candidates;
    float f_hz, freq_max;
    float* eval_freqs;
    float* coherent_mag[3];
    float* incoherent_mag[3];
    int*   is_candidate;
    float max_fft_rss_sq, fft_gate_sq, fft_promo_sq;
    int has_fft_data;
    int n_tr_harmonics;
    float* cand_freqs;
    float* cand_amps[3];
    float* cand_grad_amps;
    float* cand_grad_amps_ax[3];
    int*   cand_violations;
    float* ana_peak_freqs;
    float* ana_peak_amps[3];
    float* ana_peak_phases[3];
    float* ana_peak_widths;
    float* surviving_freqs_hz;

    /* Component-level export arrays */
    int    max_component_terms;
    int    num_component_terms;
    float* component_freqs_hz;
    float* component_amps;
    float* component_phases_rad;
    float* component_widths_hz;
    int*   component_axes;
    int*   component_def_ids;
    int*   component_contrib_ids;
    int*   component_run_ids;

    (void)peak_log10_threshold;
    (void)peak_norm_scale;
    (void)peak_eps;
    (void)peak_prominence;

    memset(&se, 0, sizeof(se));
    eval_freqs = NULL;
    is_candidate = NULL;
    cand_freqs = NULL;
    n_tr_harmonics = 0;
    cand_grad_amps = NULL;
    cand_violations = NULL;
    ana_peak_freqs = NULL;
    ana_peak_widths = NULL;
    surviving_freqs_hz = NULL;
    num_component_terms = 0;
    max_component_terms = 0;
    component_freqs_hz = NULL;
    component_amps = NULL;
    component_phases_rad = NULL;
    component_widths_hz = NULL;
    component_axes = NULL;
    component_def_ids = NULL;
    component_contrib_ids = NULL;
    component_run_ids = NULL;
    for (ax = 0; ax < 3; ++ax) {
        coherent_mag[ax] = NULL;
        incoherent_mag[ax] = NULL;
        cand_amps[ax] = NULL;
        cand_grad_amps_ax[ax] = NULL;
        ana_peak_amps[ax] = NULL;
        ana_peak_phases[ax] = NULL;
    }

    result = sa_build_structural_events(&se, desc, start_block, block_count);
    if (PULSEQLIB_FAILED(result)) return result;

    /* Compute evaluation frequency upper bound */
    if (spectra->num_freq_bins <= 0) {
        sa_free_structural_events(&se);
        return PULSEQLIB_SUCCESS;
    }
    freq_max = spectra->freq_min_hz +
               (float)(spectra->num_freq_bins - 1) * spectra->freq_spacing_hz;

    /* --- Compute FFT power gate threshold ---
     * Find peak RSS magnitude of the dense FFT and set the gate at
     * SA_FFT_GATE_FRACTION of the peak.  Work in squared magnitudes to
     * avoid sqrt in the inner loops. */
    max_fft_rss_sq = 0.0f;
    has_fft_data = (spectra->spectrum_full_gx && spectra->spectrum_full_gy &&
                    spectra->spectrum_full_gz && spectra->num_freq_bins > 0 &&
                    spectra->freq_spacing_hz > 0.0f);
    if (has_fft_data) {
        for (i = 0; i < spectra->num_freq_bins; ++i) {
            float sx = spectra->spectrum_full_gx[i];
            float sy = spectra->spectrum_full_gy[i];
            float sz = spectra->spectrum_full_gz[i];
            float rss_sq = sx * sx + sy * sy + sz * sz;
            if (rss_sq > max_fft_rss_sq) max_fft_rss_sq = rss_sq;
        }
    }
    fft_gate_sq = SA_FFT_GATE_FRACTION * SA_FFT_GATE_FRACTION * max_fft_rss_sq;
    {   /* Promotion threshold (higher than gate) */
        float promo_frac = SA_FFT_PROMOTION_FRACTION;
        fft_promo_sq = promo_frac * promo_frac * max_fft_rss_sq;
    }

    /* --- Build evaluation frequency list ---
     * TR harmonics form the primary grid.  Additionally, any local maximum
     * of the dense FFT RSS spectrum that exceeds the gate threshold and
     * falls between TR harmonics is promoted into the evaluation set.
     * This catches spectral features (e.g. the bipolar-readout fundamental
     * in bSSFP) that land between adjacent TR harmonics.
     *
     * For K=1 we pretend K=2 so the TR-harmonic comb naturally produces
     * the correct spectrum via destructive/constructive interference. */
    if (num_instances < 2) num_instances = 2;
    {
        float tr_f1 = 1.0e6f / tr_duration_us;
        int m_max = (int)(freq_max / tr_f1);
        int max_eval, m, fi;
        if (m_max < 1) m_max = 0;

        /* Upper bound on eval_freqs size */
        max_eval = m_max + (has_fft_data ? spectra->num_freq_bins : 0);
        if (max_eval < 1) {
            sa_free_structural_events(&se);
            return PULSEQLIB_SUCCESS;
        }
        eval_freqs = (float*)PULSEQLIB_ALLOC((size_t)max_eval * sizeof(float));
        if (!eval_freqs) goto alloc_fail;

        /* 1) TR harmonics */
        num_freqs = 0;
        for (m = 0; m < m_max; ++m)
            eval_freqs[num_freqs++] = (float)(m + 1) * tr_f1;
        n_tr_harmonics = num_freqs;

        /* 2) FFT local maxima between TR harmonics (peak promotion) */
        if (has_fft_data && m_max > 0) {
            for (fi = 1; fi < spectra->num_freq_bins - 1; ++fi) {
                float sx = spectra->spectrum_full_gx[fi];
                float sy = spectra->spectrum_full_gy[fi];
                float sz = spectra->spectrum_full_gz[fi];
                float rss_sq = sx * sx + sy * sy + sz * sz;
                float f_peak, nearest, dist;
                float spx, spy, spz, snx, sny, snz;
                float rss_prev, rss_next;

                if (rss_sq < fft_promo_sq) continue;

                /* Local max check on RSS² */
                spx = spectra->spectrum_full_gx[fi - 1];
                spy = spectra->spectrum_full_gy[fi - 1];
                spz = spectra->spectrum_full_gz[fi - 1];
                snx = spectra->spectrum_full_gx[fi + 1];
                sny = spectra->spectrum_full_gy[fi + 1];
                snz = spectra->spectrum_full_gz[fi + 1];
                rss_prev = spx * spx + spy * spy + spz * spz;
                rss_next = snx * snx + sny * sny + snz * snz;
                if (rss_sq <= rss_prev || rss_sq <= rss_next) continue;

                /* Not near a TR harmonic? (> 25% of harmonic spacing) */
                f_peak = spectra->freq_min_hz
                       + (float)fi * spectra->freq_spacing_hz;
                if (f_peak <= 0.0f || f_peak > freq_max) continue;
                nearest = (float)floor((double)(f_peak / tr_f1) + 0.5) * tr_f1;
                dist = f_peak - nearest;
                if (dist < 0.0f) dist = -dist;
                if (dist <= tr_f1 * 0.25f) continue;

                eval_freqs[num_freqs++] = f_peak;
            }
        }

        if (num_freqs == 0) {
            PULSEQLIB_FREE(eval_freqs);
            eval_freqs = NULL;
            sa_free_structural_events(&se);
            return PULSEQLIB_SUCCESS;
        }
    }

    /* --- Evaluate inner spectrum at each frequency for each axis --- */
    for (ax = 0; ax < 3; ++ax) {
        coherent_mag[ax] = (float*)PULSEQLIB_ALLOC((size_t)num_freqs * sizeof(float));
        incoherent_mag[ax] = (float*)PULSEQLIB_ALLOC((size_t)num_freqs * sizeof(float));
        if (!coherent_mag[ax] || !incoherent_mag[ax]) goto alloc_fail;
    }

    for (i = 0; i < num_freqs; ++i) {
        f_hz = eval_freqs[i];
        for (ax = 0; ax < 3; ++ax) {
            if (se.axes[ax].num_events == 0) {
                coherent_mag[ax][i] = 0.0f;
                incoherent_mag[ax][i] = 0.0f;
            } else {
                sa_eval_inner_spectrum(f_hz, &se.axes[ax],
                    &coherent_mag[ax][i], &incoherent_mag[ax][i]);
            }
        }
    }

    /* --- Two-tier candidate selection --- */
    is_candidate = (int*)PULSEQLIB_ALLOC((size_t)num_freqs * sizeof(int));
    if (!is_candidate) goto alloc_fail;

    for (i = 0; i < num_freqs; ++i) {
        is_candidate[i] = 0;

        /* FFT power gate: reject frequencies where the actual waveform FFT
         * shows negligible energy, regardless of analytical coherence. */
        if (has_fft_data) {
            int fft_bin = (int)((eval_freqs[i] - spectra->freq_min_hz)
                                / spectra->freq_spacing_hz + 0.5f);
            if (fft_bin >= 0 && fft_bin < spectra->num_freq_bins) {
                float sx = spectra->spectrum_full_gx[fft_bin];
                float sy = spectra->spectrum_full_gy[fft_bin];
                float sz = spectra->spectrum_full_gz[fft_bin];
                float rss_sq = sx * sx + sy * sy + sz * sz;
                if (rss_sq < fft_gate_sq) continue;
            }
        }

        /* Tier 3: FFT-promoted peaks (indices >= n_tr_harmonics) are
         * automatically candidates — the dense FFT already confirmed
         * they are significant spectral features. */
        if (i >= n_tr_harmonics) {
            is_candidate[i] = 1;
            continue;
        }

        for (ax = 0; ax < 3; ++ax) {
            float coh = coherent_mag[ax][i];
            float inc = incoherent_mag[ax][i];

            if (se.axes[ax].num_events == 0)
                continue;

            /* Tier 1: coherence-detected (N_events >= 2) */
            if (se.axes[ax].num_events >= 2 && inc > 1.0e-20f) {
                float rho = coh / inc;
                if (rho > SA_COHERENCE_RATIO_THRESHOLD && coh > SA_AMPLITUDE_FLOOR) {
                    is_candidate[i] = 1;
                    break;
                }
            }

            /* Tier 2: absolute amplitude threshold (single-event axes) */
            if (se.axes[ax].num_events == 1 && coh > SA_AMPLITUDE_SAFETY) {
                is_candidate[i] = 1;
                break;
            }
        }
    }

    /* Count candidates */
    num_candidates = 0;
    for (i = 0; i < num_freqs; ++i)
        if (is_candidate[i]) num_candidates++;

    /* --- Estimate max component terms for export --- */
    max_component_terms = 0;
    for (ax = 0; ax < 3; ++ax)
        max_component_terms += se.axes[ax].num_events;
    /* Each event contributes one component term per candidate frequency,
     * but we export terms only for candidate frequencies.  Upper bound: */
    if (num_candidates > 0 && max_component_terms > 0)
        max_component_terms = max_component_terms * num_candidates;
    if (max_component_terms > 100000) max_component_terms = 100000;

    if (max_component_terms > 0) {
        component_freqs_hz = (float*)PULSEQLIB_ALLOC((size_t)max_component_terms * sizeof(float));
        component_amps = (float*)PULSEQLIB_ALLOC((size_t)max_component_terms * sizeof(float));
        component_phases_rad = (float*)PULSEQLIB_ALLOC((size_t)max_component_terms * sizeof(float));
        component_widths_hz = (float*)PULSEQLIB_ALLOC((size_t)max_component_terms * sizeof(float));
        component_axes = (int*)PULSEQLIB_ALLOC((size_t)max_component_terms * sizeof(int));
        component_def_ids = (int*)PULSEQLIB_ALLOC((size_t)max_component_terms * sizeof(int));
        component_contrib_ids = (int*)PULSEQLIB_ALLOC((size_t)max_component_terms * sizeof(int));
        component_run_ids = (int*)PULSEQLIB_ALLOC((size_t)max_component_terms * sizeof(int));
        if (!component_freqs_hz || !component_amps || !component_phases_rad ||
            !component_widths_hz || !component_axes || !component_def_ids ||
            !component_contrib_ids || !component_run_ids)
            goto alloc_fail;
    }

    /* --- Allocate candidate arrays --- */
    if (num_candidates > 0) {
        cand_freqs = (float*)PULSEQLIB_ALLOC((size_t)num_candidates * sizeof(float));
        cand_grad_amps = (float*)PULSEQLIB_ALLOC((size_t)num_candidates * sizeof(float));
        cand_violations = (int*)PULSEQLIB_ALLOC((size_t)num_candidates * sizeof(int));
        ana_peak_freqs = (float*)PULSEQLIB_ALLOC((size_t)num_candidates * sizeof(float));
        ana_peak_widths = (float*)PULSEQLIB_ALLOC((size_t)num_candidates * sizeof(float));
        surviving_freqs_hz = (float*)PULSEQLIB_ALLOC((size_t)num_candidates * sizeof(float));
        if (!cand_freqs || !cand_grad_amps || !cand_violations ||
            !ana_peak_freqs || !ana_peak_widths || !surviving_freqs_hz)
            goto alloc_fail;
        for (ax = 0; ax < 3; ++ax) {
            cand_amps[ax] = (float*)PULSEQLIB_ALLOC((size_t)num_candidates * sizeof(float));
            cand_grad_amps_ax[ax] = (float*)PULSEQLIB_ALLOC((size_t)num_candidates * sizeof(float));
            ana_peak_amps[ax] = (float*)PULSEQLIB_ALLOC((size_t)num_candidates * sizeof(float));
            ana_peak_phases[ax] = (float*)PULSEQLIB_ALLOC((size_t)num_candidates * sizeof(float));
            if (!cand_amps[ax] || !cand_grad_amps_ax[ax] ||
                !ana_peak_amps[ax] || !ana_peak_phases[ax])
                goto alloc_fail;
        }

        /* --- Fill candidate arrays --- */
        ci = 0;
        for (i = 0; i < num_freqs; ++i) {
            float max_ga;
            if (!is_candidate[i]) continue;
            f_hz = eval_freqs[i];
            cand_freqs[ci] = f_hz;
            surviving_freqs_hz[ci] = f_hz;
            ana_peak_freqs[ci] = f_hz;

            /* FWHM estimate */
            if (num_instances > 1 && tr_duration_us > 0.0f) {
                float f_rep = 1.0e6f / tr_duration_us;
                ana_peak_widths[ci] = 1.2067091288032284f * (f_rep / (float)num_instances);
            } else {
                ana_peak_widths[ci] = 1.0f;  /* Placeholder for K=1 */
            }

            /* Per-axis analytical amplitudes and phases */
            max_ga = 0.0f;
            for (ax = 0; ax < 3; ++ax) {
                float coh = coherent_mag[ax][i];
                cand_amps[ax][ci] = coh;
                ana_peak_amps[ax][ci] = coh;

                /* Phase: compute from coherent sum */
                {
                    float sum_re = 0.0f, sum_im = 0.0f;
                    int k;
                    for (k = 0; k < se.axes[ax].num_events; ++k) {
                        float d_re, d_im;
                        sa_eval_event_spectrum(f_hz, &se.axes[ax].events[k], &d_re, &d_im);
                        sum_re += d_re;
                        sum_im += d_im;
                    }
                    ana_peak_phases[ax][ci] = (float)atan2((double)sum_im, (double)sum_re);
                }

                /* G_eff per axis */
                if (se.axes[ax].num_events > 0) {
                    float ga = sa_eval_geff(f_hz, &se.axes[ax]);
                    cand_grad_amps_ax[ax][ci] = ga;
                    if (ga > max_ga) max_ga = ga;
                } else {
                    cand_grad_amps_ax[ax][ci] = 0.0f;
                }
            }

            /* Cross-axis G_eff: spectral-weighted blend across axes */
            {
                float aw_sum = 0.0f, aw_amp_sum = 0.0f;
                for (ax = 0; ax < 3; ++ax) {
                    float s_ax = cand_amps[ax][ci];
                    float a_ax = cand_grad_amps_ax[ax][ci];
                    if (s_ax > 0.0f && a_ax > 0.0f) {
                        aw_sum += s_ax;
                        aw_amp_sum += s_ax * a_ax;
                    }
                }
                if (aw_sum > 0.0f) max_ga = aw_amp_sum / aw_sum;
            }
            cand_grad_amps[ci] = max_ga;

            /* Check against forbidden bands */
            cand_violations[ci] = 0;
            for (b = 0; b < num_forbidden_bands; ++b) {
                if (f_hz >= forbidden_bands[b].freq_min_hz &&
                    f_hz <= forbidden_bands[b].freq_max_hz) {
                    if (max_ga > forbidden_bands[b].max_amplitude_hz_per_m)
                        cand_violations[ci] = 1;
                }
            }

            /* Export component terms for this candidate frequency */
            for (ax = 0; ax < 3; ++ax) {
                int k;
                for (k = 0; k < se.axes[ax].num_events; ++k) {
                    float d_re, d_im, d_mag;
                    if (num_component_terms >= max_component_terms) break;
                    sa_eval_event_spectrum(f_hz, &se.axes[ax].events[k], &d_re, &d_im);
                    d_mag = (float)sqrt((double)(d_re * d_re + d_im * d_im));
                    component_freqs_hz[num_component_terms] = f_hz;
                    component_amps[num_component_terms] = d_mag;
                    component_phases_rad[num_component_terms] = (float)atan2((double)d_im, (double)d_re);
                    component_widths_hz[num_component_terms] = ana_peak_widths[ci];
                    component_axes[num_component_terms] = ax;
                    component_def_ids[num_component_terms] = se.axes[ax].events[k].def_id;
                    component_contrib_ids[num_component_terms] = k;
                    component_run_ids[num_component_terms] = 0;
                    num_component_terms++;
                }
            }

            ci++;
        }
    }

    /* --- Assign to output spectra struct (ownership transfer) --- */
    spectra->num_analytical_peaks   = num_candidates;
    spectra->analytical_peak_freqs  = ana_peak_freqs;     ana_peak_freqs = NULL;
    spectra->analytical_peak_amp_gx = ana_peak_amps[0];   ana_peak_amps[0] = NULL;
    spectra->analytical_peak_amp_gy = ana_peak_amps[1];   ana_peak_amps[1] = NULL;
    spectra->analytical_peak_amp_gz = ana_peak_amps[2];   ana_peak_amps[2] = NULL;
    spectra->analytical_peak_phase_gx = ana_peak_phases[0]; ana_peak_phases[0] = NULL;
    spectra->analytical_peak_phase_gy = ana_peak_phases[1]; ana_peak_phases[1] = NULL;
    spectra->analytical_peak_phase_gz = ana_peak_phases[2]; ana_peak_phases[2] = NULL;
    spectra->analytical_peak_widths_hz = ana_peak_widths; ana_peak_widths = NULL;

    spectra->num_candidates            = num_candidates;
    spectra->candidate_freqs           = cand_freqs;          cand_freqs = NULL;
    spectra->candidate_amps_gx         = cand_amps[0];        cand_amps[0] = NULL;
    spectra->candidate_amps_gy         = cand_amps[1];        cand_amps[1] = NULL;
    spectra->candidate_amps_gz         = cand_amps[2];        cand_amps[2] = NULL;
    spectra->candidate_grad_amps       = cand_grad_amps;      cand_grad_amps = NULL;
    spectra->candidate_grad_amps_gx    = cand_grad_amps_ax[0]; cand_grad_amps_ax[0] = NULL;
    spectra->candidate_grad_amps_gy    = cand_grad_amps_ax[1]; cand_grad_amps_ax[1] = NULL;
    spectra->candidate_grad_amps_gz    = cand_grad_amps_ax[2]; cand_grad_amps_ax[2] = NULL;
    spectra->candidate_violations      = cand_violations;     cand_violations = NULL;
    spectra->num_component_terms       = num_component_terms;
    spectra->component_freqs_hz        = component_freqs_hz; component_freqs_hz = NULL;
    spectra->component_amps            = component_amps; component_amps = NULL;
    spectra->component_phases_rad      = component_phases_rad; component_phases_rad = NULL;
    spectra->component_widths_hz       = component_widths_hz; component_widths_hz = NULL;
    spectra->component_axes            = component_axes; component_axes = NULL;
    spectra->component_def_ids         = component_def_ids; component_def_ids = NULL;
    spectra->component_contrib_ids     = component_contrib_ids; component_contrib_ids = NULL;
    spectra->component_run_ids         = component_run_ids; component_run_ids = NULL;
    spectra->num_surviving_freqs       = num_candidates;
    spectra->surviving_freqs_hz        = surviving_freqs_hz; surviving_freqs_hz = NULL;

    /* Cleanup */
    if (eval_freqs)    PULSEQLIB_FREE(eval_freqs);
    if (is_candidate)  PULSEQLIB_FREE(is_candidate);
    for (ax = 0; ax < 3; ++ax) {
        if (coherent_mag[ax])    PULSEQLIB_FREE(coherent_mag[ax]);
        if (incoherent_mag[ax])  PULSEQLIB_FREE(incoherent_mag[ax]);
    }
    sa_free_structural_events(&se);
    return PULSEQLIB_SUCCESS;

alloc_fail:
    for (ax = 0; ax < 3; ++ax) {
        if (coherent_mag[ax])     PULSEQLIB_FREE(coherent_mag[ax]);
        if (incoherent_mag[ax])   PULSEQLIB_FREE(incoherent_mag[ax]);
        if (cand_amps[ax])        PULSEQLIB_FREE(cand_amps[ax]);
        if (cand_grad_amps_ax[ax]) PULSEQLIB_FREE(cand_grad_amps_ax[ax]);
        if (ana_peak_amps[ax])    PULSEQLIB_FREE(ana_peak_amps[ax]);
        if (ana_peak_phases[ax])  PULSEQLIB_FREE(ana_peak_phases[ax]);
    }
    if (ana_peak_freqs)       PULSEQLIB_FREE(ana_peak_freqs);
    if (ana_peak_widths)      PULSEQLIB_FREE(ana_peak_widths);
    if (eval_freqs)           PULSEQLIB_FREE(eval_freqs);
    if (is_candidate)         PULSEQLIB_FREE(is_candidate);
    if (cand_freqs)           PULSEQLIB_FREE(cand_freqs);
    if (cand_grad_amps)       PULSEQLIB_FREE(cand_grad_amps);
    if (cand_violations)      PULSEQLIB_FREE(cand_violations);
    if (component_freqs_hz)   PULSEQLIB_FREE(component_freqs_hz);
    if (component_amps)       PULSEQLIB_FREE(component_amps);
    if (component_phases_rad) PULSEQLIB_FREE(component_phases_rad);
    if (component_widths_hz)  PULSEQLIB_FREE(component_widths_hz);
    if (component_axes)       PULSEQLIB_FREE(component_axes);
    if (component_def_ids)    PULSEQLIB_FREE(component_def_ids);
    if (component_contrib_ids) PULSEQLIB_FREE(component_contrib_ids);
    if (component_run_ids)    PULSEQLIB_FREE(component_run_ids);
    if (surviving_freqs_hz)   PULSEQLIB_FREE(surviving_freqs_hz);
    sa_free_structural_events(&se);
    return PULSEQLIB_ERR_ALLOC_FAILED;
}

/* ================================================================== */
/*  Acoustic violation check                                          */
/* ================================================================== */

static int check_acoustic_violations(
    int** out_peaks,
    const float* spectrum, const float* frequencies, int num_freq_bins,
    float max_envelope,
    const pulseqlib_forbidden_band* bands, int num_bands,
    float peak_log10_threshold,
    float peak_norm_scale,
    float peak_eps,
    float peak_prominence)
{
    int* peaks;
    int i, b;
    float freq;

    (void)max_envelope;

    if (num_bands <= 0 || !bands) {
        if (out_peaks) *out_peaks = NULL;
        return PULSEQLIB_SUCCESS;
    }

    peaks = (int*)PULSEQLIB_ALLOC((size_t)num_freq_bins * sizeof(int));
    if (!peaks) return PULSEQLIB_ERR_ALLOC_FAILED;

    detect_resonances(
        peaks,
        spectrum,
        num_freq_bins,
        peak_log10_threshold,
        peak_norm_scale,
        peak_eps,
        peak_prominence);

    /* Check peaks against forbidden bands (informational) */
    for (i = 0; i < num_freq_bins; ++i) {
        if (!peaks[i]) continue;
        if (!frequencies) continue;
        freq = frequencies[i];
        for (b = 0; b < num_bands; ++b) {
            if (freq >= bands[b].freq_min_hz && freq <= bands[b].freq_max_hz)
                break;
        }
    }

    if (out_peaks) *out_peaks = peaks;
    else           PULSEQLIB_FREE(peaks);

    return PULSEQLIB_SUCCESS;
}

/* ================================================================== */
/*  Sliding window spectra                                            */
/* ================================================================== */

/* FIX: outputs grouped before inputs, out_peaks + out_max_envelope moved up */
static int compute_sliding_window_spectra(
    acoustic_support* sup,
    float* spectra_out,
    float* out_max_envelope,
    int* out_peaks,
    const float* waveform, const float* frequencies,
    int num_samples, int padded_len, int combined,
    const pulseqlib_forbidden_band* bands, int num_bands,
    float peak_log10_threshold,
    float peak_norm_scale,
    float peak_eps,
    float peak_prominence)
{
    acoustic_waveform aw;
    int w, i, result, start_idx, window_len;
    float* win_spectrum;
    float max_env_win, abs_val, max_env_overall;
    int* win_peaks;

    win_spectrum = NULL;
    memset(&aw, 0, sizeof(aw));
    max_env_overall = 0.0f;

    if (combined) {
        win_spectrum = (float*)PULSEQLIB_ALLOC(sup->output_freq_bins * sizeof(float));
        if (!win_spectrum) return PULSEQLIB_ERR_ALLOC_FAILED;
        for (i = 0; i < sup->output_freq_bins; ++i) spectra_out[i] = 0.0f;
    }

    result = acoustic_waveform_init(&aw, sup, waveform, num_samples, padded_len);
    if (PULSEQLIB_FAILED(result)) {
        if (win_spectrum) PULSEQLIB_FREE(win_spectrum);
        return result;
    }

    for (w = 0; w < sup->num_windows; ++w) {
        start_idx  = w * sup->hop_size;
        window_len = sup->nwin;
        if (start_idx + window_len > aw.num_samples)
            window_len = aw.num_samples - start_idx;

        max_env_win = 0.0f;
        for (i = start_idx; i < start_idx + window_len; ++i) {
            abs_val = (aw.samples[i] >= 0.0f) ? aw.samples[i] : -aw.samples[i];
            if (abs_val > max_env_win) max_env_win = abs_val;
        }

        if (combined) {
            if (max_env_win > max_env_overall) max_env_overall = max_env_win;
            result = compute_window_spectrum(sup, win_spectrum, &aw, w);
            if (PULSEQLIB_FAILED(result)) {
                PULSEQLIB_FREE(win_spectrum); acoustic_waveform_free(&aw); return result;
            }
            for (i = 0; i < sup->output_freq_bins; ++i)
                if (win_spectrum[i] > spectra_out[i]) spectra_out[i] = win_spectrum[i];
        } else {
            if (out_max_envelope) out_max_envelope[w] = max_env_win;
            result = compute_window_spectrum(sup,
                &spectra_out[w * sup->output_freq_bins], &aw, w);
            if (PULSEQLIB_FAILED(result)) {
                if (win_spectrum) PULSEQLIB_FREE(win_spectrum);
                acoustic_waveform_free(&aw); return result;
            }

            if (num_bands > 0 && bands) {
                win_peaks = NULL;
                result = check_acoustic_violations(out_peaks ? &win_peaks : NULL,
                    &spectra_out[w * sup->output_freq_bins], frequencies,
                    sup->output_freq_bins, max_env_win, bands, num_bands,
                    peak_log10_threshold, peak_norm_scale, peak_eps,
                    peak_prominence);
                if (PULSEQLIB_FAILED(result)) {
                    if (win_spectrum) PULSEQLIB_FREE(win_spectrum);
                    acoustic_waveform_free(&aw); return result;
                }
                if (win_peaks && out_peaks) {
                    memcpy(&out_peaks[w * sup->output_freq_bins], win_peaks,
                           sup->output_freq_bins * sizeof(int));
                    PULSEQLIB_FREE(win_peaks);
                }
            } else if (out_peaks) {
                detect_resonances(
                    &out_peaks[w * sup->output_freq_bins],
                    &spectra_out[w * sup->output_freq_bins],
                    sup->output_freq_bins,
                    peak_log10_threshold,
                    peak_norm_scale,
                    peak_eps,
                    peak_prominence);
            }
        }
    }

    if (combined && out_max_envelope) out_max_envelope[0] = max_env_overall;

    if (win_spectrum) PULSEQLIB_FREE(win_spectrum);
    acoustic_waveform_free(&aw);
    return PULSEQLIB_SUCCESS;
}

/* ================================================================== */
/*  Sequence spectrum (full-TR FFT + harmonic sampling)               */
/* ================================================================== */

/* FIX: outputs grouped before inputs */
static int compute_sequence_spectrum(
    float* full_spectrum,
    float** seq_spectrum, float** seq_frequencies,
    int* out_num_freq_bins_full, float* out_freq_res_full,
    int* out_num_picked, float* out_max_envelope,
    int** out_seq_peaks,
    const float* waveform, int num_samples,
    float grad_raster_us, float target_spectral_res_hz,
    float max_frequency, float fundamental_freq,
    const pulseqlib_forbidden_band* bands, int num_bands,
    float peak_log10_threshold,
    float peak_norm_scale,
    float peak_eps,
    float peak_prominence)
{
    int nfft, nfreq, output_bins_full, num_picked;
    int min_nfft, max_idx;
    float freq_res, max_freq, max_env, abs_val;
    float* work;
    float* cos_win;
    kiss_fft_cpx* fft_out;
    kiss_fftr_cfg cfg;
    float* picked_mag;
    float* picked_freq;
    float mean, fft_norm;
    int i, k, freq_idx;
    float freq, freq_low, freq_high, t;
    float re_low, im_low, re_high, im_high, re_interp, im_interp;
    int* seq_peaks;
    int result;

    result      = PULSEQLIB_SUCCESS;
    work        = NULL;
    cos_win     = NULL;
    fft_out     = NULL;
    cfg         = NULL;
    picked_mag  = NULL;
    picked_freq = NULL;

    if (!waveform || num_samples <= 0) return PULSEQLIB_ERR_INVALID_ARGUMENT;

    max_env = 0.0f;
    for (i = 0; i < num_samples; ++i) {
        abs_val = (waveform[i] >= 0.0f) ? waveform[i] : -waveform[i];
        if (abs_val > max_env) max_env = abs_val;
    }
    if (out_max_envelope) *out_max_envelope = max_env;

    min_nfft = (int)ceil((double)(1.0e6 / (grad_raster_us * target_spectral_res_hz)));
    nfft = (min_nfft < num_samples) ? (int)pulseqlib__next_pow2((size_t)num_samples)
                                    : (int)pulseqlib__next_pow2((size_t)min_nfft);
    nfreq = nfft / 2 + 1;
    freq_res = (float)(1.0e6 / (grad_raster_us * (double)nfft));
    max_freq = (max_frequency > 0.0f) ? max_frequency : (float)(5.0e5 / grad_raster_us);
    max_idx  = (int)(max_freq / freq_res + 0.5f);
    if (max_idx >= nfreq) output_bins_full = nfreq;
    else if (max_idx < 1) output_bins_full = 1;
    else                  output_bins_full = max_idx + 1;

    work    = (float*)PULSEQLIB_ALLOC((size_t)nfft * sizeof(float));
    cos_win = (float*)PULSEQLIB_ALLOC((size_t)num_samples * sizeof(float));
    fft_out = (kiss_fft_cpx*)PULSEQLIB_ALLOC((size_t)nfreq * sizeof(kiss_fft_cpx));
    if (!work || !cos_win || !fft_out) { result = PULSEQLIB_ERR_ALLOC_FAILED; goto fail; }

    cfg = kiss_fftr_alloc(nfft, 0, NULL, NULL);
    if (!cfg) { result = PULSEQLIB_ERR_ALLOC_FAILED; goto fail; }

    for (i = 0; i < num_samples; ++i)
        cos_win[i] = (float)(0.5 * (1.0 - cos(2.0 * M_PI * (double)(i + 1) / (double)num_samples)));

    for (i = 0; i < num_samples; ++i) work[i] = waveform[i];
    for (i = num_samples; i < nfft; ++i) work[i] = 0.0f;

    mean = 0.0f;
    for (i = 0; i < num_samples; ++i) mean += work[i];
    mean /= (float)num_samples;
    for (i = 0; i < num_samples; ++i) work[i] -= mean;
    for (i = 0; i < num_samples; ++i) work[i] *= cos_win[i];

    kiss_fftr(cfg, work, fft_out);
    fft_norm = 1.0f / (float)nfft;
    for (i = 0; i < nfreq; ++i) { fft_out[i].r *= fft_norm; fft_out[i].i *= fft_norm; }

    if (full_spectrum) {
        for (i = 0; i < output_bins_full; ++i)
            full_spectrum[i] = (float)sqrt((double)(fft_out[i].r * fft_out[i].r +
                                                    fft_out[i].i * fft_out[i].i));
    }

    if (fundamental_freq > 0.0f && seq_spectrum) {
        num_picked = (int)(max_freq / fundamental_freq) + 1;
        picked_mag = (float*)PULSEQLIB_ALLOC((size_t)num_picked * sizeof(float));
        if (!picked_mag) { result = PULSEQLIB_ERR_ALLOC_FAILED; goto fail; }

        if (seq_frequencies) {
            picked_freq = (float*)PULSEQLIB_ALLOC((size_t)num_picked * sizeof(float));
            if (!picked_freq) { result = PULSEQLIB_ERR_ALLOC_FAILED; goto fail; }
        }

        for (k = 0; k < num_picked; ++k) {
            freq = (float)k * fundamental_freq;
            if (picked_freq) picked_freq[k] = freq;

            freq_idx = (int)(freq / freq_res);
            if (freq_idx >= nfreq - 1) {
                picked_mag[k] = 0.0f;
            } else if (freq_idx == 0) {
                picked_mag[k] = (float)sqrt((double)(fft_out[0].r * fft_out[0].r +
                                                     fft_out[0].i * fft_out[0].i));
            } else {
                freq_low  = (float)freq_idx * freq_res;
                freq_high = (float)(freq_idx + 1) * freq_res;
                re_low  = fft_out[freq_idx].r;     im_low  = fft_out[freq_idx].i;
                re_high = fft_out[freq_idx + 1].r;  im_high = fft_out[freq_idx + 1].i;
                t = (freq - freq_low) / (freq_high - freq_low);
                re_interp = re_low * (1.0f - t) + re_high * t;
                im_interp = im_low * (1.0f - t) + im_high * t;
                picked_mag[k] = (float)sqrt((double)(re_interp * re_interp +
                                                     im_interp * im_interp));
            }
        }

        *seq_spectrum = picked_mag;   picked_mag = NULL;
        if (seq_frequencies) { *seq_frequencies = picked_freq; picked_freq = NULL; }
        if (out_num_picked) *out_num_picked = num_picked;
    }

    if (out_num_freq_bins_full) *out_num_freq_bins_full = output_bins_full;
    if (out_freq_res_full)      *out_freq_res_full      = freq_res;

    if (num_bands > 0 && bands && seq_spectrum && *seq_spectrum) {
        num_picked = out_num_picked ? *out_num_picked : 0;
        result = check_acoustic_violations(out_seq_peaks,
            *seq_spectrum, (seq_frequencies && *seq_frequencies) ? *seq_frequencies : NULL,
            num_picked, max_env, bands, num_bands,
            peak_log10_threshold, peak_norm_scale, peak_eps,
            peak_prominence);
        if (PULSEQLIB_FAILED(result)) goto fail;
    } else if (out_seq_peaks && seq_spectrum && *seq_spectrum) {
        num_picked = out_num_picked ? *out_num_picked : 0;
        seq_peaks = (int*)PULSEQLIB_ALLOC((size_t)num_picked * sizeof(int));
        if (!seq_peaks) { result = PULSEQLIB_ERR_ALLOC_FAILED; goto fail; }
        detect_resonances(
            seq_peaks,
            *seq_spectrum,
            num_picked,
            peak_log10_threshold,
            peak_norm_scale,
            peak_eps,
            peak_prominence);
        *out_seq_peaks = seq_peaks;
    }

fail:
    if (work)        PULSEQLIB_FREE(work);
    if (cos_win)     PULSEQLIB_FREE(cos_win);
    if (fft_out)     PULSEQLIB_FREE(fft_out);
    if (cfg)         kiss_fftr_free(cfg);
    if (picked_mag)  PULSEQLIB_FREE(picked_mag);
    if (picked_freq) PULSEQLIB_FREE(picked_freq);
    return result;
}

/* ================================================================== */
/*  Acoustic spectra (static helper from uniform waveforms)           */
/* ================================================================== */

static int calc_mech_resonances_from_uniform(
    pulseqlib_mech_resonances_spectra* spectra,
    pulseqlib_diagnostic* diag,
    const pulseqlib__uniform_grad_waveforms* waveforms,
    int target_window_size,
    float target_spectral_resolution_hz,
    float max_frequency_hz,
    int num_trs,
    float tr_duration_us,
    int num_forbidden_bands,
    const pulseqlib_forbidden_band* forbidden_bands,
    float peak_log10_threshold,
    float peak_norm_scale,
    float peak_eps,
    float peak_prominence,
    const struct pulseqlib_sequence_descriptor* desc,
    int start_block,
    int block_count)
{
    acoustic_support sup;
    pulseqlib_diagnostic local_diag;
    int max_samples, result, output_size, padded_len;
    float fundamental_freq;
    int num_freq_bins_full, num_freq_bins_seq;
    float freq_res_full;
    float* seq_spec_gx;
    float* seq_spec_gy;
    float* seq_spec_gz;
    float* seq_freqs;

    seq_spec_gx = NULL;
    seq_spec_gy = NULL;
    seq_spec_gz = NULL;
    seq_freqs   = NULL;

    if (!diag) { pulseqlib_diagnostic_init(&local_diag); diag = &local_diag; }
    else       { pulseqlib_diagnostic_init(diag); }

    if (!waveforms || !spectra) {
        diag->code = PULSEQLIB_ERR_NULL_POINTER; return diag->code;
    }

    memset(spectra, 0, sizeof(*spectra));
    memset(&sup, 0, sizeof(sup));

    max_samples = waveforms->num_samples;
    if (max_samples <= 0) { diag->code = PULSEQLIB_ERR_MECH_RESONANCES_NO_WAVEFORM; return diag->code; }

    result = acoustic_support_init(&sup, max_samples, target_window_size,
                                   target_spectral_resolution_hz,
                                   waveforms->raster_us, max_frequency_hz);
    if (PULSEQLIB_FAILED(result)) { diag->code = result; return result; }

    spectra->freq_min_hz    = 0.0f;
    spectra->freq_spacing_hz = sup.freq_resolution;
    spectra->num_freq_bins  = sup.output_freq_bins;
    spectra->num_windows    = sup.num_windows;
    spectra->num_instances  = num_trs;

    output_size = sup.num_windows * sup.output_freq_bins;

    spectra->spectrogram_gx = (float*)PULSEQLIB_ALLOC((size_t)output_size * sizeof(float));
    spectra->spectrogram_gy = (float*)PULSEQLIB_ALLOC((size_t)output_size * sizeof(float));
    spectra->spectrogram_gz = (float*)PULSEQLIB_ALLOC((size_t)output_size * sizeof(float));
    spectra->peaks_gx = (int*)PULSEQLIB_ALLOC((size_t)output_size * sizeof(int));
    spectra->peaks_gy = (int*)PULSEQLIB_ALLOC((size_t)output_size * sizeof(int));
    spectra->peaks_gz = (int*)PULSEQLIB_ALLOC((size_t)output_size * sizeof(int));

    if (!spectra->spectrogram_gx || !spectra->spectrogram_gy || !spectra->spectrogram_gz ||
        !spectra->peaks_gx || !spectra->peaks_gy || !spectra->peaks_gz) {
        pulseqlib_mech_resonances_spectra_free(spectra);
        acoustic_support_free(&sup);
        diag->code = PULSEQLIB_ERR_ALLOC_FAILED; return diag->code;
    }

    memset(spectra->spectrogram_gx, 0, (size_t)output_size * sizeof(float));
    memset(spectra->spectrogram_gy, 0, (size_t)output_size * sizeof(float));
    memset(spectra->spectrogram_gz, 0, (size_t)output_size * sizeof(float));
    memset(spectra->peaks_gx, 0, (size_t)output_size * sizeof(int));
    memset(spectra->peaks_gy, 0, (size_t)output_size * sizeof(int));
    memset(spectra->peaks_gz, 0, (size_t)output_size * sizeof(int));

    /* padded length */
    padded_len = (max_samples <= sup.nwin) ? max_samples
                 : ((max_samples + sup.nwin - 1) / sup.nwin) * sup.nwin;

    /* sliding window: Gx */
    if (waveforms->num_samples > 0) {
        result = compute_sliding_window_spectra(&sup,
            spectra->spectrogram_gx, NULL, spectra->peaks_gx,
            waveforms->gx, NULL,
            waveforms->num_samples, padded_len, 0,
            forbidden_bands, num_forbidden_bands,
            peak_log10_threshold, peak_norm_scale, peak_eps,
            peak_prominence);
        if (PULSEQLIB_FAILED(result)) {
            pulseqlib_mech_resonances_spectra_free(spectra);
            acoustic_support_free(&sup);
            diag->code = result; return result;
        }
    }
    /* sliding window: Gy */
    if (waveforms->num_samples > 0) {
        result = compute_sliding_window_spectra(&sup,
            spectra->spectrogram_gy, NULL, spectra->peaks_gy,
            waveforms->gy, NULL,
            waveforms->num_samples, padded_len, 0,
            forbidden_bands, num_forbidden_bands,
            peak_log10_threshold, peak_norm_scale, peak_eps,
            peak_prominence);
        if (PULSEQLIB_FAILED(result)) {
            pulseqlib_mech_resonances_spectra_free(spectra);
            acoustic_support_free(&sup);
            diag->code = result; return result;
        }
    }
    /* sliding window: Gz */
    if (waveforms->num_samples > 0) {
        result = compute_sliding_window_spectra(&sup,
            spectra->spectrogram_gz, NULL, spectra->peaks_gz,
            waveforms->gz, NULL,
            waveforms->num_samples, padded_len, 0,
            forbidden_bands, num_forbidden_bands,
            peak_log10_threshold, peak_norm_scale, peak_eps,
            peak_prominence);
        if (PULSEQLIB_FAILED(result)) {
            pulseqlib_mech_resonances_spectra_free(spectra);
            acoustic_support_free(&sup);
            diag->code = result; return result;
        }
    }

    acoustic_support_free(&sup);

    /* ---- sequence spectra ---- */
    if (num_trs > 1 && tr_duration_us > 0.0f) {
        fundamental_freq = 1.0e6f / tr_duration_us;
    } else {
        fundamental_freq = 0.0f;
    }

    /* Gx (first call: determines sizes) */
    if (waveforms->num_samples > 0) {
        result = compute_sequence_spectrum(NULL,
            (fundamental_freq > 0.0f) ? &seq_spec_gx : NULL,
            (fundamental_freq > 0.0f) ? &seq_freqs   : NULL,
            &num_freq_bins_full, &freq_res_full, &num_freq_bins_seq, NULL,
            (fundamental_freq > 0.0f) ? &spectra->peaks_seq_gx : NULL,
            waveforms->gx, waveforms->num_samples,
            waveforms->raster_us, target_spectral_resolution_hz,
            max_frequency_hz, fundamental_freq,
            forbidden_bands, num_forbidden_bands,
            peak_log10_threshold, peak_norm_scale, peak_eps,
            peak_prominence);
        if (PULSEQLIB_FAILED(result)) {
            pulseqlib_mech_resonances_spectra_free(spectra);
            diag->code = result; return result;
        }
    } else {
        result = compute_sequence_spectrum(NULL, NULL, NULL,
            &num_freq_bins_full, &freq_res_full, NULL, NULL, NULL,
            waveforms->gy ? waveforms->gy : waveforms->gz,
            max_samples, waveforms->raster_us, target_spectral_resolution_hz,
            max_frequency_hz, 0.0f,
            forbidden_bands, num_forbidden_bands,
            peak_log10_threshold, peak_norm_scale, peak_eps,
            peak_prominence);
        if (PULSEQLIB_FAILED(result)) {
            pulseqlib_mech_resonances_spectra_free(spectra);
            diag->code = result; return result;
        }
    }

    /* full TR spectra allocation */
    spectra->spectrum_full_gx = (float*)PULSEQLIB_ALLOC((size_t)num_freq_bins_full * sizeof(float));
    spectra->spectrum_full_gy = (float*)PULSEQLIB_ALLOC((size_t)num_freq_bins_full * sizeof(float));
    spectra->spectrum_full_gz = (float*)PULSEQLIB_ALLOC((size_t)num_freq_bins_full * sizeof(float));
    if (!spectra->spectrum_full_gx || !spectra->spectrum_full_gy || !spectra->spectrum_full_gz) {
        if (seq_spec_gx) PULSEQLIB_FREE(seq_spec_gx);
        if (seq_freqs)   PULSEQLIB_FREE(seq_freqs);
        pulseqlib_mech_resonances_spectra_free(spectra);
        diag->code = PULSEQLIB_ERR_ALLOC_FAILED; return diag->code;
    }

    if (fundamental_freq > 0.0f && num_freq_bins_seq > 0) {
        spectra->freq_spacing_seq_hz = fundamental_freq;
        spectra->num_freq_bins_seq = num_freq_bins_seq;
        spectra->spectrum_seq_gx = (float*)PULSEQLIB_ALLOC((size_t)num_freq_bins_seq * sizeof(float));
        spectra->spectrum_seq_gy = (float*)PULSEQLIB_ALLOC((size_t)num_freq_bins_seq * sizeof(float));
        spectra->spectrum_seq_gz = (float*)PULSEQLIB_ALLOC((size_t)num_freq_bins_seq * sizeof(float));
        if (!spectra->spectrum_seq_gx || !spectra->spectrum_seq_gy || !spectra->spectrum_seq_gz) {
            if (seq_spec_gx) PULSEQLIB_FREE(seq_spec_gx);
            if (seq_freqs)   PULSEQLIB_FREE(seq_freqs);
            pulseqlib_mech_resonances_spectra_free(spectra);
            diag->code = PULSEQLIB_ERR_ALLOC_FAILED; return diag->code;
        }
        if (seq_spec_gx && waveforms->num_samples > 0) {
            memcpy(spectra->spectrum_seq_gx, seq_spec_gx, (size_t)num_freq_bins_seq * sizeof(float));
            PULSEQLIB_FREE(seq_spec_gx); seq_spec_gx = NULL;
        } else {
            memset(spectra->spectrum_seq_gx, 0, (size_t)num_freq_bins_seq * sizeof(float));
        }
    }
    if (seq_spec_gx) { PULSEQLIB_FREE(seq_spec_gx); seq_spec_gx = NULL; }
    if (seq_freqs)   { PULSEQLIB_FREE(seq_freqs);   seq_freqs = NULL; }

    /* full-TR spectrum: Gx */
    if (waveforms->num_samples > 0) {
        result = compute_sequence_spectrum(
            spectra->spectrum_full_gx, NULL, NULL,
            NULL, NULL, NULL, NULL, NULL,
            waveforms->gx, waveforms->num_samples,
            waveforms->raster_us, target_spectral_resolution_hz,
            max_frequency_hz, 0.0f,
            forbidden_bands, num_forbidden_bands,
            peak_log10_threshold, peak_norm_scale, peak_eps,
            peak_prominence);
        if (PULSEQLIB_FAILED(result)) {
            pulseqlib_mech_resonances_spectra_free(spectra);
            diag->code = result; return result;
        }
    } else {
        memset(spectra->spectrum_full_gx, 0, (size_t)num_freq_bins_full * sizeof(float));
    }

    /* Gy - full TR + seq */
    if (waveforms->num_samples > 0) {
        result = compute_sequence_spectrum(
            spectra->spectrum_full_gy,
            (fundamental_freq > 0.0f) ? &seq_spec_gy : NULL,
            NULL, NULL, NULL, NULL, NULL,
            (fundamental_freq > 0.0f) ? &spectra->peaks_seq_gy : NULL,
            waveforms->gy, waveforms->num_samples,
            waveforms->raster_us, target_spectral_resolution_hz,
            max_frequency_hz, fundamental_freq,
            forbidden_bands, num_forbidden_bands,
            peak_log10_threshold, peak_norm_scale, peak_eps,
            peak_prominence);
        if (PULSEQLIB_FAILED(result)) {
            pulseqlib_mech_resonances_spectra_free(spectra);
            diag->code = result; return result;
        }
        if (seq_spec_gy && spectra->spectrum_seq_gy) {
            memcpy(spectra->spectrum_seq_gy, seq_spec_gy, (size_t)num_freq_bins_seq * sizeof(float));
            PULSEQLIB_FREE(seq_spec_gy); seq_spec_gy = NULL;
        } else if (seq_spec_gy) { PULSEQLIB_FREE(seq_spec_gy); seq_spec_gy = NULL; }
    } else {
        memset(spectra->spectrum_full_gy, 0, (size_t)num_freq_bins_full * sizeof(float));
        if (spectra->spectrum_seq_gy)
            memset(spectra->spectrum_seq_gy, 0, (size_t)num_freq_bins_seq * sizeof(float));
        if (spectra->peaks_seq_gy)
            memset(spectra->peaks_seq_gy, 0, (size_t)num_freq_bins_seq * sizeof(int));
    }

    /* Gz - full TR + seq */
    if (waveforms->num_samples > 0) {
        result = compute_sequence_spectrum(
            spectra->spectrum_full_gz,
            (fundamental_freq > 0.0f) ? &seq_spec_gz : NULL,
            NULL, NULL, NULL, NULL, NULL,
            (fundamental_freq > 0.0f) ? &spectra->peaks_seq_gz : NULL,
            waveforms->gz, waveforms->num_samples,
            waveforms->raster_us, target_spectral_resolution_hz,
            max_frequency_hz, fundamental_freq,
            forbidden_bands, num_forbidden_bands,
            peak_log10_threshold, peak_norm_scale, peak_eps,
            peak_prominence);
        if (PULSEQLIB_FAILED(result)) {
            pulseqlib_mech_resonances_spectra_free(spectra);
            diag->code = result; return result;
        }
        if (seq_spec_gz && spectra->spectrum_seq_gz) {
            memcpy(spectra->spectrum_seq_gz, seq_spec_gz, (size_t)num_freq_bins_seq * sizeof(float));
            PULSEQLIB_FREE(seq_spec_gz); seq_spec_gz = NULL;
        } else if (seq_spec_gz) { PULSEQLIB_FREE(seq_spec_gz); seq_spec_gz = NULL; }
    } else {
        memset(spectra->spectrum_full_gz, 0, (size_t)num_freq_bins_full * sizeof(float));
        if (spectra->spectrum_seq_gz)
            memset(spectra->spectrum_seq_gz, 0, (size_t)num_freq_bins_seq * sizeof(float));
        if (spectra->peaks_seq_gz)
            memset(spectra->peaks_seq_gz, 0, (size_t)num_freq_bins_seq * sizeof(int));
    }

    /* full-TR peak detection */
    if (num_freq_bins_full > 0) {
        spectra->peaks_full_gx = (int*)PULSEQLIB_ALLOC((size_t)num_freq_bins_full * sizeof(int));
        spectra->peaks_full_gy = (int*)PULSEQLIB_ALLOC((size_t)num_freq_bins_full * sizeof(int));
        spectra->peaks_full_gz = (int*)PULSEQLIB_ALLOC((size_t)num_freq_bins_full * sizeof(int));
        if (spectra->peaks_full_gx) detect_resonances(spectra->peaks_full_gx, spectra->spectrum_full_gx, num_freq_bins_full, peak_log10_threshold, peak_norm_scale, peak_eps, peak_prominence);
        if (spectra->peaks_full_gy) detect_resonances(spectra->peaks_full_gy, spectra->spectrum_full_gy, num_freq_bins_full, peak_log10_threshold, peak_norm_scale, peak_eps, peak_prominence);
        if (spectra->peaks_full_gz) detect_resonances(spectra->peaks_full_gz, spectra->spectrum_full_gz, num_freq_bins_full, peak_log10_threshold, peak_norm_scale, peak_eps, peak_prominence);
    }

    /* ---- structural acoustic analysis ---- */
    if (desc && block_count > 0) {
        result = sa_check_structural_violations(
            spectra, desc, start_block, block_count,
            num_trs, tr_duration_us,
            num_forbidden_bands, forbidden_bands,
            peak_log10_threshold, peak_norm_scale,
            peak_eps, peak_prominence);
        if (PULSEQLIB_FAILED(result)) {
            pulseqlib_mech_resonances_spectra_free(spectra);
            diag->code = result; return result;
        }
    }

    diag->code = PULSEQLIB_SUCCESS;
    return PULSEQLIB_SUCCESS;
}

/* Select canonical TR window for safety/plotting wrappers.
 * Canonical geometry is always extracted with AMP_MAX_POS.
 * Non-degenerate prep/cooldown: full-pass canonical TR (pass-expanded).
 * Degenerate prep/cooldown: imaging TR canonical window (no pass expansion). */
static void pulseqlib__select_canonical_tr_window(
    const pulseqlib_sequence_descriptor* desc,
    int* start_block,
    int* block_count,
    int* amplitude_mode,
    int* num_instances,
    float* tr_duration_us)
{
    const pulseqlib_tr_descriptor* trd;
    int has_nd_prep, has_nd_cool, n;

    trd = &desc->tr_descriptor;
    has_nd_prep = (trd->num_prep_blocks > 0 && !trd->degenerate_prep);
    has_nd_cool = (trd->num_cooldown_blocks > 0 && !trd->degenerate_cooldown);

    if (has_nd_prep || has_nd_cool) {
        *start_block = 0;
        *block_count = desc->pass_len;
        *amplitude_mode = PULSEQLIB_AMP_MAX_POS;
        *num_instances = (desc->num_passes > 1) ? desc->num_passes : 1;

        *tr_duration_us = 0.0f;
        for (n = 0; n < desc->pass_len; ++n) {
            const pulseqlib_block_table_element* bte = &desc->block_table[n];
            const pulseqlib_block_definition* bdef = &desc->block_definitions[bte->id];
            *tr_duration_us += (bte->duration_us >= 0)
                ? (float)bte->duration_us
                : (float)bdef->duration_us;
        }
        return;
    }

    *start_block = trd->num_prep_blocks + trd->imaging_tr_start;
    *block_count = trd->tr_size;
    *amplitude_mode = PULSEQLIB_AMP_MAX_POS;
    *num_instances = trd->num_trs;
    *tr_duration_us = trd->tr_duration_us;
}

/* Find unique shot-index patterns across pass-expanded canonical windows.
 * Returns number of unique pass patterns, 0 on allocation/shape failure.
 * Caller must free out arrays when return > 0. */
int pulseqlib__find_unique_shot_passes(
    const pulseqlib_sequence_descriptor* desc,
    int** out_unique_pass_indices,
    int** out_pass_group_labels)
{
    int num_passes, pass_size, num_cols;
    int* rows;
    int* unique_defs;
    int* event_table;
    int* result_indices;
    int p, pos, col, bt_pos, block_idx, raw_id;
    int num_unique, g, match, c;
    const pulseqlib_block_table_element* bte;
    const pulseqlib_grad_table_element* gte;

    *out_unique_pass_indices = NULL;
    *out_pass_group_labels = NULL;

    num_passes = (desc->num_passes > 1) ? desc->num_passes : 1;
    pass_size = (num_passes > 0) ? (desc->scan_table_len / num_passes) : 0;
    if (num_passes <= 0 || pass_size <= 0 || !desc->scan_table_block_idx)
        return 0;

    num_cols = pass_size * 3;
    rows = (int*)PULSEQLIB_ALLOC((size_t)num_passes * (size_t)num_cols * sizeof(int));
    unique_defs = (int*)PULSEQLIB_ALLOC((size_t)num_passes * sizeof(int));
    event_table = (int*)PULSEQLIB_ALLOC((size_t)num_passes * sizeof(int));
    if (!rows || !unique_defs || !event_table) {
        if (rows) PULSEQLIB_FREE(rows);
        if (unique_defs) PULSEQLIB_FREE(unique_defs);
        if (event_table) PULSEQLIB_FREE(event_table);
        return 0;
    }

    for (p = 0; p < num_passes; ++p) {
        col = 0;
        for (pos = 0; pos < pass_size; ++pos) {
            bt_pos = p * pass_size + pos;
            if (bt_pos < 0 || bt_pos >= desc->scan_table_len) {
                rows[p * num_cols + col++] = -1;
                rows[p * num_cols + col++] = -1;
                rows[p * num_cols + col++] = -1;
                continue;
            }
            block_idx = desc->scan_table_block_idx[bt_pos];
            if (block_idx < 0 || block_idx >= desc->num_blocks) {
                rows[p * num_cols + col++] = -1;
                rows[p * num_cols + col++] = -1;
                rows[p * num_cols + col++] = -1;
                continue;
            }

            bte = &desc->block_table[block_idx];

            raw_id = bte->gx_id;
            if (raw_id >= 0 && raw_id < desc->grad_table_size) {
                gte = &desc->grad_table[raw_id];
                rows[p * num_cols + col] = gte->shot_index;
            } else rows[p * num_cols + col] = -1;
            col++;

            raw_id = bte->gy_id;
            if (raw_id >= 0 && raw_id < desc->grad_table_size) {
                gte = &desc->grad_table[raw_id];
                rows[p * num_cols + col] = gte->shot_index;
            } else rows[p * num_cols + col] = -1;
            col++;

            raw_id = bte->gz_id;
            if (raw_id >= 0 && raw_id < desc->grad_table_size) {
                gte = &desc->grad_table[raw_id];
                rows[p * num_cols + col] = gte->shot_index;
            } else rows[p * num_cols + col] = -1;
            col++;
        }
    }

    num_unique = 0;
    for (p = 0; p < num_passes; ++p) {
        match = -1;
        for (g = 0; g < num_unique; ++g) {
            int rep = unique_defs[g];
            int equal = 1;
            for (c = 0; c < num_cols; ++c) {
                if (rows[p * num_cols + c] != rows[rep * num_cols + c]) {
                    equal = 0;
                    break;
                }
            }
            if (equal) {
                match = g;
                break;
            }
        }
        if (match >= 0) {
            event_table[p] = match;
        } else {
            unique_defs[num_unique] = p;
            event_table[p] = num_unique;
            num_unique++;
        }
    }

    result_indices = (int*)PULSEQLIB_ALLOC((size_t)num_unique * sizeof(int));
    if (!result_indices) {
        PULSEQLIB_FREE(rows);
        PULSEQLIB_FREE(unique_defs);
        PULSEQLIB_FREE(event_table);
        return 0;
    }
    for (g = 0; g < num_unique; ++g) result_indices[g] = unique_defs[g];

    *out_unique_pass_indices = result_indices;
    *out_pass_group_labels = event_table;

    PULSEQLIB_FREE(rows);
    PULSEQLIB_FREE(unique_defs);
    return num_unique;
}

/* ================================================================== */
/*  Acoustic spectra (public wrapper)                                 */
/* ================================================================== */

int pulseqlib_calc_mech_resonances(
    pulseqlib_mech_resonances_spectra* spectra,
    pulseqlib_diagnostic* diag,
    const pulseqlib_collection* coll,
    int subseq_idx,
    int canonical_tr_idx,
    const pulseqlib_opts* opts,
    int target_window_size,
    float target_resolution_hz,
    float max_freq_hz,
    int num_forbidden_bands,
    const pulseqlib_forbidden_band* forbidden_bands)
{
    const pulseqlib_sequence_descriptor* desc;
    const pulseqlib_tr_descriptor* trd;
    pulseqlib__uniform_grad_waveforms uw;
    pulseqlib_diagnostic local_diag;
    int rc, start_block, block_count, amplitude_mode, num_instances;
    int sa_start_block, sa_block_count;
    int* block_order;
    int has_nd_prep, has_nd_cool;
    float tr_duration_us;
    float peak_log10_threshold;
    float peak_norm_scale;
    float peak_eps;
    float peak_prominence;

    memset(&uw, 0, sizeof(uw));
    block_order = NULL;
    if (!diag) { pulseqlib_diagnostic_init(&local_diag); diag = &local_diag; }
    else       { pulseqlib_diagnostic_init(diag); }
    if (!coll || !spectra) { diag->code = PULSEQLIB_ERR_NULL_POINTER; return diag->code; }
    if (subseq_idx < 0 || subseq_idx >= coll->num_subsequences) {
        diag->code = PULSEQLIB_ERR_INVALID_ARGUMENT; return diag->code;
    }

    peak_log10_threshold = PULSEQLIB_PEAK_LOG10_THRESHOLD_DEFAULT;
    peak_norm_scale = PULSEQLIB_PEAK_NORM_SCALE_DEFAULT;
    peak_eps = PULSEQLIB_PEAK_EPS_DEFAULT;
    peak_prominence = PULSEQLIB_PEAK_PROMINENCE_DEFAULT;
    if (opts) {
        peak_log10_threshold = opts->peak_log10_threshold;
        peak_norm_scale = opts->peak_norm_scale;
        peak_eps = opts->peak_eps;
        peak_prominence = opts->peak_prominence;
    }

    desc = &coll->descriptors[subseq_idx];
    trd = &desc->tr_descriptor;
    /* Select the canonical TR window for the given canonical_tr_idx */
    if (canonical_tr_idx < 0 || canonical_tr_idx >= desc->tr_descriptor.num_trs) {
        diag->code = PULSEQLIB_ERR_INVALID_ARGUMENT; return diag->code;
    }
    /* Use a new helper to select the correct window for the canonical_tr_idx */
    pulseqlib__select_canonical_tr_window_idx(
        desc,
        canonical_tr_idx,
        &start_block,
        &block_count,
        &amplitude_mode,
        &num_instances,
        &tr_duration_us);

    sa_start_block = start_block;
    sa_block_count = block_count;

    has_nd_prep = (trd->num_prep_blocks > 0 && !trd->degenerate_prep);
    has_nd_cool = (trd->num_cooldown_blocks > 0 && !trd->degenerate_cooldown);
    if (has_nd_prep || has_nd_cool) {
        rc = pulseqlib__build_pass_expanded_block_order(
            desc,
            start_block,
            &block_order,
            &block_count,
            &tr_duration_us);
        if (PULSEQLIB_FAILED(rc)) {
            diag->code = rc;
            return rc;
        }
        start_block = 0;
    }

    rc = pulseqlib__get_gradient_waveforms_range(desc, &uw, diag,
        start_block,
        block_count,
        amplitude_mode,
        NULL,
        0,
        block_order);
    if (PULSEQLIB_FAILED(rc)) {
        if (block_order) PULSEQLIB_FREE(block_order);
        return rc;
    }
    rc = calc_mech_resonances_from_uniform(spectra, diag, &uw,
        target_window_size, target_resolution_hz, max_freq_hz,
        num_instances,
        tr_duration_us,
        num_forbidden_bands, forbidden_bands,
        peak_log10_threshold, peak_norm_scale, peak_eps,
        peak_prominence,
        desc, sa_start_block, sa_block_count);
    pulseqlib__uniform_grad_waveforms_free(&uw);
    if (block_order) PULSEQLIB_FREE(block_order);
    return rc;
}

/* ================================================================== */
/*  PNS                                                               */
/* ================================================================== */

/* FIX: outputs before inputs */
static int build_pns_kernel(
    float** kernel, int* kernel_len,
    float dt_us, const pulseqlib_pns_params* params)
{
    int n, i;
    float c_s, dt_s, s_min, tau, denom;
    float* k;

    if (params->chronaxie_us <= 0.0f) return PULSEQLIB_ERR_PNS_INVALID_CHRONAXIE;
    if (params->rheobase_hz_per_m_per_s     <= 0.0f) return PULSEQLIB_ERR_PNS_INVALID_RHEOBASE;
    if (params->alpha        <= 0.0f) return PULSEQLIB_ERR_PNS_INVALID_PARAMS;

    c_s  = params->chronaxie_us * 1e-6f;
    dt_s = dt_us * 1e-6f;
    s_min = params->rheobase_hz_per_m_per_s / params->alpha;

    n = (int)(PNS_KERNEL_DURATION_FACTOR * c_s / dt_s) + 1;
    k = (float*)PULSEQLIB_ALLOC((size_t)n * sizeof(float));
    if (!k) return PULSEQLIB_ERR_ALLOC_FAILED;

    for (i = 0; i < n; ++i) {
        tau   = (float)i * dt_s;
        denom = (c_s + tau) * (c_s + tau);
        k[i]  = (dt_s / s_min) * (c_s / denom);
    }

    *kernel     = k;
    *kernel_len = n;
    return PULSEQLIB_SUCCESS;
}

/* FIX: output before inputs, C89-compliant declarations */
static void compute_slew_rate(
    float* slew_out,
    const float* waveform, int num_samples,
    float dt_us, float gamma_hz_per_tesla)
{
    int i;
    float dt_s;
    float inv_g;

    dt_s = dt_us * 1e-6f;
    inv_g = 1.0f / gamma_hz_per_tesla;

    for (i = 0; i < num_samples - 1; ++i)
        slew_out[i] = ((waveform[i + 1] - waveform[i]) * inv_g) / dt_s;
}

/* FIX: outputs then in-out then scratch then inputs */
static int process_pns_axis_circular(
    float* pns_axis, float* pns_total,
    float* pns_store,
    float* padded_waveform, float* slew_rate, float* pns_conv,
    const float* waveform, int num_samples,
    const float* kernel, int kernel_len,
    float grad_raster_us, float gamma_hz_per_tesla,
    int full_output_len)
{
    int i, padded_len, slew_len, rc;

    (void)full_output_len;

    if (num_samples <= 0 || !waveform) return PULSEQLIB_SUCCESS;

    padded_len = num_samples + kernel_len;
    slew_len   = padded_len - 1;

    for (i = 0; i < num_samples; ++i)
        padded_waveform[i] = waveform[i];
    for (i = 0; i < kernel_len; ++i)
        padded_waveform[num_samples + i] = waveform[i % num_samples];

    compute_slew_rate(slew_rate, padded_waveform, padded_len, grad_raster_us, gamma_hz_per_tesla);

    rc = pulseqlib__calc_convolution_fft(pns_conv, slew_rate, slew_len, kernel, kernel_len);
    if (PULSEQLIB_FAILED(rc)) return rc;

    for (i = 0; i < slew_len; ++i) {
        pns_axis[i]   = pns_conv[i] * 100.0f;
        pns_total[i] += pns_conv[i] * pns_conv[i];
    }

    if (pns_store) {
        for (i = 0; i < slew_len; ++i) pns_store[i] = pns_axis[i];
    }
    return PULSEQLIB_SUCCESS;
}

static int calc_pns_from_uniform(
    pulseqlib_pns_result* result,
    pulseqlib_diagnostic* diag,
    float gamma_hz_per_tesla,
    const pulseqlib__uniform_grad_waveforms* waveforms,
    const pulseqlib_pns_params* params)
{
    pulseqlib_diagnostic local_diag;
    int max_samples, padded_len, slew_len, full_output_len;
    int kernel_len, i;
    float* kernel;
    float* padded;
    float* slew;
    float* conv;
    float* axis;
    float* pns_x;
    float* pns_y;
    float* pns_z;
    float* pns_tot;
    int rc;

    kernel  = NULL;
    padded  = NULL;
    slew    = NULL;
    conv    = NULL;
    axis    = NULL;
    pns_x   = NULL;
    pns_y   = NULL;
    pns_z   = NULL;
    pns_tot = NULL;
    rc      = PULSEQLIB_SUCCESS;

    if (!diag) { pulseqlib_diagnostic_init(&local_diag); diag = &local_diag; }
    else       { pulseqlib_diagnostic_init(diag); }

    if (!waveforms || !params || !result) {
        diag->code = PULSEQLIB_ERR_NULL_POINTER; return diag->code;
    }

    memset(result, 0, sizeof(*result));

    max_samples = waveforms->num_samples;
    if (max_samples <= 1) { diag->code = PULSEQLIB_ERR_PNS_NO_WAVEFORM; return diag->code; }

    rc = build_pns_kernel(&kernel, &kernel_len, waveforms->raster_us, params);
    if (PULSEQLIB_FAILED(rc)) { diag->code = rc; return rc; }

    padded_len      = max_samples + kernel_len;
    slew_len        = padded_len - 1;
    full_output_len = slew_len;

    padded  = (float*)PULSEQLIB_ALLOC((size_t)padded_len * sizeof(float));
    slew    = (float*)PULSEQLIB_ALLOC((size_t)slew_len * sizeof(float));
    conv    = (float*)PULSEQLIB_ALLOC((size_t)slew_len * sizeof(float));
    axis    = (float*)PULSEQLIB_ALLOC((size_t)full_output_len * sizeof(float));
    pns_tot = (float*)PULSEQLIB_ALLOC((size_t)full_output_len * sizeof(float));
    if (!padded || !slew || !conv || !axis || !pns_tot) {
        rc = PULSEQLIB_ERR_ALLOC_FAILED; goto fail;
    }
    for (i = 0; i < full_output_len; ++i) pns_tot[i] = 0.0f;

    pns_x = (float*)PULSEQLIB_ALLOC((size_t)full_output_len * sizeof(float));
    pns_y = (float*)PULSEQLIB_ALLOC((size_t)full_output_len * sizeof(float));
    pns_z = (float*)PULSEQLIB_ALLOC((size_t)full_output_len * sizeof(float));
    if (!pns_x || !pns_y || !pns_z) { rc = PULSEQLIB_ERR_ALLOC_FAILED; goto fail; }
    for (i = 0; i < full_output_len; ++i) { pns_x[i] = 0.0f; pns_y[i] = 0.0f; pns_z[i] = 0.0f; }

    /* X */
    rc = process_pns_axis_circular(axis, pns_tot,
        pns_x,
        padded, slew, conv,
        waveforms->gx, waveforms->num_samples,
        kernel, kernel_len,
        waveforms->raster_us, gamma_hz_per_tesla,
        full_output_len);
    if (PULSEQLIB_FAILED(rc)) goto fail;

    /* Y */
    rc = process_pns_axis_circular(axis, pns_tot,
        pns_y,
        padded, slew, conv,
        waveforms->gy, waveforms->num_samples,
        kernel, kernel_len,
        waveforms->raster_us, gamma_hz_per_tesla,
        full_output_len);
    if (PULSEQLIB_FAILED(rc)) goto fail;

    /* Z */
    rc = process_pns_axis_circular(axis, pns_tot,
        pns_z,
        padded, slew, conv,
        waveforms->gz, waveforms->num_samples,
        kernel, kernel_len,
        waveforms->raster_us, gamma_hz_per_tesla,
        full_output_len);
    if (PULSEQLIB_FAILED(rc)) goto fail;

    for (i = 0; i < full_output_len; ++i) {
        pns_tot[i] = 100.0f * (float)sqrt((double)pns_tot[i]);
    }

    result->num_samples    = full_output_len;
    result->slew_x_hz_per_m_per_s = pns_x; pns_x = NULL;
    result->slew_y_hz_per_m_per_s = pns_y; pns_y = NULL;
    result->slew_z_hz_per_m_per_s = pns_z; pns_z = NULL;

    rc = PULSEQLIB_SUCCESS;
    diag->code = rc;

fail:
    if (kernel)  PULSEQLIB_FREE(kernel);
    if (padded)  PULSEQLIB_FREE(padded);
    if (slew)    PULSEQLIB_FREE(slew);
    if (conv)    PULSEQLIB_FREE(conv);
    if (axis)    PULSEQLIB_FREE(axis);
    if (pns_tot) PULSEQLIB_FREE(pns_tot);
    if (pns_x)   PULSEQLIB_FREE(pns_x);
    if (pns_y)   PULSEQLIB_FREE(pns_y);
    if (pns_z)   PULSEQLIB_FREE(pns_z);
    return rc;
}

/* ================================================================== */
/*  PNS (public wrapper)                                              */
/* ================================================================== */

int pulseqlib_calc_pns(
    pulseqlib_pns_result* result,
    pulseqlib_diagnostic* diag,
    const pulseqlib_collection* coll,
    int subseq_idx,
    int canonical_tr_idx,
    const pulseqlib_opts* opts,
    const pulseqlib_pns_params* params)
{
    const pulseqlib_sequence_descriptor* desc;
    const pulseqlib_tr_descriptor* trd;
    pulseqlib__uniform_grad_waveforms uw;
    pulseqlib_diagnostic local_diag;
    int rc, start_block, block_count, amplitude_mode, num_instances;
    int* block_order;
    int has_nd_prep, has_nd_cool;
    float tr_duration_us;

    (void)opts;
    memset(&uw, 0, sizeof(uw));
    block_order = NULL;
    if (!diag) { pulseqlib_diagnostic_init(&local_diag); diag = &local_diag; }
    else       { pulseqlib_diagnostic_init(diag); }
    if (!coll || !result || !params) { diag->code = PULSEQLIB_ERR_NULL_POINTER; return diag->code; }
    if (subseq_idx < 0 || subseq_idx >= coll->num_subsequences) {
        diag->code = PULSEQLIB_ERR_INVALID_ARGUMENT; return diag->code;
    }
    desc = &coll->descriptors[subseq_idx];
    trd = &desc->tr_descriptor;
    if (canonical_tr_idx < 0 || canonical_tr_idx >= desc->tr_descriptor.num_trs) {
        diag->code = PULSEQLIB_ERR_INVALID_ARGUMENT; return diag->code;
    }
    pulseqlib__select_canonical_tr_window_idx(
        desc,
        canonical_tr_idx,
        &start_block,
        &block_count,
        &amplitude_mode,
        &num_instances,
        &tr_duration_us);
    (void)num_instances;
    (void)tr_duration_us;

    has_nd_prep = (trd->num_prep_blocks > 0 && !trd->degenerate_prep);
    has_nd_cool = (trd->num_cooldown_blocks > 0 && !trd->degenerate_cooldown);
    if (has_nd_prep || has_nd_cool) {
        rc = pulseqlib__build_pass_expanded_block_order(
            desc,
            start_block,
            &block_order,
            &block_count,
            NULL);
        if (PULSEQLIB_FAILED(rc)) {
            diag->code = rc;
            return rc;
        }
        start_block = 0;
    }

    rc = pulseqlib__get_gradient_waveforms_range(desc, &uw, diag,
        start_block,
        block_count,
        amplitude_mode,
        NULL,
        0,
        block_order);
    if (PULSEQLIB_FAILED(rc)) {
        if (block_order) PULSEQLIB_FREE(block_order);
        return rc;
    }
    rc = calc_pns_from_uniform(result, diag, opts->gamma_hz_per_t, &uw, params);
    pulseqlib__uniform_grad_waveforms_free(&uw);
    if (block_order) PULSEQLIB_FREE(block_order);
    return rc;
}

/* ================================================================== */
/*  PNS result free                                                   */
/* ================================================================== */

void pulseqlib_pns_result_free(pulseqlib_pns_result* r)
{
    if (!r) return;
    if (r->slew_x_hz_per_m_per_s)     { PULSEQLIB_FREE(r->slew_x_hz_per_m_per_s);     r->slew_x_hz_per_m_per_s = NULL; }
    if (r->slew_y_hz_per_m_per_s)     { PULSEQLIB_FREE(r->slew_y_hz_per_m_per_s);     r->slew_y_hz_per_m_per_s = NULL; }
    if (r->slew_z_hz_per_m_per_s)     { PULSEQLIB_FREE(r->slew_z_hz_per_m_per_s);     r->slew_z_hz_per_m_per_s = NULL; }
    r->num_samples = 0;
}

/* ================================================================== */
/*  Collection-level safety check                                     */
/* ================================================================== */

int check_max_grad(
    const pulseqlib_collection* coll,
    pulseqlib_diagnostic* diag,
    const pulseqlib_opts* opts
) {
    int s, b, raw_id;
    int worst_subseq, worst_block;
    float gx_amp, gy_amp, gz_amp, gsos, gsos_max, limit_sq, hz_per_mt;
    const pulseqlib_sequence_descriptor* desc;
    const pulseqlib_block_table_element* bte;

    if (!coll || !opts) {
        if (diag) { pulseqlib_diagnostic_init(diag); diag->code = PULSEQLIB_ERR_NULL_POINTER; }
        return PULSEQLIB_ERR_NULL_POINTER;
    }
    if (diag) pulseqlib_diagnostic_init(diag);

    /* ---- max gradient amplitude (GSOS) check ---- */
    gsos_max     = 0.0f;
    limit_sq     = opts->max_grad_hz_per_m * opts->max_grad_hz_per_m;
    worst_subseq = 0;
    worst_block  = 0;

    for (s = 0; s < coll->num_subsequences; ++s) {
        desc = &coll->descriptors[s];
        for (b = 0; b < desc->num_blocks; ++b) {
            bte = &desc->block_table[b];

            gx_amp = 0.0f;
            raw_id = bte->gx_id;
            if (raw_id >= 0 && raw_id < desc->grad_table_size)
                gx_amp = desc->grad_table[raw_id].amplitude;

            gy_amp = 0.0f;
            raw_id = bte->gy_id;
            if (raw_id >= 0 && raw_id < desc->grad_table_size)
                gy_amp = desc->grad_table[raw_id].amplitude;

            gz_amp = 0.0f;
            raw_id = bte->gz_id;
            if (raw_id >= 0 && raw_id < desc->grad_table_size)
                gz_amp = desc->grad_table[raw_id].amplitude;

            gsos = gx_amp * gx_amp + gy_amp * gy_amp + gz_amp * gz_amp;
            if (gsos > gsos_max) {
                gsos_max     = gsos;
                worst_subseq = s;
                worst_block  = b;
            }
        }
    }

    if (limit_sq > 0.0f && gsos_max > limit_sq) {
        hz_per_mt = opts->gamma_hz_per_t * 0.001f;
        if (diag) {
            diag->code = PULSEQLIB_ERR_MAX_GRAD_EXCEEDED;
            pulseqlib__diag_printf(diag,
                "max grad exceeded: subseq=%d block=%d amp=%.4f limit=%.4f mT/m",
                worst_subseq, worst_block,
                (double)((float)sqrt((double)gsos_max) / hz_per_mt),
                (double)(opts->max_grad_hz_per_m / hz_per_mt));
        }
        return PULSEQLIB_ERR_MAX_GRAD_EXCEEDED;
    }

    return PULSEQLIB_SUCCESS;
}

/* ================================================================== */
/*  Gradient continuity check (cursor dry-run over n repetitions)     */
/* ================================================================== */

int check_grad_continuity(
    pulseqlib_collection* coll,
    pulseqlib_diagnostic* diag,
    const pulseqlib_opts* opts)
{
    pulseqlib_block_cursor saved_cursor;
    const pulseqlib_sequence_descriptor* desc;
    const pulseqlib_block_table_element* bte;
    const pulseqlib_block_definition* bdef;
    const pulseqlib_grad_definition* gdef;
    const pulseqlib_grad_table_element* gte;
    int n, raw_id, rot_id, status, cur_seq;
    int grad_def_ids[3];
    int shot_idx[3];
    float amp[3], first_val[3], last_val[3];
    float first_phys[3], last_phys[3], prev_phys[3];
    float max_allowed, grad_raster_s, step, hz_per_mt;

    if (!coll || !opts) {
        if (diag) { pulseqlib_diagnostic_init(diag); diag->code = PULSEQLIB_ERR_NULL_POINTER; }
        return PULSEQLIB_ERR_NULL_POINTER;
    }
    if (diag) pulseqlib_diagnostic_init(diag);

    /* save cursor state */
    saved_cursor = coll->block_cursor;
    coll->block_cursor.sequence_index = 0;
    coll->block_cursor.scan_table_position = 0;
    coll->block_cursor.from_last_reset = 0;

    prev_phys[0] = 0.0f;
    prev_phys[1] = 0.0f;
    prev_phys[2] = 0.0f;
    cur_seq = 0;
    status  = PULSEQLIB_CURSOR_BLOCK;

    desc = &coll->descriptors[0];
    grad_raster_s = desc->grad_raster_us * 1e-6f;
    max_allowed   = opts->max_slew_hz_per_m_per_s * grad_raster_s;

    while (status != PULSEQLIB_CURSOR_DONE) {
        /* detect subsequence change */
        if (coll->block_cursor.sequence_index != cur_seq) {
            /* end-of-subsequence: prev must ramp to zero */
            for (n = 0; n < 3; ++n) {
                step = prev_phys[n];
                if (step < 0.0f) step = -step;
                if (step > max_allowed) {
                    hz_per_mt = opts->gamma_hz_per_t * 0.001f;
                    if (diag) {
                        diag->code = PULSEQLIB_ERR_GRAD_DISCONTINUITY;
                        pulseqlib__diag_printf(diag,
                            "grad discontinuity at subseq boundary: axis=%d step=%.4f limit=%.4f mT/m",
                            n, (double)(step / hz_per_mt), (double)(max_allowed / hz_per_mt));
                    }
                    coll->block_cursor = saved_cursor;
                    return PULSEQLIB_ERR_GRAD_DISCONTINUITY;
                }
            }

            cur_seq = coll->block_cursor.sequence_index;
            prev_phys[0] = 0.0f;
            prev_phys[1] = 0.0f;
            prev_phys[2] = 0.0f;

            desc = &coll->descriptors[cur_seq];
            grad_raster_s = desc->grad_raster_us * 1e-6f;
            max_allowed   = opts->max_slew_hz_per_m_per_s * grad_raster_s;
        }

        /* read current block */
        {
            int bt_idx = desc->scan_table_block_idx[coll->block_cursor.scan_table_position];
            bte  = &desc->block_table[bt_idx];
        }
        bdef = &desc->block_definitions[bte->id];

        /* grad table: amplitude + shot_index */
        grad_def_ids[0] = bdef->gx_id;
        grad_def_ids[1] = bdef->gy_id;
        grad_def_ids[2] = bdef->gz_id;

        raw_id = bte->gx_id;
        if (raw_id >= 0 && raw_id < desc->grad_table_size) {
            gte = &desc->grad_table[raw_id];
            amp[0] = gte->amplitude; shot_idx[0] = gte->shot_index;
        } else { amp[0] = 0.0f; shot_idx[0] = 0; }

        raw_id = bte->gy_id;
        if (raw_id >= 0 && raw_id < desc->grad_table_size) {
            gte = &desc->grad_table[raw_id];
            amp[1] = gte->amplitude; shot_idx[1] = gte->shot_index;
        } else { amp[1] = 0.0f; shot_idx[1] = 0; }

        raw_id = bte->gz_id;
        if (raw_id >= 0 && raw_id < desc->grad_table_size) {
            gte = &desc->grad_table[raw_id];
            amp[2] = gte->amplitude; shot_idx[2] = gte->shot_index;
        } else { amp[2] = 0.0f; shot_idx[2] = 0; }

        /* first_value / last_value from grad definitions, scaled by amplitude */
        for (n = 0; n < 3; ++n) {
            if (grad_def_ids[n] >= 0 && grad_def_ids[n] < desc->num_unique_grads) {
                gdef = &desc->grad_definitions[grad_def_ids[n]];
                first_val[n] = gdef->first_value[shot_idx[n]] * amp[n];
                last_val[n]  = gdef->last_value[shot_idx[n]]  * amp[n];
            } else {
                first_val[n] = 0.0f;
                last_val[n]  = 0.0f;
            }
        }

        /* transform logical -> physical */
        rot_id = bte->rotation_id;
        if (rot_id >= 0 && rot_id < desc->num_rotations) {
            pulseqlib__apply_rotation(first_phys, desc->rotation_matrices[rot_id], first_val, 1);
            pulseqlib__apply_rotation(last_phys,  desc->rotation_matrices[rot_id], last_val,  1);
        } else {
            first_phys[0] = first_val[0]; first_phys[1] = first_val[1]; first_phys[2] = first_val[2];
            last_phys[0]  = last_val[0];  last_phys[1]  = last_val[1];  last_phys[2]  = last_val[2];
        }

        /* continuity check */
        for (n = 0; n < 3; ++n) {
            step = first_phys[n] - prev_phys[n];
            if (step < 0.0f) step = -step;
            if (step > max_allowed) {
                hz_per_mt = opts->gamma_hz_per_t * 0.001f;
                if (diag) {
                    diag->code = PULSEQLIB_ERR_GRAD_DISCONTINUITY;
                    pulseqlib__diag_printf(diag,
                        "grad discontinuity: axis=%d scan_pos=%d step=%.4f limit=%.4f mT/m",
                        n, coll->block_cursor.scan_table_position,
                        (double)(step / hz_per_mt), (double)(max_allowed / hz_per_mt));
                }
                coll->block_cursor = saved_cursor;
                return PULSEQLIB_ERR_GRAD_DISCONTINUITY;
            }
        }

        prev_phys[0] = last_phys[0];
        prev_phys[1] = last_phys[1];
        prev_phys[2] = last_phys[2];

        /* advance cursor */
        status = pulseqlib_cursor_next(coll);
    }

    /* final subsequence trailing edge */
    for (n = 0; n < 3; ++n) {
        step = prev_phys[n];
        if (step < 0.0f) step = -step;
        if (step > max_allowed) {
            hz_per_mt = opts->gamma_hz_per_t * 0.001f;
            if (diag) {
                diag->code = PULSEQLIB_ERR_GRAD_DISCONTINUITY;
                pulseqlib__diag_printf(diag,
                    "grad discontinuity at trailing edge: axis=%d step=%.4f limit=%.4f mT/m",
                    n, (double)(step / hz_per_mt), (double)(max_allowed / hz_per_mt));
            }
            coll->block_cursor = saved_cursor;
            return PULSEQLIB_ERR_GRAD_DISCONTINUITY;
        }
    }

    coll->block_cursor = saved_cursor;
    return PULSEQLIB_SUCCESS;
}

/* ================================================================== */
/*  Max slew rate check (per unique block definition)                 */
/* ================================================================== */

int check_max_slew(
    const pulseqlib_collection* coll,
    pulseqlib_diagnostic* diag,
    const pulseqlib_opts* opts)
{
    int s, d, n, grad_id, shot;
    float slew_limit, slew_phys;
    const pulseqlib_sequence_descriptor* desc;
    const pulseqlib_block_definition* bdef;
    const pulseqlib_grad_definition* gdef;
    int grad_ids[3];

    if (!coll || !opts) {
        if (diag) { pulseqlib_diagnostic_init(diag); diag->code = PULSEQLIB_ERR_NULL_POINTER; }
        return PULSEQLIB_ERR_NULL_POINTER;
    }
    if (diag) pulseqlib_diagnostic_init(diag);

    slew_limit = opts->max_slew_hz_per_m_per_s / (float)sqrt(3.0);

    for (s = 0; s < coll->num_subsequences; ++s) {
        desc = &coll->descriptors[s];

        for (d = 0; d < desc->num_unique_blocks; ++d) {
            bdef = &desc->block_definitions[d];
            grad_ids[0] = bdef->gx_id;
            grad_ids[1] = bdef->gy_id;
            grad_ids[2] = bdef->gz_id;

            for (n = 0; n < 3; ++n) {
                grad_id = grad_ids[n];
                if (grad_id < 0 || grad_id >= desc->num_unique_grads)
                    continue;

                gdef = &desc->grad_definitions[grad_id];
                for (shot = 0; shot < gdef->num_shots; ++shot) {
                    slew_phys = gdef->slew_rate[shot] * gdef->max_amplitude[shot];
                    if (slew_phys > slew_limit) {
                        if (diag) {
                            diag->code = PULSEQLIB_ERR_MAX_SLEW_EXCEEDED;
                            pulseqlib__diag_printf(diag,
                                "max slew exceeded: axis=%d def=%d slew=%.4f limit=%.4f Hz/m/s",
                                n, d, (double)slew_phys, (double)slew_limit);
                        }
                        return PULSEQLIB_ERR_MAX_SLEW_EXCEEDED;
                    }
                }
            }
        }
    }

    return PULSEQLIB_SUCCESS;
}

/* ================================================================== */
/*  Safety check                                                      */
/* ================================================================== */
int pulseqlib_check_safety(
    pulseqlib_collection* coll,
    pulseqlib_diagnostic* diag,
    const pulseqlib_opts* opts,
    int num_forbidden_bands,
    const pulseqlib_forbidden_band* forbidden_bands,
    const pulseqlib_pns_params* pns_params,
    float pns_threshold_percent)
{
    int rc, s, u, i;
    int ci, b;
    int num_unique_trs;
    int* unique_tr_indices;
    int* tr_group_labels;
    const pulseqlib_sequence_descriptor* desc;
    const pulseqlib_tr_descriptor* trd;
    pulseqlib__uniform_grad_waveforms uw;
    pulseqlib_mech_resonances_spectra spectra;
    pulseqlib_pns_result pns_result;
    int start_block, block_count, amplitude_mode, num_instances;
    int sa_start_block, sa_block_count;
    int* block_order;
    int has_nd_prep, has_nd_cool;
    float tr_duration_us;
    float pns_combined, max_pns;
    float cf_hz, ca_hz_per_m;

    if (!coll || !opts) {
        if (diag) { pulseqlib_diagnostic_init(diag); diag->code = PULSEQLIB_ERR_NULL_POINTER; }
        return PULSEQLIB_ERR_NULL_POINTER;
    }
    if (diag) pulseqlib_diagnostic_init(diag);

    /* ---- 1. max gradient amplitude ---- */
    rc = check_max_grad(coll, diag, opts);
    if (PULSEQLIB_FAILED(rc)) return rc;

    /* ---- 2. gradient continuity ---- */
    rc = check_grad_continuity(coll, diag, opts);
    if (PULSEQLIB_FAILED(rc)) return rc;

    /* ---- 3. max slew rate ---- */
    rc = check_max_slew(coll, diag, opts);
    if (PULSEQLIB_FAILED(rc)) return rc;

    /* ---- 4. per-subsequence canonical-TR acoustic + PNS ---- */
    for (s = 0; s < coll->num_subsequences; ++s) {
        desc = &coll->descriptors[s];
        trd = &desc->tr_descriptor;
        unique_tr_indices = NULL;
        tr_group_labels = NULL;

        pulseqlib__select_canonical_tr_window(
            desc,
            &start_block,
            &block_count,
            &amplitude_mode,
            &num_instances,
            &tr_duration_us);
        (void)amplitude_mode;
        sa_start_block = start_block;
        sa_block_count = block_count;
        block_order = NULL;
        has_nd_prep = (trd->num_prep_blocks > 0 && !trd->degenerate_prep);
        has_nd_cool = (trd->num_cooldown_blocks > 0 && !trd->degenerate_cooldown);
        if (has_nd_prep || has_nd_cool) {
            rc = pulseqlib__build_pass_expanded_block_order(
                desc,
                start_block,
                &block_order,
                &block_count,
                &tr_duration_us);
            if (PULSEQLIB_FAILED(rc)) {
                if (unique_tr_indices) PULSEQLIB_FREE(unique_tr_indices);
                if (tr_group_labels)   PULSEQLIB_FREE(tr_group_labels);
                if (block_order) PULSEQLIB_FREE(block_order);
                return rc;
            }
            start_block = 0;
        }

        /* Evaluate one canonical TR per shot-ID combination. */
        unique_tr_indices = NULL;
        tr_group_labels   = NULL;
        if ((trd->num_prep_blocks > 0 && !trd->degenerate_prep) ||
            (trd->num_cooldown_blocks > 0 && !trd->degenerate_cooldown)) {
            num_unique_trs = pulseqlib__find_unique_shot_passes(
                desc, &unique_tr_indices, &tr_group_labels);
        } else {
            num_unique_trs = pulseqlib__find_unique_shot_trs(
                desc, &unique_tr_indices, &tr_group_labels);
        }
        if (num_unique_trs <= 0) num_unique_trs = 1;

        for (u = 0; u < num_unique_trs; ++u) {
            memset(&uw, 0, sizeof(uw));
            rc = pulseqlib__get_gradient_waveforms_range(desc, &uw, diag,
                start_block,
                block_count,
                PULSEQLIB_AMP_MAX_POS,
                tr_group_labels, u, block_order);
            if (PULSEQLIB_FAILED(rc)) {
                if (unique_tr_indices) PULSEQLIB_FREE(unique_tr_indices);
                if (tr_group_labels)   PULSEQLIB_FREE(tr_group_labels);
                if (block_order) PULSEQLIB_FREE(block_order);
                return rc;
            }

            if (num_forbidden_bands > 0) {
                memset(&spectra, 0, sizeof(spectra));
                rc = calc_mech_resonances_from_uniform(
                    &spectra, diag, &uw,
                    0, 0.0f, 0.0f,
                    num_instances,
                    tr_duration_us,
                    num_forbidden_bands, forbidden_bands,
                    opts->peak_log10_threshold,
                    opts->peak_norm_scale,
                    opts->peak_eps,
                    opts->peak_prominence,
                    desc, sa_start_block, sa_block_count);

                /* Safety path: fail fast on first violating candidate.
                 * Pattern: for each candidate, scan union of all bands. */
                if (!PULSEQLIB_FAILED(rc) && spectra.num_candidates > 0 &&
                    spectra.candidate_freqs &&
                    spectra.candidate_grad_amps) {
                    for (ci = 0; ci < spectra.num_candidates; ++ci) {
                        cf_hz = spectra.candidate_freqs[ci];
                        ca_hz_per_m = spectra.candidate_grad_amps[ci];

                        for (b = 0; b < num_forbidden_bands; ++b) {
                            if (cf_hz < forbidden_bands[b].freq_min_hz ||
                                cf_hz > forbidden_bands[b].freq_max_hz)
                                continue;
                            if (ca_hz_per_m > forbidden_bands[b].max_amplitude_hz_per_m) {
                                rc = PULSEQLIB_ERR_MECH_RESONANCES_VIOLATION;
                                if (diag) {
                                    diag->code = rc;
                                    pulseqlib__diag_printf(
                                        diag,
                                        "mechanical resonance violation: ss=%d ctr=%d f=%.3fHz amp=%.3fHz/m band=[%.3f,%.3f]Hz limit=%.3fHz/m",
                                        s, u,
                                        (double)cf_hz, (double)ca_hz_per_m,
                                        (double)forbidden_bands[b].freq_min_hz,
                                        (double)forbidden_bands[b].freq_max_hz,
                                        (double)forbidden_bands[b].max_amplitude_hz_per_m);
                                }
                                break;
                            }
                        }
                        if (PULSEQLIB_FAILED(rc)) break;
                    }
                }

                pulseqlib_mech_resonances_spectra_free(&spectra);
                if (PULSEQLIB_FAILED(rc)) {
                    pulseqlib__uniform_grad_waveforms_free(&uw);
                    if (unique_tr_indices) PULSEQLIB_FREE(unique_tr_indices);
                    if (tr_group_labels)   PULSEQLIB_FREE(tr_group_labels);
                    if (block_order) PULSEQLIB_FREE(block_order);
                    return rc;
                }
            }

            if (pns_params) {
                memset(&pns_result, 0, sizeof(pns_result));
                rc = calc_pns_from_uniform(
                    &pns_result, diag, opts->gamma_hz_per_t,
                    &uw, pns_params);
                if (!PULSEQLIB_FAILED(rc) && pns_result.num_samples > 0) {
                    max_pns = 0.0f;
                    for (i = 0; i < pns_result.num_samples; ++i) {
                        pns_combined = 0.0f;
                        if (pns_result.slew_x_hz_per_m_per_s)
                            pns_combined += pns_result.slew_x_hz_per_m_per_s[i] * pns_result.slew_x_hz_per_m_per_s[i];
                        if (pns_result.slew_y_hz_per_m_per_s)
                            pns_combined += pns_result.slew_y_hz_per_m_per_s[i] * pns_result.slew_y_hz_per_m_per_s[i];
                        if (pns_result.slew_z_hz_per_m_per_s)
                            pns_combined += pns_result.slew_z_hz_per_m_per_s[i] * pns_result.slew_z_hz_per_m_per_s[i];
                        pns_combined = (float)sqrt((double)pns_combined);
                        if (pns_combined > max_pns) max_pns = pns_combined;
                    }
                    if (max_pns > pns_threshold_percent) rc = PULSEQLIB_ERR_PNS_THRESHOLD_EXCEEDED;
                }
                pulseqlib_pns_result_free(&pns_result);
                if (PULSEQLIB_FAILED(rc)) {
                    pulseqlib__uniform_grad_waveforms_free(&uw);
                    if (unique_tr_indices) PULSEQLIB_FREE(unique_tr_indices);
                    if (tr_group_labels)   PULSEQLIB_FREE(tr_group_labels);
                    if (block_order) PULSEQLIB_FREE(block_order);
                    return rc;
                }
            }

            pulseqlib__uniform_grad_waveforms_free(&uw);
        }

        if (unique_tr_indices) PULSEQLIB_FREE(unique_tr_indices);
        if (tr_group_labels)   PULSEQLIB_FREE(tr_group_labels);
        if (block_order) PULSEQLIB_FREE(block_order);
    }

    return PULSEQLIB_SUCCESS;
}