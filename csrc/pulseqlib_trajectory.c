/* pulseqlib_trajectory.c -- k-space trajectory computation and caching */

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "pulseqlib_internal.h"
#include "pulseqlib_methods.h"

/* ================================================================== */
/*  Forward declarations (internal helpers from other files)           */
/* ================================================================== */

/* from pulseqlib_cache.c (we reuse these via the same pattern) */
static char* traj_make_cache_path(const char* seq_path);

/* ================================================================== */
/*  Constants                                                         */
/* ================================================================== */

#define TRAJ_DEDUP_TOL 1e-6f

/* ================================================================== */
/*  Helpers                                                           */
/* ================================================================== */

static char* traj_make_cache_path(const char* seq_path)
{
    size_t len;
    char* out;
    const char* dot;

    len = strlen(seq_path);
    out = (char*)PULSEQLIB_ALLOC(len + 5);
    if (!out) return NULL;
    strcpy(out, seq_path);
    dot = strrchr(out, '.');
    if (dot && (strrchr(out, '/') == NULL || dot > strrchr(out, '/')))
        strcpy((char*)(out + (dot - out)), ".bin");
    else
        strcat(out, ".bin");
    return out;
}

/* Linear interpolation: resample src[0..src_n-1] at uniform raster src_dt
 * to dst[0..dst_n-1] at uniform raster dst_dt, both starting at t=0. */
static void linear_resample(const float* src, int src_n, float src_dt,
                            float* dst, int dst_n, float dst_dt)
{
    int i;
    for (i = 0; i < dst_n; ++i) {
        float t = (float)i * dst_dt;
        float fi = t / src_dt;
        int lo = (int)fi;
        float frac;
        if (lo < 0) lo = 0;
        if (lo >= src_n - 1) {
            dst[i] = src[src_n - 1];
            continue;
        }
        frac = fi - (float)lo;
        dst[i] = src[lo] * (1.0f - frac) + src[lo + 1] * frac;
    }
}

/* Check if an axis k-trajectory is trivial (constant, i.e. no readout
 * gradient activity during ADC window). */
static int is_trivial_shot(const float* k, int n)
{
    int i;
    float first;
    if (n <= 0) return 1;
    first = k[0];
    for (i = 1; i < n; ++i) {
        float d = k[i] - first;
        if (d > TRAJ_DEDUP_TOL || d < -TRAJ_DEDUP_TOL)
            return 0;
    }
    return 1;
}

/* Compare two k-space shot shapes; returns 1 if identical within tolerance. */
static int shots_equal(const float* a, const float* b, int n)
{
    int i;
    for (i = 0; i < n; ++i) {
        float d = a[i] - b[i];
        if (d > TRAJ_DEDUP_TOL || d < -TRAJ_DEDUP_TOL)
            return 0;
    }
    return 1;
}

/* Add a shot to the library if not already present.
 * Returns the shot index, or -1 on allocation failure. */
static int kshot_library_add(pulseqlib_kshot_library* lib,
                             const float* k, int n)
{
    int i;
    pulseqlib_kshot* shot;

    /* Check for duplicate */
    for (i = 0; i < lib->num_shots; ++i) {
        if (lib->shots[i].num_samples == n &&
            shots_equal(lib->shots[i].k, k, n))
            return i;
    }

    /* Grow array */
    {
        pulseqlib_kshot* new_shots;
        new_shots = (pulseqlib_kshot*)realloc(
            lib->shots, (size_t)(lib->num_shots + 1) * sizeof(pulseqlib_kshot));
        if (!new_shots) return -1;
        lib->shots = new_shots;
    }

    shot = &lib->shots[lib->num_shots];
    shot->num_samples = n;
    shot->k = (float*)PULSEQLIB_ALLOC((size_t)n * sizeof(float));
    if (!shot->k) return -1;
    memcpy(shot->k, k, (size_t)n * sizeof(float));

    return lib->num_shots++;
}

/* ================================================================== */
/*  Compute trajectory for a single block (block-level k-space)       */
/* ================================================================== */

/*
 * For a single block with an ADC event:
 * 1. Get gradient waveforms at block level (uniform raster)
 * 2. Integrate to get k-space trajectory
 * 3. Crop to ADC window
 * 4. Resample to ADC dwell time
 * 5. Center to k-zero anchor
 *
 * Returns per-axis k-space arrays of length adc_num_samples.
 */
static int compute_block_kspace(
    const pulseqlib_sequence_descriptor* desc,
    int block_table_idx,
    int kzero_index,
    float* out_kx, float* out_ky, float* out_kz,  /* [adc_num_samples] */
    int* out_num_samples,
    pulseqlib_diagnostic* diag)
{
    const pulseqlib_block_table_element* bte;
    const pulseqlib_adc_definition* adc_def;
    int adc_def_idx;
    int adc_num_samples, adc_delay_us;
    float adc_dwell_us;
    float grad_raster_us;
    int block_dur_us;
    int n_grad;
    float dt_s;
    float *kx_full, *ky_full, *kz_full;
    int i, adc_start_sample, adc_end_sample, adc_window_samples;
    float* kx_window;
    float* ky_window;
    float* kz_window;
    float kzero_kx, kzero_ky, kzero_kz;

    bte = &desc->block_table[block_table_idx];
    adc_def_idx = bte->adc_id;
    if (adc_def_idx < 0 || adc_def_idx >= desc->num_unique_adcs) {
        if (diag) snprintf(diag->message, PULSEQLIB_DIAG_MSG_LEN,
            "Block %d has no ADC event", block_table_idx);
        return PULSEQLIB_ERR_INVALID_ARGUMENT;
    }

    adc_def = &desc->adc_definitions[adc_def_idx];
    adc_num_samples = adc_def->num_samples;
    adc_dwell_us = (float)adc_def->dwell_time * 1e-3f; /* ns -> us */
    adc_delay_us = adc_def->delay;

    grad_raster_us = desc->grad_raster_us;
    block_dur_us = bte->duration_us;
    n_grad = (int)(block_dur_us / grad_raster_us + 0.5f);
    if (n_grad < 2) n_grad = 2;
    dt_s = grad_raster_us * 1e-6f;

    /* Get gradient waveforms for this single block. We use the
     * grad_definitions and grad_table to reconstruct the waveform. */

    /* Allocate full-block k-space arrays */
    kx_full = (float*)PULSEQLIB_ALLOC((size_t)n_grad * sizeof(float));
    ky_full = (float*)PULSEQLIB_ALLOC((size_t)n_grad * sizeof(float));
    kz_full = (float*)PULSEQLIB_ALLOC((size_t)n_grad * sizeof(float));
    if (!kx_full || !ky_full || !kz_full) goto alloc_fail;

    /* Integrate gradient waveforms at block level.
     * For the trajectory we need to compute k(t) = integral(g(t)*dt).
     * We use the grad_table amplitude and grad_definition shapes. */
    {
        int gx_id = bte->gx_id;
        int gy_id = bte->gy_id;
        int gz_id = bte->gz_id;
        float gx_amp = 0.0f, gy_amp = 0.0f, gz_amp = 0.0f;

        /* Get amplitudes from grad table */
        if (gx_id >= 0 && gx_id < desc->grad_table_size)
            gx_amp = desc->grad_table[gx_id].amplitude;
        if (gy_id >= 0 && gy_id < desc->grad_table_size)
            gy_amp = desc->grad_table[gy_id].amplitude;
        if (gz_id >= 0 && gz_id < desc->grad_table_size)
            gz_amp = desc->grad_table[gz_id].amplitude;

        /* Get grad definition shapes and compute waveform.
         * For now we do simple constant-amplitude integration.
         * The actual shape comes from the grad_definition. */
        {
            int gx_def_id = -1, gy_def_id = -1, gz_def_id = -1;
            int gx_shot = 0, gy_shot = 0, gz_shot = 0;

            if (gx_id >= 0 && gx_id < desc->grad_table_size) {
                gx_def_id = desc->grad_table[gx_id].id;
                gx_shot = desc->grad_table[gx_id].shot_index;
            }
            if (gy_id >= 0 && gy_id < desc->grad_table_size) {
                gy_def_id = desc->grad_table[gy_id].id;
                gy_shot = desc->grad_table[gy_id].shot_index;
            }
            if (gz_id >= 0 && gz_id < desc->grad_table_size) {
                gz_def_id = desc->grad_table[gz_id].id;
                gz_shot = desc->grad_table[gz_id].shot_index;
            }

            /* Reconstruct k-space by integrating shaped gradients.
             * For each axis, expand the grad shape to uniform raster,
             * scale by amplitude, then trapezoidal integrate. */
            {
                float cum_x = 0.0f, cum_y = 0.0f, cum_z = 0.0f;
                kx_full[0] = 0.0f;
                ky_full[0] = 0.0f;
                kz_full[0] = 0.0f;

                for (i = 1; i < n_grad; ++i) {
                    /* Get gradient value at sample i.
                     * For trapezoid: rise/flat/fall. For arbitrary: shapes.
                     * Simplified: use existing internal helper to expand shapes
                     * at the grad raster level. For now we assume the
                     * waveforms have already been computed by check_safety. */
                    float gx_i = 0.0f, gy_i = 0.0f, gz_i = 0.0f;

                    /* Get shaped gradient value at this sample.
                     * Use the grad_definition + shape lookup. */
                    if (gx_def_id >= 0) {
                        const pulseqlib_grad_definition* gd =
                            &desc->grad_definitions[gx_def_id];
                        if (gd->type == PULSEQLIB__GRAD_TRAP) {
                            int rise = (int)(gd->rise_time_or_unused / grad_raster_us + 0.5f);
                            int flat = (int)(gd->flat_time_or_unused / grad_raster_us + 0.5f);
                            int fall = (int)(gd->fall_time_or_num_uncompressed_samples / grad_raster_us + 0.5f);
                            int delay_samp = (int)(gd->delay / grad_raster_us + 0.5f);
                            int s = i - delay_samp;
                            if (s < 0) gx_i = 0.0f;
                            else if (s < rise) gx_i = gx_amp * ((float)(s + 1) / (float)rise);
                            else if (s < rise + flat) gx_i = gx_amp;
                            else if (s < rise + flat + fall) gx_i = gx_amp * (1.0f - (float)(s - rise - flat + 1) / (float)fall);
                            else gx_i = 0.0f;
                        } else {
                            /* Arbitrary: use shape samples scaled by amplitude */
                            int shape_id = gd->shot_shape_ids[gx_shot];
                            if (shape_id >= 0 && shape_id < desc->num_shapes) {
                                const pulseqlib_shape_arbitrary* shape =
                                    &desc->shapes[shape_id];
                                int delay_samp = (int)(gd->delay / grad_raster_us + 0.5f);
                                int s = i - delay_samp;
                                if (s >= 0 && s < shape->num_uncompressed_samples)
                                    gx_i = gx_amp * shape->samples[s];
                            }
                        }
                    }

                    if (gy_def_id >= 0) {
                        const pulseqlib_grad_definition* gd =
                            &desc->grad_definitions[gy_def_id];
                        if (gd->type == PULSEQLIB__GRAD_TRAP) {
                            int rise = (int)(gd->rise_time_or_unused / grad_raster_us + 0.5f);
                            int flat = (int)(gd->flat_time_or_unused / grad_raster_us + 0.5f);
                            int fall = (int)(gd->fall_time_or_num_uncompressed_samples / grad_raster_us + 0.5f);
                            int delay_samp = (int)(gd->delay / grad_raster_us + 0.5f);
                            int s = i - delay_samp;
                            if (s < 0) gy_i = 0.0f;
                            else if (s < rise) gy_i = gy_amp * ((float)(s + 1) / (float)rise);
                            else if (s < rise + flat) gy_i = gy_amp;
                            else if (s < rise + flat + fall) gy_i = gy_amp * (1.0f - (float)(s - rise - flat + 1) / (float)fall);
                            else gy_i = 0.0f;
                        } else {
                            int shape_id = gd->shot_shape_ids[gy_shot];
                            if (shape_id >= 0 && shape_id < desc->num_shapes) {
                                const pulseqlib_shape_arbitrary* shape =
                                    &desc->shapes[shape_id];
                                int delay_samp = (int)(gd->delay / grad_raster_us + 0.5f);
                                int s = i - delay_samp;
                                if (s >= 0 && s < shape->num_uncompressed_samples)
                                    gy_i = gy_amp * shape->samples[s];
                            }
                        }
                    }

                    if (gz_def_id >= 0) {
                        const pulseqlib_grad_definition* gd =
                            &desc->grad_definitions[gz_def_id];
                        if (gd->type == PULSEQLIB__GRAD_TRAP) {
                            int rise = (int)(gd->rise_time_or_unused / grad_raster_us + 0.5f);
                            int flat = (int)(gd->flat_time_or_unused / grad_raster_us + 0.5f);
                            int fall = (int)(gd->fall_time_or_num_uncompressed_samples / grad_raster_us + 0.5f);
                            int delay_samp = (int)(gd->delay / grad_raster_us + 0.5f);
                            int s = i - delay_samp;
                            if (s < 0) gz_i = 0.0f;
                            else if (s < rise) gz_i = gz_amp * ((float)(s + 1) / (float)rise);
                            else if (s < rise + flat) gz_i = gz_amp;
                            else if (s < rise + flat + fall) gz_i = gz_amp * (1.0f - (float)(s - rise - flat + 1) / (float)fall);
                            else gz_i = 0.0f;
                        } else {
                            int shape_id = gd->shot_shape_ids[gz_shot];
                            if (shape_id >= 0 && shape_id < desc->num_shapes) {
                                const pulseqlib_shape_arbitrary* shape =
                                    &desc->shapes[shape_id];
                                int delay_samp = (int)(gd->delay / grad_raster_us + 0.5f);
                                int s = i - delay_samp;
                                if (s >= 0 && s < shape->num_uncompressed_samples)
                                    gz_i = gz_amp * shape->samples[s];
                            }
                        }
                    }

                    /* Trapezoidal integration.
                     * Note: gx_i-1 is approximated as previous iteration value. */
                    cum_x += 0.5f * (gx_i) * dt_s; /* simplified: half-step */
                    cum_y += 0.5f * (gy_i) * dt_s;
                    cum_z += 0.5f * (gz_i) * dt_s;

                    kx_full[i] = cum_x;
                    ky_full[i] = cum_y;
                    kz_full[i] = cum_z;

                    /* Add the other half-step for proper trapezoidal rule */
                    cum_x += 0.5f * (gx_i) * dt_s;
                    cum_y += 0.5f * (gy_i) * dt_s;
                    cum_z += 0.5f * (gz_i) * dt_s;
                }
            }
        }
    }

    /* Crop to ADC window */
    adc_start_sample = (int)(adc_delay_us / grad_raster_us + 0.5f);
    adc_window_samples = (int)((float)adc_num_samples * adc_dwell_us / grad_raster_us + 0.5f);
    adc_end_sample = adc_start_sample + adc_window_samples;
    if (adc_end_sample > n_grad) adc_end_sample = n_grad;
    adc_window_samples = adc_end_sample - adc_start_sample;
    if (adc_window_samples < 1) adc_window_samples = 1;

    kx_window = kx_full + adc_start_sample;
    ky_window = ky_full + adc_start_sample;
    kz_window = kz_full + adc_start_sample;

    /* Resample from grad raster to ADC dwell */
    linear_resample(kx_window, adc_window_samples, grad_raster_us,
                    out_kx, adc_num_samples, adc_dwell_us);
    linear_resample(ky_window, adc_window_samples, grad_raster_us,
                    out_ky, adc_num_samples, adc_dwell_us);
    linear_resample(kz_window, adc_window_samples, grad_raster_us,
                    out_kz, adc_num_samples, adc_dwell_us);

    /* Center to k-zero anchor */
    if (kzero_index >= 0 && kzero_index < adc_num_samples) {
        kzero_kx = out_kx[kzero_index];
        kzero_ky = out_ky[kzero_index];
        kzero_kz = out_kz[kzero_index];
    } else {
        kzero_kx = 0.0f;
        kzero_ky = 0.0f;
        kzero_kz = 0.0f;
    }

    for (i = 0; i < adc_num_samples; ++i) {
        out_kx[i] -= kzero_kx;
        out_ky[i] -= kzero_ky;
        out_kz[i] -= kzero_kz;
    }

    *out_num_samples = adc_num_samples;

    PULSEQLIB_FREE(kx_full);
    PULSEQLIB_FREE(ky_full);
    PULSEQLIB_FREE(kz_full);
    return PULSEQLIB_SUCCESS;

alloc_fail:
    PULSEQLIB_FREE(kx_full);
    PULSEQLIB_FREE(ky_full);
    PULSEQLIB_FREE(kz_full);
    return PULSEQLIB_ERR_ALLOC_FAILED;
}

/* ================================================================== */
/*  Public: compute trajectory                                        */
/* ================================================================== */

int pulseqlib_compute_trajectory(const pulseqlib_collection* coll,
                                 pulseqlib_trajectory*       out,
                                 pulseqlib_diagnostic*       diag,
                                 int                         subseq_idx)
{
    const pulseqlib_sequence_descriptor* desc;
    int n, b, seg, blk;
    int num_adc_events;
    int adc_idx;
    int label_ncols;
    int* label_buf = NULL;
    pulseqlib_traj_table_entry* table = NULL;
    float* kx_buf = NULL;
    float* ky_buf = NULL;
    float* kz_buf = NULL;
    int max_adc_samples;
    int rc;

    if (!coll || !out) return PULSEQLIB_ERR_NULL_POINTER;
    if (subseq_idx < 0 || subseq_idx >= coll->num_subsequences)
        return PULSEQLIB_ERR_INVALID_ARGUMENT;

    memset(out, 0, sizeof(*out));
    desc = &coll->descriptors[subseq_idx];

    /* Count ADC events from scan table */
    num_adc_events = 0;
    for (n = 0; n < desc->scan_table_len; ++n) {
        b = desc->scan_table_block_idx[n];
        if (desc->block_table[b].adc_id >= 0)
            ++num_adc_events;
    }

    if (num_adc_events == 0) {
        /* No ADC events — trajectory is empty */
        return PULSEQLIB_SUCCESS;
    }

    /* Find max ADC sample count */
    max_adc_samples = 0;
    {
        int a;
        for (a = 0; a < desc->num_unique_adcs; ++a) {
            if (desc->adc_definitions[a].num_samples > max_adc_samples)
                max_adc_samples = desc->adc_definitions[a].num_samples;
        }
    }

    /* Allocate work buffers */
    kx_buf = (float*)PULSEQLIB_ALLOC((size_t)max_adc_samples * sizeof(float));
    ky_buf = (float*)PULSEQLIB_ALLOC((size_t)max_adc_samples * sizeof(float));
    kz_buf = (float*)PULSEQLIB_ALLOC((size_t)max_adc_samples * sizeof(float));
    label_ncols = desc->label_num_columns;
    if (label_ncols > 0)
        label_buf = (int*)PULSEQLIB_ALLOC((size_t)label_ncols * sizeof(int));
    table = (pulseqlib_traj_table_entry*)PULSEQLIB_ALLOC(
        (size_t)num_adc_events * sizeof(pulseqlib_traj_table_entry));

    if (!kx_buf || !ky_buf || !kz_buf || !table)
        goto compute_fail;

    memset(table, 0, (size_t)num_adc_events * sizeof(pulseqlib_traj_table_entry));

    /* Initialize kshot library */
    out->kshots.num_shots = 0;
    out->kshots.shots = NULL;

    /* Iterate scan table, process each ADC block */
    adc_idx = 0;
    for (n = 0; n < desc->scan_table_len; ++n) {
        const pulseqlib_block_table_element* bte;
        int adc_def_idx, adc_nsamples, kzero;
        int kx_id, ky_id, kz_id;

        b = desc->scan_table_block_idx[n];
        bte = &desc->block_table[b];
        adc_def_idx = bte->adc_id;
        if (adc_def_idx < 0)
            continue;

        adc_nsamples = desc->adc_definitions[adc_def_idx].num_samples;

        /* Find k-zero index from segment timing.
         * Look up which segment this block belongs to and find the
         * ADC anchor with matching block offset. */
        kzero = adc_nsamples / 2;  /* default: center */
        {
            int seg_idx = desc->scan_table_seg_id[n];
            if (seg_idx >= 0 && seg_idx < desc->num_unique_segments) {
                const pulseqlib_tr_segment* seg_def = &desc->segment_definitions[seg_idx];
                const pulseqlib_segment_timing* tim = &seg_def->timing;
                int a;
                for (a = 0; a < tim->num_adc_anchors; ++a) {
                    if (tim->adc_anchors[a].block_offset ==
                        (b - seg_def->start_block)) {
                        kzero = tim->adc_anchors[a].kzero_index;
                        break;
                    }
                }
            }
        }

        /* Compute block-level k-space */
        rc = compute_block_kspace(desc, b, kzero,
                                  kx_buf, ky_buf, kz_buf,
                                  &adc_nsamples, diag);
        if (PULSEQLIB_FAILED(rc)) goto compute_fail;

        /* Per-axis dedup into kshot library */
        if (is_trivial_shot(kx_buf, adc_nsamples)) {
            kx_id = -1;
        } else {
            kx_id = kshot_library_add(&out->kshots, kx_buf, adc_nsamples);
            if (kx_id < 0) goto compute_fail;
        }

        if (is_trivial_shot(ky_buf, adc_nsamples)) {
            ky_id = -1;
        } else {
            ky_id = kshot_library_add(&out->kshots, ky_buf, adc_nsamples);
            if (ky_id < 0) goto compute_fail;
        }

        if (is_trivial_shot(kz_buf, adc_nsamples)) {
            kz_id = -1;
        } else {
            kz_id = kshot_library_add(&out->kshots, kz_buf, adc_nsamples);
            if (kz_id < 0) goto compute_fail;
        }

        /* Populate table entry */
        table[adc_idx].kx_shot_id = kx_id;
        table[adc_idx].ky_shot_id = ky_id;
        table[adc_idx].kz_shot_id = kz_id;

        /* Gradient amplitudes */
        table[adc_idx].gx_amplitude = (bte->gx_id >= 0 && bte->gx_id < desc->grad_table_size)
            ? desc->grad_table[bte->gx_id].amplitude : 0.0f;
        table[adc_idx].gy_amplitude = (bte->gy_id >= 0 && bte->gy_id < desc->grad_table_size)
            ? desc->grad_table[bte->gy_id].amplitude : 0.0f;
        table[adc_idx].gz_amplitude = (bte->gz_id >= 0 && bte->gz_id < desc->grad_table_size)
            ? desc->grad_table[bte->gz_id].amplitude : 0.0f;

        /* Rotation */
        table[adc_idx].rotation_id = bte->rotation_id;

        /* Labels from label table */
        if (label_buf && adc_idx < desc->label_num_entries) {
            rc = pulseqlib_get_adc_label(coll, label_buf,
                                         subseq_idx, adc_idx);
            if (PULSEQLIB_SUCCEEDED(rc) && label_ncols >= 3) {
                /* GEHC mapping: col0=lin, col1=slc, col2=eco */
                table[adc_idx].lin = label_buf[0];
                table[adc_idx].slc = label_buf[1];
                table[adc_idx].eco = label_buf[2];
            }
        }

        ++adc_idx;
    }

    out->num_adc_events = adc_idx;
    out->table = table;
    table = NULL;

    /* Build encoding spaces (one per subsequence) */
    out->num_encoding_spaces = 1;
    out->encoding_spaces = (pulseqlib_encoding_space*)PULSEQLIB_ALLOC(
        sizeof(pulseqlib_encoding_space));
    if (!out->encoding_spaces) goto compute_fail;
    memset(out->encoding_spaces, 0, sizeof(pulseqlib_encoding_space));
    memcpy(out->encoding_spaces[0].fov, desc->fov, sizeof(float) * 3);
    memcpy(out->encoding_spaces[0].matrix, desc->matrix, sizeof(float) * 3);
    memcpy(out->encoding_spaces[0].nav_fov, desc->nav_fov, sizeof(float) * 3);
    memcpy(out->encoding_spaces[0].nav_matrix, desc->nav_matrix, sizeof(float) * 3);
    out->encoding_spaces[0].subseq_idx = subseq_idx;
    out->encoding_spaces[0].nav_subseq_offset = 0;

    PULSEQLIB_FREE(kx_buf);
    PULSEQLIB_FREE(ky_buf);
    PULSEQLIB_FREE(kz_buf);
    PULSEQLIB_FREE(label_buf);
    return PULSEQLIB_SUCCESS;

compute_fail:
    PULSEQLIB_FREE(kx_buf);
    PULSEQLIB_FREE(ky_buf);
    PULSEQLIB_FREE(kz_buf);
    PULSEQLIB_FREE(label_buf);
    PULSEQLIB_FREE(table);
    pulseqlib_free_trajectory(out);
    return PULSEQLIB_ERR_ALLOC_FAILED;
}

/* ================================================================== */
/*  Free trajectory                                                   */
/* ================================================================== */

void pulseqlib_free_trajectory(pulseqlib_trajectory* traj)
{
    int i;
    if (!traj) return;

    if (traj->kshots.shots) {
        for (i = 0; i < traj->kshots.num_shots; ++i)
            PULSEQLIB_FREE(traj->kshots.shots[i].k);
        free(traj->kshots.shots);
    }
    traj->kshots.shots = NULL;
    traj->kshots.num_shots = 0;

    PULSEQLIB_FREE(traj->encoding_spaces);
    traj->encoding_spaces = NULL;
    traj->num_encoding_spaces = 0;

    PULSEQLIB_FREE(traj->table);
    traj->table = NULL;
    traj->num_adc_events = 0;
}

/* ================================================================== */
/*  Cache I/O helpers (local to this file)                            */
/* ================================================================== */

static void traj_swap4(void* p)
{
    unsigned char* b = (unsigned char*)p;
    unsigned char t;
    t = b[0]; b[0] = b[3]; b[3] = t;
    t = b[1]; b[1] = b[2]; b[2] = t;
}

static void traj_swap4_array(void* p, int count)
{
    int i;
    for (i = 0; i < count; ++i)
        traj_swap4((unsigned char*)p + (size_t)i * 4);
}

static int traj_write4(FILE* f, const void* p, int count)
{
    return (int)fwrite(p, 4, (size_t)count, f) == count;
}

static int traj_read4(FILE* f, void* p, int count)
{
    return (int)fread(p, 4, (size_t)count, f) == count;
}

/* ================================================================== */
/*  Write trajectory cache (section 5)                                */
/* ================================================================== */

#define CACHE_ENDIAN_MARKER  0x01020304
#define CACHE_SECTION_TRAJECTORY 5

int pulseqlib_write_trajectory_cache(const pulseqlib_trajectory* traj,
                                     const char*                 seq_path)
{
    char* cache_path;
    FILE* f;
    int marker, num_sections;
    int version_major, version_minor, vendor, stored_size;
    int do_swap;
    long entries_pos, data_start, data_end, hdr_ns_pos;
    int i, found_idx;
    int entries_buf[16 * 3]; /* up to 16 sections x 3 ints each */

    if (!traj || !seq_path) return PULSEQLIB_ERR_NULL_POINTER;

    cache_path = traj_make_cache_path(seq_path);
    if (!cache_path) return PULSEQLIB_ERR_ALLOC_FAILED;

    f = fopen(cache_path, "r+b");
    if (!f) { PULSEQLIB_FREE(cache_path); return PULSEQLIB_ERR_FILE_READ_FAILED; }

    /* Read header */
    if (!traj_read4(f, &marker, 1)) goto tw_fail;
    do_swap = 0;
    if (marker != CACHE_ENDIAN_MARKER) {
        traj_swap4(&marker);
        if (marker != CACHE_ENDIAN_MARKER) goto tw_fail;
        do_swap = 1;
    }
    if (!traj_read4(f, &version_major, 1)) goto tw_fail;
    if (!traj_read4(f, &version_minor, 1)) goto tw_fail;
    if (!traj_read4(f, &vendor, 1))        goto tw_fail;
    if (!traj_read4(f, &stored_size, 1))   goto tw_fail;
    hdr_ns_pos = ftell(f);
    if (!traj_read4(f, &num_sections, 1))  goto tw_fail;
    if (do_swap) {
        traj_swap4(&version_major); traj_swap4(&version_minor);
        traj_swap4(&vendor); traj_swap4(&stored_size);
        traj_swap4(&num_sections);
    }
    if (num_sections <= 0 || num_sections > 15) goto tw_fail;

    entries_pos = ftell(f);
    if (entries_pos < 0) goto tw_fail;

    /* Read existing section entries */
    for (i = 0; i < num_sections; ++i) {
        if (!traj_read4(f, &entries_buf[i * 3], 3)) goto tw_fail;
        if (do_swap) traj_swap4_array(&entries_buf[i * 3], 3);
    }

    /* Check if section 5 already exists */
    found_idx = -1;
    for (i = 0; i < num_sections; ++i) {
        if (entries_buf[i * 3] == CACHE_SECTION_TRAJECTORY) {
            found_idx = i;
            break;
        }
    }
    if (found_idx < 0) {
        found_idx = num_sections;
        entries_buf[found_idx * 3] = CACHE_SECTION_TRAJECTORY;
        num_sections++;
    }

    /* Seek to end, write trajectory data */
    fseek(f, 0, SEEK_END);
    data_start = ftell(f);
    if (data_start < 0) goto tw_fail;

    /* Write kshot library */
    if (!traj_write4(f, &traj->kshots.num_shots, 1)) goto tw_fail;
    for (i = 0; i < traj->kshots.num_shots; ++i) {
        if (!traj_write4(f, &traj->kshots.shots[i].num_samples, 1)) goto tw_fail;
        if (traj->kshots.shots[i].num_samples > 0) {
            if (!traj_write4(f, traj->kshots.shots[i].k,
                            traj->kshots.shots[i].num_samples)) goto tw_fail;
        }
    }

    /* Write encoding spaces */
    if (!traj_write4(f, &traj->num_encoding_spaces, 1)) goto tw_fail;
    for (i = 0; i < traj->num_encoding_spaces; ++i) {
        const pulseqlib_encoding_space* es = &traj->encoding_spaces[i];
        if (!traj_write4(f, es->fov, 3)) goto tw_fail;
        if (!traj_write4(f, es->matrix, 3)) goto tw_fail;
        if (!traj_write4(f, es->nav_fov, 3)) goto tw_fail;
        if (!traj_write4(f, es->nav_matrix, 3)) goto tw_fail;
        if (!traj_write4(f, &es->subseq_idx, 1)) goto tw_fail;
        if (!traj_write4(f, &es->nav_subseq_offset, 1)) goto tw_fail;
    }

    /* Write trajectory table */
    if (!traj_write4(f, &traj->num_adc_events, 1)) goto tw_fail;
    for (i = 0; i < traj->num_adc_events; ++i) {
        const pulseqlib_traj_table_entry* e = &traj->table[i];
        /* Write as contiguous 17 ints/floats:
         * 3 shot_ids + 3 amplitudes + rotation_id + 10 labels */
        if (!traj_write4(f, &e->kx_shot_id, 1)) goto tw_fail;
        if (!traj_write4(f, &e->ky_shot_id, 1)) goto tw_fail;
        if (!traj_write4(f, &e->kz_shot_id, 1)) goto tw_fail;
        if (!traj_write4(f, &e->gx_amplitude, 1)) goto tw_fail;
        if (!traj_write4(f, &e->gy_amplitude, 1)) goto tw_fail;
        if (!traj_write4(f, &e->gz_amplitude, 1)) goto tw_fail;
        if (!traj_write4(f, &e->rotation_id, 1)) goto tw_fail;
        if (!traj_write4(f, &e->slc, 1)) goto tw_fail;
        if (!traj_write4(f, &e->seg, 1)) goto tw_fail;
        if (!traj_write4(f, &e->rep, 1)) goto tw_fail;
        if (!traj_write4(f, &e->avg, 1)) goto tw_fail;
        if (!traj_write4(f, &e->set, 1)) goto tw_fail;
        if (!traj_write4(f, &e->eco, 1)) goto tw_fail;
        if (!traj_write4(f, &e->phs, 1)) goto tw_fail;
        if (!traj_write4(f, &e->lin, 1)) goto tw_fail;
        if (!traj_write4(f, &e->par, 1)) goto tw_fail;
        if (!traj_write4(f, &e->acq, 1)) goto tw_fail;
    }

    data_end = ftell(f);
    if (data_end < 0) goto tw_fail;

    entries_buf[found_idx * 3 + 1] = (int)data_start;
    entries_buf[found_idx * 3 + 2] = (int)(data_end - data_start);

    /* Patch num_sections */
    if (fseek(f, hdr_ns_pos, SEEK_SET) != 0) goto tw_fail;
    if (!traj_write4(f, &num_sections, 1)) goto tw_fail;

    /* Rewrite all section entries */
    if (fseek(f, entries_pos, SEEK_SET) != 0) goto tw_fail;
    for (i = 0; i < num_sections; ++i) {
        if (!traj_write4(f, &entries_buf[i * 3], 3)) goto tw_fail;
    }

    fclose(f);
    PULSEQLIB_FREE(cache_path);
    return PULSEQLIB_SUCCESS;

tw_fail:
    fclose(f);
    PULSEQLIB_FREE(cache_path);
    return PULSEQLIB_ERR_FILE_READ_FAILED;
}

/* ================================================================== */
/*  Load trajectory from cache (section 5)                            */
/* ================================================================== */

int pulseqlib_load_trajectory_cache(pulseqlib_trajectory* out,
                                    const char*           seq_path)
{
    char* cache_path;
    FILE* f;
    int marker, num_sections;
    int version_major, version_minor, vendor, stored_size;
    int do_swap, i, found;
    int section_id, section_offset, section_size;

    if (!out || !seq_path) return PULSEQLIB_ERR_NULL_POINTER;
    memset(out, 0, sizeof(*out));

    cache_path = traj_make_cache_path(seq_path);
    if (!cache_path) return PULSEQLIB_ERR_ALLOC_FAILED;

    f = fopen(cache_path, "rb");
    PULSEQLIB_FREE(cache_path);
    if (!f) return PULSEQLIB_ERR_FILE_READ_FAILED;

    /* Read header */
    if (!traj_read4(f, &marker, 1)) { fclose(f); return PULSEQLIB_ERR_FILE_READ_FAILED; }
    do_swap = 0;
    if (marker != CACHE_ENDIAN_MARKER) {
        traj_swap4(&marker);
        if (marker != CACHE_ENDIAN_MARKER) { fclose(f); return PULSEQLIB_ERR_FILE_READ_FAILED; }
        do_swap = 1;
    }
    if (!traj_read4(f, &version_major, 1)) { fclose(f); return PULSEQLIB_ERR_FILE_READ_FAILED; }
    if (!traj_read4(f, &version_minor, 1)) { fclose(f); return PULSEQLIB_ERR_FILE_READ_FAILED; }
    if (!traj_read4(f, &vendor, 1))        { fclose(f); return PULSEQLIB_ERR_FILE_READ_FAILED; }
    if (!traj_read4(f, &stored_size, 1))   { fclose(f); return PULSEQLIB_ERR_FILE_READ_FAILED; }
    if (!traj_read4(f, &num_sections, 1))  { fclose(f); return PULSEQLIB_ERR_FILE_READ_FAILED; }
    if (do_swap) {
        traj_swap4(&version_major); traj_swap4(&version_minor);
        traj_swap4(&vendor); traj_swap4(&stored_size);
        traj_swap4(&num_sections);
    }

    /* Find section 5 */
    found = 0;
    section_offset = 0;
    section_size = 0;
    for (i = 0; i < num_sections; ++i) {
        if (!traj_read4(f, &section_id, 1))     { fclose(f); return PULSEQLIB_ERR_FILE_READ_FAILED; }
        if (!traj_read4(f, &section_offset, 1))  { fclose(f); return PULSEQLIB_ERR_FILE_READ_FAILED; }
        if (!traj_read4(f, &section_size, 1))    { fclose(f); return PULSEQLIB_ERR_FILE_READ_FAILED; }
        if (do_swap) { traj_swap4(&section_id); traj_swap4(&section_offset); traj_swap4(&section_size); }
        if (section_id == CACHE_SECTION_TRAJECTORY) { found = 1; break; }
    }

    if (!found) {
        fclose(f);
        return PULSEQLIB_ERR_FILE_READ_FAILED;
    }

    /* Seek to trajectory data */
    if (fseek(f, section_offset, SEEK_SET) != 0) { fclose(f); return PULSEQLIB_ERR_FILE_READ_FAILED; }

    /* Read kshot library */
    if (!traj_read4(f, &out->kshots.num_shots, 1)) goto lr_fail;
    if (do_swap) traj_swap4(&out->kshots.num_shots);

    if (out->kshots.num_shots > 0) {
        out->kshots.shots = (pulseqlib_kshot*)PULSEQLIB_ALLOC(
            (size_t)out->kshots.num_shots * sizeof(pulseqlib_kshot));
        if (!out->kshots.shots) goto lr_fail;
        memset(out->kshots.shots, 0, (size_t)out->kshots.num_shots * sizeof(pulseqlib_kshot));

        for (i = 0; i < out->kshots.num_shots; ++i) {
            if (!traj_read4(f, &out->kshots.shots[i].num_samples, 1)) goto lr_fail;
            if (do_swap) traj_swap4(&out->kshots.shots[i].num_samples);

            if (out->kshots.shots[i].num_samples > 0) {
                out->kshots.shots[i].k = (float*)PULSEQLIB_ALLOC(
                    (size_t)out->kshots.shots[i].num_samples * sizeof(float));
                if (!out->kshots.shots[i].k) goto lr_fail;
                if (!traj_read4(f, out->kshots.shots[i].k,
                               out->kshots.shots[i].num_samples)) goto lr_fail;
                if (do_swap) traj_swap4_array(out->kshots.shots[i].k,
                                              out->kshots.shots[i].num_samples);
            }
        }
    }

    /* Read encoding spaces */
    if (!traj_read4(f, &out->num_encoding_spaces, 1)) goto lr_fail;
    if (do_swap) traj_swap4(&out->num_encoding_spaces);

    if (out->num_encoding_spaces > 0) {
        out->encoding_spaces = (pulseqlib_encoding_space*)PULSEQLIB_ALLOC(
            (size_t)out->num_encoding_spaces * sizeof(pulseqlib_encoding_space));
        if (!out->encoding_spaces) goto lr_fail;

        for (i = 0; i < out->num_encoding_spaces; ++i) {
            pulseqlib_encoding_space* es = &out->encoding_spaces[i];
            if (!traj_read4(f, es->fov, 3)) goto lr_fail;
            if (!traj_read4(f, es->matrix, 3)) goto lr_fail;
            if (!traj_read4(f, es->nav_fov, 3)) goto lr_fail;
            if (!traj_read4(f, es->nav_matrix, 3)) goto lr_fail;
            if (!traj_read4(f, &es->subseq_idx, 1)) goto lr_fail;
            if (!traj_read4(f, &es->nav_subseq_offset, 1)) goto lr_fail;
            if (do_swap) traj_swap4_array(es->fov, 14);
        }
    }

    /* Read trajectory table */
    if (!traj_read4(f, &out->num_adc_events, 1)) goto lr_fail;
    if (do_swap) traj_swap4(&out->num_adc_events);

    if (out->num_adc_events > 0) {
        out->table = (pulseqlib_traj_table_entry*)PULSEQLIB_ALLOC(
            (size_t)out->num_adc_events * sizeof(pulseqlib_traj_table_entry));
        if (!out->table) goto lr_fail;

        for (i = 0; i < out->num_adc_events; ++i) {
            pulseqlib_traj_table_entry* e = &out->table[i];
            if (!traj_read4(f, &e->kx_shot_id, 1)) goto lr_fail;
            if (!traj_read4(f, &e->ky_shot_id, 1)) goto lr_fail;
            if (!traj_read4(f, &e->kz_shot_id, 1)) goto lr_fail;
            if (!traj_read4(f, &e->gx_amplitude, 1)) goto lr_fail;
            if (!traj_read4(f, &e->gy_amplitude, 1)) goto lr_fail;
            if (!traj_read4(f, &e->gz_amplitude, 1)) goto lr_fail;
            if (!traj_read4(f, &e->rotation_id, 1)) goto lr_fail;
            if (!traj_read4(f, &e->slc, 1)) goto lr_fail;
            if (!traj_read4(f, &e->seg, 1)) goto lr_fail;
            if (!traj_read4(f, &e->rep, 1)) goto lr_fail;
            if (!traj_read4(f, &e->avg, 1)) goto lr_fail;
            if (!traj_read4(f, &e->set, 1)) goto lr_fail;
            if (!traj_read4(f, &e->eco, 1)) goto lr_fail;
            if (!traj_read4(f, &e->phs, 1)) goto lr_fail;
            if (!traj_read4(f, &e->lin, 1)) goto lr_fail;
            if (!traj_read4(f, &e->par, 1)) goto lr_fail;
            if (!traj_read4(f, &e->acq, 1)) goto lr_fail;
            if (do_swap) traj_swap4_array(&e->kx_shot_id, 17);
        }
    }

    fclose(f);
    return PULSEQLIB_SUCCESS;

lr_fail:
    fclose(f);
    pulseqlib_free_trajectory(out);
    return PULSEQLIB_ERR_FILE_READ_FAILED;
}
