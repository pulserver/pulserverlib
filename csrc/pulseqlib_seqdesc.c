/* pulseqlib_seqdesc.c -- Sequence description walker (Section 6)
 *
 * Implements:
 *   pulseqlib_get_sequence_description()
 *   pulseqlib_free_sequence_description()
 *   pulseqlib_get_sequence_parameters()
 *
 * The sequence description provides a per-TR event list (WAIT / RF / ADC)
 * plus shape/shim libraries and composite RF group annotations, intended
 * as input for a state-machine Bloch simulator on the reconstruction server.
 */

#include "pulseqlib_internal.h"
#include <string.h>

/* ================================================================== */
/*  Internal helpers                                                  */
/* ================================================================== */

/*
 * Compute TR-relative start time (us) for each block in the imaging TR.
 * Returns a malloc'd array of length tr_size; caller must free.
 */
static float* seqdesc__build_block_start_us(
    const pulseqlib_sequence_descriptor* desc)
{
    const pulseqlib_tr_descriptor* trd = &desc->tr_descriptor;
    int    n = trd->tr_size;
    float* t;
    float  acc;
    int    i, blk;

    if (n <= 0) return NULL;
    t = (float*)PULSEQLIB_ALLOC((size_t)n * sizeof(float));
    if (!t) return NULL;

    acc = 0.0f;
    for (i = 0; i < n; ++i) {
        t[i] = acc;
        blk  = trd->imaging_tr_start + i;
        if (blk >= 0 && blk < desc->num_blocks)
            acc += (float)desc->block_table[blk].duration_us;
    }
    return t;
}

/*
 * Build a TR-relative kzero_us lookup indexed by TR block position.
 * Uses the segment timing anchors; values are -1 where no ADC anchor exists.
 * Returns a malloc'd array of length tr_size; caller must free.
 */
static float* seqdesc__build_adc_kzero_us(
    const pulseqlib_sequence_descriptor* desc,
    const float*                         blk_start_us,
    int                                  tr_size)
{
    const pulseqlib_segment_table_result* stbl = &desc->segment_table;
    float*                       kzero;
    int                          si, seg_id, ai;
    const pulseqlib_tr_segment*  seg;

    kzero = (float*)PULSEQLIB_ALLOC((size_t)tr_size * sizeof(float));
    if (!kzero) return NULL;

    for (si = 0; si < tr_size; ++si) kzero[si] = -1.0f;

    if (!stbl->main_segment_table || stbl->num_main_segments <= 0)
        return kzero;   /* fallback to midpoint will apply everywhere */

    for (si = 0; si < stbl->num_main_segments; ++si) {
        int tr_pos;
        seg_id = stbl->main_segment_table[si];
        if (seg_id < 0 || seg_id >= desc->num_unique_segments) continue;
        seg = &desc->segment_definitions[seg_id];
        if (!seg->timing.adc_anchors) continue;

        for (ai = 0; ai < seg->timing.num_adc_anchors; ++ai) {
            tr_pos = seg->start_block + seg->timing.adc_anchors[ai].block_offset;
            if (tr_pos < 0 || tr_pos >= tr_size) continue;
            if (seg->start_block < 0 || seg->start_block >= tr_size) continue;
            kzero[tr_pos] = blk_start_us[seg->start_block]
                          + seg->timing.adc_anchors[ai].kzero_us;
        }
    }
    return kzero;
}

/* ================================================================== */
/*  pulseqlib_free_sequence_description                               */
/* ================================================================== */

void pulseqlib_free_sequence_description(pulseqlib_sequence_description* desc)
{
    int i;
    if (!desc) return;

    if (desc->rf_shape_tuples) {
        for (i = 0; i < desc->num_tuples; ++i) {
            if (desc->rf_shape_tuples[i].mag)   PULSEQLIB_FREE(desc->rf_shape_tuples[i].mag);
            if (desc->rf_shape_tuples[i].phase) PULSEQLIB_FREE(desc->rf_shape_tuples[i].phase);
            if (desc->rf_shape_tuples[i].time)  PULSEQLIB_FREE(desc->rf_shape_tuples[i].time);
        }
        PULSEQLIB_FREE(desc->rf_shape_tuples);
    }
    if (desc->shim_defs)           PULSEQLIB_FREE(desc->shim_defs);
    if (desc->events)              PULSEQLIB_FREE(desc->events);
    if (desc->composite_rf_groups) PULSEQLIB_FREE(desc->composite_rf_groups);

    memset(desc, 0, sizeof(*desc));
}

/* ================================================================== */
/*  pulseqlib_get_sequence_description                                */
/* ================================================================== */

int pulseqlib_get_sequence_description(
    pulseqlib_sequence_description* out,
    const pulseqlib_collection*     coll,
    int                             subseq_idx)
{
    const pulseqlib_sequence_descriptor* desc;
    const pulseqlib_tr_descriptor*       trd;
    const pulseqlib_block_table_element* bte;
    const pulseqlib_rf_definition*       rfdef;
    const pulseqlib_adc_definition*      adcdef;

    float*                  blk_start   = NULL;
    float*                  adc_kzero   = NULL;
    int*                    rf_def_used = NULL;  /* [num_unique_rfs] -> tuple_id */
    pulseqlib_seq_event*    tmp_events  = NULL;
    int                     max_events, n_events, n_tuples, n_shims;
    int                     tr_size, i, j, rf_def_id, adc_def_id;
    float                   t_cursor, blk_s;
    int                     ret = PULSEQLIB_SUCCESS;
    pulseqlib_shape_arbitrary decomp;

    if (!out || !coll) return PULSEQLIB_ERR_NULL_POINTER;
    if (subseq_idx < 0 || subseq_idx >= coll->num_subsequences)
        return PULSEQLIB_ERR_INVALID_ARGUMENT;

    memset(out, 0, sizeof(*out));
    out->subseq_idx = subseq_idx;

    desc    = &coll->descriptors[subseq_idx];
    trd     = &desc->tr_descriptor;
    tr_size = trd->tr_size;
    out->tr_duration_us = trd->tr_duration_us;

    if (tr_size <= 0)
        return PULSEQLIB_SUCCESS;   /* empty TR — legal, return blank desc */

    /* ===========================================================
     * Step 1: Block start times within canonical TR
     * =========================================================== */
    blk_start = seqdesc__build_block_start_us(desc);
    if (!blk_start) { ret = PULSEQLIB_ERR_ALLOC_FAILED; goto cleanup; }

    /* ===========================================================
     * Step 2: ADC kzero lookup (TR-relative)
     * =========================================================== */
    adc_kzero = seqdesc__build_adc_kzero_us(desc, blk_start, tr_size);
    if (!adc_kzero) { ret = PULSEQLIB_ERR_ALLOC_FAILED; goto cleanup; }

    /* ===========================================================
     * Step 3: Collect unique RF definition IDs used in TR
     * =========================================================== */
    if (desc->num_unique_rfs > 0) {
        rf_def_used = (int*)PULSEQLIB_ALLOC(
            (size_t)desc->num_unique_rfs * sizeof(int));
        if (!rf_def_used) { ret = PULSEQLIB_ERR_ALLOC_FAILED; goto cleanup; }
        for (i = 0; i < desc->num_unique_rfs; ++i) rf_def_used[i] = -1;
    }

    n_tuples = 0;
    for (i = 0; i < tr_size; ++i) {
        int blk = trd->imaging_tr_start + i;
        if (blk < 0 || blk >= desc->num_blocks) continue;
        bte = &desc->block_table[blk];
        if (bte->rf_id < 0 || bte->rf_id >= desc->rf_table_size) continue;
        rf_def_id = desc->rf_table[bte->rf_id].id;
        if (rf_def_id < 0 || rf_def_id >= desc->num_unique_rfs) continue;
        if (rf_def_used && rf_def_used[rf_def_id] < 0)
            rf_def_used[rf_def_id] = n_tuples++;
    }

    /* ===========================================================
     * Step 4: Decompress RF shapes into shape tuples
     * =========================================================== */
    if (n_tuples > 0) {
        out->rf_shape_tuples = (pulseqlib_rf_shape_tuple*)PULSEQLIB_ALLOC(
            (size_t)n_tuples * sizeof(pulseqlib_rf_shape_tuple));
        if (!out->rf_shape_tuples) { ret = PULSEQLIB_ERR_ALLOC_FAILED; goto cleanup; }
        memset(out->rf_shape_tuples, 0,
               (size_t)n_tuples * sizeof(pulseqlib_rf_shape_tuple));

        for (rf_def_id = 0; rf_def_id < desc->num_unique_rfs; ++rf_def_id) {
            int    tuple_id;
            pulseqlib_rf_shape_tuple* tup;
            int    mag_id, phase_id, time_id;
            int    nch, npts;

            if (!rf_def_used || rf_def_used[rf_def_id] < 0) continue;
            tuple_id = rf_def_used[rf_def_id];
            tup      = &out->rf_shape_tuples[tuple_id];
            rfdef    = &desc->rf_definitions[rf_def_id];

            nch = (rfdef->num_channels > 1) ? rfdef->num_channels : 1;
            tup->tuple_id   = tuple_id;
            tup->N_tx       = nch;
            tup->rf_raster_us = desc->rf_raster_us;

            mag_id   = rfdef->mag_shape_id;
            phase_id = rfdef->phase_shape_id;
            time_id  = rfdef->time_shape_id;

            /* Magnitude */
            if (mag_id > 0 && mag_id <= desc->num_shapes) {
                decomp.num_samples = 0; decomp.num_uncompressed_samples = 0;
                decomp.samples = NULL;
                if (!pulseqlib__decompress_shape(&decomp,
                        &desc->shapes[mag_id - 1], 1.0f)) {
                    ret = PULSEQLIB_ERR_ALLOC_FAILED; goto cleanup;
                }
                npts = decomp.num_uncompressed_samples / nch;
                tup->N_samples = npts;
                tup->mag = (float*)PULSEQLIB_ALLOC(
                    (size_t)decomp.num_uncompressed_samples * sizeof(float));
                if (!tup->mag) {
                    PULSEQLIB_FREE(decomp.samples);
                    ret = PULSEQLIB_ERR_ALLOC_FAILED; goto cleanup;
                }
                memcpy(tup->mag, decomp.samples,
                       (size_t)decomp.num_uncompressed_samples * sizeof(float));
                PULSEQLIB_FREE(decomp.samples);
            }

            /* Phase (optional) */
            if (phase_id > 0 && phase_id <= desc->num_shapes) {
                decomp.num_samples = 0; decomp.num_uncompressed_samples = 0;
                decomp.samples = NULL;
                if (!pulseqlib__decompress_shape(&decomp,
                        &desc->shapes[phase_id - 1], 1.0f)) {
                    ret = PULSEQLIB_ERR_ALLOC_FAILED; goto cleanup;
                }
                if (decomp.num_uncompressed_samples > 0) {
                    tup->phase = (float*)PULSEQLIB_ALLOC(
                        (size_t)decomp.num_uncompressed_samples * sizeof(float));
                    if (!tup->phase) {
                        PULSEQLIB_FREE(decomp.samples);
                        ret = PULSEQLIB_ERR_ALLOC_FAILED; goto cleanup;
                    }
                    memcpy(tup->phase, decomp.samples,
                           (size_t)decomp.num_uncompressed_samples * sizeof(float));
                }
                PULSEQLIB_FREE(decomp.samples);
            }

            /* Time (optional) */
            if (time_id > 0 && time_id <= desc->num_shapes) {
                decomp.num_samples = 0; decomp.num_uncompressed_samples = 0;
                decomp.samples = NULL;
                if (!pulseqlib__decompress_shape(&decomp,
                        &desc->shapes[time_id - 1], desc->rf_raster_us)) {
                    ret = PULSEQLIB_ERR_ALLOC_FAILED; goto cleanup;
                }
                if (decomp.num_uncompressed_samples > 0) {
                    tup->time = (float*)PULSEQLIB_ALLOC(
                        (size_t)decomp.num_uncompressed_samples * sizeof(float));
                    if (!tup->time) {
                        PULSEQLIB_FREE(decomp.samples);
                        ret = PULSEQLIB_ERR_ALLOC_FAILED; goto cleanup;
                    }
                    memcpy(tup->time, decomp.samples,
                           (size_t)decomp.num_uncompressed_samples * sizeof(float));
                }
                PULSEQLIB_FREE(decomp.samples);
            }

            /* Multiband annotation from rf_stats */
            tup->num_bands        = rfdef->stats.num_bands;
            tup->band_bandwidth_hz = rfdef->stats.band_bandwidth_hz;
            tup->total_b1sq_power  = rfdef->stats.total_b1sq_power;
            for (j = 0; j < PULSEQLIB_MAX_BANDS; ++j)
                tup->band_freq_offsets_hz[j] = rfdef->stats.band_freq_offsets_hz[j];
        }
    }
    out->num_tuples = n_tuples;

    /* ===========================================================
     * Step 5: Build per-subsequence shim library
     * =========================================================== */
    n_shims = desc->num_rf_shims;
    if (n_shims > 0) {
        out->shim_defs = (pulseqlib_shim_def_local*)PULSEQLIB_ALLOC(
            (size_t)n_shims * sizeof(pulseqlib_shim_def_local));
        if (!out->shim_defs) { ret = PULSEQLIB_ERR_ALLOC_FAILED; goto cleanup; }
        memset(out->shim_defs, 0, (size_t)n_shims * sizeof(pulseqlib_shim_def_local));

        for (i = 0; i < n_shims; ++i) {
            const pulseqlib_rf_shim_definition* src = &desc->rf_shim_definitions[i];
            out->shim_defs[i].shim_id_local = i;
            out->shim_defs[i].N_ch          = src->num_channels;
            for (j = 0; j < src->num_channels && j < PULSEQLIB_MAX_RF_SHIM_CHANNELS; ++j) {
                out->shim_defs[i].magnitudes[j] = src->magnitudes[j];
                out->shim_defs[i].phases[j]     = src->phases[j];
            }
        }
        out->num_shims = n_shims;
    }

    /* ===========================================================
     * Step 6: Walk TR blocks — emit WAIT / RF / ADC events
     * =========================================================== */
    max_events = 3 * tr_size + 4;
    tmp_events = (pulseqlib_seq_event*)PULSEQLIB_ALLOC(
        (size_t)max_events * sizeof(pulseqlib_seq_event));
    if (!tmp_events) { ret = PULSEQLIB_ERR_ALLOC_FAILED; goto cleanup; }
    memset(tmp_events, 0, (size_t)max_events * sizeof(pulseqlib_seq_event));
    n_events = 0;
    t_cursor = 0.0f;

    for (i = 0; i < tr_size; ++i) {
        int blk_abs = trd->imaging_tr_start + i;
        int has_rf, has_adc;

        if (blk_abs < 0 || blk_abs >= desc->num_blocks) continue;
        bte   = &desc->block_table[blk_abs];
        blk_s = blk_start[i];

        has_rf  = (bte->rf_id  >= 0 && bte->rf_id  < desc->rf_table_size);
        has_adc = (bte->adc_id >= 0 && bte->adc_id < desc->adc_table_size);

        /* ---- RF event ---- */
        if (has_rf) {
            float rf_start_us, rf_end_us, rf_center_us, rf_dur_us;
            int   tuple_id_local;

            rf_def_id = desc->rf_table[bte->rf_id].id;
            if (rf_def_id >= 0 && rf_def_id < desc->num_unique_rfs) {
                rfdef = &desc->rf_definitions[rf_def_id];

                rf_dur_us    = rfdef->stats.duration_us;
                rf_start_us  = blk_s + (float)rfdef->delay;
                rf_end_us    = rf_start_us + rf_dur_us;
                rf_center_us = rf_start_us + rf_dur_us
                             - (float)rfdef->stats.isodelay_us;

                /* WAIT before RF */
                if (rf_start_us > t_cursor + 0.5f && n_events < max_events) {
                    tmp_events[n_events].type       = PULSEQLIB_SEQ_EVENT_WAIT;
                    tmp_events[n_events].params[0]  = rf_start_us - t_cursor;
                    n_events++;
                }
                t_cursor = rf_start_us;

                tuple_id_local = (rf_def_used && rf_def_id < desc->num_unique_rfs)
                               ? rf_def_used[rf_def_id] : -1;

                if (n_events < max_events) {
                    tmp_events[n_events].type       = PULSEQLIB_SEQ_EVENT_RF;
                    tmp_events[n_events].params[0]  = rf_center_us;
                    tmp_events[n_events].params[1]  = rf_dur_us;
                    tmp_events[n_events].params[2]  = rfdef->stats.flip_angle_deg;
                    tmp_events[n_events].params[3]  = (float)tuple_id_local;
                    tmp_events[n_events].params[4]  = (float)bte->rf_shim_id;
                    tmp_events[n_events].params[5]  = (bte->gz_id >= 0) ? 1.0f : 0.0f;
                    tmp_events[n_events].params[6]  = 0.0f;
                    n_events++;
                }
                t_cursor = rf_end_us;
            }
        }

        /* ---- ADC event ---- */
        if (has_adc) {
            float adc_start_us, adc_end_us, adc_center_us, adc_dur_us, dwell_us;
            int   n_samps, adc_role;

            adc_def_id = desc->adc_table[bte->adc_id].id;
            if (adc_def_id >= 0 && adc_def_id < desc->num_unique_adcs) {
                adcdef   = &desc->adc_definitions[adc_def_id];
                n_samps  = adcdef->num_samples;
                dwell_us = (float)adcdef->dwell_time * 1e-3f;   /* ns → us */
                adc_dur_us  = (float)n_samps * dwell_us;
                adc_start_us = blk_s + (float)adcdef->delay;
                adc_end_us   = adc_start_us + adc_dur_us;

                /* TR-relative kzero from segment anchor, or midpoint fallback */
                adc_center_us = (adc_kzero[i] >= 0.0f)
                              ? adc_kzero[i]
                              : (adc_start_us + 0.5f * adc_dur_us);

                /* WAIT before ADC */
                if (adc_start_us > t_cursor + 0.5f && n_events < max_events) {
                    tmp_events[n_events].type      = PULSEQLIB_SEQ_EVENT_WAIT;
                    tmp_events[n_events].params[0] = adc_start_us - t_cursor;
                    n_events++;
                }
                t_cursor = adc_start_us;

                adc_role = (bte->pmc_flag || bte->nav_flag)
                         ? PULSEQLIB_ADC_ROLE_NON_ACQUIRED
                         : PULSEQLIB_ADC_ROLE_SINGLE;  /* refined below */

                if (n_events < max_events) {
                    tmp_events[n_events].type      = PULSEQLIB_SEQ_EVENT_ADC;
                    tmp_events[n_events].params[0] = adc_center_us;
                    tmp_events[n_events].params[1] = adc_dur_us;
                    tmp_events[n_events].params[2] = (float)n_samps;
                    tmp_events[n_events].params[3] = dwell_us;
                    tmp_events[n_events].params[4] = (float)adc_role;
                    tmp_events[n_events].params[5] = 0.0f;
                    tmp_events[n_events].params[6] = 0.0f;
                    n_events++;
                }
                t_cursor = adc_end_us;
            }
        }
    }

    /* Trailing WAIT to fill TR */
    if (out->tr_duration_us > t_cursor + 0.5f && n_events < max_events) {
        tmp_events[n_events].type      = PULSEQLIB_SEQ_EVENT_WAIT;
        tmp_events[n_events].params[0] = out->tr_duration_us - t_cursor;
        n_events++;
    }

    /* Copy trimmed event list */
    if (n_events > 0) {
        out->events = (pulseqlib_seq_event*)PULSEQLIB_ALLOC(
            (size_t)n_events * sizeof(pulseqlib_seq_event));
        if (!out->events) { ret = PULSEQLIB_ERR_ALLOC_FAILED; goto cleanup; }
        memcpy(out->events, tmp_events,
               (size_t)n_events * sizeof(pulseqlib_seq_event));
    }
    out->num_events = n_events;

    /* ===========================================================
     * Step 7: ADC echo-role refinement
     *   If multiple acquired ADC events exist, mark the first as
     *   ECHO_CENTER and the rest as NON_CENTER.
     *   (Full trajectory-based detection is deferred to recon server.)
     * =========================================================== */
    {
        int n_acq = 0, first_acq = -1;
        for (i = 0; i < n_events; ++i) {
            if (out->events[i].type == PULSEQLIB_SEQ_EVENT_ADC &&
                (int)out->events[i].params[4] != PULSEQLIB_ADC_ROLE_NON_ACQUIRED) {
                n_acq++;
                if (first_acq < 0) first_acq = i;
            }
        }
        if (n_acq > 1) {
            for (i = 0; i < n_events; ++i) {
                if (out->events[i].type == PULSEQLIB_SEQ_EVENT_ADC &&
                    (int)out->events[i].params[4] != PULSEQLIB_ADC_ROLE_NON_ACQUIRED) {
                    out->events[i].params[4] = (i == first_acq)
                        ? (float)PULSEQLIB_ADC_ROLE_ECHO_CENTER
                        : (float)PULSEQLIB_ADC_ROLE_NON_CENTER;
                }
            }
        }
    }

    /* ===========================================================
     * Step 8: Composite RF group annotation
     *   Emit a group for every run of >= 2 RF events uninterrupted
     *   by an ADC event (WAIT events are transparent).
     * =========================================================== */
    {
        int n_groups = 0, g, group_start, rf_run, last_rf;

        /* Count groups */
        group_start = -1; rf_run = 0;
        for (i = 0; i < n_events; ++i) {
            int etype = out->events[i].type;
            if (etype == PULSEQLIB_SEQ_EVENT_RF) {
                if (group_start < 0) group_start = i;
                rf_run++;
            } else if (etype == PULSEQLIB_SEQ_EVENT_ADC) {
                if (rf_run >= 2) n_groups++;
                group_start = -1; rf_run = 0;
            }
            /* WAIT is transparent — does not break an RF group */
        }
        if (rf_run >= 2) n_groups++;

        if (n_groups > 0) {
            out->composite_rf_groups = (pulseqlib_composite_rf_group*)PULSEQLIB_ALLOC(
                (size_t)n_groups * sizeof(pulseqlib_composite_rf_group));
            if (!out->composite_rf_groups) {
                ret = PULSEQLIB_ERR_ALLOC_FAILED; goto cleanup;
            }

            g = 0; group_start = -1; rf_run = 0; last_rf = -1;
            for (i = 0; i <= n_events; ++i) {
                /* Sentinel: flush remaining group at loop end */
                int etype = (i < n_events)
                          ? out->events[i].type
                          : PULSEQLIB_SEQ_EVENT_ADC;
                if (etype == PULSEQLIB_SEQ_EVENT_RF) {
                    if (group_start < 0) group_start = i;
                    last_rf = i;
                    rf_run++;
                } else if (etype == PULSEQLIB_SEQ_EVENT_ADC) {
                    if (rf_run >= 2 && g < n_groups) {
                        out->composite_rf_groups[g].group_id        = g;
                        out->composite_rf_groups[g].first_event_idx = group_start;
                        out->composite_rf_groups[g].last_event_idx  = last_rf;
                        out->composite_rf_groups[g].num_pulses      = rf_run;
                        out->composite_rf_groups[g].eff_te_us       = 0.0f;
                        g++;
                    }
                    group_start = -1; rf_run = 0; last_rf = -1;
                }
            }
        }
        out->num_composite_rf_groups = n_groups;
    }

    /* Fall through to cleanup (success path) */

cleanup:
    if (blk_start)   PULSEQLIB_FREE(blk_start);
    if (adc_kzero)   PULSEQLIB_FREE(adc_kzero);
    if (rf_def_used) PULSEQLIB_FREE(rf_def_used);
    if (tmp_events)  PULSEQLIB_FREE(tmp_events);

    if (ret != PULSEQLIB_SUCCESS)
        pulseqlib_free_sequence_description(out);

    return ret;
}

/* ================================================================== */
/*  pulseqlib_get_sequence_parameters                                 */
/* ================================================================== */

int pulseqlib_get_sequence_parameters(
    pulseqlib_sequence_parameters* out,
    const pulseqlib_collection*    coll)
{
    const pulseqlib_sequence_descriptor* desc;
    const pulseqlib_tr_descriptor*       trd;
    int ss, i;
    float fa_max = 0.0f;

    if (!out || !coll) return PULSEQLIB_ERR_NULL_POINTER;
    memset(out, 0, sizeof(*out));

    out->min_te_us  = 1e30f;
    out->min_tr_us  = 1e30f;
    out->num_subseqs = coll->num_subsequences;

    for (ss = 0; ss < coll->num_subsequences; ++ss) {
        desc = &coll->descriptors[ss];
        trd  = &desc->tr_descriptor;

        /* TR bounds */
        if (trd->tr_duration_us > 0.0f) {
            if (trd->tr_duration_us < out->min_tr_us) out->min_tr_us = trd->tr_duration_us;
            if (trd->tr_duration_us > out->max_tr_us) out->max_tr_us = trd->tr_duration_us;
        }

        /* TE: earliest acquired ADC center in TR */
        for (i = 0; i < trd->tr_size; ++i) {
            int blk = trd->imaging_tr_start + i;
            if (blk < 0 || blk >= desc->num_blocks) continue;
            {
                const pulseqlib_block_table_element* bte = &desc->block_table[blk];
                float t_sum = 0.0f;
                int k;
                for (k = 0; k < i; ++k) {
                    int bk = trd->imaging_tr_start + k;
                    if (bk >= 0 && bk < desc->num_blocks)
                        t_sum += (float)desc->block_table[bk].duration_us;
                }
                if (bte->adc_id >= 0 && bte->adc_id < desc->adc_table_size
                    && !bte->pmc_flag && !bte->nav_flag) {
                    int adc_def_id = desc->adc_table[bte->adc_id].id;
                    if (adc_def_id >= 0 && adc_def_id < desc->num_unique_adcs) {
                        const pulseqlib_adc_definition* ad = &desc->adc_definitions[adc_def_id];
                        float te_est = t_sum + (float)ad->delay
                                     + 0.5f * (float)ad->num_samples
                                       * (float)ad->dwell_time * 1e-3f;
                        if (te_est < out->min_te_us) out->min_te_us = te_est;
                    }
                }
                /* Flip angle max */
                if (bte->rf_id >= 0 && bte->rf_id < desc->rf_table_size) {
                    int rf_def_id = desc->rf_table[bte->rf_id].id;
                    if (rf_def_id >= 0 && rf_def_id < desc->num_unique_rfs) {
                        float fa = desc->rf_definitions[rf_def_id].stats.flip_angle_deg;
                        if (fa > fa_max) fa_max = fa;
                    }
                }
            }
        }

        /* Total scan time approximation */
        out->total_scan_time_us += trd->tr_duration_us
                                  * (float)(trd->num_trs > 0 ? trd->num_trs : 1);
    }

    if (out->min_te_us > 1e29f) out->min_te_us = 0.0f;
    if (out->min_tr_us > 1e29f) out->min_tr_us = 0.0f;
    out->max_flip_angle_deg = fa_max;

    return PULSEQLIB_SUCCESS;
}
