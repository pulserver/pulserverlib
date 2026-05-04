/* pulseqlib_cache_seqdesc.c -- Section 5 cache writer for sequence description
 *
 * Implements:
 *   pulseqlib_write_sequence_description_cache()
 *
 * Section 5 serialization format (all values are 4 bytes, little-endian by
 * default — same convention as all other cache sections):
 *
 * [sequence parameters]
 *   min_te_us, min_tr_us, max_tr_us, max_flip_angle_deg, total_scan_time_us  (float x5)
 *   num_subseqs                                                               (int x1)
 *   reserved[3]                                                               (int x3)
 *
 * [per-subsequence block, for ss = 0 .. num_subseqs-1]
 *   subseq_idx          (int)
 *   tr_duration_us      (float)
 *   num_rows            (int)     — number of blocks in the pass
 *   [rows, for r = 0 .. num_rows-1]
 *     type              (int)     — PULSEQLIB_SEQ_EVENT_{OTHER,RF,ADC}
 *     timestamp_us      (float)   — pass-relative anchor time (us)
 *     params[6]         (float x6) — type-specific payload (zero-padded)
 *       RF:    params[0]=rf_def_id, [1]=rf_use, [2]=act_amplitude_hz,
 *              [3]=phase_offset_rad, [4]=freq_offset_hz, [5]=rf_shim_id
 *       ADC:   params[0]=adc_role, [1]=phase_offset_rad, [2..5]=0
 *       OTHER: params[0..5]=0
 */

#include <string.h>
#include <stdio.h>
#include <errno.h>

#include "pulseqlib_internal.h"
#include "pulseqlib_methods.h"

/* ------------------------------------------------------------------ */
/*  I/O helpers (duplicated from pulseqlib_trajectory.c pattern)       */
/* ------------------------------------------------------------------ */

static int sd_write4(FILE *f, const void *p, int count)
{
    return (int)fwrite(p, 4, (size_t)count, f) == count;
}

static int sd_read4(FILE *f, void *p, int count)
{
    return (int)fread(p, 4, (size_t)count, f) == count;
}

static void sd_swap4(void *p)
{
    unsigned char *b = (unsigned char *)p;
    unsigned char t;
    t = b[0];
    b[0] = b[3];
    b[3] = t;
    t = b[1];
    b[1] = b[2];
    b[2] = t;
}

static void sd_swap4_array(void *p, int count)
{
    int i;
    for (i = 0; i < count; ++i)
        sd_swap4((unsigned char *)p + (size_t)i * 4);
}

/* ------------------------------------------------------------------ */
/*  Path helper: .seq → .bin (same as in trajectory.c)                */
/* ------------------------------------------------------------------ */
static char *sd_make_cache_path(const char *seq_path)
{
    size_t len = strlen(seq_path);
    char *p;
    if (len < 4 || strcmp(seq_path + len - 4, ".seq") != 0)
        return NULL;
    p = (char *)PULSEQLIB_ALLOC(len + 1);
    if (!p)
        return NULL;
    memcpy(p, seq_path, len - 4);
    memcpy(p + len - 4, ".bin", 5);
    return p;
}

/* ------------------------------------------------------------------ */
/*  Write one subsequence's sequence description block                 */
/* ------------------------------------------------------------------ */
static int sd_write_subseq(FILE *f, const pulseqlib_sequence_description *sd)
{
    int i;

    if (!sd_write4(f, &sd->subseq_idx, 1))
        return 0;
    if (!sd_write4(f, &sd->tr_duration_us, 1))
        return 0;

    /* Compact row table: one row per pass block */
    if (!sd_write4(f, &sd->num_rows, 1))
        return 0;
    for (i = 0; i < sd->num_rows; ++i)
    {
        const pulseqlib_seq_event *row = &sd->rows[i];
        if (!sd_write4(f, &row->type, 1))
            return 0;
        if (!sd_write4(f, &row->timestamp_us, 1))
            return 0;
        if (!sd_write4(f, row->params, PULSEQLIB_SEQ_EVENT_PARAMS))
            return 0;
    }

    return 1;
}

/* ================================================================== */
/*  pulseqlib_write_sequence_description_cache                        */
/* ================================================================== */

#define SD_CACHE_ENDIAN_MARKER 0x01020304
#define SD_CACHE_SECTION_SEQDESC 5

int pulseqlib_write_sequence_description_cache(
    const pulseqlib_collection *coll,
    const char *seq_path)
{
    char *cache_path = NULL;
    FILE *f = NULL;
    int marker, num_sections;
    int version_major, version_minor, vendor, stored_size;
    int do_swap;
    long entries_pos, data_start, data_end, hdr_ns_pos;
    int i, found_idx;
    int entries_buf[18 * 3]; /* up to 18 sections x 3 ints */
    pulseqlib_sequence_description sd;
    pulseqlib_sequence_parameters sp;
    int ret = PULSEQLIB_SUCCESS;

    if (!coll || !seq_path)
        return PULSEQLIB_ERR_NULL_POINTER;

    cache_path = sd_make_cache_path(seq_path);
    if (!cache_path)
        return PULSEQLIB_ERR_ALLOC_FAILED;

    f = fopen(cache_path, "r+b");
    if (!f)
    {
        ret = PULSEQLIB_ERR_FILE_READ_FAILED;
        goto done;
    }

    /* Read file header */
    if (!sd_read4(f, &marker, 1))
    {
        ret = PULSEQLIB_ERR_FILE_READ_FAILED;
        goto done;
    }
    do_swap = 0;
    if (marker != SD_CACHE_ENDIAN_MARKER)
    {
        sd_swap4(&marker);
        if (marker != SD_CACHE_ENDIAN_MARKER)
        {
            ret = PULSEQLIB_ERR_FILE_READ_FAILED;
            goto done;
        }
        do_swap = 1;
    }
    if (!sd_read4(f, &version_major, 1))
    {
        ret = PULSEQLIB_ERR_FILE_READ_FAILED;
        goto done;
    }
    if (!sd_read4(f, &version_minor, 1))
    {
        ret = PULSEQLIB_ERR_FILE_READ_FAILED;
        goto done;
    }
    if (!sd_read4(f, &vendor, 1))
    {
        ret = PULSEQLIB_ERR_FILE_READ_FAILED;
        goto done;
    }
    if (!sd_read4(f, &stored_size, 1))
    {
        ret = PULSEQLIB_ERR_FILE_READ_FAILED;
        goto done;
    }
    hdr_ns_pos = ftell(f);
    if (hdr_ns_pos < 0)
    {
        ret = PULSEQLIB_ERR_FILE_READ_FAILED;
        goto done;
    }
    if (!sd_read4(f, &num_sections, 1))
    {
        ret = PULSEQLIB_ERR_FILE_READ_FAILED;
        goto done;
    }
    if (do_swap)
    {
        sd_swap4(&version_major);
        sd_swap4(&version_minor);
        sd_swap4(&vendor);
        sd_swap4(&stored_size);
        sd_swap4(&num_sections);
    }
    if (num_sections <= 0 || num_sections > 17)
    {
        ret = PULSEQLIB_ERR_FILE_READ_FAILED;
        goto done;
    }

    entries_pos = ftell(f);
    if (entries_pos < 0)
    {
        ret = PULSEQLIB_ERR_FILE_READ_FAILED;
        goto done;
    }

    /* Read existing section entries */
    for (i = 0; i < num_sections; ++i)
    {
        if (!sd_read4(f, &entries_buf[i * 3], 3))
        {
            ret = PULSEQLIB_ERR_FILE_READ_FAILED;
            goto done;
        }
        if (do_swap)
            sd_swap4_array(&entries_buf[i * 3], 3);
    }

    /* Find or allocate slot for section 5 */
    found_idx = -1;
    for (i = 0; i < num_sections; ++i)
    {
        if (entries_buf[i * 3] == SD_CACHE_SECTION_SEQDESC)
        {
            found_idx = i;
            break;
        }
    }
    if (found_idx < 0)
    {
        found_idx = num_sections;
        entries_buf[found_idx * 3] = SD_CACHE_SECTION_SEQDESC;
        num_sections++;
    }

    /* Seek to end and start writing */
    if (fseek(f, 0, SEEK_END) != 0)
    {
        ret = PULSEQLIB_ERR_FILE_READ_FAILED;
        goto done;
    }
    data_start = ftell(f);
    if (data_start < 0)
    {
        ret = PULSEQLIB_ERR_FILE_READ_FAILED;
        goto done;
    }

    /* Sequence parameters */
    ret = pulseqlib_get_sequence_parameters(&sp, coll);
    if (ret != PULSEQLIB_SUCCESS)
        goto done;

    if (!sd_write4(f, &sp.min_te_us, 1) ||
        !sd_write4(f, &sp.min_tr_us, 1) ||
        !sd_write4(f, &sp.max_tr_us, 1) ||
        !sd_write4(f, &sp.max_flip_angle_deg, 1) ||
        !sd_write4(f, &sp.total_scan_time_us, 1) ||
        !sd_write4(f, &sp.num_subseqs, 1) ||
        !sd_write4(f, sp.reserved, 3))
    {
        ret = PULSEQLIB_ERR_FILE_READ_FAILED;
        goto done;
    }

    /* Per-subsequence blocks */
    for (i = 0; i < coll->num_subsequences; ++i)
    {
        memset(&sd, 0, sizeof(sd));
        ret = pulseqlib_get_sequence_description(&sd, coll, i);
        if (ret != PULSEQLIB_SUCCESS)
        {
            pulseqlib_free_sequence_description(&sd);
            goto done;
        }
        if (!sd_write_subseq(f, &sd))
        {
            pulseqlib_free_sequence_description(&sd);
            ret = PULSEQLIB_ERR_FILE_READ_FAILED;
            goto done;
        }
        pulseqlib_free_sequence_description(&sd);
    }

    data_end = ftell(f);
    if (data_end < 0)
    {
        ret = PULSEQLIB_ERR_FILE_READ_FAILED;
        goto done;
    }

    entries_buf[found_idx * 3 + 1] = (int)data_start;
    entries_buf[found_idx * 3 + 2] = (int)(data_end - data_start);

    /* Patch num_sections in header */
    if (fseek(f, hdr_ns_pos, SEEK_SET) != 0)
    {
        ret = PULSEQLIB_ERR_FILE_READ_FAILED;
        goto done;
    }
    if (!sd_write4(f, &num_sections, 1))
    {
        ret = PULSEQLIB_ERR_FILE_READ_FAILED;
        goto done;
    }

    /* Rewrite all section entries */
    if (fseek(f, entries_pos, SEEK_SET) != 0)
    {
        ret = PULSEQLIB_ERR_FILE_READ_FAILED;
        goto done;
    }
    for (i = 0; i < num_sections; ++i)
    {
        if (!sd_write4(f, &entries_buf[i * 3], 3))
        {
            ret = PULSEQLIB_ERR_FILE_READ_FAILED;
            goto done;
        }
    }

done:
    if (f)
        fclose(f);
    if (cache_path)
        PULSEQLIB_FREE(cache_path);
    return ret;
}
