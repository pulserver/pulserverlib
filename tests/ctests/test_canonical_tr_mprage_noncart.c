

#include "test_helpers.h"


// test_canonical_tr_mprage_noncart.c
// Minimal test for canonical TR extraction on noncartesian MPRAGE/userotext0

#include <stdio.h>
#include <assert.h>
#include <stdio.h>
#include <assert.h>




#include "test_helpers.h"

TEST_MAYBE_UNUSED int test_canonical_tr_mprage_noncart_main(void) {
    int failed = 0;
    static const char *seq_names[] = {
        "mprage_noncart_3d_1sl_1avg_userotext0",
        "mprage_noncart_3d_1sl_3avg_userotext0",
        "mprage_noncart_3d_3sl_1avg_userotext0",
        "mprage_noncart_3d_3sl_3avg_userotext0"
    };
    static const char *seq_files[] = {
        "mprage_noncart_3d_1sl_1avg_userotext0.seq",
        "mprage_noncart_3d_1sl_3avg_userotext0.seq",
        "mprage_noncart_3d_3sl_1avg_userotext0.seq",
        "mprage_noncart_3d_3sl_3avg_userotext0.seq"
    };
    for (int case_idx = 0; case_idx < 4; ++case_idx) {
        const char *seq_name = seq_names[case_idx];
        const char *seq_file = seq_files[case_idx];
        pulseqlib_collection *coll = NULL;
        pulseqlib_opts opts;
        default_opts_init(&opts);
        int rc = load_seq(&coll, seq_file, &opts);
        if (rc < 0 || !coll) {
            printf("[CANON_TR] %s: FAILED to load %s\n", seq_name, seq_file);
            ++failed;
            continue;
        }
        int n_canon = pulseqlib_get_canonical_segment_sequence(coll, 0, NULL);
        printf("[CANON_TR] %s: num_canonical_tr = %d\n", seq_name, n_canon);
        if (n_canon != 1) {
            printf("[CANON_TR] %s: WARNING: expected 1 canonical TR, got %d\n", seq_name, n_canon);
            ++failed;
        }
        printf("[CANON_TR] %s: extracting canonical TR\n", seq_name);
        pulseqlib_tr_gradient_waveforms wf = PULSEQLIB_TR_GRADIENT_WAVEFORMS_INIT;
        rc = pulseqlib_get_tr_gradient_waveforms(coll, 0, &wf, NULL);
        if (rc < 0) {
            printf("[CANON_TR] %s: canonical TR extraction FAILED\n", seq_name);
            ++failed;
        } else {
            printf("[CANON_TR] %s: canonical TR extraction OK\n", seq_name);
            pulseqlib_tr_gradient_waveforms_free(&wf);
        }
        pulseqlib_collection_free(coll);
    }
    return failed;
}
