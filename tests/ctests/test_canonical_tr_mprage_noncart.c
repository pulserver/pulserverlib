
/*
 * test_canonical_tr_mprage_noncart.c
 * Focused test for canonical TR extraction on noncartesian MPRAGE, rotext=0.
 *
 * For each userotext0 sequence:
 *   1. The canonical segment sequence must have exactly 1 element (degenerate
 *      topology: single main segment, both prep and cooldown degenerate).
 *   2. get_tr_gradient_waveforms must succeed for the single subsequence.
 */
#include "test_helpers.h"

static void run_canon_tr_case(const char *seq_name, const char *seq_file,
                               int num_averages)
{
    pulseqlib_collection *coll = NULL;
    pulseqlib_opts opts;
    int n_canon, rc;
    pulseqlib_tr_gradient_waveforms wf = PULSEQLIB_TR_GRADIENT_WAVEFORMS_INIT;

    gre_opts_init(&opts);
    rc = load_seq_with_averages(&coll, seq_file, &opts, num_averages);
    mu_assert(PULSEQLIB_SUCCEEDED(rc) && coll != NULL, "failed to load sequence");

    n_canon = pulseqlib_get_canonical_segment_sequence(coll, 0, NULL);
    fprintf(stderr, "[CANON_TR] %s: num_canonical_segments=%d\n", seq_name, n_canon);
    mu_assert_int_eq(1, n_canon);

    rc = pulseqlib_get_tr_gradient_waveforms(coll, 0, &wf, NULL);
    mu_assert(PULSEQLIB_SUCCEEDED(rc), "canonical TR waveform extraction failed");
    pulseqlib_tr_gradient_waveforms_free(&wf);

    pulseqlib_collection_free(coll);
}

MU_TEST(test_canon_tr_mprage_nc_1sl_1avg_userotext0)
{
    run_canon_tr_case("mprage_noncart_3d_1sl_1avg_userotext0",
                      "mprage_noncart_3d_1sl_1avg_userotext0.seq", 1);
}

MU_TEST(test_canon_tr_mprage_nc_1sl_3avg_userotext0)
{
    run_canon_tr_case("mprage_noncart_3d_1sl_3avg_userotext0",
                      "mprage_noncart_3d_1sl_3avg_userotext0.seq", 3);
}

MU_TEST(test_canon_tr_mprage_nc_3sl_1avg_userotext0)
{
    run_canon_tr_case("mprage_noncart_3d_3sl_1avg_userotext0",
                      "mprage_noncart_3d_3sl_1avg_userotext0.seq", 1);
}

MU_TEST(test_canon_tr_mprage_nc_3sl_3avg_userotext0)
{
    run_canon_tr_case("mprage_noncart_3d_3sl_3avg_userotext0",
                      "mprage_noncart_3d_3sl_3avg_userotext0.seq", 3);
}

MU_TEST_SUITE(suite_canon_tr_mprage_noncart)
{
    MU_RUN_TEST(test_canon_tr_mprage_nc_1sl_1avg_userotext0);
    MU_RUN_TEST(test_canon_tr_mprage_nc_1sl_3avg_userotext0);
    MU_RUN_TEST(test_canon_tr_mprage_nc_3sl_1avg_userotext0);
    MU_RUN_TEST(test_canon_tr_mprage_nc_3sl_3avg_userotext0);
}

int test_canonical_tr_mprage_noncart_main(void)
{
    MU_RUN_SUITE(suite_canon_tr_mprage_noncart);
    MU_REPORT();
    return MU_EXIT_CODE;
}
