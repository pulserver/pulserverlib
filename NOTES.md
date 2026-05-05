## 2026-05-05 — Fix validate ADC anchors and restore rf_shape_tuples
- `get_tr_waveforms` failed on multi-block repeated segments (e.g., mprage) because segment lookup used `find_segment_for_block_pos(n)` against `segment.start_block` (first-occurrence only), causing missing anchors on later occurrences.
- Correct fix in `pulseqlib_waveforms.c`: use `scan_table_seg_id[block_idx]` only when `block_order == NULL` (degenerate/single-pass path where block_idx maps to scan-table position).
- For ADC anchor `blk_in_seg`, use precomputed occurrence offsets (`scan_blk_in_occ`) only when `block_order == NULL`; keep fallback `n - seg->start_block` for non-degenerate average-expanded passes.
- Gotcha: in non-degenerate average-expanded passes, `block_order[n]` is a block-table index, not a scan-table position; indexing `scan_table_seg_id` with it breaks anchor lookup (seen in bssfp navg=3).
- `_get_sequence_description` now emits `rf_shape_tuples` by scanning segment blocks per unique RF, retrieving stats via `pulseqlib_get_rf_stats`, and waveform arrays via `pulseqlib_get_rf_magnitude/phase/time_us`.
- Wrapper implementation uses `pulseqlib_internal.h` for read-only descriptor access (`rf_definitions`, `segment_definitions`) and frees all C-allocated RF arrays with `PULSEQLIB_FREE`.
