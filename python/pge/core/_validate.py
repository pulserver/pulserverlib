"""Validation: compare pulserver waveforms against a reference source."""

__all__ = ['validate']

from pathlib import Path

# Fallback stub for _validate_tr if not imported
def _validate_tr(seq, sidx, tidx):
    raise NotImplementedError("_validate_tr must be implemented or imported.")

import numpy as np

from ._sequence import SequenceCollection
from ._waveforms import get_tr_waveforms

# ── helpers ──────────────────────────────────────────────────────────


def _interp_to_ref(t_ref, a_ref, t_test, a_test):
    """Interpolate *test* waveform onto *ref* time base, return aligned pair."""
    t_ref = np.asarray(t_ref).ravel()
    a_ref = np.asarray(a_ref).ravel()
    t_test = np.asarray(t_test).ravel()
    a_test = np.asarray(a_test).ravel()

    if np.iscomplexobj(t_ref):
        t_ref = np.real(t_ref)
    if np.iscomplexobj(t_test):
        t_test = np.real(t_test)

    if t_ref.size == 0 or t_test.size == 0:
        return np.empty(0), np.empty(0)

    n_ref = min(t_ref.size, a_ref.size)
    n_test = min(t_test.size, a_test.size)
    if n_ref == 0 or n_test == 0:
        return np.empty(0), np.empty(0)

    t_ref = t_ref[:n_ref]
    a_ref = a_ref[:n_ref]
    t_test = t_test[:n_test]
    a_test = a_test[:n_test]

    # np.interp requires an ordered x-grid.
    order = np.argsort(t_test)
    t_test = t_test[order]
    a_test = a_test[order]

    # Remove duplicate x values so interpolation remains well-defined.
    t_test, uniq_idx = np.unique(t_test, return_index=True)
    a_test = a_test[uniq_idx]

    a_interp = np.interp(t_ref, t_test, a_test, left=0.0, right=0.0)
    return a_ref, a_interp


def _rms_error(ref, test):
    """Percentage RMS error relative to ref (0 if ref is silent)."""
    norm = np.sqrt(np.mean(ref**2))
    if norm < 1e-30:
        return 0.0
    return 100.0 * np.sqrt(np.mean((ref - test) ** 2)) / norm


def _is_non_degenerate(tr_info):
    """True for bSSFP-like layouts where each ACTUAL-mode TR spans a full pass.

    In that case the C library returns one waveform per *pass* (slice), with
    imaging blocks replicated ``num_averages`` times inside each pass.
    """
    has_nd_prep = (not tr_info['degenerate_prep'] and tr_info['num_prep_blocks'] > 0)
    has_nd_cool = (not tr_info['degenerate_cooldown'] and tr_info['num_cooldown_blocks'] > 0)
    return has_nd_prep or has_nd_cool


def _total_addressable_trs(tr_info, num_averages=1):
    """Total TRs addressable in ACTUAL mode (matches C-side logic).

    Non-degenerate (bSSFP): ``num_passes`` — each pass is a single TR that
    embeds all averages.

    Degenerate (GRE-like): prep/cooldown play once; imaging repeats
    ``num_averages`` times, giving ``num_prep_trs + num_averages * num_trs
    + num_cooldown_trs``.
    """
    if _is_non_degenerate(tr_info):
        return tr_info['num_passes']
    return (tr_info['num_prep_trs']
            + num_averages * tr_info['num_trs']
            + tr_info['num_cooldown_trs'])


def _canonical_pypulseq_block(tr_idx, tr_info, num_averages=1):
    """Map a validate-loop *tr_idx* to the 0-based block offset used to
    locate the corresponding time window in the pypulseq Sequence object.

    Non-degenerate: tr_idx selects a pass; block offset = pass_idx * pass_len.
    Degenerate: imaging indices wrap modulo ``num_trs`` to canonical blocks.
    """
    if _is_non_degenerate(tr_info):
        pass_len = (tr_info['num_prep_blocks']
                    + tr_info['num_trs'] * tr_info['tr_size']
                    + tr_info['num_cooldown_blocks'])
        return tr_idx * pass_len
    num_prep = tr_info['num_prep_trs']
    num_trs = tr_info['num_trs']
    if tr_idx < num_prep:
        canonical_tr = tr_idx
    elif tr_idx < num_prep + num_averages * num_trs:
        canonical_tr = num_prep + (tr_idx - num_prep) % num_trs
    else:
        cool_idx = tr_idx - num_prep - num_averages * num_trs
        canonical_tr = num_prep + num_trs + cool_idx
    return canonical_tr * tr_info['tr_size']


def _abs_tr_start_s(
    src_seq,
    num_prep_blocks: int,
    tr_idx: int,
    tr_duration_us: float,
    imaging_tr_start: int | None = None,
) -> float:
    """Compute the absolute start time (in seconds) of a TR.

    Sums block durations over the preparation region and adds
    ``tr_idx * tr_duration``.
    """
    prep_s = 0.0
    bd = src_seq.block_durations

    # Prefer explicit imaging-tr anchor from C when available.
    start_blocks = num_prep_blocks if imaging_tr_start is None else int(imaging_tr_start)
    if start_blocks < 0:
        start_blocks = 0

    for k in range(1, start_blocks + 1):
        prep_s += bd.get(k, 0.0) if isinstance(bd, dict) else bd[k - 1]
    return prep_s + tr_idx * (tr_duration_us * 1e-6)


def _abs_block_start_s(src_seq, block_start: int) -> float:
    """Compute absolute time in seconds for a 0-based block index."""
    if block_start <= 0:
        return 0.0

    prep_s = 0.0
    bd = src_seq.block_durations
    for k in range(1, block_start + 1):
        prep_s += bd.get(k, 0.0) if isinstance(bd, dict) else bd[k - 1]
    return prep_s


def _pypulseq_reference_window(seq, sequence_idx, abs_t0_s, window_duration_us):
    """Extract reference waveforms from pypulseq for an exact time window."""
    src_seq = seq._seqs[sequence_idx]
    gamma = seq.system.gamma

    # Expand window slightly at the start to ensure first block is included
    epsilon = 1e-9  # 1 ns
    abs_t1_s = abs_t0_s + window_duration_us * 1e-6
    result = src_seq.waveforms_and_times(append_RF=True, time_range=[abs_t0_s - epsilon, abs_t1_s])
    channels = result[0]

    hz_to_mT_per_m = 1.0 / (gamma * 1e-3)

    def _clip_window(time_us, amp):
        t = np.asarray(time_us)
        a = np.asarray(amp)
        if t.size == 0:
            # Pad with zeros for the full window
            return np.array([0.0, window_duration_us]), np.zeros(2)
        mask = (t >= -1e-6) & (t <= window_duration_us + 1e-6)
        t_clipped = t[mask]
        a_clipped = a[mask]
        # Pad start if needed
        if len(t_clipped) == 0 or t_clipped[0] > 0:
            t_clipped = np.insert(t_clipped, 0, 0.0)
            a_clipped = np.insert(a_clipped, 0, 0.0)
        # Pad end if needed
        if t_clipped[-1] < window_duration_us:
            t_clipped = np.append(t_clipped, window_duration_us)
            a_clipped = np.append(a_clipped, 0.0)
        return t_clipped, a_clipped

    ref: dict[str, tuple] = {}
    for ch_idx, ch_name in enumerate(('gx', 'gy', 'gz')):
        if ch_idx < len(channels):
            arr = channels[ch_idx]
            if arr.shape[1] > 0:
                t_rel = np.real((arr[0] - abs_t0_s) * 1e6)
                a_rel = np.real(arr[1]) * hz_to_mT_per_m
                ref[ch_name] = _clip_window(t_rel, a_rel)
                continue
        ref[ch_name] = (np.empty(0), np.empty(0))

    # RF magnitude and phase
    if len(channels) > 3:
        arr = channels[3]
        if arr.shape[1] > 0:
            t_rel = np.real((arr[0] - abs_t0_s) * 1e6)
            hz_to_uT = 1e6 / gamma
            rf_mag = np.abs(arr[1]) * hz_to_uT
            rf_phase = np.angle(arr[1])
            ref['rf_mag'] = _clip_window(t_rel, rf_mag)
            ref['rf_phase'] = _clip_window(t_rel, rf_phase)

            # Multichannel RF detection: in pTx pypulseq the channels are
            # concatenated so the first time value repeats once per channel.
            nch = 1
            n_total = arr.shape[1]
            if n_total > 1:
                n_repeat = int(np.sum(arr[0] == arr[0, 0]))
                if n_repeat > 1 and n_total % n_repeat == 0:
                    nch = n_repeat
            if nch > 1:
                npts_ch = n_total // nch
                rf_mag_ch_list = []
                rf_phase_ch_list = []
                for c in range(nch):
                    sl = slice(c * npts_ch, (c + 1) * npts_ch)
                    t_c = np.real((arr[0, sl] - abs_t0_s) * 1e6)
                    sig_c = arr[1, sl]
                    rf_mag_ch_list.append(_clip_window(t_c, np.abs(sig_c) * hz_to_uT))
                    rf_phase_ch_list.append(_clip_window(t_c, np.angle(sig_c)))
                ref['rf_mag_channels'] = rf_mag_ch_list
                ref['rf_phase_channels'] = rf_phase_ch_list
            # RF shim detection (pTx calibration mode): TODO when supported.
        else:
            ref['rf_mag'] = (np.empty(0), np.empty(0))
            ref['rf_phase'] = (np.empty(0), np.empty(0))
    else:
        ref['rf_mag'] = (np.empty(0), np.empty(0))
        ref['rf_phase'] = (np.empty(0), np.empty(0))

    # ADC per-event data from result[3] (t_adc) and result[4] (fp_adc)
    t_adc = np.asarray(result[3])          # absolute sample times (s)
    fp_adc = np.asarray(result[4])         # (num_events, 2): [freq_offset, phase_offset]
    ref['adc_events'] = []
    ref['adc_phase'] = (np.empty(0), np.empty(0))
    if t_adc.size > 0 and fp_adc.size > 0:
        fp_adc = np.atleast_2d(fp_adc)
        num_events = fp_adc.shape[0]
        t_rel_adc = (t_adc - abs_t0_s) * 1e6          # µs, TR-relative
        phase_all = np.zeros_like(t_rel_adc)
        if num_events == 1:
            ref['adc_events'] = [
                (float(t_rel_adc[0]), float(fp_adc[0, 0]), float(fp_adc[0, 1]))
            ]
            phase_all = fp_adc[0, 1] + 2 * np.pi * fp_adc[0, 0] * (t_rel_adc * 1e-6)
            ref['adc_phase'] = (t_rel_adc, phase_all)
        else:
            # Split samples across events via gap detection
            dt = np.diff(t_adc)
            median_dt = np.median(dt) if len(dt) > 0 else 1e-6
            boundaries = np.where(dt > 3 * median_dt)[0] + 1
            ev_starts = np.concatenate([[0], boundaries])
            ev_ends = np.concatenate([boundaries, [len(t_adc)]])

            phase_t_parts = []
            phase_a_parts = []
            for ev_idx in range(min(len(ev_starts), num_events)):
                onset_us = float(t_rel_adc[ev_starts[ev_idx]])
                ref['adc_events'].append(
                    (onset_us, float(fp_adc[ev_idx, 0]), float(fp_adc[ev_idx, 1]))
                )
                sl = slice(ev_starts[ev_idx], ev_ends[ev_idx])
                t_ev = t_rel_adc[sl]
                a_ev = (fp_adc[ev_idx, 1]
                        + 2 * np.pi * fp_adc[ev_idx, 0]
                        * (t_ev * 1e-6))
                phase_all[sl] = a_ev

                if t_ev.size > 0:
                    phase_t_parts.append(t_ev)
                    phase_a_parts.append(a_ev)

                # Represent non-ADC gaps explicitly as zero phase so plotting
                # does not linearly connect adjacent ADC events across gaps.
                if ev_idx + 1 < min(len(ev_starts), num_events):
                    t_gap_start = float(t_rel_adc[ev_ends[ev_idx] - 1])
                    t_gap_end = float(t_rel_adc[ev_starts[ev_idx + 1]])
                    if t_gap_end > t_gap_start:
                        phase_t_parts.append(np.array([t_gap_start, t_gap_end], dtype=float))
                        phase_a_parts.append(np.array([0.0, 0.0], dtype=float))

            if phase_t_parts:
                ref['adc_phase'] = (np.concatenate(phase_t_parts), np.concatenate(phase_a_parts))
            else:
                ref['adc_phase'] = (t_rel_adc, phase_all)

    return ref


# ── reference extraction ─────────────────────────────────────────────


def _concat_refs(segments):
    """Concatenate reference dicts from multiple (offset_us, ref) segments.

    Each *ref* is a dict returned by :func:`_pypulseq_reference_window`.
    Times in each ref are relative to that segment's start; *offset_us*
    shifts them into the combined timeline.
    """
    result = {}

    for key in ('gx', 'gy', 'gz', 'rf_mag', 'rf_phase', 'adc_phase'):
        parts_t, parts_a = [], []
        prev_t_end = None
        for off, ref in segments:
            if key in ref:
                t, a = ref[key]
                t = np.asarray(t, dtype=float)
                a = np.asarray(a, dtype=float)
                if len(t) > 0:
                    t_shift = t + off
                    # For ADC phase specifically, represent gaps between
                    # concatenated windows as zero phase to avoid spurious
                    # linear ramps across regions without ADC activity.
                    if key == 'adc_phase' and prev_t_end is not None:
                        t_start = float(t_shift[0])
                        if t_start > prev_t_end:
                            parts_t.append(np.array([prev_t_end, t_start], dtype=float))
                            parts_a.append(np.array([0.0, 0.0], dtype=float))

                    parts_t.append(t_shift)
                    parts_a.append(a)
                    prev_t_end = float(t_shift[-1])
        if parts_t:
            result[key] = (np.concatenate(parts_t), np.concatenate(parts_a))
        else:
            result[key] = (np.empty(0), np.empty(0))

    # Multichannel RF
    for key in ('rf_mag_channels', 'rf_phase_channels'):
        nch = 0
        for _, ref in segments:
            if key in ref:
                nch = max(nch, len(ref[key]))
        if nch == 0:
            continue
        ch_list = []
        for c in range(nch):
            parts_t, parts_a = [], []
            for off, ref in segments:
                if key in ref and c < len(ref[key]):
                    t, a = ref[key][c]
                    t = np.asarray(t)
                    a = np.asarray(a)
                    if len(t) > 0:
                        parts_t.append(t + off)
                        parts_a.append(a)
            if parts_t:
                ch_list.append((np.concatenate(parts_t), np.concatenate(parts_a)))
            else:
                ch_list.append((np.empty(0), np.empty(0)))
        result[key] = ch_list

    # ADC events
    adc_events = []
    for off, ref in segments:
        if 'adc_events' in ref:
            for onset, freq, phase in ref['adc_events']:
                adc_events.append((onset + off, freq, phase))
    result['adc_events'] = adc_events

    return result


def _pypulseq_reference(seq, sequence_idx, tr_idx, tr_info, num_averages=1,
                        duration_us=None):
    """Extract reference waveforms from pypulseq for a single TR.

    Parameters
    ----------
    duration_us : float | None
        Override window duration.  When *None*, uses ``tr_duration_us`` from
        *tr_info* for degenerate TRs, or the full pass duration (computed from
        block_durations) for non-degenerate.

    Returns
    -------
    dict
        Mapping of channel name to ``(time_us, amplitude)`` tuples.
        Times are relative to TR start; amplitudes in mT/m (grad) or
        µT (RF magnitude).
    """
    blk_offset = _canonical_pypulseq_block(tr_idx, tr_info, num_averages)
    t0 = _abs_block_start_s(
        seq._seqs[sequence_idx],
        blk_offset,
    )

    # Non-degenerate (bSSFP-like): the C library returns a full pass with
    # prep + eff_num_averages*imaging + cooldown.  Build a tiled reference
    # by extracting prep / imaging / cooldown separately and repeating the
    # imaging portion.  This applies even for num_averages <= 1 so the
    # overlay always spans the complete pass.
    if _is_non_degenerate(tr_info):
        eff_navg = max(num_averages, 1)
        src_seq = seq._seqs[sequence_idx]
        n_prep = tr_info['num_prep_blocks']
        n_img = tr_info['num_trs'] * tr_info['tr_size']
        n_cool = tr_info['num_cooldown_blocks']

        t_img_s = _abs_block_start_s(src_seq, blk_offset + n_prep)
        t_cool_s = _abs_block_start_s(src_seq, blk_offset + n_prep + n_img)
        t_end_s = _abs_block_start_s(src_seq, blk_offset + n_prep + n_img + n_cool)

        prep_dur_us = (t_img_s - t0) * 1e6
        img_dur_us = (t_cool_s - t_img_s) * 1e6
        cool_dur_us = (t_end_s - t_cool_s) * 1e6

        segments = []
        if prep_dur_us > 0:
            segments.append((0.0, _pypulseq_reference_window(
                seq, sequence_idx, t0, prep_dur_us)))
        ref_img = _pypulseq_reference_window(
            seq, sequence_idx, t_img_s, img_dur_us)
        for k in range(eff_navg):
            segments.append((prep_dur_us + k * img_dur_us, ref_img))
        if cool_dur_us > 0:
            segments.append((
                prep_dur_us + eff_navg * img_dur_us,
                _pypulseq_reference_window(seq, sequence_idx, t_cool_s, cool_dur_us),
            ))
        return _concat_refs(segments)

    if duration_us is None:
        duration_us = tr_info['tr_duration_us']
    return _pypulseq_reference_window(seq, sequence_idx, t0, duration_us)


def _xml_reference(xml_path):
    """Extract reference waveforms from an XML file.

    Returns
    -------
    dict
        Same format as :func:`_pypulseq_reference`.
    """
    import xml.etree.ElementTree as ET

    tree = ET.parse(str(xml_path))
    root = tree.getroot()

    g_per_cm_to_mT_per_m = 10.0
    g_to_uT = 100.0

    ref: dict[str, tuple] = {}
    for ch_name in ('gx', 'gy', 'gz'):
        elem = root.find(f'.//{ch_name}')
        if elem is not None and elem.text:
            data = np.fromstring(elem.text, sep=' ')
            if data.size >= 2:
                n = data.size // 2
                ref[ch_name] = (data[:n], data[n:] * g_per_cm_to_mT_per_m)
                continue
        ref[ch_name] = (np.empty(0), np.empty(0))

    rf_elem = root.find('.//rf')
    if rf_elem is not None and rf_elem.text:
        data = np.fromstring(rf_elem.text, sep=' ')
        if data.size >= 2:
            n = data.size // 2
            ref['rf_mag'] = (data[:n], np.abs(data[n:]) * g_to_uT)
        else:
            ref['rf_mag'] = (np.empty(0), np.empty(0))
    else:
        ref['rf_mag'] = (np.empty(0), np.empty(0))

    return ref


# ── public entry point ───────────────────────────────────────────────

def validate(
    seq: SequenceCollection,
    *,
    subsequence_idx: int | None = None,
    tr_instance: int | None = None,
    xml_path: str | Path | None = None,
    do_plot: bool = False,
    grad_atol: float | None = None,
    rf_rms_percent: float = 10.0,
) -> None:
    """
    Validate subsequences and TRs against pypulseq or an XML reference.

    Behaviour
    ---------
    - ``do_plot=False``:
      - ``subsequence_idx=None`` validates all subsequences.
      - ``tr_instance=None`` validates all TRs for each selected subsequence.
    - ``do_plot=True``:
      - A concrete (subsequence, TR) target is required.
      - ``subsequence_idx=None`` auto-resolves to 0 only if there is exactly
        one subsequence; otherwise raises ``ValueError``.
      - ``tr_instance=None`` auto-resolves to 0 only if the selected
        subsequence has exactly one addressable TR; otherwise raises
        ``ValueError``.
    """
    from ._extension._pulseqlib_wrapper import _find_tr

    num_subseq = seq.num_sequences
    num_averages = seq.num_averages

    # ── resolve validation scope ──────────────────────────────────────
    if do_plot:
        targeted = True
        if subsequence_idx is None:
            if num_subseq == 1:
                subsequence_idx = 0
            else:
                raise ValueError(
                    'subsequence_idx must be specified when do_plot=True and '
                    f'the sequence has {num_subseq} subsequences.'
                )
        if subsequence_idx < 0 or subsequence_idx >= num_subseq:
            raise ValueError(
                f'subsequence_idx={subsequence_idx} out of range '
                f'(valid: 0..{num_subseq - 1})'
            )

        tr_info_target = _find_tr(seq._cseq, subsequence_idx=subsequence_idx)
        num_total = _total_addressable_trs(tr_info_target, num_averages)
        if tr_instance is None:
            if num_total == 1:
                tr_instance = 0
            else:
                raise ValueError(
                    'tr_instance must be specified when do_plot=True and '
                    f'subsequence {subsequence_idx} has {num_total} TRs.'
                )
        elif tr_instance < 0:
            tr_instance = num_total + tr_instance
        if tr_instance < 0 or tr_instance >= num_total:
            raise ValueError(
                f'tr_instance={tr_instance} out of range for {num_total} TRs '
                f'(valid: 0..{num_total - 1} or -{num_total}..-1)'
            )

        ss_range = [subsequence_idx]
        tr_range = {subsequence_idx: [tr_instance]}
    else:
        targeted = False
        if subsequence_idx is None:
            ss_range = list(range(num_subseq))
        else:
            if subsequence_idx < 0 or subsequence_idx >= num_subseq:
                raise ValueError(
                    f'subsequence_idx={subsequence_idx} out of range '
                    f'(valid: 0..{num_subseq - 1})'
                )
            ss_range = [subsequence_idx]

        tr_range = {}
        for ss in ss_range:
            tr_info = _find_tr(seq._cseq, subsequence_idx=ss)
            num_total = _total_addressable_trs(tr_info, num_averages)
            if tr_instance is None:
                tr_range[ss] = list(range(num_total))
            else:
                tr_idx = tr_instance
                if tr_idx < 0:
                    tr_idx = num_total + tr_idx
                if tr_idx < 0 or tr_idx >= num_total:
                    raise ValueError(
                        f'tr_instance={tr_instance} out of range for subsequence {ss} '
                        f'with {num_total} TRs (valid: 0..{num_total - 1} or '
                        f'-{num_total}..-1)'
                    )
                tr_range[ss] = [tr_idx]

    # ── common derived settings ───────────────────────────────────────
    sys = seq.system
    if grad_atol is None:
        grad_atol = 3.0 * sys.max_slew * sys.grad_raster_time * 1e3
    _max_grad_plot = sys.max_grad / sys.gamma * 1e3

    # ── validation loop ───────────────────────────────────────────────
    fail = False
    fail_msg = []
    fail_subseq_idx = None
    fail_tr_idx = None

    for ss_idx in ss_range:
        tr_info = _find_tr(seq._cseq, subsequence_idx=ss_idx)
        num_trs_ss = _total_addressable_trs(tr_info, num_averages)
        multi_tr = num_trs_ss > 1
        for tr_idx in tr_range[ss_idx]:
            messages = []
            wf = get_tr_waveforms(
                seq,
                subsequence_idx=ss_idx,
                amplitude_mode='actual',
                tr_index=tr_idx,
                num_averages=num_averages,
            )
            # Reference extraction
            if xml_path is None:
                ref = _pypulseq_reference(
                    seq, ss_idx, tr_idx, tr_info, num_averages,
                    duration_us=wf.total_duration_us,
                )
            else:
                ref = _xml_reference(xml_path)
            # Sort every (t, a) entry so time axes are monotone.
            for _k in list(ref.keys()):
                if isinstance(ref[_k], tuple):
                    _t, _a = np.asarray(ref[_k][0]), np.asarray(ref[_k][1])
                    if _t.ndim == 1 and len(_t) > 1:
                        _ord = np.argsort(_t)
                        ref[_k] = (_t[_ord], _a[_ord])
            ref_t, ref_a = ref['rf_mag']

            # Compare gradients
            for ch in ('gx', 'gy', 'gz'):
                ref_t_ch, ref_a_ch = ref[ch]
                test_ch = getattr(wf, ch)
                r, t = _interp_to_ref(ref_t_ch, ref_a_ch, test_ch.time_us, test_ch.amplitude)
                err = float(np.max(np.abs(r - t))) if len(r) > 0 else 0.0
                if err > grad_atol:
                    pfx = f'TR {tr_idx}: ' if multi_tr else ''
                    messages.append(
                        f'{pfx}{ch} mismatch: max diff {err:.4f} mT/m '
                        f'(tol {grad_atol:.4f} mT/m)'
                    )

            # RF envelope (magnitude)
            ref_t_real = np.real(ref_t)
            ref_a_real = np.real(ref_a)
            test_t = wf.rf_mag.time_us
            test_a = wf.rf_mag.amplitude
            t_min = max(ref_t_real[0], test_t[0]) if len(ref_t_real) and len(test_t) else 0.0
            t_max = min(ref_t_real[-1], test_t[-1]) if len(ref_t_real) and len(test_t) else 0.0
            if t_max <= t_min:
                messages.append(f"RF time overlap is empty: ref=[{ref_t_real[0] if len(ref_t_real) else 'NA'}, {ref_t_real[-1] if len(ref_t_real) else 'NA'}], test=[{test_t[0] if len(test_t) else 'NA'}, {test_t[-1] if len(test_t) else 'NA'}]")
            else:
                t_common = np.arange(np.ceil(t_min), np.floor(t_max) + 1)
                if len(t_common) < 3:
                    messages.append(f"RF overlap region too small for robust comparison (len={len(t_common)})")
                else:
                    t_common = t_common[1:-1]
                    if len(t_common) == 0:
                        messages.append("RF overlap region empty after ignoring endpoints.")
                    else:
                        ref_interp = np.abs(np.interp(t_common, ref_t_real, ref_a_real, left=0.0, right=0.0))
                        test_interp = np.abs(np.interp(t_common, test_t, test_a, left=0.0, right=0.0))
                        rf_err = _rms_error(ref_interp, test_interp)
                        if rf_err > rf_rms_percent:
                            pfx = f'TR {tr_idx}: ' if multi_tr else ''
                            messages.append(
                                f'{pfx}RF mismatch: {rf_err:.1f}% RMS (tol {rf_rms_percent:.1f}%)'
                            )

            # RF phase validation — use complex RF comparison so that
            # sign-convention differences (π phase flip at RF edges) are
            # handled naturally instead of inflating the phase RMS.
            test_phase_t = np.real(wf.rf_phase.time_us)
            test_phase_a = np.real(np.copy(wf.rf_phase.amplitude))
            if 'rf_phase' in ref and len(test_phase_t) > 0:
                ref_phase_t, ref_phase_a = ref['rf_phase']
                ref_phase_t = np.real(ref_phase_t)
                ref_phase_a = np.real(ref_phase_a)
                t_min_p = max(ref_phase_t[0], test_phase_t[0]) if len(ref_phase_t) and len(test_phase_t) else 0.0
                t_max_p = min(ref_phase_t[-1], test_phase_t[-1]) if len(ref_phase_t) and len(test_phase_t) else 0.0
                if t_max_p > t_min_p:
                    t_common_p = np.arange(np.ceil(t_min_p), np.floor(t_max_p) + 1)
                    t_common_p = t_common_p[1:-1] if len(t_common_p) > 2 else t_common_p
                    if len(t_common_p) > 0:
                        ref_interp_p = np.interp(t_common_p, ref_phase_t, ref_phase_a, left=0.0, right=0.0)
                        test_interp_p = np.interp(t_common_p, test_phase_t, test_phase_a, left=0.0, right=0.0)
                        ref_mag_t, ref_mag_a = ref['rf_mag']
                        ref_mag_interp = np.interp(t_common_p, np.real(ref_mag_t), ref_mag_a, left=0.0, right=0.0)
                        test_mag_interp = np.interp(t_common_p, np.real(wf.rf_mag.time_us), wf.rf_mag.amplitude, left=0.0, right=0.0)
                        mag_peak = max(ref_mag_interp.max(), test_mag_interp.max(), 1e-30)
                        active = (np.maximum(ref_mag_interp, test_mag_interp) > 0.05 * mag_peak)
                        if active.sum() > 0:
                            ref_c = ref_mag_interp[active] * np.exp(1j * ref_interp_p[active])
                            test_c = test_mag_interp[active] * np.exp(1j * test_interp_p[active])
                            phase_err = float(np.sqrt(np.mean(np.abs(ref_c - test_c) ** 2))) / mag_peak
                            phase_tol = 0.1  # fractional complex RF tolerance
                            if phase_err > phase_tol:
                                pfx = f'TR {tr_idx}: ' if multi_tr else ''
                                messages.append(
                                    f'{pfx}RF phase mismatch: {phase_err:.3f} (tol {phase_tol:.3f})'
                                )

            # ADC phase validation — event-level comparison of freq/phase
            # offsets rather than interpolated waveform comparison (avoids
            # cross-gap interpolation artifacts for multi-event passes).
            if 'adc_events' in ref and len(ref['adc_events']) > 0:
                ref_adc_list = ref['adc_events']
                for adc in wf.adc_events:
                    if adc.num_samples <= 0 or adc.duration_us <= 0:
                        continue
                    # Find closest reference event by onset time
                    best = min(ref_adc_list,
                               key=lambda ev: abs(ev[0] - adc.onset_us))
                    if abs(best[0] - adc.onset_us) > 100.0:
                        continue  # no matching ref event within 100 µs
                    freq_err = abs(adc.freq_offset_hz - best[1])
                    phase_err = abs(np.angle(
                        np.exp(1j * (adc.phase_offset_rad - best[2]))))
                    adc_tol_freq = 1.0   # Hz
                    adc_tol_phase = 0.01  # rad
                    if freq_err > adc_tol_freq or phase_err > adc_tol_phase:
                        pfx = f'TR {tr_idx}: ' if multi_tr else ''
                        messages.append(
                            f'{pfx}ADC offset mismatch: freq_err={freq_err:.2f} Hz, '
                            f'phase_err={phase_err:.4f} rad'
                        )
            if messages:
                fail = True
                fail_msg = messages
                fail_subseq_idx = ss_idx
                fail_tr_idx = tr_idx
                break
        if fail:
            break

    # ── report ────────────────────────────────────────────────────────
    if fail:
        print('Validation failed:')
        for msg in fail_msg:
            print(f'  - {msg}')
        # Always plot the failing (subseq, TR) when do_plot=True
        if do_plot:
            _validation_plot(
                seq,
                sequence_idx=fail_subseq_idx,
                xml_path=xml_path,
                tr_idx=fail_tr_idx,
                max_grad_mT_per_m=_max_grad_plot,
                ok=False,
                messages=fail_msg,
                num_averages=num_averages,
            )
        raise RuntimeError('Validation failed:\n' + '\n'.join(fail_msg))

    print('Validation passed.')
    # Plot only when a specific (subseq, TR) was requested
    if do_plot and targeted:
        _validation_plot(
            seq,
            sequence_idx=subsequence_idx,
            xml_path=xml_path,
            tr_idx=tr_instance,
            max_grad_mT_per_m=_max_grad_plot,
            ok=True,
            messages=[],
            num_averages=num_averages,
        )

# ── validation plot (private) ────────────────────────────────────────

def _validation_plot(
    seq,
    *,
    sequence_idx,
    xml_path,
    tr_idx,
    max_grad_mT_per_m,
    ok,
    messages,
    num_averages=0,
):
    """Create the comparison overlay using :func:`_plot.plot`."""
    from ._plot import plot as _plot_impl

    # Branch 1: C-backend base figure
    handle = _plot_impl(
        seq,
        subsequence_idx=sequence_idx,
        tr_idx=tr_idx,
        collapse_delays=False,
        show_segments=False,
        show_blocks=False,
        show_slew=False,
        show_rf_centers=False,
        show_echoes=False,
        max_grad_mT_per_m=max_grad_mT_per_m,
        num_averages=num_averages,
        max_slew_T_per_m_per_s=None,
    )

    # Branch 2 (pypulseq) or Branch 3 (XML) overlay
    label = 'XML' if xml_path is not None else 'pypulseq'
    overlay_source = str(xml_path) if xml_path is not None else seq._seqs[sequence_idx]
    _plot_impl(overlay_source, fig=handle, label=label)

    # Annotate pass / fail
    status = 'PASS' if ok else 'FAIL'
    colour = 'green' if ok else 'red'
    handle.fig.suptitle(
        f'validate()  \u2014  {status}  (subseq={sequence_idx}, TR={tr_idx})',
        fontweight='bold', color=colour,
    )
    # Only show error messages in the plot if you want detailed debugging; suppress for user-facing plots
    # (per user request, do not print per-TR error messages on the plot)
    for gname in ('gx', 'gy', 'gz'):
        ax = handle.axes[gname]

        def _line_peak_abs(line):
            y = np.asarray(line.get_ydata(), dtype=float)
            if y.size == 0:
                return 0.0
            y = y[np.isfinite(y)]
            if y.size == 0:
                return 0.0
            # Ignore flat horizontal reference lines (e.g. +/- max_grad),
            # otherwise all gradient panels inherit the same y-scale.
            if np.ptp(y) == 0.0:
                return 0.0
            return float(np.max(np.abs(y)))

        ymax = max(
            max((_line_peak_abs(line) for line in ax.get_lines()), default=0.0),
            1.0,
        )
        ax.set_ylim(-ymax * 1.1, ymax * 1.1)
    handle.fig.tight_layout(rect=[0, 0, 1, 0.95])
