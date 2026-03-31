"""Validation: compare pulserver waveforms against a reference source."""

__all__ = ['validate']

from pathlib import Path
from typing import Union

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

    abs_t1_s = abs_t0_s + window_duration_us * 1e-6

    result = src_seq.waveforms_and_times(append_RF=True, time_range=[abs_t0_s, abs_t1_s])
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
                t_rel = (arr[0] - abs_t0_s) * 1e6
                a_rel = arr[1] * hz_to_mT_per_m
                ref[ch_name] = _clip_window(t_rel, a_rel)
                continue
        ref[ch_name] = (np.empty(0), np.empty(0))

    if len(channels) > 3:
        arr = channels[3]
        if arr.shape[1] > 0:
            t_rel = (arr[0] - abs_t0_s) * 1e6
            a_rel = np.abs(arr[1])  # Hz, match C-backend
            ref['rf_mag'] = _clip_window(t_rel, a_rel)
        else:
            ref['rf_mag'] = (np.empty(0), np.empty(0))
    else:
        ref['rf_mag'] = (np.empty(0), np.empty(0))

    return ref


# ── reference extraction ─────────────────────────────────────────────


def _pypulseq_reference(seq, sequence_idx, tr_idx, tr_info):
    """Extract reference waveforms from pypulseq for a single TR.

    Returns
    -------
    dict
        Mapping of channel name to ``(time_us, amplitude)`` tuples.
        Times are relative to TR start; amplitudes in mT/m (grad) or
        µT (RF magnitude).
    """
    t0 = _abs_tr_start_s(
        seq._seqs[sequence_idx],
        tr_info['num_prep_blocks'],
        tr_idx,
        tr_info['tr_duration_us'],
        tr_info.get('imaging_tr_start'),
    )
    return _pypulseq_reference_window(seq, sequence_idx, t0, tr_info['tr_duration_us'])


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
    sequence_idx: int = 0,
    xml_path: Union[str, Path, None] = None,
    do_plot: bool = False,
    tr_range: tuple[int, int] = (0, 1),
    show_rf_centers: bool = False,
    show_echoes: bool = False,
    show_segments: bool = True,
    show_blocks: bool = False,
    max_grad_mT_per_m: Union[float, bool, None] = True,
    grad_atol: float | None = None,
    rf_rms_percent: float = 10.0,
) -> None:
    """Compare pulserver C-backend waveforms against a reference.

    The reference is either the source pypulseq ``Sequence`` (default,
    ``xml_path=None``) or an exported XML waveform file.  Gradient
    channels are compared by maximum absolute error; RF magnitude is
    compared by RMS percentage error.

    When *do_plot* is ``True`` the base TR waveform plot is created from
    the C backend (using :func:`_plot.plot` Branch 1) and the reference
    is overlaid (Branch 2 for pypulseq, Branch 3 for XML).

    # Print first 20 values of RF time and amplitude for reference and test
    print("\n[DEBUG] Reference RF time_us (first 20):", ref_t[:20])
    print("[DEBUG] Reference RF amplitude (first 20):", ref_a[:20])
    print("[DEBUG] Test RF time_us (first 20):", wf.rf_mag.time_us[:20])
    print("[DEBUG] Test RF amplitude (first 20):", wf.rf_mag.amplitude[:20])

    Parameters
    ----------
    seq : SequenceCollection
        Sequence collection to validate.
    sequence_idx : int
        Subsequence index (0-based, default 0).
    xml_path : str, Path or None
        Path to an XML waveform file.  ``None`` (default) validates
        against the stored pypulseq ``Sequence``.
    do_plot : bool
        Show a comparison overlay plot (default ``False``).
    tr_range : (int, int)
        Half-open TR index range ``[start, stop)`` to validate.
        When nondegenerate prep/cooldown blocks are present, the range is
        adjusted to ``(0, 1)`` to validate the full pass as a single TR.
    show_rf_centers : bool
        Mark RF iso-centres on the plot (default ``False``).
    show_echoes : bool
        Mark echo (ADC centre) on the plot (default ``False``).
    show_segments : bool
        Colour-code gradients by segment (default ``True``).
    show_blocks : bool
        Block-boundary lines on the plot (default ``False``).
    max_grad_mT_per_m : float, bool or None
        Gradient-limit reference line.  ``True`` (default) derives the
        value from system limits; a float is used directly; ``False``
        or ``None`` disables the line.
    grad_atol : float or None
        Absolute gradient error tolerance in mT/m.  ``None`` (default)
        uses ``3 × max_slew × grad_raster_time``.
    rf_rms_percent : float
        RF magnitude percent-RMS error threshold (default 10).

    """
    from ._extension._pulseqlib_wrapper import _find_tr

    sys = seq.system

    # Default gradient tolerance: 3 slew steps (mT/m)
    if grad_atol is None:
        grad_atol = 3.0 * sys.max_slew * sys.grad_raster_time * 1e3

    # Resolve max-grad line for plotting
    if max_grad_mT_per_m is True:
        _max_grad_plot = sys.max_grad * 1e3  # T/m → mT/m
    elif isinstance(max_grad_mT_per_m, (int, float)) and max_grad_mT_per_m is not False:
        _max_grad_plot = float(max_grad_mT_per_m)
    else:
        _max_grad_plot = None

    # TR metadata
    tr_info = _find_tr(seq._cseq, subsequence_idx=sequence_idx)
    tr_start, tr_stop = tr_range
    if tr_start < 0 or tr_stop < tr_start or tr_stop > tr_info['num_trs']:
        raise ValueError(
            f'tr_range {tr_range} out of bounds for subsequence with '
            f'{tr_info["num_trs"]} TRs'
        )

    has_prep = tr_info['num_prep_blocks'] > 0 and not tr_info['degenerate_prep']
    has_cooldown = tr_info['num_cooldown_blocks'] > 0 and not tr_info['degenerate_cooldown']

    if has_prep or has_cooldown:
        # Treat subsequence as single full-pass TR
        tr_start = 0
        tr_stop = 1
        # (prep/cooldown always included by backend)
        multi_tr = False
    else:
        # (prep/cooldown always included by backend)
        multi_tr = (tr_stop - tr_start) > 1

    # ── Validate each TR in range ────────────────────────────
    errors_per_tr: list[dict[str, float]] = []
    messages: list[str] = []

    for tr_idx in range(tr_start, tr_stop):
        # C-backend waveforms (actual amplitudes for this TR)
        wf = get_tr_waveforms(
            seq,
            subsequence_idx=sequence_idx,
            amplitude_mode='actual',  # use actual amplitudes for validation
            tr_index=tr_idx,
        )

        # Reference waveforms
        if xml_path is None:
            # Always use full-pass window: start from block 0
            block_start = 0
            abs_t0 = _abs_block_start_s(seq._seqs[sequence_idx], block_start)
            ref = _pypulseq_reference_window(
                seq,
                sequence_idx,
                abs_t0,
                wf.total_duration_us,
            )
            # Print windowing info for reference and test
            ref_t, ref_a = ref['rf_mag']
            print(f"[DEBUG] TR {tr_idx} abs_t0: {abs_t0}")
            print(f"[DEBUG] TR {tr_idx} wf.total_duration_us: {wf.total_duration_us}")
            print(f"[DEBUG] TR {tr_idx} REF RF time_us: first={ref_t[:5]}, last={ref_t[-5:] if len(ref_t) > 5 else ref_t}")
            print(f"[DEBUG] TR {tr_idx} TEST RF time_us: first={wf.rf_mag.time_us[:5]}, last={wf.rf_mag.time_us[-5:] if len(wf.rf_mag.time_us) > 5 else wf.rf_mag.time_us}")
        else:
            ref = _xml_reference(xml_path)

        # Compare gradients
        tr_err: dict[str, float] = {}
        for ch in ('gx', 'gy', 'gz'):
            ref_t, ref_a = ref[ch]
            test_ch = getattr(wf, ch)
            r, t = _interp_to_ref(ref_t, ref_a, test_ch.time_us, test_ch.amplitude)
            err = float(np.max(np.abs(r - t))) if len(r) > 0 else 0.0
            tr_err[ch] = err
            if err > grad_atol:
                pfx = f'TR {tr_idx}: ' if multi_tr else ''
                messages.append(
                    f'{pfx}{ch} mismatch: max diff {err:.4f} mT/m '
                    f'(tol {grad_atol:.4f} mT/m)'
                )


        # Force reference RF arrays to be real
        ref_t_real = np.real(ref_t)
        ref_a_real = np.real(ref_a)

        test_t = wf.rf_mag.time_us
        test_a = wf.rf_mag.amplitude


        # Interpolate both reference and test RF onto a common grid over the overlapping region
        t_min = max(ref_t_real[0], test_t[0]) if len(ref_t_real) and len(test_t) else 0.0
        t_max = min(ref_t_real[-1], test_t[-1]) if len(ref_t_real) and len(test_t) else 0.0
        if t_max <= t_min:
            messages.append(f"RF time overlap is empty: ref=[{ref_t_real[0] if len(ref_t_real) else 'NA'}, {ref_t_real[-1] if len(ref_t_real) else 'NA'}], test=[{test_t[0] if len(test_t) else 'NA'}, {test_t[-1] if len(test_t) else 'NA'}]")
            errors_per_tr.append(tr_err)
            continue
        t_common = np.arange(np.ceil(t_min), np.floor(t_max) + 1)
        if len(t_common) < 3:
            messages.append(f"RF overlap region too small for robust comparison (len={len(t_common)})")
            errors_per_tr.append(tr_err)
            continue
        # Ignore endpoints (first and last sample) as in MATLAB reference
        t_common = t_common[1:-1]
        if len(t_common) == 0:
            messages.append(f"RF overlap region empty after ignoring endpoints.")
            errors_per_tr.append(tr_err)
            continue
        # Interpolate absolute values
        ref_interp = np.abs(np.interp(t_common, ref_t_real, ref_a_real, left=0.0, right=0.0))
        test_interp = np.abs(np.interp(t_common, test_t, test_a, left=0.0, right=0.0))

        # Print diagnostics
        print(f"[DEBUG] TR {tr_idx} REF RF (interp): len={len(ref_interp)}, min={ref_interp.min() if len(ref_interp) else 'NA'}, max={ref_interp.max() if len(ref_interp) else 'NA'}")
        print(f"[DEBUG] TR {tr_idx} TEST RF (interp): len={len(test_interp)}, min={test_interp.min() if len(test_interp) else 'NA'}, max={test_interp.max() if len(test_interp) else 'NA'}")

        rf_err = _rms_error(ref_interp, test_interp)
        tr_err['rf_mag'] = rf_err
        if rf_err > rf_rms_percent:
            pfx = f'TR {tr_idx}: ' if multi_tr else ''
            messages.append(
                f'{pfx}RF mismatch: {rf_err:.1f}% RMS ' f'(tol {rf_rms_percent:.1f}%)'
            )

        errors_per_tr.append(tr_err)

    ok = len(messages) == 0

    # ── Optional overlay plot ────────────────────────────────
    if do_plot:
        _validation_plot(
            seq,
            sequence_idx=sequence_idx,
            xml_path=xml_path,
            tr_idx=tr_range[0],
            show_rf_centers=show_rf_centers,
            show_echoes=show_echoes,
            show_segments=show_segments,
            show_blocks=show_blocks,
            max_grad_mT_per_m=_max_grad_plot,
            ok=ok,
            messages=messages,
        )

    if ok:
        print('Validation passed.')
        return

    print('Validation failed:')
    for msg in messages:
        print(f'  - {msg}')

    raise RuntimeError('Validation failed:\n' + '\n'.join(messages))


# ── validation plot (private) ────────────────────────────────────────


def _validation_plot(
    seq,
    *,
    sequence_idx,
    xml_path,
    tr_idx,
    show_rf_centers,
    show_echoes,
    show_segments,
    show_blocks,
    max_grad_mT_per_m,
    ok,
    messages,
):
    """Create the comparison overlay using :func:`_plot.plot`."""
    from ._plot import plot as _plot_impl

    # Branch 1: C-backend base figure
    handle = _plot_impl(
        seq,
        subsequence_idx=sequence_idx,
        tr_idx=tr_idx,
        collapse_delays=False,
        show_segments=show_segments,
        show_blocks=show_blocks,
        show_slew=False,
        show_rf_centers=show_rf_centers,
        show_echoes=show_echoes,
        max_grad_mT_per_m=max_grad_mT_per_m,
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
        f'validate()  \u2014  {status}', fontweight='bold', color=colour
    )
    if messages:
        handle.fig.text(
            0.5, 0.01, '\n'.join(messages), ha='center', fontsize=8, color='red'
        )
    handle.fig.tight_layout(rect=[0, 0.05 if messages else 0, 1, 0.95])
