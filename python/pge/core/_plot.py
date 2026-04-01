"""TR waveform plotting for SequenceCollection."""

__all__ = ['plot']

from pathlib import Path

import numpy as np

from ._sequence import SequenceCollection
from ._waveforms import ChannelWaveform, get_tr_waveforms

# Matplotlib is imported lazily to avoid hard dependency at import time.

# ── MATLAB lines-like palette for style parity ───────────────────────
_MATLAB_LINES = [
    '#0072BD',
    '#D95319',
    '#EDB120',
    '#7E2F8E',
    '#77AC30',
    '#4DBEEE',
    '#A2142F',
]

_MAIN_LINEWIDTH = 1.1
_OVERLAY_LINEWIDTH = 1.0

# Physical-unit conversion constants
_GAMMA_DEFAULT = 42.576e6  # Hz/T


def _seg_color(idx: int) -> str:
    """Return a colour for segment index *idx* (cycling)."""
    if idx < 0:
        return '#aaaaaa'  # prep/cooldown
    return _MATLAB_LINES[idx % len(_MATLAB_LINES)]


def _segment_idx_for_time(blocks, t_us: float) -> int:
    """Return segment index for a time sample based on block spans."""
    for blk in blocks:
        blk_start = blk.start_us
        blk_end = blk_start + blk.duration_us
        if (t_us >= blk_start - 0.5) and (t_us <= blk_end + 0.5):
            return blk.segment_idx
    return -1


def _slew_rate(ch: ChannelWaveform) -> ChannelWaveform:
    """Compute numerical derivative dA/dt -> T/m/s from mT/m over us."""
    if ch.time_us.size < 2:
        return ChannelWaveform(
            time_us=np.empty(0, dtype=np.float32),
            amplitude=np.empty(0, dtype=np.float32),
        )
    dt = np.diff(ch.time_us) * 1e-3  # us -> ms
    da = np.diff(ch.amplitude)
    dt_safe = np.where(np.abs(dt) < 1e-12, 1e-12, dt)
    slew = da / dt_safe  # mT/m per us = 1e3 T/m/s
    slew_T_per_m_per_s = slew * 1e3
    mid_t = 0.5 * (ch.time_us[:-1] + ch.time_us[1:])
    return ChannelWaveform(
        time_us=mid_t.astype(np.float32),
        amplitude=slew_T_per_m_per_s.astype(np.float32),
    )


class PlotHandle:
    """Lightweight handle returned by :func:`plot` carrying figure + metadata.

    Attributes
    ----------
    fig : matplotlib.figure.Figure
    axes : dict[str, matplotlib.axes.Axes]
    tr_duration_us : float
    num_trs : int
    first_tr_start_us : float
    _tr_start_abs_s : float
        Absolute start time of the displayed TR within the full
        pypulseq sequence (seconds).  Used by Branch 2 overlay to
        compute the correct ``time_range`` for
        ``waveforms_and_times``.
    """

    def __init__(
        self,
        fig,
        axes,
        tr_duration_us,
        num_trs,
        first_tr_start_us,
        *,
        _tr_start_abs_s=0.0,
    ):
        self.fig = fig
        self.axes = axes
        self.tr_duration_us = tr_duration_us
        self.num_trs = num_trs
        self.first_tr_start_us = first_tr_start_us
        self._tr_start_abs_s = _tr_start_abs_s


def plot(
    source,
    *,
    subsequence_idx: int = 0,
    tr_idx=0,
    collapse_delays: bool = True,
    show_segments: bool = True,
    show_blocks: bool = False,
    show_slew: bool = False,
    show_rf_centers: bool = False,
    show_echoes: bool = False,
    max_grad_mT_per_m=None,
    max_slew_T_per_m_per_s=None,
    time_unit: str = 'ms',
    figsize=None,
    fig=None,
    label=None,
):
    """Plot native-timing TR waveforms.

    Layout is (3, 2): first column has RF magnitude (top), RF phase
    (middle), ADC mask (bottom); second column has Gx, Gy, Gz.

    Parameters
    ----------
    source : SequenceCollection or pp.Sequence or str/Path
        Primary data source. A :class:`SequenceCollection` creates a
        fresh figure. A pypulseq ``Sequence`` or XML file path can be
        plotted directly (auto-creates a figure) or overlaid on an
        existing figure via ``fig``.
    subsequence_idx : int
        Subsequence index (default 0).
    tr_idx : int or {'max_pos', 'zero_var'}
        TR index (0-based) or amplitude-mode string.
    collapse_delays : bool
        Shrink pure-delay blocks to 0.1 ms at C level (default True).
    show_segments : bool
        Colour-code gradient waveforms by segment index.
    show_blocks : bool
        Draw vertical dotted lines at block boundaries (default False).
    show_slew : bool
        Overlay slew rate on gradient panels.
    show_rf_centers : bool
        Mark RF iso-centres on the RF magnitude subplot (default False).
    show_echoes : bool
        Mark echo (ADC centre) on the ADC panel (default False).
    max_grad_mT_per_m : float or None
        Horizontal reference line for max gradient amplitude.
    max_slew_T_per_m_per_s : float or None
        Horizontal reference for max slew rate (only if show_slew).
    time_unit : str
        ``'ms'`` (default) or ``'us'``.
    figsize : tuple or None
        Figure size.  Default ``(14, 8)``.
    fig : PlotHandle or None
        Existing plot handle for overlay.  Forbidden when *source* is
        a SequenceCollection. Optional for pypulseq Sequence/XML; when
        omitted, a new figure is created automatically.
    label : str or None
        Legend label for overlay traces.

    Returns
    -------
    PlotHandle
        Handle containing figure, axes, and TR metadata.
    """
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle
    
    is_collection = isinstance(source, SequenceCollection)

    # Lazy check for pypulseq
    try:
        import pypulseq as pp

        # SequenceCollection inherits from pp.Sequence; keep it in Branch 1.
        is_pulseq = (not is_collection) and isinstance(source, pp.Sequence)
    except ImportError:
        is_pulseq = False

    is_xml = isinstance(source, (str, Path))

    # ── Validate fig argument ──
    if is_collection and fig is not None:
        raise ValueError(
            "Cannot pass 'fig' handle when source is a SequenceCollection. "
            "Overlay is only supported for pypulseq Sequence or XML sources."
        )
    if (is_pulseq or is_xml) and fig is None:
        if is_pulseq:
            source = SequenceCollection(source)
            is_collection = True
            is_pulseq = False
            is_xml = False
        else:
            # If a .seq path is provided, use the same base path as collections.
            try:
                source = SequenceCollection(str(source))
                is_collection = True
                is_pulseq = False
                is_xml = False
            except Exception:
                fig = _new_empty_plot_handle(figsize=figsize)

    # ── Determine amplitude mode and tr_index ──
    if isinstance(tr_idx, str):
        amplitude_mode = tr_idx  # 'max_pos' or 'zero_var'
        tr_index = 0
    else:
        # Support negative indexing for tr_index
        tr_info = None
        if 'tr_info' not in locals():
            from ._extension._pulseqlib_wrapper import _find_tr
            if is_collection:
                tr_info = _find_tr(source._cseq, subsequence_idx=subsequence_idx)
            else:
                tr_info = None
        num_trs = tr_info['num_trs'] if tr_info else None
        idx = int(tr_idx)
        if num_trs is not None and idx < 0:
            idx = num_trs + idx
        amplitude_mode = 'actual' if idx > 0 else 'max_pos'
        tr_index = idx

    # Always plot the canonical TR window(s) as defined by the C backend
    from ._extension._pulseqlib_wrapper import _find_tr
    tr_info = None
    if is_collection:
        tr_info = _find_tr(source._cseq, subsequence_idx=subsequence_idx)
        # For 'max_pos' or 'zero_var', plot all canonical TRs (one per unique shot pattern)
        # For 'actual', plot only the requested TR
        if amplitude_mode in ('max_pos', 'zero_var'):
            wf = get_tr_waveforms(
                source,
                subsequence_idx=subsequence_idx,
                amplitude_mode=amplitude_mode,
                tr_index=0,
                collapse_delays=collapse_delays,
            )
            # Diagnostic: print block info for debugging canonical TR window
            print(f"[DEBUG] Number of waveform blocks: {len(wf.blocks)}")
            if wf.blocks:
                print(f"[DEBUG] Block start_us: {[b.start_us for b in wf.blocks]}")
                print(f"[DEBUG] Block duration_us: {[b.duration_us for b in wf.blocks]}")
            tr_dur = tr_info['tr_duration_us']
            num_trs = tr_info['num_trs']
            first_tr_start_us = wf.blocks[0].start_us if wf.blocks else 0.0
            from ._validate import _abs_tr_start_s
            _abs_s = _abs_tr_start_s(
                source._seqs[subsequence_idx],
                tr_info['num_prep_blocks'],
                0,
                tr_info['tr_duration_us'],
                tr_info.get('imaging_tr_start'),
            )
        else:
            wf = get_tr_waveforms(
                source,
                subsequence_idx=subsequence_idx,
                amplitude_mode=amplitude_mode,
                tr_index=tr_index,
                collapse_delays=collapse_delays,
            )
            # Mask to the TR window for the requested TR
            tr_dur = tr_info['tr_duration_us']
            t0 = wf.blocks[0].start_us if wf.blocks else 0.0
            t1 = t0 + tr_dur
            def _mask_tr(t, a):
                t = np.asarray(t)
                a = np.asarray(a)
                mask = (t >= t0 - 1e-3) & (t <= t1 + 1e-3)
                return t[mask], a[mask]
            wf.rf_mag.time_us, wf.rf_mag.amplitude = _mask_tr(wf.rf_mag.time_us, wf.rf_mag.amplitude)
            wf.rf_phase.time_us, wf.rf_phase.amplitude = _mask_tr(wf.rf_phase.time_us, wf.rf_phase.amplitude)
            wf.gx.time_us, wf.gx.amplitude = _mask_tr(wf.gx.time_us, wf.gx.amplitude)
            wf.gy.time_us, wf.gy.amplitude = _mask_tr(wf.gy.time_us, wf.gy.amplitude)
            wf.gz.time_us, wf.gz.amplitude = _mask_tr(wf.gz.time_us, wf.gz.amplitude)
            wf.adc_events = [adc for adc in wf.adc_events if t0 <= adc.onset_us <= t1]
            wf.blocks = [blk for blk in wf.blocks if t0 <= blk.start_us < t1]
            num_trs = 1
            first_tr_start_us = t0
            from ._validate import _abs_tr_start_s
            _abs_s = _abs_tr_start_s(
                source._seqs[subsequence_idx],
                tr_info['num_prep_blocks'],
                tr_index,
                tr_info['tr_duration_us'],
                tr_info.get('imaging_tr_start'),
            )

    # ── Time helpers ──
    t_scale = 1e-3 if time_unit == 'ms' else 1.0
    t_label = 'Time (ms)' if time_unit == 'ms' else 'Time (us)'

    # ================================================================
    #  Branch 1: SequenceCollection -> fresh figure
    # ================================================================
    if is_collection:
        # Always include non-degenerate prep/cooldown
        from ._extension._pulseqlib_wrapper import _find_tr

        tr_info = _find_tr(source._cseq, subsequence_idx=subsequence_idx)

        wf = get_tr_waveforms(
            source,
            subsequence_idx=subsequence_idx,
            amplitude_mode=amplitude_mode,
            tr_index=tr_index,
            collapse_delays=collapse_delays,
        )
        tr_dur = tr_info['tr_duration_us']
        num_trs = tr_info['num_trs']
        # Compute actual start time of first TR from waveform blocks
        first_tr_start_us = 0.0
        for blk in wf.blocks:
            if blk.segment_idx >= 0:
                first_tr_start_us = blk.start_us
                break

        # Absolute start time of displayed TR within the pypulseq seq
        # (needed for Branch 2 overlay alignment).
        from ._validate import _abs_tr_start_s

        _abs_s = _abs_tr_start_s(
            source._seqs[subsequence_idx],
            tr_info['num_prep_blocks'],
            tr_index,
            tr_info['tr_duration_us'],
            tr_info.get('imaging_tr_start'),
        )

        if figsize is None:
            figsize = (14, 8)

        fig_obj, axes_grid = plt.subplots(3, 2, figsize=figsize, sharex=True)
        fig_obj.patch.set_facecolor('white')
        axes = {
            'rf_mag': axes_grid[0, 0],
            'rf_phase': axes_grid[1, 0],
            'adc': axes_grid[2, 0],
            'gx': axes_grid[0, 1],
            'gy': axes_grid[1, 1],
            'gz': axes_grid[2, 1],
        }

        def _t(t_us):
            return np.asarray(t_us) * t_scale

        # ── RF magnitude ──
        ax = axes['rf_mag']
        ch = wf.rf_mag
        if ch.time_us.size > 0:
            ax.plot(
                _t(ch.time_us),
                ch.amplitude,
                color='gray',
                linewidth=_MAIN_LINEWIDTH,
                linestyle='-',
                label='pulserver',
            )
        ax.set_ylabel('|RF| (uT)')

        # ── RF phase ──
        ax = axes['rf_phase']
        ch = wf.rf_phase
        # Plot RF phase directly from C backend (already includes all offsets)
        if ch.time_us.size > 0:
            rf_phase_wrapped = (ch.amplitude + np.pi) % (2 * np.pi) - np.pi
            ax.plot(
                _t(ch.time_us),
                rf_phase_wrapped,
                color='gray',
                linewidth=_MAIN_LINEWIDTH,
                linestyle='-',
                label='pulserver',
            )
        ax.set_ylabel('RF phase (rad)')
        ax.set_yticks([-np.pi, 0, np.pi])
        ax.set_yticklabels(['-pi', '0', 'pi'])

        # ── ADC phase envelope ──
        ax = axes['adc']
        # Only plot for real ADC events (not dummy TRs)
        for adc in wf.adc_events:
            if adc.num_samples > 0 and adc.duration_us > 0:
                dwell_time = adc.duration_us / adc.num_samples  # us
                t_adc = adc.onset_us + np.arange(adc.num_samples) * dwell_time
                t_adc_s = t_adc * 1e-6
                adc_phase = adc.phase_offset_rad + 2 * np.pi * adc.freq_offset_hz * t_adc_s
                # Wrap phase to [-pi, pi]
                adc_phase_wrapped = (adc_phase + np.pi) % (2 * np.pi) - np.pi
                ax.plot(_t(t_adc), adc_phase_wrapped, color='purple', linewidth=_MAIN_LINEWIDTH, label='pulserver')
                # Shadowed region
                ax.fill_between(_t(t_adc), 0, np.maximum(adc_phase_wrapped, 1.0), color='purple', alpha=0.15)
        ax.set_ylabel('ADC phase (rad)')
        ax.set_yticks([-np.pi, 0, np.pi])
        ax.set_yticklabels(['-pi', '0', 'pi'])


        # ── Gradients (Gx, Gy, Gz) ──
        axis_color = {
            'gx': _MATLAB_LINES[0],
            'gy': _MATLAB_LINES[1],
            'gz': _MATLAB_LINES[2],
        }
        for gname in ('gx', 'gy', 'gz'):
            ax = axes[gname]
            ch = getattr(wf, gname)
            if ch.time_us.size > 0:
                if show_segments and len(wf.blocks) > 0:
                    _plot_segmented(
                        ax, _t, ch, wf.blocks, linewidth=_MAIN_LINEWIDTH, alpha=1.0
                    )
                else:
                    ax.plot(
                        _t(ch.time_us),
                        ch.amplitude,
                        color='gray',
                        linewidth=_MAIN_LINEWIDTH,
                        linestyle='-',
                        label='pulserver',
                    )

            ax.set_ylabel(f'{gname.upper()} (mT/m)')

            if max_grad_mT_per_m is not None:
                ax.axhline(max_grad_mT_per_m, color='gray', ls='--', lw=0.6, alpha=0.6)
                ax.axhline(-max_grad_mT_per_m, color='gray', ls='--', lw=0.6, alpha=0.6)

            # Slew rate overlay
            if show_slew and ch.time_us.size > 1:
                slew = _slew_rate(ch)
                ax2 = ax.twinx()
                ax2.plot(
                    _t(slew.time_us),
                    slew.amplitude,
                    color='orange',
                    linewidth=0.5,
                    alpha=0.5,
                )
                ax2.set_ylabel('Slew (T/m/s)', color='orange', fontsize=8)
                ax2.tick_params(axis='y', labelcolor='orange', labelsize=7)
                if max_slew_T_per_m_per_s is not None:
                    ax2.axhline(
                        max_slew_T_per_m_per_s,
                        color='orange',
                        ls=':',
                        lw=0.5,
                        alpha=0.5,
                    )
                    ax2.axhline(
                        -max_slew_T_per_m_per_s,
                        color='orange',
                        ls=':',
                        lw=0.5,
                        alpha=0.5,
                    )

        # ── Block boundaries ──
        if show_blocks and len(wf.blocks) > 0:
            for blk in wf.blocks:
                t_start = _t(np.array([blk.start_us]))[0]
                for ax_item in axes.values():
                    ax_item.axvline(t_start, color='k', ls=':', lw=0.3, alpha=0.3)

        # ── RF centres ──
        if show_rf_centers and wf.rf_mag.time_us.size > 0:
            _plot_rf_centers(axes['rf_mag'], wf, _t)

        # ── Echo markers ──
        if show_echoes:
            _plot_echo_markers(axes['adc'], wf, _t)

        # ── Finalise ──
        axes['adc'].set_xlabel(t_label)
        axes['gz'].set_xlabel(t_label)
        for ax_item in axes.values():
            ax_item.grid(True)

        # Add legend to all axes for clarity
        for ax_item in axes.values():
            handles, labels = ax_item.get_legend_handles_labels()
            if handles:
                ax_item.legend(fontsize=8, loc='upper right')

        # --- Overlay for pypulseq (Branch 2) ---
        if not is_collection and is_pulseq and fig is not None:
            # Use the same windowing as validation for the overlay
            from ._validate import _pypulseq_reference_window
            if not hasattr(source, '_seqs'):
                seq_coll = SequenceCollection(source)
            else:
                seq_coll = source
            abs_t0 = fig._tr_start_abs_s
            window_duration_us = fig.tr_duration_us
            ref = _pypulseq_reference_window(seq_coll, subsequence_idx, abs_t0, window_duration_us)
            # Overlay ADC phase for pypulseq
            ax_adc = fig.axes['adc']
            # For pypulseq, try to extract ADC events from the sequence if possible
            # (Assume only one ADC event per TR for now)
            # If not available, skip overlay
            # This is a placeholder: user may need to adapt for their pypulseq structure
            # Here, we just plot zeros for demonstration
            # TODO: Replace with real pypulseq ADC phase extraction if available
            # ax_adc.plot([], [], color='black', linestyle='--', label='pypulseq ADC phase')

        fig_obj.tight_layout()

        return PlotHandle(
            fig=fig_obj,
            axes=axes,
            tr_duration_us=tr_dur,
            num_trs=num_trs,
            first_tr_start_us=first_tr_start_us,
            _tr_start_abs_s=_abs_s,
        )

    # ================================================================
    #  Branch 2: pypulseq Sequence overlay
    # ================================================================
    if is_pulseq:
        assert fig is not None
        # Use the same windowing as validation for the overlay

        from ._validate import _pypulseq_reference_window
        # Ensure we have a SequenceCollection for reference extraction
        if not hasattr(source, '_seqs'):
            seq_coll = SequenceCollection(source)
        else:
            seq_coll = source
        abs_t0 = fig._tr_start_abs_s
        window_duration_us = fig.tr_duration_us
        ref = _pypulseq_reference_window(seq_coll, subsequence_idx, abs_t0, window_duration_us)

        alpha = 1.0
        for ch_name in ('gx', 'gy', 'gz'):
            t_us, amp = ref[ch_name]
            if len(t_us) > 0:
                fig.axes[ch_name].plot(
                    t_us * t_scale,
                    amp,
                    linewidth=_OVERLAY_LINEWIDTH,
                    alpha=alpha,
                    color='black',
                    linestyle='--',
                    label=f'{label or "pypulseq"}',
                )
        # RF magnitude
        t_us, amp = ref['rf_mag']
        if len(t_us) > 0:
            fig.axes['rf_mag'].plot(
                t_us * t_scale,
                amp,
                linewidth=_OVERLAY_LINEWIDTH,
                alpha=alpha,
                color='black',
                linestyle='--',
                label=f'{label or "pypulseq"}',
            )
        # RF phase overlay if available
        if 'rf_phase' in ref:
            t_us, amp = ref['rf_phase']
            if len(t_us) > 0:
                fig.axes['rf_phase'].plot(
                    t_us * t_scale,
                    amp,
                    linewidth=_OVERLAY_LINEWIDTH,
                    alpha=alpha,
                    color='black',
                    linestyle='--',
                    label=f'{label or "pypulseq"}',
                )
        # ADC phase overlay if available
        if 'adc_phase' in ref:
            t_us, amp = ref['adc_phase']
            if len(t_us) > 0:
                fig.axes['adc'].plot(
                    t_us * t_scale,
                    amp,
                    linewidth=_OVERLAY_LINEWIDTH,
                    alpha=alpha,
                    color='black',
                    linestyle='--',
                    label=f'{label or "pypulseq"}',
                )
        # Add legends to all axes
        for ax_item in fig.axes.values():
            handles, labels = ax_item.get_legend_handles_labels()
            if handles:
                ax_item.legend(fontsize=8, loc='upper right')
        return fig

    # ================================================================
    #  Branch 3: XML file overlay
    # ================================================================
    if is_xml:
        assert fig is not None
        _overlay_xml(fig, source, t_scale=t_scale, label=label)
        return fig

    raise TypeError(f'Unsupported source type: {type(source)}')


def _overlay_xml(handle, xml_path, *, t_scale, label):
    """Overlay waveforms from an XML file onto an existing plot."""
    import xml.etree.ElementTree as ET

    tree = ET.parse(str(xml_path))
    root = tree.getroot()

    # G/cm -> mT/m: 1 G/cm = 10 mT/m
    # G -> uT: 1 G = 100 uT
    g_per_cm_to_mT_per_m = 10.0
    g_to_uT = 100.0

    alpha = 0.7

    for ch_name in ('gx', 'gy', 'gz'):
        elem = root.find(f'.//{ch_name}')
        if elem is not None and elem.text:
            data = np.fromstring(elem.text, sep=' ')
            if data.size >= 2:
                n = data.size // 2
                t_us = data[:n]
                amp = data[n:] * g_per_cm_to_mT_per_m
                handle.axes[ch_name].plot(
                    t_us * t_scale,
                    amp,
                    linewidth=_OVERLAY_LINEWIDTH,
                    alpha=alpha,
                    label=label,
                )

    rf_elem = root.find('.//rf')
    if rf_elem is not None and rf_elem.text:
        data = np.fromstring(rf_elem.text, sep=' ')
        if data.size >= 2:
            n = data.size // 2
            t_us = data[:n]
            amp = data[n:] * g_to_uT
            handle.axes['rf_mag'].plot(
                t_us * t_scale,
                np.abs(amp),
                linewidth=_OVERLAY_LINEWIDTH,
                alpha=alpha,
                label=label,
            )

    if label:
        handle.axes['rf_mag'].legend(fontsize=8, loc='upper right')


def _new_empty_plot_handle(*, figsize=None):
    """Create an empty plot handle suitable for direct XML overlays."""
    import matplotlib.pyplot as plt

    if figsize is None:
        figsize = (14, 8)

    fig_obj, axes_grid = plt.subplots(3, 2, figsize=figsize, sharex=True)
    fig_obj.patch.set_facecolor('white')
    axes = {
        'rf_mag': axes_grid[0, 0],
        'rf_phase': axes_grid[1, 0],
        'adc': axes_grid[2, 0],
        'gx': axes_grid[0, 1],
        'gy': axes_grid[1, 1],
        'gz': axes_grid[2, 1],
    }

    axes['rf_mag'].set_ylabel('|RF| (uT)')
    axes['rf_phase'].set_ylabel('RF phase (rad)')
    axes['rf_phase'].set_yticks([-np.pi, 0, np.pi])
    axes['rf_phase'].set_yticklabels(['-pi', '0', 'pi'])
    axes['adc'].set_ylabel('ADC')
    axes['adc'].set_ylim(0, 1)
    axes['adc'].set_yticks([])
    axes['gx'].set_ylabel('GX (mT/m)')
    axes['gy'].set_ylabel('GY (mT/m)')
    axes['gz'].set_ylabel('GZ (mT/m)')
    axes['adc'].set_xlabel('Time (ms)')
    axes['gz'].set_xlabel('Time (ms)')
    for ax in axes.values():
        ax.grid(True)

    fig_obj.tight_layout()
    return PlotHandle(
        fig=fig_obj,
        axes=axes,
        tr_duration_us=0.0,
        num_trs=1,
        first_tr_start_us=0.0,
        _tr_start_abs_s=0.0,
    )


def _plot_rf_centers(ax, wf, t_fn):
    """Plot RF iso-centre markers on the RF magnitude panel.

    Approximates RF centre as the amplitude-weighted mean time for each
    contiguous RF burst.
    """
    if wf.rf_mag.time_us.size == 0:
        return
    t = wf.rf_mag.time_us
    a = np.abs(wf.rf_mag.amplitude)

    # Split into contiguous bursts (gaps > 10 us)
    if len(t) < 2:
        return
    gaps = np.where(np.diff(t) > 10.0)[0]
    starts = np.concatenate([[0], gaps + 1])
    ends = np.concatenate([gaps + 1, [len(t)]])

    for s, e in zip(starts, ends, strict=False):
        seg_t = t[s:e]
        seg_a = a[s:e]
        total = np.sum(seg_a)
        if total > 0:
            center_us = np.sum(seg_t * seg_a) / total
            center_display = t_fn(np.array([center_us]))[0]
            ax.axvline(center_display, color='r', ls='--', lw=1.0, alpha=0.7)
            ax.plot(center_display, np.max(seg_a) * 1.05, 'rv', markersize=6, alpha=0.7)


def _plot_echo_markers(ax, wf, t_fn):
    """Plot echo (ADC midpoint) markers on the ADC panel."""
    for adc in wf.adc_events:
        mid_us = adc.onset_us + adc.duration_us / 2.0
        mid_display = t_fn(np.array([mid_us]))[0]
        ax.axvline(mid_display, color='b', ls='--', lw=1.0, alpha=0.7)
        ax.plot(mid_display, 0.5, 'b^', markersize=6, alpha=0.7)


def _plot_segmented(ax, t_fn, ch, blocks, **kwargs):
    """Plot a channel waveform with colour-coded segments."""
    t = ch.time_us
    a = ch.amplitude

    # Build per-sample segment index
    seg_idx = np.full(len(t), -1, dtype=int)
    for blk in blocks:
        blk_start = blk.start_us
        blk_end = blk_start + blk.duration_us
        mask = (t >= blk_start - 0.5) & (t <= blk_end + 0.5)
        seg_idx[mask] = blk.segment_idx

    if len(t) == 0:
        return
    run_start = 0
    for i in range(1, len(t)):
        if seg_idx[i] != seg_idx[run_start]:
            ax.plot(
                t_fn(t[run_start:i]),
                a[run_start:i],
                color=_seg_color(seg_idx[run_start]),
                **kwargs,
            )
            run_start = i
    ax.plot(
        t_fn(t[run_start:]),
        a[run_start:],
        color=_seg_color(seg_idx[run_start]),
        **kwargs,
    )
