"""Gradient spectrum (mechanical resonances) visualisation for SequenceCollection (plotting only, no pass/fail)."""

__all__ = ['grad_spectrum']

import warnings

import numpy as np

from ._extension._pulseqlib_wrapper import _calc_mech_resonances, _find_tr
from ._helpers import _add_echo_spacing_axis
from ._sequence import SequenceCollection


def _parse_band_axis(tag: object) -> str | None:
    """Normalize a user/system channel tag to one of gx/gy/gz."""
    if isinstance(tag, str):
        t = tag.strip().lower()
        if t in ('gx', 'x'):
            return 'gx'
        if t in ('gy', 'y'):
            return 'gy'
        if t in ('gz', 'z'):
            return 'gz'
    if isinstance(tag, (int, np.integer)):
        idx = int(tag)
        if idx == 0:
            return 'gx'
        if idx == 1:
            return 'gy'
        if idx == 2:
            return 'gz'
    return None


def _plot_grad_spectrum_for_subsequence(
    seq: SequenceCollection,
    *,
    subsequence_idx: int,
    canonical_tr_indices: list[int],
    forbidden_bands: list[tuple[float, float, float] | tuple[float, float, float, str]],
    target_window_size: int,
    spectral_resolution: float,
    max_frequency: float,
    peak_log10_threshold: float | None,
    peak_norm_scale: float | None,
    peak_eps: float | None,
    peak_prominence: float | None,
) -> None:
    """Plot gradient spectra for all canonical TRs of one subsequence, overlaid.

    One figure with 3 axes (Gx, Gy, Gz) is produced per subsequence.  Each
    canonical TR is drawn with a different alpha level so all waveforms are
    visible.  The analytical (Dirichlet) structural spectra are shown as
    dashed overlays.  Forbidden bands are shaded.  Violation candidates are
    annotated with their peak gradient amplitude (mT/m).
    """
    import matplotlib.pyplot as plt

    gamma_hz_per_t = float(seq.system.gamma)
    hz_per_m_to_mT_per_m = 1e3 / gamma_hz_per_t
    n_ctrs = len(canonical_tr_indices)

    axis_names = ('gx', 'gy', 'gz')

    # Normalize forbidden-band input: keep channel metadata for plotting,
    # strip to plain triples for C backend calls.
    plot_bands: list[tuple[float, float, float, str]] = []
    c_forbidden_bands: list[tuple[float, float, float]] = []
    for bi, band in enumerate(forbidden_bands):
        if len(band) == 4:
            b0, b1, b2, b3 = band
            b_axis = _parse_band_axis(b3)
            if b_axis is None:
                b_axis = axis_names[bi % 3]
        else:
            b0, b1, b2 = band
            b_axis = axis_names[bi % 3]

        fmin = float(b0)
        fmax = float(b1)
        alim = float(b2)
        c_forbidden_bands.append((fmin, fmax, alim))
        plot_bands.append((fmin, fmax, alim, b_axis))

    # Collect per-canonical-TR result dicts
    all_rd = []
    for ctr_idx in canonical_tr_indices:
        rd = _calc_mech_resonances(
            seq._cseq,
            subsequence_idx=subsequence_idx,
            canonical_tr_idx=ctr_idx,
            target_window_size=target_window_size,
            target_resolution_hz=spectral_resolution,
            max_freq_hz=max_frequency,
            forbidden_bands=c_forbidden_bands,
            peak_log10_threshold=peak_log10_threshold,
            peak_norm_scale=peak_norm_scale,
            peak_eps=peak_eps,
            peak_prominence=peak_prominence,
        )
        all_rd.append(rd)

    # Use the frequency grid of the first result (all share the same grid)
    rd0 = all_rd[0]
    num_freq_bins = rd0['num_freq_bins']
    frequencies = rd0['freq_min_hz'] + np.arange(num_freq_bins) * rd0['freq_spacing_hz']
    freq_min = 0.0
    freq_max = float(frequencies[-1]) if num_freq_bins > 0 else max_frequency

    axis_labels = {'gx': 'Gx', 'gy': 'Gy', 'gz': 'Gz'}
    colors = {'gx': 'C0', 'gy': 'C1', 'gz': 'C2'}
    resonance_styles = {
        'gx': ((0.0, (1.0, 1.8)), 2.0),
        'gy': ((1.1, (1.0, 2.3)), 1.8),
        'gz': ((2.2, (1.0, 2.8)), 1.6),
    }

    # Alpha levels: first CTR is most opaque; subsequent ones are lighter
    alphas = [max(0.2, 0.5 - 0.15 * k) for k in range(n_ctrs)]

    title_suffix = (
        f'CTR{canonical_tr_indices[0]} - CTR{canonical_tr_indices[-1]}'
        if n_ctrs > 1
        else f'CTR{canonical_tr_indices[0]}'
    )
    title = f'[SS{subsequence_idx}, {title_suffix}] Mechanical Resonances'

    fig, axes = plt.subplots(3, 1, figsize=(12, 8), sharex=True)
    global_ymax = 0.0
    traces_by_axis: dict[
        str, list[tuple[np.ndarray, np.ndarray, float, str | None]]
    ] = {
        'gx': [],
        'gy': [],
        'gz': [],
    }
    axis_to_idx = {'gx': 0, 'gy': 1, 'gz': 2}

    for ax_name in axis_names:
        for k, (ctr_idx, rd) in enumerate(zip(canonical_tr_indices, all_rd, strict=False)):
            label = f'CTR{ctr_idx}' if n_ctrs > 1 else None

            comp_freqs = np.asarray(rd['component_freqs_hz'], dtype=np.float64)
            comp_amps = np.asarray(rd['component_amps'], dtype=np.float64)
            comp_phases = np.asarray(rd['component_phases_rad'], dtype=np.float64)
            comp_widths = np.asarray(rd['component_widths_hz'], dtype=np.float64)
            comp_axes = np.asarray(rd['component_axes'], dtype=np.int32)
            comp_def_ids = np.asarray(rd['component_def_ids'], dtype=np.int32)
            comp_contrib_ids = np.asarray(rd['component_contrib_ids'], dtype=np.int32)

            n_comp = min(
                len(comp_freqs),
                len(comp_amps),
                len(comp_phases),
                len(comp_widths),
                len(comp_axes),
                len(comp_def_ids),
                len(comp_contrib_ids),
            )
            if n_comp == 0:
                f_dense = np.linspace(freq_min, freq_max, 2000)
                rebuilt = np.zeros_like(f_dense)
                traces_by_axis[ax_name].append((f_dense, rebuilt, alphas[k], label))
                continue

            f_dense = np.linspace(freq_min, freq_max, 2000)
            x_half = 0.6033545644016142
            idx_ax = axis_to_idx[ax_name]

            # Stage 1: sum spacing/duration terms within each contribution.
            contrib_curves: dict[int, np.ndarray] = {}
            contrib_def: dict[int, int] = {}
            for ci in range(n_comp):
                if int(comp_axes[ci]) != idx_ax:
                    continue
                fpk = float(comp_freqs[ci])
                apk = float(comp_amps[ci])
                ppk = float(comp_phases[ci])
                wpk = float(comp_widths[ci])
                contrib_id = int(comp_contrib_ids[ci])
                def_id = int(comp_def_ids[ci])
                if (
                    not np.isfinite(fpk)
                    or not np.isfinite(apk)
                    or not np.isfinite(ppk)
                    or not np.isfinite(wpk)
                ):
                    continue
                if wpk <= 0.0:
                    raise RuntimeError(
                        'Invalid non-positive component width in analytical export'
                    )
                sinc_scale = wpk / (2.0 * x_half)
                term = apk * np.exp(1j * ppk) * np.sinc((f_dense - fpk) / sinc_scale)
                if contrib_id in contrib_curves:
                    contrib_curves[contrib_id] += term
                else:
                    contrib_curves[contrib_id] = term
                    contrib_def[contrib_id] = def_id

            if len(contrib_curves) == 0:
                rebuilt = np.zeros_like(f_dense)
                traces_by_axis[ax_name].append((f_dense, rebuilt, alphas[k], label))
                continue

            # Stage 2: sum contributions by grad def, then across defs.
            def_curves: dict[int, np.ndarray] = {}
            for contrib_id, curve in contrib_curves.items():
                def_id = contrib_def[contrib_id]
                if def_id in def_curves:
                    def_curves[def_id] += curve
                else:
                    def_curves[def_id] = curve.copy()

            rebuilt_complex = np.zeros_like(f_dense, dtype=np.complex128)
            for curve in def_curves.values():
                rebuilt_complex += curve

            # Plot the inner-run spectral envelope (sum of sinc terms per
            # contribution/def) without applying the outer-TR Dirichlet D_K.
            # D_K has FWHM = 1.2067 * seq_df / K, which is typically sub-Hz
            # (e.g. 0.12 Hz for K=200, TR=50ms) — far below the f_dense grid
            # spacing — so multiplying by D_K produces sub-pixel spikes that
            # make the traces invisible.  The outer-TR modulation is already
            # encoded in the candidate positions (found by C at outer harmonics):
            # the vertical candidate lines show exactly which peaks C selected.
            rebuilt = np.abs(rebuilt_complex)

            local_max = float(np.max(rebuilt)) if rebuilt.size > 0 else 0.0
            if np.isfinite(local_max) and local_max > global_ymax:
                global_ymax = local_max
            traces_by_axis[ax_name].append((f_dense, rebuilt, alphas[k], label))

    band_colors = [colors[b[3]] for b in plot_bands]

    def _candidate_is_viol(rd: dict, ci: int, cf: float) -> bool:
        arr = np.asarray(rd.get('candidate_grad_amps', []), dtype=np.float64)
        if ci >= len(arr):
            return False
        amp_hz = float(arr[ci])
        return any(blo <= cf <= bhi and amp_hz > blim for blo, bhi, blim, _ in plot_bands)

    for panel_idx, ax_name in enumerate(axis_names):
        ax = axes[panel_idx]
        color = colors[ax_name]

        for f_dense, rebuilt, alpha_k, label in traces_by_axis[ax_name]:
            y_plot = rebuilt / global_ymax if global_ymax > 0.0 else rebuilt
            ax.plot(
                f_dense,
                y_plot,
                color=color,
                lw=0.9,
                alpha=alpha_k,
                label=label,
            )

        for bi, band in enumerate(plot_bands):
            ax.axvspan(
                band[0],
                band[1],
                alpha=0.10,
                color=band_colors[bi],
                zorder=0,
            )

        # Draw channel-colored resonance markers using per-channel candidate
        # amplitudes, then propagate selected channel markers to all subplots.
        # A channel is shown for candidate i only if its amplitude is at least
        # rel_axis_thresh of the candidate max across x/y/z.
        rel_axis_thresh = 0.25
        drawn_freqs: dict[tuple[float, str], int] = {}
        for rd in all_rd:
            cf_arr = np.asarray(rd.get('candidate_freqs', []), dtype=np.float64)
            a_gx = np.asarray(rd.get('candidate_amps_gx', []), dtype=np.float64)
            a_gy = np.asarray(rd.get('candidate_amps_gy', []), dtype=np.float64)
            a_gz = np.asarray(rd.get('candidate_amps_gz', []), dtype=np.float64)
            n_c = len(cf_arr)
            for ci in range(n_c):
                cf = float(cf_arr[ci])
                if cf < freq_min or cf > freq_max:
                    continue

                ax_amp = {
                    'gx': float(a_gx[ci]) if ci < len(a_gx) else 0.0,
                    'gy': float(a_gy[ci]) if ci < len(a_gy) else 0.0,
                    'gz': float(a_gz[ci]) if ci < len(a_gz) else 0.0,
                }
                amax = max(ax_amp.values())
                if not np.isfinite(amax) or amax <= 0.0:
                    continue

                for src_axis in axis_names:
                    if ax_amp[src_axis] >= rel_axis_thresh * amax:
                        key = (cf, src_axis)
                        iv = 1 if _candidate_is_viol(rd, ci, cf) else 0
                        drawn_freqs[key] = max(drawn_freqs.get(key, 0), iv)

        for (cf, src_axis), is_viol in sorted(
            drawn_freqs.items(), key=lambda kv: (kv[0][0], kv[0][1])
        ):
            ls, lw_base = resonance_styles[src_axis]
            ax.axvline(
                cf,
                color=colors[src_axis],
                linestyle=ls,
                linewidth=lw_base + (0.4 if is_viol else 0.0),
                alpha=0.9,
                zorder=1,
            )

        ax.set_xlim(freq_min, freq_max)
        ax.set_ylabel(f'{axis_labels[ax_name]} (norm.)')
        ax.grid(True, alpha=0.3)
    for ax in axes:
        ax.set_ylim(0.0, 1.0)

    # --- Annotations on top panel only ---
    ax_top = axes[0]

    # Per-candidate grad-amp labels (mT/m) — worst case across all CTRs,
    # only for candidates that violate any forbidden band.
    viol_amps: dict[float, float] = {}
    for rd in all_rd:
        cf_arr = np.asarray(rd.get('candidate_freqs', []), dtype=np.float64)
        amp_arr = np.asarray(rd.get('candidate_grad_amps', []), dtype=np.float64)
        for ci in range(len(cf_arr)):
            cf = float(cf_arr[ci])
            if cf < freq_min or cf > freq_max:
                continue
            if ci >= len(amp_arr):
                continue
            gamp_hz = float(amp_arr[ci])
            if gamp_hz <= 0.0:
                continue
            if not any(
                blo <= cf <= bhi and gamp_hz > blim for blo, bhi, blim, _ in plot_bands
            ):
                continue
            gamp = gamp_hz * hz_per_m_to_mT_per_m
            viol_amps[cf] = max(viol_amps.get(cf, 0.0), gamp)
    for cf, gamp_mT in viol_amps.items():
        ax_top.annotate(
            f'{gamp_mT:.1f} mT/m',
            xy=(cf, 1.02),
            xycoords=('data', 'axes fraction'),
            ha='center',
            va='bottom',
            fontsize=7,
            fontweight='bold',
            color='#333333',
            alpha=0.9,
            clip_on=False,
        )

    # Per-forbidden-band: max candidate amplitude within band (mT/m),
    # worst case across all canonical TRs.
    band_label_base_y = {'gx': 0.95, 'gy': 0.86, 'gz': 0.77}
    band_label_step = 0.09
    band_axis_counts: dict[str, int] = {'gx': 0, 'gy': 0, 'gz': 0}
    for bi, band in enumerate(plot_bands):
        band_lo, band_hi, band_limit = band[0], band[1], band[2]
        band_axis = band[3]
        band_limit_mT = band_limit * hz_per_m_to_mT_per_m
        max_cand_mT = 0.0
        for rd in all_rd:
            cf_arr = np.asarray(rd.get('candidate_freqs', []), dtype=np.float64)
            amp_arr = np.asarray(rd.get('candidate_grad_amps', []), dtype=np.float64)
            for ci in range(len(cf_arr)):
                cf = float(cf_arr[ci])
                if band_lo <= cf <= band_hi and ci < len(amp_arr):
                    amp_hz = float(amp_arr[ci])
                    if amp_hz > band_limit:
                        val = amp_hz * hz_per_m_to_mT_per_m
                        if val > max_cand_mT:
                            max_cand_mT = val
        band_center = 0.5 * (band_lo + band_hi)
        label = f'limit: {band_limit_mT:.1f} mT/m'
        if max_cand_mT > 0.0:
            label += f'\nmax: {max_cand_mT:.1f} mT/m'

        y_anchor = (
            band_label_base_y[band_axis] - band_axis_counts[band_axis] * band_label_step
        )
        y_anchor = max(0.20, y_anchor)
        band_axis_counts[band_axis] += 1
        ax_top.annotate(
            label,
            xy=(band_center, y_anchor),
            xycoords=('data', 'axes fraction'),
            ha='center',
            va='top',
            fontsize=7,
            color=band_colors[bi],
            alpha=0.85,
            bbox={
                'boxstyle': 'round,pad=0.2',
                'fc': 'white',
                'ec': band_colors[bi],
                'alpha': 0.6,
            },
        )

    # Title and axis labels
    axes[0].set_title(title)
    axes[2].set_xlabel('Frequency (Hz)')

    # Single bottom legend: explicit channel mapping for spectra, resonance,
    # and forbidden-band shading.
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch

    legend_handles = [
        Line2D([0], [0], color=colors['gx'], lw=1.6, label='Gx spectrum'),
        Line2D([0], [0], color=colors['gy'], lw=1.6, label='Gy spectrum'),
        Line2D([0], [0], color=colors['gz'], lw=1.6, label='Gz spectrum'),
        Line2D(
            [0],
            [0],
            color=colors['gx'],
            lw=resonance_styles['gx'][1],
            linestyle=resonance_styles['gx'][0],
            label='Gx resonance',
        ),
        Line2D(
            [0],
            [0],
            color=colors['gy'],
            lw=resonance_styles['gy'][1],
            linestyle=resonance_styles['gy'][0],
            label='Gy resonance',
        ),
        Line2D(
            [0],
            [0],
            color=colors['gz'],
            lw=resonance_styles['gz'][1],
            linestyle=resonance_styles['gz'][0],
            label='Gz resonance',
        ),
        Patch(facecolor=colors['gx'], alpha=0.10, edgecolor='none', label='Gx bands'),
        Patch(facecolor=colors['gy'], alpha=0.10, edgecolor='none', label='Gy bands'),
        Patch(facecolor=colors['gz'], alpha=0.10, edgecolor='none', label='Gz bands'),
    ]
    fig.legend(
        handles=legend_handles,
        loc='lower center',
        ncol=3,
        frameon=True,
        fontsize=8,
        bbox_to_anchor=(0.5, 0.01),
    )

    # Echo spacing secondary axis on top panel only
    _add_echo_spacing_axis(axes[0], freq_min, freq_max)

    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        fig.tight_layout(rect=(0.0, 0.07, 1.0, 1.0))


def grad_spectrum(
    seq: SequenceCollection,
    *,
    sequence_idx: int | None = None,
    forbidden_bands: (
        list[tuple[float, float, float] | tuple[float, float, float, str]] | None
    ) = None,
    window_duration: float = 25.0e-3,
    spectral_resolution: float = 5.0,
    max_frequency: float = 3000.0,
    peak_log10_threshold: float | None = None,
    peak_norm_scale: float | None = None,
    peak_eps: float | None = None,
    peak_prominence: float | None = None,
) -> None:
    """Plot gradient spectrum (mechanical resonance candidates) for the canonical TR.

    Creates a single-panel figure with three overlaid spectra (Gx, Gy, Gz) showing
    surviving peak candidates as vertical lines, forbidden bands as shaded regions,
    and max gradient amplitude annotations (mT/m) per candidate and per band.

    Two frequency axes are provided: frequency in Hz (bottom) and the corresponding
    echo spacing in ms (top). The vertical axis is in arbitrary units (a.u.).

    No pass/fail check is performed — use
    :meth:`SequenceCollection.check` for that.

    Parameters
    ----------
    seq : SequenceCollection
        The sequence to analyse.
    sequence_idx : int, optional
        Subsequence index (0-based, default all).
    forbidden_bands : list of (freq_min, freq_max, max_amplitude), optional
        Forbidden frequency bands.  Each tuple gives
        ``(freq_min_Hz, freq_max_Hz, max_allowed_amplitude_Hz_per_m)``.
        Drawn as shaded regions; the max allowed amplitude is labelled in mT/m.
    window_duration : float
        Sliding-window size in seconds (default 25 ms, used for spectral resolution).
    spectral_resolution : float
        Target frequency resolution in Hz (default 5 Hz).
    max_frequency : float
        Upper frequency limit in Hz (default 3000 Hz).
    peak_log10_threshold : float, optional
        Resonance detector threshold in log10 space.
    peak_norm_scale : float, optional
        Normalization scale used before the log transform in resonance detection.
    peak_eps : float, optional
        Positive epsilon added before log transform for numerical stability.
    peak_prominence : float, optional
        Minimum prominence (in log10 units) for a peak to be retained.
    """
    if forbidden_bands is None:
        forbidden_bands = []

    grad_raster_time = seq.system.grad_raster_time
    target_window_size = int(2.0 * window_duration / grad_raster_time)

    if sequence_idx is None:
        subsequence_indices = list(range(seq.num_sequences))
    else:
        if sequence_idx < 0 or sequence_idx >= seq.num_sequences:
            raise ValueError(
                f'sequence_idx={sequence_idx} out of range for {seq.num_sequences} subsequences'
            )
        subsequence_indices = [int(sequence_idx)]

    for ss_idx in subsequence_indices:
        tr_info = _find_tr(seq._cseq, subsequence_idx=ss_idx)
        num_canonical = int(tr_info.get('num_canonical_trs', 1))
        _plot_grad_spectrum_for_subsequence(
            seq,
            subsequence_idx=ss_idx,
            canonical_tr_indices=list(range(num_canonical)),
            forbidden_bands=forbidden_bands,
            target_window_size=target_window_size,
            spectral_resolution=spectral_resolution,
            max_frequency=max_frequency,
            peak_log10_threshold=peak_log10_threshold,
            peak_norm_scale=peak_norm_scale,
            peak_eps=peak_eps,
            peak_prominence=peak_prominence,
        )


# Alias
mechanical_resonances = grad_spectrum
