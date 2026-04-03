"""Gradient spectrum (mechanical resonances) visualisation for SequenceCollection (plotting only, no pass/fail)."""

__all__ = ['grad_spectrum']

import warnings

import numpy as np

from ._extension._pulseqlib_wrapper import _calc_mech_resonances, _find_tr
from ._helpers import _add_echo_spacing_axis
from ._sequence import SequenceCollection


def _plot_grad_spectrum_for_subsequence(
    seq: SequenceCollection,
    *,
    subsequence_idx: int,
    canonical_tr_indices: list[int],
    forbidden_bands: list[tuple[float, float, float]],
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
            forbidden_bands=forbidden_bands,
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

    axis_names = ('gx', 'gy', 'gz')
    axis_labels = {'gx': 'Gx', 'gy': 'Gy', 'gz': 'Gz'}
    colors = {'gx': 'C0', 'gy': 'C1', 'gz': 'C2'}

    # Alpha levels: first CTR is most opaque; subsequent ones are lighter
    alphas = [max(0.2, 0.5 - 0.15 * k) for k in range(n_ctrs)]

    title_suffix = (
        f'CTR{canonical_tr_indices[0]}–CTR{canonical_tr_indices[-1]}'
        if n_ctrs > 1
        else f'CTR{canonical_tr_indices[0]}'
    )
    title = f'[SS{subsequence_idx}, {title_suffix}] Mechanical Resonances'

    fig, axes = plt.subplots(3, 1, figsize=(12, 8), sharex=True)

    for panel_idx, ax_name in enumerate(axis_names):
        ax = axes[panel_idx]
        color = colors[ax_name]

        for k, (ctr_idx, rd) in enumerate(zip(canonical_tr_indices, all_rd)):
            num_bins_k = rd['num_freq_bins']
            freqs_k = rd['freq_min_hz'] + np.arange(num_bins_k) * rd['freq_spacing_hz']
            spec = np.asarray(rd[f'spectrum_full_{ax_name}'], dtype=np.float32)
            label = f'CTR{ctr_idx}' if n_ctrs > 1 else None

            # FFT spectrum
            ax.plot(
                freqs_k,
                spec,
                color=color,
                lw=0.8,
                alpha=alphas[k],
                label=label,
            )

            # Dirichlet kernel envelope
            num_instances = int(rd.get('num_instances', 0))
            f0 = float(rd.get('freq_spacing_seq_hz', 0.0))
            if num_instances > 1 and f0 > 0.0:
                N = num_instances
                arg = np.pi * freqs_k / f0
                with np.errstate(divide='ignore', invalid='ignore'):
                    dkern = np.abs(np.sin(N * arg) / (N * np.sin(arg)))
                dkern[~np.isfinite(dkern)] = 1.0
                ax.plot(
                    freqs_k,
                    spec * dkern,
                    color=color,
                    lw=0.6,
                    alpha=max(0.15, alphas[k] * 0.7),
                    linestyle=':',
                )

            # Analytical structural spectrum overlay
            raw_analytical = rd.get(f'analytical_{ax_name}', [])
            if len(raw_analytical) > 0:
                analytical = np.asarray(raw_analytical, dtype=np.float32)
                ax.plot(
                    freqs_k,
                    analytical,
                    color='k',
                    lw=0.7,
                    alpha=max(0.2, alphas[k] * 0.9),
                    linestyle='--',
                )

        # Forbidden bands — shaded regions (drawn once, shared for all CTRs)
        for band in forbidden_bands:
            ax.axvspan(band[0], band[1], alpha=0.15, color='red', zorder=0)

        # Candidate frequency lines — union over all canonical TRs
        # Violations from any CTR are shown in red; others in dimgray.
        drawn_freqs: dict[float, int] = {}  # freq -> max is_viol flag
        for rd in all_rd:
            cf_arr = np.asarray(rd.get('candidate_freqs', []), dtype=np.float64)
            viol_arr = np.asarray(rd.get('candidate_violations', []), dtype=np.int32)
            for ci in range(len(cf_arr)):
                cf = float(cf_arr[ci])
                if cf < freq_min or cf > freq_max:
                    continue
                iv = int(viol_arr[ci]) if ci < len(viol_arr) else 0
                drawn_freqs[cf] = max(drawn_freqs.get(cf, 0), iv)
        for cf, is_viol in drawn_freqs.items():
            ax.axvline(
                cf,
                color='red' if is_viol else 'dimgray',
                linestyle='-' if is_viol else '--',
                linewidth=1.0 if is_viol else 0.6,
                alpha=0.7 if is_viol else 0.4,
                zorder=1,
            )

        ax.set_xlim(freq_min, freq_max)
        ax.set_ylabel(f'{axis_labels[ax_name]} (a.u.)')
        ax.grid(True, alpha=0.3)
        if n_ctrs > 1 and panel_idx == 0:
            ax.legend(fontsize=7, loc='upper right', framealpha=0.6)

    # --- Annotations on top panel only ---
    ax_top = axes[0]

    # Per-candidate grad-amp labels (mT/m) — worst case across all CTRs,
    # only for candidates that violate a forbidden band in any CTR.
    viol_amps: dict[float, float] = {}
    for rd in all_rd:
        cf_arr = np.asarray(rd.get('candidate_freqs', []), dtype=np.float64)
        viol_arr = np.asarray(rd.get('candidate_violations', []), dtype=np.int32)
        amp_arr = np.asarray(rd.get('candidate_grad_amps', []), dtype=np.float64)
        for ci in range(len(cf_arr)):
            cf = float(cf_arr[ci])
            if cf < freq_min or cf > freq_max:
                continue
            is_viol = int(viol_arr[ci]) if ci < len(viol_arr) else 0
            if is_viol and ci < len(amp_arr) and amp_arr[ci] > 0.0:
                gamp = float(amp_arr[ci]) * hz_per_m_to_mT_per_m
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
            color='red',
            alpha=0.9,
            clip_on=False,
        )

    # Per-forbidden-band: max candidate amplitude within band (mT/m),
    # worst case across all canonical TRs.
    for band in forbidden_bands:
        band_lo, band_hi, band_limit = band[0], band[1], band[2]
        band_limit_mT = band_limit * hz_per_m_to_mT_per_m
        max_cand_mT = 0.0
        for rd in all_rd:
            cf_arr = np.asarray(rd.get('candidate_freqs', []), dtype=np.float64)
            amp_arr = np.asarray(rd.get('candidate_grad_amps', []), dtype=np.float64)
            for ci in range(len(cf_arr)):
                cf = float(cf_arr[ci])
                if band_lo <= cf <= band_hi and ci < len(amp_arr):
                    val = float(amp_arr[ci]) * hz_per_m_to_mT_per_m
                    if val > max_cand_mT:
                        max_cand_mT = val
        band_center = 0.5 * (band_lo + band_hi)
        label = f'limit: {band_limit_mT:.1f} mT/m'
        if max_cand_mT > 0.0:
            label += f'\nmax: {max_cand_mT:.1f} mT/m'
        ax_top.annotate(
            label,
            xy=(band_center, 0.95),
            xycoords=('data', 'axes fraction'),
            ha='center',
            va='top',
            fontsize=7,
            color='red',
            alpha=0.85,
            bbox={'boxstyle': 'round,pad=0.2', 'fc': 'white', 'ec': 'red', 'alpha': 0.6},
        )

    # Title and axis labels
    axes[0].set_title(title)
    axes[2].set_xlabel('Frequency (Hz)')

    # Echo spacing secondary axis on top panel only
    _add_echo_spacing_axis(axes[0], freq_min, freq_max)

    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        fig.tight_layout()


def grad_spectrum(
    seq: SequenceCollection,
    *,
    sequence_idx: int | None = None,
    forbidden_bands: list[tuple[float, float, float]] | None = None,
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
