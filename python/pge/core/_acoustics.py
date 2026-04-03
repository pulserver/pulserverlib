"""Gradient spectrum (mechanical resonances) visualisation for SequenceCollection (plotting only, no pass/fail)."""

__all__ = ['grad_spectrum']

import warnings

import numpy as np

from ._extension._pulseqlib_wrapper import _calc_mech_resonances, _find_tr
from ._helpers import _add_echo_spacing_axis
from ._sequence import SequenceCollection


def _plot_grad_spectrum_single(
    seq: SequenceCollection,
    *,
    subsequence_idx: int,
    canonical_tr_idx: int,
    forbidden_bands: list[tuple[float, float, float]],
    target_window_size: int,
    spectral_resolution: float,
    max_frequency: float,
    peak_log10_threshold: float | None,
    peak_norm_scale: float | None,
    peak_eps: float | None,
    peak_prominence: float | None,
) -> None:
    import matplotlib.pyplot as plt

    gamma_hz_per_t = float(seq.system.gamma)
    hz_per_m_to_mT_per_m = 1e3 / gamma_hz_per_t

    rd = _calc_mech_resonances(
        seq._cseq,
        subsequence_idx=subsequence_idx,
        canonical_tr_idx=canonical_tr_idx,
        target_window_size=target_window_size,
        target_resolution_hz=spectral_resolution,
        max_freq_hz=max_frequency,
        forbidden_bands=forbidden_bands,
        peak_log10_threshold=peak_log10_threshold,
        peak_norm_scale=peak_norm_scale,
        peak_eps=peak_eps,
        peak_prominence=peak_prominence,
    )

    num_freq_bins = rd['num_freq_bins']
    frequencies = rd['freq_min_hz'] + np.arange(num_freq_bins) * rd['freq_spacing_hz']

    spectrum_full = {}
    for ax_name in ('gx', 'gy', 'gz'):
        spectrum_full[ax_name] = np.asarray(
            rd[f'spectrum_full_{ax_name}'],
            dtype=np.float32,
        )

    # Dirichlet kernel envelope for multi-TR visualization
    num_instances = int(rd.get('num_instances', 0))
    f0 = float(rd.get('freq_spacing_seq_hz', 0.0))
    dirichlet_env = {}
    if num_instances > 1 and f0 > 0.0:
        N = num_instances
        arg = np.pi * frequencies / f0
        with np.errstate(divide='ignore', invalid='ignore'):
            dkern = np.abs(np.sin(N * arg) / (N * np.sin(arg)))
        # Replace 0/0 at exact multiples of f0 with the correct limit (1.0)
        dkern[~np.isfinite(dkern)] = 1.0
        for ax_name in ('gx', 'gy', 'gz'):
            dirichlet_env[ax_name] = spectrum_full[ax_name] * dkern

    freq_min = 0.0
    freq_max = float(frequencies[-1])

    axis_labels = {'gx': 'Gx', 'gy': 'Gy', 'gz': 'Gz'}
    colors = {'gx': 'C0', 'gy': 'C1', 'gz': 'C2'}
    title = f'[SS{subsequence_idx}, CTR{canonical_tr_idx}] Mechanical Resonances'

    # Candidate data (shared across axes)
    cand_freqs = np.asarray(rd.get('candidate_freqs', []), dtype=np.float64)
    cand_viols = np.asarray(rd.get('candidate_violations', []), dtype=np.int32)
    cand_grad_amps = np.asarray(rd.get('candidate_grad_amps', []), dtype=np.float64)

    # Analytical structural spectrum (dense grid, one per axis)
    analytical = {}
    for ax_name in ('gx', 'gy', 'gz'):
        raw = rd.get(f'analytical_{ax_name}', [])
        analytical[ax_name] = (
            np.asarray(raw, dtype=np.float32) if len(raw) > 0 else None
        )

    fig, axes = plt.subplots(3, 1, figsize=(12, 8), sharex=True)

    axis_names = ('gx', 'gy', 'gz')
    axis_labels = {'gx': 'Gx', 'gy': 'Gy', 'gz': 'Gz'}
    colors = {'gx': 'C0', 'gy': 'C1', 'gz': 'C2'}

    for panel_idx, ax_name in enumerate(axis_names):
        ax = axes[panel_idx]
        color = colors[ax_name]

        # FFT spectrum
        ax.plot(
            frequencies,
            spectrum_full[ax_name],
            color=color,
            lw=0.8,
            alpha=0.5,
        )

        # Dirichlet kernel envelope
        if ax_name in dirichlet_env:
            ax.plot(
                frequencies,
                dirichlet_env[ax_name],
                color=color,
                lw=0.6,
                alpha=0.35,
                linestyle=':',
            )

        # Analytical structural spectrum overlay
        if analytical.get(ax_name) is not None:
            ax.plot(
                frequencies,
                analytical[ax_name],
                color='k',
                lw=0.7,
                alpha=0.45,
            )

        # Forbidden bands — shaded regions
        for band in forbidden_bands:
            ax.axvspan(band[0], band[1], alpha=0.15, color='red', zorder=0)

        # Candidate frequency lines
        for ci in range(len(cand_freqs)):
            cf = cand_freqs[ci]
            if cf < freq_min or cf > freq_max:
                continue
            is_viol = int(cand_viols[ci]) if ci < len(cand_viols) else 0
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

    # --- Annotations on top panel only ---

    # Per-candidate grad-amp labels (horizontal, mT/m, on top panel)
    # Only show annotation for candidates that violate a forbidden band.
    ax_top = axes[0]
    for ci in range(len(cand_freqs)):
        cf = cand_freqs[ci]
        if cf < freq_min or cf > freq_max:
            continue
        is_viol = int(cand_viols[ci]) if ci < len(cand_viols) else 0
        if is_viol and ci < len(cand_grad_amps) and cand_grad_amps[ci] > 0.0:
            gamp_mT = float(cand_grad_amps[ci]) * hz_per_m_to_mT_per_m
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

    # Per-forbidden-band: max candidate amplitude within band (mT/m)
    for band in forbidden_bands:
        band_lo, band_hi, band_limit = band[0], band[1], band[2]
        band_limit_mT = band_limit * hz_per_m_to_mT_per_m
        # Find max candidate grad_amp within this band
        max_cand_mT = 0.0
        for ci in range(len(cand_freqs)):
            cf = cand_freqs[ci]
            if band_lo <= cf <= band_hi and ci < len(cand_grad_amps):
                val = float(cand_grad_amps[ci]) * hz_per_m_to_mT_per_m
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
        for canonical_tr_idx in range(num_canonical):
            _plot_grad_spectrum_single(
                seq,
                subsequence_idx=ss_idx,
                canonical_tr_idx=canonical_tr_idx,
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
