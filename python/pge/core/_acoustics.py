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
    cand_grad_amps = np.asarray(rd.get('candidate_grad_amps', []), dtype=np.float64)

    fig, ax = plt.subplots(1, 1, figsize=(12, 5))

    for ax_name in ('gx', 'gy', 'gz'):
        color = colors[ax_name]
        label = axis_labels[ax_name]

        ax.plot(
            frequencies,
            spectrum_full[ax_name],
            color=color,
            lw=0.8,
            alpha=0.5,
            label=label,
        )

        if ax_name in dirichlet_env:
            ax.plot(
                frequencies,
                dirichlet_env[ax_name],
                color=color,
                lw=0.6,
                alpha=0.35,
                linestyle=':',
            )

    # Forbidden bands — shaded with max-allowed label
    y_top = ax.get_ylim()[1]
    for band in forbidden_bands:
        ax.axvspan(
            band[0],
            band[1],
            alpha=0.15,
            color='red',
            zorder=0,
        )
        max_allowed_mT = band[2] * hz_per_m_to_mT_per_m
        band_center = 0.5 * (band[0] + band[1])
        ax.annotate(
            f'{max_allowed_mT:.2f} mT/m',
            xy=(band_center, 1.0),
            xycoords=('data', 'axes fraction'),
            ha='center',
            va='bottom',
            fontsize=7,
            color='red',
            alpha=0.85,
            rotation=90,
        )

    # Candidate frequency lines with max grad-amp annotation
    cand_freqs_gx = np.asarray(rd.get('candidate_freqs_gx', []), dtype=np.float64)
    cand_viols_gx = np.asarray(rd.get('candidate_violations_gx', []), dtype=np.int32)
    for ci in range(len(cand_freqs_gx)):
        cf = cand_freqs_gx[ci]
        if cf < freq_min or cf > freq_max:
            continue
        is_viol = int(cand_viols_gx[ci]) if ci < len(cand_viols_gx) else 0
        ax.axvline(
            cf,
            color='red' if is_viol else 'dimgray',
            linestyle='-' if is_viol else '--',
            linewidth=1.0 if is_viol else 0.6,
            alpha=0.7 if is_viol else 0.4,
            zorder=1,
        )
        if ci < len(cand_grad_amps) and cand_grad_amps[ci] > 0.0:
            gamp_mT = float(cand_grad_amps[ci]) * hz_per_m_to_mT_per_m
            ax.annotate(
                f'{gamp_mT:.2f}',
                xy=(cf, 1.0),
                xycoords=('data', 'axes fraction'),
                ha='left',
                va='top',
                fontsize=6,
                color='red' if is_viol else 'dimgray',
                alpha=0.9,
                rotation=90,
                clip_on=True,
            )

    ax.set_xlim(freq_min, freq_max)
    ax.set_xlabel('Frequency (Hz)')
    ax.set_ylabel('Magnitude (a.u.)')
    ax.set_title(title)
    ax.legend(fontsize=8, loc='upper right')
    ax.grid(True, alpha=0.3)
    _add_echo_spacing_axis(ax, freq_min, freq_max)

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

