"""Gradient spectrum visualisation for SequenceCollection (plotting only, no pass/fail)."""

__all__ = ['grad_spectrum']

import warnings

import numpy as np
from matplotlib.colors import Normalize

from ._extension._pulseqlib_wrapper import _calc_acoustic_spectra, _find_tr
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
    thresholds: list[float],
    peak_log10_threshold: float | None,
    peak_norm_scale: float | None,
    peak_eps: float | None,
    peak_prominence: float | None,
) -> None:
    import matplotlib.pyplot as plt

    rd = _calc_acoustic_spectra(
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

    num_windows = rd['num_windows']
    num_freq_bins = rd['num_freq_bins']
    frequencies = rd['freq_min_hz'] + np.arange(num_freq_bins) * rd['freq_spacing_hz']

    spectrograms = {}
    peaks = {}
    for ax_name in ('gx', 'gy', 'gz'):
        spectrograms[ax_name] = np.asarray(
            rd[f'spectrogram_{ax_name}'],
            dtype=np.float32,
        ).reshape(num_windows, num_freq_bins)
        peaks[ax_name] = np.asarray(
            rd[f'peaks_{ax_name}'],
            dtype=np.int32,
        ).reshape(num_windows, num_freq_bins)

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
        # Fix 0/0 at multiples of f0
        dkern[~np.isfinite(dkern)] = 1.0
        for ax_name in ('gx', 'gy', 'gz'):
            dirichlet_env[ax_name] = spectrum_full[ax_name] * dkern

    freq_min = 0.0
    freq_max = float(frequencies[-1])

    axis_labels = {'gx': 'Gx', 'gy': 'Gy', 'gz': 'Gz'}
    colors = {'gx': 'C0', 'gy': 'C1', 'gz': 'C2'}
    title_prefix = f'[SS{subsequence_idx}, CTR{canonical_tr_idx}] '

    use_sliding_windows = num_windows > 1
    if use_sliding_windows:
        fig, axes = plt.subplots(
            2,
            3,
            figsize=(16, 8),
            gridspec_kw={'height_ratios': [1, 1]},
        )

        for col, ax_name in enumerate(('gx', 'gy', 'gz')):
            ax = axes[0, col]
            sg = spectrograms[ax_name]
            pk = peaks[ax_name]

            ax.pcolormesh(
                frequencies,
                np.arange(num_windows),
                sg,
                cmap='viridis',
                shading='auto',
                norm=Normalize(vmin=sg.min(), vmax=sg.max()),
            )

            peak_coords = np.where(pk > 0)
            if len(peak_coords[0]) > 0:
                ax.plot(
                    frequencies[peak_coords[1]],
                    peak_coords[0],
                    marker='*',
                    color='red',
                    linestyle='none',
                    markersize=10,
                    markeredgewidth=0,
                )

            for band in forbidden_bands:
                ax.axvline(band[0], color='white', ls='--', lw=1.2, alpha=0.8)
                ax.axvline(band[1], color='white', ls='--', lw=1.2, alpha=0.8)

            ax.set_xlim(freq_min, freq_max)
            ax.set_ylim(-0.5, num_windows - 0.5)
            ax.set_xlabel('Frequency (Hz)')
            ax.set_ylabel('Window Index')
            ax.set_title(f'{title_prefix}{axis_labels[ax_name]} Spectrogram')
            _add_echo_spacing_axis(ax, freq_min, freq_max)
    else:
        fig, axes = plt.subplots(1, 3, figsize=(16, 4.5))

    for col, ax_name in enumerate(('gx', 'gy', 'gz')):
        ax = axes[1, col] if use_sliding_windows else axes[col]
        color = colors[ax_name]

        ax.plot(
            frequencies,
            spectrum_full[ax_name],
            color=color,
            lw=0.8,
            alpha=0.5,
            label='Canonical TR',
        )

        if ax_name in dirichlet_env:
            ax.plot(
                frequencies,
                dirichlet_env[ax_name],
                color=color,
                lw=0.6,
                alpha=0.35,
                linestyle=':',
                label=f'Dirichlet (N={num_instances})',
            )

        for band in forbidden_bands:
            ax.axvspan(
                band[0],
                band[1],
                alpha=0.15,
                color='red',
                zorder=0,
            )

        # Structural candidate frequency overlays
        cand_freqs = np.asarray(rd.get(f'candidate_freqs_{ax_name}', []), dtype=np.float64)
        cand_viols = np.asarray(rd.get(f'candidate_violations_{ax_name}', []), dtype=np.int32)
        for ci in range(len(cand_freqs)):
            cf = cand_freqs[ci]
            if cf < freq_min or cf > freq_max:
                continue
            is_viol = int(cand_viols[ci]) if ci < len(cand_viols) else 0
            ax.axvline(
                cf,
                color='red' if is_viol else 'green',
                linestyle='-' if is_viol else '--',
                linewidth=1.0 if is_viol else 0.6,
                alpha=0.7 if is_viol else 0.4,
                zorder=1,
            )

        for thr in thresholds:
            ax.axhline(
                thr,
                color='orange',
                linestyle='--',
                linewidth=1.0,
                alpha=0.8,
                label=f'threshold {thr:g}%',
            )

        ax.set_xlim(freq_min, freq_max)
        ax.set_xlabel('Frequency (Hz)')
        ax.set_ylabel('Magnitude (a.u.)')
        ax.set_title(f'{title_prefix}{axis_labels[ax_name]} Harmonic Spectrum')
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
    threshold_percent: float | list[float] | tuple[float, ...] | None = None,
    peak_log10_threshold: float | None = None,
    peak_norm_scale: float | None = None,
    peak_eps: float | None = None,
    peak_prominence: float | None = None,
) -> None:
    """Plot acoustic spectra for gradient waveforms in the canonical TR.

    Creates a two-row figure:

    * **Top row** — three sliding-window spectrograms (Gx, Gy, Gz).
        * **Bottom row** — canonical-TR harmonic spectrum (max envelope across
      windows) with forbidden-band overlays.

        Canonical-TR selection follows the C safety backend per shot-ID
        combination, using shot-filtered ``AMP_MAX_POS`` for each group.

    No pass/fail check is performed — use
    :meth:`SequenceCollection.check` for that.

    Parameters
    ----------
    seq : SequenceCollection
        The sequence to analyse.
    sequence_idx : int
        Subsequence index (0-based, default 0).
    forbidden_bands : list of (freq_min, freq_max, max_amplitude), optional
        Forbidden frequency bands.  Each tuple gives
        ``(freq_min_Hz, freq_max_Hz, max_allowed_amplitude_Hz_per_m)``.
        Drawn as shaded regions on the harmonic plot.
    window_duration : float
        Sliding-window size in seconds (default 25 ms).
    spectral_resolution : float
        Target frequency resolution in Hz (default 5 Hz).
    max_frequency : float
        Upper frequency limit in Hz (default 3000 Hz).
    threshold_percent : float or sequence of float, optional
        Extra horizontal threshold guide(s) drawn on the harmonic plots.
        Accepts a single value or a list/tuple, e.g. ``80.0`` or
        ``[80.0, 100.0]``.
    peak_log10_threshold : float, optional
        Resonance detector threshold in log10 space. Higher values detect
        fewer peaks.
    peak_norm_scale : float, optional
        Normalization scale used before the log transform in resonance
        detection.
    peak_eps : float, optional
        Positive epsilon added before log transform for numerical stability.
    peak_prominence : float, optional
        Minimum prominence (in log10 units) for a peak to be retained.
        Peaks with prominence below this value are discarded.
    """
    if forbidden_bands is None:
        forbidden_bands = []

    if threshold_percent is None:
        thresholds: list[float] = []
    elif isinstance(threshold_percent, (list, tuple)):
        thresholds = [float(v) for v in threshold_percent]
    else:
        thresholds = [float(threshold_percent)]

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
                thresholds=thresholds,
                peak_log10_threshold=peak_log10_threshold,
                peak_norm_scale=peak_norm_scale,
                peak_eps=peak_eps,
                peak_prominence=peak_prominence,
            )
