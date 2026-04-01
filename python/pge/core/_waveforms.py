def calc_acoustic_spectra(
    seq: 'SequenceCollection',
    subsequence_idx: int = 0,
    canonical_tr_idx: int = 0,
    target_window_size: int = 0,
    target_resolution_hz: float = 0.0,
    max_freq_hz: float = 0.0,
    forbidden_bands=None,
) -> dict:
    """Compute acoustic spectra for a specific canonical TR of a subsequence."""
    if forbidden_bands is None:
        forbidden_bands = []
    return seq._cseq._calc_acoustic_spectra(
        subsequence_idx,
        canonical_tr_idx,
        target_window_size,
        target_resolution_hz,
        max_freq_hz,
        forbidden_bands,
    )

def calc_pns(
    seq: 'SequenceCollection',
    subsequence_idx: int = 0,
    canonical_tr_idx: int = 0,
    chronaxie_us: float = 360.0,
    rheobase: float = 20.0,
    alpha: float = 0.333,
) -> dict:
    """Compute PNS slew-rate waveforms for a specific canonical TR of a subsequence."""
    return seq._cseq._calc_pns(
        subsequence_idx,
        canonical_tr_idx,
        chronaxie_us,
        rheobase,
        alpha,
    )
"""Native-timing TR waveform extraction for SequenceCollection."""

__all__ = ['TrWaveforms', 'get_tr_waveforms']

from dataclasses import dataclass, field
from typing import Literal

import numpy as np


# Import SequenceCollection before any type annotations use it
from ._sequence import SequenceCollection
from ._extension._pulseqlib_wrapper import _find_tr, _get_tr_waveforms

# Physical-unit conversion constants
_GAMMA_DEFAULT = 42.576e6  # Hz/T


@dataclass
class ChannelWaveform:
    """Single-channel waveform with native (non-uniform) timing.

    Attributes
    ----------
    time_us : np.ndarray
        Time points in microseconds, shape ``(N,)``.
    amplitude : np.ndarray
        Amplitude values, shape ``(N,)``.  Units depend on channel:
        mT/m for gradients, µT for RF magnitude, rad for RF phase.
    """

    time_us: np.ndarray = field(default_factory=lambda: np.empty(0, dtype=np.float32))
    amplitude: np.ndarray = field(default_factory=lambda: np.empty(0, dtype=np.float32))


@dataclass
class AdcEvent:
    """ADC event descriptor."""

    onset_us: float = 0.0
    duration_us: float = 0.0
    num_samples: int = 0
    freq_offset_hz: float = 0.0
    phase_offset_rad: float = 0.0


@dataclass
class BlockDescriptor:
    """Per-block metadata within a TR."""

    start_us: float = 0.0
    duration_us: float = 0.0
    segment_idx: int = -1
    rf_freq_offset_hz: float = 0.0
    rf_phase_offset_rad: float = 0.0
    adc_freq_offset_hz: float = 0.0
    adc_phase_offset_rad: float = 0.0
    rotation_matrix: list = field(default_factory=list)


@dataclass
class TrWaveforms:
    """Complete native-timing TR waveforms for plotting.

    Gradient amplitudes are in mT/m,  RF magnitude in µT,
    RF phase in radians.  Times are in microseconds.

    Attributes
    ----------
    gx, gy, gz : ChannelWaveform
        Gradient waveforms (mT/m).
    rf_mag : ChannelWaveform
        RF magnitude envelope (µT).
    rf_phase : ChannelWaveform
        RF phase (rad).
    adc_events : list[AdcEvent]
        ADC event descriptors.
    blocks : list[BlockDescriptor]
        Per-block metadata (timing, segment assignment).
    total_duration_us : float
        Total duration of the extracted block range (µs).
    gamma_hz_per_t : float
        Gyromagnetic ratio used for unit conversion.
    tr_duration_us : float
        Duration of a single main TR (µs), for pypulseq overlay.
    prep_duration_us : float
        Total duration of prep blocks (µs), for pypulseq overlay.
    """

    gx: ChannelWaveform = field(default_factory=ChannelWaveform)
    gy: ChannelWaveform = field(default_factory=ChannelWaveform)
    gz: ChannelWaveform = field(default_factory=ChannelWaveform)
    rf_mag: ChannelWaveform = field(default_factory=ChannelWaveform)
    rf_phase: ChannelWaveform = field(default_factory=ChannelWaveform)
    adc_events: list = field(default_factory=list)
    blocks: list = field(default_factory=list)
    total_duration_us: float = 0.0
    gamma_hz_per_t: float = _GAMMA_DEFAULT
    tr_duration_us: float = 0.0
    prep_duration_us: float = 0.0


# Amplitude mode constants (must match C defines)


# Amplitude mode constants (must match C defines)
AMP_MAX_POS = 0
AMP_ZERO_VAR = 1
AMP_ACTUAL = 2


def get_tr_waveforms(
    seq: SequenceCollection,
    subsequence_idx: int = 0,
    amplitude_mode: Literal['max_pos', 'zero_var', 'actual'] = 'max_pos',
    tr_index: int = 0,
    collapse_delays: bool = False,
) -> TrWaveforms:
    """Extract native-timing TR waveforms from the segmented representation.

    Returns gradient, RF, and ADC waveforms with per-channel native
    timing (NOT interpolated to a uniform raster).  Amplitudes are
    converted to physical units: mT/m for gradients, µT for RF.

    Parameters
    ----------
    seq : SequenceCollection
        The sequence to analyse.
    subsequence_idx : int
        Subsequence index (0-based, default 0).
    amplitude_mode : {'max_pos', 'zero_var', 'actual'}
        How to resolve multi-shot gradient amplitudes:

        - ``'max_pos'`` — position-max across all TRs (safety view).
        - ``'zero_var'`` — zero variable grads, keep constant (k-space view).
        - ``'actual'``  — signed amplitude for a specific TR instance.
    tr_index : int
        TR instance index (0-based).  Only used for ``'actual'`` mode.
    (Preparation and cooldown blocks are always included if present.)
    collapse_delays : bool
        Shrink pure-delay blocks (no RF, no grad, no ADC) to 0.1 ms
        at the C level, producing a compact waveform timeline.

    Returns
    -------
    TrWaveforms
        Dataclass with per-channel waveforms and metadata.
    """
    # Accept 'canonical' as alias for 'max_pos' (structural canonical TR)

    if amplitude_mode == 'canonical':
        amplitude_mode = 'max_pos'
    mode_map = {'max_pos': AMP_MAX_POS, 'zero_var': AMP_ZERO_VAR, 'actual': AMP_ACTUAL}
    c_mode = mode_map[amplitude_mode]

    gamma = seq.system.gamma  # Hz/T

    # Call C extension
    raw = _get_tr_waveforms(
        seq._cseq,
        subsequence_idx,
        c_mode,
        tr_index,
        collapse_delays,
    )

    # Unpack raw waveform arrays
    hz_per_m_to_mT_per_m = 1e3 / gamma
    gx = ChannelWaveform(np.array(raw['gx']['time_us']), np.array(raw['gx']['amplitude']) * hz_per_m_to_mT_per_m)
    gy = ChannelWaveform(np.array(raw['gy']['time_us']), np.array(raw['gy']['amplitude']) * hz_per_m_to_mT_per_m)
    gz = ChannelWaveform(np.array(raw['gz']['time_us']), np.array(raw['gz']['amplitude']) * hz_per_m_to_mT_per_m)
    # Convert RF magnitude from Hz to µT for validation (µT = 1e6 / gamma * Hz).
    # The C backend stores signed magnitude (negative for negative lobes) with
    # phase separate.  Convert to unsigned magnitude + adjusted phase so the
    # phase of negative lobes is π, matching the pypulseq convention.
    rf_mag_hz = np.array(raw['rf_mag']['amplitude'])
    rf_phase_raw = np.array(raw['rf_phase']['amplitude'])
    neg = rf_mag_hz < 0
    rf_mag_hz = np.abs(rf_mag_hz)
    rf_phase_raw = rf_phase_raw.copy()
    rf_phase_raw[neg] += np.pi
    rf_mag_uT = rf_mag_hz * (1e6 / gamma)
    rf_mag = ChannelWaveform(np.array(raw['rf_mag']['time_us']), rf_mag_uT)
    rf_phase = ChannelWaveform(np.array(raw['rf_phase']['time_us']), rf_phase_raw)
    adc_events = [AdcEvent(**adc) for adc in raw['adc_events']]
    # Extract per-block scan table parameters if present
    blocks = []
    for blk in raw['blocks']:
        block_kwargs = dict(
            start_us=blk.get('start_us', 0.0),
            duration_us=blk.get('duration_us', 0.0),
            segment_idx=blk.get('segment_idx', -1),
            rf_freq_offset_hz=blk.get('rf_freq_offset_hz', 0.0),
            rf_phase_offset_rad=blk.get('rf_phase_offset_rad', 0.0),
            adc_freq_offset_hz=blk.get('adc_freq_offset_hz', 0.0),
            adc_phase_offset_rad=blk.get('adc_phase_offset_rad', 0.0),
            rotation_matrix=blk.get('rotation_matrix', []),
        )
        blocks.append(BlockDescriptor(**block_kwargs))

    # Get TR timing info for pypulseq overlay
    tr_info = _find_tr(seq._cseq, subsequence_idx=subsequence_idx)
    tr_dur = tr_info['tr_duration_us']
    prep_dur = sum(
        b['duration_us']
        for b in raw['blocks']
        if b['segment_idx'] < 0  # prep/cooldown heuristic — refine if needed
    )

    return TrWaveforms(
        gx=gx,
        gy=gy,
        gz=gz,
        rf_mag=rf_mag,
        rf_phase=rf_phase,
        adc_events=adc_events,
        blocks=blocks,
        total_duration_us=raw['total_duration_us'],
        gamma_hz_per_t=gamma,
        tr_duration_us=tr_dur,
        prep_duration_us=prep_dur,
    )
