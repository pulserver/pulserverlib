"""Pure-Python reader for pulseqlib binary cache sections and Python dataclasses.

Provides:
  * Dataclasses for SequenceDescriptionInfo and TrajectoryInfo.
  * ``build_sequence_description_info`` — build SequenceDescriptionInfo from
    the pybind ``_get_sequence_parameters`` / ``_get_sequence_description``
    output dicts (always available).
  * ``read_trajectory_info`` — parse section 5 (TRAJECTORY) from a co-located
    .bin cache file (present when the PSD predownload phase has been run).
"""

from __future__ import annotations

__all__ = [
    'EncodingSpace',
    'LabelLimits',
    'SeqRow',
    'SequenceDescription',
    'SequenceDescriptionInfo',
    'SequenceParameters',
    'TrajTableEntry',
    'TrajectoryInfo',
    'build_sequence_description_info',
    'read_trajectory_info',
]

import struct
from dataclasses import dataclass, field
from pathlib import Path

# ── Constants (mirror pulseqlib C sources) ────────────────────────────

_CACHE_ENDIAN_MARKER    = 0x01020304
_SECTION_TRAJECTORY     = 4
_SEQ_EVENT_TYPE_WAIT    = 0
_SEQ_EVENT_TYPE_RF      = 1
_SEQ_EVENT_TYPE_ADC     = 2

# ── Dataclasses — SequenceDescription (built via pybind) ─────────────

@dataclass
class SequenceParameters:
    """Scan-global sequence parameters (mirrors pulseqlib_sequence_parameters)."""
    min_te_us:          float
    min_tr_us:          float
    max_tr_us:          float
    max_flip_angle_deg: float
    total_scan_time_us: float
    num_subseqs:        int


@dataclass
class SeqRow:
    """Single pass-block row in the compact sequence description.

    Param layout (mirrors pulseqlib_seq_event):

    RF    — params[0]=rf_def_id, [1]=rf_use, [2]=act_amplitude_hz,
             [3]=phase_offset_rad, [4]=freq_offset_hz, [5]=rf_shim_id
    ADC   — params[0]=adc_role, [1]=phase_offset_rad, [2..5]=0
    OTHER — params[0..5]=0
    """
    type:         int          # 0=OTHER, 1=RF, 2=ADC
    timestamp_us: float        # pass-relative anchor time (us)
    params:       list[float]  # 6 floats, type-specific

    @property
    def is_other(self) -> bool:
        return self.type == _SEQ_EVENT_TYPE_WAIT  # alias for OTHER

    @property
    def is_rf(self) -> bool:
        return self.type == _SEQ_EVENT_TYPE_RF

    @property
    def is_adc(self) -> bool:
        return self.type == _SEQ_EVENT_TYPE_ADC

    # ── RF accessors ─────────────────────────────────────────
    @property
    def rf_def_id(self) -> int:
        return int(self.params[0])

    @property
    def rf_use(self) -> int:
        return int(self.params[1])

    @property
    def rf_act_amplitude_hz(self) -> float:
        return self.params[2]

    @property
    def rf_phase_offset_rad(self) -> float:
        return self.params[3]

    @property
    def rf_freq_offset_hz(self) -> float:
        return self.params[4]

    @property
    def rf_shim_id(self) -> int:
        return int(self.params[5])

    # ── ADC accessors ────────────────────────────────────────
    @property
    def adc_role(self) -> int:
        return int(self.params[0])

    @property
    def adc_phase_offset_rad(self) -> float:
        return self.params[1]


@dataclass
class SequenceDescription:
    """Per-subsequence compact row table (one row per pass block)."""
    subseq_idx:    int
    tr_duration_us: float
    rows:          list[SeqRow]

    @property
    def num_rf_rows(self) -> int:
        return sum(1 for r in self.rows if r.is_rf)

    @property
    def num_adc_rows(self) -> int:
        return sum(1 for r in self.rows if r.is_adc)

    @property
    def flip_angles_deg(self) -> list[float]:
        """Nominal flip angles (degrees) — requires rf_stats lookup; placeholder."""
        return []


@dataclass
class SequenceDescriptionInfo:
    """Full sequence description: global parameters + per-subsequence data."""
    seq_params: SequenceParameters
    subseqs:    list[SequenceDescription]


# ── Dataclasses — TrajectoryInfo (parsed from .bin section 5) ────────

@dataclass
class LabelLimits:
    """Min/max label range for one dimension."""
    min: int
    max: int


@dataclass
class EncodingSpace:
    """One encoding space (mirrors mrdserver::EncodingSpace)."""
    fov:               tuple[float, float, float]
    matrix:            tuple[float, float, float]
    nav_fov:           tuple[float, float, float]
    nav_matrix:        tuple[float, float, float]
    subseq_idx:        int
    nav_subseq_offset: int
    label_limits:      dict[str, LabelLimits]  # keys: slc phs rep avg seg set eco par lin acq


@dataclass
class TrajTableEntry:
    """One ADC event in the trajectory table."""
    kx_shot_id:         int
    ky_shot_id:         int
    kz_shot_id:         int
    gx_amplitude:       float
    gy_amplitude:       float
    gz_amplitude:       float
    rotation_id:        int
    slc: int; seg: int; rep: int; avg: int
    set: int; eco: int; phs: int; lin: int; par: int; acq: int
    flags:              int   # 64-bit bitmask stored as Python int
    center_sample:      int
    sample_time_us:     float
    encoding_space_ref: int
    off:                int = 0  # Pulseq LABELSET OFF flag (1 = discard)


@dataclass
class TrajectoryInfo:
    """Parsed trajectory data from cache section 4."""
    kshots:          list[list[float]]      # [num_shots][num_samples]
    encoding_spaces: list[EncodingSpace]
    table:           list[TrajTableEntry]


# ── Builder from pybind dicts ─────────────────────────────────────────

def build_sequence_description_info(
    params_dict: dict,
    desc_dicts: list[dict],
) -> SequenceDescriptionInfo:
    """Build a ``SequenceDescriptionInfo`` from pybind output dicts.

    Parameters
    ----------
    params_dict : dict
        Output of ``_get_sequence_parameters(cseq)``.
    desc_dicts : list of dict
        Output of ``_get_sequence_description(cseq, i)`` for each subsequence.
    """
    sp = SequenceParameters(
        min_te_us          = float(params_dict['min_te_us']),
        min_tr_us          = float(params_dict['min_tr_us']),
        max_tr_us          = float(params_dict['max_tr_us']),
        max_flip_angle_deg = float(params_dict['max_flip_angle_deg']),
        total_scan_time_us = float(params_dict['total_scan_time_us']),
        num_subseqs        = int(params_dict['num_subseqs']),
    )

    subseqs = []
    for d in desc_dicts:
        rows = [
            SeqRow(
                type         = int(r['type']),
                timestamp_us = float(r['timestamp_us']),
                params       = [float(v) for v in r['params']],
            )
            for r in d['rows']
        ]
        subseqs.append(SequenceDescription(
            subseq_idx     = int(d['subseq_idx']),
            tr_duration_us = float(d['tr_duration_us']),
            rows           = rows,
        ))

    return SequenceDescriptionInfo(seq_params=sp, subseqs=subseqs)


# ── Pure-Python .bin reader for section 4 (TRAJECTORY) ───────────────

_LL_NAMES = ('slc', 'phs', 'rep', 'avg', 'seg', 'set', 'eco', 'par', 'lin', 'acq')


def _find_section(data: bytes, section_id: int):
    """Scan the cache file header; return (offset, size, do_swap) or None."""
    if len(data) < 28:
        return None

    marker_le = struct.unpack_from('<I', data, 0)[0]
    do_swap = False
    if marker_le != _CACHE_ENDIAN_MARKER:
        marker_be = struct.unpack_from('>I', data, 0)[0]
        if marker_be != _CACHE_ENDIAN_MARKER:
            return None
        do_swap = True

    endian = '>' if do_swap else '<'
    # offset 4: version_major, version_minor, vendor, stored_size, num_sections (5 × int32)
    num_sections = struct.unpack_from(f'{endian}5i', data, 4)[4]
    if not 1 <= num_sections <= 16:
        return None

    hdr_end = 4 + 5 * 4  # 24 bytes
    for i in range(num_sections):
        sid, soff, ssz = struct.unpack_from(f'{endian}3i', data, hdr_end + i * 12)
        if sid == section_id:
            return (soff, ssz, do_swap)

    return None


def read_trajectory_info(seq_path: str | Path) -> TrajectoryInfo | None:
    """Read trajectory data (cache section 4) from a co-located ``.bin`` file.

    Returns ``None`` when the ``.bin`` file does not exist or section 4 is
    absent (e.g. for Cartesian sequences or when the cache was not written
    by the PSD predownload phase).

    Parameters
    ----------
    seq_path : str or Path
        Path to the ``.seq`` sequence file.  The function looks for a file
        with the same base name but ``.bin`` extension.
    """
    bin_path = Path(seq_path).with_suffix('.bin')
    if not bin_path.exists():
        return None

    try:
        data = bin_path.read_bytes()
    except OSError:
        return None

    found = _find_section(data, _SECTION_TRAJECTORY)
    if found is None:
        return None

    offset, _size, do_swap = found
    endian = '>' if do_swap else '<'
    fmt_i = f'{endian}i'
    fmt_f = f'{endian}f'
    pos = offset

    def _ri() -> int:
        nonlocal pos
        v = struct.unpack_from(fmt_i, data, pos)[0]; pos += 4
        return v

    def _rf() -> float:
        nonlocal pos
        v = struct.unpack_from(fmt_f, data, pos)[0]; pos += 4
        return v

    def _ri_n(n: int) -> tuple:
        nonlocal pos
        v = struct.unpack_from(f'{endian}{n}i', data, pos); pos += 4 * n
        return v

    def _rf_n(n: int) -> tuple:
        nonlocal pos
        v = struct.unpack_from(f'{endian}{n}f', data, pos); pos += 4 * n
        return v

    try:
        # ── kshot library ────────────────────────────────────
        num_shots = _ri()
        kshots: list[list[float]] = []
        for _ in range(num_shots):
            ns = _ri()
            k = list(_rf_n(ns)) if ns > 0 else []
            kshots.append(k)

        # ── encoding spaces ──────────────────────────────────
        num_es = _ri()
        encoding_spaces: list[EncodingSpace] = []
        for _ in range(num_es):
            fov        = _rf_n(3)
            matrix     = _rf_n(3)
            nav_fov    = _rf_n(3)
            nav_matrix = _rf_n(3)
            subseq_idx        = _ri()
            nav_subseq_offset = _ri()
            # 10 × {min, max} ints  =  20 ints
            ll_raw = _ri_n(20)
            label_limits = {
                _LL_NAMES[i]: LabelLimits(ll_raw[i * 2], ll_raw[i * 2 + 1])
                for i in range(10)
            }
            encoding_spaces.append(EncodingSpace(
                fov=fov, matrix=matrix, nav_fov=nav_fov, nav_matrix=nav_matrix,
                subseq_idx=subseq_idx, nav_subseq_offset=nav_subseq_offset,
                label_limits=label_limits,
            ))

        # ── trajectory table ─────────────────────────────────
        num_entries = _ri()
        table: list[TrajTableEntry] = []
        for _ in range(num_entries):
            kx_id = _ri(); ky_id = _ri(); kz_id = _ri()
            gx = _rf();  gy = _rf();  gz = _rf()
            rot_id = _ri()
            slc, seg, rep, avg, set_, eco, phs, lin, par, acq = _ri_n(10)
            flags_lo = _ri(); flags_hi = _ri()
            flags = (flags_hi << 32) | (flags_lo & 0xFFFFFFFF)
            center_sample = _ri()
            sample_time_us = _rf()
            encoding_space_ref = _ri()
            off = _ri()
            table.append(TrajTableEntry(
                kx_shot_id=kx_id, ky_shot_id=ky_id, kz_shot_id=kz_id,
                gx_amplitude=gx, gy_amplitude=gy, gz_amplitude=gz,
                rotation_id=rot_id,
                slc=slc, seg=seg, rep=rep, avg=avg,
                set=set_, eco=eco, phs=phs, lin=lin, par=par, acq=acq,
                flags=flags,
                center_sample=center_sample,
                sample_time_us=sample_time_us,
                encoding_space_ref=encoding_space_ref,
                off=off,
            ))

        return TrajectoryInfo(kshots=kshots, encoding_spaces=encoding_spaces, table=table)

    except struct.error:
        # Malformed / truncated section — degrade gracefully
        return None
