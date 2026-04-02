"""Tests for SequenceCollection construction and basic accessors."""

import numpy as np
import pypulseq as pp
import pytest

from pge import Opts, SequenceCollection


def _safe_delay_collection() -> SequenceCollection:
    """Return a trivially safe collection for check() smoke tests."""
    seq = pp.Sequence(
        system=pp.Opts(
            max_grad=200,
            grad_unit="mT/m",
            max_slew=10000,
            slew_unit="T/m/s",
        )
    )
    seq.add_block(pp.make_delay(1e-3))
    seq.add_block(pp.make_delay(1e-3))
    return SequenceCollection(seq)


def test_from_sequence_object(simple_gre_seq):
    sc = SequenceCollection(simple_gre_seq)
    assert sc.num_sequences == 1


def test_from_sequence_list(simple_gre_seq):
    sc = SequenceCollection([simple_gre_seq])
    assert sc.num_sequences == 1


def test_from_file(simple_gre_seq, tmp_path):
    seq_file = tmp_path / "test.seq"
    simple_gre_seq.write(str(seq_file))
    sc = SequenceCollection(str(seq_file))
    assert sc.num_sequences == 1


def test_invalid_type_raises():
    with pytest.raises(TypeError):
        SequenceCollection(42)


def test_num_blocks(simple_gre_seq):
    sc = SequenceCollection(simple_gre_seq)
    n = sc.num_blocks(0)
    assert n > 0


def test_get_block(simple_gre_seq):
    sc = SequenceCollection(simple_gre_seq)
    n = sc.num_blocks(0)

    block = None
    for idx in range(n):
        try:
            block = sc.get_block(0, idx)
            break
        except KeyError:
            continue

    assert block is not None


def test_get_sequence_copy(simple_gre_seq):
    sc = SequenceCollection(simple_gre_seq)
    seq_copy = sc.get_sequence(0)
    assert isinstance(seq_copy, pp.Sequence)


def test_get_sequence_out_of_range(simple_gre_seq):
    sc = SequenceCollection(simple_gre_seq)
    with pytest.raises(IndexError):
        sc.get_sequence(99)


def test_num_segments(simple_gre_seq):
    sc = SequenceCollection(simple_gre_seq)
    ns = sc.num_segments(0)
    assert ns >= 1


def test_tr_size(simple_gre_seq):
    sc = SequenceCollection(simple_gre_seq)
    ts = sc.tr_size(0)
    assert ts >= 1


def test_segment_size(simple_gre_seq):
    sc = SequenceCollection(simple_gre_seq)
    ns = sc.num_segments(0)
    for i in range(ns):
        sz = sc.segment_size(0, i)
        assert sz >= 1


def test_report_returns_list(simple_gre_seq):
    sc = SequenceCollection(simple_gre_seq)
    reps = sc.report()
    assert isinstance(reps, list)
    assert len(reps) >= 1


def test_report_fields(simple_gre_seq):
    sc = SequenceCollection(simple_gre_seq)
    reps = sc.report()
    r = reps[0]
    assert hasattr(r, "num_blocks")
    assert hasattr(r, "segments")
    assert hasattr(r, "tr_size")
    assert hasattr(r, "tr_duration_s")
    assert r.tr_size >= 1
    assert r.tr_duration_s > 0


def test_report_print(simple_gre_seq):
    sc = SequenceCollection(simple_gre_seq)
    text = sc.report(do_print=True)
    assert isinstance(text, str)
    assert "Subsequence" in text


def test_check_passes():
    sc = _safe_delay_collection()
    sc.check()


def test_check_with_pns():
    sc = _safe_delay_collection()
    sc.check(
        stim_threshold=1e9,
        decay_constant_us=360.0,
        pns_threshold_percent=1000.0,
    )


def test_sequence_collection_opts_override_system(simple_gre_seq):
    opts = Opts(
        gamma=float(simple_gre_seq.system.gamma),
        B0=float(simple_gre_seq.system.B0),
        max_grad=40.0,
        max_slew=150.0,
        b1_max_uT=10.0,
        rf_raster_time=float(simple_gre_seq.system.rf_raster_time),
        grad_raster_time=float(simple_gre_seq.system.grad_raster_time),
        adc_raster_time=float(simple_gre_seq.system.adc_raster_time),
        block_duration_raster=float(simple_gre_seq.system.block_duration_raster),
    )
    sc = SequenceCollection(simple_gre_seq, system=opts)
    assert sc.system is opts


def test_geopts_from_coil_model_defaults():
    opts = Opts.from_coil_model(
        'hrmw',
        gamma=42.576e6,
        B0=3.0,
        b1_max_uT=10.0,
    )
    assert np.isclose(opts.grad_raster_time, 4e-6)
    assert np.isclose(opts.rf_raster_time, 2e-6)
    assert np.isclose(opts.adc_raster_time, 2e-6)
    assert np.isclose(opts.chronaxie_us, 642.4)
    assert np.isclose(opts.default_stim_threshold(), 17.9 / 0.310)
    assert len(opts.forbidden_bands) > 0


def test_check_uses_opts_forbidden_bands_default(simple_gre_seq, monkeypatch):
    opts = Opts(
        gamma=float(simple_gre_seq.system.gamma),
        B0=float(simple_gre_seq.system.B0),
        max_grad=40.0,
        max_slew=150.0,
        forbidden_bands=[(100.0, 200.0, 1.0)],
        rf_raster_time=float(simple_gre_seq.system.rf_raster_time),
        grad_raster_time=float(simple_gre_seq.system.grad_raster_time),
        adc_raster_time=float(simple_gre_seq.system.adc_raster_time),
        block_duration_raster=float(simple_gre_seq.system.block_duration_raster),
    )
    sc = SequenceCollection(simple_gre_seq, system=opts)

    captured = {}

    def _fake_consistency(_cseq):
        return None

    def _fake_safety(_cseq, **kwargs):
        captured.update(kwargs)

    import pge.core._sequence as _sequence_mod

    monkeypatch.setattr(_sequence_mod, '_check_consistency', _fake_consistency)
    monkeypatch.setattr(_sequence_mod, '_check_safety', _fake_safety)

    sc.check()

    bands = captured['forbidden_bands']
    assert len(bands) == 1
    assert np.isclose(bands[0][0], 100.0)
    assert np.isclose(bands[0][1], 200.0)
    assert np.isclose(bands[0][2], 1.0e-3 * float(opts.gamma))
