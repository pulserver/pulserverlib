from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pytest

from pge import SequenceCollection
from pge.core._plot import plot as _plot_impl


def test_generated_sequences_wrapper_smoke(generated_seq_path):
    sc = SequenceCollection(str(generated_seq_path))
    report = sc.report()

    assert sc.num_sequences >= 1
    assert isinstance(report, list)
    assert len(report) >= 1
    assert report[0].num_blocks > 0


@pytest.mark.parametrize("threshold_percent", [80.0, [80.0, 100.0], (80.0, 100.0), []])
def test_grad_spectrum_threshold_percent_variants(
    generated_seq_path, threshold_percent
):
    sc = SequenceCollection(str(generated_seq_path))

    sc.grad_spectrum(
        sequence_idx=0,
        forbidden_bands=[(500.0, 600.0, 2000.0)],
        threshold_percent=threshold_percent,
        max_frequency=1200.0,
    )

    fig = plt.gcf()
    assert fig is not None
    # Layout can be either two-row (spectrogram + harmonics) or
    # harmonics-only when sliding windows are not meaningful.
    assert len(fig.axes) >= 3

    harmonic_axes = [
        ax
        for ax in fig.axes
        if isinstance(ax.get_title(), str) and ax.get_title().endswith('Harmonic Spectrum')
    ]
    assert len(harmonic_axes) == 3

    expected_count = (
        1 if isinstance(threshold_percent, (int, float)) else len(threshold_percent)
    )
    if expected_count > 0:
        for ax in harmonic_axes:
            labels = [ln.get_label() for ln in ax.lines]
            threshold_labels = [
                lb
                for lb in labels
                if isinstance(lb, str) and lb.startswith("threshold ")
            ]
            assert len(threshold_labels) >= expected_count

    plt.close(fig)


@pytest.mark.parametrize("threshold_percent", [80.0, [80.0, 100.0], (80.0, 100.0), []])
def test_pns_threshold_percent_variants(representative_generated_seq_path, threshold_percent):
    sc = SequenceCollection(str(representative_generated_seq_path))

    sc.pns(
        sequence_idx=0,
        stim_threshold=23.4 / 0.333,
        decay_constant_us=360.0,
        threshold_percent=threshold_percent,
    )

    fig = plt.gcf()
    assert fig is not None
    assert len(fig.axes) >= 1

    expected_count = 1 if isinstance(threshold_percent, (int, float)) else len(threshold_percent)
    labels = [ln.get_label() for ln in fig.axes[0].lines]
    threshold_labels = [lb for lb in labels if isinstance(lb, str) and lb.endswith('% threshold')]
    assert len(threshold_labels) >= expected_count

    plt.close(fig)


def test_plot_accepts_pypulseq_without_fig(simple_gre_seq):
    handle = _plot_impl(simple_gre_seq, tr_idx=0, time_unit="ms")

    assert handle is not None
    assert handle.fig is not None
    assert "gx" in handle.axes

    plt.close(handle.fig)


def test_plot_sequence_collection_without_fig(simple_gre_seq):
    sc = SequenceCollection(simple_gre_seq)
    handle = _plot_impl(sc, tr_idx=0, time_unit='ms')

    assert handle is not None
    assert handle.fig is not None
    assert 'gx' in handle.axes

    plt.close(handle.fig)


def test_plot_accepts_xml_without_fig(tmp_path):
    xml_path = tmp_path / "wf.xml"
    xml_path.write_text(
        "<root>"
        "<gx>0 10 20 0.0 1.0 0.0</gx>"
        "<gy>0 10 20 0.0 0.5 0.0</gy>"
        "<gz>0 10 20 0.0 0.2 0.0</gz>"
        "<rf>0 10 20 0.0 0.1 0.0</rf>"
        "</root>",
        encoding="utf-8",
    )

    handle = _plot_impl(str(xml_path), time_unit="us", label="xml")

    assert handle is not None
    assert handle.fig is not None
    assert "rf_mag" in handle.axes

    plt.close(handle.fig)


def test_validate_handles_unsorted_reference_times(simple_gre_seq, monkeypatch):
    import pge.core._validate as _validate_mod

    sc = SequenceCollection(simple_gre_seq)
    orig = _validate_mod._pypulseq_reference

    def _unsorted_reference(seq, sequence_idx, tr_idx, tr_info):
        ref = orig(seq, sequence_idx, tr_idx, tr_info)
        out = {}
        for key, value in ref.items():
            t, a = value
            t_arr = np.asarray(t)
            a_arr = np.asarray(a)
            if t_arr.size > 1:
                idx = np.arange(t_arr.size)[::-1]
                out[key] = (t_arr[idx], a_arr[idx])
            else:
                out[key] = (t_arr, a_arr)
        return out

    monkeypatch.setattr(_validate_mod, "_pypulseq_reference", _unsorted_reference)

    result = sc.validate(
        sequence_idx=0,
        do_plot=False,
        tr_range=(0, 1),
    )

    assert result is None


def test_validate_passes_for_supported_generated_sequences(validate_pass_seq_path, capsys):
    sc = SequenceCollection(str(validate_pass_seq_path))

    result = sc.validate(sequence_idx=0, do_plot=False, tr_range=(0, 1))

    captured = capsys.readouterr()
    assert result is None
    assert 'Validation passed.' in captured.out


def test_validate_bssfp_full_range_includes_prep_and_cooldown(capsys):
    sc = SequenceCollection(
        str(
            Path('/home/mcencini/pulserverlib/tests/utils/expected')
            / 'bssfp_2d_1sl_1avg.seq'
        )
    )

    result = sc.validate(sequence_idx=0, do_plot=False, tr_range=(0, 8))

    captured = capsys.readouterr()
    assert result is None
    assert 'Validation passed.' in captured.out


def test_validate_known_generated_failures_raise(known_validate_failure_seq_path, capsys):
    sc = SequenceCollection(str(known_validate_failure_seq_path))

    with pytest.raises(RuntimeError, match='Validation failed'):
        sc.validate(sequence_idx=0, do_plot=False, tr_range=(0, 1))

    captured = capsys.readouterr()
    assert 'Validation failed:' in captured.out


def test_validate_raises_on_failure(simple_gre_seq, monkeypatch):
    import pge.core._validate as _validate_mod

    sc = SequenceCollection(simple_gre_seq)
    orig = _validate_mod._pypulseq_reference

    def _bad_reference(seq, sequence_idx, tr_idx, tr_info):
        ref = orig(seq, sequence_idx, tr_idx, tr_info)
        out = {}
        for key, value in ref.items():
            t, a = value
            t_arr = np.asarray(t)
            a_arr = np.asarray(a)
            if t_arr.size > 0:
                out[key] = (t_arr + 500.0, a_arr)
            else:
                out[key] = (t_arr, a_arr)
        return out

    monkeypatch.setattr(_validate_mod, '_pypulseq_reference', _bad_reference)

    with pytest.raises(RuntimeError, match='Validation failed'):
        sc.validate(
            sequence_idx=0,
            do_plot=False,
            tr_range=(0, 1),
            grad_atol=1e-9,
            rf_rms_percent=1e-9,
        )
