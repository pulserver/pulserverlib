# Mechanical Resonance Safety: Structural Acoustic Analysis

## Overview

MRI gradient systems produce acoustic noise whose spectral content is
determined by the temporal pattern of gradient pulses within and across
repetition periods.  Certain frequencies excite mechanical resonances of
the scanner bore, cryostat, or gradient coil assembly, potentially causing
hardware damage or exceeding acoustic safety limits.  The structural
acoustic analysis module computes a physics-informed spectral model of the
gradient waveform, identifies candidate resonance frequencies, and checks
each candidate against user-defined forbidden frequency bands.

The analysis is performed **per gradient axis independently** (Gx, Gy, Gz)
and **per subsequence** in a sequence collection.  Forbidden bands from
all axes are collected into a single union set: a violation on *any* axis
flags the candidate.

## Analytical spectral model

The spectrum is built from the *structural timing* of gradient events
within one canonical repetition period, without relying on a windowed FFT.
The pipeline has four stages.

### Stage 1 — event extraction

For each axis, every block in the canonical TR window is inspected.
If a block contains a gradient event, an *occurrence* is recorded:

- **def_id** — gradient definition index.
- **start_time** — cumulative time from the start of the TR (µs).
- **amplitude** — worst-case positional amplitude
  $|\text{gte.amplitude}| \times |\text{gdef.max\_amplitude}|$ (Hz/m).

### Stage 2 — waveform response model

Each gradient definition is classified into one of three response models.

**Trapezoid gradients** are represented as a 4-vertex piecewise-linear
(PWL) waveform:

$$(0, 0) \;\to\; (\tau_r, 1) \;\to\; (\tau_r + \tau_f, 1) \;\to\; (\tau_r + \tau_f + \tau_d, 0)$$

The normalised Fourier transform is computed analytically from the
vertices.  Because the waveform is continuous and starts and ends at
zero, the response decays as $\sim 1/f^2$.

**Arbitrary waveforms with many samples** (≥ 10) are tested for internal
sub-periodicity via normalised autocorrelation.  If a repeating sub-period
is detected (autocorrelation > 0.5, minimum 3 repetitions), the waveform
is *decomposed*: each sub-period becomes a separate event whose PWL model
is the sub-period shape, and consecutive sub-periods emit separate events
at spaced start times.  If no sub-period is detected, the full waveform is
approximated by a uniformly-sampled PWL with up to 16 vertices.

**Arbitrary waveforms with few samples** (< 10) are treated as PWL
directly from their sample values.

### Stage 3 — repetition tagging (multi-average)

For sequences with non-degenerate prep/cooldown blocks (i.e. where the
prep or cooldown varies between passes), the canonical TR is the full
*pass* (prep + imaging + cooldown).  When `num_averages > 1`, the imaging
section repeats `num_averages` times within each pass.  Rather than
explicitly expanding all copies, each imaging event is tagged with:

- **num_reps** = `num_averages`
- **rep_period_us** = total imaging section duration

Cooldown events are shifted in time to their position in the expanded
pass.  Prep events are left at their original position with
`num_reps = 1`.

This lets the Dirichlet kernel (see below) handle repetition analytically
in $O(1)$ per event per frequency, instead of enumerating all copies.

### Stage 4 — spectral evaluation

At a given frequency $f$, the complex contribution of a single event $k$
is:

$$a_k(f) \;=\; A_k \; W_k(f) \; e^{-j\,2\pi f\, t_k}$$

where $A_k$ is the signed amplitude (Hz/m), $W_k(f)$ the normalised PWL
response, and $t_k$ the start time within the TR.

When an event has `num_reps` $= N > 1$, its contribution is multiplied by
the Dirichlet kernel:

$$D_N(f, T) \;=\; \frac{\sin(N\,\pi\, f\, T)}{\sin(\pi\, f\, T)} \;\cdot\; e^{-j\,(N-1)\,\pi\, f\, T}$$

where $T$ is `rep_period_us` in seconds.  This preserves complex phase and
produces sharp peaks at multiples of $1/T$.

The **coherent magnitude** on one axis is the magnitude of the complex sum
over all events:

$$M_\text{coh}(f) \;=\; \Bigl|\,\sum_k a_k(f)\,D_{N_k}(f,T_k)\,\Bigr|$$

The **incoherent magnitude** is the root-sum-of-powers:

$$M_\text{inc}(f) \;=\; \sqrt{\sum_k N_k \,|a_k(f)|^2}$$

The coherence ratio $\rho = M_\text{coh} / M_\text{inc}$ measures how
strongly events interfere constructively at $f$.

## Evaluation frequency grid

The coherent and incoherent magnitudes are evaluated on a composite
frequency grid:

1. **TR harmonics**: $f_m = m \cdot f_1$ for
   $m = 1, \ldots, \lfloor f_\text{max}/f_1 \rfloor$, where
   $f_1 = 1/T_\text{TR}$ is the fundamental of the canonical TR.
2. **FFT-promoted peaks**: any local maximum of the dense FFT RSS spectrum
   that exceeds 15% of the global FFT peak and falls more than 25% of the
   harmonic spacing away from the nearest TR harmonic is added as an extra
   evaluation frequency.  This catches spectral features that land between
   adjacent TR harmonics (e.g. bipolar-readout fundamentals in bSSFP).

## Candidate selection

After computing magnitudes at every evaluation frequency, candidates are
selected through a two-stage filter.

### Analytical power gate

The cross-axis RSS² of the coherent magnitudes,

$$P(f) \;=\; \sum_\text{ax} M_\text{coh,ax}^2(f)$$

is computed.  Any frequency where $P(f) < 0.05 \times \max_f P(f)$ is
rejected, unless it is present in the **FFT bypass set** (see below).

### FFT bypass set

The analytical PWL response $W_k(f)$ can underestimate power at high-order
harmonics where the true waveform differs from the PWL approximation.  To
compensate, the dense FFT is scanned for local maxima whose RSS exceeds
15% of the global FFT peak.  Each such FFT peak is snapped to its nearest
TR harmonic, and that harmonic is marked to bypass the analytical power
gate.  FFT-bypassed TR harmonics also use a relaxed tier check (amplitude
floor only, no coherence ratio requirement), since the FFT already confirms
a real spectral peak.

### Per-axis tier checks

Frequencies that survive the power gate are evaluated per axis through
three tiers:

| Tier | Condition |
|------|-----------|
| **Tier 1** — multi-event coherence | ≥ 2 events on the axis, coherence ratio $\rho > 2$, and $M_\text{coh} > 1$ Hz/m |
| **Tier 2** — single-event amplitude | Exactly 1 event and $M_\text{coh} > 10^6$ Hz/m |
| **Tier 3** — FFT-promoted peaks | Between TR harmonics; $M_\text{coh} > 1$ Hz/m |

A frequency becomes a **candidate** if any axis qualifies through any
tier.  The candidate list is shared across axes.

## Forbidden-band check and effective gradient amplitude

For each candidate frequency $f_c$, the effective gradient amplitude
$G_\text{eff}$ is computed **per axis**:

$$G_{\text{eff},\text{ax}}(f_c) \;=\; \frac{\displaystyle\sum_k |A_k|\;|a_k(f_c)|\;|D_{N_k}(f_c,T_k)|}{\displaystyle\sum_k |a_k(f_c)|\;|D_{N_k}(f_c,T_k)|}$$

This is a spectral-weight–averaged amplitude: events that contribute more
to the spectral peak at $f_c$ dominate the average.  The result is in the
same units as the time-domain gradient amplitude (Hz/m).

Each candidate is then checked against the union set of forbidden bands:

$$f_\text{lo} \;\le\; f_c \;\le\; f_\text{hi} \quad\text{and}\quad G_{\text{eff},\text{ax}}(f_c) \;>\; A_\text{limit}$$

A **violation** is flagged if *any single axis* exceeds the band limit.
This per-axis check ensures that a resonance driven by one gradient channel
is caught even if the other channels are quiet at that frequency.

## Per-subsequence processing

A sequence collection may contain multiple subsequences (e.g. a localiser
followed by a scan).  The analysis iterates over each subsequence
independently:

1. Identify the canonical TR window (or pass window for non-degenerate
   sequences) and enumerate unique TR variants by shot ID.
2. For each unique variant, extract the gradient waveforms, run the
   structural analysis, and check candidates against the forbidden bands.
3. In the safety-check path, the first violation triggers an immediate
   failure with a diagnostic message identifying the subsequence, TR
   variant, axis, frequency, amplitude, and band that was exceeded.
