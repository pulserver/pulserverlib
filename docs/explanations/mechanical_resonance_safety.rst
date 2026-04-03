Mechanical Resonance Safety: Structural Acoustic Analysis
==========================================================

Overview
--------

MRI gradient systems produce acoustic noise whose spectral content is
determined by the temporal pattern of gradient pulses within and across
repetition periods.  Certain frequencies excite mechanical resonances of
the scanner bore, cryostat, or gradient coil assembly, potentially causing
hardware damage or exceeding acoustic safety limits.  The structural
acoustic analysis module computes a physics-informed spectral model of the
gradient waveform, identifies candidate resonance frequencies, and
checks each candidate against user-defined forbidden frequency bands.

Problem statement
~~~~~~~~~~~~~~~~~

Given a pulse sequence whose canonical TR is repeated :math:`K` times, we
need to:

1. Predict the acoustic spectrum without relying on a measured or
   windowed FFT — i.e.\ analytically, from the structural timing of the
   gradient definitions.
2. Identify the discrete set of candidate resonance frequencies where the
   spectrum concentrates its energy.
3. For each candidate, determine the worst-case gradient amplitude among
   the gradient definitions that contribute to that peak.
4. Flag candidates that fall inside forbidden bands and exceed the
   amplitude limit of the band.


Analytical spectral model
-------------------------

Definitions and notation
~~~~~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :widths: 20 80

   * - :math:`K`
     - Number of TR repetitions (``num_instances``).
   * - :math:`T_\mathrm{TR}`
     - Repetition time (seconds).
   * - :math:`d`
     - Index over gradient definitions on a single axis.
   * - :math:`A_d`
     - Worst-case time-domain amplitude of definition :math:`d` (Hz/m).
   * - :math:`H_d(f)`
     - Intra-pulse waveform response of definition :math:`d`.
   * - :math:`R_d`
     - Number of spacing runs for definition :math:`d` within the canonical TR.
   * - :math:`T_r`, :math:`N_r`, :math:`t_r`
     - Spacing, occurrence count, and start time of run :math:`r`.
   * - :math:`D_N(x)`
     - Dirichlet kernel :math:`\sin(N\pi x)/\sin(\pi x)`.


Per-definition spectrum
~~~~~~~~~~~~~~~~~~~~~~~

Each gradient definition :math:`d` produces a complex spectral
contribution:

.. math::

   S_d(f) \;=\; A_d \, H_d(f) \;
              \sum_{r=0}^{R_d - 1}
                 e^{-j\,2\pi f \, t_r}\;
                 D_{N_r}\!\bigl(f \, T_r\bigr)

The Dirichlet kernel is evaluated with its complex phase:

.. math::

   e^{-j\,2\pi f \, t_r}\, D_{N_r}(f\,T_r)
   \;=\;
   \frac{\sin(N_r \,\pi\, f\, T_r)}{\sin(\pi\, f\, T_r)}
   \;\cdot\;
   e^{-j\,2\pi f \, t_r}

Each run within a definition contributes a frequency comb whose peak
positions are at multiples of :math:`1/T_r` and whose lobe width narrows
as :math:`N_r` increases.  The complex exponential encodes the absolute
start time :math:`t_r` of the run within the TR, so that interference
between runs from the same or different definitions is correctly modelled.

When a definition appears only once in the TR (:math:`R_d = 0` runs), the
inner sum is unity and the spectrum reduces to :math:`A_d \, H_d(f)`.


Per-axis spectrum
~~~~~~~~~~~~~~~~~

On each gradient axis, the contributions from all definitions are summed
coherently (complex sum) before taking the magnitude:

.. math::

   \tilde{S}_\mathrm{axis}(f)
   \;=\;
   \sum_{d} S_d(f)

The outer TR Dirichlet kernel accounts for the :math:`K`-fold repetition
of the canonical TR:

.. math::

   S_\mathrm{axis}(f)
   \;=\;
   \bigl|\,
      D_K\!\bigl(f\,T_\mathrm{TR}\bigr)
      \;\cdot\;
      \tilde{S}_\mathrm{axis}(f)
   \,\bigr|

where the product is a complex multiplication and the final magnitude is
taken.  This concentrates the intra-TR spectral content at harmonics of
:math:`1/T_\mathrm{TR}`, with lobe width inversely proportional to
:math:`K`.


Cross-axis combined spectrum
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The three axes are combined by taking the element-wise maximum:

.. math::

   S_\mathrm{combined}(f)
   \;=\;
   \max\bigl(
      S_\mathrm{gx}(f),\;
      S_\mathrm{gy}(f),\;
      S_\mathrm{gz}(f)
   \bigr)

This is conservative: a forbidden-band violation on any axis is detected.


Waveform response :math:`H_d(f)`
---------------------------------

The intra-pulse waveform response :math:`H_d(f)` captures the spectral
content of a single gradient pulse.  It is normalised so that
:math:`H_d(0) = 1` and characterises how rapidly acoustic energy decays
at higher frequencies.  Four cases are handled:

Case 1: Trapezoid gradients — piecewise-linear analytical FT
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A trapezoid gradient with rise time :math:`\tau_r`, flat time
:math:`\tau_f`, and fall time :math:`\tau_d` is a piecewise-linear
waveform with four vertices:

.. math::

   (0,\,0) \;\to\; (\tau_r,\,1) \;\to\; (\tau_r + \tau_f,\,1)
   \;\to\; (\tau_r + \tau_f + \tau_d,\,0)

The Fourier transform of any piecewise-linear waveform
:math:`g(t)` defined by :math:`n` vertices :math:`\{(t_k, v_k)\}` is:

.. math::

   G(f) = \sum_{k=0}^{n-2}
            e^{-j\omega t_k}
            \bigl[a_k \, I_0(\omega, T_k) + b_k \, I_1(\omega, T_k)\bigr]

where :math:`\omega = 2\pi f`, :math:`T_k = t_{k+1} - t_k`,
:math:`a_k = v_k`, :math:`b_k = (v_{k+1} - v_k)/T_k`, and:

.. math::

   I_0(\omega, T) &= \frac{1 - e^{-j\omega T}}{j\omega} \\[4pt]
   I_1(\omega, T) &= \frac{I_0(\omega, T) - T\,e^{-j\omega T}}{j\omega}

The normalised response is:

.. math::

   H_d(f) = \frac{|G(f)|}{|G(0)|}

Because :math:`g(t)` is continuous and starts and ends at zero, the
leading-order behaviour at high frequency is:

.. math::

   H_d(f) \;\sim\; \frac{1}{(2\pi f)^2}
   \qquad (f \to \infty)

This :math:`1/f^2` decay correctly attenuates high-order harmonics and
avoids the generation of spurious high-frequency candidates.

Case 2: Arbitrary waveforms with many samples and uniform time
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For arbitrary (free-form) gradient waveforms with at least 10 samples on
a uniform time raster, the spectral content is determined via FFT-based
peak detection:

1. Decompress the shape samples from the sequence file.
2. Remove the DC component: :math:`g_s \leftarrow g_s - \bar{g}`.
3. Apply a Hann taper:
   :math:`w_s = \tfrac{1}{2}\bigl(1 - \cos\!\bigl(\tfrac{2\pi(s+1)}{N_\mathrm{wf}}\bigr)\bigr)`.
4. Zero-pad to the next power of two and compute the real FFT.
5. Run the resonance peak detector on the magnitude spectrum.

If a peak is detected, the waveform is classified as *resonant* at
:math:`f_\mathrm{peak}` (the frequency of the absolute spectral maximum)
with effective cycle count:

.. math::

   N_\mathrm{eff} = \mathrm{round}(T_\mathrm{wf} \times f_\mathrm{peak})

The response is a normalised Dirichlet kernel:

.. math::

   H_d(f) = \frac{1}{N_\mathrm{eff}}
   \;\Bigl|\, D_{N_\mathrm{eff}}\!\Bigl(\frac{f}{f_\mathrm{peak}}\Bigr) \,\Bigr|

If no peak is detected, the waveform is classified as broadband:
:math:`H_d(f) = 1`.

Case 3: Arbitrary waveforms with non-uniform time shape
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Some Pulseq arbitrary gradients carry an explicit non-uniform time
array (``time_shape_id > 0``).  In this case the waveform samples are
not on a regular raster and cannot be directly FFT'd.

The procedure is:

1. Decompress both the amplitude shape (scale = 1.0) and the time shape
   (scale = ``grad_raster_us``) to obtain the non-uniform
   :math:`(t_i, g_i)` pairs.
2. Linearly interpolate onto a uniform grid at the gradient raster
   interval.
3. Apply the same FFT-based resonance detection as Case 2.

Case 4: Arbitrary waveforms with few samples (< 10)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Arbitrary waveforms with fewer than 10 samples typically represent
simple piecewise-linear shapes (ramps, flat-top segments, spoiler
gradients) stored in the Pulseq arbitrary format rather than as a
trapezoid.  Using FFT on so few points would be numerically unstable.

Instead, these are treated identically to trapezoids: the sample
values and their time coordinates (from the time shape or uniform raster)
define a piecewise-linear waveform whose Fourier transform is computed
analytically via the same :math:`I_0 / I_1` formulation.  This produces
the correct :math:`\sim 1/f^2` decay.


Hierarchical spacing extraction
---------------------------------

The timing structure of the canonical TR is extracted in four stages.

Stage 1: occurrence extraction
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For each gradient axis, every block in the canonical TR window is
inspected.  If a block contains a gradient event on the axis, an
*occurrence* is recorded with:

- **def_id**: gradient definition index
- **start_time**: cumulative time from the start of the TR (µs)
- **amplitude**: :math:`|\text{gte.amplitude}| \times |\text{gdef.max\_amplitude}|`
  (Hz/m)

The result is an ordered list of all gradient pulse firings on one axis
within one TR.

Stage 2: deduplication
~~~~~~~~~~~~~~~~~~~~~~

Occurrences are grouped by ``def_id``.  For each unique gradient
definition, the list of occurrence times is extracted and passed to the
spacing clusterer.

Stage 3: hierarchical spacing clustering
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Given :math:`M` occurrence times :math:`\{t_0, t_1, \ldots, t_{M-1}\}`
for one gradient definition, the consecutive spacings are:

.. math::

   \Delta_i = t_{i+1} - t_i, \qquad i = 0, \ldots, M-2

These spacings are clustered into *runs* of approximately equal duration.

**Level-1 clustering (consecutive runs):**

1. Initialise a run with spacing :math:`\Delta_0`.
2. For each subsequent :math:`\Delta_i`, compute a tolerance:
   :math:`\epsilon = \max(\epsilon_\mathrm{abs},\; \epsilon_\mathrm{rel} \times \bar{\Delta}_\mathrm{run})`.
3. If :math:`|\Delta_i - \bar{\Delta}_\mathrm{run}| \le \epsilon`:
   extend the current run and update the running average.
4. Otherwise: close the current run and start a new one.

Each run stores its average spacing :math:`T_r`, its occurrence count
:math:`N_r`, and the time of its first member :math:`t_r`.

**Level-2 clustering (hierarchical):**

If the level-1 step produces :math:`\ge 2` runs, their start times are
treated as a new set of occurrence times and the clustering is applied
recursively.  This captures timing structure at multiple scales.

**Example: multi-slice multi-echo EPI.**
A readout gradient in such a sequence fires at times governed by (at
least) three scales:

- **ESP** (echo spacing, ~0.5 ms): consecutive echoes within one echo
  train.
- **Δ TE** or **T_slice** (~15 ms): gap between echo trains of
  successive slices.
- **Longer pauses**: gaps between the last echo of one slice group and
  the first echo of the next.

Level-1 clustering captures the ESP runs.  Level-2 clustering, applied to
the run start times, captures the slice-level periodicity.  Each scale
becomes a separate run entry in the output, so that the Dirichlet kernel
sum in :math:`S_d(f)` naturally produces interference patterns at all
relevant harmonic families.

Stage 4: arbitrary waveform analysis
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

After clustering, each contribution from an arbitrary-waveform gradient
definition undergoes the FFT-based resonance detection described above to
determine :math:`f_\mathrm{peak}` and :math:`N_\mathrm{eff}`.


Candidate detection
-------------------

Dense analytical evaluation
~~~~~~~~~~~~~~~~~~~~~~~~~~~

The three per-axis spectra :math:`S_\mathrm{gx}(f)`,
:math:`S_\mathrm{gy}(f)`, :math:`S_\mathrm{gz}(f)` and the combined
spectrum :math:`S_\mathrm{combined}(f)` are evaluated on every bin of the
frequency grid used by the sliding-window FFT analysis:

.. math::

   f_i = f_\mathrm{min} + i \, \Delta f, \qquad i = 0, \ldots, N_\mathrm{bins} - 1

This reuses the same grid populated by the FFT-based spectrogram, so
that the analytical and FFT spectra share a common frequency axis for
overlay plotting.

Prominence-based peak detection
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Candidates are the local maxima of :math:`S_\mathrm{combined}(f)` that
satisfy:

1. **Amplitude threshold:** :math:`S_\mathrm{combined}(f_i) > 0.01 \times S_\mathrm{max}`.
2. **Prominence threshold:** the peak prominence exceeds
   :math:`0.20 \times S_\mathrm{max}`.

The *prominence* of a local maximum is the height it rises above the
deeper of its two flanking valleys, where a valley extends until a higher
peak is encountered.  This metric naturally distinguishes main Dirichlet
lobes (high prominence, spanning a full null-to-null interval) from
inter-lobe sidelobes (low prominence, sitting between closely-spaced
nulls).

For a simple GRE with :math:`K = 8` repetitions at
:math:`T_\mathrm{TR} = 6.7` ms, the main peaks at harmonics of
:math:`1/T_\mathrm{TR} \approx 149` Hz each have prominence close to 100%
of the global maximum, while the six sidelobes between each pair of main
peaks have prominence of 10–15%.  The 20% threshold cleanly separates
these two populations.


Per-peak gradient amplitude attribution
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For each candidate frequency :math:`f_c`, the worst-case gradient amplitude
is determined among the gradient definitions that actually contribute to
the peak:

1. On each axis, evaluate the complex per-definition spectra:

   .. math::

      S_d(f_c) = A_d \, H_d(f_c) \, \sum_r e^{-j2\pi f_c t_r}\, D_{N_r}(f_c T_r)

2. Compute the axis total:

   .. math::

      S_\mathrm{axis}^\mathrm{inner}(f_c)
      = \Bigl|\sum_d S_d(f_c)\Bigr|

3. A definition :math:`d` is considered a *contributor* if:

   .. math::

      |S_d(f_c)| > 0.01 \times S_\mathrm{axis}^\mathrm{inner}(f_c)

4. The attributed gradient amplitude for the candidate is:

   .. math::

      A_\mathrm{max}(f_c) = \max_{d \in \text{contributors}} A_d

Note that the 1% threshold on the inner (pre-TR-Dirichlet) sum means
the outer :math:`D_K` kernel cancels in the ratio, so the contributor
set does not depend on :math:`K`.


Forbidden-band violation check
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Each candidate :math:`(f_c, A_\mathrm{max})` is checked against the list
of forbidden bands :math:`\{[f_\mathrm{lo},\, f_\mathrm{hi},\,
A_\mathrm{limit}]\}`.  A violation is flagged when:

.. math::

   f_\mathrm{lo} \le f_c \le f_\mathrm{hi}
   \quad\text{and}\quad
   A_\mathrm{max}(f_c) > A_\mathrm{limit}


Output
------

The analysis produces the following outputs:

.. list-table::
   :widths: 30 70

   * - ``analytical_gx/gy/gz``
     - Dense analytical spectrum arrays :math:`S_\mathrm{axis}(f_i)`,
       one value per frequency bin, for plotting as an overlay on the
       FFT-based spectrum.
   * - ``candidate_freqs``
     - Frequencies of the detected prominent peaks in
       :math:`S_\mathrm{combined}(f)`.
   * - ``candidate_amps_gx/gy/gz``
     - Per-axis analytical amplitudes at each candidate frequency.
   * - ``candidate_grad_amps``
     - Worst-case time-domain gradient amplitude (Hz/m) among
       contributing definitions, per candidate.
   * - ``candidate_violations``
     - Binary flag per candidate: 1 if the candidate falls inside a
       forbidden band and its attributed amplitude exceeds the band
       limit.


Discussion
----------

Harmonics and spectral decay
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For a sequence with :math:`K` identical TR repetitions, the outer
Dirichlet kernel :math:`D_K(f \, T_\mathrm{TR})` produces peaks at every
integer multiple of :math:`1/T_\mathrm{TR}`.  At exact multiples,
:math:`|D_K| = K` regardless of harmonic order.  The *envelope* across
harmonics is then determined by the intra-TR spectral content, which is
the product :math:`H_d(f) \cdot |\sum_r D_{N_r}(f\,T_r)|`.

For trapezoid and piecewise-linear waveforms, :math:`H_d(f) \sim 1/f^2`
(since the waveform is continuous with piecewise-constant derivative),
so the harmonic amplitudes decay quadratically.  A 10× higher harmonic
thus has :math:`\sim 1\%` of the fundamental's amplitude.  This is
important in practice: a GRE with :math:`T_\mathrm{TR} = 7` ms
produces ~20 harmonics below 3 kHz, but only the lowest ~7 survive the
envelope decay and the prominence filter.

For oscillating arbitrary waveforms (e.g.\ spiral readouts), the
Dirichlet-kernel model for :math:`H_d(f)` concentrates energy near
multiples of the oscillation frequency :math:`f_\mathrm{peak}`, which
may amplify certain harmonics and suppress others.

Coherent summation and interference
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The use of complex (coherent) summation across gradient definitions and
timing runs means that constructive and destructive interference is
faithfully modelled.  Two readout lobes separated by half a period will
partially cancel at the fundamental; two lobes separated by exactly one
period will reinforce.  This is essential for sequences with multiple
echo trains per TR, where the relative timing of gradient events
determines which harmonics are amplified and which are suppressed.
