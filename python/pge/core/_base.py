"""Abstract base class for pulserver sequence plugins."""

from __future__ import annotations

from abc import ABC, abstractmethod

import pypulseq as pp

from ._params import Protocol


class PulseqSequence(ABC):
    """Abstract base class for pulserver sequence plugin implementations.

    Subclass this abstract base class to create a custom sequence plugin
    that integrates with the pulserver framework. Implement all three
    abstract methods to satisfy the contract. The plugin discovery bridge
    automatically locates subclasses via Python's ``inspect`` module.

    Notes
    -----
    Subclasses should be placed in modules that the plugin loader can
    discover (typically in a plugins directory or registered via entry points).
    The bridge will instantiate the subclass and call its methods in order:
    1. ``get_default_protocol()`` (once, at load time)
    2. ``validate_protocol()`` (repeatedly, on parameter changes)
    3. ``make_sequence()`` (on sequence generation request)

    Examples
    --------
    >>> from python.pge.core import PulseqSequence, Protocol
    >>> import pypulseq as pp
    >>>
    >>> class MyGRESequence(PulseqSequence):
    ...     def get_default_protocol(self, opts: pp.Opts) -> Protocol:
    ...         # Return default protocol dict
    ...         from python.pge.core._params import FloatParam, IntParam
    ...         return {
    ...             'TE': FloatParam(value=10.0, min=5.0, max=50.0, incr=1.0, unit='ms'),
    ...             'TR': FloatParam(value=500.0, min=100.0, max=5000.0, incr=10.0, unit='ms'),
    ...         }
    ...
    ...     def validate_protocol(self, opts: pp.Opts, protocol: Protocol) -> dict:
    ...         # Validate and return extended info
    ...         te = protocol['TE'].value
    ...         tr = protocol['TR'].value
    ...         valid = te < tr and te > 0
    ...         duration_s = 1.0  # placeholder
    ...         return {'valid': valid, 'duration': duration_s, 'info': 'GRE sequence'}
    ...
    ...     def make_sequence(self, opts: pp.Opts, protocol: Protocol, output_path: str) -> None:
    ...         # Build sequence and write to output_path
    ...         seq = pp.Sequence(system=pp.Opts())
    ...         seq.write(output_path)
    """

    @abstractmethod
    def get_default_protocol(self, opts: pp.Opts) -> Protocol:
        """Return the default protocol for this sequence.

        Called once when the plugin is loaded. Should return a fully
        instantiated protocol dict with sensible defaults.

        Parameters
        ----------
        opts : pp.Opts
            System hardware specification (gradients, RF, timing constraints).

        Returns
        -------
        Protocol
            Dictionary mapping UIParam (or str) keys to ProtocolValue dataclasses
            (FloatParam, IntParam, BoolParam, StringListParam, Description).
            Should include all parameters that the sequence supports.

        Notes
        -----
        The returned protocol is used as the baseline for all subsequent
        validation and sequence generation calls. Must be valid (i.e.,
        would pass validation against the provided opts).

        Examples
        --------
        >>> from python.pge.core._params import FloatParam, UIParam
        >>> def get_default_protocol(self, opts):
        ...     return {
        ...         UIParam.TE: FloatParam(10.0, 5.0, 50.0, 1.0, 'ms'),
        ...         UIParam.TR: FloatParam(500.0, 100.0, 5000.0, 10.0, 'ms'),
        ...         UIParam.FLIP_ANGLE: FloatParam(90.0, 1.0, 180.0, 1.0, 'deg'),
        ...     }
        """
        ...

    @abstractmethod
    def validate_protocol(self, opts: pp.Opts, protocol: Protocol) -> dict:
        """Validate *protocol* against hardware *opts*.

        Called on every parameter change — must be fast with no file I/O
        or external communication. Typically checks parameter bounds,
        interdependencies, and hardware constraints.

        Parameters
        ----------
        opts : pp.Opts
            System hardware specification (gradients, RF, timing constraints).
        protocol : Protocol
            Current parameter dictionary (may include out-of-bounds values).

        Returns
        -------
        dict
            Dictionary with keys:

            - ``'valid'`` (bool): Whether the protocol is valid for the
              given hardware and parameter values.
            - ``'duration'`` (float | None): Estimated scan time in seconds
              (for valid protocols; ``None`` if invalid).
            - ``'info'`` (str | None): Human-readable message explaining
              the validation result (e.g., error reason or duration estimate).

        Notes
        -----
        This method is called frequently during interactive parameter adjustment.
        Performance is critical — avoid blocking operations, file I/O, or
        external calls. Keep validation logic deterministic and fast.

        Examples
        --------
        >>> def validate_protocol(self, opts, protocol):
        ...     te = protocol['TE'].value
        ...     tr = protocol['TR'].value
        ...     if te >= tr:
        ...         return {'valid': False, 'duration': None, 'info': 'TE >= TR'}
        ...     if tr > opts.rf_dead_time:
        ...         duration = 1.0  # placeholder calculation
        ...         return {'valid': True, 'duration': duration, 'info': f'Scan time: {duration}s'}
        ...     return {'valid': False, 'duration': None, 'info': 'TR too short'}
        """
        ...

    @abstractmethod
    def make_sequence(
        self, opts: pp.Opts, protocol: Protocol, output_path: str
    ) -> None:
        """Build the full sequence and write the ``.seq`` file to *output_path*.

        Called when the user requests sequence generation. Expects *protocol*
        to be valid (i.e., a prior call to ``validate_protocol()`` returned
        ``'valid': True``). Responsibility of the bridge to ensure this precondition.

        Parameters
        ----------
        opts : pp.Opts
            System hardware specification used for sequence construction.
        protocol : Protocol
            Validated parameter dictionary (guaranteed to pass validation).
        output_path : str
            File system path where the ``.seq`` file should be written.
            Path is guaranteed to be writable.

        Raises
        ------
        RuntimeError
            If sequence construction fails (e.g., inconsistent parameters,
            unsupported hardware constraints). The error message should be
            human-readable.

        Notes
        -----
        This method may take several seconds (sequence generation can be
        computationally expensive). It runs asynchronously in the bridge
        and the user is notified upon completion or failure.

        File I/O is allowed (and expected) in this method. Implementations
        may write temporary files or perform other setup as needed.

        Examples
        --------
        >>> def make_sequence(self, opts, protocol, output_path):
        ...     seq = pp.Sequence(system=opts)
        ...     te_ms = protocol['TE'].value
        ...     tr_ms = protocol['TR'].value
        ...     # Build sequence using te_ms and tr_ms...
        ...     seq.write(output_path)
        """
        ...
