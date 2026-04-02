"""Typed protocol keys and values for pulserver sequence plugins."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from enum import StrEnum


class Validate(StrEnum):
    """Per-parameter validation strategy for handling out-of-bounds values.

    Determines how the system responds when a user sets a parameter value
    outside its specified ``[min, max]`` range. Maps to vendor UI framework
    validation concepts (nimpulseqgui PropertyValidate in GE).

    Members
    -------
    SEARCH : 'search'
        Binary-search for the nearest valid value inside [min, max].
        Useful for parameters where any valid value is acceptable.
    CLIP : 'clip'
        Clamp (clip) the value to the range [min, max].
        Useful for hard physical limits (e.g., max gradient).
    NONE : 'none'
        No automatic validation; accept any value (custom validation elsewhere).
        Useful when complex interdependencies override simple bounds.
    """

    SEARCH = 'search'  # binary-search for valid min/max
    CLIP = 'clip'  # clamp to [min, max]
    NONE = 'none'  # no auto-validation


class UIParam(StrEnum):
    """Standard MR protocol parameter keys for UI frameworks.

    Defines a canonical set of parameter names for MR pulse sequences,
    recognized by vendor UI bridges (e.g., GE nimpulseqgui). StrEnum members
    are strings themselves (zero conversion cost at language boundaries).
    Plugin authors may extend these with subclassing or raw strings for
    custom, vendor-specific parameters.

    Timing Parameters
    -----------------
    TE : 'TE'
        Echo time (ms).
    TR : 'TR'
        Repetition time (ms).
    TI : 'TI'
        Inversion time (ms, for IR pulses).

    Spatial Parameters
    ------------------
    FOV : 'FOV'
        Field of view (mm, typically).
    SLICE_THICKNESS : 'SliceThickness'
        Slice thickness (mm).
    NSLICES : 'NSlices'
        Number of slices.
    MATRIX : 'Matrix'
        Acquisition matrix (e.g., 256 or 256x256).
    NECHOES : 'NEchoes'
        Number of echoes per TR.

    Contrast Parameters
    -------------------
    FLIP_ANGLE : 'FlipAngle'
        RF excitation flip angle (degrees).
    BANDWIDTH : 'Bandwidth'
        Receiver bandwidth (Hz or Hz/pixel).

    Sequence Feature Flags
    ----------------------
    FAT_SAT : 'FatSat'
        Fat saturation enabled.
    SPOILER : 'Spoiler'
        Spoiler gradient enabled.
    RF_SPOILING : 'RFSpoiling'
        RF spoiling enabled.

    Scan Timing
    -----------
    TA : 'TA'
        Total acquisition time (seconds).
    """

    # Timing
    TE = 'TE'
    TR = 'TR'
    TI = 'TI'
    # Spatial
    FOV = 'FOV'
    SLICE_THICKNESS = 'SliceThickness'
    NSLICES = 'NSlices'
    MATRIX = 'Matrix'
    NECHOES = 'NEchoes'
    # Contrast
    FLIP_ANGLE = 'FlipAngle'
    BANDWIDTH = 'Bandwidth'
    # Flags
    FAT_SAT = 'FatSat'
    SPOILER = 'Spoiler'
    RF_SPOILING = 'RFSpoiling'
    # Description row
    TA = 'TA'

    @staticmethod
    def user(n: int) -> str:
        """GE user CV slot: UIParam.user(0) -> 'User0'."""
        return f'User{n}'


# ---------------------------------------------------------------------------
# Protocol value dataclasses
# ---------------------------------------------------------------------------


@dataclass
class FloatParam:
    """Floating-point protocol parameter with bounds and validation step.

    Represents a numeric slider parameter that can be adjusted between
    ``min`` and ``max`` with specified increment ``incr``. Useful for
    timing parameters (TE, TR), spatial parameters (FOV), and other
    floating-point MR protocol settings.

    Attributes
    ----------
    value : float
        Current parameter value.
    min : float
        Minimum allowed value.
    max : float
        Maximum allowed value.
    incr : float
        Step size for slider increments or binary searches.
    unit : str, optional
        Unit string for display (e.g., ``'ms'``, ``'mm'``, ``'Hz'``).
        Empty string if unitless (default).
    validate : Validate, default Validate.SEARCH
        Validation strategy when user sets value outside bounds:
        - SEARCH: binary-search for valid min/max
        - CLIP: clamp to range [min, max]
        - NONE: accept any value
    type : str, default 'float'
        Serialization type tag; reserved for protocol encoding.
    """

    value: float
    min: float
    max: float
    incr: float
    unit: str = ''
    validate: Validate = Validate.SEARCH
    type: str = 'float'


@dataclass
class IntParam:
    """Integer protocol parameter with bounds and validation step.

    Represents an integer-only parameter (e.g., number of slices, matrix size,
    echo count) that can be changed via spinner or slider with specified
    increment ``incr``. Validation strategies ensure the value remains
    within bounds.

    Attributes
    ----------
    value : int
        Current parameter value (integer).
    min : int
        Minimum allowed value.
    max : int
        Maximum allowed value.
    incr : int
        Step size for spinner/slider increments or binary searches.
    unit : str, optional
        Unit string for display (e.g., ``'slices'``, ``'lines'``).
        Empty string if unitless (default).
    validate : Validate, default Validate.SEARCH
        Validation strategy when user sets value outside bounds:
        - SEARCH: binary-search for valid min/max
        - CLIP: clamp to range [min, max]
        - NONE: accept any value
    type : str, default 'int'
        Serialization type tag; reserved for protocol encoding.
    """

    value: int
    min: int
    max: int
    incr: int
    unit: str = ''
    validate: Validate = Validate.SEARCH
    type: str = 'int'


@dataclass
class BoolParam:
    """Boolean toggle (on/off) protocol parameter.

    Represents a simple on/off switch parameter (e.g., FatSat, RF Spoiling,
    Spoiler Enable) with no bounds or validation complexity.

    Attributes
    ----------
    value : bool
        Current toggle state (``True`` = enabled, ``False`` = disabled).
    type : str, default 'bool'
        Serialization type tag; reserved for protocol encoding.
    """

    value: bool
    type: str = 'bool'


@dataclass
class StringListParam:
    """Dropdown (enumerated choice) protocol parameter.

    Represents a selector with a fixed list of string options, where
    only one option is active at a time (identified by ``index``).
    Useful for mode selectors (Pulse Type, Excitation Mode, etc.).

    Attributes
    ----------
    options : list[str]
        Available option strings (e.g., ``['SMS', 'Multi-Band', 'Single-Slice']``).
    index : int
        Current selection index (0-based) into the ``options`` list.
        Must be 0 ≤ ``index`` < ``len(options)``.
    type : str, default 'stringlist'
        Serialization type tag; reserved for protocol encoding.
    """

    options: list[str]
    index: int
    type: str = 'stringlist'


@dataclass
class Description:
    """Read-only description or section header row in protocol.

    Represents non-interactive informational text displayed in the
    protocol UI (e.g., section headers, explanatory notes, warnings).
    Descriptions are not modifiable parameters but provide user guidance.

    Attributes
    ----------
    text : str
        Display text (e.g., ``'--- Contrast Settings ---'``,
        ``'Sequence supports SMS up to factor 4'``).
    type : str, default 'description'
        Serialization type tag; reserved for protocol encoding.
    """

    text: str
    type: str = 'description'


ProtocolValue = FloatParam | IntParam | BoolParam | StringListParam | Description
Protocol = dict[UIParam | str, ProtocolValue]


# ---------------------------------------------------------------------------
# Bridge serialization helpers (used at Nim <-> Python boundary)
# ---------------------------------------------------------------------------

_TYPE_MAP: dict[str, type] = {
    'float': FloatParam,
    'int': IntParam,
    'bool': BoolParam,
    'stringlist': StringListParam,
    'description': Description,
}


def param_to_dict(p: ProtocolValue) -> dict:
    """Convert a protocol value dataclass to a plain dict.

    Serializes a typed protocol parameter (FloatParam, IntParam, etc.)
    into a dictionary suitable for JSON encoding or transmission across
    language boundaries (Nim/Python bridge).

    Parameters
    ----------
    p : ProtocolValue
        A protocol value dataclass (FloatParam, IntParam, BoolParam,
        StringListParam, or Description).

    Returns
    -------
    dict
        Plain dictionary with all dataclass fields and a ``'type'`` tag
        for deserialization.

    Examples
    --------
    >>> f = FloatParam(value=1.5, min=0.5, max=5.0, incr=0.1, unit='ms')
    >>> d = param_to_dict(f)
    >>> print(d)
    {'value': 1.5, 'min': 0.5, 'max': 5.0, 'incr': 0.1, 'unit': 'ms',
     'validate': 'search', 'type': 'float'}
    """
    return asdict(p)


def dict_to_param(d: dict) -> ProtocolValue:
    """Reconstruct a protocol value dataclass from a plain dict.

    Deserializes a dictionary back into a typed protocol parameter using
    the ``'type'`` field to identify the correct dataclass. The reverse
    operation of :func:`param_to_dict`.

    Parameters
    ----------
    d : dict
        Dictionary with a ``'type'`` field and corresponding value fields.
        Common types: ``'float'``, ``'int'``, ``'bool'``, ``'stringlist'``,
        ``'description'``.

    Returns
    -------
    ProtocolValue
        Reconstructed dataclass instance (FloatParam, IntParam, BoolParam,
        StringListParam, or Description).

    Raises
    ------
    KeyError
        If the ``'type'`` field is missing or does not map to a known class.

    Examples
    --------
    >>> d = {'value': 2.0, 'min': 1.0, 'max': 3.0, 'incr': 0.1,
    ...      'unit': 'ms', 'validate': 'search', 'type': 'float'}
    >>> f = dict_to_param(d)
    >>> assert isinstance(f, FloatParam)
    >>> assert f.value == 2.0
    """
    d = dict(d)  # shallow copy to avoid mutating caller's dict
    tag = d.pop('type')
    cls = _TYPE_MAP[tag]
    return cls(**d)


def protocol_to_dict(protocol: Protocol) -> dict[str, dict]:
    """Serialize an entire protocol to nested plain dicts.

    Converts a complete protocol dictionary (all parameters and descriptions)
    to a JSON-serializable nested dict structure. Useful for persisting
    protocol state or transmitting across process/language boundaries.

    Parameters
    ----------
    protocol : Protocol
        Dictionary mapping parameter keys (UIParam or str) to ProtocolValue
        dataclasses.

    Returns
    -------
    dict[str, dict]
        Nested dict with string keys and dict values. Both keys and values
        are plain Python types suitable for JSON encoding.

    Examples
    --------
    >>> protocol = {
    ...     'TE': FloatParam(value=30.0, min=5.0, max=100.0, incr=1.0, unit='ms'),
    ...     'TR': FloatParam(value=2000.0, min=500.0, max=5000.0, incr=10.0, unit='ms'),
    ...     'FatSat': BoolParam(value=True),
    ... }
    >>> d = protocol_to_dict(protocol)
    >>> # d can now be JSON-serialized: json.dumps(d)
    """
    return {str(k): param_to_dict(v) for k, v in protocol.items()}


def dict_to_protocol(d: dict[str, dict]) -> Protocol:
    """Deserialize nested plain dicts back to a typed Protocol.

    Reconstructs a complete protocol from a JSON-decoded nested dict
    structure. The reverse operation of :func:`protocol_to_dict`.

    Parameters
    ----------
    d : dict[str, dict]
        Nested dict with string keys (parameter names) and dict values
        (parameter specifications with ``'type'`` field).

    Returns
    -------
    Protocol
        Dictionary mapping string keys to ProtocolValue dataclasses.

    Examples
    --------
    >>> nested_dict = {
    ...     'TE': {'value': 30.0, 'min': 5.0, 'max': 100.0, 'incr': 1.0,
    ...            'unit': 'ms', 'validate': 'search', 'type': 'float'},
    ...     'FatSat': {'value': True, 'type': 'bool'},
    ... }
    >>> protocol = dict_to_protocol(nested_dict)
    >>> protocol['TE'].value
    30.0
    """
    return {k: dict_to_param(v) for k, v in d.items()}
