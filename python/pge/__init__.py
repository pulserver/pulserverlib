"""
Pulserver — offline validation and debugging for pulseq MRI sequences.

Public API:

    SequenceCollection    Wraps a pypulseq Sequence with C-backed analysis.
    serialize             Save collection to linked .seq file chain.
    deserialize           Restore collection from linked .seq file chain.
"""

__all__ = [
    'SequenceCollection',
    'deserialize',
    'serialize',
]

from .core import (
    SequenceCollection,
    deserialize,
    serialize,
)
