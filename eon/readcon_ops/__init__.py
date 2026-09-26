"""Operations on structures. Not the CON codec.

``readcon`` parses and writes ``ConFrame``. ``readcon-db`` stores a corpus
of those frames. This package is the third piece: comparisons and rigid
alignment of two frames. Neighbor lists stay in :mod:`eon.geometry` until
the same functions can take a ``ConFrame`` without an eOn ``Structure``.

eOn keeps :mod:`eon.atoms` as the compatibility surface. New callers should
import from here.
"""

from eon.readcon_ops.match import identical, internal_motion, rotate

__all__ = ["identical", "internal_motion", "rotate"]
