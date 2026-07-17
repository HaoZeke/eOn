Release the GIL around Matter energy/force accessors and NEB construction,
``update_forces``, and per-image energy reads so metatomic/torch autograd can
run from C++ without deadlocking. Optional ``torch/cuda.h`` and ``torch/mps.h``
includes (``__has_include``) so CPU-only pip torch builds compile without
CUDA/MPS headers.
