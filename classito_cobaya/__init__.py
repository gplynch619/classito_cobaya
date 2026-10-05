"""
Cobaya components for classito.

``classito``      -- a drop-in replacement for cobaya's ``classy`` theory that additionally
                    accepts a tabulated sigma_e^2(z) and Delta_c(z) from another Theory.
``ClumpingModel`` -- a template for that other Theory: implement your clumping model in one
                    method and expose its parameters to the sampler.

See ``example.yaml`` for how the two are wired together.
"""

from .classito import classito
from .clumping_template import ClumpingModel

__all__ = ["classito", "ClumpingModel"]
