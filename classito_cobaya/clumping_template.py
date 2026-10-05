"""
Template for a small-scale clumping model.

Copy this file, edit :meth:`ClumpingModel.clumping_functions`, and list your parameters in
:meth:`get_can_support_params`. Everything else is wiring you should not need to touch.

The one thing to be careful about: cobaya's Boltzmann components are *parameter-agnostic*, which
means classito will silently absorb any sampled parameter that no other component claims, and CLASS
will then reject it as unknown. So every parameter your model uses must be claimed here, either by
returning it from ``get_can_support_params`` or by declaring a ``params`` block.
"""

import numpy as np

from cobaya.theory import Theory


class ClumpingModel(Theory):
    """Supplies sigma_e^2(z) and Delta_c(z) to classito."""

    def initialize(self):
        # Filled in by must_provide() below, from the grid classito asks for.
        self._z = None

    def must_provide(self, **requirements):
        # Chain to the base, which clears any cached states when the requirements change.
        super().must_provide(**requirements)
        if "ito_clumping" in requirements:
            self._z = np.atleast_1d(requirements["ito_clumping"]["z"])
        return {}

    # ------------------------------------------------------------------ edit me

    def get_can_support_params(self):
        """The sampled parameters your model reads. They must be listed here."""
        return ["clumping_amplitude", "clumping_scale", "clumping_tilt"]

    def clumping_functions(self, z, **params):
        """Return (sigma_e^2(z), Delta_c(z)) on the redshift grid ``z``.

        ``Delta_c`` is the comoving coherence length in Mpc.
        classito forms tau_c(z) = Delta_c(z) * kappa'(z) itself, once thermodynamics has run. This wrapper
        assumes your model's Delta_c can be calculated without the thermodynamic state of the system. A coupled
        evolution requires directly modifying the code.

        Both functions are held at their endpoint values outside the requested grid, so if your
        model should switch off there, return zeros at the ends.

        Mind the dynamic range. The grid spans ito_z_min to ito_z_max, three decades in 1+z by
        default, so an unbounded power law pivoted at recombination varies over many decades across
        it and will leave the tabulated sigma_e^2 range of [0, 2]. So keep sigma_e^2 bounded, or
        narrow ito_z_min/ito_z_max to the range your model is really about.

        The placeholder below is a bounded lognormal bump in 1+z centred on recombination, with a
        power-law coherence length. Replace it with your model.
        """
        amplitude = params["clumping_amplitude"]
        scale = params["clumping_scale"]
        tilt = params["clumping_tilt"]

        x = np.log((1. + z) / 1101.)
        sigma_e2 = amplitude * np.exp(-0.5 * (x / 1.5) ** 2)
        Delta_c = scale * np.exp(tilt * x)
        return sigma_e2, Delta_c

    # ------------------------------------------------------------------ wiring

    def calculate(self, state, want_derived=True, **params_values_dict):
        sigma_e2, Delta_c = self.clumping_functions(self._z, **params_values_dict)
        state["ito_clumping"] = {"z": self._z,
                                 "sigma_e2": np.asarray(sigma_e2, dtype=float),
                                 "Delta_c": np.asarray(Delta_c, dtype=float)}

    def get_ito_clumping(self):
        # Defining this method is what registers "ito_clumping" as a product this component can
        # provide: cobaya discovers providable quantities by introspecting get_* methods.
        return self.current_state["ito_clumping"]
