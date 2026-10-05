"""
A cobaya Theory wrapping classito.

Cobaya's ``classy`` component is a pure producer: it never reads from ``self.provider``, so it
cannot accept a quantity computed by another Theory. This subclass adds exactly that for the one
product classito needs, following the pattern CAMB uses for ``external_primordial_pk``.

This subclass inherits from Cobaya's ``classy`` component, and accepts a z,sigma_e^2(z),Delta_c(z) table
produced by a companion theory code.

The dependency runs

    likelihood  --requires-->  C_l
    classito    --requires-->  ito_clumping   (this class)
    your model  --provides-->  ito_clumping   (see clumping_template.py)

so the sampler's parameters reach CLASS through your model rather than through classito.
"""

import hashlib
import os
import shutil
import tempfile

import numpy as np

from cobaya.theories.classy import classy
from cobaya.log import LoggedError

# Every CLASS parameter value is copied into a fixed-size buffer (_ARGUMENT_LENGTH_MAX_ = 1024
# bytes). Stay clear of it, and fall back to a file when a table does not fit.
_MAX_VALUE_BYTES = 1000

# The product this component asks another Theory for.
_CLUMPING = "ito_clumping"


class classito(classy):
    """classy plus a tabulated sigma_e^2(z) and Delta_c(z) supplied by another Theory."""

    # Nodes on which the clumping functions are requested. They are log-spaced in 1+z, and the
    # defaults span recombination with room on either side. Both functions are held at their
    # endpoint values outside this range, so widen it if your model is still active there.
    ito_nodes: int = 48
    ito_z_min: float = 10.
    ito_z_max: float = 2.e4

    # "auto"   -- inline lists, falling back to a file if they would overflow
    # "inline" -- inline lists only, error if they would overflow
    # "file"   -- always write a file
    ito_transport: str = "auto"

    def initialize(self):
        super().initialize()
        if self.ito_transport not in ("auto", "inline", "file"):
            raise LoggedError(self.log, "ito_transport must be 'auto', 'inline' or 'file', not %r",
                              self.ito_transport)
        if self.ito_nodes < 4:
            raise LoggedError(self.log, "ito_nodes must be at least 4 to spline, got %d",
                              self.ito_nodes)
        self._z_nodes = np.geomspace(1. + self.ito_z_min,
                                     1. + self.ito_z_max, self.ito_nodes) - 1.
        self._wants_clumping = self._detect_clumping()
        self._ito_args = {}
        # One directory per instance keeps concurrent MPI ranks from writing over each other.
        self._table_dir = None

    def close(self):
        if self._table_dir is not None:
            shutil.rmtree(self._table_dir, ignore_errors=True)
            self._table_dir = None
        super().close()

    def _detect_clumping(self):
        """Whether the ito module is switched on and in table mode.

        Decided once, at initialisation, and cached. Reading it out of extra_args on every call
        would be fragile: set() below temporarily writes the table into extra_args, so the flag
        would depend on bookkeeping that has nothing to do with the user's configuration.
        """
        return (str(self.extra_args.get("has_ito", "no")).lower().startswith("y") and
                str(self.extra_args.get("ito_input_mode", "scalar")) == "table")

    def must_provide(self, **requirements):
        # BoltzmannBase raises on products it does not recognise, so let it see the request first
        # and only then add our own conditional requirement.
        super().must_provide(**requirements)
        if self._wants_clumping:
            return {_CLUMPING: {"z": self._z_nodes}}
        return None

    def calculate(self, state, want_derived=True, **params_values_dict):
        if self._wants_clumping:
            self._ito_args = self._encode(self.provider.get_ito_clumping())
        return super().calculate(state, want_derived, **params_values_dict)

    def set(self, params_values_dict):
        # classy.set is the single funnel into CLASS, so the table is merged here. Only the keys we
        # added are removed afterwards, rather than restoring a snapshot: classy.set mutates
        # extra_args itself, and those changes should survive. Removing them at all matters because
        # the transport can switch between inline lists and a file from one call to the next, and a
        # stale key from the other route would otherwise linger.
        previous = {k: self.extra_args[k] for k in self._ito_args if k in self.extra_args}
        introduced = [k for k in self._ito_args if k not in self.extra_args]
        self.extra_args.update(self._ito_args)
        try:
            super().set(params_values_dict)
        finally:
            # Put back what was there (ito_input_mode comes from the user's yaml) and drop only the
            # keys we introduced, so a switch between the inline and file routes leaves nothing stale.
            self.extra_args.update(previous)
            for key in introduced:
                self.extra_args.pop(key, None)

    # ------------------------------------------------------------------ helpers

    def _encode(self, clumping):
        """Turn the provider's arrays into CLASS input parameters."""
        z, var_e, Delta_c = self._validate(clumping)

        if self.ito_transport != "file":
            args = {"ito_input_mode": "table",
                    "ito_z_table": _as_list(z),
                    "ito_var_e_table": _as_list(var_e),
                    "ito_Delta_c_table": _as_list(Delta_c)}
            longest = max(len(v.encode()) for v in list(args.values())[1:])
            if longest <= _MAX_VALUE_BYTES:
                return args
            if self.ito_transport == "inline":
                raise LoggedError(
                    self.log,
                    "The clumping table needs %d bytes per parameter but CLASS accepts at most "
                    "%d. Reduce ito_nodes (currently %d) or set ito_transport to 'auto' or "
                    "'file'.", longest, _MAX_VALUE_BYTES, self.ito_nodes)

        return {"ito_input_mode": "table", "ito_table_file": self._write(z, var_e, Delta_c)}

    def _validate(self, clumping):
        try:
            z = np.atleast_1d(np.asarray(clumping["z"], dtype=float))
            var_e = np.atleast_1d(np.asarray(clumping["sigma_e2"], dtype=float))
            Delta_c = np.atleast_1d(np.asarray(clumping["Delta_c"], dtype=float))
        except (TypeError, KeyError, ValueError) as excpt:
            raise LoggedError(
                self.log, "The component providing %r must return a dict with keys 'z', "
                          "'sigma_e2' and 'Delta_c'. Got %r (%s)", _CLUMPING, clumping, excpt)

        if not (z.shape == var_e.shape == Delta_c.shape):
            raise LoggedError(
                self.log, "'z', 'sigma_e2' and 'Delta_c' must have the same shape, got %r, %r, %r",
                z.shape, var_e.shape, Delta_c.shape)
        if not np.all(np.isfinite(var_e)) or not np.all(np.isfinite(Delta_c)):
            raise LoggedError(self.log, "Your clumping model returned a non-finite sigma_e^2 or "
                                        "Delta_c. CLASS cannot use it.")
        if np.any(var_e < 0) or np.any(Delta_c < 0):
            raise LoggedError(self.log, "sigma_e^2 and Delta_c must be non-negative; got minima "
                                        "%g and %g.", var_e.min(), Delta_c.min())
        return z, var_e, Delta_c

    def _write(self, z, var_e, Delta_c):
        """Write the table and return its path, named by a hash of its contents.

        The name has to change when the numbers do. CLASS's python wrapper skips recomputing when
        the parameter dictionary it is handed is unchanged, and a fixed filename looks unchanged
        however different the file behind it is -- which would silently freeze the theory at the
        first point of a chain. Hashing the contents also makes the converse true: an identical
        table reuses the same path and is correctly served from cache.
        """
        if self._table_dir is None:
            self._table_dir = tempfile.mkdtemp(prefix="classito_")

        body = "%d\n" % len(z) + "".join(
            "%s %s %s\n" % tuple(_fmt(v) for v in row) for row in zip(z, var_e, Delta_c))
        digest = hashlib.sha1(body.encode()).hexdigest()[:16]
        path = os.path.join(self._table_dir, "ito_table_%s.dat" % digest)

        if not os.path.exists(path):
            with open(path, "w") as table_file:
                table_file.write("# z  sigma_e^2  Delta_c[Mpc], from the classito cobaya wrapper\n")
                table_file.write(body)
            self._prune(keep=path)
        return path

    def _prune(self, keep, limit=8):
        """Keep the directory from growing without bound over a long chain."""
        try:
            files = [os.path.join(self._table_dir, f) for f in os.listdir(self._table_dir)]
            if len(files) <= limit:
                return
            for stale in sorted(files, key=os.path.getmtime)[:len(files) - limit]:
                if stale != keep:
                    os.remove(stale)
        except OSError:
            pass  # a full cache directory is not worth failing a likelihood evaluation over


def _fmt(value):
    """Format one number so that it parses back to exactly the same double.

    repr() gives the shortest string that round-trips, which keeps the inline lists compact
    without losing anything. Both transports use it, so that whether a table travels inline or
    through a file - a choice ito_transport: auto makes silently, based on length - cannot change
    the numbers CLASS receives.
    """
    return repr(float(value))


def _as_list(values):
    """Format an array as the comma-separated list CLASS's parser expects."""
    return ",".join(_fmt(v) for v in values)
