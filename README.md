# classito-cobaya

Cobaya interface for [classito](https://github.com/gplynch619/classito), the CLASS fork that
models small-scale electron density fluctuations (Chluba, Vasil & Battye 2025,
[arXiv:2505.22242](https://arxiv.org/abs/2505.22242)).

If you have a physical clumping model (e.g. primordial magnetic fields, a phase transition,
enhanced small-scale power) that predicts a variance $\sigma_e^2(z)$ and a comoving coherence
length $\Delta_c(z)$, this package lets you sample its parameters in an MCMC by writing one
Python method.

## How it fits together

This wrapper is used with classito:

- **`classito_cobaya.classito`** — a `classy` subclass that additionally requires a product
  called `ito_clumping` and forwards it into classito's table mode
  (`ito_input_mode = table`).
- **`classito_cobaya.clumping_template.ClumpingModel`** — a template cobaya
  `Theory` that provides `ito_clumping`.

```
likelihood   --requires-->  C_l
classito     --requires-->  ito_clumping
your model   --provides-->  ito_clumping
```

## Install

classito and this package are separate products; you need both.

```bash
# 1. classito itself (builds the classy python wrapper from the fork)
git clone https://github.com/gplynch619/classito
cd classito && make -j

# 2. this package
git clone https://github.com/gplynch619/classito-cobaya
cd classito-cobaya && pip install .
```

## Quick start

Copy `classito_cobaya/clumping_template.py`, rename the class, and edit two methods:

```python
def get_can_support_params(self):
    return ["my_amplitude", "my_scale"]        # every parameter your model reads

def clumping_functions(self, z, **params):
    sigma_e2 = ...                             # your model, evaluated on z
    Delta_c  = ...                             # comoving Mpc
    return sigma_e2, Delta_c
```

Then wire it up in yaml (full example in `examples/example.yaml`):

```yaml
theory:
  classito_cobaya.classito:
    extra_args: {has_ito: yes, ito_input_mode: table, electron_dist: lognormal,
                 gauge: newtonian, radiation_streaming_approximation: 3,
                 tight_coupling_trigger_tau_c_over_tau_h: 0.005,
                 tight_coupling_trigger_tau_c_over_tau_k: 0.008}
    ito_nodes: 48          # grid the clumping functions are requested on,
    ito_z_min: 10.         # log-spaced in 1+z
    ito_z_max: 2.e4
  my_module.MyClumpingModel:
    provides: [ito_clumping]
```

The `extra_args` block above needs to use the required classito settings. This is discussed in Appendix B of the classito paper.

For a worked example without a sampler, `examples/clumping_model_example.ipynb` defines a model,
evaluates it at a point through Cobaya, and plots the resulting $C_\ell$ against a fiducial ΛCDM.

## Options

| option | default | meaning |
|---|---|---|
| `ito_nodes` | 48 | number of redshift nodes requested from your model |
| `ito_z_min`, `ito_z_max` | 10, 2e4 | range of that grid |
| `ito_transport` | `auto` | `inline` passes the table as CLASS parameters; `file` writes it to a per-process temp file; `auto` uses inline when it fits (CLASS caps each parameter at 1024 bytes) and falls back to the file. All three give identical results. |

## Citation

Please cite the CLASS papers, Chluba, Vasil & Battye (arXiv:2505.22242), and the classito
release paper.
