# Stage 1 PEGDA solver implementation

##### [**Undermind**](https://undermind.ai)

---


## Table of Contents

- [Stage 1 PEGDA solver implementation](#stage-1-pegda-solver-implementation)
  - [Purpose](#purpose)
  - [Solver boundary](#solver-boundary)
  - [Suggested implementation layout](#suggested-implementation-layout)
  - [State and parameter conventions](#state-and-parameter-conventions)
  - [Core equations](#core-equations)
  - [Kinetic closures](#kinetic-closures)
  - [Parameter file template](#parameter-file-template)
  - [Function interfaces](#function-interfaces)
    - [Intensity protocol](#intensity-protocol)
    - [Closure functions](#closure-functions)
    - [ODE right-hand side](#ode-right-hand-side)
  - [Positivity and diagnostics](#positivity-and-diagnostics)
  - [Parameter-identification workflow](#parameter-identification-workflow)
    - [1. Photoinitiator block](#photoinitiator-block)
    - [2. Oxygen block](#oxygen-block)
    - [3. Intensity-history block](#intensity-history-block)
    - [4. Thermal block](#thermal-block)
    - [5. Network handoff](#network-handoff)
  - [Test protocols](#test-protocols)
  - [Output schema](#output-schema)
  - [Scope boundary](#scope-boundary)
  - [References](#references)

# Stage 1 PEGDA solver implementation

## Purpose

This is a clean-room implementation specification for the Stage 1 local PEGDA kinetics model. It is not recovered author code: the papers provide equations and parameter information but no public repository was found. The implementation should reproduce local conversion, radical, oxygen, and temperature histories before any spatial or flow coupling is added.

The model combines the reduced spatial chemistry of \[Mon22\] with the free-volume, reaction-diffusion, and thermal ideas in \[Dob24\]. PEGDA network closure is deferred to \[Zhu20\].

## Solver boundary

Stage 1 solves a local initial-value problem driven by an arbitrary intensity history:

``` math
\mathbf{y}(t)=\left(C_I,C_R,C_O,C_M,T\right).
```

It does not yet include:

- spatial diffusion or advection,
- DMD pixel overlap or measured PSFs,
- moving-mask kinematics,
- network topology,
- swelling or mechanics.

The later spatial solver will call the same local right-hand side at every grid cell or material trajectory.

## Suggested implementation layout

``` text
stage1_pegda/
├── README.md
├── parameters.yaml
├── src/
│   ├── state.py          # state vector and units
│   ├── closures.py       # kp, kt, temperature, diffusivity closures
│   ├── protocols.py      # constant, pulse, ramp, and measured I(t)
│   ├── kinetics.py       # ODE right-hand side
│   ├── solve.py          # BDF/Radau integration and diagnostics
│   ├── fit.py            # staged parameter identification
│   └── observables.py    # conversion, rate, oxygen, gel-input outputs
├── tests/
│   ├── test_limits.py
│   ├── test_reciprocity.py
│   └── test_mass_balance.py
├── data/
│   ├── raw/
│   └── processed/
└── results/
```

The implementation language can be Python, MATLAB, Julia, or another environment with a stiff implicit ODE solver. The equations below are language-independent.

## State and parameter conventions

Use SI units internally:

- concentration: $`\mathrm{mol\,m^{-3}}`$,
- time: $`\mathrm{s}`$,
- intensity: $`\mathrm{W\,m^{-2}}`$,
- diffusivity: $`\mathrm{m^2\,s^{-1}}`$,
- temperature: $`\mathrm{K}`$.

The photolysis coefficient $`\beta`$ from \[Mon22\] is tied to that paper’s intensity convention. It must be recalibrated if intensity is converted to SI units.

Use `conversion` for the degree of conversion and reserve `A` for optical attenuation.

## Core equations

Given state $`y=(C_I,C_R,C_O,C_M,T)`$ and prescribed $`I(t)`$, calculate

``` math
p=1-\frac{C_M}{C_{M,0}}.
```

Photoinitiator consumption:

``` math
\dot C_I=-\beta I C_I.
```

Radical generation:

``` math
R_{\mathrm{init}}=m\beta I C_I.
```

Radical balance:

``` math
\dot C_R=R_{\mathrm{init}}-2k_t(p,T)C_R^2-k_OC_OC_R.
```

Oxygen balance:

``` math
\dot C_O=-k_OC_OC_R.
```

Acrylate consumption:

``` math
\dot C_M=-k_p(p,T)C_RC_M.
```

For the initial implementation, use $`T=T_0`$. The optional thermal equation is

``` math
\rho c_p\dot T=\Delta H_p k_pC_RC_M-h_aa_s(T-T_a)+q_{\mathrm{abs}}.
```

## Kinetic closures

Use the \[Mon22\] conversion-dependent propagation law:

``` math
k_p(p,T_0)=
\frac{k_{p0}}
{1+\dfrac{k_{p0}}{k_{p,D0}}\exp(c_pp)}.
```

Use the combined translational and reaction-diffusion termination law:

``` math
k_t(p,T_0)=
\left[
\frac{1}{k_{t,SD}}+\frac{\exp(c_pp)}{k_{t,TD0}}
\right]^{-1}
+\frac{C_{RD}(1-p)k_{p0}}
{1+\dfrac{k_{p0}}{k_{p,D0}}\exp(c_pp)}.
```

The \[Dob24\] high-fidelity model should be used to test whether these reduced closures remain valid over the intended intensity, conversion, and temperature range. If not, replace the closures with free-volume or chain-length-dependent forms rather than refitting a single global exponent.

Optional temperature corrections are

``` math
k_j(p,T)=k_j(p,T_0)
\exp\left[-\frac{E_j}{R}
\left(\frac{1}{T}-\frac{1}{T_0}\right)\right],
\qquad j\in\{p,t\}.
```

## Parameter file template

``` yaml
chemistry:
  m: 2.0
  beta: 2.70e-3       # source intensity convention; recalibrate
  k_oxygen: 3.50e3   # m^3 mol^-1 s^-1
  kp0: 1.86           # m^3 mol^-1 s^-1
  kp_diffusion: 8.99e8
  conversion_sensitivity: 34.15
  kt_translational: 4.39e3
  kt_transport: 1.00e4
  reaction_diffusion_factor: 1.01

formulation:
  monomer_concentration: null       # calculate from PEGDA recipe
  photoinitiator_concentration: null
  oxygen_concentration: null         # measure after equilibration

thermal:
  enabled: false
  initial_temperature: null
  ambient_temperature: null
  volumetric_heat_capacity: null
  heat_of_polymerization: null
  lumped_heat_loss: null
  propagation_activation_energy: null
  termination_activation_energy: null

solver:
  method: BDF
  relative_tolerance: 1.0e-7
  absolute_tolerance: 1.0e-10
  max_step: null
```

The numerical values for the kinetic block are starting values reported or calibrated in \[Mon22\]. They are not expected to transfer unchanged to another PEGDA molecular weight, photoinitiator, wavelength, oxygen level, or solvent composition.

## Function interfaces

### Intensity protocol

``` text
intensity(t, protocol) -> I(t)
```

Required protocol types:

- constant intensity,
- square pulse,
- pulse train,
- ramp,
- measured photodiode trace,
- sampled DMD/material-element history.

### Closure functions

``` text
conversion(CM, CM0) -> p
kp(p, T, params) -> kp_value
kt(p, T, params) -> kt_value
reaction_rate(CR, CM, kp_value) -> Rp
```

### ODE right-hand side

``` text
rhs(t, y, protocol, params):
    CI, CR, CO, CM, T = y
    I = intensity(t, protocol)
    p = 1 - CM / CM0
    kp_value = kp(p, T, params)
    kt_value = kt(p, T, params)

    Rinit = m * beta * I * CI
    dCI = -beta * I * CI
    dCR = Rinit - 2 * kt_value * CR**2 - k_oxygen * CO * CR
    dCO = -k_oxygen * CO * CR
    dCM = -kp_value * CR * CM

    if thermal_enabled:
        dT = thermal_rhs(...)
    else:
        dT = 0

    return [dCI, dCR, dCO, dCM, dT]
```

Use an implicit stiff integrator. Direct explicit integration will generally fail or require impractically small time steps during the radical transient.

## Positivity and diagnostics

The solver should monitor:

- $`C_I,C_R,C_O,C_M\ge0`$,
- $`0\le p\le1`$,
- monotonic decrease of $`C_I`$ and $`C_M`$,
- oxygen consumption when $`I>0`$,
- radical balance residual,
- time-integrated reaction rate.

If the integrator produces small negative concentrations from numerical error, use a positivity-preserving tolerance strategy or integrate logarithmic concentrations. Do not silently clip large negative values.

## Parameter-identification workflow

### 1. Photoinitiator block

Fit $`\beta C_{I,0}`$ from photoinitiator absorbance decay or low-conversion experiments. Avoid fitting $`\beta`$ and $`C_{I,0}`$ independently unless the initial initiator concentration is independently known.

### 2. Oxygen block

Fit $`k_O`$ and validate $`C_{O,0}`$ using induction time, inhibition-layer thickness, or oxygen-controlled curing experiments. Flow experiments should be reserved for the later \[Slu19\] benchmark.

### 3. Intensity-history block

Use several pairs of exposures with equal nominal dose but different intensity. Fit $`k_{p0}`$, $`k_{t,SD}`$, $`k_{t,TD0}`$, and $`C_{RD}`$ against conversion-versus-time data. A single dose curve cannot identify these parameters reliably.

### 4. Thermal block

Activate temperature only after the isothermal model fails systematically or measured temperature rises are significant. Fit heat generation and heat loss against temperature traces and late-time conversion.

### 5. Network handoff

Pass $`p(t)`$ and, later, local functional-group histories to the \[Zhu20\] network closure. Do not fit modulus directly into the Stage 1 reaction rate unless the network model is intentionally being replaced by an empirical surrogate.

## Test protocols

| Test | Input | Expected diagnostic |
|:---|:---|:---|
| Dark control | $`I(t)=0`$ | No photoinitiator-driven conversion |
| Oxygen-free limit | $`C_{O,0}=0`$ | No oxygen induction period |
| Initiator sweep | Vary $`C_{I,0}`$ | Predictable change in initiation rate |
| Equal-dose sweep | $`I_1t_1=I_2t_2`$ | Conversion curves do not collapse in general \[Wyd14, Lin19\] |
| High-intensity pulse | Short, bright exposure | Strong radical termination and history dependence |
| Low-intensity exposure | Long, dim exposure | Different conversion at the same dose |
| Thermal comparison | Isothermal versus thermal mode | Divergence only when heat matters |
| Parameter transfer | New PEGDA formulation | Exposes which parameters are formulation-specific |

## Output schema

Each simulation should export a time series with at least:

``` text
time
intensity
photoinitiator_concentration
radical_concentration
oxygen_concentration
monomer_concentration
conversion
propagation_rate
termination_rate
oxygen_inhibition_rate
temperature
```

These outputs become the local constitutive response used later by the DMD projection, reaction–diffusion, and Eulerian–Lagrangian flow modules.

## Scope boundary

This implementation is deliberately smaller than the full \[Dob24\] model. Its role is to establish a stable and testable local constitutive law. The next upgrades should be made in this order:

1.  add spatial diffusion and dynamic Beer–Lambert attenuation,
2.  insert the measured moving-DMD projection field from \[Mee20\],
3.  add prescribed flow and material trajectories using \[Slu19\] scaling,
4.  add the \[Zhu20\] network closure,
5.  add thermal feedback, conversion-dependent viscosity, and mechanics.

---

## References

\[Mon22\] S. Montgomery, C. M. Hamel, J. Skovran, and H. Q, “A reaction-diffusion model for grayscale digital light processing 3D printing,” Mar. 01, 2022. doi: [10.1016/j.eml.2022.101714](https://doi.org/10.1016/j.eml.2022.101714).

\[Dob24\] A. L. Dobson and C. N. Bowman, “A Comprehensive, Multidimensional First‐Principles Model for Free‐Radical Photopolymerizations in Bulk and Thin Films,” *Advanced Functional Materials*, vol. 34, Feb. 2024, doi: [10.1002/adfm.202312607](https://doi.org/10.1002/adfm.202312607).

\[Zhu20\] H. Zhu, X. Yang, G. Genin, T. Lu, F. Xu, and M. Lin, “Modeling the mechanics, kinetics, and network evolution of photopolymerized hydrogels,” Sep. 01, 2020. doi: [10.1016/j.jmps.2020.104041](https://doi.org/10.1016/j.jmps.2020.104041).

\[Slu19\] M. Slutzky, H. Stone, and J. Nunes, “A quantitative study of the effect of flow on the photopolymerization of fibers.” *Soft matter*, Nov. 2019, doi: [10.1039/c9sm01485c](https://doi.org/10.1039/c9sm01485c).

\[Wyd14\] J. W. Wydra, N. Cramer, J. Stansbury, and C. Bowman, “The reciprocity law concerning light dose–relationships applied to BisGMA/TEGDMA photopolymers: Theoretical analysis and experimental characterization,” *Dental materials : official publication of the Academy of Dental Materials*, vol. 30, pp. 605–612, Mar. 2014, doi: [10.1016/j.dental.2014.02.021](https://doi.org/10.1016/j.dental.2014.02.021).

\[Lin19\] J.-T. Lin, H. Liu, K.-T. Chen, and D.-C. Cheng, “Modeling the Kinetics, Curing Depth, and Efficacy of Radical-Mediated Photopolymerization: The Role of Oxygen Inhibition, Viscosity, and Dynamic Light Intensity,” *Frontiers in Chemistry*, vol. 7, Nov. 2019, doi: [10.3389/fchem.2019.00760](https://doi.org/10.3389/fchem.2019.00760).

\[Mee20\] V. Meenakshisundaram, L. Sturm, and C. Williams, “Modeling A Scanning-Mask Projection Vat Photopolymerization System For Multiscale Additive Manufacturing,” May 01, 2020. doi: [10.1016/j.jmatprotec.2019.116546](https://doi.org/10.1016/j.jmatprotec.2019.116546).
