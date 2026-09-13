# Photopolymerization

Staged Python solver for **PEGDA free-radical photopolymerization**: a known intensity history (later a DMD field and a flow field) is mapped onto acrylate conversion \(p\), then onto loop-aware gelation and modulus.

This repository is the implementation of the campaign roadmap in `docs/roadmap_research_paper.md`. It is cited from the DFG proposal *Continuous-Flow Digital Lithography of Multi-Porosity Hydrogel Patches* (AVT.CVT, RWTH Aachen).

**Status (2026-09-13):** the module layout and gates are specified; Stage 1 (local stiff ODE) is under implementation. No conversion field has been tagged `validated` yet. Engineering success (a green integrator) is not scientific validity.

## Why modules, not one mixed 3D code

Five papers close different layers of the same process. Mixing them in a single PDE on day one would not tell you which layer is wrong.

| Module | What it computes | Literature counterpart |
|--------|------------------|------------------------|
| Kinetics | Four-species local RHS vs \(I(t)\): PI, radicals, O₂, monomer; \(k_p(p)\), \(k_t(p)\) | Montgomery et al., *Extreme Mech. Lett.* (2022); Dobson & Bowman, *Adv. Funct. Mater.* (2024) as high-fidelity check |
| Projection | Surface irradiance \(I_\perp(x,y,t)\) from bitmaps, pixel PSF, and scan/frame timing | Meenakshisundaram, Sturm & Williams, *J. Mater. Process. Technol.* (2020) |
| Optics / RD | Beer–Lambert attenuation and conversion-dependent diffusion | Montgomery Eqs. 4–23 |
| Flow | Plug-flow oxygen replenishment, then prescribed \(\mathbf{u}\) | Slutzky, Stone & Nunes, *Soft Matter* (2019) |
| Network | Loop-aware \(P_{\mathrm{gel}}\), \(\eta\), \(G\) **after** conversion is trusted | Zhu et al., *J. Mech. Phys. Solids* (2020) |

**Pattern velocity \(\mathbf{v}_{\mathrm{pattern}}\) and resin velocity \(\mathbf{u}\) are independent.** A fluid element sees \(\mathbf{v}_{\mathrm{pattern}}-\mathbf{u}\). Equal dose \(I\times t\) is not equal conversion: bimolecular termination and oxygen inhibition make \(p\) history-dependent.

Dobson’s multi-species free-volume model is a **reference**, not the production right-hand side. Zhu’s eosin-Y bleaching ODEs are not mixed into an Irgacure Type-I rate law. A Jacobs working curve is a diagnostic overlay, not the polymer state.

## Stage ladder

Each stage has a frozen config, a forward run, named pass/fail gates, and one evidence tag (`proposed` / `smoke-tested` / `validated` / `matched`). Do not open the next stage until the gates pass.

0. Parameter book (YAML with paper provenance; PEGDA 250 ≠ 575 ≠ HDDA).
1. **Local ODE** (this is the first code): Montgomery 16–19 with \(D_i=0\), `scipy.integrate.solve_ivp` BDF/Radau.
2. Analytic limits (constant-rate Slutzky/Zhu identities).
3. Projection engine.
4. 1-D then 2-D reaction–diffusion + Beer–Lambert.
5. Slutzky plug-flow benchmark.
6. Zhu network post-processor.

Parked: Navier–Stokes–gel coupling, scattering, inverse DMD design, PNIPAM / thiol–ene closures, Dobson 3-D species set as everyday solver.

## Intended package layout

```text
src/photopolymerization/
  kinetics/     # Stage 1 RHS; same function later used in every spatial cell
  optics/
  projection/
  transport/
  network/
```

Libraries: NumPy and SciPy only for Stages 1–5.

## Cite

```
Wessling, M. Photopolymerization: staged PEGDA reaction–transport solver.
https://github.com/MatthiasWessling/photopolymerization (2026).
```

Publisher PDFs are **not** in this repository (see `Literature/README.md`).

## License

MIT. See `LICENSE`.
