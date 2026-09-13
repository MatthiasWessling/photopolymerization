# Stage 4: 1-D reaction–diffusion and Beer–Lambert

**Date:** 2026-09-13  
**Evidence:** `validated` (internal gates). Not `matched` to Montgomery printed widths or Dendukuri PDMS films.

## Motivation

Single claim: the Stage 1 four-species rates, evaluated on a 1-D method-of-lines grid with Montgomery Beer–Lambert \(A(p,C_I)\) and harmonic \(D_i(p)\), produce (i) a conversion front when the film is optically thick, (ii) an oxygen dead zone at a Dirichlet wall, and (iii) a G0|G80 conversion transition wider than the optical \(\sigma\). Flow is off. An \(x\)–\(z\) print grid is not required for these gates.

## What was done

- Parameter book: Montgomery Table 1 optics and Table 3 diffusivities (`alpha_I`, \(A_{\mathrm{mon}}\), \(A_{\mathrm{poly}}\), \(A_{\mathrm{abs}}\), \(D_i^{\mathrm{liquid/solid}}\)).
- `src/photopolymerization/rd/`: Beer–Lambert (Eqs. 4–6), harmonic \(D(p)\) (Eq. 28), conservative 1-D diffusion, MOL RHS that calls `reaction_rates` from Stage 1, cosine MMS, BDF integrator.
- `config.pegda_stage4.json`, `run_pegda_stage4.py`, `validate_pegda_stage4.py`, `tests/test_pegda_stage4.py`.
- Figure: `docs/figures/pegda_stage4_overview.pdf`.
- Geometry: \(H=200\,\mu\mathrm{m}\) labelled film (not an FTIR path length). Thick-vat absorber \(w=0.06\) wt% as in Montgomery’s multilayer example. Sealed Neumann vs \(C_O\) Dirichlet at the illuminated node.

## What improved

| Gate | Result |
|------|--------|
| Explicit Fourier \(dt\) documented | \(0.125\,\mathrm{s}\) at \(n_z=41\); integrator is implicit BDF |
| Dark sealed \(\int C_M\) | relative residual \(0\) |
| Cosine MMS (constant \(D_I^{\mathrm{liq}}\)) | rel. max error \(6.9\times 10^{-3}\) |
| Surface cell \(D_i=0\) vs Stage 1 ODE | \(\lvert\Delta p\rvert\sim 4\times 10^{-6}\) at \(8\,\mathrm{s}\), \(I=64\,\mathrm{W\,m^{-2}}\) |
| Thick-vat front | \(\beta_{\mathrm{opt}}=1.50\), \(p(0)-p(H)=0.52\) |
| O2 Dirichlet dead zone | \(p_{\mathrm{wall}}=0.0069\), \(p_{\mathrm{mid}}=0.53\) (sealed face \(0.55\)) |
| G0\|G80 10–90% width | \(53\,\mu\mathrm{m} > \sigma=17.7\,\mu\mathrm{m}\) |
| Table 3 \(D_i\) increment | \(w_{\mathrm{on}}-w_{\mathrm{off}}=0.21\,\mu\mathrm{m}\) at \(8\,\mathrm{s}\) (logged, not a fail) |

**Reasonability:** radical/monomer diffusion at Table 3 values does not add a Montgomery-scale extra blur on an \(8\,\mathrm{s}\) strip. The 2–3 pixel interface in their Fig. 7 is consistent with the Gaussian PSF (Stage 3). \(C_{O,0}\) remains `still_required`. \(w\) in Eq. 6 is treated as weight percent following the paper’s wording. This stage is **not** an \(x\)–\(z\) DLP build.

## Suggestions to proceed

### Immediate
- [ ] Stage 5 plug flow only (Slutzky). Do not add Zhu \(P_{\mathrm{gel}}\).
- [ ] Do not cite Stage 4 grayscale widths as evidence that \(D_R\) sets feature size.

### Medium
- [ ] Coupled \(x\)–\(z\) strip once a print layer height is fixed.
- [ ] Refine near the oxygen wall if inhibition-layer thickness is a `matched` target.

### Long-term
- [ ] Dobson free-volume \(D_i\) only if Table 3 plus FTIR depth profiles fail.

## Verification

```
python3 campaigns/pegda/drivers/validate_pegda_stage4.py --smoke
python3 campaigns/pegda/drivers/validate_pegda_stage4.py
python3 campaigns/pegda/drivers/run_pegda_stage4.py
python3 -m pytest campaigns/pegda/tests -q
```

OVERALL PASS (smoke and full). 19 tests passed.

## References

- Montgomery et al. (2022) Eqs. 4–6, 15–19, 28; Tables 1 and 3.
- Dendukuri et al. (2008) for the oxygen-wall variant (not a geometry match).
- `docs/campaign_roadmap.pdf` Stage 4.
- Prior: `journal/entries/2026-09-13-pegda-stage3.md`.
