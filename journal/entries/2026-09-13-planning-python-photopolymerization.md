# Planning note: PEGDA photopolymerization Python ladder

**Date:** 2026-09-13  
**Evidence:** planning decision (not a validation campaign)

## Motivation

The repository contained only PDFs. The five unique papers in `Literature/` are complementary, not interchangeable. Without a scope lock, the natural failure mode is a single 3D “Dobson+DLP+flow+mechanics” script that cannot be validated.

## Locked operating mode

- Hydraulics / flow: **off** until Stage 5; then **prescribed plug flow** only.
- Control: **forward** simulation of $p$ (and later $\eta$); no inverse grayscale optimizer.
- Fouling / RNG: **out of scope**.
- Primary claim: staged Python modules can turn the Literature corpus into a testable PEGDA conversion solver.

## Explicit non-goals (this phase)

1. Full Dobson SI chemistry as the default RHS.
2. Coupled flow–mechanics.
3. Scattering, inverse design, PNIPAM, thiol–ene.

## Parked for later (paths only, no implementation)

| Item | When to reopen | Needed inputs |
|------|----------------|---------------|
| Dobson 3D / free volume | Reduced $k_t(p)$ fails two-intensity data | Kinetic series |
| Zhu eosin Y ODEs | Lab uses Type II visible initiation | EY/TEOA spectra |
| Inverse DMD | Stage 4 grayscale `matched` | Target maps |

## What was done

- Read all unique PDFs in `Literature/` (Dobson AFM 2024 + SI, Montgomery EML 2022, Meenakshisundaram JMPT 2020, Slutzky *Soft Matter* 2019, Zhu JMPS 2020).
- Mapped equations to a Python package and campaign stages 0–6.
- Wrote `docs/roadmap_research_paper.md` and `docs/CAMPAIGN_INDEX.md`.
- Aligned with existing Undermind notes; did not execute code.

## What improved

A PhD-readable path from papers → modules → gates, with evidence tags still `proposed`.

## Suggestions to proceed

### Immediate

- [ ] Implement Stage 1 ODE (`solve_ivp` BDF) with the validator table in the roadmap.
- [ ] Parameter YAML with paper provenance.

### Medium

- [ ] Projection engine (Stage 3) then 1D RD (Stage 4).
- [ ] Slutzky plug-flow (Stage 5).

### Long-term

- [ ] Zhu network (Stage 6) and lab `matched` tags.

## Verification

- `Literature/` listing: 7 files, of which two are identical JMPT PDFs and one is Dobson SI.
- No `src/` or tests existed at planning time.

## References

- `docs/roadmap_research_paper.md`
- `docs/CAMPAIGN_INDEX.md`
- https://doi.org/10.1002/adfm.202312607
- https://doi.org/10.1016/j.eml.2022.101714
- https://doi.org/10.1016/j.jmatprotec.2019.116546
- https://doi.org/10.1039/c9sm01485c
- https://doi.org/10.1016/j.jmps.2020.104041
