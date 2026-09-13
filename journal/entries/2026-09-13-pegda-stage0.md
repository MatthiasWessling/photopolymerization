# Stage 0: Montgomery 2022 parameter book

**Date:** 2026-09-13  
**Evidence:** `validated` (bookkeeping claim only)

## Motivation

Stage 0 asks whether every coefficient that will enter the Type-I PEGDA-250 ODE has a paper source, a recipe calculation, or an explicit `still_required` flag — before any kinetic fit.

## What was done

- `campaigns/pegda/configs/parameter_book.montgomery_2022.json` (20 entries, Type I, do-not-mix list).
- `config.pegda_stage0.json`, `run_pegda_stage0.py`, `validate_pegda_stage0.py`, `tests/test_pegda_stage0.py`.
- Outputs: `outputs/pegda_stage0_out/provenance_table.csv`, `docs/figures/pegda_stage0_overview.{pdf,png}`.

## What improved

Provenance is machine-checkable. $C_{O,0}$ is flagged `still_required` (assumed $1\,\mathrm{mol\,m^{-3}}$). $\beta$ units are `s^2/kg` for $I$ in W/m².

## Suggestions to proceed

### Immediate
- [ ] Do not use this book for Zhu PEGDA-575 / eosin Y or Slutzky jet numbers.

### Medium
- [ ] Measure dissolved oxygen for in-house resin; retag $C_{O,0}$.

### Long-term
- [ ] Parallel books for Zhu and Slutzky with separate `chemistry_key`s.

## Verification

```
python3 campaigns/pegda/drivers/validate_pegda_stage0.py --smoke
python3 campaigns/pegda/drivers/validate_pegda_stage0.py
python3 campaigns/pegda/drivers/run_pegda_stage0.py
python3 -m pytest campaigns/pegda/tests/test_pegda_stage0.py -q
```

OVERALL PASS (full). pytest included in 7 passed with Stage 1.

## References

- Montgomery et al., *Extreme Mech. Lett.* **53**, 101714 (2022), Table 2 / Table 1.
- `docs/campaign_roadmap.pdf` Stage 0.
