# Stage 1: local Montgomery four-species ODE

**Date:** 2026-09-13  
**Evidence:** `validated` (internal gates). Not `matched` to FTIR.

## Motivation

Single claim: for prescribed $I(t)$, the Montgomery 2022 four-species balances with $k_p(p)$ and $k_t(p)$ (Eqs. 16–19, 22–27; diffusion off) produce a physically admissible $p(t)$ with oxygen induction and failed equal-dose reciprocity.

## What was done

- Package `src/photopolymerization/` (parameters, closures, RHS, BDF `solve_ivp`).
- `config.pegda_stage1.json`, `run_pegda_stage1.py`, `validate_pegda_stage1.py`, `tests/test_pegda_stage1.py`.
- Outputs: `outputs/pegda_stage1_out/*.csv`, `docs/figures/pegda_stage1_overview.{pdf,png}`.

## What improved

Named gates pass: dark $p=0$; positivity $\min y \approx -3\times 10^{-14}\,\mathrm{mol\,m^{-3}}$; $0\le p\le 1$; PI monotone; acrylate $\int R_p$ residual $\sim 5\times 10^{-5}$; equal dose $I=64\,\mathrm{W\,m^{-2}}$ for $20\,\mathrm{s}$ gives $p=0.740$ vs $I=16$ for $80\,\mathrm{s}$ gives $p=0.800$; $t(p=0.02)$ is $2.11\,\mathrm{s}$ with O₂ vs $0.90\,\mathrm{s}$ without.

This is **not** a match to Montgomery’s FTIR figure. $C_{O,0}$ is assumed.

## Suggestions to proceed

### Immediate
- [ ] Stage 2 constant-rate analytic identities (Slutzky Eq. 3, QSS radicals).
- [ ] Do not start Stage 3–6.

### Medium
- [ ] In-house FTIR equal-dose pair to seek `matched`.

### Long-term
- [ ] Reopen Dobson closures only if equal-dose behavior is wrong after $C_{O,0}$ is measured.

## Verification

```
python3 campaigns/pegda/drivers/validate_pegda_stage1.py --smoke
python3 campaigns/pegda/drivers/validate_pegda_stage1.py
python3 campaigns/pegda/drivers/run_pegda_stage1.py
python3 -m pytest campaigns/pegda/tests -q
```

`--smoke` OVERALL PASS; full OVERALL PASS; 7 passed.

## References

- Montgomery et al. (2022) Eqs. 16–27, Table 2.
- Wydra et al., *Dent. Mater.* **30**, 605 (2014) (reciprocity phenomenology).
- `docs/campaign_roadmap.pdf` Stage 1.
