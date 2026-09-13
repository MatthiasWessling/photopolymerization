# Stage 2: constant-rate analytic identities

**Date:** 2026-09-13  
**Evidence:** `validated` (algebraic identities on the constant-rate fork). Not `matched`.

## Motivation

Single claim: the transcribed RHS contains the limits the papers actually solved — Slutzky Eq. 3 for constant \(k_p,k_O\), and oxygen-free quasi-steady radicals — rather than only an S-shaped \(p(t)\).

## What was done

- Constant-rate fork: `rhs(..., constant_rate=True)` freezes \(k_p=k_{p0}\) and \(k_t=k_t(p=0)\).
- `src/photopolymerization/kinetics/limits.py`
- `config.pegda_stage2.json`, `run_pegda_stage2.py`, `validate_pegda_stage2.py`, `tests/test_pegda_stage2.py`
- Figure: `docs/figures/pegda_stage2_overview.pdf`

## What improved

Slutzky \(C_O/C_{O,0}=(C_M/C_{M,0})^{k_O/k_p}\) holds to \(1.12\times 10^{-6}\) relative residual while oxygen remains (induction window). QSS \(C_R\approx\sqrt{R_i/2k_t}\) residual \(2.1\times 10^{-3}\) after \(0.5\,\mathrm{s}\) with \(C_{O,0}=0\).

**Reasonability:** \(k_O/k_p\approx 1.9\times 10^3\), so oxygen is gone at tiny conversion. Eq. 3 is only tested while \(C_O/C_{O,0}>0.05\). In that window \(p\approx 0\), so the *Montgomery* \(k_p(p)\) model also satisfies Eq. 3 (residual \(1.10\times 10^{-6}\)). That is expected, not a reason to drop the fork. The fork is real at high conversion: \(k_p(0.8)/k_p(0)=7\times 10^{-4}\).

## Suggestions to proceed

### Immediate
- [ ] Stage 3 projection engine only (no chemistry).

### Medium
- [ ] Do not cite Stage 2 as experimental oxygen–monomer data.

### Long-term
- [ ] Revisit Eq. 3 in plug flow (Stage 5) along \(x\) rather than \(t\).

## Verification

```
python3 campaigns/pegda/drivers/validate_pegda_stage2.py --smoke
python3 campaigns/pegda/drivers/validate_pegda_stage2.py
python3 campaigns/pegda/drivers/run_pegda_stage2.py
python3 -m pytest campaigns/pegda/tests -q
```

`--smoke` OVERALL PASS; full OVERALL PASS; 10 passed.

## References

- Slutzky, Stone & Nunes, *Soft Matter* **15**, 9553 (2019), Eq. 3.
- Zhu et al., *J. Mech. Phys. Solids* **142**, 104041 (2020), Eq. 15 analogue.
- Montgomery et al. (2022) conversion-dependent \(k_p(p)\) as the production law, not the identity fork.
- `docs/campaign_roadmap.pdf` Stage 2.
