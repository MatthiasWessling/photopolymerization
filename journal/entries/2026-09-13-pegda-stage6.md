# Stage 6: Zhu network post-processor

**Date:** 2026-09-13  
**Evidence:** `validated` (algebraic identities). Not `matched` to Zhu shear-modulus or indentation data.

## Motivation

Single claim: given acrylate conversion \(p\) and \([\mathrm{PEGDA}]_0\), Zhu’s loop-aware Macosko–Miller map yields \(P_{\mathrm{gel}}\), \(\eta\le 1\), and \(G<G_{\mathrm{ideal}}\), with no gel when the precursor is at or below \(C_0\). This stage does **not** integrate species and does **not** run eosin Y bleaching ODEs. It does not retune Montgomery \(k_p\).

## What was done

- Book: `campaigns/pegda/configs/parameter_book.zhu_2020.json` (\(C_0=3.45\) vol%, Table A2; Type II flag; bleaching ODEs off).
- `src/photopolymerization/network/macosko.py`: Eqs. 21, 26, 32–38; no-loop Eq. 36 \(\eta=p^2\).
- Applied \(\eta(p)\) to a Stage 1 Montgomery trajectory **labelled hypothetical**.
- `config.pegda_stage6.json`, `run_pegda_stage6.py`, `validate_pegda_stage6.py`, `tests/test_pegda_stage6.py`.
- Figure: `docs/figures/pegda_stage6_overview.pdf`.

## What improved

| Gate | Result |
|------|--------|
| \(P_{\mathrm{gel}}=C_0/c\) at 20 vol% | \(0.1725\) |
| \(c\le C_0\) | \(\eta(p=1)=0\) |
| \(\theta\to 1\) | \(P_{\mathrm{gel}}\to 0\) |
| No-loop \(\eta=p^2\) | identity |
| Loops at 20 vol%, \(p=1\) | \(\eta=0.378<1\), \(G/G_{\mathrm{ideal}}=0.378\) |
| Pregel | \(\eta=0\) for \(p<P_{\mathrm{gel}}\) |
| Not Slutzky 2% | \(0.173\) vs \(0.02\) |

**Reasonability:** \(C_0\) is Zhu’s **fitted** Table A2 value, not a universal PEGDA constant. \(G\) uses an assumed PEGDA density to convert vol% to mol/m³; ratios \(G/G_{\mathrm{ideal}}=\eta\) do not depend on that density. Mapping Irgacure Type I \(p(t)\) onto a Type II hydrogel \(\eta\) is a thought experiment.

## Suggestions to proceed

### Immediate
- [ ] Do not open inverse DMD or swelling FEM without a new request.
- [ ] Do not cite Stage 6 \(G\) as Montgomery FTIR evidence.

### Medium
- [ ] Overlay Zhu \(G(p)\) / indentation for a `matched` retag on **this** eosin Y / PEGDA 575 water system.
- [ ] If Type II initiation is the lab chemistry, implement bleaching ODEs as a **separate** kinetics branch.

### Long-term
- [ ] Spatial \(\eta(x)\) on a Stage 4 field, still as a post-processor.

## Verification

```
python3 campaigns/pegda/drivers/validate_pegda_stage6.py --smoke
python3 campaigns/pegda/drivers/validate_pegda_stage6.py
python3 campaigns/pegda/drivers/run_pegda_stage6.py
python3 -m pytest campaigns/pegda/tests -q
```

OVERALL PASS (smoke and full). 29 tests passed.

## References

- Zhu et al., *J. Mech. Phys. Solids* **142**, 104041 (2020), Eqs. 21, 26, 32–38, Table A2.
- Macosko & Miller (1976) for the recursive looking-out construction.
- Prior: `journal/entries/2026-09-13-pegda-stage5.md`.
- `docs/campaign_roadmap.pdf` Stage 6.
