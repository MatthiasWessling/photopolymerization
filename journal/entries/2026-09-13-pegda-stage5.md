# Stage 5: Slutzky plug flow

**Date:** 2026-09-13  
**Evidence:** `validated` (internal ODE gates on Table 2). Not `matched` to their experimental \(U_c\) map.

## Motivation

Single claim: a uniform jet of speed \(U\) through a fixed UV window of length \(L\) obeys Slutzky Eqs. 1a–d, so fast flow does not gel, slow/high-\(k_d\) flow crosses their operational \(M/M_0=0.98\) threshold, Eq. 3 holds along \(x\), and \(U_c\) increases with \(k_d\) and \([\mathrm{PI}]_0\). This is **not** Montgomery PEGDA-250 kinetics and **not** Zhu \(P_{\mathrm{gel}}\).

## What was done

- Separate book: `campaigns/pegda/configs/parameter_book.slutzky_2019.json` (PEGDA 575 / Darocur jet).
- `src/photopolymerization/flow/plug.py`: \(U\,\mathrm{d}y/\mathrm{d}x=R\), \(k_d=0\) for \(x>L\), termination \(k_t R^2\) as written, one radical per photolysis.
- Fig. 2 / Fig. 3 parameter cases from the paper; gel predicate documented as 2% double-bond conversion.
- `config.pegda_stage5.json`, `run_pegda_stage5.py`, `validate_pegda_stage5.py`, `tests/test_pegda_stage5.py`.
- Figure: `docs/figures/pegda_stage5_overview.pdf`.
- Wall oxygen diffusion off (their \(\sqrt{Dt}\ll w_1\) argument).

## What improved

| Gate | Result |
|------|--------|
| Fig. 2 high O₂, \(U=0.003\,\mathrm{m\,s^{-1}}\) | \(\min M/M_0=1.00\), no gel |
| Fig. 3 low O₂ | \(\min M/M_0=0.863\) (paper \(\simeq 0.86\)) |
| Eq. 3 along \(x\) | residual \(7.9\times 10^{-7}\) |
| \(U=1\,\mathrm{m\,s^{-1}}\) | \(\min M/M_0=1\) |
| \(U_c(k_d)\) | \(0.0018\to 0.0072\,\mathrm{m\,s^{-1}}\) when \(k_d\) \(\times 4\) |
| \(U_c([\mathrm{PI}]_0)\) | \(0.0023\to 0.0045\,\mathrm{m\,s^{-1}}\) when \(\mathrm{Pi}_0\) \(\times 2\) |
| Pulse length \(U t_{\mathrm{uv}}\) | increases with \(t_{\mathrm{uv}}\) |

**Reasonability:** numerical Fig. 3 conversion matches the paper’s quoted \(M\simeq 0.86\). Scaling \(U_c\sim L k_d \mathrm{Pi}_0/y_0\) is recovered. \(M_0\) is calculated from 54 vol% PEGDA-575 with an assumed density (does not affect \(M/M_0\)). Pulse lengths are **kinematic** (\(U t_{\mathrm{uv}}\)), not a measured fiber-length campaign.

## Suggestions to proceed

### Immediate
- [ ] Stage 6 Zhu network post-processor only. Do not retune \(k_p\) from \(G\).
- [ ] Do not plot \(M/M_0=0.98\) as \(P_{\mathrm{gel}}\).

### Medium
- [ ] Overlay their experimental \(U_c(I,[\mathrm{PI}]_0)\) for a `matched` retag.
- [ ] Lagrangian check \(t=x/U\) documented as equivalent (already the ODE).

### Long-term
- [ ] Prescribed \(\mathbf{u}(\mathbf{x})\) beyond plug flow after Stage 6.

## Verification

```
python3 campaigns/pegda/drivers/validate_pegda_stage5.py --smoke
python3 campaigns/pegda/drivers/validate_pegda_stage5.py
python3 campaigns/pegda/drivers/run_pegda_stage5.py
python3 -m pytest campaigns/pegda/tests -q
```

OVERALL PASS (smoke and full). 24 tests passed.

## References

- Slutzky, Stone & Nunes, *Soft Matter* **15**, 9553 (2019), Eqs. 1a–d, 3, Table 2, Figs. 2–3.
- Dendukuri et al. for the 2% gel convention they cite.
- Prior: `journal/entries/2026-09-13-pegda-stage4.md`.
- `docs/campaign_roadmap.pdf` Stage 5.
