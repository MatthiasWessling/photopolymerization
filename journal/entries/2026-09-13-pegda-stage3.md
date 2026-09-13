# Stage 3: Gaussian projection engine

**Date:** 2026-09-13  
**Evidence:** `validated` (Montgomery Eq. 7 fallback). Not `matched` to a camera PSF.

## Motivation

Single claim: a bitmap plus Montgomery’s Gaussian pixel kernel produces a surface irradiance \(I_\perp\) that conserves the analytic pixel integral, bleeds across commanded edges, and yields scan dose \(\propto 1/v\). This stage has **no chemistry** and is not a polymer state.

## What was done

- `src/photopolymerization/projection/psf.py` (Eq. 7 superposition, 1-D/2-D integrals, linear G0–G100 map).
- `config.pegda_stage3.json`, `run_pegda_stage3.py`, `validate_pegda_stage3.py`, `tests/test_pegda_stage3.py`.
- Figure: `docs/figures/pegda_stage3_overview.pdf` (Fig. 3-style 1-D panels).

## What improved

| Gate | Result |
|------|--------|
| 1-D \(\int I\,dx\) vs \(I_p\,2\sigma\sqrt{\pi}\) | rel. error \(0\) at 1201 points |
| 2-D \(\iint I\,dA\) vs \(I_p\,4\pi\sigma^2\) | rel. error \(\sim 10^{-16}\) |
| Dark bitmap | \(I=0\) |
| Ten-pixel bleed | \(I(\mathrm{edge})=25.4\), \(I(\mathrm{interior})=50.8\,\mathrm{W\,m^{-2}}\) |
| \(\max(E)\,v\) at two speeds | identical to \(2\times 10^{-2}\) |

**Reasonability:** \(\int I\,dA\) is **not** \(I_{\mathrm{pixel}}\times(50\,\mu\mathrm{m})^2\). Peak \(I_{\mathrm{pixel}}=41.08\,\mathrm{W\,m^{-2}}\) is an irradiance, not a power. Grayscale is a linear placeholder; Montgomery SI RGB calibration was not transcribed. Meenakshisundaram measured non-Gaussian kernels — still required for `matched`.

## Suggestions to proceed

### Immediate
- [ ] Stage 4 1-D RD + Beer–Lambert only.
- [ ] Do not feed Jacobs \(C_d\) as conversion.

### Medium
- [ ] Replace Gaussian with a measured single-pixel map.

### Long-term
- [ ] Discrete DMD frame cycling when \(\chi=\Delta t_{\mathrm{frame}}/t_{\mathrm{chem}}\) is not small.

## Verification

```
python3 campaigns/pegda/drivers/validate_pegda_stage3.py --smoke
python3 campaigns/pegda/drivers/validate_pegda_stage3.py
python3 campaigns/pegda/drivers/run_pegda_stage3.py
python3 -m pytest campaigns/pegda/tests -q
```

OVERALL PASS (smoke and full). 14 tests passed.

## References

- Montgomery et al. (2022) Eq. 7, Table 1, Fig. 3.
- Meenakshisundaram et al. (2020) for measured PSF (not implemented).
- `docs/campaign_roadmap.pdf` Stage 3.
