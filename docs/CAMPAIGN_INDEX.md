# Campaign index

**Updated:** 2026-09-13  
**Active scientific object:** PEGDA free-radical photopolymerization in Python — local conversion first, then light, space, flow, then Zhu network.  
**Parked:** full Dobson 3D species set; NS–gel coupling; scattering; inverse DMD; PNIPAM/thiol–ene; eosin Y kinetics mixed into Irgacure RHS.  
**Thesis-style report:** `docs/thesis_photopolymerization.tex` / `.pdf` (Stages 0–6 internally validated; none matched).

## Active campaign

| Stage | Claim (one line) | Evidence | Entry points |
|-------|------------------|----------|--------------|
| 0 Parameter book | Rates/optics traced to Literature PDFs | validated | `campaigns/pegda/configs/parameter_book.montgomery_2022.json` · `validate_pegda_stage0.py` · journal 2026-09-13 |
| 1 Local ODE | Montgomery four-species + $k_p(p),k_t(p)$ vs $I(t)$ | validated (internal gates; not matched) | `run_pegda_stage1.py` · `validate_pegda_stage1.py` · `docs/figures/pegda_stage1_overview.pdf` |
| 2 Analytic limits | Constant-rate Slutzky/Zhu QSS identities | validated (analytic; not matched) | `validate_pegda_stage2.py` · `docs/figures/pegda_stage2_overview.pdf` |
| 3 Projection | Bitmap + PSF + scan → $I_\perp,E_\perp$ | validated (Gaussian fallback; PSF not measured) | `validate_pegda_stage3.py` · `docs/figures/pegda_stage3_overview.pdf` |
| 4 RD + Beer–Lambert | 1D species + front, O₂ wall, grayscale vs σ | validated (internal; D increment unresolved; not matched) | `validate_pegda_stage4.py` · `docs/figures/pegda_stage4_overview.pdf` |
| 5 Plug flow | Slutzky $U_c$, $t_\mathrm{uv}/t_L$, 2% gel predicate | validated (Table 2 ODE; not matched to $U_c$ data) | `validate_pegda_stage5.py` · `docs/figures/pegda_stage5_overview.pdf` |
| 6 Network | Zhu $P_\mathrm{gel}$, $\eta$, $G$ post-process | validated (algebraic; not matched to $G(p)$ data) | `validate_pegda_stage6.py` · `docs/figures/pegda_stage6_overview.pdf` |

## Archived / contrast campaigns

| Campaign | Status | What it answered | Do not reuse as evidence for |
|----------|--------|------------------|------------------------------|
| none | — | empty repo except Literature | — |

## Open gates

- [x] Stage 0 parameter provenance (C_O0 still `still_required`)
- [x] Stage 1 BDF integrator + positivity / reciprocity / induction tests
- [x] Stage 2 analytic identities
- [x] Stage 6 Zhu network post-processor (eosin Y ODEs off; not matched)
- [ ] FTIR / oxygen probe to retag Stage 1 `matched`
- [ ] Measured DMD PSF (Meenakshisundaram) to retag Stage 3 `matched`
- [ ] Stage 5 experimental $U_c$ map to retag `matched`
- [ ] Zhu $G(p)$ rheology to retag Stage 6 `matched`

## Invalid / do-not-cite

| Artifact / claim | Why | Replacement |
|------------------|-----|-------------|
| Jacobs $C_d=D_p\ln(E/E_c)$ as the polymer state | Dose threshold is not species conversion | Keep as diagnostic overlay only |
| Montgomery $k_p,k_t$ as universal PEGDA constants | Fit to one thin-film PEGDA 250 recipe | Recalibrate; tag `matched` only vs that paper’s FTIR |
| Duplicate JMPT PDF | Same file twice | Cite `1-s2.0-S0924013619305199-main.pdf` only |

## Documentation entry points

| Doc | Role |
|-----|------|
| `docs/roadmap_research_paper.md` | Short code-facing stage list |
| `docs/campaign_roadmap.pdf` | Full staged roadmap: actions, gates, literature agreement |
| `docs/background_photopolymerization.pdf` | Chemistry/process ontology (masters primer) |
| `journal/entries/2026-09-13-planning-python-photopolymerization.md` | Scope lock |
| Undermind modelling notes | Physics synthesis (not executed code) |

## Process note

Literature-only solver bootstrap. Stages 0–6 of the planned ladder are internally `validated`. Next: retag `matched` on same-formulation data, or reopen a parked item.
