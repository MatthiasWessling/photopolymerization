# Campaign index

**Updated:** 2026-09-13  
**Active scientific object:** PEGDA free-radical photopolymerization in Python — local conversion first, then light, space, flow, then Zhu network.  
**Parked:** full Dobson 3D species set; NS–gel coupling; scattering; inverse DMD; PNIPAM/thiol–ene; eosin Y kinetics mixed into Irgacure RHS.

## Active campaign

| Stage | Claim (one line) | Evidence | Entry points |
|-------|------------------|----------|--------------|
| 0 Parameter book | Rates/optics traced to Literature PDFs | proposed | roadmap |
| 1 Local ODE | Montgomery four-species + $k_p(p),k_t(p)$ vs $I(t)$ | proposed | roadmap §Stage 1 |
| 2 Analytic limits | Constant-rate Slutzky/Zhu QSS identities | proposed | roadmap §Stage 2 |
| 3 Projection | Bitmap + PSF + scan → $I_\perp,E_\perp$ | proposed | Meenakshisundaram / Montgomery Eq. 7 |
| 4 RD + Beer–Lambert | 1D/2D species + grayscale blur | proposed | Montgomery Eqs. 4–23 |
| 5 Plug flow | Slutzky $U_c$, $t_\mathrm{uv}/t_L$ | proposed | Soft Matter Eqs. 1a–d |
| 6 Network | Zhu $P_\mathrm{gel}$, $\eta$, $G$ post-process | proposed | JMPS recursive model |

## Archived / contrast campaigns

| Campaign | Status | What it answered | Do not reuse as evidence for |
|----------|--------|------------------|------------------------------|
| none | — | empty repo except Literature | — |

## Open gates

- [ ] Stage 1 BDF integrator + positivity / reciprocity tests
- [ ] YAML parameter provenance (do not mix PEGDA 250 vs 575 vs HDDA)
- [ ] Choose lab initiator class (Irgacure 819 vs eosin Y) before mixing Zhu bleaching with Montgomery RHS

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

Literature-only repository. Next process skill: **campaign-stage** for Stage 1 ODE only.
