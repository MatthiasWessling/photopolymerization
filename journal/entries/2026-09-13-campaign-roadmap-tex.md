# Lab note: campaign roadmap TeX/PDF

**Date:** 2026-09-13  
**Evidence:** planning document (`proposed`); no solver run

## Motivation

The markdown roadmap listed stages but did not freeze, in a citable lab-note form, what is *done* at each stage, what success means, which gates block the next stage, and how literature agreement is judged (phenomenology vs analytic identity vs same-formulation match).

## What was done

- Wrote `docs/campaign_roadmap.tex` (companion to `docs/background_photopolymerization.tex`, shared `.bib`).
- Added `docs/figures/fig_roadmap_stage_ladder.pdf`.
- Compiled with `latexmk -pdf` to `docs/campaign_roadmap.pdf` (17 pages, citations resolved).
- Linked from `docs/CAMPAIGN_INDEX.md`.

## What improved

Stage cards 0–6 now have: scientific question, actions, success factors, named gates, reasonability/literature protocol, tag-and-stop. Cross-cutting reflection (five journal questions) and three kinds of literature agreement are explicit. Parked work is listed with reopen conditions.

## Suggestions to proceed

### Immediate

- [ ] Stage 0 parameter book with provenance columns.
- [ ] Stage 1 local constitutive response and the Stage 1 gate table.

### Medium

- [ ] Stage 2 analytic identities before any spatial PDE.
- [ ] Stage 3 projection with chemistry *off*.

### Long-term

- [ ] Retag stages `matched` only on same-formulation figures (Montgomery FTIR, Slutzky $U_c$, Zhu $G$).

## Verification

```
latexmk -pdf -interaction=nonstopmode docs/campaign_roadmap.tex
```

Final log: 17 pages; no remaining undefined citations.

## References

- `docs/campaign_roadmap.pdf`
- `docs/background_photopolymerization.pdf`
- `docs/roadmap_research_paper.md`
- Journal `2026-09-13-planning-python-photopolymerization.md`
