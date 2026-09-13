# Lab note: masters background on photopolymerization ontology

**Date:** 2026-09-13  
**Evidence:** documentation (`proposed` ontology; not a validation campaign)

## Motivation

Before coding Stage 1, freeze a chemistry-and-process background that a new masters student can read: what objects exist, which equivalences are forbidden, and how the five campaign papers partition those objects. Programming jargon was explicitly out of scope.

## What was done

- Wrote `docs/background_photopolymerization.tex` and `.bib`.
- Generated teaching cartoons (not data) with `docs/figures/make_background_schematics.py`.
- Compiled with BasicTeX `latexmk -pdf` to `docs/background_photopolymerization.pdf` (15 pages).
- Linked the PDF from `docs/CAMPAIGN_INDEX.md`.

## What improved

A citable, equation-numbered ontology: Type I vs II, conversion vs gelation vs modulus, two velocities, oxygen as inventory vs wall flux, failed reciprocity. Aligns the Python ladder with the Undermind Stage 1 notes without describing software.

## Suggestions to proceed

### Immediate

- [ ] Implement Stage 1 local kinetics against this ontology (BDF ODE, equal-dose test allowed to fail).
- [ ] Do not mix Irgacure and eosin Y rate laws in one constitutive object.

### Medium

- [ ] Add Den08-style oxygen-wall thought experiment as a Stage 4 variant in the roadmap.
- [ ] Replace schematic Fig. 2 with a replot of a named Wydra/Montgomery figure if copyright-clean.

### Long-term

- [ ] Formulation-specific parameter table (PEGDA 250 vs 575 vs HDDA) as a one-page addendum.

## Verification

```
latexmk -pdf -interaction=nonstopmode docs/background_photopolymerization.tex
```

Final log: PDF written, 15 pages; no remaining `undefined citations` in `.log`.

## References

- `docs/background_photopolymerization.pdf`
- `docs/roadmap_research_paper.md`
- Campaign papers as in the `.bib` file (Dobson 2024; Montgomery 2022; Meenakshisundaram 2020; Slutzky 2019; Zhu 2020; Wydra 2014; Dendukuri 2008).
