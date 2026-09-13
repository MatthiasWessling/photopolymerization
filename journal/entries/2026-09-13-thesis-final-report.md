# Thesis-style final report (TeX + PDF)

**Date:** 2026-09-13  
**Evidence:** laboratory documentation (`validated` campaign stages 0–6 as cited). The report itself is not a `matched` materials paper.

## Motivation

Stages 0–6 are internally gated. The group still needed a single, citable, master’s-thesis-length synthesis: chemistry ontology, literature partition, research gap, software, numerical methods, stage results, and an honest matching chapter. Companion notes (`background_photopolymerization.tex`, `campaign_roadmap.tex`) stay as primers; they are not replaced.

## What was done

- Driver: `docs/thesis_photopolymerization.tex` (report class, BasicTeX `latexmk`).
- Chapters: `docs/thesis_ch_intro.tex`, `thesis_ch_background.tex` (body extracted from the background note), `thesis_ch_literature.tex`, `thesis_ch_methods.tex`, `thesis_ch_results.tex`, `thesis_ch_matching.tex` (includes conclusions), `thesis_ch_appendix.tex`.
- Bibliography: `docs/background_photopolymerization.bib`.
- PDF: `docs/thesis_photopolymerization.pdf`.
- Campaign index points at the report.

No new kinetic runs. No parked physics started. Publisher PDFs were not copied into the report.

## What improved

A reader can follow species → rates → light → flow → network, then see which Python module owns each claim and which gate passed, without upgrading `validated` to `matched`. Forbidden equivalences (dose ≢ p ≢ gel ≢ η ≢ G; Type I ≢ Type II; recipe mixing) are stated in the background chapter and enforced in the results chapter.

## Suggestions to proceed

### Immediate
- [ ] Treat FTIR overlay (Stage 1 `--match`) as the next *scientific* step, not further TeX expansion.
- [ ] Do not cite the thesis PDF as experimental PEGDA properties.

### Medium
- [ ] When a stage is retagged `matched`, add a short erratum paragraph to the thesis rather than silently rewriting history.

### Long-term
- [ ] If the document is submitted as an actual master’s thesis, add an examination cover, sworn statement, and supervisor names; those are not in this lab version.

## Verification

```
PATH="/Library/TeX/texbin:$PATH" latexmk -pdf -interaction=nonstopmode -cd docs/thesis_photopolymerization.tex
```

`latexmk` wrote `docs/thesis_photopolymerization.pdf` (**47 pages**).

## References

- Stage journals `journal/entries/2026-09-13-pegda-stage0.md` … `stage6.md`.
- `docs/background_photopolymerization.pdf`, `docs/campaign_roadmap.pdf`.
- Montgomery 2022; Slutzky 2019; Zhu 2020; Dobson 2024; Meenakshisundaram 2020 (unique JMPT file only).
