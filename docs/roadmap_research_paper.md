# Roadmap: from photopolymerization literature to Python

**Date:** 2026-09-13  
**Evidence:** planning decision (`proposed`) — no solver has been run in this repository  
**Active scientific object:** a **PEGDA free-radical photopolymerization solver** that maps a known intensity history (later a DMD field and a flow field) onto conversion, then onto loop-aware gelation and modulus.

This document is the **code-facing** plan for the PDFs in `Literature/`. A companion physics synthesis already exists in the Undermind workspace [DFG DMD hydrogel patches](https://app.undermind.ai/projects/c73c9d01-f906-4382-8811-3c89a1330326?path=/modelling/Unified%20photopolymerization%20modelling%20roadmap) (`Unified photopolymerization modelling roadmap`, 2026-09-13). Do not treat that note as executed code.

---

## Why a staged Python code, not one paper’s model

The five unique papers in `Literature/` solve **different pieces** of the same process. None of them is a complete “print a grayscale hydrogel in flow” code. Copying any single paper into a 3D PDE on day one would mix untested optics, stiff kinetics, oxygen transport, and network statistics — so a green run would not tell you *which* physics is wrong.

| File in `Literature/` | Cite | What it actually solves | What it should become in Python |
|-----------------------|------|-------------------------|----------------------------------|
| `Adv Funct Materials - 2024 - Dobson - …pdf` + `adfm202312607-sup-0001-suppmat.pdf` | Dobson & Bowman, *Adv. Funct. Mater.* **34**, 2312607 (2024), [DOI](https://doi.org/10.1002/adfm.202312607) | First-principles **multi-species** FRP: free-volume $D_C$, short/long radicals, oxygen, heat, 1D/3D balances (their Eq. 1–3 and SI scheme) | **High-fidelity 0D/1D reference**, not the first production solver. ~27 parameters; HDDA, not PEGDA. |
| `1-s2.0-S2352431622000645-main.pdf` | Montgomery et al., *Extreme Mech. Lett.* **53**, 101714 (2022), [DOI](https://doi.org/10.1016/j.eml.2022.101714) | **Four-species** reaction–diffusion + Beer–Lambert $A(p,C_I)$ + Gaussian DMD superposition + conversion-dependent $k_p$, $k_t$ (Eqs. 4–23) | **Production spatial RHS** after the local ODE is stable. PEGDA 250 / Irgacure 819. Finite differences. |
| `1-s2.0-S0924013619305199-main.pdf` (duplicate: `… (1).pdf`) | Meenakshisundaram et al., *J. Mater. Process. Tech.* **279**, 116546 (2020), [DOI](https://doi.org/10.1016/j.jmatprotec.2019.116546) | **Scan + mask** dose from **measured** (non-Gaussian) pixel irradiance; Jacobs-style cure profile | **Projection engine** $I(x,y,0,t)$ only. Do not use a scalar $E_c$/$D_p$ threshold as the chemistry. |
| `c9sm01485c.pdf` | Slutzky, Stone & Nunes, *Soft Matter* **15**, 9553 (2019), [DOI](https://doi.org/10.1039/c9sm01485c) | Steady **plug-flow** $U\,\mathrm{d}c/\mathrm{d}x = R$ for PI, $R^\bullet$, $M$, $\mathrm{O}_2$ (Eqs. 1a–d); critical speed and $t_\mathrm{uv}/t_L$ | **Mandatory flow benchmark** before 2D advection. PEGDA 575 jet. |
| `1-s2.0-S0022509620302763-main.pdf` | Zhu et al., *J. Mech. Phys. Solids* **142**, 104041 (2020), [DOI](https://doi.org/10.1016/j.jmps.2020.104041) | Type-II **eosin Y photobleaching** + acrylate conversion + **Macosko/Miller loops** → $P_\mathrm{gel}$, $\eta$, $G$ (Eqs. 1–34+) | **Network post-processor** on $p(\mathbf{x},t)$. Different initiator class than Irgacure 819. |

**One-sentence takeaway:** kinetics (Montgomery/Dobson) × optics (Meenakshisundaram + Beer–Lambert) × transport (Slutzky) × network (Zhu) must be **separate Python modules** with separate tests.

---

## Locked operating mode (this repository)

- **Chemistry:** PEGDA chain-growth free-radical photopolymerization (acrylate conversion $p$). First numerical rates from Montgomery; Dobson used to *stress-test* closures, not as the default RHS.
- **Illumination:** prescribed $I(t)$ first; then a DMD-generated surface field; then depth-resolved attenuation. Pattern velocity $\mathbf{v}_\mathrm{pattern}$ and resin velocity $\mathbf{u}$ are **independent**.
- **Mechanics:** Zhu’s loop-aware $\eta(p,[\mathrm{PEGDA}]_0)$ **after** conversion is trusted. No hyperelastic FEM in this ladder.
- **Primary claim of the ladder:** a well-tested Python stack can predict $p(\mathbf{x},t)$ for known light and (later) known plug flow, and then $\eta$ and $G$ for PEGDA — tagged `smoke-tested` / `validated` / `matched` per stage.

### Explicit non-goals (this phase)

1. Full Dobson 3D voxel + chain-length composite model as the everyday solver.
2. Coupled Navier–Stokes / gelation (viscosity feedback on $\mathbf{u}$).
3. Scattering / Maxwell optics / ceramic-filled resins.
4. Inverse grayscale design, shrinkage stress, poroelastic swelling, PNIPAM / thiol–ene closures.
5. Treating Jacobs working-curve dimensions as a substitute for species balances.

### Parked for later

| Item | When to reopen | Needed inputs |
|------|----------------|---------------|
| Dobson full SI species set | After reduced $k_p(p),k_t(p)$ fail equal-dose tests | HDDA FTIR or own two-intensity kinetic series |
| Thermal energy equation | After isothermal conversion is wrong vs measured $T$ | $\Delta H_p$, $h$, $c_p$ |
| Zhu eosin Y bleaching ODE (Eqs. 3–8) | If Type-II visible initiation is the lab chemistry | EY/TEOA UV–vis |
| Inverse control of DMD bitmaps | After Stage 4 grayscale widths are `matched` | Target $\eta$ maps |

---

## Conceptual stack (what the code must compute)

```text
bitmap / scan  ─►  I_⊥(x,y,t)     [Meenakshisundaram]
                         │
                         ▼
               I(x,y,z,t) = I_⊥ exp(−∫ A dz)   [Montgomery Eq. 4–6]
                         │
                         ▼
        ∂c_i/∂t + ∇·(u c_i) = ∇·(D_i ∇c_i) + R_i
        i ∈ {PI, R•, O2, M}                     [Montgomery 16–19; Slutzky 1a–d]
                         │
                         ▼
                    p = 1 − C_M / C_{M,0}
                         │
                         ▼
              P_gel, η, G(p, [PEGDA]_0)         [Zhu]
```

**Why $I\times t$ is not enough:** bimolecular termination $\propto C_R^2$ and oxygen inhibition $\propto C_O C_R$ make conversion **history-dependent**. Equal dose at two intensities must *not* collapse (Montgomery/Dobson phenomenology; literature on reciprocity failure). That is a Stage 1 test, not a later surprise.

**Why projector motion ≠ resin motion:** Meenakshisundaram’s scan speed sets how the **pattern** dwells on the vat. Slutzky’s $U$ is how **material** crosses a fixed UV window. The relative velocity $\mathbf{v}_\mathrm{pattern}-\mathbf{u}$ is the quantity a fluid element sees.

---

## Python package layout (implement in this order)

Keep the local ODE RHS identical to the spatial cell RHS so Stage 4 cannot silently rewrite chemistry.

```text
photopolymerization/
  pyproject.toml
  src/photopolymerization/
    units.py              # SI internally: mol m^-3, W m^-2, s, K
    state.py              # CI, CR, CO, CM, T
    kinetics/
      rhs.py              # Montgomery four-species ODE (Stage 1)
      closures.py         # kp(p), kt(p) = kt,D + kt,RD
      protocols.py        # I(t): constant, pulse, equal-dose pairs
    optics/
      beer_lambert.py     # dI/dz = −A I ; A from Montgomery Eq. 6
    projection/
      psf.py              # measured kernel or Gaussian fallback
      dose.py             # I_⊥ = Σ B_ij(t) K_ij ; scan frames
    transport/
      reaction_diffusion.py
      plug_flow.py        # Slutzky 1D
    network/
      zhu_macosko.py      # P_gel, η, G — post-process only
    io/
      parameters.py       # YAML → typed params
  campaigns/pegda/
    configs/
    drivers/              # run_*, validate_*
    tests/
  docs/figures/
  journal/entries/
```

**Libraries (deliberate, small):** NumPy; SciPy `solve_ivp(..., method="BDF")` for Stage 1 (radicals are stiff); SciPy `solve_bvp` or method-of-lines + BDF for 1D; later `scipy.sparse` for 2D implicit diffusion. Do not start with FEniCS/MOOSE.

**Units trap:** Montgomery’s $\beta$ is tied to *their* intensity convention. Recalibrate $\beta$ whenever $I$ is converted (mW cm$^{-2}$ vs W m$^{-2}$).

---

## Campaign stages (one quartet per stage)

Each stage: `config` · `run` · `validate` · `pytest` · figure · journal. **Stop** after tagging evidence. Do not auto-start the next stage.

### Stage 0 — Parameter book (`proposed`)

**Claim:** every rate, diffusivity, and optical coefficient is traced to a paper table or an explicit “unknown — measure” flag.

- Transcribe Montgomery Table 1 (optics) and kinetic closures; flag PEGDA MW mismatch (250 vs Zhu 575 vs Slutzky 575).
- Transcribe Slutzky Table 2 rates and $k_d$ range.
- No fitting until Stage 1 tests exist.

**Gates:** YAML schema; no silent `NaN` defaults; units in comments.

### Stage 1 — Local stiff ODE (`proposed` → then `smoke-tested`)

**Claim:** for prescribed $I(t)$, the four-species Montgomery RHS plus $k_p(p)$, $k_t(p)$ integrates without negative concentrations and shows an oxygen induction period.

Equations (Montgomery 16–19 with $D_i=0$):

\[
\dot C_I = -\beta I C_I,\quad
\dot C_R = m\beta I C_I - 2k_t(p)C_R^2 - k_O C_O C_R,
\]
\[
\dot C_O = -k_O C_O C_R,\quad
\dot C_M = -k_p(p) C_R C_M,\quad
p=1-C_M/C_{M,0}.
\]

**Gates (`validate_pegda_stage1.py`):**

| Gate | Pass criterion |
|------|----------------|
| Dark | $I=0$ ⇒ $p$ unchanged at solver tolerance |
| Positivity | $C_i\ge -10^{-12}$ mol m$^{-3}$ |
| Bounds | $0\le p\le 1$ |
| PI decay | $C_I$ monotone when $I>0$ |
| Reciprocity | two equal-dose protocols; $\lvert p_1-p_2\rvert$ at $t_\mathrm{end}$ **not** both ~0 (history dependence visible) |
| Mass | $\int R_p\,\mathrm{d}t$ matches $\Delta C_M$ |

**Code:** `scipy.integrate.solve_ivp` BDF/Radau, `rtol=1e-7`. Optional thermal $\dot T$ **off**.

### Stage 2 — Analytic / reduced limits (`proposed`)

**Claim:** the ODE recovers known limits, not just “looks like a conversion curve.”

- Slutzky $C_O/C_{O,0}=(C_M/C_{M,0})^{k_O/k_p}$ when $k_p,k_t$ constant and no diffusion (their Eq. 3) — use a **constant-rate fork** of the RHS for this MMS-style check.
- Quasi-steady radical $C_R\sim\sqrt{R_i/2k_t}$ when $C_O=0$ and $k_t$ constant (Zhu Eq. 15 analogue).
- Dobson used as **qualitative** reference: $R_p$ vs $I$ exponent should move from ~0.5 (bimolecular) toward ~1 when oxygen dominates — **do not** require HDDA numbers.

**Gates:** `np.testing.assert_allclose` on the constant-rate analytic relations.

### Stage 3 — Projection engine (`proposed`)

**Claim:** a bitmap + PSF + scan produces $I_\perp(x,y,t)$ and $E_\perp=\int I_\perp\mathrm{d}t$ that conserve energy (sum of pixel powers) and, for a Gaussian PSF, match Montgomery Eq. 7.

**Gates:** single-pixel integral ≈ $I_\mathrm{pixel}\times$ area; ten overlapping Gaussians reproduce the paper’s qualitative bleed (their Fig. 3); Meenakshisundaram scan: exposure independent of $v$ only if dwell $\propto 1/v$ (document if not).

**Non-goal:** cure depth from Jacobs $C_d=D_p\ln(E/E_c)$ except as a **diagnostic overlay**.

### Stage 4 — 1D then 2D reaction–diffusion (`proposed`)

**Claim:** method-of-lines + same RHS + $D_i(p)$ (Montgomery harmonic liquid–solid) + Beer–Lambert reproduces (i) oxygen dead-zone at an $\mathrm{O}_2$ Dirichlet wall, (ii) lateral grayscale blur wider than the optical PSF.

**Gates:** Neumann dark box conserves $\int C_M$; manufactured diffusion-only; CFL/Fourier documented; grayscale step $G_0$|$G_{80}$ transition width $> \sigma_\mathrm{pixel}$ when $D_R$ or $D_M$ on.

**Solver:** 1D first (depth $z$ only), then $x$–$z$ grayscale strip. Not 3D.

### Stage 5 — Plug-flow benchmark (`proposed`)

**Claim:** Slutzky Eqs. 1a–d as $\mathrm{d}y/\mathrm{d}x = R(y)/U$ recover $U_c$ scaling with $I$ and $C_{\mathrm{PI},0}/C_{\mathrm{O}_2,0}$ and the $t_\mathrm{uv}/t_L$ collapse for fiber length **qualitatively**. Quantitative `matched` only if we use their Table 2 numbers and gel criterion $M/M_0=0.98$.

**Gates:** $U\to\infty$ ⇒ no conversion; $U\to 0$ with long $t_\mathrm{uv}$ ⇒ gel; oxygen diffusion optional (they argue $\sqrt{Dt}\ll w_1$).

### Stage 6 — Network post-process (`proposed`)

**Claim:** given $p$ and $[\mathrm{PEGDA}]_0$, Zhu $P_\mathrm{gel}=(1-\theta)/(2\theta)$ with $\theta=[\mathrm{PEGDA}]_0/([\mathrm{PEGDA}]_0+2C_0)$ predicts no gel below $C_0$, and $G=\eta N_A [\mathrm{PEGDA}]_0 k_B T$ is below affine ideal-network $G$.

**Gates:** $\theta\to 1$ ⇒ $P_\mathrm{gel}\to 0$; $[\mathrm{PEGDA}]_0\le C_0$ ⇒ never gel; $\eta\le 1$.

**Do not** mix eosin Y bleaching kinetics into the Irgacure RHS without a chemistry flag.

### Stage 7+ (parked)

Moving pattern $\mathbf{v}_\mathrm{pattern}\ne 0$; Eulerian–Lagrangian markers; Dobson free-volume $D_C$; thermal; inverse DMD.

---

## Dimensionless tests (use these to choose the next PDE)

| Group | Meaning | If large / O(1) |
|-------|---------|-----------------|
| $\beta_\mathrm{opt}=A H$ | optical depth | resolve $z$ (Stage 4) |
| $\mathrm{Da}_O\sim \phi\varepsilon C_I I H^2/(D_{\mathrm{O}_2} C_{\mathrm{O}_2})$ | oxygen reaction vs diffusion | wall dead-zone |
| $\Pi_\mathrm{flow}=L k_d C_{I,0}/(U C_{\mathrm{O}_2,0})$ | flow oxygen depletion | Stage 5 gel/no-gel |
| $\tau_\mathrm{exp}=t_\mathrm{uv}/t_L$ | pulse vs residence | Slutzky fiber length |
| $\Lambda_v=\lvert\mathbf{v}_\mathrm{pattern}-\mathbf{u}\rvert/U_*$ | relative pattern motion | do not collapse scan and flow |
| $\mathrm{Pe}_i=U H/D_i$ | advection vs diffusion | which species need $\nabla\cdot D\nabla c$ |
| $\chi=\Delta t_\mathrm{frame}/t_\mathrm{chem}$ | DMD frame vs chemistry | keep discrete frames |

---

## Evidence rules

- **Engineering success** (finite, no NaN, pytest green) is not **scientific validity**.
- `validated` requires the named gates above, not a pretty conversion plot.
- `matched` is reserved for numbers taken from a paper figure (e.g. Montgomery FTIR $p(t)$, Slutzky $U_c$) with the **same formulation**.
- Montgomery kinetic numbers **do not transfer** to Zhu’s eosin Y / PEGDA 575 water gel without refit.
- Duplicate JMPT PDF: cite one file only.

---

## Immediate / medium / long-term

### Immediate

- [ ] Scaffold `src/photopolymerization` + `pyproject.toml` + Stage 1 YAML.
- [ ] Implement Stage 1 ODE + `validate` gates + pytest.
- [ ] Deduplicate `Literature/` (keep one JMPT PDF).
- [ ] Build a parameter table in YAML with paper provenance.

### Medium

- [ ] Stage 3 Gaussian PSF (then measured PSF if a camera map exists).
- [ ] Stage 4 1D Beer–Lambert + oxygen wall.
- [ ] Stage 5 Slutzky integrator vs Table 2.

### Long-term

- [ ] Stage 6 Zhu network on spatial $p$.
- [ ] Lab FTIR / radiometry to retag stages `matched`.
- [ ] Reopen Dobson SI only if reduced $k_t(p)$ fails two-intensity data.

---

## Verification of *this* planning document

- Corpus: five unique PDFs + Dobson SI in `Literature/` (listed above).
- Equations cited by paper number, not invented.
- No Python solver executed in this repository as of 2026-09-13 (`proposed`).

## References

- Local PDFs as in the table.  
- Undermind: [Unified photopolymerization modelling roadmap](https://app.undermind.ai/projects/c73c9d01-f906-4382-8811-3c89a1330326?path=/modelling/Unified%20photopolymerization%20modelling%20roadmap); [Stage 1 PEGDA solver implementation](https://app.undermind.ai/projects/c73c9d01-f906-4382-8811-3c89a1330326?path=/modelling/Stage%201%20PEGDA%20solver%20implementation).  
- Journal: `journal/entries/2026-09-13-planning-python-photopolymerization.md`.  
- Index: `docs/CAMPAIGN_INDEX.md`.
