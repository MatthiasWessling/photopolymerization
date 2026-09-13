# Unified photopolymerization modelling roadmap

##### [**Undermind**](https://undermind.ai)

---


## Table of Contents

- [Unified photopolymerization modelling roadmap](#unified-photopolymerization-modelling-roadmap)
  - [Aim](#aim)
  - [Two independent motions must be separated](#two-independent-motions-must-be-separated)
  - [Revised model architecture](#revised-model-architecture)
    - [1. Projection-kinematics engine](#projection-kinematics-engine)
    - [2. High-fidelity chemistry reference](#high-fidelity-chemistry-reference)
    - [3. Reduced Eulerian reaction–diffusion–advection model](#reduced-eulerian-reactiondiffusionadvection-model)
    - [4. Dynamic optical transport](#dynamic-optical-transport)
    - [5. Eulerian–Lagrangian flow coupling](#eulerianlagrangian-flow-coupling)
    - [6. Mechanistic crosslink-density closure](#mechanistic-crosslink-density-closure)
  - [What the full texts establish](#what-the-full-texts-establish)
  - [Staged modelling roadmap](#staged-modelling-roadmap)
  - [Dimensionless groups and model-selection tests](#dimensionless-groups-and-model-selection-tests)
  - [Validation ladder](#validation-ladder)
  - [Recommended first implementation](#recommended-first-implementation)
  - [Remaining gap](#remaining-gap)
  - [References](#references)

# Unified photopolymerization modelling roadmap

## Aim

The target is a mechanistic model that predicts local conversion and crosslink density when illumination varies in space and time while the resin is stationary or flowing. A single dose field is insufficient: radical termination, oxygen transport, photoinitiator bleaching, diffusion, and viscosity evolution make the result depend on the full intensity history \[Wyd14, Lin19, Mon22\].

The current literature supports four complementary anchors:

- \[Dob24\] provides a high-fidelity reaction–diffusion–thermal reference model for multifunctional acrylates.
- \[Mon22\] provides the clearest spatial grayscale-DLP implementation, including Gaussian pixel overlap and conversion-dependent diffusivity.
- \[Slu19\] provides an experimentally validated flow-through, residence-time formulation with a critical curing velocity.
- \[Zhu20\] provides a mechanistic conversion-to-gelation-to-crosslink-density closure for PEGDA hydrogels.
- \[Mee20\] adds the missing projection-kinematics layer: measured DMD point-spread functions, discrete frame cycling, scan synchronization, and dose accumulation.

No paper combines all of these with a moving DMD, flowing resin, and a full hydrogel network model. The roadmap should therefore use a projection engine, a high-fidelity reference chemistry, a reduced production solver, and staged validation rather than one monolithic model.

## Two independent motions must be separated

The central conceptual refinement from \[Mee20\] is that **projector motion and resin motion are different velocities**.

- The projected pattern may move because the optical head or mask moves.
- The projected pattern may be deliberately scrolled on the DMD to cancel that motion.
- The resin may independently move through the reactor with velocity $`\mathbf{u}(\mathbf{x},t)`$.

Let $`\mathbf{v}_{\mathrm{stage}}`$ be the projected stage velocity and $`\mathbf{v}_{\mathrm{scroll}}`$ the DMD pattern-scroll velocity. The pattern velocity in resin coordinates is approximately

``` math
\mathbf{v}_{\mathrm{pattern}}
=\mathbf{v}_{\mathrm{stage}}+Q\mathbf{v}_{\mathrm{scroll}},
```

where $`Q`$ is optical magnification. Synchronized scanning imposes $`\mathbf{v}_{\mathrm{pattern}}\approx0`$, producing a stationary projected feature in the vat frame even though the projector is moving. A moving-mask or moving-source process without compensation has $`\mathbf{v}_{\mathrm{pattern}}\ne0`$.

The material velocity $`\mathbf{u}`$ is independent of both. The final model must therefore distinguish

``` math
\mathbf{v}_{\mathrm{pattern}},
\qquad
\mathbf{u},
\qquad
\mathbf{v}_{\mathrm{pattern}}-\mathbf{u}.
```

The last quantity is the relative motion between the illumination pattern and a material element.

## Revised model architecture

### 1. Projection-kinematics engine

\[Mee20\] should be used as an optical boundary-condition generator, not as the polymerization model. Its three useful components are:

1.  **Measured pixel point-spread function.** The irradiance from one micromirror is measured experimentally and used as a spatial kernel $`K_{ij}(x,y)`$. This is preferable to assuming a Gaussian PSF: the measured profile can have a flat top and substantially more pixel bleeding than a Gaussian approximation.
2.  **Frame cycling and scan synchronization.** A bitmap is updated column by column at a frame rate set by projected pixel pitch and scan speed:

``` math
T_{\mathrm{frame}}=\frac{p_{\mathrm{pix}}}{v_{\mathrm{resin}}},
\qquad
F_{\mathrm{frame}}=\frac{v_{\mathrm{resin}}}{p_{\mathrm{pix}}},
```

where $`p_{\mathrm{pix}}=Q(2w)`$ is the projected pixel pitch and $`v_{\mathrm{resin}}=Qv_{\mathrm{proj}}`$. 3. **Dose accumulation.** The surface irradiance is the superposition

``` math
I_\perp(x,y,t)=\sum_{i,j}B_{ij}(t)K_{ij}(x,y),
```

and the local incident exposure is

``` math
E_\perp(x,y)=\int I_\perp(x,y,t)\,dt.
```

The discrete sums in \[Mee20\] can be replaced by a continuous convolution when the frame rate is high enough. The discrete form should be retained when frame cycling, pixel dwell time, or aliasing may affect feature fidelity.

This projection engine should be validated separately against surface irradiance maps and dose uniformity before it is coupled to reaction kinetics.

### 2. High-fidelity chemistry reference

\[Dob24\] tracks photoinitiator, primary radicals, short and network-bound macroradicals, unreacted and partially reacted multifunctional monomers, oxygen, and terminated species. It includes species diffusion, free-volume-dependent propagation and termination, chain-length-dependent radical mobility, optical attenuation, and reaction heat.

This is the best reference model for testing reduced kinetics. It should not be the first full flow solver: the model has many coupled species, stiff radical time scales, and roughly 27 fitted parameters. The practical strategy is to use it in 0D or 1D to generate and validate reduced closures for $`k_p(\alpha,T)`$, $`k_t(\alpha,T)`$, radical efficiency, and oxygen inhibition.

### 3. Reduced Eulerian reaction–diffusion–advection model

For the first production model, retain the state variables

``` math
 c_{\mathrm{PI}},\quad c_R,\quad c_{\mathrm{O_2}},\quad c_M,\quad \alpha,\quad T,
```

with the conservative balance

``` math
\frac{\partial c_i}{\partial t}+\nabla\cdot(\mathbf{u}c_i)
=\nabla\cdot(D_i(\alpha,T)\nabla c_i)+R_i.
```

A reduced chain-growth chemistry is

``` math
\frac{D c_{\mathrm{PI}}}{D t}=-k_b I c_{\mathrm{PI}},
```

``` math
\frac{D c_R}{D t}=2f k_b I c_{\mathrm{PI}}
-2k_t(\alpha,T)c_R^2-k_Oc_Rc_{\mathrm{O_2}},
```

``` math
\frac{D c_M}{D t}=-k_p(\alpha,T)c_Rc_M,
```

``` math
\frac{\partial c_{\mathrm{O_2}}}{\partial t}
+\nabla\cdot(\mathbf{u}c_{\mathrm{O_2}})
=\nabla\cdot(D_{\mathrm{O_2}}(\alpha,T)\nabla c_{\mathrm{O_2}})-k_Oc_Rc_{\mathrm{O_2}}.
```

Conversion is

``` math
\alpha=1-\frac{c_M}{c_{M,0}}.
```

The diffusion coefficients should evolve strongly with conversion. \[Mon22\] uses a liquid-to-solid harmonic transition spanning roughly four orders of magnitude, while \[Dob24\] derives mobility from free volume. The first implementation can use the calibrated \[Mon22\] form and later replace it with the \[Dob24\] free-volume closure.

If thermal acceleration matters, add

``` math
\rho c_p\frac{D T}{D t}
=\nabla\cdot(k_T\nabla T)+\Delta H_pR_p+q_{\mathrm{abs}}.
```

The thermal equation is part of \[Dob24\] but should initially be activated only after isothermal predictions are validated.

### 4. Dynamic optical transport

The projection engine supplies the surface field. The chemistry-dependent optical solver then propagates it through the resin:

``` math
\frac{\partial I}{\partial z}+A(\mathbf{x},t)I=0,
```

with

``` math
A=\alpha_Ic_{\mathrm{PI}}+\alpha_{\mathrm{poly}}\alpha
+\alpha_{\mathrm{mon}}(1-\alpha)+A_{\mathrm{abs}}.
```

This formulation follows \[Mon22\] and explicitly allows initiator depletion, polymer–monomer optical differences, and fixed absorbers. The incident field is therefore

``` math
I(x,y,0,t)=\mathcal{P}\left[B_{ij}(t),K_{ij},\mathbf{v}_{\mathrm{pattern}}\right],
```

where $`\mathcal{P}`$ is the measured projection operator from \[Mee20\].

This separates three optical approximations:

- **Level 1:** empirical measured PSF and linear superposition \[Mee20\].
- **Level 2:** dynamic Beer–Lambert attenuation coupled to photochemistry \[Mon22\].
- **Level 3:** scattering, reflection, numerical-aperture divergence, and radiative-transfer corrections.

A scalar dose or cure threshold should be retained only as a diagnostic. Gelation should be identified from the evolving conversion and network state.

### 5. Eulerian–Lagrangian flow coupling

The material trajectories are

``` math
\frac{d\mathbf{x}_p}{dt}=\mathbf{u}(\mathbf{x}_p,t),
\qquad
I_p(t)=I(\mathbf{x}_p(t),t).
```

Eulerian fields are still required for oxygen and other diffusing species; Lagrangian markers provide the illumination and reaction histories of material elements.

\[Slu19\] supplies the essential reduced benchmark. In a plug-flow illumination zone of length $`L`$,

``` math
t_L=\frac{L}{U},
\qquad
\frac{t_{\mathrm{uv}}}{t_L}=\frac{Ut_{\mathrm{uv}}}{L}.
```

The critical processing speed scales as

``` math
U_c\approx\frac{L\,\phi\,\varepsilon I c_{\mathrm{PI},0}}{c_{\mathrm{O_2},0}}.
```

This expresses a key physical point: flow continuously supplies oxygenated resin, so a fluid element must consume the local oxygen inventory before appreciable polymerization occurs. The scaling and the fiber-length collapse against $`t_{\mathrm{uv}}/t_L`$ should become mandatory tests for the flow module.

The \[Slu19\] equations are a steady 1D Eulerian advection–reaction model and are equivalent to a Lagrangian material-element model under $`t=x/U`$. They are therefore the correct analytical bridge between the static reaction–diffusion model and the full flowing-resin formulation.

### 6. Mechanistic crosslink-density closure

\[Zhu20\] substantially strengthens the network stage for PEGDA. It couples spatial conversion to loop-aware gelation using a modified Macosko recursive model. The gel point depends on precursor concentration:

``` math
P_{\mathrm{gel}}=\frac{C_0}{[\mathrm{PEGDA}]_0}.
```

After gelation, the fractional elastically active network density $`\eta`$ is obtained from the probability that reacted arms connect to the infinite network, including intramolecular loops. The shear modulus is then

``` math
G(\mathbf{x},t)=\eta(\mathbf{x},t)N_A[\mathrm{PEGDA}]_0k_BT.
```

This replaces the earlier generic closure $`\nu_x=\mathcal{F}(\alpha)`$ for PEGDA. It prevents the model from treating every reacted acrylate as an elastically effective crosslink and predicts a modulus ceiling below the ideal network value.

The closure is not universal. It does not include advection, species diffusion, thermal transport, swelling, or poroelasticity. PNIPAM and thiol–ene systems will require separate functional-group and network closures.

## What the full texts establish

| Physical block | Strongest result | Consequence for the roadmap |
|:---|:---|:---|
| Projection kinematics | \[Mee20\] gives measured PSFs, frame cycling, synchronization, and scan-wise dose convolution | Add a projection engine before the optical–chemical PDE; do not represent a moving DMD only as a generic time-dependent intensity |
| Bulk kinetics | \[Dob24\] shows that radical species, chain length, free volume, and reaction-diffusion-controlled termination matter across widely different initiation rates | Use \[Dob24\] as the high-fidelity reference; do not fit one fixed $`k_t`$ across all intensities |
| Grayscale DLP | \[Mon22\] couples four species, conversion-dependent diffusion, Gaussian pixel overlap, dynamic attenuation, and a moving multilayer boundary | Replace scalar threshold models with a transient reaction–diffusion solver |
| Flow | \[Slu19\] shows that residence time, oxygen replenishment, and $`t_{\mathrm{uv}}/t_L`$ control whether fibers form | Add a plug-flow benchmark before attempting general 3D advection |
| Hydrogel network | \[Zhu20\] predicts composition-dependent gelation, loop defects, crosslink density, and modulus | Use a PEGDA-specific network closure rather than a direct empirical modulus–conversion fit |
| Moving illumination | \[Bri21, Guv22, Uzc20\] provide moving-source and grayscale property-gradient formulations | Use these for validation, but supply their illumination fields with \[Mee20\] and their kinetics with \[Mon22\] |
| Oxygen inhibition | \[Den08, Tak20, Kru16\] establish wall, interface, and dead-zone effects | Retain oxygen as a transported field, not a local correction factor |

## Staged modelling roadmap

| Stage | Model scope | Main evidence | Deliverable |
|:---|:---|:---|:---|
| 0\. Projection and material calibration | Measure PSF, pixel pitch, magnification, scan speed, frame timing, absorption, bleaching, oxygen, viscosity, and velocity | \[Mee20, Mon22, Ste23\] | Calibrated optical and material parameter set |
| 1\. High-fidelity reference | Implement the concrete local model and solver structure in Stage 1 PEGDA reaction kinetics and Stage 1 PEGDA solver implementation, using \[Dob24\] as the high-fidelity reference | \[Dob24\] | Reference conversion, radical, oxygen, and temperature responses |
| 2\. Reduced local chemistry | Derive a stable reduced model from Stage 1; test equal-dose and unequal-intensity histories | \[Dob24, Wyd14, Lin19\] | Validated reduced $`k_p`$, $`k_t`$, bleaching, and oxygen closures |
| 3\. Projection engine | Generate $`I(x,y,0,t)`$ from measured PSFs, bitmap frames, pattern motion, and pixel cycling | \[Mee20\] | Surface intensity and exposure fields with verified scan synchronization |
| 4\. Static grayscale reaction–diffusion | Couple the projection engine to the \[Mon22\] four-species PDE with dynamic attenuation | \[Mon22, Den08, Tak20\] | Grayscale conversion fields, gel fronts, and lateral transition widths |
| 5\. Moving illumination | Allow nonzero pattern velocity for moving masks or sources; retain zero-relative-velocity scrolling as a special case | \[Mee20, Bri21, Guv22, Uzc20\] | Conversion and crosslink-density maps under arbitrary exposure histories |
| 6\. Flow benchmark | Solve the \[Slu19\] plug-flow model and recover critical velocity and fiber-length scalings | \[Slu19\] | Verified residence-time and oxygen-replenishment module |
| 7\. Eulerian–Lagrangian flow | Combine the projection field, spatial PDEs, prescribed $`\mathbf{u}(\mathbf{x},t)`$, oxygen diffusion, and material trajectories | \[Slu19, Wan25h, Kru16\] | Inlet-to-outlet conversion and crosslink-density histories |
| 8\. Network mechanics | Apply \[Zhu20\] for PEGDA; add chemistry-specific closures for PNIPAM and thiol–ene systems | \[Zhu20, Zhu18, Hig20\] | Gelation, crosslink density, modulus, and swelling inputs |
| 9\. Feedback and inverse control | Optimize intensity, pattern, scan speed, and residence time for a target $`\alpha`$ or $`\eta`$ field | \[Guv22, Ste23\] | Illumination trajectories with transport compensation |

## Dimensionless groups and model-selection tests

- **Optical depth:** $`\beta=AH`$. Large values require depth-resolved optical transport.
- **Diffusive oxygen Damköhler number:** use the scaling from \[Den08\], $`\mathrm{Da}_{O}\sim \phi\varepsilon c_{\mathrm{PI}}IH^2/(D_{\mathrm{O_2}}c_{\mathrm{O_2}})`$.
- **Flow oxygen-depletion number:** $`\Pi_{\mathrm{flow}}=Lk_dc_{\mathrm{PI},0}/(Uc_{\mathrm{O_2},0})`$. Values near unity mark the transition between non-gelling and curing flow.
- **Exposure-to-residence ratio:** $`\tau_{\mathrm{exp}}=t_{\mathrm{uv}}/t_L`$. This is the main scaling variable for pulsed flow exposure \[Slu19\].
- **Pattern-to-material velocity ratio:** $`\Lambda_v=\lvert\mathbf{v}_{\mathrm{pattern}}-\mathbf{u}\rvert/U_*`$. It distinguishes a stationary projected feature, a swept exposure, and a pattern moving relative to the resin.
- **Species Péclet number:** $`\mathrm{Pe}_i=UH/D_i`$. This determines whether oxygen, monomer, or radicals are diffusion- or advection-dominated.
- **Photobleaching number:** $`\mathrm{Da}_{\mathrm{bleach}}=k_bIt`$. Large values require dynamic optical attenuation.
- **Frame-resolution ratio:** $`\chi=\Delta t_{\mathrm{frame}}/t_{\mathrm{chem}}`$. If $`\chi`$ is not small, discrete DMD frame timing must be retained rather than replaced by a continuous scan approximation.
- **Thermal number:** compare reaction heat generation with thermal conduction and boundary cooling before activating the energy equation.
- **Gelation or viscosity number:** track proximity to the network transition where fluid flow can no longer be treated independently of the evolving solid.

## Validation ladder

| Validation target | Measurement | Model component tested |
|:---|:---|:---|
| Projection field | Camera or radiometer maps of single-pixel PSF, frame timing, and scan exposure | Projection engine and synchronization \[Mee20\] |
| Local kinetics | FTIR conversion versus intensity, time, and temperature | Reduced chemistry against \[Dob24\], \[Wyd14\], and \[Lin19\] |
| Optical evolution | UV–visible bleaching and cure-depth profiles | Dynamic attenuation and photoinitiator depletion \[Mon22, Ste23\] |
| Grayscale interface | Conversion transition width and feature growth | Pixel overlap, species diffusion, and gel-front evolution \[Mon22\] |
| Oxygen field | Induction time, inhibition thickness, and wall-adjacent conversion | Oxygen reaction–diffusion \[Den08, Tak20, Kru16\] |
| Flow threshold | Critical velocity, fiber onset, and fiber length | Residence-time and oxygen-replenishment scaling \[Slu19\] |
| Moving illumination | Conversion or modulus versus pattern velocity, scan position, and scan speed | Relative pattern/material kinematics \[Mee20, Bri21, Uzc20\] |
| Network state | Gel point, modulus, loop density, and swelling | PEGDA network closure \[Zhu20\] |

## Recommended first implementation

The most defensible first solver is now a two-dimensional operator-split PEGDA model with an explicit projection engine:

1.  Measure the single-pixel PSF, pixel pitch, magnification, scan speed, and frame timing.
2.  Generate $`I(x,y,0,t)`$ from bitmap frames, measured PSFs, and the desired pattern velocity. Use synchronized frame cycling to impose $`\mathbf{v}_{\mathrm{pattern}}=0`$ when appropriate.
3.  Propagate that field with dynamic attenuation from \[Mon22\].
4.  Advance photoinitiator, radicals, oxygen, and monomer with reaction–diffusion equations and conversion-dependent diffusivity.
5.  Benchmark the local chemistry against the reduced form of \[Dob24\].
6.  Add prescribed flow and verify the plug-flow limits from \[Slu19\].
7.  Advect material markers through the Eulerian velocity field to record individual intensity histories.
8.  Apply the loop-aware \[Zhu20\] closure to obtain gelation and crosslink density.

PEGDA should be the first chemistry because \[Mon22\] and \[Zhu20\] directly cover its spatial kinetics and network evolution. PNIPAM can reuse the transport architecture but requires different monomer diffusivity, solvent, temperature, and network parameters. Thiol–ene systems should remain a separate chemistry branch because oxygen tolerance and step-growth network formation change the kinetic closure.

The full \[Dob24\] chemistry, Navier–Stokes feedback, scattering, shrinkage, swelling, and inverse control should be added only after projection fields, conversion fields, and oxygen-inhibited boundaries are correct. This keeps the model interpretable and makes each added physical effect testable.

## Remaining gap

The central unresolved problem is now precise: coupling the measured projection operator from \[Mee20\], transient grayscale reaction–diffusion from \[Mon22\], residence-time oxygen depletion from \[Slu19\], and loop-aware network evolution from \[Zhu20\] under arbitrary pattern motion and a spatially varying velocity field. The roadmap no longer lacks the moving-mask kinematic foundation; the remaining gap is the fully coupled chemistry–transport implementation.

---

## References

\[Wyd14\] J. W. Wydra, N. Cramer, J. Stansbury, and C. Bowman, “The reciprocity law concerning light dose–relationships applied to BisGMA/TEGDMA photopolymers: Theoretical analysis and experimental characterization,” *Dental materials : official publication of the Academy of Dental Materials*, vol. 30, pp. 605–612, Mar. 2014, doi: [10.1016/j.dental.2014.02.021](https://doi.org/10.1016/j.dental.2014.02.021).

\[Lin19\] J.-T. Lin, H. Liu, K.-T. Chen, and D.-C. Cheng, “Modeling the Kinetics, Curing Depth, and Efficacy of Radical-Mediated Photopolymerization: The Role of Oxygen Inhibition, Viscosity, and Dynamic Light Intensity,” *Frontiers in Chemistry*, vol. 7, Nov. 2019, doi: [10.3389/fchem.2019.00760](https://doi.org/10.3389/fchem.2019.00760).

\[Mon22\] S. Montgomery, C. M. Hamel, J. Skovran, and H. Q, “A reaction-diffusion model for grayscale digital light processing 3D printing,” Mar. 01, 2022. doi: [10.1016/j.eml.2022.101714](https://doi.org/10.1016/j.eml.2022.101714).

\[Dob24\] A. L. Dobson and C. N. Bowman, “A Comprehensive, Multidimensional First‐Principles Model for Free‐Radical Photopolymerizations in Bulk and Thin Films,” *Advanced Functional Materials*, vol. 34, Feb. 2024, doi: [10.1002/adfm.202312607](https://doi.org/10.1002/adfm.202312607).

\[Slu19\] M. Slutzky, H. Stone, and J. Nunes, “A quantitative study of the effect of flow on the photopolymerization of fibers.” *Soft matter*, Nov. 2019, doi: [10.1039/c9sm01485c](https://doi.org/10.1039/c9sm01485c).

\[Zhu20\] H. Zhu, X. Yang, G. Genin, T. Lu, F. Xu, and M. Lin, “Modeling the mechanics, kinetics, and network evolution of photopolymerized hydrogels,” Sep. 01, 2020. doi: [10.1016/j.jmps.2020.104041](https://doi.org/10.1016/j.jmps.2020.104041).

\[Mee20\] V. Meenakshisundaram, L. Sturm, and C. Williams, “Modeling A Scanning-Mask Projection Vat Photopolymerization System For Multiscale Additive Manufacturing,” May 01, 2020. doi: [10.1016/j.jmatprotec.2019.116546](https://doi.org/10.1016/j.jmatprotec.2019.116546).

\[Bri21\] R. Brighenti, M. P. Cosma, L. Marșavina, A. Spagnoli, and M. Terzano, “Multiphysics modelling of the mechanical properties in polymers obtained via photo-induced polymerization,” Jul. 30, 2021. doi: [10.1007/s00170-021-07273-2](https://doi.org/10.1007/s00170-021-07273-2).

\[Guv22\] E. Guven, Y. Karpat, and M. Çakmakci, “Improving the Dimensional Accuracy of Micro Parts 3d Printed with Projection-Based Continuous Vat Photopolymerization Using a Model-Based Grayscale Optimization Method,” 2022. doi: [10.2139/ssrn.4055141](https://doi.org/10.2139/ssrn.4055141).

\[Uzc20\] A. C. Uzcategui *et al.*, “Microscale Photopatterning of Through‐Thickness Modulus in a Monolithic and Functionally Graded 3D‐Printed Part,” *Small Science*, vol. 1, Dec. 2020, doi: [10.1002/smsc.202000017](https://doi.org/10.1002/smsc.202000017).

\[Den08\] D. Dendukuri, P. Panda, R. Haghgooie, J. M. Kim, T. Hatton, and P. Doyle, “Modeling of Oxygen-Inhibited Free Radical Photopolymerization in a PDMS Microfluidic Device,” Oct. 22, 2008. doi: [10.1021/MA801219W](https://doi.org/10.1021/MA801219W).

\[Tak20\] K. Taki, “A Simplified 2D Numerical Simulation of Photopolymerization Kinetics and Oxygen Diffusion–Reaction for the Continuous Liquid Interface Production (CLIP) System,” *Polymers*, vol. 12, Apr. 2020, doi: [10.3390/polym12040875](https://doi.org/10.3390/polym12040875).

\[Kru16\] K. Krutkramelis, B. Xia, and J. Oakey, “Monodisperse Polyethylene Glycol Diacrylate Hydrogel Microsphere Formation by Oxygen-Controlled Photopolymerization in a Microfluidic Device,” *Lab on a chip*, vol. 16, pp. 1457–1465, Apr. 2016, doi: [10.1039/c6lc00254d](https://doi.org/10.1039/c6lc00254d).

\[Ste23\] L. M. Stevens *et al.*, “Counting All Photons: Efficient Optimization of Visible Light 3D Printing,” *Advanced Materials Technologies*, vol. 8, Apr. 2023, doi: [10.1002/admt.202300052](https://doi.org/10.1002/admt.202300052).

\[Wan25h\] X. Wang, A. S. K. Kho, J. Liu, T. Mao, M. Gilchrist, and N. Zhang, “Mechanistic Modelling of Coupled UV Energy Penetration and Resin Flow Dynamics in Digital Light Processing (DLP)-Based Microfluidic Chip Printing,” *Micromachines*, vol. 16, Jan. 2025, doi: [10.3390/mi16020115](https://doi.org/10.3390/mi16020115).

\[Zhu18\] H. Zhu, X. Yang, G. Genin, T. Lu, F. Xu, and M. Lin, “The Relationship between Thiol-acrylate Photopolymerization Kinetics and Hydrogel Mechanics: An Improved Model Incorporating Photobleaching and Thiol-Michael Addition,” *Journal of the mechanical behavior of biomedical materials*, vol. 88, pp. 160–169, Aug. 2018, doi: [10.1016/j.jmbbm.2018.08.013](https://doi.org/10.1016/j.jmbbm.2018.08.013).

\[Hig20\] C. I. Higgins, J. Killgore, F. DelRio, S. Bryant, and R. McLeod, “Photo-tunable hydrogel mechanical heterogeneity informed by predictive transport kinetics model,” *Soft matter*, vol. 16, pp. 4131–4141, Mar. 2020, doi: [10.1039/d0sm00052c](https://doi.org/10.1039/d0sm00052c).
