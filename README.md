# Lithium Diffusion and Stress in Silicon Battery Electrodes

Silicon can store a large amount of lithium, making it attractive for rechargeable batteries. However, it swells as lithium enters. Uneven swelling creates internal stress, which can contribute to damage and also change how quickly lithium moves.

This project aims to reproduce a published model of that feedback in silicon nanowire electrodes. It uses a phase-field formulation, which describes material behaviour through continuous fields and free energy, and the finite element method (FEM), which divides the electrode into small elements to solve the equations numerically.

**Status:** Agreement with the paper's results has not yet been established. The model studies lithium insertion into an individual electrode, rather than a complete battery.

## Questions

1. When does stress peak at a location, and how does this relate to the arrival of the lithium concentration front?
2. Which measurements best describe how quickly lithium spreads through the electrode?
3. How much does stress change lithium transport compared with ordinary diffusion, and why?
4. Does the charging process have distinct stages, and how do these depend on electrode size?
5. How does the fastest local concentration change relate to the overall state-of-charge rate?

## Main assumptions

- **Geometry:** a solid cylindrical silicon nanowire with circular symmetry; a quarter of its cross-section can represent the full cross-section.
- **Deformation:** small strains and linear, isotropic elasticity; no permanent plastic deformation or crack growth.
- **Axial constraint:** plane strain, so expansion along the wire is prevented, although axial stress can develop.
- **Mechanical balance:** deformation adjusts much faster than diffusion; inertia and body forces such as gravity are neglected.
- **Constant properties:** temperature, elastic properties, partial molar volume and diffusivity remain constant during lithiation. Stress affects the transport driving force even though diffusivity itself is constant.
- **Lithium entry:** surface reactions are sufficiently fast that the outer surface immediately reaches a fixed concentration, $c_{\max}$. Transport inside silicon limits the insertion rate.
- **Initial and mechanical conditions:** the interior starts at uniform concentration $c_0$; the curved outer surface is mechanically free, with symmetry conditions on the cut edges.
- **Numerical simplification:** following the paper, gradient-energy contributions are omitted from the simulated transport equation.

## Governing equations

The unknowns are lithium concentration $c(\mathbf{x},t)$ and displacement $\mathbf{u}(\mathbf{x},t)$.

**Free-energy starting point** combines chemical, elastic and concentration-gradient energy:

$$
F=\int_V\left(f_{\mathrm{chem}}+f_{\mathrm{el}}+\frac{\kappa}{2}|\nabla c|^2\right)\,dV.
$$

The gradient contribution motivates the full phase-field formulation, but its transport terms are dropped for reproduction of the paper's simulations.

**Deformation and swelling** determine elastic stress:

$$
\boldsymbol\varepsilon=\frac12\left[\nabla\mathbf u+(\nabla\mathbf u)^T\right],\qquad
\boldsymbol\varepsilon^{\mathrm{ch}}=\frac{\Omega(c-c_{\mathrm{ref}})}{3}\mathbf I,
$$

$$
\boldsymbol\sigma=\mathbb C:(\boldsymbol\varepsilon-\boldsymbol\varepsilon^{\mathrm{ch}}).
$$

**Mechanical equilibrium** balances internal forces. Hydrostatic stress is the average normal stress:

$$
\nabla\cdot\boldsymbol\sigma=\mathbf0,\qquad
\sigma_h=\frac{\sigma_{xx}+\sigma_{yy}+\sigma_{zz}}3.
$$

**Lithium flux** includes concentration and stress gradients:

$$
\mathbf J=-D\left[\nabla c-\frac{\Omega c}{R_gT}\nabla\sigma_h\right].
$$

**Lithium conservation** updates concentration as lithium moves:

$$
\frac{\partial c}{\partial t}+\nabla\cdot\mathbf J=0.
$$

Here, $\boldsymbol\varepsilon$ is strain, $\mathbb C$ is elastic stiffness, $\Omega$ is partial molar volume, $c_{\mathrm{ref}}$ is reference concentration, $\mathbf I$ is the identity tensor, $D$ is diffusivity, $R_g$ is the gas constant and $T$ is temperature. The colon denotes tensor contraction. Physical concentration is expressed in $\mathrm{mol/m^3}$.

**Symmetric-cylinder check:** under the paper's assumptions, eliminating mechanics gives

$$
\frac{\partial c}{\partial t}
=\nabla\cdot\left[D(1+\Theta c)\nabla c\right],\qquad
\Theta=\frac{2\Omega^2E}{9R_gT(1-\nu)},
$$

where $E$ is Young's modulus and $\nu$ is Poisson's ratio. The paper relates diffusivity to molecular mobility $M$ through $D=Mk_BT$, where $k_B$ is Boltzmann's constant.

Comparisons will examine concentration and stress profiles, front movement, local concentration rate $\partial c/\partial t$, its spatial maximum, and state of charge and its rate. Verification will include the ordinary-diffusion limit and checks that refining the mesh and time step no longer materially changes the results.

## Reference

Bernhard Eidel and Masoud Ahmadi, *Insights into the Interplay of Diffusion and Stress in the Intercalation Dynamics of a Li-Ion Battery Electrode by Phase-Field Finite Element Modeling*. Manuscript dated 13 January 2020. Motivation and questions follow Section 1; assumptions and equations follow Sections 2-4.
