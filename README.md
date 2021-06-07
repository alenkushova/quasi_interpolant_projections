# Quasi interpolant projections
Implementation of local quasi interpolant projection .

## Description 
This files provide 1D and 2D projections of functions into spline spaces according
with the De Rham spline complexes commutative diagrams. The projections are local and built with the 
Algorithm of Lee, Lyche, Morken. We deal both with Dirichlet or periodic boundary conditions.

### Content
For Dirichlet boundary conditions:
- **Lyche_1D.m**  is the 1D quasi interpolant projection defined by Lee Lyche Morken. (we use open knot vectors)
- **Lyche_c_1D.m** is the commutative projection of the previous one based on Dirichlet boundary conditions (H^1_0 and L2 De Rham complex). 
- **Lyche_c_1D_new.m** is the commutative projection in 1D De Rham complex between H1 and L2 (it projects with interpolation at the boundary).
- **Lyche_2D_Drch.m** is the tensor product of 2 univariate Lyche_1D projections.
- **Lyche_c_2D_Drch.m** is the tensor product between Lyche_1D and Lyche_c_1D in order to project in subspace of H(div, Omega)
- **Lyche_c_2D_Drch_new.m** is the tensor product as before wih Lyche_c_1D_new file.

For periodic boundary conditions:
- **Lyche_1D_Periodic.m**  is the 1D quasi interpolant projection for periodic functions. (we use closed knot vectors)
- **Lyche_c_1D_Periodic.m** is the commutative projection of the previous one based on periodic boundary conditions. 
- **Lyche_2D_Periodic.m** is the tensor product of 2 univariate Lyche_1D_Periodic projections.
- **Lyche_c_2D_Periodic.m** as previous but between Lyche_1D_Periodic and Lyche_c_1D_Periodic.

## Important
This lybrary requires [GeoPDEs periodic](https://github.com/rafavzqz/geopdes/tree/periodic) software.

### Examples 
Two examples are available. Run the Matlab file **approximation_examples.m** to give a look. 
