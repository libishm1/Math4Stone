# The discrete Laplace-Beltrami operator and Poisson equations

Prereqs: [[17_discrete-curvature]], [[14_linear-least-squares]] => this => [[20_sdf-signed-heat]], [[27_gnn-message-passing]]

## Why it matters for the stone scan->CAD->IFC / CADFit workflow

A raw stone scan arrives as a noisy, dense triangle mesh with high-frequency sensor noise that no NURBS fitter should ever see directly. The cotangent Laplacian is the operator you use to smooth and fair that mesh while preserving its true low-frequency geometry, so downstream `GeomAPI_PointsToBSplineSurface` fits a clean surface instead of chasing scanner speckle. The same operator computes the mean-curvature field that drives feature-aware segmentation, telling you where to split one stone block into separate B-spline patches versus where a single chart is valid. When you build a UV parameterization to resample the irregular cloud onto the ordered grid OCCT needs, you solve a Poisson or harmonic system whose matrix is exactly this Laplacian. Finally, the Signed Heat method that recovers a clean signed distance field from broken scan geometry is built on the very same cotan stiffness and mass matrices, so mastering this one operator unlocks smoothing, curvature, parameterization, and reconstruction at once.

## Intuition

Think of the value at each vertex as the temperature of a tiny metal plate welded to the surface. The Laplacian at a vertex measures how much that vertex differs from the weighted average of its neighbors. If a vertex is hotter than its surroundings, the Laplacian is negative and heat wants to flow out; if colder, the Laplacian is positive and heat flows in. Smoothing a mesh is letting heat diffuse: bumps cool down, dents fill in, and the surface relaxes toward something fair.

The cotangent weights are the one twist that makes this work on an irregular mesh instead of a uniform grid. A naive average using equal weights ignores that triangles are stretched and skewed. The cotangents of the two angles facing each edge are precisely the correction that makes the discrete operator agree with the smooth Laplace-Beltrami operator as the mesh refines. That is why this single matrix is called the Swiss army knife of geometry processing: smoothing, curvature, distance, and parameterization are all just different right-hand sides for the same operator.

## The math (first principles)

### Setup and conventions

Work on a triangle mesh $M = (V, E, F)$ embedding a surface in $\mathbb{R}^3$. Let $f : V \to \mathbb{R}$ be a scalar function sampled at vertices, stored as a column vector $\mathbf{f} \in \mathbb{R}^{|V|}$. For a vector-valued quantity such as positions, stack them as a matrix $\mathbf{V} \in \mathbb{R}^{|V| \times 3}$ and apply operators column-wise. Units are whatever the mesh uses; for stone scans keep meters at the quarry scale and millimeters at the shop scale, and never mix them inside one solve. The mesh must be manifold and consistently oriented; the cotan formula assumes triangular faces.

The smooth Laplace-Beltrami operator $\Delta$ generalizes the flat Laplacian $\nabla^2$ to curved surfaces. Convention warning: geometers usually define $\Delta = \operatorname{div}\operatorname{grad}$, which is negative semidefinite, while analysts often use $-\Delta$. This document and the libigl library use the geometry convention, where the discrete stiffness matrix has a negative diagonal.

### The cotangent formula

For an interior vertex $i$ with one-ring neighbors $N(i)$, the cotangent discretization of the Laplace-Beltrami operator is

$$
(\Delta f)_i = \frac{1}{2 A_i} \sum_{j \in N(i)} \left( \cot \alpha_{ij} + \cot \beta_{ij} \right) (f_j - f_i).
$$

Here $\alpha_{ij}$ and $\beta_{ij}$ are the two interior angles opposite the edge $(i,j)$, one in each triangle sharing that edge, and $A_i$ is the area of the dual cell around vertex $i$. Each edge weight is

$$
w_{ij} = \frac{1}{2}\left( \cot \alpha_{ij} + \cot \beta_{ij} \right).
$$

This formula was introduced by Pinkall and Polthier and refined by Meyer, Desbrun, Schroder and Barr. It converges to the smooth operator as the mesh refines, which is why it is preferred over the uniform (graph) Laplacian.

### Matrix form: stiffness and mass

Split the operator into two sparse matrices. The cotangent stiffness (or weak Laplacian) matrix $C \in \mathbb{R}^{|V| \times |V|}$ is

$$
C_{ij} = \begin{cases}
-w_{ij} = -\tfrac{1}{2}(\cot \alpha_{ij} + \cot \beta_{ij}) & j \in N(i), \\[4pt]
\displaystyle \sum_{k \in N(i)} w_{ik} & j = i, \\[4pt]
0 & \text{otherwise.}
\end{cases}
$$

Wait: in the libigl/geometry sign convention the signs flip so that $C$ is negative semidefinite. Concretely `igl::cotmatrix` builds $L$ with positive off-diagonal weights $+w_{ij}$ and a negative diagonal $-\sum_k w_{ik}$, so each row sums to zero and constant functions lie in the null space. We adopt that convention: call it $L$. Then $L \mathbf{f}$ is the integrated (area-weighted) Laplacian, not the pointwise one.

The mass matrix $M \in \mathbb{R}^{|V| \times |V|}$ is diagonal with $M_{ii} = A_i$, the dual area of vertex $i$. Two common choices:

- Barycentric: $A_i = \tfrac{1}{3}\sum_{t \ni i} \operatorname{area}(t)$, one third of each incident triangle. Always positive, simplest.
- Voronoi (circumcentric, Meyer): $A_i$ is the Voronoi cell area, using the mixed formula that falls back to barycentric on obtuse triangles to stay positive. More accurate.

The pointwise Laplace-Beltrami operator is then recovered as

$$
(\Delta f) = M^{-1} L \mathbf{f}, \qquad \text{i.e.} \quad (\Delta f)_i = \frac{1}{A_i}(L\mathbf{f})_i,
$$

which is exactly the cotangent formula above. Because $M$ is diagonal, $M^{-1}$ is trivial. Keep $L$ and $M$ separate: $L$ is symmetric, $M$ encodes the inner product, and most solves want the symmetric $L$ on the left.

### The Poisson equation

The Poisson problem is: find $u$ with $\Delta u = \rho$ for a given source $\rho$. Discretized with the cotan operator and mass matrix,

$$
L \mathbf{u} = M \boldsymbol{\rho}.
$$

The factor $M$ on the right turns the pointwise source $\rho$ into an integrated source over each dual cell, matching the integrated left side $L\mathbf{u}$. Because $L$ has a one-dimensional null space (constants), the system is solvable only when $\sum_i A_i \rho_i = 0$ (the source integrates to zero), and the solution is unique up to an additive constant. Pin one vertex or add the constraint $\mathbf{1}^\top M \mathbf{u} = 0$ to remove the freedom.

### Mean curvature from the Laplacian

Applying the Laplacian to the embedding $\mathbf{V}$ (the position function) gives the mean curvature normal:

$$
\Delta \mathbf{x} = -2 H \mathbf{n}, \qquad \text{discretely} \quad \mathbf{H}\mathbf{n} = -M^{-1} L \mathbf{V}.
$$

The row norm $\lVert (M^{-1}L\mathbf{V})_i \rVert$ gives $2|H_i|$, a per-vertex mean-curvature magnitude you can threshold for segmentation.

### Implicit smoothing (curvature flow)

Mean curvature flow moves each point along $-H\mathbf{n}$. Explicit Euler $\mathbf{V}^{k+1} = \mathbf{V}^k + \delta\, M^{-1} L \mathbf{V}^k$ is unstable for large $\delta$. The backward-Euler (implicit) step of Desbrun, Meyer, Schroder, Barr is stable for any $\delta$:

$$
(M - \delta L)\, \mathbf{V}^{k+1} = M \mathbf{V}^{k}.
$$

This is one sparse symmetric positive-definite solve per step. Larger $\delta$ smooths more aggressively. Recompute $M$ (and optionally $L$) each step if you want the geometrically correct flow, or freeze $L$ for cheap fairing.

## Worked example

Take the smallest mesh where a cotangent weight is interesting: a single edge $(i,j)$ shared by two triangles, vertex $i$ at the center of a small fan. To keep arithmetic checkable, use one flat triangle and compute one edge weight by hand.

Triangle with vertices $P_i = (0,0)$, $P_j = (4,0)$, $P_k = (0,3)$. The angle opposite edge $(i,j)$ is the angle at $P_k$. Form the two edge vectors from $P_k$:

$$
\mathbf{u} = P_i - P_k = (0,-3), \qquad \mathbf{v} = P_j - P_k = (4,-3).
$$

The cotangent of the angle between two vectors is $\cot\theta = \dfrac{\mathbf{u}\cdot\mathbf{v}}{\lVert \mathbf{u}\times\mathbf{v}\rVert}$.

Dot product: $\mathbf{u}\cdot\mathbf{v} = (0)(4) + (-3)(-3) = 9.$

Cross magnitude (2D scalar cross): $\lVert \mathbf{u}\times\mathbf{v}\rVert = |(0)(-3) - (-3)(4)| = |12| = 12.$

So $\cot\alpha_{ij} = 9/12 = 0.75$. With a single triangle the edge $(i,j)$ is a boundary edge, so $\beta$ is absent and the half-weight is $w_{ij} = \tfrac12(0.75) = 0.375$.

Now add the mirror triangle $P_{k'} = (0,-3)$ below the edge, giving the same $\cot\beta_{ij} = 0.75$ by symmetry. The full edge weight is

$$
w_{ij} = \tfrac12(\cot\alpha_{ij} + \cot\beta_{ij}) = \tfrac12(0.75 + 0.75) = 0.75.
$$

Sanity check on the stiffness diagonal: if $i$ had only this one neighbor, $L_{ii} = -w_{ij} \cdot (\text{wait, sum convention})$. In libigl convention $L_{ij} = +0.75$ and $L_{ii} = -0.75$, so the row $(L\mathbf{f})_i = 0.75(f_j - f_i)$ sums to zero. Good: a constant $f$ gives $(L\mathbf{f})_i = 0$.

Barycentric area at $i$: triangle area is $\tfrac12 \cdot 4 \cdot 3 = 6$, so each incident triangle contributes $6/3 = 2$; with two triangles $A_i = 4$. The pointwise Laplacian contribution is $(\Delta f)_i = \tfrac{1}{A_i}(L\mathbf f)_i = \tfrac{1}{4}\cdot 0.75 (f_j - f_i) = 0.1875(f_j - f_i)$. Every number here is reproducible with a calculator.

## Code

Idiomatic C++ using Eigen and libigl. libigl ships the standard cotan operators, so production code assembles them rather than hand-rolling cotangents. The snippet builds $L$ and $M$, computes mean curvature, runs one implicit smoothing step, and solves a Poisson system.

```cpp
// laplace_poisson.cpp
// Build the cotan Laplacian, do curvature, smoothing, and a Poisson solve.
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <igl/cotmatrix.h>       // C: cotan stiffness, libigl sign (neg diagonal)
#include <igl/massmatrix.h>      // M: lumped (diagonal) mass matrix
#include <igl/invert_diag.h>
#include <iostream>

int main() {
    Eigen::MatrixXd V;   // |V| x 3 vertex positions  (the embedding x)
    Eigen::MatrixXi F;   // |F| x 3 triangle indices
    // ... load V, F from your scan mesh here ...

    // --- assemble operators -------------------------------------------------
    Eigen::SparseMatrix<double> L, M, Minv;
    igl::cotmatrix(V, F, L);                                 // L  (= the matrix C)
    igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_VORONOI, M);  // M_ii = Voronoi area A_i
    igl::invert_diag(M, Minv);                               // M^{-1}, cheap (diagonal)

    // --- mean curvature normal:  H n = -M^{-1} L x --------------------------
    Eigen::MatrixXd HN = -(Minv * (L * V));     // row i = mean-curvature normal at i
    Eigen::VectorXd H  = 0.5 * HN.rowwise().norm(); // |H_i| per vertex (segmentation cue)

    // --- one implicit smoothing step:  (M - delta L) V' = M V ---------------
    const double delta = 1e-3;                  // larger delta => more smoothing
    Eigen::SparseMatrix<double> A = (M - delta * L); // SPD: solvable with Cholesky
    Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> chol(A);
    Eigen::MatrixXd Vs = chol.solve(M * V);     // smoothed positions, column-wise

    // --- Poisson solve:  L u = M rho  (one scalar field) --------------------
    // L is rank-deficient (constants in null space). Pin vertex 0 to fix the DOF.
    Eigen::VectorXd rho = Eigen::VectorXd::Zero(V.rows()); // your source term
    // ... fill rho so that 1^T M rho == 0 (zero net source) ...
    Eigen::VectorXd b = M * rho;
    Eigen::SparseMatrix<double> Lp = L;         // copy to pin a vertex
    // Dirichlet pin: force u[0] = 0 by replacing row/col 0 with identity.
    for (int k = 0; k < Lp.outerSize(); ++k)
        for (Eigen::SparseMatrix<double>::InnerIterator it(Lp, k); it; ++it)
            if (it.row() == 0 || it.col() == 0) it.valueRef() = (it.row()==it.col()) ? 1.0 : 0.0;
    b(0) = 0.0;
    Eigen::SparseLU<Eigen::SparseMatrix<double>> lu(Lp);
    Eigen::VectorXd u = lu.solve(b);            // u defined up to the pinned constant

    std::cout << "max |H| = " << H.maxCoeff() << "\n";
    return 0;
}
```

Symbol map: `L` is the stiffness matrix $L$ (libigl negative-diagonal convention), `M` is the mass matrix $M$ with Voronoi areas $A_i$, `Minv` is $M^{-1}$, `HN` is $-M^{-1}L\mathbf{V} = H\mathbf{n}$, the smoothing line solves $(M-\delta L)\mathbf{V}' = M\mathbf{V}$, and the Poisson block solves $L\mathbf{u} = M\boldsymbol\rho$ with one pinned vertex.

A short Python version clarifies the same call sequence with libigl bindings.

```python
import igl
import numpy as np
import scipy.sparse.linalg as spla

# v: (n,3) float, f: (m,3) int  -- your scan mesh
L = igl.cotmatrix(v, f)                                  # stiffness, neg diagonal
M = igl.massmatrix(v, f, igl.MASSMATRIX_TYPE_VORONOI)    # diagonal mass, A_i
Minv = igl.invert_diag(M)

HN = -Minv.dot(L.dot(v))          # mean-curvature normal H n
H  = 0.5 * np.linalg.norm(HN, axis=1)

delta = 1e-3
vs = spla.spsolve(M - delta * L, M.dot(v))   # implicit smoothing step
```

## Connections

- [[17_discrete-curvature]]: mean curvature is literally $\mathbf{H}\mathbf{n} = -M^{-1}L\mathbf{x}$, so the Laplacian is the operator that produces the curvature you learned to measure; Gaussian curvature uses the angle-defect, but mean curvature is pure Laplacian.
- [[14_linear-least-squares]]: every use here ends in a sparse linear solve $A\mathbf{x}=\mathbf{b}$ with $A$ symmetric positive (semi)definite, and fairing-with-constraints is a regularized least-squares problem; the prereq gives you the normal-equation and Cholesky machinery.
- [[20_sdf-signed-heat]]: the Signed Heat method of Feng and Crane diffuses with the heat operator and integrates a field by solving a Poisson equation $L\mathbf{u}=M\boldsymbol\rho$, reusing exactly the matrices built here to get a signed distance field from broken scans.
- [[27_gnn-message-passing]]: the cotan Laplacian is a weighted message-passing operator on the mesh graph; $(L\mathbf{f})_i$ aggregates neighbor differences with learned-or-fixed weights, which is the bridge from classical geometry processing to graph neural networks.

## Pitfalls and failure modes

- Sign convention. The single most common bug. libigl `cotmatrix` is negative semidefinite, so mean curvature needs the minus sign in $-M^{-1}L\mathbf{V}$, and implicit smoothing uses $M-\delta L$ (plus, because $-L$ is positive). If you flip the convention, your surface inflates and explodes instead of smoothing. Check that a constant function gives $L\mathbf{1}=\mathbf{0}$.
- Obtuse and degenerate triangles. A very obtuse angle makes its cotangent large and negative, so edge weights can go negative; the operator loses the maximum principle and smoothing can overshoot. Use the Voronoi mass matrix with the mixed (obtuse-safe) area, and consider an intrinsic Delaunay flip to restore positive weights.
- Near-degenerate triangles. Slivers from raw scans give cotangents that blow up numerically (division by a near-zero cross product). Clean or remesh first; never feed raw scanner triangulation straight into the cotan formula.
- Rank deficiency. $L$ always has constants in its null space, so $L\mathbf{u}=\mathbf{b}$ is singular. A direct solver fails or returns garbage unless you pin a vertex or add a mean-zero constraint, and the source must satisfy $\mathbf{1}^\top M\boldsymbol\rho=0$.
- Units and scale. $L$ scales like $1/\text{length}^2$ and $M$ like $\text{length}^2$. Mixing meters and millimeters in one mesh makes $M-\delta L$ ill-conditioned and picks a meaningless $\delta$. Keep one unit per solve and choose $\delta$ relative to mean edge length squared.
- Lumped vs full mass. The diagonal (lumped) $M$ is standard and keeps $M^{-1}$ trivial, but it is an approximation of the full FEM mass matrix. For high-accuracy spectral work the full $M$ matters; for smoothing and curvature the lumped one is fine.
- Boundary vertices. On an open mesh (a scanned slab patch), boundary edges have only one opposite angle. Handle boundaries explicitly or set Neumann/Dirichlet conditions; ignoring them corrupts curvature near the rim where your fitting tolerance matters most.

## Practice

1. By hand, build the $3\times 3$ stiffness matrix $L$ for a single triangle with vertices $(0,0)$, $(1,0)$, $(0,1)$. Verify each row sums to zero and the matrix is symmetric. (Inspectable: print the $3\times 3$ matrix.)
2. Implement the cotangent weight $w_{ij}=\tfrac12(\cot\alpha+\cot\beta)$ from scratch in C++ for one edge given the four vertex coordinates, and confirm it matches `igl::cotmatrix` on a two-triangle quad. (Output: the scalar weight and the libigl entry side by side.)
3. Load a noisy mesh, run 1, 5, and 20 implicit smoothing steps with fixed $\delta$, and save each result. Observe that high-frequency noise vanishes first while overall shape is preserved. (Output: three meshes.)
4. Compute the per-vertex mean curvature $|H_i|$ on a sphere mesh and check it concentrates near the analytic value $1/r$. (Output: a histogram of $|H_i|$.)
5. Solve $L\mathbf{u}=M\boldsymbol\rho$ with $\rho$ a single positive spike at one vertex balanced by a uniform negative background so the source is mean-zero, pin one vertex, and visualize $u$ as a smooth bump (a discrete Green's function). (Output: colored mesh.)
6. Build a harmonic UV parameterization by solving $L\mathbf{u}=\mathbf{0}$ with the boundary pinned to a circle, then verify the interior map has no flipped triangles on a simple disk mesh. (Output: the UV layout.)

## References

- Keenan Crane, "Discrete Differential Geometry: An Applied Introduction" (CMU DDG course notes and the "Laplacian: Swiss Army Knife" assignment). https://www.cs.cmu.edu/~kmcrane/Projects/DDG/
- Keenan Crane, "The n-dimensional cotangent formula." https://www.cs.cmu.edu/~kmcrane/Projects/Other/nDCotanFormula.pdf
- Meyer, Desbrun, Schroder, Barr, "Discrete Differential-Geometry Operators for Triangulated 2-Manifolds" (Voronoi/mixed area, cotan operator), 2003.
- Pinkall and Polthier, "Computing Discrete Minimal Surfaces and Their Conjugates" (origin of cotan weights), 1993.
- Bobenko and Springborn, "A discrete Laplace-Beltrami operator for simplicial surfaces," Discrete & Computational Geometry, 2007.
- Desbrun, Meyer, Schroder, Barr, "Implicit Fairing of Irregular Meshes using Diffusion and Curvature Flow," SIGGRAPH 1999 (implicit smoothing $(M-\delta L)\mathbf{V}'=M\mathbf{V}$).
- libigl tutorial, Chapter 1 (cotmatrix, massmatrix, curvature, smoothing). https://libigl.github.io/libigl-python-bindings/tut-chapter1/ and https://libigl.github.io/tutorial/
- Feng and Crane, "A Heat Method for Generalized Signed Distance," ACM TOG 43(4), SIGGRAPH 2024. https://nzfeng.github.io/research/SignedHeatMethod/
- OpenCASCADE Reference Manual, `GeomAPI_PointsToBSplineSurface` (downstream surface fitting). https://dev.opencascade.org/doc/overview/html/
- Geogram documentation (Voronoi/remeshing infrastructure). https://github.com/BrunoLevy/geogram

## Flashcard seeds

- Q: Derive the cotangent of the angle $\alpha$ at the vertex opposite an edge, written using only the two edge vectors $\mathbf{u},\mathbf{v}$ from that vertex (derive, step 1 of 2). <br> <b>Formula:</b> $\cot\alpha=\dfrac{\cos\alpha}{\sin\alpha}$, $\ \mathbf{u}\cdot\mathbf{v}=\|\mathbf{u}\|\|\mathbf{v}\|\cos\alpha$, $\ \|\mathbf{u}\times\mathbf{v}\|=\|\mathbf{u}\|\|\mathbf{v}\|\sin\alpha$; target $\cot\alpha=\dfrac{\mathbf{u}\cdot\mathbf{v}}{\|\mathbf{u}\times\mathbf{v}\|}$ <br> <i>Hint:</i> divide the dot-product identity by the cross-product identity so $\|\mathbf{u}\|\|\mathbf{v}\|$ cancels. :: A: $\cot\alpha=\dfrac{\cos\alpha}{\sin\alpha}$ <br> divide the two identities so $\|\mathbf{u}\|\|\mathbf{v}\|$ cancels <br> $\dfrac{\mathbf{u}\cdot\mathbf{v}}{\|\mathbf{u}\times\mathbf{v}\|}=\dfrac{\|\mathbf{u}\|\|\mathbf{v}\|\cos\alpha}{\|\mathbf{u}\|\|\mathbf{v}\|\sin\alpha}$ <br> $\boxed{\cot\alpha=\dfrac{\mathbf{u}\cdot\mathbf{v}}{\|\mathbf{u}\times\mathbf{v}\|}}$
- Q: Derive the cotan edge weight $w_{ij}$ for an interior edge $(i,j)$ shared by two triangles (derive, step 2 of 2). <br> <b>Formula:</b> each triangle contributes a half-cotangent $\tfrac12\cot(\cdot)$; target $w_{ij}=\tfrac12(\cot\alpha_{ij}+\cot\beta_{ij})$ <br> <i>Hint:</i> add the half-cotangent from the $\alpha$-side triangle to the half-cotangent from the $\beta$-side triangle. :: A: triangle on the $\alpha$ side gives $\tfrac12\cot\alpha_{ij}$ <br> triangle on the $\beta$ side gives $\tfrac12\cot\beta_{ij}$ <br> sum the two faces sharing the edge <br> $\boxed{w_{ij}=\tfrac12(\cot\alpha_{ij}+\cot\beta_{ij})}$
- Q: Derive the cotangent Laplacian at vertex $i$ in fully expanded form (derive, step 1 of 2). <br> <b>Formula:</b> $(L\mathbf{f})_i=\sum_{j\in N(i)} w_{ij}(f_j-f_i)$ with $w_{ij}=\tfrac12(\cot\alpha_{ij}+\cot\beta_{ij})$ <br> <i>Hint:</i> substitute the edge weight into the sum and pull the $\tfrac12$ out front. :: A: substitute the edge weight into the sum <br> $(L\mathbf{f})_i=\sum_{j\in N(i)}\tfrac12(\cot\alpha_{ij}+\cot\beta_{ij})(f_j-f_i)$ <br> $\boxed{(L\mathbf{f})_i=\tfrac12\sum_{j\in N(i)}(\cot\alpha_{ij}+\cot\beta_{ij})(f_j-f_i)}$
- Q: Derive the pointwise Laplace-Beltrami value $(\Delta f)_i$ from the integrated form (derive, step 2 of 2). <br> <b>Formula:</b> $(L\mathbf{f})_i=\tfrac12\sum_{j}(\cot\alpha_{ij}+\cot\beta_{ij})(f_j-f_i)$ and $\Delta f=M^{-1}L\mathbf{f}$ with $M_{ii}=A_i$ <br> <i>Hint:</i> $M$ is diagonal, so $M^{-1}$ just divides row $i$ by the dual area $A_i$. :: A: $M$ is diagonal so $(\Delta f)_i=\tfrac{1}{A_i}(L\mathbf{f})_i$ <br> divide the integrated sum by the dual area $A_i$ <br> $\boxed{(\Delta f)_i=\dfrac{1}{2A_i}\sum_{j\in N(i)}(\cot\alpha_{ij}+\cot\beta_{ij})(f_j-f_i)}$
- Q: Derive the diagonal entry $L_{ii}$ of the stiffness matrix from the requirement that constant functions vanish. <br> <b>Formula:</b> $(L\mathbf{f})_i=L_{ii}f_i+\sum_{j\in N(i)} L_{ij}f_j$ with off-diagonals $L_{ij}=+w_{ij}$, and the constraint $(L\mathbf{1})_i=0$ <br> <i>Hint:</i> set $f_i=f_j=1$ and solve the resulting equation for $L_{ii}$. :: A: $(L\mathbf{f})_i=L_{ii}f_i+\sum_{j\in N(i)} w_{ij} f_j$ <br> a constant $f\equiv1$ must give $0$: $L_{ii}+\sum_j w_{ij}=0$ <br> $\boxed{L_{ii}=-\sum_{j\in N(i)} w_{ij}}$ (negative diagonal, rows sum to zero)
- Q: Derive the discrete Poisson system from the continuous equation $\Delta u=\rho$. <br> <b>Formula:</b> $(\Delta u)_i=\tfrac{1}{A_i}(L\mathbf{u})_i$, $\ A_i=M_{ii}$; target $L\mathbf{u}=M\boldsymbol\rho$ <br> <i>Hint:</i> set the pointwise Laplacian equal to $\rho_i$, then clear the $\tfrac{1}{A_i}$ by multiplying through by $A_i$. :: A: $\tfrac{1}{A_i}(L\mathbf{u})_i=\rho_i$ <br> multiply both sides by $A_i=M_{ii}$ <br> $(L\mathbf{u})_i=A_i\rho_i=(M\boldsymbol\rho)_i$ <br> $\boxed{L\mathbf{u}=M\boldsymbol\rho}$
- Q: Derive the solvability (compatibility) condition for the Poisson system $L\mathbf{u}=M\boldsymbol\rho$. <br> <b>Formula:</b> $L\mathbf{1}=\mathbf{0}$, $L=L^\top$; target $\sum_i A_i\rho_i=0$ <br> <i>Hint:</i> left-multiply both sides by $\mathbf{1}^\top$ and use $\mathbf{1}^\top L=(L\mathbf{1})^\top=\mathbf{0}$ to kill the left side. :: A: $\mathbf{1}^\top L\mathbf{u}=\mathbf{1}^\top M\boldsymbol\rho$ <br> $L$ symmetric and $L\mathbf{1}=\mathbf{0}$ give $\mathbf{1}^\top L=\mathbf{0}$ <br> so the left side is $0$ <br> $\boxed{\mathbf{1}^\top M\boldsymbol\rho=\sum_i A_i\rho_i=0}$ (source must integrate to zero)
- Q: Derive the implicit (backward-Euler) smoothing step (derive, step 1 of 2). <br> <b>Formula:</b> flow $\dot{\mathbf{V}}=M^{-1}L\mathbf{V}$; backward Euler evaluates the right side at the NEW state: $\dfrac{\mathbf{V}^{k+1}-\mathbf{V}^k}{\delta}=M^{-1}L\mathbf{V}^{k+1}$ <br> <i>Hint:</i> multiply through by $\delta$ and gather both $\mathbf{V}^{k+1}$ terms on the left. :: A: $\dfrac{\mathbf{V}^{k+1}-\mathbf{V}^k}{\delta}=M^{-1}L\mathbf{V}^{k+1}$ <br> multiply through by $\delta$: $\mathbf{V}^{k+1}-\mathbf{V}^k=\delta M^{-1}L\mathbf{V}^{k+1}$ <br> $\boxed{\mathbf{V}^{k+1}-\delta M^{-1}L\mathbf{V}^{k+1}=\mathbf{V}^k}$
- Q: Derive the symmetric SPD smoothing system (derive, step 2 of 2). <br> <b>Formula:</b> $\mathbf{V}^{k+1}-\delta M^{-1}L\mathbf{V}^{k+1}=\mathbf{V}^k$; target $(M-\delta L)\mathbf{V}^{k+1}=M\mathbf{V}^k$ <br> <i>Hint:</i> left-multiply the whole equation by $M$ to clear $M^{-1}$, then factor $\mathbf{V}^{k+1}$. :: A: $M\mathbf{V}^{k+1}-\delta L\mathbf{V}^{k+1}=M\mathbf{V}^k$ <br> factor the unknown on the left <br> $\boxed{(M-\delta L)\,\mathbf{V}^{k+1}=M\mathbf{V}^k}$ ($-L$ is positive, so $M-\delta L$ is SPD for any $\delta$)
- Q: Worked example (small, 2D). For the triangle $P_i=(0,0)$, $P_j=(4,0)$, $P_k=(0,3)$, compute the half-weight from the angle at $P_k$ opposite edge $(i,j)$. <br> <b>Formula:</b> $\mathbf{u}=P_i-P_k$, $\mathbf{v}=P_j-P_k$, $\ \cot\alpha=\dfrac{\mathbf{u}\cdot\mathbf{v}}{\|\mathbf{u}\times\mathbf{v}\|}$, $\ w^{\text{half}}=\tfrac12\cot\alpha$; 2D cross $=|u_xv_y-u_yv_x|$ <br> <i>Hint:</i> first form the two edge vectors from $P_k$, then plug into the dot and cross. :: A: $\mathbf{u}=P_i-P_k=(0,-3)$, $\mathbf{v}=P_j-P_k=(4,-3)$ <br> $\mathbf{u}\cdot\mathbf{v}=0\cdot4+(-3)(-3)=9$ <br> $\|\mathbf{u}\times\mathbf{v}\|=|0\cdot(-3)-(-3)\cdot4|=12$ <br> $\cot\alpha_{ij}=9/12=0.75$ <br> $\boxed{w_{ij}^{\text{half}}=\tfrac12(0.75)=0.375}$ (boundary edge, one face)
- Q: Worked example (two-triangle quad). Mirror the triangle above with $P_{k'}=(0,-3)$ so edge $(i,j)$ is interior with $\cot\beta_{ij}=0.75$ by symmetry. Find the full weight and $(\Delta f)_i$. <br> <b>Formula:</b> $w_{ij}=\tfrac12(\cot\alpha_{ij}+\cot\beta_{ij})$; libigl $L_{ij}=+w_{ij}$, $L_{ii}=-w_{ij}$; $\ (\Delta f)_i=\tfrac{1}{A_i}(L\mathbf{f})_i$ with $A_i=4$ <br> <i>Hint:</i> compute $w_{ij}$ first, then assemble the single-neighbor row $(L\mathbf{f})_i=w_{ij}(f_j-f_i)$ before dividing by $A_i$. :: A: $w_{ij}=\tfrac12(0.75+0.75)=0.75$ <br> libigl: $L_{ij}=+0.75$, $L_{ii}=-0.75$, so $(L\mathbf{f})_i=0.75(f_j-f_i)$ <br> $(\Delta f)_i=\tfrac{1}{A_i}(L\mathbf{f})_i=\tfrac14\cdot0.75(f_j-f_i)$ <br> $\boxed{(\Delta f)_i=0.1875\,(f_j-f_i)}$
- Q: Worked example (unit right triangle, Practice #1). For vertices $(0,0)$, $(1,0)$, $(0,1)$, give the single-triangle stiffness diagonals $L_{00},L_{11},L_{22}$. <br> <b>Formula:</b> half-weight $w=\tfrac12\cot(\text{opposite angle})$; $\cot90^\circ=0$, $\cot45^\circ=1$; diagonal $L_{ii}=-\sum_{j}w_{ij}$ over edges at $i$ <br> <i>Hint:</i> find the three corner angles first, convert each to its half-weight, then each diagonal is minus the sum of its two incident half-weights. :: A: angle at $(0,0)$ is $90^\circ$ so $\cot=0$; the other two are $45^\circ$ so $\cot=1$ <br> half-weights: $w_{01}=\tfrac12(1)=0.5$, $w_{02}=\tfrac12(1)=0.5$, $w_{12}=\tfrac12(0)=0$ <br> $L_{00}=-(w_{01}+w_{02})$, $L_{11}=-(w_{01}+w_{12})$, $L_{22}=-(w_{02}+w_{12})$ <br> $\boxed{L_{00}=-1,\;L_{11}=-0.5,\;L_{22}=-0.5}$ (each row sums to zero)
- Q: Worked example (equilateral, regular regime). On a mesh of equilateral triangles every angle opposite an interior edge is $60^\circ$; find the interior edge weight. <br> <b>Formula:</b> $\cot\theta=\dfrac{\cos\theta}{\sin\theta}$, $\ w_{ij}=\tfrac12(\cot\alpha_{ij}+\cot\beta_{ij})$; $\cos60^\circ=\tfrac12$, $\sin60^\circ=\tfrac{\sqrt3}{2}$ <br> <i>Hint:</i> evaluate $\cot60^\circ$ first, then both opposite angles are equal so just double the half-cotangent. :: A: $\cot 60^\circ=\dfrac{\cos60^\circ}{\sin60^\circ}=\dfrac{1/2}{\sqrt3/2}=\dfrac{1}{\sqrt3}$ <br> both opposite angles equal $60^\circ$ <br> $w_{ij}=\tfrac12(\tfrac{1}{\sqrt3}+\tfrac{1}{\sqrt3})=\tfrac{1}{\sqrt3}$ <br> $\boxed{w_{ij}=\tfrac{1}{\sqrt3}\approx0.577}$ (positive, well-conditioned)
- Q: Worked example (large scale, units). A quarry slab mesh in meters has mean edge length $h\approx0.2$ m; pick the smoothing step and rescale it to millimeters. <br> <b>Formula:</b> $\delta=\tfrac12 h^2$ (choose $\delta$ relative to $h^2$); same physical mesh, $1\text{ m}=1000\text{ mm}$ <br> <i>Hint:</i> compute $\delta$ in $\text{m}^2$ first, then redo with $h=200$ mm to see $\delta$ jump by $1000^2$. :: A: meters: $\delta=\tfrac12(0.2)^2=0.02\ \text{m}^2$ <br> in mm, $h=200$ mm so $\delta=\tfrac12(200)^2=2.0\times10^4\ \text{mm}^2$ <br> $\boxed{\delta\propto h^2\text{; choose it from mean edge length, never mix units}}$
- Q: Test yourself (close the loop). For the unit right triangle stiffness row $L_{00}=-1,L_{01}=0.5,L_{02}=0.5$, verify the row by applying it to the constant $\mathbf{f}=(1,1,1)$. <br> <b>Formula:</b> $(L\mathbf{f})_0=L_{00}f_0+L_{01}f_1+L_{02}f_2$; a correct stiffness row satisfies $(L\mathbf{1})_0=0$ <br> <i>Hint:</i> a pass means the three entries sum to zero. :: A: $(L\mathbf{f})_0=L_{00}\cdot1+L_{01}\cdot1+L_{02}\cdot1$ <br> $=-1+0.5+0.5$ <br> $=0$ <br> $\boxed{(L\mathbf{1})_0=0\ \checkmark}$ constants are in the null space, as required
- Q: Test yourself (limiting case + units). Verify the pointwise Laplacian $\Delta f=M^{-1}L\mathbf{f}$ has units of $f$ per length$^2$. <br> <b>Formula:</b> cotangents are dimensionless so $(L\mathbf{f})_i\sim\text{area}\cdot[f]\sim[f]\cdot\text{length}^2$; $M_{ii}=A_i\sim\text{length}^2$ <br> <i>Hint:</i> a pass means dividing $(L\mathbf{f})$ by $M$ leaves exactly $[f]/\text{length}^2$, matching the smooth $\nabla^2 f$. :: A: $L$ entries (cotangents) are dimensionless, but $(L\mathbf{f})_i\approx \text{area}\cdot\Delta f$ carries one factor of area <br> $M_{ii}=A_i\sim\text{length}^2$ <br> $M^{-1}L\mathbf{f}\sim\dfrac{[f]\cdot\text{length}^2}{\text{length}^2\cdot\text{length}^2}=\dfrac{[f]}{\text{length}^2}$ <br> $\boxed{[\Delta f]=[f]/\text{length}^2\ \checkmark}$ matching the smooth $\nabla^2 f$
- Q: Test yourself (sphere, analytic check). On a sphere of radius $r=2$ m, check the discrete mean curvature against the closed form. <br> <b>Formula:</b> $\mathbf{H}\mathbf{n}=-M^{-1}L\mathbf{V}$, row norm $\|(M^{-1}L\mathbf{V})_i\|=2|H_i|$; analytic sphere $|H|=1/r$ <br> <i>Hint:</i> a pass means the row norm clusters at $2/r$, i.e. $|H_i|\to1/r$, as the mesh refines. :: A: analytic mean curvature of a sphere is $|H|=1/r$ <br> $|H|=1/2=0.5\ \text{m}^{-1}$ <br> the row norm $\|(M^{-1}L\mathbf{V})_i\|=2|H_i|$ should concentrate at $1.0\ \text{m}^{-1}$ <br> $\boxed{|H_i|\to0.5\ \text{m}^{-1}\ \checkmark}$ as the mesh refines
- Q: When should you reach for the cotan operator in the scan->CAD->IFC pipeline? <br> <i>Hint:</i> think of the four jobs that share one matrix. :: A: To smooth/fair a noisy scan before NURBS fitting, to compute the mean-curvature field for feature-aware segmentation, to solve harmonic/Poisson UV parameterizations, and as the stiffness/mass for the Signed Heat SDF.
- Q: Why does the Poisson right-hand side carry the mass matrix $M$ (recall)? <br> <b>Formula:</b> $L\mathbf{u}=M\boldsymbol\rho$ <br> <i>Hint:</i> compare the units of the two sides — $L\mathbf{u}$ is already integrated over each dual cell. :: A: $M$ integrates the pointwise source $\rho$ over each dual cell, matching the integrated left side $L\mathbf{u}$; without it the units and scaling do not balance.
- Q: Pitfall: what is the single most common cotan-Laplacian bug, and the one-line check for it? <br> <b>Formula:</b> libigl uses $-M^{-1}L\mathbf{V}$ and $M-\delta L$; check $L\mathbf{1}=\mathbf{0}$ <br> <i>Hint:</i> a pass means a constant function maps to the zero vector. :: A: A flipped sign convention. libigl $L$ is negative semidefinite, so use $-M^{-1}L\mathbf{V}$ and $M-\delta L$; check that $L\mathbf{1}=\mathbf{0}$ or the surface inflates and explodes.
- Q: Pitfall: what goes wrong on a very obtuse triangle, and the fix? <br> <b>Formula:</b> $w_{ij}=\tfrac12(\cot\alpha_{ij}+\cot\beta_{ij})$ <br> <i>Hint:</i> ask what happens to $\cot\theta$ as $\theta\to180^\circ$. :: A: Its cotangent is large and negative, so an edge weight can go negative and smoothing loses the maximum principle; use the Voronoi/mixed (obtuse-safe) mass area or an intrinsic Delaunay flip.
- CLOZE: The cotan Laplacian splits into a {{c1::stiffness}} matrix $L$ and a {{c2::mass}} matrix $M$, with the pointwise operator equal to $M^{-1}L$.
- CLOZE: In the libigl convention the cotmatrix is {{c1::negative semidefinite}}, so mean curvature carries a leading minus sign in $-M^{-1}L\mathbf{V}$.
- CLOZE: The Poisson system $L\mathbf{u}=M\boldsymbol\rho$ is singular because constants lie in the {{c1::null space}} of $L$, so you pin a vertex or add a mean-zero constraint. <br> <i>hint:</i> $L\mathbf{1}=\mathbf{0}$.
- CLOZE: The cotangent of an angle from its two edge vectors is $\cot\theta={{c1::\mathbf{u}\cdot\mathbf{v}}}/{{c2::\|\mathbf{u}\times\mathbf{v}\|}}$.
- CLOZE: Implicit fairing is stable for any time step because backward Euler solves the SPD system $({{c1::M-\delta L}})\mathbf{V}'=M\mathbf{V}$.
- CLOZE: To assemble the operators in libigl you call {{c1::igl::cotmatrix}} for $L$ and {{c2::igl::massmatrix}} for $M$, then evaluate $-M^{-1}LV$.
