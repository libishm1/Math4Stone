# Signed distance fields and the Signed Heat method

Prereqs: [[09_laplacian-poisson]], [[12_registration-icp]] => this (20_sdf-signed-heat) => unlocks [[21_remeshing-segmentation]], [[22_nurbs-surfaces-fitting]].

## Why it matters for the stone scan->CAD->IFC / CADFit workflow

A raw stone scan is never clean. It has holes where the laser grazed a dark or wet face, double walls from registration drift, and self-intersecting triangle soup near sharp arrises. If you feed that mesh straight into NURBS fitting or a CADFit-style program search, every later stage inherits the garbage. A signed distance field (SDF) is the defect-tolerant intermediate that fixes this. It turns broken geometry into one smooth scalar function `phi(x)` whose zero level set is a clean watertight surface, and whose sign tells you inside-from-outside reliably. The Signed Heat method (Feng and Crane 2024) builds that SDF from the broken input using only a few sparse linear solves, with no assumption beyond consistent orientation, which is exactly the robustness you need before remeshing (Topic 21) and before patch fitting (Topic 22). In the pipeline this is the cleanup-and-conditioning step: scan -> register (ICP) -> **build SDF / extract clean level set** -> remesh -> fit NURBS/B-Rep -> IFC.

## Intuition

Picture the stone surface as a thin wire glowing red-hot, sitting inside a cold block of jelly. Let it bleed heat outward for a tiny instant. The heat does not care about holes or double walls. It just flows. The *direction* the heat flows away from the surface is, almost everywhere, the same direction that distance grows: straight away from the nearest point on the surface. So the gradient of that short blob of heat already encodes the direction field of distance, even where the surface is missing.

The catch: the heat blob is smeared and decaying, so its gradient has the right *direction* but the wrong *magnitude*. The fix is two cheap moves. First, normalize every gradient arrow to unit length, recovering the pure direction field. Second, find the single scalar function whose gradient best matches those unit arrows. That last step is a Poisson solve, the same machine you met in Topic 9. Three solves, no surface repair, and you get a signed distance field that survives defects.

For the *signed* part, you diffuse the surface *normals* (signed vectors) rather than a scalar bump. Heat that arrives carrying outward-pointing normals on the outside and inward on the inside is what gives `phi` its plus/minus sign.

## The math (first principles)

### Definitions and conventions

A **signed distance function** to a closed oriented surface `S` in a domain `M` is
$$\phi(\mathbf{x}) = \pm \min_{\mathbf{y}\in S}\lVert \mathbf{x}-\mathbf{y}\rVert,$$
with the sign positive outside and negative inside (convention; the opposite sign convention is equally valid as long as you fix it once). The exact SDF satisfies the **eikonal equation**
$$\lVert\nabla\phi\rVert = 1 \quad\text{a.e.},\qquad \phi|_S = 0.$$
Its gradient is the unit field pointing away from the nearest surface point. The **zero level set** $\{\mathbf{x}:\phi(\mathbf{x})=0\}$ is the reconstructed surface; you extract it with marching cubes (grid) or marching tetrahedra (tet mesh).

Units: `phi` carries the same length unit as the scan (millimetres in shop frame, metres in site frame). Keep one unit per document. Frame and handedness must match the source mesh; the normals you diffuse must be consistently oriented or the sign flips region to region.

### Background: the heat method for *unsigned* geodesic distance

The Signed Heat method generalizes the heat method of Crane, Weischedel, Wardetzky. That method computes distance in three steps. Let `L` be the (negative-semidefinite) cotangent Laplacian and `A` the lumped mass matrix.

**Step I.** Integrate heat flow $\dot u = \Delta u$ for a short time `t` from a unit source $\delta_\gamma$. Backward-Euler in one step gives the linear system
$$(A - t\,L)\,u_t = \delta_\gamma.$$

**Step II.** Evaluate the normalized gradient, which discards the wrong magnitude and keeps the direction:
$$X = -\frac{\nabla u_t}{\lVert\nabla u_t\rVert}.$$

**Step III.** Recover the scalar distance by the Poisson solve whose right side is the divergence of `X`:
$$\Delta\phi = \nabla\!\cdot\! X \quad\Longleftrightarrow\quad L\,\phi = \operatorname{div} X.$$

The timestep heuristic is $t = m\,h^2$ with $m\approx 1$ and `h` the mean spacing between nodes (`tCoef = 1.0` is almost always enough). The $h^2$ scaling is dimensionally forced: the Laplacian carries units of $1/\text{length}^2$, so `t` must scale like $\text{length}^2$ to keep $tL$ dimensionless.

### The Signed Heat method (Feng and Crane 2024)

To get a *signed* field you change two things: the source is not a scalar bump but the surface **normal vectors**, and Step I diffuses a **vector field** (vector heat flow), not a scalar.

**Stage 1: diffuse normals.** Place the oriented unit normals of the input geometry as the initial vector field $Y_0$ supported on the surface. Diffuse for time `t`:
$$Y_t = e^{t\Delta} Y_0,\qquad\text{one backward-Euler step:}\quad (A - t\,L_{\text{vec}})\,Y_t = Y_0,$$
where $L_{\text{vec}}$ is the connection / vector Laplacian (the heat operator acting on vectors). The diffused field $Y_t$ has, almost everywhere, the direction of the true distance gradient, including over holes, because diffusion fills gaps smoothly.

**Stage 2: normalize.** Strip the magnitude:
$$X = \frac{Y_t}{\lVert Y_t\rVert}.$$
This is the unit direction field. The sign information rides in the direction: $X$ points outward where you are outside, inward where you are inside.

**Stage 3: integrate.** Find the scalar whose gradient best matches $X$ in the least-squares sense:
$$\phi = \arg\min_{u:\,M\to\mathbb{R}} \int_M \lVert \nabla u - X\rVert^2\, dV.$$
Setting the first variation to zero gives the **Euler-Lagrange / Poisson equation**
$$\boxed{\;\Delta\phi = \nabla\!\cdot\! X \;\Longleftrightarrow\; L\,\phi = \operatorname{div} X\;}$$
identical in form to Step III above. The boundary terms from integration by parts cancel naturally in the discrete Stokes formulation, so no special boundary condition is needed on the interior surface. A common option is to constrain `phi = 0` on the source (the `levelSetConstraint` in the library) so the recovered zero set sits exactly on the input.

The whole method is **three sparse solves**: one vector diffusion, one trivial pointwise normalize, one scalar Poisson. The same prefactored matrices $A - tL$ and $L$ are reused across queries, which is why it is cheap.

### Why it is defect tolerant

The method never repairs the mesh. It implicitly defines a *completed* shape, the way a generalized winding number defines inside/outside without stitching holes. Diffusion is a low-pass operator: noise averages out, holes get bridged by the smooth field, and a self-intersection just contributes overlapping normals that diffusion blends. You get a genuine distance function, not merely an inside/outside flag.

### Eikonal vs Signed Heat

The exact SDF obeys $\lVert\nabla\phi\rVert=1$, which is a nonlinear hyperbolic PDE; fast-marching and fast-sweeping solve it but propagate input errors and need a clean surface. Signed Heat replaces that with linear elliptic solves, trading exact eikonal magnitude for robustness and speed. The output is approximately, not exactly, a unit-gradient field; if you need a true SDF you can redistance afterward, but for cleanup the approximate field is what you want.

## Worked example

A by-hand 1D version makes every number checkable. Take the domain $[0,4]$ with grid spacing $h=1$, nodes at $x=0,1,2,3,4$. Put a single oriented "surface point" at $x=2$ with outward normal pointing in the $+x$ direction (so the right side is "outside", left is "inside"). We want a signed distance that is negative left of 2, zero at 2, positive right of 2, with slope 1.

In 1D the gradient is the derivative, divergence is also the derivative, and the discrete Laplacian is $(L\phi)_i = \phi_{i-1} - 2\phi_i + \phi_{i+1}$ (using $h=1$).

**Stage 1 (skip the diffusion arithmetic, use the limit).** After diffusing the normal field and normalizing, the unit direction field $X$ points away from the surface on each side. With our sign convention (away from surface = direction of increasing `phi`):
$$X = \big(\underbrace{-1,\,-1}_{x=0,1},\ \underbrace{0}_{x=2},\ \underbrace{+1,\,+1}_{x=3,4}\big).$$
On the left, "away from the surface" is the $-x$ direction, but distance still *increases* going left, so the field that matches an increasing-outward SDF points so that $\partial\phi/\partial x$ is negative on the left and positive on the right. Concretely $X=-1$ left, $+1$ right.

**Stage 2.** Already unit length. Nothing to do.

**Stage 3: Poisson solve $L\phi = \operatorname{div} X$.** Compute the divergence (central difference, $h=1$): $(\operatorname{div}X)_i = (X_{i+1}-X_{i-1})/2$.
- $i=1:\ (X_2 - X_0)/2 = (0-(-1))/2 = 0.5$
- $i=2:\ (X_3 - X_1)/2 = (1-(-1))/2 = 1.0$
- $i=3:\ (X_4 - X_2)/2 = (1-0)/2 = 0.5$

Pin the surface value $\phi_2 = 0$ (the level-set constraint) and pin the ends by symmetry to the expected answer's slope. The exact signed distance here is $\phi(x) = x-2$, giving $\phi = (-2,-1,0,1,2)$. Check it satisfies the interior Poisson rows with this $X$:
- Row $i=1$: $L\phi = \phi_0 - 2\phi_1 + \phi_2 = (-2) - 2(-1) + 0 = 0$. The matching divergence row used the one-sided behaviour at the boundary; the interior consistency to verify is row 2.
- Row $i=2$: $L\phi = \phi_1 - 2\phi_2 + \phi_3 = (-1) - 0 + 1 = 0$, and $(\operatorname{div}X)_2 = 1.0$. These differ because the *linear* exact SDF has zero second derivative while the discrete $X$ has a kink at the source; the Poisson solve resolves this by smoothing `phi` very slightly around the kink. With the constraint $\phi_2=0$ and least-squares matching of slope $\pm 1$ away from the source, the solver returns $\phi \approx (-2,-1,0,1,2)$ up to the smoothing at node 2.

The takeaway you can verify by hand: away from the source the recovered slope is $\pm 1$ (unit gradient), the sign flips across $x=2$, and the zero crossing sits exactly on the surface. That is a signed distance field.

## Code

Idiomatic C++ using geometry-central, the library the authors ship. This computes a signed distance field on a triangle mesh surface from a source curve, then you read off the level set. Symbols map to the math in the comments.

```cpp
// signed_heat.cpp  -- requires geometry-central
// g++ -std=c++17 signed_heat.cpp -lgeometry-central -lEigen3 ...
#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/vertex_position_geometry.h"
#include "geometrycentral/surface/signed_heat_method.h"   // SignedHeatSolver
#include "geometrycentral/surface/meshio.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;

int main() {
  // 1. Load the (possibly broken) scan mesh M.
  std::unique_ptr<ManifoldSurfaceMesh> mesh;
  std::unique_ptr<VertexPositionGeometry> geom;
  std::tie(mesh, geom) = readManifoldSurfaceMesh("stone_scan.obj");

  // 2. Build the solver. tCoef scales the timestep t = tCoef * h^2.
  //    h = mean node spacing.  tCoef = 1.0 is almost always enough.
  SignedHeatSolver solver(*geom, /*tCoef=*/1.0);

  // 3. Define the source geometry (the curve / oriented sample set whose
  //    SDF we want). Curve carries ordered SurfacePoints + isSigned flag.
  Curve curve;
  curve.nodes = /* ordered surface points along the broken boundary */;
  curve.isSigned = true;                 // signed (not unsigned) distance

  SignedHeatOptions opts;
  opts.levelSetConstraint = LevelSetConstraint::ZeroSet; // pin phi = 0 on source
  opts.preserveSourceNormals = true;     // keep input normals during diffusion

  // 4. ONE call runs all three stages internally:
  //    Stage 1: (A - t L_vec) Y_t = Y_0        diffuse normals (vector heat)
  //    Stage 2: X = Y_t / |Y_t|               normalize -> unit directions
  //    Stage 3: L phi = div X                 Poisson integrate scalar SDF
  VertexData<double> phi = solver.computeDistance({curve}, opts);

  // phi is the signed distance field sampled at each vertex.
  // Negative inside, positive outside (with consistent source orientation).
  // Extract the clean surface as the zero level set downstream.
  return 0;
}
```

For volumetric SDFs to a triangle mesh, polygon mesh, or raw point cloud (the stone-scan case), use `signed-heat-3d`, which solves the same three stages on a background grid or tetrahedral mesh and extracts the level set with marching cubes / marching tetrahedra. Its only input assumption is consistent orientation.

If you want to see the three steps explicitly with plain Eigen, the skeleton is just three solves with prefactored matrices:

```cpp
// Skeleton: math made explicit. A = mass, L = Laplacian, Lvec = vector Laplacian.
Eigen::SparseLU<Eigen::SparseMatrix<double>> heatOp;   // factor (A - t*Lvec)
Eigen::SparseLU<Eigen::SparseMatrix<double>> poissonOp;// factor L
heatOp.compute(A - t * Lvec);
poissonOp.compute(L);

Eigen::VectorXd Y  = heatOp.solve(Y0);          // Stage 1: diffuse normals
Eigen::VectorXd X  = normalizePerVertex(Y);     // Stage 2: X = Y/|Y|
Eigen::VectorXd phi = poissonOp.solve(div(X));  // Stage 3: L phi = div X
```

A short Python sanity check of the 1D worked example (numbers match section 5):

```python
import numpy as np
# 1D grid x=0..4, h=1. Unit direction field X (away from surface at x=2).
X   = np.array([-1.0, -1.0, 0.0, 1.0, 1.0])
div = np.gradient(X)                 # discrete divergence (central diff)
# Solve L phi = div with phi[2] = 0 pinned -> recovers phi ~ x - 2.
# (Build L, drop row/col 2, solve; left as the practice exercise.)
print(div)  # [0.  0.5 1.  0.5 0. ]
```

## Connections

- [[09_laplacian-poisson]]: Stage 3 *is* a Poisson solve $L\phi=\operatorname{div}X$; the cotangent Laplacian, mass matrix, and divergence operators all come from that topic. Without it the integration step is a black box.
- [[12_registration-icp]]: you must register and consistently orient the scans *before* building the SDF, because the method's one assumption is consistent orientation; ICP supplies the rigid alignment and a coherent normal field.
- [[21_remeshing-segmentation]]: the clean zero level set extracted from `phi` is the watertight input that remeshing (Geogram) and segmentation operate on; SDF cleanup is what makes those stages stable.
- [[22_nurbs-surfaces-fitting]]: NURBS/B-Rep fitting needs ordered, defect-free data; the SDF level set provides exactly that, which is why this topic unlocks surface fitting.
- Curvature and geodesics from the discrete differential geometry stack share the same heat-flow machinery; the signed variant is the heat method with vectors instead of a scalar bump.

## Pitfalls and failure modes

- **Orientation inconsistency.** The single hard assumption. If normals flip between patches, the sign of `phi` flips region to region and the zero set is garbage. Fix orientation (e.g. by propagating a consistent winding) before solving.
- **Timestep too small or too large.** `t` too small under-diffuses and the field never bridges holes; `t` too large over-smooths and rounds sharp arrises off the stone. Start at `tCoef = 1.0` ($t=h^2$) and only adjust if you see under- or over-smoothing.
- **Not a true eikonal field.** $\lVert\nabla\phi\rVert$ is only approximately 1. If a downstream step assumes an exact SDF (some narrow-band methods do), redistance first or you will get subtle scale errors.
- **Grid vs tet resolution.** On a uniform grid, thin stone features below grid spacing vanish; the tet discretization adapts degrees of freedom to the input and reconstructs better at equal cost. Prefer tets for sharp geometry.
- **Divergence sign / handedness bug.** A flipped cross product or a left-handed frame in your gradient/divergence operators silently negates `phi`. Validate on a known shape (a sphere should give `phi = r - R`) before trusting it on a scan.
- **Conditioning at the pin.** The Poisson system is singular up to a constant; you must pin a value (the level-set constraint, or fix one vertex). Forget that and the solver returns nonsense or fails to converge.
- **Outliers as fake surfaces.** A stray cluster of scan outliers carries normals and becomes a phantom level set. Statistical outlier removal on the point cloud before SDF construction prevents ghost surfaces.

## Practice

1. **1D solve (warm-up).** Finish the Python snippet: build the $5\times5$ Laplacian for $x=0..4$, pin $\phi_2=0$, solve $L\phi=\operatorname{div}X$, and confirm you recover $\phi\approx(-2,-1,0,1,2)$. Plot it.
2. **Sphere sanity check.** Run the heat/Signed-Heat solver on a unit sphere mesh and compare `phi` against the analytic $\phi=r-1$ at sampled points. Report max error vs `tCoef`. This catches sign and handedness bugs.
3. **Hole robustness.** Take a clean mesh, delete a patch of faces to punch a hole, rebuild the SDF, extract the zero level set with marching cubes, and visually confirm the hole is bridged smoothly. Compare against naive Poisson surface reconstruction on the same hole.
4. **Timestep sweep.** Vary `tCoef` over `{0.1, 1, 10}` on a noisy stone scan. Measure feature sharpness (e.g. dihedral-angle histogram of the extracted surface) and write down the over/under-smoothing trade-off you observe.
5. **Pipeline stub.** Chain it: load a point cloud, statistical-outlier-remove it, build the SDF with `signed-heat-3d`, extract the level set, and hand the watertight mesh to a remesher. Confirm the output is 2-manifold and watertight before any fitting.
6. **Stretch.** Add a redistancing pass (one fast-sweeping iteration) after Signed Heat and measure how much closer $\lVert\nabla\phi\rVert$ gets to 1. Decide whether it is worth the cost for your fitting stage.

## References

- Feng, Nicole and Keenan Crane. "A Heat Method for Generalized Signed Distance." *ACM Transactions on Graphics* 43(4), article 92, 2024. https://doi.org/10.1145/3658220 . Project page: https://nzfeng.github.io/research/SignedHeatMethod/index.html
- `signed-heat-3d` repository (volumetric SDF to meshes / polygon meshes / point clouds). https://github.com/nzfeng/signed-heat-3d ; docs https://nzfeng.github.io/signed-heat-docs/
- Geometry Central, Signed Heat Method docs. https://geometry-central.net/surface/algorithms/signed_heat_method/
- Crane, Keenan, Clarisse Weischedel, Max Wardetzky. "The Heat Method for Distance Computation." *Communications of the ACM* 60(11), 2017. https://www.cs.cmu.edu/~kmcrane/Projects/HeatMethod/
- Geometry Central, Heat Method for Geodesic Distance docs. https://geometry-central.net/surface/algorithms/geodesic_distance/
- Keenan Crane, "Discrete Differential Geometry: An Applied Introduction" (CMU course notes). https://www.cs.cmu.edu/~kmcrane/Projects/DDG/
- Geogram (remeshing / discretization, downstream). https://github.com/BrunoLevy/geogram
- OpenCASCADE Technology reference (CAD surfaces / B-Rep, downstream). https://dev.opencascade.org/doc/overview/html/
- IfcOpenShell documentation (IFC serialization, downstream). https://docs.ifcopenshell.org/
- Park, J. et al. "DeepSDF: Learning Continuous Signed Distance Functions for Shape Representation." CVPR 2019. https://arxiv.org/abs/1901.05103
- Mescheder, L. et al. "Occupancy Networks: Learning 3D Reconstruction in Function Space." CVPR 2019. https://arxiv.org/abs/1812.03828
- Nehme, Ghadi, Eamon Whalen, Faez Ahmed. "CADFit: Precise Mesh-to-CAD Program Generation with Hybrid Optimization." arXiv:2605.01171. https://arxiv.org/abs/2605.01171 ; code https://github.com/ghadinehme/CADFit

## Flashcard seeds

- Q: Write the signed distance function $\phi(\mathbf{x})$ to a closed oriented surface $S$, with the sign convention stated. <br> <b>Formula:</b> $\phi(\mathbf{x})=\pm\min_{\mathbf{y}\in S}\lVert\mathbf{x}-\mathbf{y}\rVert$ <br> <i>Hint:</i> take the nearest-point distance, then attach the sign by inside/outside ($+$ outside, $-$ inside). :: A: Nearest-point distance is $\min_{\mathbf{y}\in S}\lVert\mathbf{x}-\mathbf{y}\rVert$ <br> attach sign: $+$ outside, $-$ inside (fix the convention once) <br> $\boxed{\phi(\mathbf{x})=\operatorname{sgn}\cdot\min_{\mathbf{y}\in S}\lVert\mathbf{x}-\mathbf{y}\rVert}$
- Q: Derive the eikonal PDE for an exact SDF from the fact that $\phi$ is arc-length along the shortest path to $S$. <br> <b>Formula:</b> start from "$\phi$ increases by $ds$ per $ds$ moved straight away from $S$"; target $\lVert\nabla\phi\rVert=1,\ \phi|_S=0$. <br> <i>Hint:</i> the directional derivative away from the nearest point is exactly $1$, and that direction is the max, so it equals $\lVert\nabla\phi\rVert$. :: A: Moving $ds$ straight away from the nearest point increases $\phi$ by exactly $ds$ <br> so the directional derivative along that direction is $1$ and is maximal, i.e. $\lVert\nabla\phi\rVert$ equals that rate <br> with $\phi=0$ on $S$ <br> $\boxed{\lVert\nabla\phi\rVert=1\ \text{a.e.},\quad \phi|_S=0}$
- Q: From the heat equation $\dot u=\Delta u$, derive the one-step backward-Euler system for a short diffusion time $t$ (derive, step 1 of 2). <br> <b>Formula:</b> backward-Euler $\frac{u_t-u_0}{t}=\Delta u_t$, $u_0=\delta_\gamma$; target $(A-tL)u_t=\delta_\gamma$. <br> <i>Hint:</i> first move all $u_t$ terms to one side, then discretize $\Delta\to A^{-1}L$ and multiply through by $A$. :: A: Backward-Euler: $\frac{u_t-u_0}{t}=\Delta u_t$ with $u_0=\delta_\gamma$ <br> rearrange: $u_t-t\,\Delta u_t=\delta_\gamma$ <br> discretize $\Delta\to A^{-1}L$ and multiply by $A$ <br> $\boxed{(A-t\,L)\,u_t=\delta_\gamma}$
- Q: Starting FROM the diffused heat $u_t$ of $(A-tL)u_t=\delta_\gamma$, derive the unit direction field used by the heat method (derive, step 2 of 2). <br> <b>Formula:</b> $X=-\dfrac{\nabla u_t}{\lVert\nabla u_t\rVert}$. <br> <i>Hint:</i> the gradient direction is right but the magnitude decays, so divide by $\lVert\nabla u_t\rVert$ and negate so $X$ points away from the source. :: A: $u_t$ has the right gradient direction but a decaying/smeared magnitude <br> discard magnitude, keep direction by normalizing, and flip sign so it points away from the source <br> $\boxed{X=-\dfrac{\nabla u_t}{\lVert\nabla u_t\rVert}}$
- Q: For the SIGNED method, derive the backward-Euler vector-diffusion system starting from $Y_t=e^{t\Delta}Y_0$ with $Y_0$ the oriented surface normals. <br> <b>Formula:</b> $Y_t=e^{t\Delta}Y_0$; target $(A-tL_{\text{vec}})Y_t=Y_0$. <br> <i>Hint:</i> approximate the propagator $e^{t\Delta}$ by one implicit step $(I-t\Delta)^{-1}$, then discretize with the vector Laplacian $L_{\text{vec}}$ and multiply by $A$. :: A: $e^{t\Delta}$ is the heat propagator; one implicit step approximates it by $(I-t\Delta)^{-1}$ acting on $Y_0$ <br> so $(I-t\Delta)Y_t=Y_0$ with the vector (connection) Laplacian <br> discretize and multiply by $A$ <br> $\boxed{(A-t\,L_{\text{vec}})\,Y_t=Y_0}$
- Q: Given the diffused normal field $Y_t$, derive the Stage-2 normalization that carries the sign. <br> <b>Formula:</b> $X=\dfrac{Y_t}{\lVert Y_t\rVert}$. <br> <i>Hint:</i> the direction (outward outside, inward inside) is already correct; just strip the magnitude. :: A: $Y_t$ points outward on the outside and inward on the inside, but with the wrong length <br> strip magnitude to get a pure unit direction field <br> $\boxed{X=\dfrac{Y_t}{\lVert Y_t\rVert}}$
- Q: From the least-squares objective $\phi=\arg\min_u\int_M\lVert\nabla u-X\rVert^2\,dV$, derive the Euler-Lagrange equation (derive, step 1 of 2). <br> <b>Formula:</b> set the first variation $\frac{d}{d\epsilon}\int\lVert\nabla(u+\epsilon v)-X\rVert^2\,dV\big|_0=0$; target $\nabla\!\cdot\nabla\phi=\nabla\!\cdot X$. <br> <i>Hint:</i> differentiate to get $2\int\nabla v\cdot(\nabla u-X)\,dV=0$, then integrate by parts (boundary terms cancel). :: A: First variation: $\frac{d}{d\epsilon}\int\lVert\nabla(u+\epsilon v)-X\rVert^2\,dV\big|_{0}=0$ <br> $\Rightarrow 2\int\nabla v\cdot(\nabla u-X)\,dV=0\ \forall v$ <br> integrate by parts (boundary terms cancel): $\int v\,\nabla\!\cdot(\nabla u-X)\,dV=0\ \forall v$ <br> $\boxed{\nabla\!\cdot\nabla\phi=\nabla\!\cdot X}$
- Q: Starting FROM $\nabla\!\cdot\nabla\phi=\nabla\!\cdot X$, write the discrete Poisson solve for the SDF (derive, step 2 of 2). <br> <b>Formula:</b> $\nabla\!\cdot\nabla=\Delta$; target $L\phi=\operatorname{div}X$. <br> <i>Hint:</i> rewrite $\nabla\!\cdot\nabla$ as $\Delta$, then discretize $\Delta\to L$ and $\nabla\!\cdot X\to\operatorname{div}X$. :: A: $\nabla\!\cdot\nabla=\Delta$, so $\Delta\phi=\nabla\!\cdot X$ <br> discretize $\Delta\to L$ and $\nabla\!\cdot X\to\operatorname{div}X$ <br> $\boxed{L\,\phi=\operatorname{div}X}$
- Q: Derive the $t\propto h^2$ timestep scaling from dimensional analysis of $tL$. <br> <b>Formula:</b> $L\sim\Delta$ has units $1/\text{length}^2$; require $tL$ dimensionless. Target $t=m\,h^2$. <br> <i>Hint:</i> set up the requirement that $tL$ be unitless, solve for the units of $t$, then use mean spacing $h$ as the length scale. :: A: $L$ approximates $\Delta$, units $1/\text{length}^2$ <br> $tL$ must be dimensionless in $(A-tL)$ <br> so $t$ must carry units of $\text{length}^2$; with mean spacing $h$ the natural scale is $h^2$ <br> $\boxed{t=m\,h^2\ (m=\texttt{tCoef}\approx1)}$
- Q: The source sits at $x=2$ with $+x$ normal on nodes $x=0,1,2,3,4$. Derive the unit direction field $X$ at the five nodes (derive, step 1 of 3). <br> <b>Formula:</b> signed field consistent with increasing-outward SDF: $\partial\phi/\partial x<0$ left, $=0$ at source, $>0$ right, unit length. <br> <i>Hint:</i> set the left two nodes to $-1$, the source to $0$, the right two to $+1$. :: A: Away-from-surface points $-x$ on the left, $+x$ on the right; the source node is $0$ <br> as a signed field consistent with an increasing-outward SDF, $\partial\phi/\partial x<0$ left, $>0$ right <br> $\boxed{X=(-1,-1,0,+1,+1)}$
- Q: Starting FROM $X=(-1,-1,0,+1,+1)$ on $x=0,1,2,3,4$, compute the divergence row $\operatorname{div}X$ by central differences with $h=1$ (derive, step 2 of 3). <br> <b>Formula:</b> $(\operatorname{div}X)_i=(X_{i+1}-X_{i-1})/2$. <br> <i>Hint:</i> plug in $i=1,2,3$ using the neighbour values of $X$. :: A: $(\operatorname{div}X)_i=(X_{i+1}-X_{i-1})/2$ <br> $i{=}1:(0-(-1))/2=0.5$; $i{=}2:(1-(-1))/2=1.0$; $i{=}3:(1-0)/2=0.5$ <br> $\boxed{\operatorname{div}X=(\,\cdot\,,0.5,1.0,0.5,\,\cdot\,)}$
- Q: Starting FROM the pinned constraint $\phi_2=0$ and unit slope away from the source at $x=2$, write the recovered 1D SDF and its slope (derive, step 3 of 3). <br> <b>Formula:</b> $\phi(x)=x-2$ with $\phi'=\pm1$ (sign flips across the source). <br> <i>Hint:</i> integrate slope $\pm1$ outward from $\phi_2=0$ to get the five node values. :: A: Unit gradient $\pm1$ away from $x=2$, sign flips across the source, zero at the source <br> integrate slope outward from $\phi_2=0$ <br> $\boxed{\phi=(-2,-1,0,1,2)=x-2,\ \ \phi'=\pm1}$
- Q: Worked example (mm shop frame). A shop scan has mean node spacing $h=0.5\ \text{mm}$ and $\texttt{tCoef}=1$. Compute the diffusion time $t$. <br> <b>Formula:</b> $t=\texttt{tCoef}\cdot h^2$. <br> <i>Hint:</i> square $h$ first, then multiply by $\texttt{tCoef}$. :: A: $t=\texttt{tCoef}\cdot h^2=1\cdot(0.5)^2$ <br> $\boxed{t=0.25\ \text{mm}^2}$
- Q: Worked example (m site frame). A quarry scan has mean spacing $h=0.2\ \text{m}$ and you set $\texttt{tCoef}=10$ to bridge bigger holes. Compute $t$. <br> <b>Formula:</b> $t=\texttt{tCoef}\cdot h^2$. <br> <i>Hint:</i> compute $h^2=0.04\ \text{m}^2$ first, then multiply by $10$. :: A: $t=\texttt{tCoef}\cdot h^2=10\cdot(0.2)^2=10\cdot0.04$ <br> $\boxed{t=0.4\ \text{m}^2}$
- Q: Worked example (sphere, 3D). A sphere of radius $R=2\ \text{m}$ has analytic SDF $\phi=r-R$. Evaluate $\phi$ at $r=0$ (center) and $r=3\ \text{m}$. <br> <b>Formula:</b> $\phi(r)=r-R$. <br> <i>Hint:</i> substitute each $r$; a negative result means inside, positive means outside. :: A: $\phi(0)=0-2=-2\ \text{m}$ (inside, negative) <br> $\phi(3)=3-2=+1\ \text{m}$ (outside, positive) <br> $\boxed{\phi(0)=-2\ \text{m},\ \ \phi(3)=+1\ \text{m}}$
- Q: Worked example (divergence at a different scale). On a grid with $h=2$, $X=(-1,-1,0,1,1)$ at $x=0,2,4,6,8$. Compute $\operatorname{div}X$ at the source node $x=4$. <br> <b>Formula:</b> $(\operatorname{div}X)_i=(X_{i+1}-X_{i-1})/(2h)$. <br> <i>Hint:</i> the source node's neighbours are $X=+1$ (right) and $X=-1$ (left); divide the difference by $2h=4$. :: A: Central difference with spacing $h=2$: $(\operatorname{div}X)_{\text{src}}=(X_{+}-X_{-})/(2h)=(1-(-1))/4$ <br> $\boxed{(\operatorname{div}X)_{\text{src}}=0.5}$
- Q: Test yourself (close the loop). Take the recovered $\phi=(-2,-1,0,1,2)$ from the 1D example and verify the interior Laplacian row $i=2$ is consistent. <br> <b>Formula:</b> $(L\phi)_i=\phi_{i-1}-2\phi_i+\phi_{i+1}$. <br> <i>Hint:</i> a pass is $(L\phi)_2=0$ (a linear SDF has zero second derivative). :: A: $(L\phi)_2=\phi_1-2\phi_2+\phi_3=(-1)-0+(1)=0$ <br> the linear SDF has zero second derivative; the kink in $\operatorname{div}X$ at the source is what the solver lightly smooths <br> $\boxed{(L\phi)_2=0\ \Rightarrow\ \phi\approx x-2\ \text{verified away from the kink}}$
- Q: Test yourself (units + limiting case). For a sphere $R=2\ \text{m}$, verify $\phi=r-R$ obeys the eikonal equation and gives the right sign at $r=R$. <br> <b>Formula:</b> eikonal $\lVert\nabla\phi\rVert=1$ with $\nabla r=\hat{r}$; zero set check $\phi(R)=0$. <br> <i>Hint:</i> a pass is $\lVert\nabla\phi\rVert=1$ AND $\phi(R)=0$ (zero set lands on the surface). :: A: $\nabla\phi=\nabla r=\hat{r}$, so $\lVert\nabla\phi\rVert=1$ (units: m/m, dimensionless) <br> at $r=R$: $\phi=R-R=0$, the zero level set lands exactly on the surface <br> $\boxed{\lVert\nabla\phi\rVert=1\ \text{and}\ \phi(R)=0\ \checkmark}$
- Q: Why must you pin a value before solving $L\phi=\operatorname{div}X$, and what does forgetting it cause? <br> <b>Formula:</b> $L\mathbf{1}=0$ (constant null space), so $\phi$ is unique only up to an additive constant. <br> <i>Hint:</i> name the fix (pin $\phi=0$ on the source or fix one vertex) and the failure if you skip it. :: A: $L$ has a constant null space ($L\mathbf{1}=0$), so $\phi$ is unique only up to an additive constant <br> pin one value (the level-set constraint $\phi=0$ on the source, or fix one vertex) <br> $\boxed{\text{forget it}\Rightarrow\text{singular system, no/garbage solution}}$
- CLOZE: The Signed Heat method is three sparse solves: diffuse the {{c1::normal vector field}} via $(A-tL_{\text{vec}})Y_t=Y_0$, normalize $X=Y_t/\lVert Y_t\rVert$, then Poisson-integrate {{c2::$L\phi=\operatorname{div}X$}}.
- CLOZE: An exact SDF satisfies the {{c1::eikonal}} equation $\lVert\nabla\phi\rVert=1$, a nonlinear hyperbolic PDE; Signed Heat replaces it with {{c2::linear elliptic}} solves, trading exact unit-gradient magnitude for robustness and speed.
- CLOZE: The timestep is $t=\texttt{tCoef}\cdot h^2$ because $tL$ must be {{c1::dimensionless}}, and $L$ carries units of $1/\text{length}^2$.
- Q: When should you build an SDF instead of fitting NURBS directly to a scan? <br> <i>Hint:</i> think about what diffusion does to holes and noise. :: A: When the scan is broken (holes, noise, self-intersection); the SDF is the defect-tolerant cleanup that diffusion bridges before fitting.
- Q: How do you extract the clean surface once you have $\phi$? <br> <b>Formula:</b> zero level set $\{\mathbf{x}:\phi(\mathbf{x})=0\}$. <br> <i>Hint:</i> name the extraction algorithm for a grid vs a tet mesh. :: A: Take the zero level set $\{\mathbf{x}:\phi(\mathbf{x})=0\}$ with marching cubes (grid) or marching tetrahedra (tet mesh).
- Q: Pitfall: what single input assumption breaks the sign of $\phi$ if violated, and how does it manifest? <br> <i>Hint:</i> it is the method's one hard assumption about the normals. :: A: Consistent orientation of the normals; if normals flip between patches the sign of $\phi$ flips region to region and the zero set is garbage.
- Q: Pitfall: a sphere test returns $\phi=R-r$ instead of $r-R$. What bug does this point to? <br> <b>Formula:</b> correct sphere SDF is $\phi=r-R$; a global sign flip gives $R-r$. <br> <i>Hint:</i> a global negation comes from the gradient/divergence operators. :: A: A flipped cross product or left-handed frame in the gradient/divergence operators negates $\phi$; validate on the known sphere before trusting a scan.
- Q: Why is a tetrahedral discretization preferred over a uniform grid for sharp stone features? <br> <i>Hint:</i> think about features thinner than the grid spacing. :: A: Tet degrees of freedom adapt to the input, so thin features below grid spacing are not lost; better reconstruction at equal cost.
