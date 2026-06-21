# NURBS surfaces and surface fitting from scans

Prereqs: [[22_bsplines-nurbs-curves]], [[06_linear-least-squares]], [[21_remeshing-segmentation]] => this => unlocks [[24_brep-boolean-occt]].

## Why it matters for the stone scan->CAD->IFC / CADFit workflow

A raw stone scan is millions of unordered 3D points. A CAD kernel like OpenCASCADE cannot consume that. It needs a smooth parametric surface defined over a clean $(u,v)$ domain, so a NURBS surface is the bridge from "cloud of points" to "trimmable B-Rep face". In the pipeline this is the step *fit a smooth field or spline*: after you align the scan, remesh it, and segment it into chart-sized patches, you resample each patch onto an ordered grid and fit a B-spline surface that downstream booleans, offsets, and IFC `IfcAdvancedBrep` faces can use. The hard, load-bearing detail is that `GeomAPI_PointsToBSplineSurface` demands an *ordered* $(u,v)$ grid, so the real engineering work is template-guided UV resampling of an unorganized cloud, not the fit call itself. Get the parameterization wrong and the surface ripples between samples, which then poisons every CAD operation built on top of it.

## Intuition

A NURBS surface is a rubber sheet pinned by a grid of control points. You do not touch the sheet directly. You move the pins, and the sheet follows smoothly, staying inside the convex hull of the pins. The grid structure is the key word *tensor product*: pick a curve basis for the $u$ direction and a separate curve basis for the $v$ direction, then the surface basis at $(u,v)$ is just the product of the two. So a surface is a curve of curves. Slide along $v$ and each row of control points sweeps out a $u$-curve; the family of those curves is the surface.

Analogy: think of a topographic map of a hillside. The control points are survey stakes on a regular grid. The NURBS surface is the smoothest tent you can pull over those stakes. If your stakes are evenly spaced you get a clean tent. If three stakes are crammed together and the next is far away, a naive tent overshoots and oscillates near the crowded stakes. Centripetal parameterization is the trick that spaces the *parameter* values to match the *physical* spacing of stakes, which tames that overshoot.

## The math (first principles)

### Conventions

Right-handed coordinates, units in meters at site scale and millimeters at shop scale (be consistent within one fit). Parameter domain is $(u,v) \in [u_0, u_{\max}] \times [v_0, v_{\max}]$, typically normalized to $[0,1]^2$. Degrees are $p$ in $u$ and $q$ in $v$; cubic ($p=q=3$) is the default for $C^2$ continuity. Tolerance `Tol3D` is the max allowed 3D distance from a sample to the fitted surface, in model units.

### B-spline basis (recap)

The degree-$p$ B-spline basis $N_{i,p}(u)$ is built from a knot vector $U = \{u_0,\dots,u_m\}$ by the Cox-de Boor recursion. The defining count relation is

$$m = n + p + 1,$$

where $n+1$ is the number of control points in $u$ and $m+1$ is the number of knots. The basis is nonnegative, has local support over $p+1$ knot spans, and forms a partition of unity, $\sum_{i=0}^{n} N_{i,p}(u) = 1$ on the domain. Same in $v$ with basis $M_{j,q}(v)$, degree $q$, knot vector $V$.

### Tensor-product NURBS surface

A NURBS surface with control net $P_{ij}$ (a grid of $(n+1)\times(m+1)$ points) and weights $w_{ij} > 0$ is

$$S(u,v) = \frac{\displaystyle\sum_{i=0}^{n}\sum_{j=0}^{m} w_{ij}\, P_{ij}\, N_{i,p}(u)\, M_{j,q}(v)}{\displaystyle\sum_{i=0}^{n}\sum_{j=0}^{m} w_{ij}\, N_{i,p}(u)\, M_{j,q}(v)}.$$

Define the rational basis function

$$R_{ij}(u,v) = \frac{w_{ij}\, N_{i,p}(u)\, M_{j,q}(v)}{\sum_{k}\sum_{l} w_{kl}\, N_{k,p}(u)\, M_{l,q}(v)},$$

so that $S(u,v) = \sum_{i}\sum_{j} R_{ij}(u,v)\, P_{ij}$ and $\sum_{i}\sum_{j} R_{ij}(u,v) = 1$. The rational basis is also a partition of unity, which gives the convex-hull and affine-invariance properties: transform the control points and the surface transforms identically.

If all weights are equal, the denominator collapses to $\sum N \sum M = 1$ and the surface reduces to a non-rational tensor-product **B-spline surface**:

$$S(u,v) = \sum_{i=0}^{n}\sum_{j=0}^{m} P_{ij}\, N_{i,p}(u)\, M_{j,q}(v).$$

For scan fitting you almost always fit a non-rational B-spline surface (all weights 1). Rational weights matter for exact conics and quadrics, which raw stone does not have. OCCT's `GeomAPI_PointsToBSplineSurface` produces a non-rational result.

### Interpolation vs approximation

Given a grid of data points $Q_{rs}$, $r=0\dots R$, $s=0\dots S$:

**Interpolation** finds control points $P_{ij}$ so the surface passes through every $Q_{rs}$ exactly: $S(\bar u_r, \bar v_s) = Q_{rs}$. The number of control points equals the number of data points. This is a square linear system. With scan data this is the wrong choice: it forces the surface to honor every noisy sample, so noise becomes geometry and the surface ripples.

**Approximation** chooses fewer control points than data points and minimizes squared distance:

$$\min_{\{P_{ij}\}} \sum_{r}\sum_{s} \big\| S(\bar u_r, \bar v_s) - Q_{rs} \big\|^2 .$$

Because $S$ is linear in the control points (for fixed knots and parameters), this is an ordinary linear least-squares problem. Stack the basis evaluations into a matrix $B$ where row $(r,s)$ holds $N_{i,p}(\bar u_r) M_{j,q}(\bar v_s)$ across all $(i,j)$. Then $P$ solves the normal equations

$$B^\top B\, P = B^\top Q.$$

This is exactly [[06_linear-least-squares]] applied per coordinate ($x$, $y$, $z$). Fewer control points act as a low-pass filter, which is why approximation is the default for scans. The fit error is bounded by `Tol3D`; OCCT raises the control-point count until that tolerance is met.

### Parameterization: the choice that makes or breaks the fit

Before you can evaluate any basis function you must assign a parameter value $\bar u_r$ to each row and $\bar v_s$ to each column of data. Three standard schemes, given a sequence of points $Q_0,\dots,Q_R$ along one direction with total length $d = \sum_{k=1}^{R} \lvert Q_k - Q_{k-1}\rvert$:

**Uniform (equally spaced):**

$$\bar u_0 = 0,\quad \bar u_R = 1,\quad \bar u_r = \frac{r}{R}.$$

Ignores geometry. Overshoots badly when sample spacing is uneven.

**Chord-length:**

$$\bar u_r = \bar u_{r-1} + \frac{\lvert Q_r - Q_{r-1}\rvert}{d}.$$

Parameter advances in proportion to physical distance. Good when spacing is moderate, but big jumps in spacing still cause loops and overshoot.

**Centripetal (Lee 1989):**

$$\bar u_r = \bar u_{r-1} + \frac{\sqrt{\lvert Q_r - Q_{r-1}\rvert}}{\sum_{k=1}^{R}\sqrt{\lvert Q_k - Q_{k-1}\rvert}}.$$

Uses the square root of the distance, equivalent to raising each chord to the power $\alpha = \tfrac12$. The square root compresses the dynamic range of the spacings, so a sample that is far away does not grab a disproportionate slice of parameter space. This is why centripetal almost always gives the best shape on irregular sample spacing, which is exactly the regime of stone scans. OCCT exposes this as `Approx_Centripetal`.

The pattern: $\alpha = 0$ recovers uniform, $\alpha = 1$ recovers chord-length, $\alpha = \tfrac12$ is centripetal. For a *surface*, you compute parameters per row and per column and then average (OCCT handles this internally given an ordered array).

### The ordered-grid requirement

`GeomAPI_PointsToBSplineSurface` takes a `TColgp_Array2OfPnt`, a 2D array indexed $(i,j)$. The fit assumes $Q_{ij}$ is topologically a grid: row $i$ varies smoothly in $u$, column $j$ varies smoothly in $v$, with no holes and no folds. A scan cloud has none of that structure. The resampling step (see Code) is what manufactures the grid: lay a template parameter grid over the patch, and for each grid node compute a robust surface point from nearby cloud samples (PCA-plane projection or weighted median along the template normal). Cells with too few samples become holes that you inpaint or split off into a separate patch.

## Worked example

Take a tiny non-rational B-spline **surface** with bidegree $(p,q)=(1,1)$, that is a bilinear patch, control net $2\times 2$, all weights 1. Domain $[0,1]^2$. Degree-1 basis over knot vector $\{0,0,1,1\}$ gives the two hat functions

$$N_{0,1}(u) = 1-u,\qquad N_{1,1}(u) = u,$$

and identically $M_{0,1}(v)=1-v$, $M_{1,1}(v)=v$. Check partition of unity: $N_0+N_1 = (1-u)+u = 1$. Good.

Control points (heights $z$ over a unit square, $x=i$, $y=j$):

$$P_{00}=(0,0,0),\ P_{10}=(1,0,2),\ P_{01}=(0,1,1),\ P_{11}=(1,1,4).$$

The surface is

$$S(u,v) = (1-u)(1-v)P_{00} + u(1-v)P_{10} + (1-u)v\,P_{01} + uv\,P_{11}.$$

Evaluate the center $(u,v)=(0.5,0.5)$. Each basis product is $0.5\cdot 0.5 = 0.25$, so

$$S(0.5,0.5) = 0.25\,(P_{00}+P_{10}+P_{01}+P_{11}) = 0.25\,(0+1+0+1,\ 0+0+1+1,\ 0+2+1+4).$$

That is $0.25\,(2,2,7) = (0.5,\ 0.5,\ 1.75)$. The center height is $1.75$, the average of the four corner heights $\{0,2,1,4\}$. Sanity check: a bilinear patch interpolates its corners and the center is the corner average. Correct.

Now the parameterization point. Suppose along one row the three sample points sit at chord distances $1, 1, 8$ (one big gap). Total chord $= 10$.

Chord-length parameters: $0,\ 0.1,\ 0.2,\ 1.0$. The first three samples are squeezed into $[0,0.2]$; the basis must bend hard there, inviting overshoot.

Centripetal: square roots are $\sqrt1, \sqrt1, \sqrt8 = 1, 1, 2.83$, total $4.83$. Parameters: $0,\ 0.207,\ 0.414,\ 1.0$. The crowded samples now occupy $[0,0.414]$, a far gentler distribution. Same points, calmer parameterization, less ripple. This is the centripetal advantage in one arithmetic step.

## Code

### C++ with OpenCASCADE (the fit) plus Eigen (the resample)

```cpp
// nurbs_fit.cpp
// Build: link OpenCASCADE (TKMath TKG3d TKGeomBase TKBRep TKTopAlgo) and Eigen.
//
// Pipeline: unorganized cloud  ->  template-guided ordered (u,v) grid
//           ->  GeomAPI_PointsToBSplineSurface  ->  Geom_BSplineSurface.
//
// Math-to-code map:
//   S(u,v) = sum_ij P_ij N_ip(u) M_jq(v)   <-- the fitted Geom_BSplineSurface
//   Tol3D                                  <-- max ||S(u,v) - Q|| at samples
//   Approx_Centripetal                     <-- bar_u_r via sqrt(|dQ|) (Lee 1989)
//   [DegMin, DegMax]                       <-- allowed degree range p,q
//   GeomAbs_C2                             <-- continuity target (cubic, C2)

#include <Eigen/Dense>
#include <vector>
#include <limits>
#include <stdexcept>

#include <gp_Pnt.hxx>
#include <TColgp_Array2OfPnt.hxx>
#include <GeomAPI_PointsToBSplineSurface.hxx>
#include <Geom_BSplineSurface.hxx>
#include <Approx_ParametrizationType.hxx>
#include <GeomAbs_Shape.hxx>

using Vec3 = Eigen::Vector3d;

// Step 1: resample an unorganized cloud onto an ordered nu x nv grid.
// Template = best-fit plane (PCA). Project all points; bin into a regular
// (u,v) lattice; each cell's height z is the mean of the points that land in it.
// Cells with no support are flagged so the caller can inpaint or split a patch.
TColgp_Array2OfPnt resampleToGrid(const std::vector<Vec3>& cloud,
                                  int nu, int nv,
                                  std::vector<std::vector<bool>>& filled)
{
    if (cloud.size() < 4) throw std::runtime_error("need >=4 points");

    // PCA frame: centroid + covariance eigenvectors (see eig-svd-pca topic).
    Vec3 c = Vec3::Zero();
    for (const auto& p : cloud) c += p;
    c /= double(cloud.size());

    Eigen::Matrix3d cov = Eigen::Matrix3d::Zero();
    for (const auto& p : cloud) { Vec3 d = p - c; cov += d * d.transpose(); }
    cov /= double(cloud.size() - 1);

    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> es(cov);
    // eigenvalues ascending: col(2),col(1) span the plane; col(0) ~ normal.
    Vec3 ax = es.eigenvectors().col(2).normalized(); // u axis (max spread)
    Vec3 ay = es.eigenvectors().col(1).normalized(); // v axis
    Vec3 az = ax.cross(ay).normalized();             // template normal
    ay = az.cross(ax).normalized();                  // re-orthogonalize

    // Project to (a,b) plane coords; find bounds.
    double amin =  1e300, amax = -1e300, bmin = 1e300, bmax = -1e300;
    std::vector<double> aa(cloud.size()), bb(cloud.size()), hh(cloud.size());
    for (size_t k = 0; k < cloud.size(); ++k) {
        Vec3 d = cloud[k] - c;
        aa[k] = d.dot(ax); bb[k] = d.dot(ay); hh[k] = d.dot(az);
        amin = std::min(amin, aa[k]); amax = std::max(amax, aa[k]);
        bmin = std::min(bmin, bb[k]); bmax = std::max(bmax, bb[k]);
    }

    // Accumulate heights per cell.
    std::vector<double> zsum(nu * nv, 0.0);
    std::vector<int>    zcnt(nu * nv, 0);
    auto cell = [&](int i, int j) { return i * nv + j; };
    for (size_t k = 0; k < cloud.size(); ++k) {
        int i = int((aa[k] - amin) / (amax - amin) * (nu - 1) + 0.5);
        int j = int((bb[k] - bmin) / (bmax - bmin) * (nv - 1) + 0.5);
        i = std::max(0, std::min(nu - 1, i));
        j = std::max(0, std::min(nv - 1, j));
        zsum[cell(i, j)] += hh[k];
        zcnt[cell(i, j)] += 1;
    }

    // Build the OCC array. Note: TColgp_Array2OfPnt is 1-based.
    TColgp_Array2OfPnt grid(1, nu, 1, nv);
    filled.assign(nu, std::vector<bool>(nv, false));
    for (int i = 0; i < nu; ++i)
        for (int j = 0; j < nv; ++j) {
            double a = amin + (amax - amin) * i / (nu - 1);
            double b = bmin + (bmax - bmin) * j / (nv - 1);
            double z = (zcnt[cell(i, j)] > 0)
                       ? zsum[cell(i, j)] / zcnt[cell(i, j)] : 0.0;
            filled[i][j] = (zcnt[cell(i, j)] > 0);
            Vec3 world = c + a * ax + b * ay + z * az;
            grid.SetValue(i + 1, j + 1, gp_Pnt(world.x(), world.y(), world.z()));
        }
    return grid; // caller must inpaint holes before fitting
}

// Step 2: fit a C2 cubic B-spline surface with centripetal parameterization.
Handle(Geom_BSplineSurface) fitSurface(const TColgp_Array2OfPnt& grid,
                                       double tol3D = 1.0e-3)
{
    GeomAPI_PointsToBSplineSurface fitter;
    fitter.Init(grid,
                Approx_Centripetal,  // ParType: sqrt-distance spacing (Lee 1989)
                3,                   // DegMin
                8,                   // DegMax (OCC default upper bound)
                GeomAbs_C2,          // continuity target
                tol3D);              // Tol3D: max surface-to-point distance
    if (!fitter.IsDone())
        throw std::runtime_error("B-spline surface fit failed");
    return fitter.Surface();         // Geom_BSplineSurface, non-rational
}
```

Two facts the comments depend on, both verified against the OCCT reference: the constructor and `Init` defaults are `DegMin=3, DegMax=8, Continuity=GeomAbs_C2, Tol3D=1.0e-3`, and `Approx_Centripetal` advances each parameter by $\sqrt{\text{distance}}$ which "can get better result for irregular distances between points."

### Python sanity check (de Boor, no CAD kernel)

A short non-rational surface evaluator clarifies the tensor-product formula and lets you unit-test a fit before trusting OCCT.

```python
import numpy as np
from scipy.interpolate import BSpline

def surface_point(u, v, ctrl, ku, kv, p, q):
    # ctrl: (nu, nv, 3). Tensor product: evaluate v-curves, then a u-curve.
    nu, nv, _ = ctrl.shape
    Nu = np.array([BSpline.basis_element(ku[i:i+p+2])(u)
                   if ku[i] <= u <= ku[i+p+1] else 0.0 for i in range(nu)])
    Nv = np.array([BSpline.basis_element(kv[j:j+q+2])(v)
                   if kv[j] <= v <= kv[j+q+1] else 0.0 for j in range(nv)])
    # S(u,v) = sum_i sum_j N_i(u) M_j(v) P_ij  (all weights = 1)
    return np.einsum('i,j,ijc->c', Nu, Nv, ctrl)
```

## Connections

- [[22_bsplines-nurbs-curves]]: the surface is a tensor product of two curve bases; every basis property (locality, partition of unity, Cox-de Boor, knot relation $m=n+p+1$) carries over unchanged.
- [[06_linear-least-squares]]: surface *approximation* is least squares per coordinate with the basis matrix $B$; the normal equations $B^\top B P = B^\top Q$ are the engine inside `GeomAPI_PointsToBSplineSurface`.
- [[21_remeshing-segmentation]]: you must remesh and segment the scan into single-valued, chart-sized patches before fitting, because one tensor-product surface cannot represent an overhang or a fold.
- [[03_eig-svd-pca]]: PCA gives the local template plane and $(u,v)$ axes used to resample the cloud into an ordered grid.
- [[24_brep-boolean-occt]]: the fitted `Geom_BSplineSurface` becomes a B-Rep face (via `BRepBuilderAPI_MakeFace`) that boolean operations and IFC `IfcAdvancedBrep` consume; bad surface quality here breaks every boolean downstream.

## Pitfalls and failure modes

- **Rippling / oscillation.** The signature failure. Causes: interpolating noisy data instead of approximating; uniform parameterization on irregular spacing; too many control points; a single patch forced across geometry that bends back on itself. Fixes: approximate not interpolate, use `Approx_Centripetal`, raise `Tol3D` so OCC uses fewer control points, split into multiple patches.
- **Unordered input.** Passing a cloud that is not a true grid produces garbage. The array index $(i,j)$ must be a real topological grid: row $i$ monotone in $u$, column $j$ monotone in $v$. Resample first.
- **Holes in the grid.** Empty cells default to plane height ($z=0$ in the code above), which punches a dent. Inpaint small holes by local averaging before fitting; for large holes split the patch.
- **Folds and overhangs.** A tensor-product surface is single-valued: each $(u,v)$ maps to exactly one point. Stone with overhangs violates this. Segment so each patch is a height field over its template plane.
- **Mixed units.** Fitting a grid in millimeters with `Tol3D = 1.0e-3` (interpreted as $0.001$ mm) over-refines and explodes the control count. Keep grid units and `Tol3D` in the same system and set `Tol3D` to a physically meaningful fraction of feature size.
- **Frame inconsistency.** If the resample plane (PCA frame) and the world frame are confused, the surface fits in the wrong orientation and the later B-Rep lands rotated. Carry one explicit frame and reproject consistently.
- **Outliers.** A few stray scan points drag the mean height of their cells. Use a robust per-cell estimator (median or PCA-plane projection), not the raw mean, and pre-filter with [[09_robust-ransac]] style rejection.
- **Degree too high.** Large `DegMax` lets the fitter chase noise. Cubic ($C^2$) is enough for fabrication; do not raise it without a reason.

## Practice

1. Implement the bilinear patch from the worked example in code. Evaluate $S$ on a $5\times5$ grid of $(u,v)$, print the heights, and confirm the center is $1.75$ and the corners match $P_{ij}$.
2. Generate 1000 points by sampling $z = \sin(2\pi x)\cos(2\pi y)$ on $[0,1]^2$ and adding Gaussian noise $\sigma=0.02$. Resample to a $20\times20$ grid, fit with OCCT once interpolating (`Tol3D` near 0) and once approximating (`Tol3D=0.05`). Mesh both and compare ripple visually.
3. Take one row of 8 points with deliberately uneven spacing. Compute uniform, chord-length, and centripetal parameters by hand and tabulate them. Confirm centripetal compresses the crowded region least aggressively.
4. Write a hole-inpainting pass for the resampled grid (iterative neighbor averaging on unfilled cells). Feed a grid with a punched-out 3x3 hole through it and through the fit; verify the surface no longer dents.
5. Fit a real stone-patch cloud (or a synthetic boulder cap) and report: control-net size, max fit error vs `Tol3D`, and whether `IsDone()` succeeded. Sweep `Tol3D` over $\{0.001, 0.01, 0.1\}$ and plot control count vs tolerance.
6. Take a fitted `Geom_BSplineSurface`, build a face with `BRepBuilderAPI_MakeFace`, and run `BRepCheck_Analyzer`. Diagnose any invalidity; this is the handoff to [[24_brep-boolean-occt]].

## References

- Les Piegl and Wayne Tiller, *The NURBS Book*, 2nd ed., Springer, 1997. Tensor-product surfaces, parameterization, interpolation and approximation (chapters 9-10). https://link.springer.com/book/10.1007/978-3-642-59223-2
- Les Piegl, "On NURBS: a survey," *IEEE Computer Graphics and Applications*, 11(1), 1991. Why rational B-splines unify CAD geometry. https://ieeexplore.ieee.org/document/75531
- E. T. Y. Lee, "Choosing nodes in parametric curve interpolation," *Computer-Aided Design*, 21(6), 1989, pp. 363-370. Origin of centripetal parameterization. https://doi.org/10.1016/0010-4485(89)90003-1
- OpenCASCADE Reference Manual, `GeomAPI_PointsToBSplineSurface` class. Constructor/Init signatures, degree, continuity, `Tol3D`, `Approx_Centripetal`. https://dev.opencascade.org/doc/refman/html/class_geom_a_p_i___points_to_b_spline_surface.html
- OpenCASCADE User Guide, *Modeling Data*. Parameterization types (chord length, centripetal, isoparametric) and approximation criteria. https://dev.opencascade.org/doc/overview/html/occt_user_guides__modeling_data.html
- C. de Boor, *A Practical Guide to Splines*, revised ed., Springer, 2001. Cox-de Boor recursion and least-squares spline fitting.
- IfcOpenShell documentation, geometry processing. Consuming OCCT B-Reps for IFC `IfcAdvancedBrep` faces. https://docs.ifcopenshell.org/
- Keenan Crane, *Discrete Differential Geometry* course notes, CMU. Surface and parameterization background. https://www.cs.cmu.edu/~kmcrane/Projects/DDG/
- ghadinehme/CADFit, CAD program synthesis. The downstream synthesis layer that consumes fitted surfaces. https://github.com/ghadinehme/CADFit (arXiv:2605.01171)

## Flashcard seeds

- Q: Construct the non-rational tensor-product B-spline surface from a $u$-curve basis $N_{i,p}(u)$ and a $v$-curve basis $M_{j,q}(v)$ over a control net $P_{ij}$. <br> <b>Formula:</b> surface basis $=N_{i,p}(u)\,M_{j,q}(v)$; target $S(u,v)=\sum_{i=0}^{n}\sum_{j=0}^{m}P_{ij}\,N_{i,p}(u)\,M_{j,q}(v)$ <br> <i>Hint:</i> "tensor product" means multiply the two bases, then weight each product by its control point $P_{ij}$ and sum over the net. :: A: Tensor product: surface basis at $(u,v)$ is the product $N_{i,p}(u)M_{j,q}(v)$ <br> weight each by control point $P_{ij}$ and sum over all $(i,j)$ <br> $\boxed{S(u,v)=\sum_{i=0}^{n}\sum_{j=0}^{m} P_{ij}\,N_{i,p}(u)\,M_{j,q}(v)}$.
- Q: Derive the rational NURBS surface (step 1 of 2): write the homogeneous numerator and the normalizing denominator from weights $w_{ij}>0$. <br> <b>Formula:</b> $\text{num}=\sum_i\sum_j w_{ij}P_{ij}N_{i,p}M_{j,q}$, $\text{den}=\sum_i\sum_j w_{ij}N_{i,p}M_{j,q}$ <br> <i>Hint:</i> the denominator is the numerator with $P_{ij}$ deleted -- it just renormalizes the weighted basis. :: A: Numerator weights each term by $w_{ij}$: $\sum_i\sum_j w_{ij}P_{ij}N_{i,p}M_{j,q}$ <br> denominator is the same sum without $P_{ij}$: $\sum_i\sum_j w_{ij}N_{i,p}M_{j,q}$ <br> $\boxed{\text{num}=\sum_i\sum_j w_{ij}P_{ij}N_{i,p}M_{j,q},\quad \text{den}=\sum_i\sum_j w_{ij}N_{i,p}M_{j,q}}$.
- Q: Derive the rational NURBS surface (step 2 of 2): form $S=\text{num}/\text{den}$, then show all $w_{ij}=1$ recovers the B-spline surface. <br> <b>Formula:</b> $S=\dfrac{\sum_i\sum_j w_{ij}P_{ij}N_{i,p}M_{j,q}}{\sum_i\sum_j w_{ij}N_{i,p}M_{j,q}}$; partition of unity $\sum_i N_{i,p}=1$, $\sum_j M_{j,q}=1$ <br> <i>Hint:</i> set every $w_{ij}=1$ and use the two partition-of-unity sums to collapse the denominator to $1\cdot1$. :: A: $S(u,v)=\dfrac{\sum_i\sum_j w_{ij}P_{ij}N_{i,p}M_{j,q}}{\sum_i\sum_j w_{ij}N_{i,p}M_{j,q}}$ <br> set $w_{ij}=1$: denominator $=\sum_iN_{i,p}\sum_jM_{j,q}=1\cdot1=1$ <br> $\boxed{S(u,v)=\sum_i\sum_j P_{ij}N_{i,p}(u)M_{j,q}(v)}$.
- Q: Derive the rational basis $R_{ij}(u,v)$ so that $S=\sum_{ij}R_{ij}P_{ij}$, and show it sums to 1. <br> <b>Formula:</b> $R_{ij}=\dfrac{w_{ij}N_{i,p}M_{j,q}}{\sum_k\sum_l w_{kl}N_{k,p}M_{l,q}}$ <br> <i>Hint:</i> factor $P_{ij}$ out of the rational $S$; to sum $R_{ij}$, notice summing the numerator over $ij$ gives exactly the denominator. :: A: Pull $P_{ij}$ out of the rational form: $R_{ij}=\dfrac{w_{ij}N_{i,p}M_{j,q}}{\sum_k\sum_l w_{kl}N_{k,p}M_{l,q}}$ <br> summing over $ij$ makes numerator $=$ denominator <br> $\boxed{\sum_i\sum_j R_{ij}(u,v)=1}$ (partition of unity $\Rightarrow$ affine invariance + convex hull).
- Q: Surface approximation (chain, step 1 of 3): write what you minimize given grid data $Q_{rs}$ and surface samples $S(\bar u_r,\bar v_s)$. <br> <b>Formula:</b> $\min_{\{P_{ij}\}}\sum_r\sum_s\big\|S(\bar u_r,\bar v_s)-Q_{rs}\big\|^2$ <br> <i>Hint:</i> the misfit at each node is the residual $S-Q$; penalize its squared 3D length and sum over the grid. :: A: For each data node the misfit is the residual $S(\bar u_r,\bar v_s)-Q_{rs}$ <br> penalize the squared 3D length, summed over the whole grid <br> $\boxed{\min_{\{P_{ij}\}}\sum_r\sum_s\big\|S(\bar u_r,\bar v_s)-Q_{rs}\big\|^2}$.
- Q: Surface fit (chain, step 2 of 3): from $S=\sum_{ij}P_{ij}N_{i,p}(\bar u_r)M_{j,q}(\bar v_s)$, write the linear model $S=BP$ for fixed knots and parameters. <br> <b>Formula:</b> $S=BP$ with $B_{(r,s),(i,j)}=N_{i,p}(\bar u_r)M_{j,q}(\bar v_s)$ <br> <i>Hint:</i> at fixed knots/params $S$ is linear in $P_{ij}$; build $B$ by putting one row per node $(r,s)$, one column per $(i,j)$. :: A: $S$ is linear in the unknown $P_{ij}$ at fixed knots/params <br> stack one row per node $(r,s)$ with entries $N_{i,p}(\bar u_r)M_{j,q}(\bar v_s)$ across all $(i,j)$ <br> $\boxed{S=BP,\ \ B_{(r,s),(i,j)}=N_{i,p}(\bar u_r)M_{j,q}(\bar v_s)}$.
- Q: Surface fit (chain, step 3 of 3): from $\min_P\|BP-Q\|^2$, derive the normal equations. <br> <b>Formula:</b> $\nabla_P\|BP-Q\|^2=2B^\top(BP-Q)$; target $B^\top B\,P=B^\top Q$ <br> <i>Hint:</i> set the gradient to zero and divide by 2, then move the $B^\top Q$ term to the right. :: A: Set the gradient to zero: $\nabla_P\|BP-Q\|^2=2B^\top(BP-Q)=0$ <br> rearrange <br> $\boxed{B^\top B\,P=B^\top Q}$ (solve once per coordinate $x,y,z$).
- Q: Derive the centripetal parameter increment (step 1 of 2): from the rule "advance by $\sqrt{\text{chord}}$, normalized to the total," write $\Delta\bar u_r$. <br> <b>Formula:</b> $\Delta\bar u_r=\dfrac{\sqrt{|Q_r-Q_{r-1}|}}{\sum_{k=1}^{R}\sqrt{|Q_k-Q_{k-1}|}}$ <br> <i>Hint:</i> the numerator is the square root of this one chord; the denominator is the sum of all square-root chords (so $\bar u$ ends at 1). :: A: Each chord $|Q_r-Q_{r-1}|$ contributes $\sqrt{|Q_r-Q_{r-1}|}$ <br> normalize by the total of all square-root chords so $\bar u$ ends at 1 <br> $\boxed{\Delta\bar u_r=\dfrac{\sqrt{|Q_r-Q_{r-1}|}}{\sum_{k=1}^{R}\sqrt{|Q_k-Q_{k-1}|}}}$.
- Q: Centripetal (step 2 of 2): from $\bar u_r=\bar u_{r-1}+\Delta\bar u_r$ with $\bar u_0=0$, write the running parameter as a partial sum. <br> <b>Formula:</b> $\bar u_r=\dfrac{\sum_{k=1}^{r}\sqrt{|Q_k-Q_{k-1}|}}{\sum_{k=1}^{R}\sqrt{|Q_k-Q_{k-1}|}}$ <br> <i>Hint:</i> unroll the recursion -- the numerator accumulates square-root chords up to index $r$, the denominator stays the full total. :: A: $\bar u_r=\bar u_{r-1}+\Delta\bar u_r$ with $\bar u_0=0$ <br> unrolling the recursion gives the partial sum <br> $\boxed{\bar u_r=\dfrac{\sum_{k=1}^{r}\sqrt{|Q_k-Q_{k-1}|}}{\sum_{k=1}^{R}\sqrt{|Q_k-Q_{k-1}|}}}$.
- Q: From the general rule $\Delta\bar u_r\propto|Q_r-Q_{r-1}|^{\alpha}$, derive uniform, chord-length, and centripetal as the three special cases. <br> <b>Formula:</b> $\alpha=0$ (uniform), $\alpha=\tfrac12$ (centripetal), $\alpha=1$ (chord-length) <br> <i>Hint:</i> substitute each $\alpha$ into the weight $|Q_r-Q_{r-1}|^{\alpha}$ -- start with $\alpha=0$, which makes every chord weigh $1$. :: A: $\alpha=0$: every chord weighted $1\Rightarrow\bar u_r=r/R$ (uniform) <br> $\alpha=1$: weight $=$ raw distance (chord-length) <br> $\alpha=\tfrac12$: weight $=\sqrt{\text{distance}}$ <br> $\boxed{\alpha\in\{0,\tfrac12,1\}\to\{\text{uniform, centripetal, chord-length}\}}$.
- Q: Worked example (bilinear, step 1 of 2): bidegree $(1,1)$ over knots $\{0,0,1,1\}$. Find the four basis products at $(u,v)=(0.5,0.5)$. <br> <b>Formula:</b> degree-1 hats $N_0(u)=1-u,\ N_1(u)=u$; same in $v$: $M_0(v)=1-v,\ M_1(v)=v$ <br> <i>Hint:</i> evaluate each hat at $0.5$ first, then multiply $N_i\cdot M_j$ for the four combinations. :: A: At $u=0.5$: $N_0=N_1=0.5$; at $v=0.5$: $M_0=M_1=0.5$ <br> each product $=0.5\cdot0.5$ <br> $\boxed{N_iM_j=0.25\ \text{for all four}\ (i,j)}$.
- Q: Worked example (bilinear, step 2 of 2): with $P_{00}=(0,0,0),P_{10}=(1,0,2),P_{01}=(0,1,1),P_{11}=(1,1,4)$ and each basis product $=0.25$, evaluate $S(0.5,0.5)$. <br> <b>Formula:</b> $S(0.5,0.5)=0.25\,(P_{00}+P_{10}+P_{01}+P_{11})$ <br> <i>Hint:</i> add the four control points component-by-component first, then multiply the sum by $0.25$. :: A: Sum $=(0{+}1{+}0{+}1,\,0{+}0{+}1{+}1,\,0{+}2{+}1{+}4)=(2,2,7)$ <br> multiply by $0.25$ <br> $\boxed{S(0.5,0.5)=(0.5,0.5,1.75)}$ (center $z$ = mean of corner heights).
- Q: Worked example (centripetal, irregular row): three chords measure $1,1,8$ (one big gap). Compute the four centripetal parameters. <br> <b>Formula:</b> $\bar u_r=\dfrac{\sum_{k=1}^{r}\sqrt{\text{chord}_k}}{\sum_{k=1}^{R}\sqrt{\text{chord}_k}}$ <br> <i>Hint:</i> take the square root of each chord ($\sqrt8\approx2.828$), then accumulate and divide by the total of the roots. :: A: Roots of chords: $1,1,2.828$, total $=4.828$ <br> cumulative over total: $0,\ 1/4.828,\ 2/4.828,\ 4.828/4.828$ <br> $\boxed{\bar u=\{0,\,0.207,\,0.414,\,1.0\}}$ (crowded samples spread over $[0,0.414]$).
- Q: Worked example (chord-length, same row): chords $1,1,8$, total $10$. Compute the chord-length parameters and compare to centripetal. <br> <b>Formula:</b> $\bar u_r=\dfrac{\sum_{k=1}^{r}|Q_k-Q_{k-1}|}{d}$, $d=\sum|Q_k-Q_{k-1}|$ <br> <i>Hint:</i> here $\alpha=1$, so use raw distances (no square root); accumulate and divide by total $10$. :: A: Cumulative raw distance over total: $0,\ 1/10,\ 2/10,\ 10/10$ <br> $\bar u=\{0,0.1,0.2,1.0\}$ <br> $\boxed{\text{crowded into }[0,0.2]\text{ vs }[0,0.414]\text{ centripetal}}$ -- chord-length squeezes harder, inviting overshoot.
- Q: Worked example (mm vs m scale): a row spans chords $2,2,16$ in mm; recompute in m ($0.002,0.002,0.016$). Do the centripetal parameters change? <br> <b>Formula:</b> $\bar u_r=\dfrac{\sum_{k=1}^{r}\sqrt{\text{chord}_k}}{\sum_{k=1}^{R}\sqrt{\text{chord}_k}}$ <br> <i>Hint:</i> compute the normalized ratio in each unit; a common scale factor cancels top and bottom. :: A: mm: $\sqrt{}=\{1.414,1.414,4\}$, total $6.828$; m: $\sqrt{}=\{0.0447,0.0447,0.1265\}$, total $0.2159$ <br> ratios identical: $1.414/6.828=0.0447/0.2159$ <br> $\boxed{\bar u=\{0,0.207,0.414,1.0\}\ \text{either unit}}$ -- params are scale-invariant.
- Q: Worked example (knot count): a cubic ($p=3$) B-spline has $n+1=10$ control points in $u$. How many knots $m+1$ are there? <br> <b>Formula:</b> $m=n+p+1$ <br> <i>Hint:</i> read off $n=9$ and $p=3$, plug into $m=n+p+1$, then the knot count is $m+1$. :: A: $n=9,\ p=3\Rightarrow m=9+3+1=13$ <br> knot count $=m+1=14$ <br> $\boxed{m+1=14\ \text{knots}}$.
- Q: Test yourself (bilinear corner check): evaluate the worked patch at $(u,v)=(1,0)$ and verify it returns $P_{10}$. <br> <b>Formula:</b> $N_0(1)=0,N_1(1)=1$; $M_0(0)=1,M_1(0)=0$; $S=\sum_{ij}N_iM_jP_{ij}$ <br> <i>Hint:</i> a pass means only one product survives -- find the $(i,j)$ where both hats equal $1$. :: A: Products: only $N_1M_0=1\cdot1=1$ is nonzero <br> so $S=1\cdot P_{10}$ <br> $\boxed{S(1,0)=(1,0,2)=P_{10}}$ -- bilinear interpolates its corners. Loop closed.
- Q: Test yourself (partition of unity at the center): for the bilinear patch at $(0.5,0.5)$, verify the four basis products sum to 1. <br> <b>Formula:</b> $\sum_{ij}N_i(u)M_j(v)=1$ <br> <i>Hint:</i> a pass means the four products ($0.25$ each) total exactly $1.00$, so $S$ is a convex average of the corners. :: A: Each product is $0.25$ <br> $0.25\times4=1.00$ <br> $\boxed{\sum_{ij}N_iM_j=1}$ -- so $S$ is an affine (convex) average of the corners, inside their hull. Verified.
- Q: Test yourself (centripetal endpoint): confirm the irregular-row centripetal parameters end at exactly $\bar u_R=1$. <br> <b>Formula:</b> $\bar u_R=\dfrac{\sum_{k=1}^{R}\sqrt{\text{chord}_k}}{\sum_{k=1}^{R}\sqrt{\text{chord}_k}}$ <br> <i>Hint:</i> a pass means the final numerator equals the denominator (full total), so the ratio is exactly $1$. :: A: Final partial sum $=$ total: $1+1+2.828=4.828$ <br> $4.828/4.828=1$ <br> $\boxed{\bar u_3=1.0}$ -- normalization guarantees the last parameter is $1$. Loop closed.
- Q: Test yourself (fit-error gate): after a fit, you measure max sample-to-surface distance $0.6$ mm against $\text{Tol3D}=1.0$ mm. Does the fit pass? <br> <b>Formula:</b> pass iff $\max_{r,s}\|S(\bar u_r,\bar v_s)-Q_{rs}\|\le\text{Tol3D}$ <br> <i>Hint:</i> a pass means the largest residual is at or below Tol3D; compare $0.6$ to $1.0$. :: A: Max residual $=0.6$ mm, Tol3D $=1.0$ mm <br> $0.6\le1.0$, so the tolerance is met <br> $\boxed{\text{PASS}}$ -- OCCT stops raising the control count once $\max\le\text{Tol3D}$.
- CLOZE: The B-spline count relation linking knots, control points, and degree is {{c1::$m=n+p+1$}}, with $n+1$ control points and $m+1$ knots in $u$.
- CLOZE: A non-rational B-spline surface is the special case of a NURBS surface where all weights are {{c1::equal (e.g. $w_{ij}=1$)}}, making the denominator collapse to 1.
- CLOZE: Surface approximation solves the normal equations {{c1::$B^\top B\,P=B^\top Q$}} once per coordinate $x,y,z$. <br> <i>hint:</i> least squares, set the gradient of $\|BP-Q\|^2$ to zero.
- CLOZE: Centripetal parameterization raises each chord to the power {{c1::$\alpha=1/2$}}, compressing the dynamic range of irregular spacings.
- CLOZE: In OCCT the centripetal parameterization type is the enum {{c1::Approx_Centripetal}}, and the surface fitter is {{c2::GeomAPI_PointsToBSplineSurface}}.
- Q: When to use approximation vs interpolation for a scan, and why? <br> <i>Hint:</i> think about what interpolation does to noisy samples. :: A: Approximate (fewer control points than data) -- interpolation forces the surface through every noisy sample, turning noise into geometry and rippling; approximation acts as a low-pass filter.
- Q: Pitfall -- you pass an unorganized cloud straight to the surface fitter. What breaks and what is the fix? <br> <i>Hint:</i> recall what the array index $(i,j)$ is required to be. :: A: The index $(i,j)$ must be a true topological grid (row monotone in $u$, column in $v$); a raw cloud yields garbage. Fix: template-guided UV resample onto an ordered grid first.
- Q: Pitfall -- a grid in mm fit with $\text{Tol3D}=1.0\text{e-}3$. What goes wrong? <br> <i>Hint:</i> what units is Tol3D in, and how small is $0.001$ in mm? :: A: Tol3D is in model units, so $0.001$ mm over-refines and explodes the control count. Keep grid units and Tol3D in one system; set Tol3D to a physical fraction of feature size.
- Q: Why must an overhang be split into separate patches before fitting? <br> <i>Hint:</i> count how many surface points a tensor-product surface allows per $(u,v)$. :: A: A tensor-product surface is single-valued: each $(u,v)$ maps to exactly one point. An overhang has two surface points over one template location, so one patch cannot represent it -- segment into height-field patches.
- CLOZE: OCCT surface-fit defaults are DegMin={{c1::3}}, DegMax={{c2::8}}, continuity {{c3::GeomAbs_C2}}, Tol3D={{c4::1.0e-3}}.
