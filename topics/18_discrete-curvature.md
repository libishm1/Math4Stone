# Discrete curvature: angle defect, mean curvature, normals

Prereqs: [[mesh-topology]], [[trig-rotations]] => this => Unlocks: [[laplacian-poisson]], [[remeshing-segmentation]]

## Why it matters for the stone scan->CAD->IFC / CADFit workflow

A raw stone scan is a triangle soup with no notion of where the sharp arrises, the flat dressed faces, or the rounded broken edges are. Discrete curvature is the cheapest reliable signal that turns vertices into those semantic regions. Angle defect (integrated Gaussian curvature) flags corners and cone points; discrete mean curvature flags ridge and valley creases; the cotangent mean-curvature normal drives the Laplacian smoothing you run before you ever fit a NURBS patch. In the pipeline this is the gate between "cleaned mesh" and "feature-aware segmentation": you threshold curvature to split the stone into chartable patches, then hand each near-developable patch to `GeomAPI_PointsToBSplineSurface` so a single spline is not forced across a fold. The same curvature scalars feed a CADFit-style synthesis layer as per-region descriptors that say "this is a planar face," "this is a fillet," "this is an apex," which is exactly the program-prediction prior a CAD-from-scan model needs.

## Intuition

Imagine standing on a surface and walking a tiny loop around your feet. On a flat floor the directions you faced sum to a full turn, $2\pi$. On the tip of a cone, the cone "is missing a wedge," so the angles of the triangles meeting at the tip add up to less than $2\pi$. That shortfall is the angle defect, and it is concentrated Gaussian curvature. A saddle has too much angle, so the defect goes negative.

Mean curvature is a different question: not "how much does the surface fail to be flat in a topological sense," but "how hard is the surface bending, and which way." Picture pinching a sheet of paper into a crease. The crease has zero Gaussian curvature (paper does not stretch) but large mean curvature (it bends sharply). The clean analogy: angle defect is about the *cone* at a vertex; mean curvature is about the *crease* along an edge. A vertex normal is just the averaged "up" direction, and how you weight the average (by triangle area, by angle, or by the cotan operator) changes how robust it is to bad triangles.

## The math (first principles)

### Assumptions and conventions

Work on a manifold triangle mesh $M = (V, E, F)$, oriented, with consistent counter-clockwise face winding so face normals point "outward." Coordinates are in the scan's working units (meters for a quarried block; keep them consistent). Curvature has units of inverse length (mean) or inverse length squared (Gaussian) once made pointwise; the *integrated* quantities below are unitless angles (Gaussian) or length times angle (mean). All angles are in radians. Right-handed frames throughout. Tolerance note: every formula has a cotangent or a divide-by-area, so degenerate (sliver or zero-area) triangles must be guarded.

### Angle defect = integrated Gaussian curvature

For interior vertex $i$, let $\phi_i^{jk}$ be the interior angle at $i$ in triangle $ijk$. The discrete Gaussian curvature (angle defect) is

$$K_i = 2\pi - \sum_{ijk \in F} \phi_i^{jk}.$$

This is the *integrated* Gaussian curvature over the dual cell of $i$, written $\int_{A_i} K \, dA$. It is a discrete 2-form: a number attached to a region, not a pointwise density. $K_i > 0$ is convex/cone-like, $K_i < 0$ is saddle-like, $K_i = 0$ is developable (flat or bendable without stretch). For a boundary vertex with total interior angle $\Theta_i$, the defect is $K_i = \pi - \Theta_i$ instead of $2\pi - \Theta_i$.

To recover a pointwise density, divide by a dual area $A_i$:

$$\kappa_i \approx \frac{K_i}{A_i} = \frac{2\pi - \sum_{ijk} \phi_i^{jk}}{A_i}.$$

### Discrete Gauss-Bonnet

Summing the defect over a closed mesh recovers topology, independent of the triangulation:

$$\sum_{i \in V} K_i = 2\pi\,\chi(M), \qquad \chi(M) = V - E + F = 2 - 2g.$$

Here $\chi$ is the Euler characteristic and $g$ the genus. A sphere ($g=0$) gives $\sum_i K_i = 4\pi$ no matter how you mesh it. This is your free correctness check: compute every angle defect, sum, and confirm you land on $2\pi\chi$. If you do not, your topology or your angle bookkeeping is wrong.

*Sketch of why:* each triangle's three interior angles sum to $\pi$, so $\sum_F \pi = \pi F$. Each vertex contributes $2\pi$ minus its angle sum, so $\sum_i K_i = 2\pi V - \pi F$. On a closed triangle mesh every face has 3 edges and every edge is shared by 2 faces, so $3F = 2E$, giving $\pi F = 2\pi(E - F)$. Substituting, $\sum_i K_i = 2\pi(V - E + F) = 2\pi\chi$.

### Dihedral angle

Curvature along an edge needs the bend between its two faces. For edge $ij$ shared by triangles $ijk$ and $ijl$ with unit face normals $N_{ijk}, N_{ijl}$, the signed dihedral angle is

$$\theta_{ij} = \operatorname{atan2}\!\left(\frac{e_{ij}}{\lVert e_{ij}\rVert}\cdot\left(N_{ijk}\times N_{ijl}\right),\; N_{ijk}\cdot N_{ijl}\right),$$

where $e_{ij}$ is the edge vector from $i$ to $j$. Using `atan2` (not `acos`) keeps the sign: convex edges and concave edges get opposite signs, which is what distinguishes a ridge from a valley.

### Discrete mean curvature

The scalar (integrated) mean curvature at vertex $i$ sums dihedral angle times edge length over incident edges:

$$H_i = \frac{1}{2}\sum_{ij \in E} \theta_{ij}\,\lVert e_{ij}\rVert.$$

Like $K_i$, this is an integrated quantity; divide by $A_i$ for a pointwise density.

### Vertex normals (several weightings)

The "true" normal at a vertex is undefined for a piecewise-flat surface, so we average face normals $N_{ijk}$. Three standard weightings:

Area-weighted (default, cheap, robust to thin triangles):
$$N_i^A = \sum_{ijk \in F} A_{ijk}\,N_{ijk}, \quad\text{then normalize.}$$

Angle-weighted (uses the tip angle $\phi_i^{jk}$ at $i$; less sensitive to long skinny triangles around the vertex):
$$N_i^\phi = \sum_{ijk \in F} \phi_i^{jk}\,N_{ijk}, \quad\text{then normalize.}$$

Sphere-inscribed (Max's formula; weights by inverse squared edge lengths, good near sharp features):
$$N_i^S = \sum_{ijk \in F} \frac{e_{ij}\times e_{ik}}{\lVert e_{ij}\rVert^2\,\lVert e_{ik}\rVert^2}, \quad\text{then normalize.}$$

All three return a direction; you must normalize. They differ only in weighting, and the differences show up exactly where the mesh is anisotropic, which is common on scans.

### The cotangent mean-curvature normal (the key operator)

The deepest identity links mean curvature to the Laplacian of position, $\Delta \mathbf{x} = 2H\mathbf{n}$. Its discrete form is the famous cotan formula. The mean-curvature normal at $i$ is

$$HN_i = \frac{1}{2}\sum_{ij \in E}\left(\cot\alpha_{ij} + \cot\beta_{ij}\right)e_{ij},$$

where $\alpha_{ij}$ and $\beta_{ij}$ are the two angles *opposite* edge $ij$ in its two incident triangles. Pointwise mean curvature magnitude is then

$$H_i = \frac{1}{2}\,\frac{\lVert HN_i\rVert}{A_i},$$

with $A_i$ the dual area. This same cotan sum is the row of the cotangent Laplacian $L$, which is why this topic unlocks [[laplacian-poisson]]: smoothing position by $\mathbf{x} \leftarrow \mathbf{x} - \tau\,L\mathbf{x}$ is mean-curvature flow.

There is also a Gaussian-curvature normal built from dihedral angles:
$$KN_i = \frac{1}{2}\sum_{ij \in E}\frac{\theta_{ij}}{\lVert e_{ij}\rVert}\,e_{ij}.$$

### Dual area: barycentric vs Voronoi vs mixed

The pointwise conversion needs an area $A_i$ attached to vertex $i$. Three choices:

Barycentric: one third of each incident triangle's area. Simple, always positive, but a crude approximation.
$$A_i^{\text{bary}} = \frac{1}{3}\sum_{ijk \in F} A_{ijk}.$$

Circumcentric (Voronoi): the region closer to $i$ than to any other vertex. For a triangle $ijk$ its contribution to $A_i$ is
$$\frac{1}{8}\left(\lVert e_{ij}\rVert^2\cot\beta_{ij} + \lVert e_{ik}\rVert^2\cot\alpha_{ik}\right),$$
where $\alpha, \beta$ are the angles opposite those edges. Summed over incident faces this gives the Crane circumcentric dual area. It matches the cotan operator exactly but can go negative on obtuse triangles, because the circumcenter falls outside the triangle.

Mixed (Meyer et al., the production choice): use the Voronoi region when the triangle is non-obtuse; when it is obtuse, fall back to triangle-area fractions. For triangle $T$ at vertex $i$:
$$A_i \mathrel{+}= \begin{cases} \text{Voronoi region of } T \text{ at } i & T \text{ non-obtuse} \\ \tfrac{1}{2}\,\operatorname{area}(T) & T \text{ obtuse at } i \\ \tfrac{1}{4}\,\operatorname{area}(T) & T \text{ obtuse, not at } i \end{cases}$$

The mixed area is always positive and stays accurate, which is why it is the default for scan meshes that contain obtuse slivers.

## Worked example

Take a unit-square-based pyramid apex. Vertex $i$ is the apex; four identical isosceles triangles meet there. Suppose each triangle's tip angle at the apex is $\phi = 80^\circ = 1.3963$ rad.

Angle defect:
$$K_i = 2\pi - 4(1.3963) = 6.2832 - 5.5851 = 0.6981 \text{ rad}.$$
Positive, as expected for a convex apex. Equivalently $360^\circ - 4(80^\circ) = 40^\circ$ of missing angle.

Now flatten the cone: if the four triangles were coplanar each tip angle would be $90^\circ$, giving $K_i = 2\pi - 4(\pi/2) = 0$. Good: a flat vertex has zero defect.

Make it a saddle instead: four triangles with tip angle $100^\circ$ each. Then
$$K_i = 2\pi - 4(1.7453) = 6.2832 - 6.9813 = -0.6981 \text{ rad},$$
negative, the saddle signature.

Gauss-Bonnet sanity on a tetrahedron (topological sphere, $\chi = 2$): each of the 4 vertices has 3 equilateral triangles meeting, tip angle $60^\circ$ each, so $K_i = 2\pi - 3(\pi/3) = 2\pi - \pi = \pi$. Sum over 4 vertices: $4\pi = 2\pi\chi = 2\pi(2)$. Confirmed.

Pointwise: if the apex dual area is $A_i = 0.5$ m$^2$, the Gaussian curvature density is $\kappa_i = 0.6981 / 0.5 = 1.40$ m$^{-2}$.

## Code

Idiomatic C++ with Eigen. A mesh is `V` (n x 3 vertices) and `F` (m x 3 indices). The functions map directly to the formulas above.

```cpp
// discrete_curvature.hpp
#pragma once
#include <Eigen/Dense>
#include <vector>
#include <cmath>

namespace dcurv {

using Vec3 = Eigen::Vector3d;

// Interior angle at vertex `a` in triangle (a,b,c). Symbol: phi_i^{jk}.
inline double cornerAngle(const Vec3& a, const Vec3& b, const Vec3& c) {
    Vec3 u = (b - a), w = (c - a);
    double d = u.normalized().dot(w.normalized());
    d = std::max(-1.0, std::min(1.0, d));   // clamp against fp drift
    return std::acos(d);
}

// cot of the angle at vertex `a` in triangle (a,b,c). Used by cotan operator.
inline double cotAngle(const Vec3& a, const Vec3& b, const Vec3& c) {
    Vec3 u = (b - a), w = (c - a);
    double cosA = u.dot(w);
    double sinA = u.cross(w).norm();        // |u||w| sin = cross magnitude
    if (sinA < 1e-12) sinA = 1e-12;         // guard sliver triangles
    return cosA / sinA;
}

// Angle defect K_i = 2*pi - sum of incident tip angles.  (integrated Gaussian)
// vertexFaces[i] = list of (a,b,c) face index triples where a == i.
inline double angleDefect(int i,
                          const Eigen::MatrixXd& V,
                          const std::vector<Eigen::Vector3i>& incident) {
    double sum = 0.0;
    for (const auto& f : incident) {
        // ensure the corner is at vertex i
        int a = f[0], b = f[1], c = f[2];
        if (b == i) std::swap(a, b);
        else if (c == i) std::swap(a, c);
        sum += cornerAngle(V.row(a), V.row(b), V.row(c));
    }
    return 2.0 * M_PI - sum;                 // K_i
}

// Area-weighted vertex normal:  N_i^A = sum A_ijk * N_ijk, normalized.
inline Vec3 areaWeightedNormal(int i,
                               const Eigen::MatrixXd& V,
                               const std::vector<Eigen::Vector3i>& incident) {
    Vec3 n = Vec3::Zero();
    for (const auto& f : incident) {
        Vec3 p0 = V.row(f[0]), p1 = V.row(f[1]), p2 = V.row(f[2]);
        Vec3 fn = (p1 - p0).cross(p2 - p0);  // 2*A * unit-normal; area folded in
        n += fn;
    }
    return n.normalized();
}

// Mixed dual area at vertex i (Meyer et al.): Voronoi if non-obtuse, else T/2 or T/4.
inline double mixedArea(int i,
                        const Eigen::MatrixXd& V,
                        const std::vector<Eigen::Vector3i>& incident) {
    double A = 0.0;
    for (const auto& f : incident) {
        int a = f[0], b = f[1], c = f[2];
        if (b == i) std::swap(a, b);
        else if (c == i) std::swap(a, c);    // a == i now
        Vec3 P = V.row(a), Q = V.row(b), R = V.row(c);
        double Tarea = 0.5 * (Q - P).cross(R - P).norm();
        double angA = cornerAngle(P, Q, R);
        double angB = cornerAngle(Q, R, P);
        double angC = cornerAngle(R, P, Q);
        bool obtuse = (angA > M_PI/2) || (angB > M_PI/2) || (angC > M_PI/2);
        if (!obtuse) {
            // circumcentric (Voronoi) contribution at i:
            // (1/8)( |ib|^2 cot(angle at c) + |ic|^2 cot(angle at b) )
            double e_ib2 = (Q - P).squaredNorm();
            double e_ic2 = (R - P).squaredNorm();
            A += (e_ib2 * (1.0/std::tan(angC)) + e_ic2 * (1.0/std::tan(angB))) / 8.0;
        } else {
            A += (angA > M_PI/2) ? (Tarea / 2.0) : (Tarea / 4.0);
        }
    }
    return A;
}

// Cotan mean-curvature normal HN_i = 1/2 sum (cot a + cot b) e_ij.
// neighbors[i] = list of (j, oppK, oppL): edge endpoint and the two opposite verts.
inline Vec3 meanCurvatureNormal(int i,
                                const Eigen::MatrixXd& V,
                                const std::vector<Eigen::Vector3i>& edges) {
    Vec3 HN = Vec3::Zero();
    Vec3 xi = V.row(i);
    for (const auto& e : edges) {
        int j = e[0], k = e[1], l = e[2];
        Vec3 xj = V.row(j), xk = V.row(k), xl = V.row(l);
        double cotA = cotAngle(xk, xi, xj);  // angle opposite ij in triangle ijk
        double cotB = cotAngle(xl, xi, xj);  // angle opposite ij in triangle ijl
        HN += 0.5 * (cotA + cotB) * (xj - xi);
    }
    return HN;   // pointwise H = 0.5 * |HN| / A_i ; direction approximates the normal
}

} // namespace dcurv
```

Notes mapping symbols to code: `cornerAngle` is $\phi_i^{jk}$, `cotAngle` supplies $\cot\alpha_{ij}$, `angleDefect` is $K_i$, `mixedArea` is $A_i$, `meanCurvatureNormal` is $HN_i$. To get the Gauss-Bonnet check, sum `angleDefect` over all vertices and compare to $2\pi\chi$.

A short Python check using libigl makes the Gauss-Bonnet test a one-liner and is worth keeping as a reference oracle:

```python
import igl
import numpy as np
V, F = igl.read_triangle_mesh("stone.obj")
K = igl.gaussian_curvature(V, F)          # integrated angle defect per vertex
chi = V.shape[0] - igl.edges(F).shape[0] + F.shape[0]
print(K.sum(), 2*np.pi*chi)               # should match for a closed mesh
H = 0.5 * np.linalg.norm(                  # mean curvature via cotan Laplacian
        igl.cotmatrix(V, F) @ V, axis=1)
M = igl.massmatrix(V, F, igl.MASSMATRIX_TYPE_VORONOI).diagonal()
H = H / M                                   # pointwise
```

## Connections

- [[mesh-topology]] is the hard prerequisite: angle defect, dual areas, and the cotan operator all iterate over the one-ring and over edge-incident faces, which only exist once you have a half-edge or face-vertex adjacency.
- [[trig-rotations]] supplies the angle, dot/cross, and `atan2` machinery; the dihedral angle and every corner angle are trig identities on edge vectors.
- [[laplacian-poisson]] is the direct unlock: the cotan mean-curvature normal *is* a row of the cotangent Laplacian, so mean-curvature flow and Poisson smoothing reuse this exact operator.
- [[remeshing-segmentation]] consumes curvature as the feature signal: threshold $K_i$ and $H_i$ to find sharp edges and corners, then cut the mesh into near-developable charts for spline fitting.
- [[normal-estimation-pca]] is the point-cloud cousin: when you have no mesh connectivity you estimate normals by PCA instead of these one-ring averages; curvature there comes from the eigenvalue ratio.

## Pitfalls and failure modes

- Frame and winding inconsistency: if face winding is not coherent, face normals flip sign per triangle and every weighted vertex normal cancels into noise. Make the mesh orientable and consistently wound before computing anything.
- `acos` blowups: clamp the dot product to $[-1,1]$ before `acos`. Floating-point drift produces values like 1.0000001 that throw NaN.
- Cotangent on slivers: $\cot$ explodes as a triangle angle approaches 0 or $\pi$. Scan meshes are full of slivers. Guard the `sin` term and remesh degenerate triangles first; otherwise the cotan Laplacian becomes ill-conditioned and smoothing diverges.
- Negative circumcentric area: the pure Voronoi dual area goes negative on obtuse triangles. Use the mixed (Meyer) area in production, never raw circumcentric, or your pointwise $\kappa_i = K_i / A_i$ flips sign spuriously.
- Confusing integrated with pointwise: $K_i$ and $H_i$ from the formulas are integrated 2-forms (units of angle, or length times angle). Dividing by $A_i$ is mandatory to get curvature densities. Comparing an unnormalized $K_i$ across meshes of different resolution is meaningless.
- Boundary vertices: the $2\pi - \sum\phi$ formula is only correct for interior vertices. Boundary vertices need $\pi - \sum\phi$, and the cotan operator has missing half-edges. Mask the boundary or handle it explicitly, or Gauss-Bonnet will be off by the boundary's geodesic curvature.
- Non-commutativity of normalization and averaging: normalize *after* summing weighted face normals, not before. Normalizing each face normal first throws away the area or angle weighting you intended.
- Outliers from holes: a hole in the scan leaves a fan of boundary triangles whose defect is huge and meaningless. Fill or mask holes before trusting curvature near them.

## Practice

1. Implement `angleDefect` and verify discrete Gauss-Bonnet on three closed meshes: a tetrahedron (expect $4\pi$), an icosahedron (expect $4\pi$), and a torus (expect $0$). Print the sum and the target $2\pi\chi$ side by side.
2. Compute area-weighted, angle-weighted, and sphere-inscribed vertex normals on a mesh with one deliberately long sliver triangle. Measure the angular difference between the three at the sliver's shared vertex and explain which you trust.
3. On a smooth analytic surface sampled as a mesh (a unit sphere of known curvature $\kappa = 1$), compute pointwise $K_i / A_i$ using barycentric, circumcentric, and mixed dual areas. Tabulate the mean and the worst-case error against the true value for each area type.
4. Color a stone-scan mesh by signed dihedral angle and confirm that convex arrises and concave grooves get opposite signs. Inspect the output in a viewer.
5. Run a few steps of mean-curvature flow $\mathbf{x} \leftarrow \mathbf{x} - \tau\,A^{-1} HN$ on a noisy mesh and watch high-curvature noise shrink. Find the largest stable $\tau$ empirically and relate it to the smallest triangle.
6. Build a feature mask: threshold pointwise $|H_i|$ and angle defect $|K_i|$, then export only the flagged vertices as a point set. Verify it isolates the corners and creases you would cut along for patch fitting.

## References

- Keenan Crane, *Discrete Differential Geometry: An Applied Introduction* (CMU 15-458/858 course notes), https://www.cs.cmu.edu/~kmcrane/Projects/DDG/
- Keenan Crane, "Assignment 2 (Coding): Investigating Curvature," CS 15-458/858 (Spring 2019), exact normal, defect, dihedral, and dual-area formulas, https://brickisland.net/DDGSpring2019/2019/02/27/assignment-2-coding-investigating-curvature-due-3-19/
- Keenan Crane, "Discrete Curvature I (Integral)," Lecture 15 slides, https://brickisland.net/DDGSpring2019/wp-content/uploads/2019/03/DDG_458_SP19_Lecture15_DiscreteCurvatureI-1.pdf
- M. Meyer, M. Desbrun, P. Schroder, A. H. Barr, "Discrete Differential-Geometry Operators for Triangulated 2-Manifolds," *Visualization and Mathematics III*, Springer, 2003 (mixed Voronoi area, cotan mean-curvature normal), https://multires.caltech.edu/pubs/diffGeoOps.pdf
- S. Upadhyay, "Gauss-Bonnet for Discrete Surfaces," University of Chicago REU 2015, https://math.uchicago.edu/~may/REU2015/REUPapers/Upadhyay.pdf
- libigl tutorial, Gaussian curvature and cotangent Laplacian, https://libigl.github.io/tutorial/
- CGAL Weights package, Mixed Voronoi Region Weight, https://doc.cgal.org/latest/Weights/group__PkgWeightsRefMixedVoronoiRegionWeights.html
- OpenCASCADE reference manual, `GeomAPI_PointsToBSplineSurface` (downstream patch fitting), https://dev.opencascade.org/doc/overview/html/
- IfcOpenShell geometry processing docs (downstream IFC serialization), https://docs.ifcopenshell.org/
- N. Nehme, CADFit: CAD program synthesis from geometry, arXiv:2605.01171; code at https://github.com/ghadinehme/CADFit

## Flashcard seeds

- Q: Derive the angle defect $K_i$ at an interior vertex, starting from the flat-case reference that the tip angles around a flat vertex make a full turn. <br> <b>Formula:</b> start from $\sum_{ijk}\phi_i^{jk}=2\pi$ (flat) and target $K_i = 2\pi - \sum_{ijk}\phi_i^{jk}$ <br> <i>Hint:</i> write the actual angle sum as the flat value minus a shortfall $K_i$, then solve for $K_i$. :: A: A flat (developable) vertex tiles a full disc: $\sum_{ijk}\phi_i^{jk}=2\pi$ <br> a cone is missing a wedge, so the actual sum falls short by the defect: $\sum_{ijk}\phi_i^{jk}=2\pi-K_i$ <br> solve for the shortfall <br> $\boxed{K_i = 2\pi - \sum_{ijk}\phi_i^{jk}}$.
- Q: Show why $K_i$ is the *integrated* Gaussian curvature and derive the pointwise density $\kappa_i$ from it. <br> <b>Formula:</b> $K_i=\int_{A_i}K\,dA$ and target $\kappa_i \approx K_i/A_i$ <br> <i>Hint:</i> treat $K$ as roughly constant over the small dual cell, so the integral becomes $\kappa_i A_i$; then divide by $A_i$. :: A: $K_i$ is a number attached to the dual cell of $i$, i.e. $K_i=\int_{A_i}K\,dA$ (a discrete 2-form, units of angle) <br> assume $K$ roughly constant over the small cell: $\int_{A_i}K\,dA \approx \kappa_i\,A_i$ <br> divide by the dual area $A_i$ <br> $\boxed{\kappa_i \approx K_i / A_i}$.
- Q: Derive the per-face angle-sum total used in Gauss-Bonnet (step 1 of 3). <br> <b>Formula:</b> each triangle's interior angles sum to $\pi$; target $\sum_{\text{faces}}(\text{angles}) = \pi F$ <br> <i>Hint:</i> just sum the per-triangle value $\pi$ over all $F$ faces. :: A: Each triangle has interior angles summing to $\pi$ <br> sum over all $F$ faces <br> $\boxed{\sum_{\text{faces}} (\text{3 interior angles}) = \pi F}$.
- Q: Continue the Gauss-Bonnet derivation (step 2 of 3): eliminate the angle sum. <br> <b>Formula:</b> $\sum_i K_i = \sum_i\big(2\pi-\sum\phi\big)$ and $\sum_F(\text{angles})=\pi F$; target $\sum_i K_i = 2\pi V - \pi F$ <br> <i>Hint:</i> the total of all interior angles is the same whether you group by vertex or by face, so substitute $\pi F$ for $\sum\phi$. :: A: Every interior angle belongs to exactly one vertex's defect sum and one face's angle sum, so the two total angle sums are equal <br> $\sum_i K_i = 2\pi V - \sum_{\text{all angles}}\phi = 2\pi V - \pi F$ <br> $\boxed{\sum_i K_i = 2\pi V - \pi F}$.
- Q: Finish Gauss-Bonnet (step 3 of 3): reach $\chi$ from the angle/face count. <br> <b>Formula:</b> $\sum_i K_i = 2\pi V - \pi F$ with closed-mesh relation $3F=2E$; target $\sum_i K_i = 2\pi\chi$, $\chi=V-E+F$ <br> <i>Hint:</i> use $3F=2E$ to rewrite $\pi F = 2\pi(E-F)$, then substitute and factor $2\pi$. :: A: Closed triangle mesh: each face has 3 edges, each edge shared by 2 faces, so $3F=2E$, i.e. $\pi F = 2\pi(E-F)$ <br> substitute: $\sum_i K_i = 2\pi V - 2\pi(E-F) = 2\pi(V-E+F)$ <br> $\boxed{\sum_i K_i = 2\pi\,\chi(M)},\quad \chi = V-E+F$.
- CLOZE: Gauss-Bonnet makes the defect sum depend only on topology: $\sum_i K_i = 2\pi\chi$ with $\chi = {{c1::2-2g}}$. <br> <i>hint:</i> $\chi = V-E+F$ written in terms of the genus $g$.
- Q: Construct the scalar (integrated) mean curvature $H_i$ at a vertex from dihedral angles and edge lengths. <br> <b>Formula:</b> per-edge bend $\theta_{ij}\lVert e_{ij}\rVert$; target $H_i = \tfrac{1}{2}\sum_{ij}\theta_{ij}\,\lVert e_{ij}\rVert$ <br> <i>Hint:</i> sum the per-edge term over incident edges, then halve because each edge is shared by two vertices. :: A: Each incident edge bends by signed dihedral angle $\theta_{ij}$ over length $\lVert e_{ij}\rVert$, contributing curvature $\theta_{ij}\lVert e_{ij}\rVert$ <br> each edge is shared by two vertices, so split its contribution in half <br> sum over incident edges <br> $\boxed{H_i = \tfrac{1}{2}\sum_{ij}\theta_{ij}\,\lVert e_{ij}\rVert}$.
- Q: From the continuous identity $\Delta\mathbf{x}=2H\mathbf{n}$, write the discrete cotan mean-curvature normal it becomes. <br> <b>Formula:</b> cotan weights $\tfrac{1}{2}(\cot\alpha_{ij}+\cot\beta_{ij})$ with $\alpha_{ij},\beta_{ij}$ opposite edge $ij$; target $HN_i = \tfrac{1}{2}\sum_{ij}(\cot\alpha_{ij}+\cot\beta_{ij})\,e_{ij}$ <br> <i>Hint:</i> replace the Laplacian of position by its cotan-weighted one-ring sum, each weight multiplying the edge vector $e_{ij}=x_j-x_i$. :: A: Discretize the Laplacian of position over the one-ring with cotan weights; $\alpha_{ij},\beta_{ij}$ are the angles opposite edge $ij$ in its two triangles <br> $\tfrac{1}{2}\Delta\mathbf{x}$ collects $\tfrac{1}{2}(\cot\alpha+\cot\beta)$ per edge times the edge vector <br> $\boxed{HN_i = \tfrac{1}{2}\sum_{ij}(\cot\alpha_{ij}+\cot\beta_{ij})\,e_{ij}}$.
- Q: Derive the pointwise mean-curvature magnitude $H_i$ from $HN_i$. <br> <b>Formula:</b> $\lVert HN_i\rVert \approx 2H_i\,A_i$ (since $HN_i$ integrates $2H\mathbf{n}$ over the dual cell); target $H_i = \tfrac{1}{2}\lVert HN_i\rVert/A_i$ <br> <i>Hint:</i> take the magnitude of $HN_i$ and divide out the $2A_i$ to isolate $H_i$. :: A: $HN_i$ integrates $2H\mathbf{n}$ over the dual cell, so $\lVert HN_i\rVert \approx 2H_i\,A_i$ <br> solve for $H_i$ <br> $\boxed{H_i = \tfrac{1}{2}\,\lVert HN_i\rVert / A_i}$.
- Q: Derive the barycentric dual area $A_i^{\text{bary}}$ from splitting each triangle's area equally among its 3 vertices. <br> <b>Formula:</b> each incident face gives $\tfrac{1}{3}A_{ijk}$; target $A_i^{\text{bary}} = \tfrac{1}{3}\sum_{ijk}A_{ijk}$ <br> <i>Hint:</i> sum the one-third share over every face incident to $i$. :: A: Triangle $ijk$ has area $A_{ijk}$; give each of its 3 vertices an equal third <br> vertex $i$ collects one third from every incident face <br> $\boxed{A_i^{\text{bary}} = \tfrac{1}{3}\sum_{ijk}A_{ijk}}$.
- Q: State the sign convention of the angle defect and the flat-reference relation it starts from. <br> <b>Formula:</b> $K_i = 2\pi - \sum\phi_i^{jk}$, flat reference $\sum\phi=2\pi$ <br> <i>Hint:</i> ask whether the actual angle sum falls short of, equals, or exceeds $2\pi$. :: A: Given (flat reference): $\sum\phi_i^{jk}=2\pi$ for a developable vertex <br> shortfall $2\pi-\sum\phi>0$ means a convex cone, $<0$ means a saddle (too much angle), $=0$ means flat <br> $\boxed{K_i>0\ \text{cone},\ K_i<0\ \text{saddle},\ K_i=0\ \text{developable}}$.
- Q: Worked example (convex apex, radians). Four isosceles triangles meet at a pyramid apex, each tip angle $\phi=80^\circ=1.3963$ rad. Compute $K_i$. <br> <b>Formula:</b> $K_i = 2\pi - \sum\phi = 2\pi - 4\phi$ <br> <i>Hint:</i> multiply the tip angle by 4, then subtract from $2\pi=6.2832$. :: A: $K_i = 2\pi - 4\phi$ <br> $= 6.2832 - 4(1.3963) = 6.2832 - 5.5851$ <br> $\boxed{K_i \approx 0.6981\ \text{rad}}$ (positive = convex), i.e. $40^\circ$ of missing angle.
- Q: Worked example (saddle regime). Same four-triangle vertex but each tip angle is now $\phi=100^\circ=1.7453$ rad. Compute $K_i$ and name the sign. <br> <b>Formula:</b> $K_i = 2\pi - 4\phi$ <br> <i>Hint:</i> compute $4\phi$ first; if it exceeds $2\pi$ the defect goes negative. :: A: $K_i = 2\pi - 4(1.7453) = 6.2832 - 6.9813$ <br> $\boxed{K_i \approx -0.6981\ \text{rad}}$ <br> negative defect = saddle: the angles overflow a full turn.
- Q: Worked example (pointwise density, m scale). For the convex apex $K_i = 0.6981$ rad with dual area $A_i = 0.5\ \text{m}^2$, give $\kappa_i$ and its units. <br> <b>Formula:</b> $\kappa_i = K_i / A_i$ <br> <i>Hint:</i> divide the defect by the area; the units are (radians are dimensionless) inverse area. :: A: $\kappa_i = K_i / A_i = 0.6981 / 0.5$ <br> $\boxed{\kappa_i \approx 1.40\ \text{m}^{-2}}$ <br> units of inverse-length-squared, the Gaussian-curvature density.
- Q: Worked example (mm vs m, mean curvature of a sphere). A fillet locally matches a sphere of radius $R$; mean curvature is $H=1/R$. Compare a coarse arris $R=2\ \text{m}$ with a fine ground edge $R=50\ \text{mm}$. <br> <b>Formula:</b> $H = 1/R$ <br> <i>Hint:</i> convert $50\ \text{mm}$ to $0.05\ \text{m}$ first so both radii are in meters, then take reciprocals. :: A: Coarse: $H = 1/2 = 0.5\ \text{m}^{-1}$ <br> fine: $R=0.05\ \text{m}\Rightarrow H = 1/0.05 = 20\ \text{m}^{-1}$ <br> $\boxed{H_{\text{fine}}/H_{\text{coarse}} = 40\times}$ — tighter rounding reads as far higher mean curvature.
- Q: Closing test (tetrahedron). Predict $\sum_i K_i$ from per-vertex defects, then verify against Gauss-Bonnet. <br> <b>Formula:</b> $K_i = 2\pi - 3(\pi/3)$ per vertex; check $2\pi\chi$, $\chi=V-E+F$ <br> <i>Hint:</i> a pass means the 4-vertex sum equals $2\pi\chi$ with $\chi=4-6+4=2$. :: A: Each vertex: 3 equilateral triangles, tip $60^\circ=\pi/3$, so $K_i = 2\pi - 3(\pi/3) = 2\pi-\pi = \pi$ <br> sum over 4 vertices: $4\pi$ <br> verify: $\chi=V-E+F=4-6+4=2$, so $2\pi\chi = 4\pi$ <br> $\boxed{4\pi = 4\pi\ \checkmark}$.
- Q: Closing test (icosahedron). 12 vertices, 5 equilateral triangles per vertex. Compute one defect, the total, and check Gauss-Bonnet. <br> <b>Formula:</b> $K_i = 2\pi - 5(\pi/3)$; total $=12K_i$; check $2\pi\chi$ with $\chi=2$ <br> <i>Hint:</i> a pass means $12K_i$ matches $4\pi$ regardless of the triangulation. :: A: $K_i = 2\pi - 5(\pi/3) = 2\pi - 5\pi/3 = \pi/3 \approx 1.0472$ <br> total $= 12(\pi/3) = 4\pi \approx 12.566$ <br> check: closed genus-0, $\chi=2$, $2\pi\chi=4\pi$ <br> $\boxed{4\pi = 4\pi\ \checkmark}$ (triangulation-independent).
- Q: Closing test (limiting case). Show the angle-defect formula gives the right answer when the four pyramid triangles are coplanar (flatten the cone). <br> <b>Formula:</b> $K_i = 2\pi - 4\phi$ with $\phi=\pi/2$ when flat <br> <i>Hint:</i> a pass means the defect collapses to exactly $0$ for the flat case. :: A: Coplanar = each tip angle $90^\circ=\pi/2$ <br> $K_i = 2\pi - 4(\pi/2) = 2\pi - 2\pi$ <br> $\boxed{K_i = 0}$ — a flat interior vertex has zero defect, the developable limit. $\checkmark$
- Q: Why use `atan2` rather than `acos` to compute the signed dihedral angle $\theta_{ij}$? <br> <b>Formula:</b> $\theta_{ij} = \operatorname{atan2}\!\big((e_{ij}/\lVert e_{ij}\rVert)\cdot(N_{ijk}\times N_{ijl}),\ N_{ijk}\cdot N_{ijl}\big)$ <br> <i>Hint:</i> compare the output ranges of the two functions and ask which keeps a sign. :: A: `acos` returns only $[0,\pi]$ and loses sign <br> `atan2` of (cross-product-along-edge, dot-of-normals) returns the full $(-\pi,\pi]$ and keeps the sign <br> $\boxed{\text{sign separates convex ridges from concave valleys}}$.
- Q: When is the mixed (Meyer) dual area required instead of raw circumcentric area, and what does it substitute? <br> <b>Formula:</b> $A_i \mathrel{+}= \tfrac{1}{2}\text{area}(T)$ if obtuse at $i$, $\tfrac{1}{4}\text{area}(T)$ if obtuse elsewhere, else the Voronoi region <br> <i>Hint:</i> recall that the circumcenter leaves the triangle when it is obtuse, which is what makes the raw area go negative. :: A: Raw Voronoi area goes negative on obtuse triangles (circumcenter falls outside), flipping the sign of $\kappa_i=K_i/A_i$ <br> mixed area uses Voronoi when non-obtuse, else $\tfrac{1}{2}\text{area}(T)$ if obtuse at $i$, $\tfrac{1}{4}\text{area}(T)$ otherwise <br> $\boxed{\text{always positive, the production choice for sliver-ridden scans}}$.
- CLOZE: $K_i$ and $H_i$ from the formulas are {{c1::integrated}} quantities; you must divide by $A_i$ to get a pointwise {{c2::density}} comparable across mesh resolutions.
- CLOZE: For a boundary vertex the angle defect is {{c1::$\pi$}} minus the sum of interior angles, not $2\pi$. <br> <i>hint:</i> a boundary vertex sees only a half-disc instead of a full turn.
- CLOZE: The cotan mean-curvature normal $HN_i$ is exactly one row of the {{c1::cotangent Laplacian}}, which is why this topic unlocks mean-curvature flow.
- CLOZE: Normalize a weighted vertex normal {{c1::after}} summing the face normals, never before, or you discard the weighting.
- Q: Pitfall recall: name two code guards that keep the curvature operators numerically safe on scan meshes. <br> <b>Formula:</b> clamp dot to $[-1,1]$ before `acos`; floor $\sin=\lVert u\times w\rVert$ at $\varepsilon$ before dividing in $\cot=\cos/\sin$ <br> <i>Hint:</i> one guard protects the `acos` domain, the other protects the `cot` divide on slivers. :: A: Clamp the dot product to $[-1,1]$ before `acos` (kills NaN from fp drift) <br> floor the `sin` (cross-product magnitude) at a small $\varepsilon$ before dividing in `cot` <br> $\boxed{\text{guard }acos\text{ domain and }cot\text{ slivers}}$.
