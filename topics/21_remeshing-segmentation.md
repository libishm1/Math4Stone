# Remeshing and segmentation (CVT, anisotropic, curvature-guided)

Prereqs: [[discrete-curvature]], [[sdf-signed-heat]], [[robust-ransac]] => this => unlocks [[nurbs-surfaces-fitting]], [[cadfit-program-synthesis]].

## Why it matters for the stone scan->CAD->IFC / CADFit workflow

A raw stone scan is a dense, ragged triangle soup: millions of skinny triangles, uneven sampling, and no notion of which face is a sawn plane and which is a natural fracture. Before you can fit anything, you must do two things. First, **remesh** the soup into well-shaped triangles whose density follows the surface, so that downstream normal estimation, parameterization, and fitting are numerically stable. Second, **segment** the mesh into regions that behave the same way: split the smooth sawn, contact, and template faces from the rough natural faces. That split is the routing decision of the whole pipeline. Smooth low-curvature regions go to NURBS or analytic primitive fitting and become editable CAD via [[cadfit-program-synthesis]]; high-curvature rough regions stay freeform meshes or become a tessellated IFC fallback. Get the remesh wrong and curvature is noise; get the segmentation wrong and you either over-fit a plane to a rough boulder face or waste a freeform patch on a flat slab.

## Intuition

Think of remeshing as re-tiling a floor. The original scan tiled the floor with thousands of random shards. You want to re-tile it with tiles that are all roughly the same nice shape, but where the floor curves sharply you want **smaller** tiles, and where it is flat you can use **big** tiles. Centroidal Voronoi tessellation (CVT) is the rule that makes the tiles regular: drop seed points, give each one the territory closest to it, then slide each seed to the center of mass of its own territory, and repeat. The seeds spread out evenly, like repelling magnets settling into a hexagonal grid.

Anisotropic remeshing adds one twist: on a long flat ridge you want **stretched** tiles aligned with the ridge, not squares. The trick is to lift each point into a 6D space where its position is glued to its surface normal. Two points that sit close in 3D but face different directions become far apart in 6D. Run ordinary even-spacing CVT in 6D, project back to 3D, and the tiles come out automatically stretched along low-curvature directions and small across sharp features.

Segmentation is then a flood-fill with a fence. Curvature builds fences along sharp creases. Water poured on a flat sawn face floods the whole face and stops at the crease where it meets the rough rock. Each puddle is one segment.

## The math (first principles)

**Conventions.** Surfaces are oriented 2-manifold triangle meshes $M=(V,E,F)$ with outward unit normals. Points are column vectors in $\mathbb{R}^3$, right-handed frame, units in millimetres at shop scale. Curvature uses the [[discrete-curvature]] conventions (cotangent Laplacian, mean curvature normal, angle defect for Gaussian curvature). Tolerances are relative to a characteristic length $L$ (object diagonal): use $\varepsilon = \max(\varepsilon_{\text{floor}},\, 10^{-3} L)$.

### Voronoi and centroidal Voronoi tessellation

Given generators (seeds) $X=\{x_1,\dots,x_k\}\subset\Omega$, the **Voronoi region** of $x_i$ is

$$V_i = \{\, y \in \Omega : \lVert y - x_i\rVert \le \lVert y - x_j\rVert \ \ \forall j \,\}.$$

With a density $\rho(y) \ge 0$, the **mass** and **centroid** of a region are

$$m_i = \int_{V_i}\rho(y)\,dy, \qquad c_i = \frac{1}{m_i}\int_{V_i} y\,\rho(y)\,dy.$$

A tessellation is **centroidal** when every generator sits at its own region's centroid, $x_i = c_i$ for all $i$. The quantity CVT minimizes is the **quantization energy** (Du, Faber, Gunzburger 1999):

$$F(X) = \sum_{i=1}^{k}\int_{V_i}\rho(y)\,\lVert y - x_i\rVert^2\,dy.$$

This is the variance of points around their nearest seed. It is the same least-squares error you have seen everywhere, summed per region.

**Gradient.** Although $V_i$ depends on $X$, the cell-boundary terms cancel (each boundary point is equidistant from two seeds), so the gradient has a clean closed form:

$$\frac{\partial F}{\partial x_i} = 2\,m_i\,(x_i - c_i).$$

Two facts fall out. First, the energy is **stationary** ($\nabla F = 0$) exactly when $x_i = c_i$, i.e. at a CVT. Second, $F$ is $C^2$ on a convex domain (Liu et al. 2009), so you can use L-BFGS rather than plain Lloyd.

**Lloyd's algorithm** is the fixed-point iteration that follows the gradient with a fixed full step:

1. Build the Voronoi diagram of $X$.
2. Compute each centroid $c_i$ by numerical integration over $V_i$.
3. Set $x_i \leftarrow c_i$.
4. Repeat until $\max_i \lVert x_i - c_i\rVert < \varepsilon$.

Each step never increases $F$, so it converges to a critical point (usually a good local minimum). It is gradient descent with step that needs no line search.

### Restricted CVT for surface remeshing

For a surface $S$ embedded in $\mathbb{R}^3$ you do not want full-3D cells, you want cells **restricted to the surface**. The **restricted Voronoi diagram** (RVD) intersects each 3D Voronoi cell with the surface:

$$V_i\big|_S = V_i \cap S.$$

Surface CVT minimizes the same energy but integrates over $V_i|_S$, and the centroids are projected back onto $S$. The output triangulation is the dual (the restricted Delaunay triangulation) of the converged seeds. With uniform density this yields near-equilateral triangles; with $\rho$ tied to curvature it yields a **curvature-adaptive** mesh.

**Sizing field.** To put more triangles where the surface bends, drive density by curvature. A standard choice uses the maximum absolute principal curvature $\kappa_{\max}$:

$$\rho(y) \;\propto\; \kappa_{\max}(y)^{\,2\beta/(2\beta+\dim)}, \qquad \text{often } \rho \propto \kappa_{\max}^{\,4/5} \ \text{for surfaces}.$$

The edge length then scales like $h(y)\propto \rho(y)^{-1/\dim}$, so flat regions get long edges and creases get short ones.

### Anisotropic remeshing by normal lifting (the 6D trick)

Isotropic CVT makes round triangles. Real stone slabs want triangles **stretched** along flat directions and **short** across sharp ones. Lévy and Bonneel (2013) and the feature-sensitive extension (Lévy 2015) do this without any UV parameterization by lifting the surface into 6D.

Define the lifting map from a point $x$ with unit normal $N(x)$:

$$\Phi(x) = \big(\,x,\; s\,N(x)\,\big) \in \mathbb{R}^6,$$

where $s \ge 0$ is the **anisotropy strength**. The 6D distance between two lifted points is

$$\lVert \Phi(x) - \Phi(x')\rVert^2 = \lVert x - x'\rVert^2 + s^2\,\lVert N(x) - N(x')\rVert^2.$$

The first term is ordinary 3D distance. The second penalizes **normal variation**. Where the surface is flat, normals barely change, so the extra term is tiny and CVT is free to space seeds far apart along that direction (big stretched triangles). Where the surface bends, normals swing fast, the second term blows up, and CVT packs seeds tightly across the bend (small triangles aligned with the feature). Anisotropy is curvature adaptivity in disguise: $\lVert N(x)-N(x')\rVert$ is a first-order proxy for curvature times arc length.

Run isotropic CVT in $\mathbb{R}^6$ (build the restricted Voronoi diagram of the lifted surface), then **project** the result to the first three coordinates to get the anisotropic 3D mesh. No charts, no seams, triangles may cross any chart boundary.

In **Geogram** the strength is set as a scaled parameter. The pipeline is: compute normals (append $n_x,n_y,n_z$ to each vertex so vertices are 6D), Laplacian-smooth the normals (critical for noisy reconstructed scans), set anisotropy via `set_anisotropy(mesh, anisotropy*0.02)` (so the effective $s$ is `anisotropy*0.02`), then `remesh_smooth`. Anisotropy of $0$ recovers isotropic CVT remeshing.

### Curvature-guided segmentation

Now split rough from smooth. Build a per-face roughness signal from curvature. Two standard signals:

- **Dihedral angle** across edge $e$ shared by faces $f,g$ with normals $n_f,n_g$:
$$\theta_e = \arccos\!\big(n_f \cdot n_g\big) \in [0,\pi].$$
A small $\theta_e$ means the two faces are nearly coplanar; a large $\theta_e$ marks a crease.

- **Discrete mean curvature** $H$ from the cotangent Laplacian (see [[discrete-curvature]]); large $|H|$ marks rough or highly curved regions.

The **segmentation rule** (the report's selectivity rule, which decides freeform vs CADFit routing):

1. Mark an edge as a **feature edge** if $\theta_e > \theta_{\text{crease}}$ (a common default is $\theta_{\text{crease}}\approx 30^\circ$ for sawn-vs-natural splits; tune per dataset).
2. **Region-grow** over faces: start a seed face, add a neighbour face only if the shared edge is **not** a feature edge and the candidate's local curvature stays below a roughness threshold $\tau$. The grow criterion is $\theta_e \le \theta_{\text{crease}}$ **and** $\overline{|H|}_{\text{region}} \le \tau$.
3. Classify each finished region: a region whose curvature statistics are low and whose RANSAC plane/cylinder fit (see [[robust-ransac]]) leaves residuals under $\varepsilon$ is a **smooth/sawn/contact/template face** -> route to analytic primitive or NURBS fitting. A region with high curvature variance and no consensus primitive is a **rough/natural face** -> keep as freeform mesh or tessellated fallback.

This is the gate: low-curvature consensus regions become editable CAD; the rest stay mesh.

## Worked example

Take a tiny 1D CVT so the arithmetic is checkable by hand. Domain $\Omega=[0,1]$, density $\rho\equiv 1$, two generators $x_1=0.2$, $x_2=0.6$.

**Voronoi split.** The boundary sits at the midpoint $b=(x_1+x_2)/2 = 0.4$. So $V_1=[0,0.4]$, $V_2=[0.4,1]$.

**Centroids.** With uniform density the centroid of an interval is its midpoint:
$$c_1 = \frac{0+0.4}{2}=0.2, \qquad c_2=\frac{0.4+1}{2}=0.7.$$

**Gradient check.** Masses are the lengths: $m_1=0.4$, $m_2=0.6$.
$$\frac{\partial F}{\partial x_1}=2 m_1 (x_1-c_1)=2(0.4)(0.2-0.2)=0,$$
$$\frac{\partial F}{\partial x_2}=2 m_2 (x_2-c_2)=2(0.6)(0.6-0.7)=-0.12.$$
Generator 1 is already at its centroid; generator 2 must move right (negative gradient means increase $x_2$).

**Lloyd step.** Set $x_i\leftarrow c_i$: new $x_1=0.2$, $x_2=0.7$. New boundary $b=0.45$, new cells $[0,0.45],[0.45,1]$, new centroids $c_1=0.225$, $c_2=0.725$. Still moving, so iterate.

**Energy drop.** Before the step:
$$F=\int_0^{0.4}(y-0.2)^2dy+\int_{0.4}^{1}(y-0.6)^2dy.$$
First integral $=\big[(y-0.2)^3/3\big]_0^{0.4}=(0.2^3+0.2^3)/3=0.016/3\cdot? $ compute: $(0.2)^3=0.008$, and at $y=0$ term is $(-0.2)^3/3=-0.008/3$, so first $=0.008/3-(-0.008/3)=0.016/3\approx 0.005333$. Second integral $=\big[(y-0.6)^3/3\big]_{0.4}^{1}=(0.4^3-(-0.2)^3)/3=(0.064+0.008)/3=0.072/3=0.024$. Total $F\approx 0.029333$.

After moving $x_2$ to $0.7$ (cells $[0,0.45],[0.45,1]$, seeds $0.2,0.7$): first $=\big[(y-0.2)^3/3\big]_0^{0.45}=((0.25)^3+(0.2)^3)/3=(0.015625+0.008)/3\approx 0.007875$; second $=\big[(y-0.7)^3/3\big]_{0.45}^{1}=((0.3)^3+(0.25)^3)/3=(0.027+0.015625)/3\approx 0.014208$. Total $F\approx 0.022083$. Energy fell from $0.0293$ to $0.0221$, as Lloyd guarantees.

**Anisotropy sanity check.** Two adjacent points on a flat region: $N(x)=N(x')$, so the 6D extra term $s^2\lVert N-N'\rVert^2=0$ and they may be far apart. Across a $90^\circ$ crease: $\lVert N-N'\rVert=\sqrt{2}$, so with $s=0.5$ the extra squared distance is $0.25\cdot 2=0.5$, pushing seeds apart in 6D and forcing extra triangles across the crease. The math matches the intuition.

## Code

C++ first. A self-contained Lloyd relaxation in 2D using a coarse grid to integrate the Voronoi cells (no external Voronoi library needed, so it compiles with Eigen only). It shows the exact math symbols: cells by nearest seed, centroid $c_i$, mass $m_i$, gradient $2m_i(x_i-c_i)$, and the monotone energy drop.

```cpp
// lloyd_cvt.cpp  -- Centroidal Voronoi Tessellation via Lloyd relaxation.
// Build: g++ -O2 -I/path/to/eigen lloyd_cvt.cpp -o lloyd_cvt
#include <Eigen/Dense>
#include <vector>
#include <iostream>
#include <limits>

using Vec2 = Eigen::Vector2d;

int main() {
    // Domain [0,1]^2 sampled on an NxN grid acts as the integration measure (rho = 1).
    const int N = 200;
    const int k = 8;                      // number of generators (seeds)
    std::vector<Vec2> x(k);               // generators x_i
    for (int i = 0; i < k; ++i)           // deterministic scatter
        x[i] = Vec2(0.1 + 0.8 * ((i * 7 % k) / double(k)),
                    0.1 + 0.8 * ((i * 3 % k) / double(k)));

    auto energy_and_step = [&](bool apply) {
        std::vector<Vec2> sum(k, Vec2::Zero());   // accumulates y over cell -> for c_i
        std::vector<double> mass(k, 0.0);         // m_i = area of V_i (count of grid cells)
        double F = 0.0;                           // quantization energy
        const double cell = 1.0 / (N * N);        // dy (uniform density)
        for (int a = 0; a < N; ++a)
            for (int b = 0; b < N; ++b) {
                Vec2 y((a + 0.5) / N, (b + 0.5) / N);
                int best = 0; double bd2 = std::numeric_limits<double>::max();
                for (int i = 0; i < k; ++i) {     // nearest seed => Voronoi region V_i
                    double d2 = (y - x[i]).squaredNorm();
                    if (d2 < bd2) { bd2 = d2; best = i; }
                }
                sum[best] += y; mass[best] += 1.0; // integrate y and 1 over V_best
                F += bd2 * cell;                   // ||y - x_i||^2 * dy
            }
        if (apply)
            for (int i = 0; i < k; ++i)
                if (mass[i] > 0) x[i] = sum[i] / mass[i];  // x_i <- c_i  (centroid step)
        return F;
    };

    double prev = energy_and_step(false);
    std::cout << "F0 = " << prev << "\n";
    for (int it = 0; it < 50; ++it) {
        energy_and_step(true);                 // Lloyd: move every seed to its centroid
        double F = energy_and_step(false);     // re-measure energy after the move
        std::cout << "it " << it << "  F = " << F << "\n";
        if (prev - F < 1e-9) break;            // converged: x_i ~= c_i, gradient ~ 0
        prev = F;
    }
    for (int i = 0; i < k; ++i)
        std::cout << "seed " << i << ": " << x[i].transpose() << "\n";
    return 0;
}
```

The energy printed must decrease every iteration. That is your unit test: if it ever rises, your centroid or cell assignment is wrong.

A short Python sketch of the **segmentation** gate (region growing over a mesh by dihedral angle), which clarifies the routing logic. Assumes `trimesh`.

```python
import numpy as np, trimesh
from collections import deque

def segment_by_dihedral(mesh, crease_deg=30.0):
    fn = mesh.face_normals                         # n_f per face
    adj = mesh.face_adjacency                      # pairs (f, g) sharing an edge
    cos_thr = np.cos(np.radians(crease_deg))
    # feature edge if angle between face normals exceeds the crease threshold
    cosang = np.einsum('ij,ij->i', fn[adj[:,0]], fn[adj[:,1]])
    smooth = cosang >= cos_thr                     # theta_e <= crease  => not a feature edge
    # build face graph keeping only non-feature edges
    nbr = {i: [] for i in range(len(mesh.faces))}
    for (f, g), ok in zip(adj, smooth):
        if ok:
            nbr[f].append(g); nbr[g].append(f)
    label = -np.ones(len(mesh.faces), dtype=int); cur = 0
    for s in range(len(mesh.faces)):               # flood fill within crease fences
        if label[s] != -1: continue
        q = deque([s]); label[s] = cur
        while q:
            f = q.popleft()
            for g in nbr[f]:
                if label[g] == -1:
                    label[g] = cur; q.append(g)
        cur += 1
    return label                                   # one integer region id per face

# Routing: a low-curvature, planar-consensus region -> CADFit/NURBS; else freeform mesh.
```

## Connections

- [[discrete-curvature]] supplies the mean curvature $H$, principal curvatures $\kappa_{\max}$, and dihedral angles that drive both the sizing field and the segmentation fences. Remeshing is only as good as the curvature estimate feeding it.
- [[sdf-signed-heat]] is the upstream defect-tolerant cleanup: you reconstruct a clean watertight surface from broken scan data first, then remesh it. Remeshing a soup with holes produces garbage; remesh the signed-distance reconstruction instead.
- [[robust-ransac]] is the per-segment primitive fitter. After segmentation labels a region "smooth", RANSAC decides whether it is truly a plane or cylinder and returns the parameters that become CAD tokens.
- [[nurbs-surfaces-fitting]] consumes a remeshed, segmented smooth patch as an ordered, well-conditioned sample set; CVT-regular triangles are what make a stable tensor-product fit possible.
- [[cadfit-program-synthesis]] is the ultimate consumer of the segmentation decision: only the smooth, primitive-consensus regions are handed to program synthesis; the rough regions are excluded so the synthesizer is not asked to fit a parametric extrusion to a fracture surface.

## Pitfalls and failure modes

- **Normal frame inconsistency.** The 6D lifting needs a **consistently oriented** normal field. If half the mesh has flipped normals, the lifted points scatter and CVT produces shredded triangles. Orient the mesh (consistent winding) before lifting, and smooth the normals.
- **Un-smoothed normals on reconstructed scans.** Raw scan normals are noisy. Feed them into anisotropy and you get speckled, over-refined triangles chasing noise. Always Laplacian-smooth normals first; this is why Geogram's workflow has an explicit smoothing step.
- **Lloyd local minima and slow tail.** Lloyd converges to a critical point, not the global minimum, and crawls near convergence. Bad seeding gives lopsided cells. Seed by a quasi-random or curvature-weighted sample, and switch to L-BFGS on the $C^2$ energy for the last mile.
- **Density-curvature feedback.** Driving density by curvature on a noisy mesh amplifies noise into a wildly non-uniform sizing field. Clamp $\rho$ to $[\rho_{\min},\rho_{\max}]$ and bound the **gradation** (ratio of adjacent edge lengths) so triangle size changes smoothly.
- **Degenerate restricted cells.** On thin plates or near self-intersections, one Voronoi cell can split into disconnected surface patches; its centroid then lands off the surface. Use a localized RVD (LRVD) or reproject centroids onto $S$.
- **Wrong crease threshold = wrong routing.** Set $\theta_{\text{crease}}$ too low and a single rough boulder fragments into hundreds of tiny regions; too high and a sawn face merges with the rock behind it. The threshold is the single most sensitive knob; calibrate it on labelled stone, not by guessing.
- **Anisotropy is not commutative with fitting.** Stretched anisotropic triangles are great for representation but can bias a naive least-squares plane fit (samples cluster along the stretch). Fit on area-weighted or uniform samples, not raw vertices.
- **Tolerance scale.** A crease angle test is unit-free, but the curvature roughness threshold $\tau$ and residual $\varepsilon$ are not. Scale them by object size $L$ or your mm-vs-m mismatch silently mislabels every face.

## Practice

1. **Lloyd energy monitor.** Run the C++ CVT above with $k=8$ on $[0,1]^2$. Print $F$ each iteration and assert it is monotonically non-increasing. Then deliberately break it: in the centroid step use the cell's bounding-box center instead of the true centroid, and observe $F$ rising. Explain why.

2. **Gradient verification.** For the 1D worked example, perturb $x_2$ by $\pm 0.01$ and estimate $\partial F/\partial x_2$ by finite difference. Confirm it matches $2 m_2 (x_2 - c_2)$ to three digits at $x_2=0.6$.

3. **Curvature sizing field.** Take a mesh with a flat region and a tight fillet. Compute $\kappa_{\max}$ per vertex, build $\rho \propto \kappa_{\max}^{4/5}$, and color the mesh by target edge length $h\propto \rho^{-1/2}$. Verify visually that the fillet gets short edges and the flat gets long ones.

4. **6D lifting experiment.** Lift a unit cube's surface to 6D with $s=0,\,0.5,\,2.0$. For each, measure the 6D distance between two adjacent vertices on the same face versus across an edge. Show the across-edge distance grows with $s$ while the same-face distance does not. Relate this to triangle stretching.

5. **Segmentation gate.** Run the Python region-grower on a scanned slab mesh at $\theta_{\text{crease}}=15^\circ,\,30^\circ,\,45^\circ$. Count regions at each threshold and overlay the labels. Identify the threshold that cleanly separates the sawn top from the split sides, and justify routing each region to NURBS vs freeform.

6. **End-to-end (stretch).** Chain it: signed-distance reconstruct a noisy scan, anisotropically remesh in Geogram, segment by dihedral angle, RANSAC-fit a plane to the largest smooth region, and report the plane normal and residual. This is the exact front half of the CADFit feed.

## References

- Du, Faber, Gunzburger. "Centroidal Voronoi Tessellations: Applications and Algorithms." SIAM Review 41(4), 1999, pp. 637-676. https://epubs.siam.org/doi/10.1137/S0036144599352836
- Liu, Wang, Lévy, Sun, Yan, Lu, Yang. "On Centroidal Voronoi Tessellation -- Energy Smoothness and Fast Computation." ACM Transactions on Graphics 28(4), 2009. https://www.microsoft.com/en-us/research/wp-content/uploads/2016/12/On-Centroidal-Voronoi-Tessellation-Energy-Smoothness-and-Fast-Computation.pdf
- Lévy, Bonneel. "Variational Anisotropic Surface Meshing with Voronoi Parallel Linear Enumeration." Proc. 21st International Meshing Roundtable, 2013. https://members.loria.fr/Bruno.Levy/papers/CVT_IMR_2012.pdf
- Lévy, "Anisotropic and feature sensitive triangular remeshing using normal lifting." Journal of Computational and Applied Mathematics, 2015. https://hal.science/hal-01202738
- Geogram documentation and Remeshing wiki (Bruno Lévy). https://github.com/BrunoLevy/geogram/wiki/Remeshing and project site http://alice.loria.fr/software/geogram/doc/html/index.html
- Keenan Crane. "Discrete Differential Geometry: An Applied Introduction" (CMU course notes, curvature and Laplacian chapters). https://www.cs.cmu.edu/~kmcrane/Projects/DDG/
- Feng, Crane. "A Heat Method for Generalized Signed Distance." ACM Transactions on Graphics, 2024. https://www.cs.cmu.edu/~kmcrane/Projects/SignedHeatMethod/
- Yan, Lévy, Liu, Sun, Wang. "Isotropic Remeshing with Fast and Exact Computation of Restricted Voronoi Diagram." Computer Graphics Forum (SGP), 2009. https://www.microsoft.com/en-us/research/wp-content/uploads/2016/12/Isotropic-Remeshing-with-Fast-and-Exact-Computation-of-Restricted-Voronoi-Diagram.pdf
- Nehme et al. "CADFit: Precise Mesh-to-CAD Program Generation with Hybrid Optimization." arXiv:2605.01171. https://arxiv.org/abs/2605.01171 and code https://github.com/ghadinehme/CADFit
- OpenCASCADE reference manual, `GeomAPI_PointsToBSplineSurface` (downstream NURBS fitting). https://dev.opencascade.org/doc/refman/html/class_geom_a_p_i___points_to_b_spline_surface.html

## Flashcard seeds

- Q: Derive the CVT quantization energy $F(X)$ (step 1 of 2). Start from the per-region density-weighted squared error of points $y$ from their seed $x_i$, namely $\int_{V_i}\rho(y)\lVert y-x_i\rVert^2\,dy$. <br> <b>Formula:</b> $F(X)=\sum_{i=1}^{k}\int_{V_i}\rho(y)\lVert y-x_i\rVert^2\,dy$ <br> <i>Hint:</i> sum the one-region error over all $k$ Voronoi regions; each point belongs to exactly its nearest seed's region. :: A: One region contributes $\int_{V_i}\rho(y)\lVert y-x_i\rVert^2\,dy$ <br> The regions $V_i$ tile the domain, so sum over $i=1,\dots,k$ <br> $\boxed{F(X)=\sum_{i=1}^{k}\int_{V_i}\rho(y)\lVert y-x_i\rVert^2\,dy}$
- Q: Derive the CVT energy gradient (step 2 of 2). Start from $F=\sum_i\int_{V_i}\rho\lVert y-x_i\rVert^2\,dy$; cell-boundary motion contributes nothing. <br> <b>Formula:</b> $\dfrac{\partial F}{\partial x_i}=2\,m_i\,(x_i-c_i)$ with $m_i=\int_{V_i}\rho\,dy,\ c_i=\frac1{m_i}\int_{V_i}\rho\,y\,dy$ <br> <i>Hint:</i> differentiate only the integrand using $\partial_{x_i}\lVert y-x_i\rVert^2=-2(y-x_i)$, then split the integral into $m_i$ and $c_i$. :: A: Boundary terms cancel (each cell wall is equidistant from two seeds), so $\partial F/\partial x_i=\int_{V_i}\rho\,\partial_{x_i}\lVert y-x_i\rVert^2\,dy=-2\int_{V_i}\rho(y)(y-x_i)\,dy$ <br> Split: $-2\big(\int_{V_i}\rho\,y\,dy-x_i\,m_i\big)$ <br> Substitute $\int_{V_i}\rho\,y\,dy=m_i c_i$: $-2(m_i c_i-m_i x_i)$ <br> $\boxed{\dfrac{\partial F}{\partial x_i}=2\,m_i\,(x_i-c_i)}$
- Q: From the gradient $\partial F/\partial x_i=2m_i(x_i-c_i)$, derive the centroidal (stationarity) condition. <br> <b>Formula:</b> set $\nabla F=0$; result $x_i=c_i\ \forall i$ <br> <i>Hint:</i> a non-empty cell has mass $m_i>0$, so decide which factor of the product must be zero. :: A: A critical point needs $2m_i(x_i-c_i)=0$ for every $i$ <br> Mass $m_i>0$ for any non-empty cell, so the vanishing factor is $(x_i-c_i)$ <br> $\boxed{x_i=c_i\ \ \forall i}$ — every seed sits at its region's centroid
- Q: Derive the centroid $c_i$ that minimizes a single cell's variance $g(x_i)=\int_{V_i}\rho\lVert y-x_i\rVert^2\,dy$. <br> <b>Formula:</b> $c_i=\dfrac{1}{m_i}\int_{V_i}y\,\rho(y)\,dy,\quad m_i=\int_{V_i}\rho\,dy$ <br> <i>Hint:</i> set $g'(x_i)=0$ using $g'(x_i)=-2\int_{V_i}\rho(y-x_i)\,dy$, then solve for $x_i$. :: A: Set $g'(x_i)=-2\int_{V_i}\rho(y-x_i)\,dy=0$ <br> So $\int_{V_i}\rho\,y\,dy=x_i\int_{V_i}\rho\,dy=x_i\,m_i$ <br> Divide by $m_i$: $x_i=\frac1{m_i}\int_{V_i}\rho\,y\,dy$ <br> $\boxed{c_i=\dfrac{1}{m_i}\int_{V_i}y\,\rho(y)\,dy}$
- Q: Find the Voronoi boundary $b$ between two 1D seeds $x_1=0.2,\ x_2=0.6$. <br> <b>Formula:</b> $b=(x_1+x_2)/2$ (point equidistant from both seeds) <br> <i>Hint:</i> set $\lvert b-x_1\rvert=\lvert b-x_2\rvert$ and solve — it is just the midpoint. :: A: Equidistance gives $b=(x_1+x_2)/2$ <br> $b=(0.2+0.6)/2=0.4$, so $V_1=[0,0.4],\ V_2=[0.4,1]$ <br> $\boxed{b=0.4}$
- Q: Worked example (1D, $\rho\equiv1$). Compute both centroids of cells $V_1=[0,0.4],\ V_2=[0.4,1]$. <br> <b>Formula:</b> uniform density $\Rightarrow c_i=(\text{lo}+\text{hi})/2$ (interval midpoint) <br> <i>Hint:</i> just average the two endpoints of each interval. :: A: $c_1=(0+0.4)/2=0.2$ <br> $c_2=(0.4+1)/2=0.7$ <br> $\boxed{c_1=0.2,\ c_2=0.7}$
- Q: Worked example. With $m_2=0.6,\ x_2=0.6,\ c_2=0.7$, evaluate $\partial F/\partial x_2$ and state which way $x_2$ moves. <br> <b>Formula:</b> $\partial F/\partial x_2=2m_2(x_2-c_2)$ <br> <i>Hint:</i> plug in the three numbers; the sign of the result tells the descent direction (move opposite the gradient). :: A: $\partial F/\partial x_2=2(0.6)(0.6-0.7)=1.2\times(-0.1)$ <br> $\boxed{\partial F/\partial x_2=-0.12}$ <br> Negative gradient $\Rightarrow$ $x_2$ increases (moves right toward $0.7$)
- Q: Closing test (close the loop). One Lloyd step sends $x_2\leftarrow c_2=0.7$. Verify the gradient at the frozen cell is now zero, and that a finite-difference estimate at $x_2=0.6$ matched $-0.12$. <br> <b>Formula:</b> $\partial F/\partial x_2=2m_2(x_2-c_2)$; FD: $\partial F/\partial x_2\approx[F(0.61)-F(0.59)]/0.02$ <br> <i>Hint:</i> a pass means the analytic value is $0$ after the move and the FD slope agrees with $2m_2(x_2-c_2)$ to ~3 digits. :: A: After the move $x_2=c_2=0.7$, so $2m_2(c_2-c_2)=0$ — the descent fixed point <br> FD at $x_2=0.6$: $[F(0.61)-F(0.59)]/0.02\approx-0.12$, matching $2(0.6)(0.6-0.7)=-0.12$ <br> $\boxed{\text{grad}=0\text{ after move; FD }\approx-0.12\ \checkmark}$
- Q: Worked example (energy before a step). Compute $F$ for seeds $0.2,0.6$, cells $[0,0.4],[0.4,1]$, $\rho\equiv1$. <br> <b>Formula:</b> $F=\int_0^{0.4}(y-0.2)^2dy+\int_{0.4}^{1}(y-0.6)^2dy$, with $\int(y-a)^2dy=(y-a)^3/3$ <br> <i>Hint:</i> integrate each cell separately as $[(y-a)^3/3]$ over its endpoints, then add. :: A: First $=[(y-0.2)^3/3]_0^{0.4}=(0.2^3-(-0.2)^3)/3=(0.008+0.008)/3\approx0.005333$ <br> Second $=[(y-0.6)^3/3]_{0.4}^{1}=(0.4^3-(-0.2)^3)/3=(0.064+0.008)/3=0.024$ <br> $\boxed{F\approx0.029333}$
- Q: Closing test (energy must fall). After Lloyd ($x_2\to0.7$, cells $[0,0.45],[0.45,1]$) recompute $F$ and confirm Lloyd's monotone-decrease guarantee. <br> <b>Formula:</b> $F=\int_0^{0.45}(y-0.2)^2dy+\int_{0.45}^{1}(y-0.7)^2dy$; pass if $F<0.029333$ <br> <i>Hint:</i> integrate each cell as $[(y-a)^3/3]$, add, then compare to the pre-step $F\approx0.029333$. :: A: First $=[(y-0.2)^3/3]_0^{0.45}=(0.25^3+0.2^3)/3\approx0.007875$ <br> Second $=[(y-0.7)^3/3]_{0.45}^{1}=(0.3^3+0.25^3)/3\approx0.014208$ <br> $F\approx0.022083<0.029333$ <br> $\boxed{\Delta F\approx-0.00725<0}$ — energy fell, as Lloyd guarantees
- Q: Derive the squared 6D distance under the normal-lifting map (step 1 of 2). Start from $\Phi(x)=(x,\,sN(x))$. <br> <b>Formula:</b> $\lVert\Phi(x)-\Phi(x')\rVert^2=\lVert x-x'\rVert^2+s^2\lVert N(x)-N(x')\rVert^2$ <br> <i>Hint:</i> the difference $\Phi(x)-\Phi(x')$ is two stacked blocks; the squared norm of stacked blocks is the sum of each block's squared norm. :: A: $\Phi(x)-\Phi(x')=\big(x-x',\ s(N(x)-N(x'))\big)$ <br> $\lVert(a,b)\rVert^2=\lVert a\rVert^2+\lVert b\rVert^2$, and the second block carries the factor $s$ <br> $\boxed{\lVert\Phi(x)-\Phi(x')\rVert^2=\lVert x-x'\rVert^2+s^2\lVert N(x)-N(x')\rVert^2}$
- Q: From the 6D squared distance (step 2 of 2), explain why flat regions get big triangles and creases get small ones. <br> <b>Formula:</b> $\lVert\Phi-\Phi'\rVert^2=\lVert x-x'\rVert^2+s^2\lVert N-N'\rVert^2$ <br> <i>Hint:</i> ask what the normal-penalty term does when $N\approx N'$ (flat) versus when $N$ swings fast (crease). :: A: Flat patch: $N(x)\approx N(x')$, so $s^2\lVert N-N'\rVert^2\approx0$; CVT may space seeds far apart $\Rightarrow$ large triangles <br> Crease: $\lVert N-N'\rVert$ large, penalty blows up; CVT packs seeds tightly $\Rightarrow$ small triangles <br> $\boxed{\text{seed spacing}\downarrow\ \text{as}\ s^2\lVert N-N'\rVert^2\uparrow}$ — anisotropy is curvature adaptivity in disguise
- Q: Worked example (anisotropy, mm-scale crease). Across a $90^\circ$ crease with unit normals and $s=0.5$, compute the extra 6D squared distance added to the 3D term. <br> <b>Formula:</b> extra $=s^2\lVert N-N'\rVert^2$, with $\lVert N-N'\rVert^2=2-2\,N\cdot N'$ for unit normals <br> <i>Hint:</i> first compute $N\cdot N'$ for perpendicular normals, then plug into $\lVert N-N'\rVert^2$, then multiply by $s^2$. :: A: Perpendicular unit normals: $N\cdot N'=0$, so $\lVert N-N'\rVert^2=2-0=2$ <br> Extra $=s^2\cdot2=0.25\times2$ <br> $\boxed{=0.5}$ added to the 3D squared distance, pushing seeds apart across the crease
- Q: Worked example (anisotropy, larger $s$). Same $90^\circ$ crease but $s=2.0$. Compute the extra term and compare to the $s=0.5$ value ($0.5$). <br> <b>Formula:</b> extra $=s^2\lVert N-N'\rVert^2$, with $\lVert N-N'\rVert^2=2$ here <br> <i>Hint:</i> the geometry (hence $\lVert N-N'\rVert^2$) is unchanged; only $s^2$ changes, so take the ratio of $s^2$ values. :: A: $\lVert N-N'\rVert^2=2$ (geometry unchanged) <br> Extra $=s^2\cdot2=4\times2=8$ <br> Ratio to $s=0.5$: $8/0.5=16$ (also $s^2$ ratio $2^2/0.5^2=16$) <br> $\boxed{=8,\ \ 16\times\text{ the }s=0.5\text{ value}}$ — stronger $s$ forces far more triangles across the same crease
- Q: From the sizing field $\rho\propto\kappa_{\max}^{4/5}$, derive how target edge length $h$ scales with curvature on a surface ($\dim=2$). <br> <b>Formula:</b> $h\propto\rho^{-1/\dim}=\rho^{-1/2}$, and $\rho\propto\kappa_{\max}^{4/5}$ <br> <i>Hint:</i> substitute the curvature law into $h\propto\rho^{-1/2}$ and multiply the exponents. :: A: $h\propto\rho^{-1/2}$ <br> Substitute $\rho\propto\kappa_{\max}^{4/5}$: $h\propto(\kappa_{\max}^{4/5})^{-1/2}=\kappa_{\max}^{-2/5}$ <br> $\boxed{h\propto\kappa_{\max}^{-2/5}}$ — flat ($\kappa_{\max}\to0$) gets long edges, creases get short ones
- Q: Worked example (sizing field). A fillet has $\kappa_{\max}$ that is $32\times$ the flat region's. By what factor does target edge length $h$ shrink? <br> <b>Formula:</b> $h\propto\kappa_{\max}^{-2/5}$, so $h_{\text{fillet}}/h_{\text{flat}}=32^{-2/5}$ <br> <i>Hint:</i> write $32=2^5$ to make the fractional power clean. :: A: $h_{\text{fillet}}/h_{\text{flat}}=(32)^{-2/5}$ <br> $32=2^5$, so $32^{-2/5}=2^{5\cdot(-2/5)}=2^{-2}=1/4$ <br> $\boxed{h\ \text{shrinks to }1/4}$ on the fillet relative to the flat
- Q: Derive the feature-edge test from two adjacent unit face normals $n_f,n_g$. <br> <b>Formula:</b> $\theta_e=\arccos(n_f\cdot n_g)$; feature edge $\iff\theta_e>\theta_{\text{crease}}$ <br> <i>Hint:</i> the dot product of unit normals is the cosine of the dihedral angle; a crease is where that angle exceeds the threshold. :: A: For unit normals, $n_f\cdot n_g=\cos\theta_e$, so $\theta_e=\arccos(n_f\cdot n_g)\in[0,\pi]$ <br> A crease is a sharp deviation, i.e. $\theta_e$ above a threshold $\theta_{\text{crease}}$ <br> $\boxed{\text{feature edge}\iff\theta_e=\arccos(n_f\cdot n_g)>\theta_{\text{crease}}}$
- Q: Worked example (segmentation). Two faces have unit normals with $n_f\cdot n_g=0.5$. Is the edge a feature edge at $\theta_{\text{crease}}=30^\circ$? <br> <b>Formula:</b> $\theta_e=\arccos(n_f\cdot n_g)$; feature edge if $\theta_e>\theta_{\text{crease}}$ <br> <i>Hint:</i> compute $\arccos(0.5)$ in degrees first, then compare to $30^\circ$. :: A: $\theta_e=\arccos(0.5)=60^\circ$ <br> $60^\circ>30^\circ$ <br> $\boxed{\text{Yes — feature edge}}$, so region growing stops here (crease fence)
- Q: Closing test (segmentation, flat-face limit). Two nearly coplanar faces have $n_f\cdot n_g=0.99$. Classify the edge at $\theta_{\text{crease}}=30^\circ$ and check the perfectly-flat limit. <br> <b>Formula:</b> $\theta_e=\arccos(n_f\cdot n_g)$; not a feature edge if $\theta_e\le\theta_{\text{crease}}$ <br> <i>Hint:</i> a pass means $0.99$ gives $\theta_e$ well under $30^\circ$ (crossed by region growing) and $n_f\cdot n_g=1$ gives $\theta_e=0$ (always merged). :: A: $\theta_e=\arccos(0.99)\approx8.1^\circ$ <br> $8.1^\circ\le30^\circ\Rightarrow$ NOT a feature edge — region growing crosses it (same smooth face) <br> Flat limit: $n_f\cdot n_g=1\Rightarrow\theta_e=\arccos(1)=0^\circ\le\theta_{\text{crease}}$, always merged <br> $\boxed{\text{not a feature edge; flat limit }\theta_e\to0\ \checkmark}$
- CLOZE: Lloyd's algorithm is gradient descent on the CVT energy with a fixed full step, because moving each seed to its centroid sets {{c1::$x_i \leftarrow c_i$}} and never increases $F$.
- CLOZE: In normal lifting, large triangles form where normals barely change ({{c1::flat/low-curvature}} regions) and small triangles form where normals swing fast (sharp features).
- CLOZE: Before anisotropic remeshing on a reconstructed scan you must {{c1::Laplacian-smooth the normals}}, or refinement chases normal noise. <br> <i>hint:</i> raw scan normals are noisy.
- CLOZE: In Geogram the effective anisotropy scale is the user value times {{c1::$0.02$}}, via `set_anisotropy(mesh, anisotropy*0.02)`.
- Q: When-to-use: which regions does the segmentation gate route to NURBS/CADFit vs freeform? <br> <i>Hint:</i> the gate keys on curvature plus whether a RANSAC primitive reaches consensus under $\varepsilon$. :: A: Low-curvature regions with a RANSAC plane/cylinder consensus (residual under $\varepsilon$) $\to$ analytic primitive or NURBS fitting <br> High-curvature regions with no consensus primitive $\to$ freeform mesh or tessellated IFC fallback <br> $\boxed{\text{smooth+consensus}\to\text{CAD};\ \text{rough/no-fit}\to\text{freeform}}$
- Q: Why prefer L-BFGS over plain Lloyd for CVT? <br> <i>Hint:</i> recall the smoothness class of the CVT energy on convex domains. :: A: The CVT energy is $C^2$ on convex domains (Liu et al. 2009) <br> Quasi-Newton (L-BFGS) uses that curvature to converge faster than Lloyd's slow linear tail near the minimum <br> $\boxed{C^2\ \text{energy}\Rightarrow\text{L-BFGS beats Lloyd's slow tail}}$
- Q: Pitfall: what does a localized RVD (LRVD) fix, and why does it matter? <br> <i>Hint:</i> think about a thin plate where one cell touches the surface in two disconnected patches. :: A: On thin plates or near self-intersections a single Voronoi cell can split into disconnected surface patches; its centroid then lands off the surface <br> LRVD (or reprojecting centroids onto $S$) keeps seeds on the surface <br> $\boxed{\text{LRVD/reproject}\Rightarrow\text{centroids stay on }S}$
- Q: Pitfall: why is $\theta_{\text{crease}}$ the most sensitive knob in stone segmentation? <br> <i>Hint:</i> consider what happens to region count as the threshold goes too low, then too high. :: A: Too low fragments one rough boulder face into hundreds of regions; too high merges a sawn face with the natural rock behind it <br> It is unit-free, so calibrate on labelled stone, not by guessing <br> $\boxed{\text{mis-set }\theta_{\text{crease}}\Rightarrow\text{wrong freeform-vs-CAD routing}}$
