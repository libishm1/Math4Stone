# Rigid registration: Kabsch-Umeyama and ICP

Prereqs: [[03_eig-svd-pca]], [[linear-least-squares]] => this => [[sdf-signed-heat]]

## Why it matters for the stone scan->CAD->IFC / CADFit workflow

When you scan an irregular block from two stations, or scan it today and re-scan it tomorrow on the cutting bed, the two clouds live in different coordinate frames. Registration finds the single rigid motion that snaps one onto the other so you can merge partial scans into one watertight cloud before reconstruction. The same machinery aligns a finished scan to a CAD template or to the machine frame, which is the step that lets a robot or saw act on the geometry you designed. Inside the pipeline the order is fixed: a coarse global alignment (or a known marker frame) first, then Kabsch-Umeyama on matched features, then ICP to polish the fit to sub-millimeter, and only then signed-distance reconstruction and B-spline fitting. A CADFit-style synthesizer also needs its primitives expressed in a canonical, well-registered frame, otherwise every fitted parameter inherits the misalignment as error.

## Intuition

You have two copies of the same rocky shape, drawn on two sheets of tracing paper, in different positions and rotations. You want to slide and spin one sheet until it lies exactly on top of the other. If you already know which dot on sheet A is which dot on sheet B (the correspondences), there is a one-shot closed-form answer: line up the two centers of mass, then find the single rotation that best aligns the leftover offsets. That one-shot answer is Kabsch-Umeyama, and it is just one SVD.

The catch is the real world rarely hands you the correspondences. So ICP guesses: for each point on A, call its nearest point on B the match, solve the one-shot rotation, move A a little, then re-guess the matches and repeat. Each round the matches get better and the fit tightens. It is hill-climbing. Analogy: parking a car by feel. You glance at the curb (correspond), straighten the wheel a bit (solve and move), glance again, straighten again. It works beautifully when you started roughly in the spot. It fails when you started facing the wrong way, because every glance only tells you the *local* correction.

## The math (first principles)

### Conventions

Points are column vectors in $\mathbb{R}^3$. Frames are right-handed. A proper rotation $R \in SO(3)$ satisfies $R^TR = I$ and $\det R = +1$; if $\det R = -1$ it is a reflection, which is a flipped, mirror-image fit and a bug here. A rigid transform is the pair $(R, t)$ acting as $x \mapsto Rx + t$, or the homogeneous $4\times4$ matrix $T = \begin{bmatrix} R & t \\ 0 & 1\end{bmatrix} \in SE(3)$. Source points $\{p_i\}$ are moved onto target points $\{q_i\}$. Units are whatever the scan uses; keep one unit per computation (meters at site scale, millimeters at shop scale). Tolerances are relative to object size: a registration "converged" when the mean residual change falls below $\varepsilon \cdot L$ with $L$ the object diameter, never a fixed absolute number across scales.

### The point-to-point objective

For known correspondences, rigid registration minimizes squared Euclidean residuals:

$$\min_{R \in SO(3),\, t \in \mathbb{R}^3} \; \sum_{i=1}^{N} \big\| R p_i + t - q_i \big\|^2.$$

This is the **orthogonal Procrustes** problem with a translation. It is *not* an ordinary linear least squares, because $R$ is constrained to the curved manifold $SO(3)$. The trick is that the constraint has a closed-form solution through the SVD.

### Step 1: translation separates out

Take the gradient with respect to $t$ and set it to zero. With $R$ fixed,

$$\frac{\partial}{\partial t} \sum_i \| R p_i + t - q_i \|^2 = 2 \sum_i (R p_i + t - q_i) = 0 \;\Longrightarrow\; t = \bar q - R \bar p,$$

where the centroids are

$$\bar p = \frac{1}{N}\sum_{i=1}^N p_i, \qquad \bar q = \frac{1}{N}\sum_{i=1}^N q_i.$$

So the optimal translation always maps the source centroid onto the target centroid. Substitute it back. Define the centered points

$$\tilde p_i = p_i - \bar p, \qquad \tilde q_i = q_i - \bar q.$$

The objective collapses to a pure rotation problem on the centered clouds:

$$\min_{R \in SO(3)} \; \sum_i \big\| R \tilde p_i - \tilde q_i \big\|^2.$$

### Step 2: rotation becomes a trace maximization

Expand one term: $\| R\tilde p_i - \tilde q_i \|^2 = \| \tilde p_i\|^2 + \|\tilde q_i\|^2 - 2\, \tilde q_i^T R \tilde p_i$, using $\|R\tilde p_i\| = \|\tilde p_i\|$ because $R$ is orthogonal. The first two terms do not depend on $R$. So minimizing the sum is the same as **maximizing** $\sum_i \tilde q_i^T R \tilde p_i$. Rewrite that with the trace identity $a^T b = \operatorname{tr}(b a^T)$:

$$\sum_i \tilde q_i^T R \tilde p_i = \operatorname{tr}\!\Big( R \sum_i \tilde p_i \tilde q_i^T \Big) = \operatorname{tr}(R H), \qquad H = \sum_{i=1}^N \tilde p_i \tilde q_i^T.$$

$H \in \mathbb{R}^{3\times3}$ is the **cross-covariance** of the two centered clouds. Note the convention: with $H = \sum \tilde p_i \tilde q_i^T$ (source on the left, target on the right) the maximization is $\max_R \operatorname{tr}(RH)$.

### Step 3: SVD gives the rotation (Kabsch-Umeyama)

Take the SVD $H = U \Sigma V^T$ with $\Sigma = \operatorname{diag}(\sigma_1 \ge \sigma_2 \ge \sigma_3 \ge 0)$. Then $\operatorname{tr}(RH) = \operatorname{tr}(R U \Sigma V^T) = \operatorname{tr}(\Sigma V^T R U)$. The matrix $M = V^T R U$ is orthogonal, so every entry satisfies $|M_{jj}| \le 1$, and $\operatorname{tr}(\Sigma M) = \sum_j \sigma_j M_{jj} \le \sum_j \sigma_j$. The bound is hit when $M = I$, i.e. $V^T R U = I$, giving $R = V U^T$.

That maximizes the trace over all orthogonal matrices, but $V U^T$ might have $\det = -1$ (a reflection). To stay in $SO(3)$, flip the sign of the term tied to the smallest singular value. The **Umeyama correction** is the determinant term:

$$\boxed{\;R = V \begin{bmatrix} 1 & 0 & 0 \\ 0 & 1 & 0 \\ 0 & 0 & d \end{bmatrix} U^T, \qquad d = \det(V U^T), \qquad t = \bar q - R \bar p.\;}$$

Here $d = \pm 1$. When the data is clean $d = +1$ and the bracket is the identity. When the best orthogonal fit is a reflection $d = -1$, and flipping the last column costs only $2\sigma_3$ (the smallest singular value), which is the cheapest possible repair. This is exactly what Arun et al. 1987 missed and Umeyama 1991 fixed: without the $d$ term, a noisy or near-degenerate cloud can return a mirror image. Umeyama also gives the optimal scale $c = \frac{1}{\sigma_p^2}\operatorname{tr}(\Sigma\, \mathrm{diag}(1,1,d))$ with $\sigma_p^2 = \frac1N\sum_i\|\tilde p_i\|^2$ if you want a similarity transform; for rigid registration set $c = 1$.

Degeneracy guard: if two singular values are equal (e.g. points all on a line or a sphere of symmetry), the rotation about that axis is unconstrained and $R$ is not unique. Check $\sigma_2 - \sigma_3$ relative to $\sigma_1$ before trusting the result.

### The ICP loop (unknown correspondences)

When correspondences are unknown, alternate:

1. **Correspond.** For each transformed source point $T p_i$, find the nearest target point $q_{c(i)} = \arg\min_q \| T p_i - q \|$. Use a kd-tree; brute force is $O(NM)$.
2. **Solve.** Run Kabsch-Umeyama on the matched pairs $\{(p_i, q_{c(i)})\}$ to get an incremental $\Delta T$.
3. **Update.** Compose $T \leftarrow \Delta T \cdot T$ and move the source.
4. **Repeat** until the mean residual stops dropping.

This is coordinate descent on the joint cost over (correspondences, transform). Each step does not increase the error, so ICP converges, but only to a **local** minimum. The fixed-point it reaches depends entirely on the starting pose (Besl and McKay 1992; Chen and Medioni 1991).

### Point-to-plane variant

Point-to-point pulls each source point toward a fixed target point, which fights the natural sliding of two surfaces along each other. Point-to-plane instead penalizes only the distance **along the target normal** $n_i$, letting the surfaces slide tangentially. It converges in far fewer iterations on structured surfaces (Chen and Medioni 1991; Low 2004):

$$\min_{R,\,t} \; \sum_i \Big( n_i^T \big( R p_i + t - q_i \big) \Big)^2.$$

This is nonlinear in $R$, but for the small per-iteration correction of ICP you linearize. Write the small rotation as $R \approx I + [r]_\times$ where $r = (\alpha, \beta, \gamma)^T$ are small angles and $[r]_\times$ is the skew-symmetric cross-product matrix. The residual becomes, to first order,

$$n_i^T(p_i - q_i) + n_i^T([r]_\times p_i) + n_i^T t = n_i^T(p_i - q_i) + (p_i \times n_i)^T r + n_i^T t,$$

using the identity $n^T([r]_\times p) = n^T(r \times p) = (p \times n)^T r$. Stack the six unknowns $x = (\alpha, \beta, \gamma, t_x, t_y, t_z)^T$. Each correspondence gives one linear equation $a_i^T x = b_i$ with

$$a_i = \begin{bmatrix} p_i \times n_i \\ n_i \end{bmatrix} \in \mathbb{R}^6, \qquad b_i = -\,n_i^T(p_i - q_i).$$

Solve the $6\times6$ normal equations $\big(\textstyle\sum_i a_i a_i^T\big)\, x = \sum_i a_i b_i$ (linear least squares, the [[linear-least-squares]] prereq), then rebuild $R$ from $(\alpha,\beta,\gamma)$ and re-orthonormalize. Because $r$ is small per step, the $I + [r]_\times$ approximation is valid; iterate to refine.

### Initialization and global registration

ICP is a local method. Its basin of convergence is roughly "less than half the feature spacing and within tens of degrees" of the true pose. Beyond that it locks onto a wrong local minimum (two flat faces snapping together backwards is the classic stone failure). The fixes:

- **Good init:** PCA axis alignment, fiducial markers, a known turntable angle, or a feature-based coarse match (FPFH + RANSAC) before ICP.
- **Go-ICP** (Yang et al. 2016): branch-and-bound over the whole motion space $SE(3)$ with provable upper and lower bounds on the L2 error, giving the **globally optimal** rigid registration with no initial guess. It nests a local ICP inside the search to speed it up. It is slower than ICP but it cannot be fooled by initialization. Use it when you have no reliable prior pose.

## Worked example

Align four 2D source points to four target points (a clean rotation + translation, no noise, so we expect an exact fit). Work in 2D so the SVD is a $2\times2$.

Source $p$: $(0,0), (2,0), (2,1), (0,1)$. Target $q$: $(1,1), (1,3), (0,3), (0,1)$.

**Centroids.** $\bar p = (1, 0.5)$, $\bar q = (0.5, 2)$.

**Centered points.**

| $i$ | $\tilde p_i$ | $\tilde q_i$ |
|---|---|---|
| 1 | $(-1,-0.5)$ | $(0.5,-1)$ |
| 2 | $(1,-0.5)$ | $(0.5,1)$ |
| 3 | $(1,0.5)$ | $(-0.5,1)$ |
| 4 | $(-1,0.5)$ | $(-0.5,-1)$ |

**Cross-covariance** $H = \sum_i \tilde p_i \tilde q_i^T$. Each outer product is a $2\times2$:

$$\tilde p_1 \tilde q_1^T = \begin{bmatrix}-1\\-0.5\end{bmatrix}\begin{bmatrix}0.5 & -1\end{bmatrix} = \begin{bmatrix}-0.5 & 1\\ -0.25 & 0.5\end{bmatrix}.$$

Summing all four (the symmetric off-diagonals cancel):

$$H = \begin{bmatrix} -0.5 - 0.5 - 0.5 - 0.5 & \; 1 + 1 + 1 + 1 \\[2pt] -0.25 + 0.25 - 0.25 + 0.25 & \; 0.5 + 0.5 + 0.5 + 0.5 \end{bmatrix} = \begin{bmatrix} -2 & 4 \\ 0 & 2 \end{bmatrix}.$$

Wait, recompute the $(1,2)$ entry carefully: contributions are $(-1)(-1)=1$, $(1)(1)=1$, $(1)(1)=1$, $(-1)(-1)=1$, sum $4$. The $(1,1)$: $(-1)(0.5)=-0.5$ four times $=-2$. The $(2,1)$: $(-0.5)(0.5)=-0.25$, $(-0.5)(0.5)=-0.25$, $(0.5)(-0.5)=-0.25$, $(0.5)(-0.5)=-0.25$, sum $-1$. The $(2,2)$: $(-0.5)(-1)=0.5$ four times $=2$. So

$$H = \begin{bmatrix} -2 & 4 \\ -1 & 2 \end{bmatrix}.$$

**Guess the answer first to check.** The target is the source rotated $90^\circ$ counterclockwise then translated. A $90^\circ$ CCW rotation is $R = \begin{bmatrix}0 & -1\\ 1 & 0\end{bmatrix}$. Verify on $\tilde p_2 = (1,-0.5)$: $R\tilde p_2 = (0.5, 1) = \tilde q_2$. Yes.

**SVD route.** We want $R = V\,\mathrm{diag}(1, d)\,U^T$ with $d = \det(VU^T)$. Rather than grind the $2\times2$ SVD by hand, confirm the known $R$ satisfies the optimality condition $\operatorname{tr}(RH)$ is maximal and $RH$ is symmetric positive semidefinite (the optimality test for Procrustes):

$$RH = \begin{bmatrix}0 & -1\\ 1 & 0\end{bmatrix}\begin{bmatrix}-2 & 4\\ -1 & 2\end{bmatrix} = \begin{bmatrix}1 & -2\\ -2 & 4\end{bmatrix}.$$

That is symmetric, and its eigenvalues are $0$ and $5$ (trace $5$, det $0$), both $\ge 0$, so $RH$ is positive semidefinite. That is exactly the condition the SVD construction guarantees, so $R$ is the optimizer. $\operatorname{tr}(RH) = 5 = \sigma_1 + \sigma_2$ confirms the maximum.

**Translation.** $t = \bar q - R\bar p = (0.5, 2) - \begin{bmatrix}0&-1\\1&0\end{bmatrix}(1, 0.5) = (0.5, 2) - (-0.5, 1) = (1, 1).$

**Check the fit.** $R p_1 + t = R(0,0) + (1,1) = (1,1) = q_1$. Exact, residual zero, as expected for clean data. The determinant term: $\det R = +1$, so $d=+1$ and no reflection repair was needed.

## Code

Idiomatic C++ with Eigen. The core solve is one SVD. Comments map each line to the math symbol.

```cpp
// kabsch_icp.hpp  -- Kabsch-Umeyama rigid solve + naive ICP, header-only.
#pragma once
#include <Eigen/Dense>
#include <vector>
#include <limits>
#include <stdexcept>

namespace reg {

using Vec3 = Eigen::Vector3d;
using Mat3 = Eigen::Matrix3d;
using Mat4 = Eigen::Matrix4d;

// Kabsch-Umeyama: best rigid (R,t) mapping src onto dst, fixed correspondences.
// Solves min_{R in SO(3), t} sum_i || R*p_i + t - q_i ||^2.
inline Mat4 bestRigidTransform(const std::vector<Vec3>& src,   // p_i
                               const std::vector<Vec3>& dst) { // q_i
    const auto N = src.size();
    if (N != dst.size() || N < 3)
        throw std::runtime_error("need >= 3 paired points of equal count");

    // Centroids: p_bar, q_bar.
    Vec3 pbar = Vec3::Zero(), qbar = Vec3::Zero();
    for (size_t i = 0; i < N; ++i) { pbar += src[i]; qbar += dst[i]; }
    pbar /= double(N); qbar /= double(N);

    // Cross-covariance H = sum_i (p_i - p_bar)(q_i - q_bar)^T  (source x target).
    Mat3 H = Mat3::Zero();
    for (size_t i = 0; i < N; ++i)
        H += (src[i] - pbar) * (dst[i] - qbar).transpose();

    // SVD: H = U * Sigma * V^T.
    Eigen::JacobiSVD<Mat3> svd(H, Eigen::ComputeFullU | Eigen::ComputeFullV);
    const Mat3 U = svd.matrixU();
    const Mat3 V = svd.matrixV();

    // Umeyama reflection guard: d = det(V U^T) is +1 or -1.
    Mat3 D = Mat3::Identity();
    D(2, 2) = (V * U.transpose()).determinant();   // flips only the last axis

    // R = V * diag(1,1,d) * U^T ;  t = q_bar - R p_bar.
    const Mat3 R = V * D * U.transpose();
    const Vec3 t = qbar - R * pbar;

    Mat4 T = Mat4::Identity();
    T.block<3,3>(0,0) = R;
    T.block<3,1>(0,3) = t;
    return T;
}

inline std::vector<Vec3> applyT(const std::vector<Vec3>& pts, const Mat4& T) {
    const Mat3 R = T.block<3,3>(0,0);
    const Vec3 t = T.block<3,1>(0,3);
    std::vector<Vec3> out; out.reserve(pts.size());
    for (const auto& p : pts) out.push_back(R * p + t);
    return out;
}

// Naive point-to-point ICP. Brute-force nearest neighbor is O(N*M):
// swap in a kd-tree (nanoflann / PCL) for production.
inline Mat4 icp(const std::vector<Vec3>& source,
                const std::vector<Vec3>& target,
                int max_iters = 50, double tol = 1e-9) {
    Mat4 T = Mat4::Identity();
    std::vector<Vec3> cur = source;            // moving copy of the source
    double prev = std::numeric_limits<double>::infinity();

    for (int it = 0; it < max_iters; ++it) {
        // 1. CORRESPOND: nearest target point for each current source point.
        std::vector<Vec3> ps, qs; ps.reserve(cur.size()); qs.reserve(cur.size());
        double err = 0.0;
        for (const auto& p : cur) {
            double best = std::numeric_limits<double>::infinity();
            Vec3 bq = target.front();
            for (const auto& q : target) {
                double d2 = (p - q).squaredNorm();
                if (d2 < best) { best = d2; bq = q; }
            }
            ps.push_back(p); qs.push_back(bq); err += best;
        }
        // 2. SOLVE incremental transform on the matches.
        Mat4 dT = bestRigidTransform(ps, qs);
        // 3. UPDATE: compose and move the source.
        T = dT * T;
        cur = applyT(cur, dT);
        // 4. STOP when the mean residual stops improving (relative tol).
        err /= double(cur.size());
        if (std::abs(prev - err) < tol) break;
        prev = err;
    }
    return T;
}

} // namespace reg
```

Library shortcut: Eigen ships the Procrustes solve as `Eigen::umeyama(src, dst, with_scaling)`, where `src` and `dst` are $d \times n$ matrices with points as **columns**. It returns the homogeneous $(d{+}1)\times(d{+}1)$ transform $T = \begin{bmatrix} cR & t \\ 0 & 1\end{bmatrix}$. Pass `with_scaling = false` for a pure rigid fit ($c = 1$). For full ICP on real clouds, use PCL:

```cpp
#include <pcl/registration/icp.h>
pcl::IterativeClosestPoint<pcl::PointXYZ, pcl::PointXYZ> icp;
icp.setInputSource(cloud_src);          // p_i
icp.setInputTarget(cloud_tgt);          // q_i
icp.setMaxCorrespondenceDistance(0.05); // reject matches farther than 5 cm (gating!)
icp.setMaximumIterations(50);
icp.setTransformationEpsilon(1e-9);     // relative stop on the transform
pcl::PointCloud<pcl::PointXYZ> aligned;
icp.align(aligned);                     // give it a good initial guess for a 2nd arg
if (icp.hasConverged())
    Eigen::Matrix4f T = icp.getFinalTransformation();
```

`pcl::GeneralizedIterativeClosestPoint` (GICP) is the production default; it is a plane-to-plane variant that is more robust than basic point-to-point on noisy stone scans.

## Connections

- [[03_eig-svd-pca]]: the rigid solve *is* an SVD of the $3\times3$ cross-covariance; PCA on each cloud also supplies the coarse axis alignment that initializes ICP.
- [[linear-least-squares]]: the point-to-plane step is a plain $6\times6$ linear least squares after small-angle linearization; the point-to-point step is its constrained (manifold) cousin solved in closed form.
- [[05_transforms-se3]]: the output is an element of $SE(3)$; composing ICP increments ($T \leftarrow \Delta T\, T$) and re-orthonormalizing $R$ are $SE(3)$ operations.
- [[08_nonlinear-least-squares]]: ICP is block coordinate descent; point-to-plane ICP is one Gauss-Newton step per iteration on the surface-normal residual.
- [[09_robust-ransac]]: distance gating and trimmed ICP are the robust-estimation defense against wrong correspondences and partial overlap.
- [[sdf-signed-heat]] (unlocks): you register and merge partial scans into one cloud *before* building a signed-distance field; misregistration shows up as a blurred or doubled surface in the SDF.

## Pitfalls and failure modes

- **Frame and direction confusion.** "Align A to B" and "align B to A" give inverse transforms. Pick one direction (source onto target) and stamp it everywhere. A swapped $H = \sum \tilde q_i \tilde p_i^T$ returns $R^T$, a silent inverse rotation.
- **Reflection without the $d$ term.** Dropping $d = \det(VU^T)$ lets noisy or planar data return a mirror image with $\det R = -1$. Always include the diagonal $\mathrm{diag}(1,1,d)$. Assert $\det R > 0$ after solving.
- **Bad initialization = wrong local minimum.** ICP only finds the nearest basin. Two flat faces can snap together flipped or rotated $180^\circ$ and ICP will happily report low residual. Coarse-align first (PCA, markers, FPFH+RANSAC) or use Go-ICP.
- **No correspondence gating.** Without `setMaxCorrespondenceDistance`, points in a non-overlapping region get matched to garbage and drag the fit. For partial scans this is the number-one bug. Gate, trim, or reweight.
- **Degenerate geometry.** A sphere, a flat plane, a straight edge, or a surface of revolution leaves a rotation unconstrained ($\sigma_2 \approx \sigma_3$, or $\sigma_3 \approx 0$). The SVD returns *an* answer but it is arbitrary along the free axis. Check singular-value gaps relative to $\sigma_1$.
- **Non-commutativity of updates.** $\Delta T \cdot T \ne T \cdot \Delta T$. The increment is in the *current* frame, so it must left-multiply. Getting the side wrong makes ICP wander.
- **Absolute tolerances across scales.** A $10^{-6}$ m stop is meaningless on a 2 m block and impossibly tight on a 100 mm cobble. Make the convergence test relative to object size; this matches the scale-invariance rule for the whole stone pipeline.
- **Outliers and partial overlap.** Plain least squares has a 0% breakdown point: one bad match shifts the fit. Use trimmed ICP (keep the best fraction of residuals) or robust weights.
- **Re-orthonormalization drift.** After many point-to-plane iterations, the rebuilt $R$ accumulates non-orthogonality. Re-project to $SO(3)$ (its own SVD, $R \leftarrow U V^T$) periodically.

## Practice

1. **By hand.** Redo the worked example but make the target a $90^\circ$ **clockwise** rotation of the source ($R = \begin{bmatrix}0&1\\-1&0\end{bmatrix}$) plus translation $(2, 0)$. Recompute $H$, verify $RH$ is symmetric PSD, and report $t$.
2. **Reflection trap.** Build two 3D point sets where the naive $R = VU^T$ has $\det = -1$ (hint: nearly coplanar points with noise). Solve once without the $d$ term and once with it. Print $\det R$ and the RMSD for both. Confirm the $d$-corrected fit is the proper rotation.
3. **ICP basin.** Take a synthetic stone cloud, apply a known transform with rotation angle $\theta$, and run `icp` from identity. Sweep $\theta$ from $0^\circ$ to $180^\circ$ and plot final RMSD vs $\theta$. Find the angle where ICP stops converging to the truth.
4. **Point-to-plane.** Estimate target normals (smallest-eigenvector PCA per neighborhood), implement the $6\times6$ point-to-plane solve, and compare its iteration count to point-to-point on a flat-faced block. Report iterations-to-converge for both.
5. **Gating.** Take two partially overlapping scans (delete 40% of the target). Run ICP with and without a max-correspondence-distance gate. Show the un-gated fit drifts and the gated one holds.
6. **Eigen vs PCL.** Solve the same paired-point problem with your `bestRigidTransform` and with `Eigen::umeyama(..., false)`. Confirm the transforms match to $10^{-12}$. Then register a real cloud pair with `pcl::IterativeClosestPoint` and compare `getFinalTransformation` to your naive ICP.

## References

- Umeyama, S. "Least-Squares Estimation of Transformation Parameters Between Two Point Patterns." *IEEE TPAMI*, 13(4), 1991, 376-380. https://doi.org/10.1109/34.88573 (the determinant correction; canonical rigid/similarity solve)
- Kabsch, W. "A solution for the best rotation to relate two sets of vectors." *Acta Cryst.* A32, 1976, 922-923. Summary: https://en.wikipedia.org/wiki/Kabsch_algorithm
- Arun, K. S., Huang, T. S., Blostein, S. D. "Least-Squares Fitting of Two 3-D Point Sets." *IEEE TPAMI*, 9(5), 1987, 698-700. (the SVD solve Umeyama later corrected)
- Besl, P. J. and McKay, N. D. "A Method for Registration of 3-D Shapes." *IEEE TPAMI*, 14(2), 1992, 239-256. https://doi.org/10.1109/34.121791 (the ICP algorithm, point-to-point)
- Chen, Y. and Medioni, G. "Object modelling by registration of multiple range images." *Image and Vision Computing*, 10(3), 1992. (point-to-plane ICP)
- Low, K.-L. "Linear Least-Squares Optimization for Point-to-Plane ICP Surface Registration." Technical Report TR04-004, UNC Chapel Hill, 2004. https://www.cs.unc.edu/techreports/04-004.pdf (the $6\times6$ small-angle linearization)
- Rusinkiewicz, S. and Levoy, M. "Efficient Variants of the ICP Algorithm." *3DIM*, 2001. https://www.cs.princeton.edu/~smr/papers/fasticp/ (correspondence, weighting, point-to-plane survey)
- Yang, J., Li, H., Campbell, D., Jia, Y. "Go-ICP: A Globally Optimal Solution to 3D ICP Point-Set Registration." *IEEE TPAMI*, 38(11), 2016. arXiv:1605.03344. https://arxiv.org/abs/1605.03344 ; code https://github.com/yangjiaolong/Go-ICP
- Eigen `umeyama` reference (Geometry module). https://eigen.tuxfamily.org/dox/group__Geometry__Module.html
- PCL `IterativeClosestPoint` class reference and tutorial. https://pointclouds.org/documentation/classpcl_1_1_iterative_closest_point.html ; https://pointclouds.org/documentation/tutorials/iterative_closest_point.html
- Lynch, K. and Park, F. *Modern Robotics*, Ch. 3 ($SE(3)$, rotations). https://hades.mech.northwestern.edu/index.php/Modern_Robotics
- CADFit: CAD program synthesis. arXiv:2605.01171 ; code https://github.com/ghadinehme/CADFit (canonical-frame primitive fitting the pipeline feeds)

## Flashcard seeds

- Q: Derive the optimal translation $t$ for fixed $R$ (step 1 of 2). <br> <b>Formula:</b> objective $\sum_i\|Rp_i+t-q_i\|^2$; centroids $\bar p=\tfrac1N\sum_i p_i,\ \bar q=\tfrac1N\sum_i q_i$; target result $t=\bar q-R\bar p$. <br> <i>Hint:</i> take $\partial/\partial t$, set it to $0$, and pull the constant $Nt$ out of the sum. :: A: $\dfrac{\partial}{\partial t}\sum_i\|Rp_i+t-q_i\|^2=2\sum_i(Rp_i+t-q_i)=0$ <br> $\sum_i(Rp_i)+Nt-\sum_i q_i=0$ <br> divide by $N$: $R\bar p+t-\bar q=0$ <br> $\boxed{t=\bar q-R\bar p}$.
- Q: Substitute $t=\bar q-R\bar p$ back in and derive the centered, pure-rotation objective (step 2 of 2). <br> <b>Formula:</b> centered points $\tilde p_i=p_i-\bar p,\ \tilde q_i=q_i-\bar q$; target result $\min_{R\in SO(3)}\sum_i\|R\tilde p_i-\tilde q_i\|^2$. <br> <i>Hint:</i> plug $t$ into $Rp_i+t-q_i$ and group the $\bar p,\bar q$ terms with their points. :: A: $Rp_i+t-q_i=Rp_i+\bar q-R\bar p-q_i=R(p_i-\bar p)-(q_i-\bar q)$ <br> rename $\tilde p_i=p_i-\bar p,\ \tilde q_i=q_i-\bar q$ <br> $\boxed{\min_{R\in SO(3)}\sum_i\|R\tilde p_i-\tilde q_i\|^2}$ — translation has dropped out.
- Q: Show the centered rotation objective reduces to a trace maximization (step 1 of 2). <br> <b>Formula:</b> $\|R\tilde p_i-\tilde q_i\|^2=\|R\tilde p_i\|^2+\|\tilde q_i\|^2-2\,\tilde q_i^TR\tilde p_i$; target result $\max_R\sum_i\tilde q_i^TR\tilde p_i$. <br> <i>Hint:</i> use $\|R\tilde p_i\|=\|\tilde p_i\|$ (orthogonal $R$) to drop the $R$-independent terms. :: A: expand one term: $\|R\tilde p_i-\tilde q_i\|^2=\|R\tilde p_i\|^2+\|\tilde q_i\|^2-2\,\tilde q_i^TR\tilde p_i$ <br> $\|R\tilde p_i\|=\|\tilde p_i\|$, so the first two terms are constant in $R$ <br> minimizing the sum $\equiv$ $\boxed{\max_R\sum_i\tilde q_i^TR\tilde p_i}$.
- Q: Turn $\max_R\sum_i\tilde q_i^TR\tilde p_i$ into the trace form and identify $H$ (step 2 of 2). <br> <b>Formula:</b> identity $a^Tb=\operatorname{tr}(ba^T)$; target result $\max_R\operatorname{tr}(RH),\ H=\sum_i\tilde p_i\tilde q_i^T$. <br> <i>Hint:</i> rewrite each scalar $\tilde q_i^TR\tilde p_i$ as a trace, then pull $R$ and the sum out. :: A: $\tilde q_i^TR\tilde p_i=\operatorname{tr}(R\tilde p_i\tilde q_i^T)$ <br> $\sum_i\operatorname{tr}(R\tilde p_i\tilde q_i^T)=\operatorname{tr}\!\big(R\sum_i\tilde p_i\tilde q_i^T\big)$ <br> $\boxed{\max_R\operatorname{tr}(RH),\quad H=\sum_i\tilde p_i\tilde q_i^T}$ (source-left, target-right convention).
- Q: From $H=U\Sigma V^T$, derive the orthogonal maximizer of $\operatorname{tr}(RH)$ before any reflection fix (step 1 of 2). <br> <b>Formula:</b> $\operatorname{tr}(RH)=\operatorname{tr}(\Sigma\,V^TRU)$, with $\Sigma=\mathrm{diag}(\sigma_1,\sigma_2,\sigma_3)$; target result $R=VU^T$. <br> <i>Hint:</i> set $M=V^TRU$ (orthogonal, so $|M_{jj}|\le1$) and bound $\operatorname{tr}(\Sigma M)$. :: A: $\operatorname{tr}(RH)=\operatorname{tr}(RU\Sigma V^T)=\operatorname{tr}(\Sigma\,V^TRU)$ <br> with $M=V^TRU$: $\operatorname{tr}(\Sigma M)=\sum_j\sigma_jM_{jj}\le\sum_j\sigma_j$ <br> bound hit when $M=I$, i.e. $V^TRU=I$ <br> $\boxed{R=VU^T}$.
- Q: From $R=VU^T$, derive the Umeyama-corrected rotation that stays in $SO(3)$ (step 2 of 2). <br> <b>Formula:</b> $d=\det(VU^T)\in\{+1,-1\}$; target result $R=V\,\mathrm{diag}(1,1,d)\,U^T$. <br> <i>Hint:</i> if $\det(VU^T)=-1$ flip the axis tied to the smallest $\sigma_3$ by inserting $\mathrm{diag}(1,1,d)$. :: A: $VU^T$ may have $\det=-1$ (a reflection); to enforce $\det R=+1$ flip the $\sigma_3$ axis, costing only $2\sigma_3$ <br> set $d=\det(VU^T)$ and insert $\mathrm{diag}(1,1,d)$ <br> $\boxed{R=V\,\mathrm{diag}(1,1,d)\,U^T,\ \ d=\det(VU^T)}$.
- Q: Once $R$ is known, derive the translation that re-attaches the clouds. <br> <b>Formula:</b> optimal $t$ sends $\bar p\mapsto\bar q$, i.e. $R\bar p+t=\bar q$; target result $t=\bar q-R\bar p$. <br> <i>Hint:</i> require the source centroid to land on the target centroid and solve for $t$. :: A: optimal $t$ maps source centroid onto target centroid <br> $R\bar p+t=\bar q$ <br> $\boxed{t=\bar q-R\bar p}$.
- Q: Derive the Umeyama optimal scale $c$ for a similarity fit. <br> <b>Formula:</b> $\min_c\sum_i\|cR\tilde p_i-\tilde q_i\|^2$ gives $c=\dfrac{\sum_i\tilde q_i^TR\tilde p_i}{\sum_i\|\tilde p_i\|^2}$, with $\sigma_p^2=\tfrac1N\sum_i\|\tilde p_i\|^2$; target result $c=\tfrac{1}{\sigma_p^2}\operatorname{tr}(\Sigma\,\mathrm{diag}(1,1,d))$. <br> <i>Hint:</i> differentiate in $c$, set to $0$, then replace the numerator by $\operatorname{tr}(\Sigma\,\mathrm{diag}(1,1,d))$ and the denominator by $N\sigma_p^2$. :: A: $\dfrac{d}{dc}\sum_i\|cR\tilde p_i-\tilde q_i\|^2=0\Rightarrow c=\dfrac{\sum_i\tilde q_i^TR\tilde p_i}{\sum_i\|\tilde p_i\|^2}$ <br> numerator $=\operatorname{tr}(\Sigma\,\mathrm{diag}(1,1,d))$, denominator $=N\sigma_p^2$ <br> $\boxed{c=\dfrac{1}{\sigma_p^2}\operatorname{tr}(\Sigma\,\mathrm{diag}(1,1,d))}$; for rigid registration set $c=1$.
- Q: Derive the per-correspondence row $a_i$ and value $b_i$ for the linearized point-to-plane step. <br> <b>Formula:</b> residual $\approx n_i^T(p_i-q_i)+n_i^T([r]_\times p_i)+n_i^Tt$, with $R\approx I+[r]_\times$ and identity $n^T([r]_\times p)=(p\times n)^Tr$; target row $a_i^Tx=b_i$. <br> <i>Hint:</i> swap $n_i^T([r]_\times p_i)$ for $(p_i\times n_i)^Tr$, then read off $a_i$ as the coefficient of $x=(r;t)$. :: A: residual $\approx n_i^T(p_i-q_i)+n_i^T([r]_\times p_i)+n_i^Tt$ <br> $n_i^T([r]_\times p_i)=(p_i\times n_i)^Tr$; stack $x=(r;t)$ <br> $\boxed{a_i=\begin{bmatrix}p_i\times n_i\\ n_i\end{bmatrix},\ b_i=-\,n_i^T(p_i-q_i)}$.
- Q: Worked example (2D, clean rotation). Source $p$: $(0,0),(2,0),(2,1),(0,1)$; target $q$: $(1,1),(1,3),(0,3),(0,1)$. Compute $H$. <br> <b>Formula:</b> $\bar p=\tfrac1N\sum p_i,\ \bar q=\tfrac1N\sum q_i$, $\tilde p_i=p_i-\bar p,\ \tilde q_i=q_i-\bar q$, $H=\sum_i\tilde p_i\tilde q_i^T$. <br> <i>Hint:</i> get the two centroids first, then center every point before summing the outer products. :: A: $\bar p=(1,0.5),\ \bar q=(0.5,2)$ <br> $\tilde p$: $(-1,-.5),(1,-.5),(1,.5),(-1,.5)$; $\tilde q$: $(.5,-1),(.5,1),(-.5,1),(-.5,-1)$ <br> $(1,1)$ entry $\sum\tilde p_{ix}\tilde q_{ix}=-.5{\times}4=-2$; $(1,2)=4$; $(2,1)=-1$; $(2,2)=2$ <br> $\boxed{H=\begin{bmatrix}-2&4\\-1&2\end{bmatrix}}$.
- Q: Worked example (continue, 2D). With $H=\begin{bmatrix}-2&4\\-1&2\end{bmatrix}$ and candidate $R=\begin{bmatrix}0&-1\\1&0\end{bmatrix}$ ($90^\circ$ CCW), verify optimality. <br> <b>Formula:</b> $R$ is the Procrustes optimizer iff $RH$ is symmetric PSD, and then $\operatorname{tr}(RH)=\sigma_1+\sigma_2$. <br> <i>Hint:</i> multiply $RH$, check symmetry, then get its eigenvalues from $\operatorname{tr}$ and $\det$. :: A: $RH=\begin{bmatrix}0&-1\\1&0\end{bmatrix}\begin{bmatrix}-2&4\\-1&2\end{bmatrix}=\begin{bmatrix}1&-2\\-2&4\end{bmatrix}$ <br> symmetric; trace $5$, det $0$ give eigenvalues $0,5\ge0$, so $RH\succeq0$ <br> $\boxed{\operatorname{tr}(RH)=5=\sigma_1+\sigma_2\Rightarrow R\text{ is the maximizer}}$.
- Q: Worked example (3D, mm shop scale, identity case). Two scans of a block are already aligned: $p_i=q_i$ for all $i$. Solve $(R,t)$. <br> <b>Formula:</b> $H=\sum_i\tilde p_i\tilde q_i^T$, $R=V\,\mathrm{diag}(1,1,d)\,U^T$, $d=\det(VU^T)$, $t=\bar q-R\bar p$. <br> <i>Hint:</i> with $p_i=q_i$ note $\tilde p_i=\tilde q_i$, so $H$ is symmetric PSD and $U=V$. :: A: $\bar p=\bar q$, so $\tilde p_i=\tilde q_i$ <br> $H=\sum\tilde p_i\tilde p_i^T\succeq0$ symmetric $\Rightarrow U=V$, $d=\det(VU^T)=+1$ <br> $R=VU^T=I$, $t=\bar q-I\bar p=0$ <br> $\boxed{(R,t)=(I,0)}$ — no motion, residual $0\,$mm.
- Q: Worked example (similarity fit, scale at site scale). Centered clouds match in shape but the target is $\times4$ bigger: $\tilde q_i=4R\tilde p_i$ with $R$ a pure rotation. Find $c$. <br> <b>Formula:</b> $c=\dfrac{\sum_i\tilde q_i^TR\tilde p_i}{\sum_i\|\tilde p_i\|^2}$. <br> <i>Hint:</i> substitute $\tilde q_i=4R\tilde p_i$ and use $\tilde q_i^TR\tilde p_i=4\|R\tilde p_i\|^2=4\|\tilde p_i\|^2$. :: A: $R$ from SVD aligns directions, $d=+1$ <br> $c=\dfrac{\sum 4\,(R\tilde p_i)^TR\tilde p_i}{\sum\|\tilde p_i\|^2}=\dfrac{4\sum\|\tilde p_i\|^2}{\sum\|\tilde p_i\|^2}$ <br> $\boxed{c=4}$ — the scale ratio is recovered exactly.
- Q: Worked example (reflection regime). The naive $R_0=VU^T$ returns $\det R_0=-1$ on a noisy near-planar cloud. Give the corrected $R$ and the cost of the fix. <br> <b>Formula:</b> $R=V\,\mathrm{diag}(1,1,d)\,U^T$, $d=\det(VU^T)$; flipping the smallest-$\sigma$ axis adds $2\sigma_3$ to the squared error. <br> <i>Hint:</i> here $d=-1$, so insert $\mathrm{diag}(1,1,-1)$ and recompute $\det R$. :: A: $d=\det(VU^T)=-1$, insert $\mathrm{diag}(1,1,-1)$ <br> $R=V\,\mathrm{diag}(1,1,-1)\,U^T$, now $\det R=(-1)\det(VU^T)=+1$ <br> flipping the $\sigma_3$ axis raises the objective by $2\sigma_3$ <br> $\boxed{\text{extra squared error}=2\sigma_3,\ \det R=+1}$.
- Q: Closing test (close the loop, 2D). Using $R=\begin{bmatrix}0&-1\\1&0\end{bmatrix}$, $\bar p=(1,0.5)$, $\bar q=(0.5,2)$, find $t$ then verify the fit on $p_1=(0,0)$ (target $q_1=(1,1)$). <br> <b>Formula:</b> $t=\bar q-R\bar p$, then check $Rp_1+t=q_1$. <br> <i>Hint:</i> a pass means the residual is $0$ and $\det R=+1$. :: A: $t=\bar q-R\bar p=(0.5,2)-\begin{bmatrix}0&-1\\1&0\end{bmatrix}(1,0.5)=(0.5,2)-(-0.5,1)=(1,1)$ <br> check: $Rp_1+t=R(0,0)+(1,1)=(1,1)=q_1$ <br> $\boxed{\text{residual}=0,\ \det R=+1}$ — exact fit on clean data.
- Q: Closing test (limiting case). Verify that when $H$ is symmetric positive-definite the solved rotation must be $R=I$. <br> <b>Formula:</b> SPD $H=U\Sigma U^T$, so $V=U$; $R=V\,\mathrm{diag}(1,1,d)\,U^T$ with $d=\det(VU^T)$. <br> <i>Hint:</i> a pass is $R=UU^T=I$; substitute $V=U$ and confirm $d=+1$. :: A: SPD $H$ has eigendecomposition $H=U\Sigma U^T$, so $V=U$ <br> $d=\det(UU^T)=+1$, so $R=U\,I\,U^T=UU^T=I$ <br> $\boxed{R=I}$ — already-aligned centered clouds need no rotation.
- Q: Closing test (convention check). A colleague computes the swapped $H'=\sum_i\tilde q_i\tilde p_i^T$. What rotation does the solve return, and how do you fix it? <br> <b>Formula:</b> swapping gives $H'=H^T$, which swaps $U\leftrightarrow V$ in the SVD, so $R'=R^T$. <br> <i>Hint:</i> a pass means recognizing $R^T$ is the inverse (target-onto-source) rotation. :: A: $H'=H^T\Rightarrow$ SVD swaps $U\leftrightarrow V$, so $R'=U\,\mathrm{diag}(1,1,d)\,V^T=R^T$ <br> $R^T$ is the inverse rotation <br> $\boxed{R'=R^T;\ \text{transpose it, or fix }H=\sum\tilde p_i\tilde q_i^T}$.
- Q: What problem does Kabsch-Umeyama solve in closed form, and what is its one expensive step? <br> <b>Formula:</b> $\min_{R\in SO(3),t}\sum_i\|Rp_i+t-q_i\|^2$ via $H=\sum_i\tilde p_i\tilde q_i^T$. <br> <i>Hint:</i> name the matrix decomposition done on $H$. :: A: optimal rigid $(R,t)$ minimizing $\sum_i\|Rp_i+t-q_i\|^2$ for known correspondences <br> $\boxed{\text{one }3\times3\text{ SVD of }H}$.
- Q: Why does ICP only reach a local minimum? <br> <i>Hint:</i> think about what it alternates between and whether either step can increase the error. :: A: it is coordinate descent on (correspondences, transform); each step is non-increasing but the fixed point depends entirely on the start pose.
- Q: When do you pick point-to-plane over point-to-point ICP? <br> <b>Formula:</b> point-to-plane cost $\sum_i\big(n_i^T(Rp_i+t-q_i)\big)^2$. <br> <i>Hint:</i> note what the normal $n_i$ lets the surfaces do tangentially. :: A: when reliable target normals exist; penalizing only the along-normal distance lets surfaces slide tangentially and converges in far fewer iterations.
- Q: What does Go-ICP guarantee that plain ICP cannot? <br> <i>Hint:</i> think about searching all of $SE(3)$ with bounds rather than from one start pose. :: A: globally optimal rigid registration via branch-and-bound over $SE(3)$ with L2 bounds, needing no initial guess.
- Q: Degeneracy guard: which geometries leave the rotation under-constrained, and what do you check? <br> <b>Formula:</b> inspect the singular values of $H$: $\sigma_2-\sigma_3$ and $\sigma_3$, both relative to $\sigma_1$. <br> <i>Hint:</i> a near-zero gap means the rotation about that axis is free. :: A: spheres, planes, straight lines, surfaces of revolution; check the gaps $\sigma_2-\sigma_3$ and $\sigma_3$ relative to $\sigma_1$.
- CLOZE: The optimal translation always maps the source centroid onto the {{c1::target centroid}}, giving $t=\bar q-R\bar p$.
- CLOZE: Without the $d=\det(VU^T)$ term the solve can return a {{c1::reflection}} with $\det R=-1$ on noisy or planar data.
- CLOZE: Point-to-plane ICP linearizes the rotation as $R\approx I+{{c1::[r]_\times}}$ with $r$ a small-angle vector. <br> <i>hint:</i> the skew-symmetric cross-product matrix of $r$.
- CLOZE: The ICP increment must {{c1::left}}-multiply ($T\leftarrow\Delta T\,T$) because the correction lives in the current frame and $SE(3)$ does not commute.
- CLOZE: ICP convergence is tested relative to {{c1::object size}} ($\varepsilon\cdot L$), never a fixed absolute tolerance across scales.
- Q: Eigen and PCL entry points for the rigid solve and full ICP. <br> <i>Hint:</i> one closed-form call for matched points, one iterative class for raw clouds. :: A: `Eigen::umeyama(src, dst, false)` for the closed-form rigid fit (points as columns); `pcl::IterativeClosestPoint` then `getFinalTransformation()` for iterative ICP.
