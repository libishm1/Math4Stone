# Eigenvalues, eigenvectors, SVD, and PCA

Prereqs: [[02_matrices]] => this => [[04_registration-icp]], [[05_remeshing-segmentation]], [[06_robust-ransac]]

## Why it matters for the stone scan->CAD->IFC / CADFit workflow

Every stage of the stone pipeline runs on this one topic. When you scan an irregular block, the raw point cloud has no orientation; you recover a local frame by taking the eigenvectors of the local covariance matrix (PCA), and the eigenvector for the smallest eigenvalue is the surface normal you need before any meshing or B-spline fitting. When you align two scans, or align a scan to a CAD template, the rigid-body solve is a single SVD of a $3\times3$ cross-covariance matrix (the Kabsch/Umeyama core inside ICP). When you fit a plane, a cylinder axis, or a least-squares solution to a noisy linear system, the condition number from the singular values tells you whether the answer is trustworthy or numerical garbage. CADFit-style program synthesis also leans on stable local frames and canonical axes, because a primitive fitted in a badly chosen frame produces a wrong CAD parameter. Master this and you can read, debug, and extend the geometric heart of the whole pipeline.

## Intuition

A square matrix is a transformation: it takes vectors and moves them. Most vectors get rotated and stretched at the same time. An **eigenvector** is a special direction that the matrix does not rotate; it only stretches it. The stretch factor is the **eigenvalue**. If you imagine the transformation acting on a sphere of vectors, the eigenvectors of a symmetric matrix are the axes of the resulting ellipsoid, and the eigenvalues are how far each axis is pushed out.

The **SVD** says any matrix, even a non-square one, is really just three simple moves in sequence: rotate, then scale along axes, then rotate again. $A = U\Sigma V^T$. $V^T$ rotates the input, $\Sigma$ stretches each axis by a singular value, $U$ rotates the result. No shearing, no skew once you pick the right axes. That is the entire geometric content.

**PCA** is the application that ties it to stone. Take a cloud of points near one spot on the surface. The covariance matrix describes the shape of that little blob. Its eigenvectors are the natural axes of the blob: two long axes lying flat in the surface, and one short axis poking out perpendicular to it. The short axis is the normal. Analogy: drop a handful of rice on a table and it spreads into a flat patch. The two directions it spreads are the tangent plane; the thin direction (table-to-ceiling) is the normal.

## The math (first principles)

### Conventions

Vectors are column vectors in $\mathbb{R}^3$ unless stated otherwise. Frames are right-handed. A rotation matrix $R$ satisfies $R^TR = I$ and $\det R = +1$; if $\det R = -1$ it is a reflection, which is a bug in a rigid-fit context. Points stack as **rows** in a data matrix $P \in \mathbb{R}^{N\times 3}$ (row $i$ is point $p_i^T$). Units are whatever the scan uses (meters at site scale, millimeters at shop scale; keep one unit per computation). Tolerances are relative: an eigenvalue is "zero" when it falls below $\varepsilon \cdot \lambda_{\max}$ with $\varepsilon$ near machine precision times problem size, not below a fixed absolute number.

### Eigenvalue problem

For a square matrix $A \in \mathbb{R}^{n\times n}$, a nonzero vector $v$ and scalar $\lambda$ form an eigenpair when

$$A v = \lambda v.$$

Geometrically $A$ does not change the direction of $v$, only its length. The eigenvalues are the roots of the characteristic polynomial $\det(A - \lambda I) = 0$. For each $\lambda$, the eigenvectors span the null space of $A - \lambda I$.

### Symmetric matrices have orthogonal eigenvectors (spectral theorem)

If $A = A^T$ (real symmetric), then all eigenvalues are real, and eigenvectors for distinct eigenvalues are orthogonal. You can choose an orthonormal eigenbasis, giving the **eigendecomposition**

$$A = Q \Lambda Q^T, \qquad Q^TQ = I, \quad \Lambda = \operatorname{diag}(\lambda_1,\dots,\lambda_n).$$

This is why covariance matrices (always symmetric, always positive semidefinite) give clean orthogonal frames. Quick proof of orthogonality: let $Av_1 = \lambda_1 v_1$ and $Av_2 = \lambda_2 v_2$ with $\lambda_1 \neq \lambda_2$. Then $\lambda_1 (v_1^Tv_2) = (Av_1)^Tv_2 = v_1^T A v_2 = \lambda_2 (v_1^Tv_2)$, using $A = A^T$. So $(\lambda_1 - \lambda_2)(v_1^Tv_2) = 0$, and since $\lambda_1 \neq \lambda_2$ we get $v_1^Tv_2 = 0$.

### Singular value decomposition

For **any** matrix $A \in \mathbb{R}^{m\times n}$ there exist orthogonal $U \in \mathbb{R}^{m\times m}$, orthogonal $V \in \mathbb{R}^{n\times n}$, and diagonal $\Sigma \in \mathbb{R}^{m\times n}$ with nonnegative entries $\sigma_1 \ge \sigma_2 \ge \dots \ge 0$ such that

$$A = U \Sigma V^T.$$

The $\sigma_i$ are the **singular values**, the columns of $U$ are **left singular vectors**, the columns of $V$ are **right singular vectors**. Geometry: $V^T$ rotates, $\Sigma$ scales along axes, $U$ rotates again. So a matrix is rotate-scale-rotate.

Link to eigendecomposition (verified): the right singular vectors are the eigenvectors of $A^TA$, and the singular values are the square roots of its eigenvalues:

$$A^TA = V \Sigma^T\Sigma V^T, \qquad \sigma_i = \sqrt{\lambda_i(A^TA)}.$$

For a symmetric positive semidefinite matrix, SVD and eigendecomposition coincide ($U = V = Q$, $\sigma_i = \lambda_i$).

### PCA = eigendecomposition of the covariance

Given points $p_1,\dots,p_N \in \mathbb{R}^3$, the centroid and covariance are

$$\bar p = \frac{1}{N}\sum_{i=1}^N p_i, \qquad C = \frac{1}{N-1}\sum_{i=1}^N (p_i - \bar p)(p_i - \bar p)^T.$$

$C$ is $3\times3$, symmetric, positive semidefinite. Its eigendecomposition $C = Q\Lambda Q^T$ orders directions by variance. With $\lambda_0 \le \lambda_1 \le \lambda_2$:

- eigenvector for $\lambda_2$ = direction of maximum spread (first principal axis),
- eigenvector for $\lambda_1$ = second tangent direction,
- eigenvector for $\lambda_0$ (smallest) = **surface normal** $n$, the least-spread direction.

The columns $[\,v_2 \;\; v_1 \;\; n\,]$ form a local frame. Force right-handedness by setting the third column to the cross product of the first two, $n = v_2 \times v_1$, then re-derive $v_1 = n \times v_2$ for orthonormality.

Surface curvature estimate (PCL convention, verified):

$$\sigma_{\text{curv}} = \frac{\lambda_0}{\lambda_0 + \lambda_1 + \lambda_2}.$$

The normal direction is only defined up to sign from PCA. Orient it consistently, for example toward the sensor viewpoint $v_p$: flip $n$ if $n^T(v_p - \bar p) < 0$.

### Rigid alignment via SVD (Kabsch / Umeyama)

Given correspondences $p_i \leftrightarrow q_i$, find rotation $R$ and translation $t$ minimizing

$$\min_{R \in SO(3),\, t} \sum_i \lVert R p_i + t - q_i \rVert^2.$$

Center both sets ($\tilde p_i = p_i - \bar p$, $\tilde q_i = q_i - \bar q$). Build the $3\times3$ cross-covariance

$$H = \sum_i \tilde p_i \tilde q_i^T \quad (=\; \tilde P^T \tilde Q \text{ in row-stacked form}).$$

Take $H = U\Sigma V^T$. The optimal rotation that maps source $p$ onto target $q$ is

$$d = \det(V U^T), \qquad R = V \begin{bmatrix} 1 & 0 & 0 \\ 0 & 1 & 0 \\ 0 & 0 & d \end{bmatrix} U^T, \qquad t = \bar q - R \bar p.$$

The $\det$ correction forces $\det R = +1$. Without it, severely corrupted or coplanar data can return a reflection ($\det R = -1$), which silently mirrors your geometry. (Convention note: this $H = \tilde P^T\tilde Q$ with $R = V\,\mathrm{diag}(1,1,d)\,U^T$ maps source onto target. Wikipedia uses $H = \tilde P^T\tilde Q$ but writes $R = U\,\mathrm{diag}(1,1,d)\,V^T$ to map target onto source; pick one and stay consistent.)

### Condition number and numerical stability

The 2-norm condition number of a matrix is the ratio of extreme singular values:

$$\kappa_2(A) = \frac{\sigma_{\max}}{\sigma_{\min}}.$$

$\kappa = 1$ is perfectly conditioned; $\kappa \to \infty$ means near-singular. A linear solve or least-squares fit can lose roughly $\log_{10}\kappa$ digits of accuracy. If $\sigma_{\min} = 0$ the matrix is rank-deficient and $\kappa = \infty$: the fit is geometrically degenerate (for plane fitting, the points are collinear; for ICP, the correspondences do not constrain all axes). Always check the singular values, not just whether the solver returned without error.

## Worked example

### PCA on five points (by hand)

Points (a flat patch in the $z=0$ plane with one bump):

$$p_1=(0,0,0),\; p_2=(2,0,0),\; p_3=(0,2,0),\; p_4=(2,2,0),\; p_5=(1,1,0).$$

Centroid: $\bar p = \frac{1}{5}(5,5,0) = (1,1,0)$.

Centered points: $(-1,-1,0),(1,-1,0),(-1,1,0),(1,1,0),(0,0,0)$.

Covariance (using $N-1 = 4$):

$$C = \frac{1}{4}\sum \tilde p_i \tilde p_i^T.$$

Sum of $\tilde x^2 = 1+1+1+1+0 = 4$. Sum of $\tilde y^2 = 4$. Sum of $\tilde x\tilde y = 1-1-1+1+0 = 0$. All $z$ terms are zero.

$$C = \frac{1}{4}\begin{bmatrix} 4 & 0 & 0 \\ 0 & 4 & 0 \\ 0 & 0 & 0 \end{bmatrix} = \begin{bmatrix} 1 & 0 & 0 \\ 0 & 1 & 0 \\ 0 & 0 & 0 \end{bmatrix}.$$

Eigenvalues: $\lambda = 1, 1, 0$. The smallest is $\lambda_0 = 0$ with eigenvector $(0,0,1)$. That is the normal: the patch lies flat, so its normal points straight up. Correct. Curvature $= 0/(0+1+1) = 0$, a perfectly flat patch.

### Rigid solve (by hand)

Source $p_1=(0,0,0), p_2=(1,0,0)$. Target is the source rotated $90°$ about $z$: $q_1=(0,0,0), q_2=(0,1,0)$.

Centroids: $\bar p = (0.5,0,0)$, $\bar q = (0,0.5,0)$. Centered source: $(-0.5,0,0),(0.5,0,0)$. Centered target: $(0,-0.5,0),(0,0.5,0)$.

$$H = \sum \tilde p_i \tilde q_i^T = (-0.5,0,0)^T(0,-0.5,0) + (0.5,0,0)^T(0,0.5,0) = \begin{bmatrix} 0 & 0.5 & 0 \\ 0 & 0 & 0 \\ 0 & 0 & 0 \end{bmatrix}.$$

By inspection $H$ maps $x$ to $y$. The rotation recovered by the SVD core is

$$R = \begin{bmatrix} 0 & -1 & 0 \\ 1 & 0 & 0 \\ 0 & 0 & 1 \end{bmatrix},$$

which is exactly the $+90°$ rotation about $z$, and $t = \bar q - R\bar p = (0,0.5,0) - (0,0.5,0) = 0$. The solve recovers the truth.

## Code

Idiomatic C++ with Eigen. PCA frame plus normal, then the Kabsch rigid solve. Both are small and compilable against Eigen 3.

```cpp
// eig_svd_pca.cpp
// build: g++ -O2 -I /path/to/eigen eig_svd_pca.cpp -o demo
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>   // SelfAdjointEigenSolver
#include <iostream>
#include <vector>

using Eigen::Vector3d;
using Eigen::Matrix3d;
using Eigen::Matrix4d;

// PCA local frame.
// Returns R whose columns are [tangent1, tangent2, normal].
// normal = eigenvector of smallest eigenvalue of the covariance C.
Matrix3d pcaFrame(const std::vector<Vector3d>& pts, Vector3d& centroid) {
    const int N = static_cast<int>(pts.size());
    centroid = Vector3d::Zero();
    for (const auto& p : pts) centroid += p;
    centroid /= static_cast<double>(N);          // p_bar

    Matrix3d C = Matrix3d::Zero();                // covariance C
    for (const auto& p : pts) {
        Vector3d d = p - centroid;                // p_i - p_bar
        C += d * d.transpose();                   // outer product
    }
    C /= static_cast<double>(N - 1);              // unbiased estimate

    // C is symmetric: use the self-adjoint solver. Eigenvalues come out
    // in ASCENDING order, so col(0) is the smallest = the normal.
    Eigen::SelfAdjointEigenSolver<Matrix3d> es(C);
    Vector3d normal = es.eigenvectors().col(0);   // lambda_0 -> n
    Vector3d t1     = es.eigenvectors().col(2);   // lambda_2 -> max spread
    Vector3d t2     = normal.cross(t1).normalized();
    t1              = t2.cross(normal).normalized(); // re-orthonormalize

    Matrix3d R;
    R.col(0) = t1; R.col(1) = t2; R.col(2) = normal;
    return R;
}

// Kabsch/Umeyama rigid solve. Returns 4x4 T mapping src onto dst: dst ~= R*src + t.
Matrix4d kabsch(const std::vector<Vector3d>& src, const std::vector<Vector3d>& dst) {
    const int N = static_cast<int>(src.size());
    Vector3d cs = Vector3d::Zero(), cd = Vector3d::Zero();
    for (int i = 0; i < N; ++i) { cs += src[i]; cd += dst[i]; }
    cs /= N; cd /= N;                              // p_bar, q_bar

    Matrix3d H = Matrix3d::Zero();                 // H = sum (p_i - p_bar)(q_i - q_bar)^T
    for (int i = 0; i < N; ++i)
        H += (src[i] - cs) * (dst[i] - cd).transpose();

    Eigen::JacobiSVD<Matrix3d> svd(H, Eigen::ComputeFullU | Eigen::ComputeFullV);
    Matrix3d U = svd.matrixU();
    Matrix3d V = svd.matrixV();

    // d = det(V U^T) forces a proper rotation (no reflection).
    double d = (V * U.transpose()).determinant();
    Matrix3d D = Matrix3d::Identity();
    D(2, 2) = (d < 0) ? -1.0 : 1.0;
    Matrix3d R = V * D * U.transpose();            // R = V diag(1,1,sign(d)) U^T
    Vector3d t = cd - R * cs;                       // t = q_bar - R p_bar

    Matrix4d T = Matrix4d::Identity();
    T.block<3,3>(0,0) = R;
    T.block<3,1>(0,3) = t;
    return T;
}

int main() {
    std::vector<Vector3d> pts = {
        {0,0,0},{2,0,0},{0,2,0},{2,2,0},{1,1,0}
    };
    Vector3d c;
    Matrix3d R = pcaFrame(pts, c);
    std::cout << "centroid: " << c.transpose() << "\n";
    std::cout << "normal (col 2):\n" << R.col(2).transpose() << "\n";

    std::vector<Vector3d> src = {{0,0,0},{1,0,0}};
    std::vector<Vector3d> dst = {{0,0,0},{0,1,0}};
    std::cout << "rigid T:\n" << kabsch(src, dst) << "\n";

    // Condition number of a matrix from its singular values.
    Matrix3d A; A << 1,0,0, 0,1e-6,0, 0,0,1;
    Eigen::JacobiSVD<Matrix3d> svd(A);
    auto sv = svd.singularValues();
    std::cout << "kappa = " << sv(0)/sv(sv.size()-1) << "\n"; // sigma_max/sigma_min
    return 0;
}
```

Key idiom: use `SelfAdjointEigenSolver` for symmetric matrices (covariances), not the general `EigenSolver`. It is faster, returns real eigenvalues in ascending order, and avoids complex arithmetic. Use `JacobiSVD` for small $3\times3$ SVD; it is accurate and the right tool for the Kabsch core.

Python sanity check (NumPy), for cross-validating the C++:

```python
import numpy as np
def kabsch(src, dst):
    cs, cd = src.mean(0), dst.mean(0)
    H = (src - cs).T @ (dst - cd)          # H = P'^T Q'
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T)) # det(V U^T)
    D = np.diag([1, 1, d])
    R = Vt.T @ D @ U.T                     # R = V diag(1,1,d) U^T
    t = cd - R @ cs
    return R, t
```

## Connections

- [[02_matrices]]: eigen/SVD operate on matrices; you need matrix multiplication, transpose, determinant, and orthogonality first.
- [[04_registration-icp]]: each ICP iteration calls the Kabsch SVD solve for the rigid update; this topic is its mathematical core.
- [[05_remeshing-segmentation]]: PCA normals and curvature from eigenvalue ratios drive feature detection, segmentation, and feature-aware remeshing.
- [[06_robust-ransac]]: RANSAC fits primitives (planes, cylinders) whose inlier sets are then refined by PCA/SVD; the condition number flags degenerate (collinear, coplanar) samples.
- [[07_least-squares]]: SVD gives the most stable least-squares solver and the condition number that predicts its accuracy.
- [[09_rigid-transforms-se3]]: the $R$ from Kabsch must satisfy $R \in SO(3)$; the $\det$ correction is what keeps it a rotation rather than a reflection.

## Pitfalls and failure modes

- **Reflection instead of rotation.** Skipping the $\det(VU^T)$ correction returns $\det R = -1$ on coplanar or noisy data. The geometry mirrors and downstream booleans/IFC fail silently. Always apply the diagonal correction.
- **Sign-flipped normals.** PCA gives the normal axis but not its orientation. Adjacent patches can flip, producing inside-out surfaces. Orient all normals consistently against a viewpoint or by propagation over a neighbor graph.
- **Eigenvalue order assumption.** `SelfAdjointEigenSolver` returns ascending eigenvalues (smallest first); many references and LAPACK wrappers differ. Pick the normal by smallest eigenvalue explicitly, never by a hardcoded column index copied from a different library.
- **Degenerate covariance.** A collinear neighborhood gives two near-zero eigenvalues and an undefined normal. Too few neighbors gives a rank-deficient $C$. Check that $\lambda_1 \gg \lambda_0$ before trusting the normal; reject or enlarge the neighborhood otherwise.
- **Ill-conditioning.** A large $\kappa_2 = \sigma_{\max}/\sigma_{\min}$ means the fit is unstable; small data noise becomes large parameter error. Do not invert near-singular matrices; use SVD with a singular-value floor (truncate $\sigma_i$ below a relative tolerance).
- **Scale and units.** Covariance entries scale with the square of your length unit; eigenvalue thresholds must be relative ($\varepsilon\,\lambda_{\max}$), never a fixed absolute number, or they break when you switch meters to millimeters.
- **Non-commutativity.** $R$ then $t$ is not $t$ then $R$. The Kabsch translation is $t = \bar q - R\bar p$, computed after $R$. Reversing the order misplaces the cloud.
- **Outliers.** A single bad correspondence corrupts $H$ and tilts the whole rigid solve, because least squares squares the error. Reject outliers (RANSAC, trimming) before the SVD solve, not after.

## Practice

1. **Hand PCA.** Take the four corners of a unit square in the $z=0$ plane plus one point lifted to $z=0.5$. Compute the covariance by hand, find the three eigenvalues, and confirm the smallest-eigenvalue eigenvector tilts away from $(0,0,1)$. Inspect: does the normal lean toward the lifted corner?
2. **Reflection bug.** Generate two paired point sets that are mirror images (negate one axis). Run the Kabsch solve with and without the $\det$ correction. Print $\det R$ in both cases and confirm one is $-1$.
3. **Normal field.** Load or synthesize a noisy hemisphere point cloud. For each point, PCA over its 20 nearest neighbors, output the normal. Orient all normals outward from the center and write them as short line segments to an OBJ for visual inspection.
4. **Condition number sweep.** Build $3\times3$ matrices with singular values $(1,1,s)$ for $s = 1, 10^{-2}, 10^{-6}, 0$. Solve $Ax=b$ for a fixed $b$, perturb $b$ slightly, and measure how much $x$ changes. Plot the change against $\kappa$ and confirm it tracks $\log_{10}\kappa$ lost digits.
5. **ICP step.** Implement one ICP iteration: nearest-neighbor correspondences, then Kabsch. Apply a known rotation+translation to a copy of a cloud, then recover it. Report the residual rotation error in degrees after one step and after ten.
6. **Curvature map.** Using the PCA eigenvalues, compute $\lambda_0/(\lambda_0+\lambda_1+\lambda_2)$ per point on a scan with both flat and sharp regions. Color the cloud by curvature and confirm edges light up.

## References

- Kabsch algorithm, Wikipedia. Cross-covariance, SVD, determinant correction, optimal rotation. https://en.wikipedia.org/wiki/Kabsch_algorithm
- Hunter Heidenreich, "Kabsch Algorithm" and "Umeyama's Method: Corrected SVD for Point Alignment." Centered covariance $H=P'^TQ'$, $d=\mathrm{sign}(\det(VU^T))$, $R=V\,\mathrm{diag}(1,1,d)\,U^T$. https://hunterheidenreich.com/posts/kabsch-algorithm/ and https://hunterheidenreich.com/notes/interdisciplinary/computational-biology/umeyama-similarity-transformation/
- S. Umeyama, "Least-Squares Estimation of Transformation Parameters Between Two Point Patterns," IEEE TPAMI, 1991. The proper-rotation correction. https://web.stanford.edu/class/cs273/refs/umeyama.pdf
- PointCloudLibrary, "Estimating Surface Normals in a PointCloud." PCA on the neighborhood covariance; normal = smallest-eigenvalue eigenvector; curvature from eigenvalue ratio. https://pointclouds.org/documentation/tutorials/normal_estimation.html
- Understanding Linear Algebra (Austin), "Singular Value Decompositions," and the spectral theorem chapter. SVD as generalized orthogonal diagonalization. https://understandinglinearalgebra.org/sec-svd-intro.html
- CS 357 (Illinois), "Singular Value Decompositions" and "Condition Numbers." $\kappa_2 = \sigma_{\max}/\sigma_{\min}$. https://cs357.cs.illinois.edu/textbook/notes/condition.html
- Eigen documentation, `SelfAdjointEigenSolver` and `JacobiSVD`. C++ idioms used above. https://eigen.tuxfamily.org/dox/
- K. Crane, "Discrete Differential Geometry: An Applied Introduction," CMU. Covariance, curvature, and local-frame thinking for meshes. https://www.cs.cmu.edu/~kmcrane/Projects/DDG/
- L. Piegl, "On NURBS: A Survey," IEEE CG&A, 1991. Affine invariance of fitted geometry (downstream of PCA frames). https://ieeexplore.ieee.org/document/75590
- OpenCASCADE reference manual, `GeomAPI_PointsToBSplineSurface`. Ordered-grid fitting fed by PCA-aligned resampling. https://dev.opencascade.org/doc/overview/html/
- CADFit: CAD program synthesis. arXiv:2605.01171; code https://github.com/ghadinehme/CADFit

## Flashcard seeds

- Q: Derive the equation whose roots are the eigenvalues of a square matrix $A$. <br> <b>Formula:</b> start from $Av=\lambda v$ ($v\neq 0$); target $\det(A-\lambda I)=0$ <br> <i>Hint:</i> move everything to one side and factor out $v$, then ask when a nonzero $v$ can exist. :: A: rewrite as $Av-\lambda v=0$ <br> factor: $(A-\lambda I)v=0$ <br> a nonzero $v$ exists only if $A-\lambda I$ is singular (not invertible) <br> $\boxed{\det(A-\lambda I)=0}$
- Q: Prove eigenvectors of distinct eigenvalues of a symmetric $A=A^T$ are orthogonal (step 1 of 2). <br> <b>Formula:</b> $Av_1=\lambda_1 v_1$, $Av_2=\lambda_2 v_2$, with $A=A^T$ <br> <i>Hint:</i> compute the scalar $v_2^TAv_1$ two different ways and set them equal. :: A: way 1: $v_2^T(Av_1)=\lambda_1(v_2^Tv_1)$ <br> way 2: $(Av_2)^Tv_1=\lambda_2(v_2^Tv_1)$, using $A=A^T$ so $v_2^TA=(Av_2)^T$ <br> the two expressions are the same scalar, so equate them <br> $\boxed{(\lambda_1-\lambda_2)(v_1^Tv_2)=0}$
- Q: Finish the orthogonality proof (step 2 of 2). <br> <b>Formula:</b> $(\lambda_1-\lambda_2)(v_1^Tv_2)=0$ with $\lambda_1\neq\lambda_2$ <br> <i>Hint:</i> divide by the factor you know is nonzero. :: A: distinct eigenvalues give $\lambda_1-\lambda_2\neq 0$ <br> divide both sides by that nonzero scalar <br> $\boxed{v_1^Tv_2=0}$ so the eigenvectors are orthogonal
- Q: From the SVD, derive the eigendecomposition of $A^TA$ (step 1 of 2). <br> <b>Formula:</b> $A=U\Sigma V^T$ with $U^TU=I$; target $A^TA=V(\Sigma^T\Sigma)V^T$ <br> <i>Hint:</i> write out $A^TA=(U\Sigma V^T)^T(U\Sigma V^T)$, then cancel $U^TU$. :: A: $A^TA=(U\Sigma V^T)^T(U\Sigma V^T)$ <br> $=V\Sigma^T U^TU\Sigma V^T$ <br> use $U^TU=I$ to drop the middle pair <br> $\boxed{A^TA=V(\Sigma^T\Sigma)V^T}$
- Q: Relate the singular values $\sigma_i$ to the eigenvalues of $A^TA$ (step 2 of 2). <br> <b>Formula:</b> $A^TA=V(\Sigma^T\Sigma)V^T$; target $\sigma_i=\sqrt{\lambda_i(A^TA)}$ <br> <i>Hint:</i> recognize this as an eigendecomposition and read off the $i$-th diagonal entry. :: A: $V$ holds eigenvectors and $\Sigma^T\Sigma$ is diagonal with the eigenvalues <br> the $i$-th diagonal entry is $\sigma_i^2=\lambda_i(A^TA)$ <br> take the nonnegative square root <br> $\boxed{\sigma_i=\sqrt{\lambda_i(A^TA)}}$
- Q: Derive the covariance matrix of a centered point neighborhood. <br> <b>Formula:</b> centroid $\bar p=\frac1N\sum_i p_i$; target $C=\frac{1}{N-1}\sum_i (p_i-\bar p)(p_i-\bar p)^T$ <br> <i>Hint:</i> first center each point, then average the outer products $\tilde p_i\tilde p_i^T$ using the unbiased divisor $N-1$. :: A: center each point: $\tilde p_i=p_i-\bar p$ <br> each contributes the outer product $\tilde p_i\tilde p_i^T$ <br> average with the unbiased divisor $N-1$ <br> $\boxed{C=\frac{1}{N-1}\sum_i (p_i-\bar p)(p_i-\bar p)^T}$
- Q: From the covariance eigendecomposition, derive which eigenvector is the surface normal. <br> <b>Formula:</b> $C=Q\Lambda Q^T$, $\lambda_0\le\lambda_1\le\lambda_2$; each $\lambda$ = variance along its eigenvector axis <br> <i>Hint:</i> the surface is locally flat, so ask along which axis the spread is smallest. :: A: each $\lambda$ is the variance (spread) along its eigenvector axis <br> a locally flat patch spreads least perpendicular to itself <br> least spread = smallest eigenvalue $\lambda_0$ <br> $\boxed{n=\text{eigenvector of }\lambda_0}$
- Q: Build a right-handed orthonormal local frame from the max-spread axis $v_2$ and the normal $n$. <br> <b>Formula:</b> $v_1=n\times v_2$; frame $[\,v_2\;\;v_1\;\;n\,]$ <br> <i>Hint:</i> the third tangent must be orthogonal to both $n$ and $v_2$ — take their cross product. :: A: the missing tangent must be orthogonal to both $n$ and $v_2$ <br> cross product gives it: $v_1=n\times v_2$ <br> $v_2$, $v_1=n\times v_2$, $n$ are mutually orthogonal and unit <br> $\boxed{[\,v_2\;\;n\times v_2\;\;n\,]\ \text{is right-handed orthonormal}}$
- Q: Derive the PCL surface-curvature estimate as a normalized least-spread ratio. <br> <b>Formula:</b> target $\sigma_{\text{curv}}=\dfrac{\lambda_0}{\lambda_0+\lambda_1+\lambda_2}$ <br> <i>Hint:</i> total spread is the sum of the three eigenvalues; put the smallest over the total so a flat patch gives $0$. :: A: total spread is the sum of variances $\lambda_0+\lambda_1+\lambda_2$ <br> flatness is captured by the smallest share $\lambda_0$ <br> normalize so a flat patch ($\lambda_0=0$) gives $0$ <br> $\boxed{\sigma_{\text{curv}}=\dfrac{\lambda_0}{\lambda_0+\lambda_1+\lambda_2}}$
- Q: Derive the rule to flip a PCA normal toward the sensor viewpoint $v_p$. <br> <b>Formula:</b> the normal should satisfy $n^T(v_p-\bar p)>0$ <br> <i>Hint:</i> the sign of the dot product tells you if $n$ points toward or away from the sensor; flip when it is negative. :: A: PCA fixes the normal axis but not its sign <br> the normal should point toward the sensor: $n^T(v_p-\bar p)>0$ <br> a negative dot product means $n$ points away <br> $\boxed{\text{flip }n\to -n\text{ when }n^T(v_p-\bar p)<0}$
- Q: From the rigid-fit objective, derive the cross-covariance $H$ used to solve for $R$ (step 1 of 2). <br> <b>Formula:</b> $\min_{R,t}\sum_i\lVert Rp_i+t-q_i\rVert^2$; target $H=\sum_i(p_i-\bar p)(q_i-\bar q)^T$ <br> <i>Hint:</i> the optimal $t$ centers both clouds; then expand the squared norm and collect the $R$-dependent term as $\mathrm{tr}(RH)$. :: A: optimal $t$ centers both clouds, leaving $\tilde p_i,\tilde q_i$ <br> the rotation part maximizes $\sum_i \tilde q_i^T R\tilde p_i=\mathrm{tr}(R\sum_i\tilde p_i\tilde q_i^T)$ <br> collect the $3\times3$ matrix inside the trace <br> $\boxed{H=\sum_i(p_i-\bar p)(q_i-\bar q)^T}$
- Q: From $H=U\Sigma V^T$, derive the optimal proper rotation (step 2 of 2). <br> <b>Formula:</b> $\mathrm{tr}(RH)$ is maximized by $R=VU^T$; correction $d=\det(VU^T)$ <br> <i>Hint:</i> the bare maximizer can be a reflection ($\det=-1$) — insert $\mathrm{diag}(1,1,d)$ to force $\det R=+1$. :: A: $\mathrm{tr}(RH)$ is maximized by $R=VU^T$ <br> but this can have $\det=-1$ (a reflection) on degenerate data <br> insert $\mathrm{diag}(1,1,d)$ with $d=\det(VU^T)$ to force $\det R=+1$ <br> $\boxed{R=V\,\mathrm{diag}(1,1,d)\,U^T}$
- Q: Derive the Kabsch translation $t$ once $R$ is known. <br> <b>Formula:</b> centroids must map exactly: $R\bar p+t=\bar q$ <br> <i>Hint:</i> isolate $t$ by subtracting $R\bar p$ from both sides. :: A: $R$ already maps the centered clouds, so the centroids must map onto each other <br> $R\bar p+t=\bar q$ <br> solve for $t$ <br> $\boxed{t=\bar q-R\bar p}$
- Q: Derive why the 2-norm condition number is the ratio of extreme singular values. <br> <b>Formula:</b> $A=U\Sigma V^T$ maps the unit sphere to an ellipsoid with axes $\sigma_i$; target $\kappa_2(A)=\sigma_{\max}/\sigma_{\min}$ <br> <i>Hint:</i> the largest stretch is $\sigma_{\max}$, the smallest is $\sigma_{\min}$; the ratio bounds relative error growth. :: A: $A=U\Sigma V^T$ stretches the unit sphere into an ellipsoid with axes $\sigma_i$ <br> worst-case amplification is $\sigma_{\max}$, worst-case shrink is $\sigma_{\min}$ <br> their ratio bounds relative error growth <br> $\boxed{\kappa_2(A)=\sigma_{\max}/\sigma_{\min}}$
- Q: Worked example (2D flat patch, m scale). Five points $(0,0),(2,0),(0,2),(2,2),(1,1)$ in the $z=0$ plane (meters). Find $C$, its eigenvalues, and the normal. <br> <b>Formula:</b> $\bar p=\frac1N\sum p_i$, $C=\frac{1}{N-1}\sum(p_i-\bar p)(p_i-\bar p)^T$, $N=5$ <br> <i>Hint:</i> compute the centroid first, then the sums $\sum\tilde x^2$, $\sum\tilde y^2$, $\sum\tilde x\tilde y$ to fill $C$. :: A: centroid $\bar p=(1,1,0)$ <br> centered: $(-1,-1,0),(1,-1,0),(-1,1,0),(1,1,0),(0,0,0)$ <br> $\sum\tilde x^2=\sum\tilde y^2=4$, $\sum\tilde x\tilde y=0$, all $z=0$ <br> $C=\frac14\,\mathrm{diag}(4,4,0)=\mathrm{diag}(1,1,0)$ <br> eigenvalues $1,1,0$ <br> $\boxed{n=(0,0,1),\ \lambda_0=0}$
- Q: Worked example (mm scale, scale check). The same five-point square is now in millimeters (side $2\text{ mm}$). What is $C$ and the curvature $\sigma_{\text{curv}}$? <br> <b>Formula:</b> $C=\frac{1}{N-1}\sum\tilde p_i\tilde p_i^T$; $\sigma_{\text{curv}}=\lambda_0/(\lambda_0+\lambda_1+\lambda_2)$ <br> <i>Hint:</i> the coordinate numbers are identical, only the unit changes; note $C$ carries units of length$^2$ but curvature is a ratio. :: A: same numbers, so $C=\mathrm{diag}(1,1,0)\ \text{mm}^2$ (entries unchanged, units now mm$^2$) <br> eigenvalues $1,1,0$ in mm$^2$ <br> curvature is a ratio, so the mm$^2$ units cancel <br> $\boxed{\sigma_{\text{curv}}=0/(0+1+1)=0\ \text{(unitless, scale-free)}}$
- Q: Worked example (tilted patch, 3D). Covariance eigenvalues are $\lambda_0=0.01,\ \lambda_1=0.9,\ \lambda_2=1.0$. Compute the curvature and judge flatness. <br> <b>Formula:</b> $\sigma_{\text{curv}}=\dfrac{\lambda_0}{\lambda_0+\lambda_1+\lambda_2}$ <br> <i>Hint:</i> sum the three eigenvalues for the denominator first, then divide. :: A: denominator $=0.01+0.9+1.0=1.91$ <br> $\sigma_{\text{curv}}=0.01/1.91$ <br> $\boxed{\sigma_{\text{curv}}\approx 0.0052}$ near-flat; normal trustworthy since $\lambda_1\gg\lambda_0$
- Q: Worked example (rigid solve, $90°$ about $z$). Source $(0,0,0),(1,0,0)$; target $(0,0,0),(0,1,0)$. Build $H$ and read off $R$. <br> <b>Formula:</b> $H=\sum_i(p_i-\bar p)(q_i-\bar q)^T$, then $H=U\Sigma V^T\Rightarrow R=V\,\mathrm{diag}(1,1,d)\,U^T$ <br> <i>Hint:</i> compute both centroids first, center the points, then form the two outer products and add. :: A: centroids $\bar p=(0.5,0,0)$, $\bar q=(0,0.5,0)$ <br> centered source $(-0.5,0,0),(0.5,0,0)$; centered target $(0,-0.5,0),(0,0.5,0)$ <br> $H=(-0.5,0,0)^T(0,-0.5,0)+(0.5,0,0)^T(0,0.5,0)=\begin{bmatrix}0&0.5&0\\0&0&0\\0&0&0\end{bmatrix}$ <br> $H$ maps $x\to y$, so $\boxed{R=\begin{bmatrix}0&-1&0\\1&0&0\\0&0&1\end{bmatrix}}$ the $+90°$ $z$-rotation
- Q: Worked example (condition number, near-singular). Singular values are $(1,\,10^{-6},\,1)$. Compute $\kappa_2$ and the digits lost. <br> <b>Formula:</b> $\kappa_2=\sigma_{\max}/\sigma_{\min}$; digits lost $\approx\log_{10}\kappa_2$ <br> <i>Hint:</i> pick out the largest and smallest singular values first, then take the ratio. :: A: $\sigma_{\max}=1$, $\sigma_{\min}=10^{-6}$ <br> $\kappa_2=1/10^{-6}=10^{6}$ <br> digits lost $\approx\log_{10}(10^6)=6$ <br> $\boxed{\kappa_2=10^{6},\ \approx 6\text{ digits lost}}$
- Q: Closing test (rigid solve). From the $90°$ example you found $t=\bar q-R\bar p$. Verify the source lands on the target. <br> <b>Formula:</b> $t=\bar q-R\bar p$, then check $Rp_2+t=q_2$ <br> <i>Hint:</i> a pass means $t=0$ here and $Rp_2+t$ equals $q_2=(0,1,0)$ exactly. :: A: $R\bar p=\begin{bmatrix}0&-1&0\\1&0&0\\0&0&1\end{bmatrix}(0.5,0,0)^T=(0,0.5,0)$ <br> $t=\bar q-R\bar p=(0,0.5,0)-(0,0.5,0)=(0,0,0)$ <br> check $p_2$: $R(1,0,0)^T+t=(0,1,0)=q_2$ <br> $\boxed{t=0\text{ and correspondences match}}$
- Q: Closing test (degenerate neighborhood). A perfectly collinear set of points: find the covariance rank and confirm the normal is undefined. <br> <b>Formula:</b> trust guard $\lambda_1\gg\lambda_0$; rank $=$ number of nonzero eigenvalues of $C$ <br> <i>Hint:</i> a pass is the guard FIRING (refusing the normal) because two eigenvalues are zero. :: A: all points lie on one line, so spread exists in one direction only <br> rank $C=1$: eigenvalues $\lambda_2>0,\ \lambda_1=\lambda_0=0$ <br> two zero eigenvalues = a whole plane of tied "least-spread" eigenvectors <br> guard check: $\lambda_1\not\gg\lambda_0$ fires <br> $\boxed{\text{normal undefined; reject or enlarge the neighborhood}}$
- Q: Closing test (reflection check). You skipped the $\det$ correction on mirrored data and got $\det R=-1$. Verify it is a reflection and fix it. <br> <b>Formula:</b> proper rotation needs $\det R=+1$; fix via $R=V\,\mathrm{diag}(1,1,d)\,U^T$, $d=\det(VU^T)$ <br> <i>Hint:</i> a pass is $\det R=+1$ after inserting $d=-1$ on the last diagonal entry. :: A: $\det R=-1\neq +1$, so $R\notin SO(3)$: it mirrors geometry <br> recompute $d=\det(VU^T)=-1$ <br> apply $R=V\,\mathrm{diag}(1,1,-1)\,U^T$ <br> now $\det R=\det V\cdot(1\cdot1\cdot d)\cdot\det U^T=+1$ <br> $\boxed{\det R=+1\text{ restores a proper rotation}}$
- Q: Recall: when do the SVD and the eigendecomposition give the same factorization? <br> <i>Hint:</i> think which matrices have $U=V$. :: A: for a symmetric positive semidefinite matrix the two coincide <br> $\boxed{U=V=Q,\ \sigma_i=\lambda_i}$
- Q: Pitfall: when is a PCA-estimated normal untrustworthy? <br> <b>Formula:</b> trust requires $\lambda_1\gg\lambda_0$ <br> <i>Hint:</i> think about collinear or coplanar neighborhoods where two eigenvalues nearly tie. :: A: when the neighborhood is degenerate (collinear/coplanar) and $\lambda_1$ is not $\gg\lambda_0$ <br> the smallest two eigenvalues nearly tie, leaving the least-spread direction ambiguous <br> $\boxed{\text{require }\lambda_1\gg\lambda_0\text{ before trusting }n}$
- CLOZE: For a symmetric positive semidefinite matrix the SVD and eigendecomposition coincide: $U=V=Q$ and $\sigma_i={{c1::\lambda_i}}$.
- CLOZE: `Eigen::SelfAdjointEigenSolver` returns eigenvalues in {{c1::ascending}} order, so the normal is column $0$. <br> <i>hint:</i> smallest first.
- CLOZE: Eigenvalue "zero" thresholds must be {{c1::relative}} ($\varepsilon\,\lambda_{\max}$), because covariance entries scale with the square of the length unit. <br> <i>hint:</i> a fixed absolute cutoff breaks when you switch m to mm.
- CLOZE: Use {{c1::SelfAdjointEigenSolver}} for symmetric covariance matrices and {{c2::JacobiSVD}} for the small $3\times3$ Kabsch SVD.
