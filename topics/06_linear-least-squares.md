# Linear least squares

Prereqs: [[01_matrices]], [[05_eig-svd-pca]] => this => [[07_nonlinear-least-squares]], [[08_calibration-pnp]], [[09_registration-icp]]

## Why it matters for the stone scan->CAD->IFC / CADFit workflow

A stone scan gives you thousands of noisy 3D points and almost no clean equations. Every step that turns those points into CAD asks the same question: which simple model best explains an overdetermined pile of measurements? Fitting a plane to a sawn face, a line to a drill axis, a circle to a coring mark, or the best rigid transform that snaps a scan onto a reference block are all linear (or linearized) least squares problems. The Kabsch-Umeyama rigid alignment that seats a raw scan into a CAD frame is a least squares solve dressed up with an SVD, and it is the closed-form core inside ICP registration. When CADFit-style program synthesis proposes a primitive (an extrusion, a cylinder, a fillet radius), the fit-and-validate inner loop measures how well that primitive matches the mesh, and the cheapest honest fit is least squares. Get this topic right and plane-fit, axis-fit, calibration, and registration all become the same three lines of Eigen with different matrices.

## Intuition

You have more equations than unknowns and they disagree, because measurements are noisy. You cannot satisfy all of them, so you pick the unknowns that come closest in total squared error.

Picture it geometrically. The matrix $A$ has $n$ columns. All vectors of the form $Ax$ live in a flat subspace, the column space of $A$. Your data vector $b$ usually sticks out of that subspace. Least squares drops a perpendicular from $b$ onto the subspace. The foot of that perpendicular is the closest reachable point $A\hat x$, and the leftover stick-out is the residual.

Analogy: you stand off a straight wall (the column space) and want the nearest point on it. You do not walk diagonally to some random spot. You walk straight at the wall along the perpendicular. The perpendicular direction is exactly the condition $A^\top(b - A\hat x) = 0$: the residual is orthogonal to every column of $A$. That one orthogonality fact generates the entire theory.

## The math (first principles)

### Setup and conventions

We solve the overdetermined system

$$A x \approx b, \qquad A \in \mathbb{R}^{m \times n}, \quad m \ge n, \quad x \in \mathbb{R}^n, \quad b \in \mathbb{R}^m.$$

Each row of $A$ and the matching entry of $b$ is one measurement equation. $x$ holds the model parameters. Conventions used throughout:

- Vectors are columns. $A^\top$ is transpose. Norm is the Euclidean 2-norm $\lVert v \rVert = \sqrt{v^\top v}$.
- Units live in $A$ and $b$, not in the math. If columns of $A$ have wildly different units (millimetres next to a constant 1), the problem is artificially ill-conditioned. Normalize or use weights.
- We assume real data. No complex conjugates.
- "Full column rank" means $\operatorname{rank}(A) = n$, equivalently the $n$ columns are linearly independent. This is the non-degenerate case.

The objective is the sum of squared residuals:

$$\hat x = \arg\min_{x} \; \lVert A x - b \rVert^2 = \arg\min_x \sum_{i=1}^{m} \big( a_i^\top x - b_i \big)^2,$$

where $a_i^\top$ is row $i$ of $A$.

### Normal equations from orthogonality

Expand the objective $f(x) = \lVert Ax - b \rVert^2 = x^\top A^\top A x - 2 x^\top A^\top b + b^\top b$. It is a convex quadratic. Set the gradient to zero:

$$\nabla f(x) = 2 A^\top A x - 2 A^\top b = 0 \;\Longrightarrow\; \boxed{A^\top A\, x = A^\top b}.$$

These are the **normal equations**. The name comes from the geometry: rearranged, they say $A^\top (b - A\hat x) = 0$, i.e. the residual $r = b - A\hat x$ is **normal (orthogonal)** to the column space of $A$. If $A$ has full column rank, $A^\top A$ is symmetric positive definite and invertible, so the solution is unique:

$$\hat x = (A^\top A)^{-1} A^\top b.$$

### Geometric meaning: orthogonal projection

Substitute $\hat x$ back. The best reachable point is

$$\hat b = A \hat x = \underbrace{A (A^\top A)^{-1} A^\top}_{P} \, b.$$

$P$ is the **orthogonal projection matrix** (the "hat matrix") onto $\operatorname{Col}(A)$. It satisfies $P = P^\top$ and $P^2 = P$. Least squares literally projects $b$ onto the column space; the residual is the rejected component in the orthogonal complement. This is why the answer is the closest point: projection minimizes distance.

### Solving via QR (numerically stable, recommended default)

Forming $A^\top A$ squares the conditioning (see below), so prefer a factorization of $A$ itself. The thin QR decomposition writes $A = QR$ with $Q \in \mathbb{R}^{m \times n}$ having orthonormal columns ($Q^\top Q = I$) and $R \in \mathbb{R}^{n \times n}$ upper triangular. Then

$$\lVert Ax - b \rVert^2 = \lVert QRx - b \rVert^2 = \lVert Rx - Q^\top b \rVert^2 + \text{const},$$

because multiplying by an orthonormal matrix preserves length. Minimize by solving the small triangular system

$$R \hat x = Q^\top b \;\Longrightarrow\; \hat x = R^{-1} Q^\top b,$$

solved cheaply by back-substitution. Column-pivoted Householder QR handles rank-deficient and near-degenerate $A$ gracefully.

### Solving via SVD and the pseudoinverse (most robust)

Let $A = U \Sigma V^\top$ be the (thin) SVD: $U \in \mathbb{R}^{m \times n}$ and $V \in \mathbb{R}^{n \times n}$ have orthonormal columns, $\Sigma = \operatorname{diag}(\sigma_1, \dots, \sigma_n)$ with $\sigma_1 \ge \dots \ge \sigma_n \ge 0$. The **Moore-Penrose pseudoinverse** is

$$A^{+} = V \Sigma^{+} U^\top, \qquad \Sigma^{+} = \operatorname{diag}\!\left( \tfrac{1}{\sigma_1}, \dots, \tfrac{1}{\sigma_r}, 0, \dots, 0 \right),$$

where only nonzero singular values (above a tolerance) are inverted; the rest map to $0$. The least squares solution is

$$\hat x = A^{+} b = \sum_{i=1}^{r} \frac{u_i^\top b}{\sigma_i}\, v_i.$$

When $A$ is rank-deficient there are infinitely many minimizers; $A^+ b$ returns the unique one of **minimum norm** $\lVert x \rVert$. SVD is the gold standard when you suspect degeneracy, because you can see the small $\sigma_i$ and truncate them.

### Weighted least squares

If measurements have different reliabilities, weight them. With a symmetric positive definite weight matrix $W$ (often diagonal $W = \operatorname{diag}(w_1, \dots, w_m)$ with $w_i = 1/\sigma_i^2$ from per-measurement variances), minimize $\lVert Ax - b \rVert_W^2 = (Ax-b)^\top W (Ax-b)$. The gradient gives the **weighted normal equations**:

$$\boxed{A^\top W A\, x = A^\top W b}, \qquad \hat x = (A^\top W A)^{-1} A^\top W b.$$

Equivalently, factor $W = C^\top C$ (Cholesky), define $\tilde A = C A$, $\tilde b = C b$, and run ordinary least squares on $(\tilde A, \tilde b)$. Under the Gauss-Markov assumptions (zero-mean, uncorrelated, equal-variance errors), ordinary least squares is the best linear unbiased estimator (BLUE); with known unequal variances, the variance-weighted solution is BLUE.

### Conditioning: when normal equations fail

The sensitivity of the solution scales with the condition number $\kappa(A) = \sigma_1 / \sigma_n$. The fatal fact about normal equations is

$$\kappa(A^\top A) = \kappa(A)^2.$$

Forming $A^\top A$ **squares** the condition number. If $\kappa(A) \approx 10^8$ in double precision (machine epsilon $\approx 10^{-16}$), then $\kappa(A^\top A) \approx 10^{16}$ and the Cholesky factorization of $A^\top A$ loses all accuracy or fails outright. QR works directly on $A$ and keeps $\kappa(A)$; SVD is even safer because it exposes the singular values. Rule: use normal equations only when $A$ is small, well-conditioned, and speed matters; otherwise use QR; use SVD when you fear rank deficiency.

## Worked example

Fit a line $y = c_0 + c_1 x$ to three points: $(1, 1)$, $(2, 2)$, $(3, 2)$. Unknowns are $x = (c_0, c_1)^\top$.

$$A = \begin{bmatrix} 1 & 1 \\ 1 & 2 \\ 1 & 3 \end{bmatrix}, \qquad b = \begin{bmatrix} 1 \\ 2 \\ 2 \end{bmatrix}.$$

Form the normal equations $A^\top A\, x = A^\top b$:

$$A^\top A = \begin{bmatrix} 3 & 6 \\ 6 & 14 \end{bmatrix}, \qquad A^\top b = \begin{bmatrix} 5 \\ 11 \end{bmatrix}.$$

Check: top-left $= 1+1+1 = 3$; off-diagonal $= 1\cdot1 + 1\cdot2 + 1\cdot3 = 6$; bottom-right $= 1+4+9 = 14$. Right side: $\sum b_i = 5$; $\sum x_i b_i = 1\cdot1 + 2\cdot2 + 3\cdot2 = 11$.

Invert the $2\times2$ system. Determinant $= 3\cdot14 - 6\cdot6 = 42 - 36 = 6$.

$$\hat x = \frac{1}{6} \begin{bmatrix} 14 & -6 \\ -6 & 3 \end{bmatrix} \begin{bmatrix} 5 \\ 11 \end{bmatrix} = \frac{1}{6} \begin{bmatrix} 70 - 66 \\ -30 + 33 \end{bmatrix} = \frac{1}{6} \begin{bmatrix} 4 \\ 3 \end{bmatrix} = \begin{bmatrix} 2/3 \\ 1/2 \end{bmatrix}.$$

So the best-fit line is $y = \tfrac{2}{3} + \tfrac{1}{2} x$. Verify the residual is orthogonal to the columns:

$$r = b - A\hat x = \begin{bmatrix} 1 \\ 2 \\ 2 \end{bmatrix} - \begin{bmatrix} 7/6 \\ 5/3 \\ 13/6 \end{bmatrix} = \begin{bmatrix} -1/6 \\ 1/3 \\ -1/6 \end{bmatrix}.$$

Then $\mathbf{1}^\top r = -\tfrac16 + \tfrac13 - \tfrac16 = 0$ and $x^\top r = -\tfrac16 + \tfrac23 - \tfrac12 = 0$, so $A^\top r = 0$ as required. The minimum squared error is $\lVert r \rVert^2 = \tfrac{1}{36}(1 + 4 + 1) = \tfrac{6}{36} = \tfrac{1}{6}$.

## Code

Idiomatic C++ with Eigen. The same matrix `A` and vector `b` solved three ways so you can compare. Symbol map: `A` = $A$, `b` = $b$, `x` = $\hat x$, `r` = residual.

```cpp
// least_squares.cpp
// Build: g++ -std=c++17 -I /path/to/eigen least_squares.cpp -o ls
#include <Eigen/Dense>
#include <iostream>

int main() {
    // Fit y = c0 + c1 x to (1,1),(2,2),(3,2). Design matrix [1, x_i].
    Eigen::MatrixXd A(3, 2);
    A << 1, 1,
         1, 2,
         1, 3;
    Eigen::VectorXd b(3);
    b << 1, 2, 2;

    // --- Method 1: QR (the recommended default: stable, fast) ---
    // colPivHouseholderQr solves min ||A x - b|| directly on A (no A^T A).
    Eigen::VectorXd x_qr = A.colPivHouseholderQr().solve(b);

    // --- Method 2: SVD via pseudoinverse (most robust; exposes rank) ---
    // bdcSvd scales to large A and falls back to JacobiSVD for small A.
    Eigen::BDCSVD<Eigen::MatrixXd> svd(A,
        Eigen::ComputeThinU | Eigen::ComputeThinV);
    Eigen::VectorXd x_svd = svd.solve(b);  // x = A^+ b, minimum-norm if rank-deficient

    // --- Method 3: normal equations A^T A x = A^T b (small, well-conditioned only) ---
    Eigen::MatrixXd AtA = A.transpose() * A;          // squares the conditioning
    Eigen::VectorXd Atb = A.transpose() * b;
    Eigen::VectorXd x_ne = AtA.ldlt().solve(Atb);     // LDL^T for symmetric pos-def

    Eigen::VectorXd r = b - A * x_qr;                  // residual
    std::cout << "x_qr  = " << x_qr.transpose()  << "\n"   // expect 0.6667 0.5
              << "x_svd = " << x_svd.transpose() << "\n"
              << "x_ne  = " << x_ne.transpose()  << "\n"
              << "||r||^2 = " << r.squaredNorm() << "\n"    // expect 0.16667
              << "cond(A) = " << svd.singularValues()(0)
                               / svd.singularValues()(svd.singularValues().size()-1)
              << "\n";

    // Weighted least squares: w_i = 1/variance_i. Scale rows by sqrt(w), then OLS.
    Eigen::VectorXd w(3); w << 4, 1, 1;                // trust point 1 four times more
    Eigen::ArrayXd s = w.array().sqrt();
    Eigen::MatrixXd Aw = s.matrix().asDiagonal() * A;   // C A with C = diag(sqrt(w))
    Eigen::VectorXd bw = s.matrix().asDiagonal() * b;   // C b
    Eigen::VectorXd x_w = Aw.colPivHouseholderQr().solve(bw);
    std::cout << "x_weighted = " << x_w.transpose() << "\n";
    return 0;
}
```

For genuinely rank-deficient `A`, prefer `A.completeOrthogonalDecomposition().solve(b)` (rank-revealing, returns the minimum-norm solution). For SVD you can set a rank tolerance with `svd.setThreshold(tol)` so tiny singular values are treated as zero.

Tiny Python mirror for quick checks. `np.linalg.lstsq` uses SVD under the hood:

```python
import numpy as np
A = np.array([[1, 1], [1, 2], [1, 3]], float)
b = np.array([1, 2, 2], float)
x, residuals, rank, sv = np.linalg.lstsq(A, b, rcond=None)
print(x)            # [0.6667 0.5]  -> y = 2/3 + x/2
print(A.T @ (b - A @ x))   # ~ [0, 0]  residual orthogonal to columns
```

## Connections

- [[01_matrices]]: least squares is built from matrix-vector products, transpose, and inversion of $A^\top A$; you must read $A$ as a stack of measurement rows.
- [[05_eig-svd-pca]]: the SVD route, the pseudoinverse, and the condition number $\kappa(A) = \sigma_1/\sigma_n$ all come straight from singular values; PCA is least squares orthogonal-distance fitting in disguise.
- [[07_nonlinear-least-squares]]: Gauss-Newton and Levenberg-Marquardt solve a **linear** least squares problem $J^\top J\, \delta = -J^\top r$ at every iteration, so this topic is the inner kernel of the nonlinear one.
- [[08_calibration-pnp]]: camera intrinsics and pose estimation start from a linear (DLT) least squares solve, then refine nonlinearly.
- [[09_registration-icp]]: the Kabsch-Umeyama rigid alignment ($H = \sum \tilde p_i \tilde q_i^\top$, $H = U\Sigma V^\top$, $R = V \operatorname{diag}(1,1,\det(VU^\top)) U^\top$) is a least squares solve over rotations, run once per ICP iteration.

## Pitfalls and failure modes

- **Normal equations on ill-conditioned $A$.** Forming $A^\top A$ squares $\kappa$. With $\kappa(A) > 10^8$ in double precision, Cholesky fails or the answer is garbage. Default to QR; reach for SVD when degeneracy is possible.
- **Mixed units / unscaled columns.** A design column of constant $1$ next to a column of millimetre values in the thousands inflates $\kappa$ artificially. Center and scale columns (or work in a sane unit), then unscale the solution.
- **Rank deficiency / degenerate geometry.** Fitting a plane to nearly collinear points, or a line to a single cluster, makes $A^\top A$ singular. Normal equations explode; SVD/COD silently return the minimum-norm answer, which may be meaningless. Always check the smallest singular value before trusting the fit.
- **Frame inconsistency.** When the rows of $A$ mix points expressed in different coordinate frames (scan frame vs CAD frame vs machine frame), the fit is nonsense even though the math runs. Transform everything into one frame first, and fix handedness (right-handed) before solving.
- **Outliers.** Least squares squares the error, so one bad measurement dominates. A single mis-registered scan point can tilt your plane fit. Use weighted least squares to down-weight, or move to RANSAC / robust M-estimators (an IRLS loop on top of weighted least squares).
- **$W$ not positive definite.** Weighted least squares assumes $W \succ 0$. A zero or negative weight breaks the Cholesky factor $W = C^\top C$ and the BLUE guarantee.
- **Confusing minimum-residual with minimum-norm.** When $A$ is rank-deficient, many $x$ achieve the same residual. The pseudoinverse picks the smallest-norm one; if you needed a specific parameterization, add a constraint instead of trusting the default.

## Practice

1. By hand, fit the constant model $y = c$ (so $A$ is a column of ones) to $b = (2, 4, 6)$. Show the normal equation reduces to $c = \operatorname{mean}(b) = 4$, and confirm the residual sums to zero. This proves least squares with a constant model is just the average.
2. In Eigen, fit a plane $z = a x + b y + c$ to at least four non-coplanar 3D points. Print the residual norm and the three singular values of the design matrix. Then make three points nearly collinear and watch the smallest singular value collapse.
3. Take a well-conditioned $5 \times 3$ system, solve it three ways (normal equations, QR, SVD), and print the answers. Then multiply one column by $10^{-9}$ and rerun. Report which method first produces a visibly wrong answer and relate it to $\kappa(A)^2$.
4. Implement weighted least squares two ways: directly via $(A^\top W A)^{-1} A^\top W b$, and via row-scaling by $\sqrt{w_i}$ followed by ordinary QR. Confirm the two solutions match to machine precision.
5. Generate a clean line fit, then corrupt one of twenty points to a gross outlier. Compare the ordinary least squares slope to a weighted fit where the outlier gets weight $0.01$. Quantify how far the outlier pulled the unweighted line.
6. (Bridge to registration.) Given two small matched point sets $\{p_i\}$ and $\{q_i\}$, build $H = \sum (p_i - \bar p)(q_i - \bar q)^\top$, take its SVD, and recover $R$ and $t$. Verify $R^\top R = I$ and $\det R = +1$. This is the least squares rigid solve you will reuse in ICP.

## References

- Gilbert Strang, "Introduction to Linear Algebra," chapters on projections and least squares. https://math.mit.edu/~gs/linearalgebra/
- Trefethen and Bau, "Numerical Linear Algebra" (SIAM, 1997), lectures 11 and 18-19 on least squares, QR, and conditioning. https://epubs.siam.org/doi/book/10.1137/1.9780898719574
- Golub and Van Loan, "Matrix Computations," 4th ed., chapter 5 (least squares, QR, SVD). https://jhupbooks.press.jhu.edu/title/matrix-computations
- Toby Driscoll, "Fundamentals of Numerical Computation," normal equations and conditioning. https://tobydriscoll.net/fnc-julia/leastsq/normaleqns.html
- Eigen documentation, "Solving linear least squares systems" (QR, SVD, normal equations idioms). https://libeigen.gitlab.io/eigen/docs-nightly/group__LeastSquares.html
- Eigen documentation, "Linear algebra and decompositions" tutorial. https://eigen.tuxfamily.org/dox/group__TutorialLinearAlgebra.html
- Cameron Musco, "Linear Least Squares, Projection, Pseudoinverses" lecture notes (UMass). https://people.cs.umass.edu/~cmusco/personal_site/pdfs/linear_least_squares.pdf
- Kabsch / Umeyama least squares rigid alignment: S. Umeyama, "Least-squares estimation of transformation parameters between two point patterns," IEEE TPAMI 13(4), 1991. https://web.stanford.edu/class/cs273/refs/umeyama.pdf
- NumPy `numpy.linalg.lstsq` reference (SVD-based least squares). https://numpy.org/doc/stable/reference/generated/numpy.linalg.lstsq.html
- Nehme, Whalen, Ahmed, "CADFit: Precise Mesh-to-CAD Program Generation with Hybrid Optimization," arXiv:2605.01171. https://arxiv.org/abs/2605.01171 and code https://github.com/ghadinehme/CADFit

## Flashcard seeds

- Q: Derive the normal equations, step 1 of 2: turn the orthogonality fact into a matrix equation. <br> <b>Formula:</b> residual $r=b-A\hat x$ obeys $a_j^\top r=0$ for every column $a_j$ of $A$; target $A^\top(b-A\hat x)=0$ <br> <i>Hint:</i> stack the $n$ scalar conditions $a_j^\top r=0$ as the rows of one matrix product. :: A: Each column gives $a_j^\top r=0$ <br> stacking all $j$ as rows is exactly $A^\top r=0$ <br> substitute $r=b-A\hat x$ <br> $\boxed{A^\top(b - A\hat x)=0}$
- Q: Derive the normal equations and closed form, step 2 of 2, from the orthogonality result. <br> <b>Formula:</b> start $A^\top(b-A\hat x)=0$; target $\hat x=(A^\top A)^{-1}A^\top b$ <br> <i>Hint:</i> distribute $A^\top$ across the bracket, move the $\hat x$ term to the other side, then left-multiply by $(A^\top A)^{-1}$. :: A: Distribute: $A^\top b - A^\top A\,\hat x = 0$ <br> rearrange: $A^\top A\,\hat x = A^\top b$ <br> full column rank $\Rightarrow A^\top A$ is SPD hence invertible <br> $\boxed{\hat x=(A^\top A)^{-1}A^\top b}$
- Q: Derive the normal equations by calculus: minimize $f(x)=\lVert Ax-b\rVert^2$. <br> <b>Formula:</b> $f(x)=x^\top A^\top A x - 2x^\top A^\top b + b^\top b$; target $A^\top A\,x = A^\top b$ <br> <i>Hint:</i> differentiate the quadratic, then set $\nabla f=0$. :: A: $\nabla f = 2A^\top A x - 2A^\top b$ <br> set $\nabla f=0$ <br> divide by 2 and rearrange <br> $\boxed{A^\top A\,x = A^\top b}$
- Q: Derive the orthogonal projection matrix $P$ that sends $b$ to $\hat b=A\hat x$. <br> <b>Formula:</b> $\hat x=(A^\top A)^{-1}A^\top b$; target $\hat b=Pb$ <br> <i>Hint:</i> substitute $\hat x$ into $\hat b=A\hat x$, then read off everything that multiplies $b$. :: A: $\hat b = A\hat x = A(A^\top A)^{-1}A^\top b$ <br> the operator acting on $b$ is the bracketed product <br> $\boxed{P=A(A^\top A)^{-1}A^\top}$
- Q: Show the projection matrix is idempotent ($P^2=P$). <br> <b>Formula:</b> $P=A(A^\top A)^{-1}A^\top$; target $P^2=P$ <br> <i>Hint:</i> write out $P\cdot P$ and collapse the central $A^\top A(A^\top A)^{-1}$. :: A: $P^2=A(A^\top A)^{-1}\,\underbrace{A^\top A(A^\top A)^{-1}}_{=\,I}\,A^\top$ <br> middle factor is $I$ <br> $P^2=A(A^\top A)^{-1}A^\top$ <br> $\boxed{P^2=P}$ (projecting twice changes nothing)
- Q: Derive the QR least squares solution, step 1 of 2: show $\lVert Ax-b\rVert^2=\lVert Rx-Q^\top b\rVert^2$ for $A=QR$. <br> <b>Formula:</b> $A=QR$ with $Q^\top Q=I$; orthonormal $Q$ preserves length, $\lVert Qv\rVert=\lVert v\rVert$ <br> <i>Hint:</i> substitute $A=QR$, then multiply the inside by $Q^\top$ (length-preserving). :: A: $\lVert Ax-b\rVert^2=\lVert QRx-b\rVert^2$ <br> apply $Q^\top$ inside (preserves length): $=\lVert Q^\top(QRx-b)\rVert^2$ <br> $Q^\top Q=I$ <br> $\boxed{\lVert Ax-b\rVert^2=\lVert Rx-Q^\top b\rVert^2}$
- Q: Finish the QR solve, step 2 of 2, from $\lVert Rx-Q^\top b\rVert^2$. <br> <b>Formula:</b> minimize $\lVert Rx-Q^\top b\rVert^2$; target $\hat x=R^{-1}Q^\top b$ <br> <i>Hint:</i> for full-rank $A$ this norm can be driven to zero, so set the bracket itself to zero and invert the triangular $R$. :: A: The norm is minimized when $Rx=Q^\top b$ (drives it to $0$) <br> $R$ is square upper-triangular, invertible for full-rank $A$ <br> solve by back-substitution <br> $\boxed{\hat x=R^{-1}Q^\top b}$
- Q: Derive the SVD / pseudoinverse solution from $A=U\Sigma V^\top$. <br> <b>Formula:</b> $A=U\Sigma V^\top$ with $U^\top U=I$; target $\hat x=A^+b=\sum_{i=1}^r \tfrac{u_i^\top b}{\sigma_i}v_i$ <br> <i>Hint:</i> substitute the SVD, set $y=V^\top x$, then minimize $\lVert\Sigma y-U^\top b\rVert^2$ one coordinate at a time. :: A: $\lVert Ax-b\rVert^2=\lVert \Sigma y-U^\top b\rVert^2$ with $y=V^\top x$ <br> per-coordinate min: $y_i=(u_i^\top b)/\sigma_i$ for $\sigma_i>0$, else $0$ <br> map back $x=Vy$ <br> $\boxed{\hat x=A^+b=\sum_{i=1}^r \tfrac{u_i^\top b}{\sigma_i}v_i}$
- Q: Derive the weighted normal equations, step 1 of 2: gradient of $f(x)=(Ax-b)^\top W(Ax-b)$. <br> <b>Formula:</b> $f=x^\top A^\top W A x - 2x^\top A^\top W b + b^\top W b$ (using $W=W^\top$); target $\nabla f=2A^\top W A x - 2A^\top W b$ <br> <i>Hint:</i> differentiate term by term, treating $A^\top W A$ like the symmetric matrix in $x^\top M x$. :: A: $\tfrac{d}{dx}(x^\top A^\top W A x)=2A^\top W A x$ <br> $\tfrac{d}{dx}(-2x^\top A^\top W b)=-2A^\top W b$ <br> constant term drops <br> $\boxed{\nabla f = 2A^\top W A x - 2A^\top W b}$
- Q: Finish the weighted solve, step 2 of 2, from $\nabla f=2A^\top W A x-2A^\top W b$. <br> <b>Formula:</b> target $\hat x=(A^\top W A)^{-1}A^\top W b$ <br> <i>Hint:</i> set $\nabla f=0$ and left-multiply by $(A^\top W A)^{-1}$, which exists because $W\succ0$. :: A: Set $\nabla f=0$: $A^\top W A\,x = A^\top W b$ <br> $W\succ0$ makes $A^\top W A$ SPD (full-rank $A$) <br> invert <br> $\boxed{\hat x=(A^\top W A)^{-1}A^\top W b}$
- Q: Show weighted least squares is ordinary LS on whitened data. <br> <b>Formula:</b> factor $W=C^\top C$ (Cholesky); target $\min\lVert \tilde A x-\tilde b\rVert^2$ with $\tilde A=CA,\ \tilde b=Cb$ <br> <i>Hint:</i> substitute $W=C^\top C$ into $(Ax-b)^\top W(Ax-b)$ and recognise $v^\top C^\top C v=\lVert Cv\rVert^2$. :: A: $(Ax-b)^\top C^\top C(Ax-b)=\lVert C(Ax-b)\rVert^2$ <br> $=\lVert \tilde A x-\tilde b\rVert^2$ with $\tilde A=CA,\ \tilde b=Cb$ <br> $\boxed{\min\lVert \tilde A x-\tilde b\rVert^2}$ (row-scale by $\sqrt{w_i}$ for diagonal $W$)
- Q: Derive that the condition number squares: $\kappa(A^\top A)=\kappa(A)^2$. <br> <b>Formula:</b> $A=U\Sigma V^\top$, $\kappa(A)=\sigma_1/\sigma_n$; target $\kappa(A^\top A)=\kappa(A)^2$ <br> <i>Hint:</i> compute $A^\top A$ from the SVD, use $U^\top U=I$ to collapse it to $V\Sigma^2 V^\top$, then read off its singular values. :: A: $A^\top A = V\Sigma U^\top U\Sigma V^\top = V\Sigma^2 V^\top$ <br> singular values of $A^\top A$ are $\sigma_i^2$ <br> $\kappa(A^\top A)=\sigma_1^2/\sigma_n^2=(\sigma_1/\sigma_n)^2$ <br> $\boxed{\kappa(A^\top A)=\kappa(A)^2}$
- Q: Worked example (small, 2D, unitless). Fit a line $y=c_0+c_1x$ to $(1,1),(2,2),(3,2)$ via the normal equations. <br> <b>Formula:</b> $\hat x=(A^\top A)^{-1}A^\top b$, with $A=\big[\begin{smallmatrix}1&1\\1&2\\1&3\end{smallmatrix}\big]$, $b=(1,2,2)^\top$ <br> <i>Hint:</i> build the $2\times2$ matrix $A^\top A$ and vector $A^\top b$ first, then invert with $\det=ad-bc$. :: A: $A^\top A=\begin{bmatrix}3&6\\6&14\end{bmatrix}$, $A^\top b=\begin{bmatrix}5\\11\end{bmatrix}$, $\det=42-36=6$ <br> $\hat x=\tfrac16\begin{bmatrix}14&-6\\-6&3\end{bmatrix}\begin{bmatrix}5\\11\end{bmatrix}=\tfrac16\begin{bmatrix}4\\3\end{bmatrix}$ <br> $\boxed{y=\tfrac23+\tfrac12 x}$
- Q: Worked example (degenerate model). Fit the constant model $y=c$ (so $A=\mathbf 1$, a column of ones) to $b=(2,4,6)$. <br> <b>Formula:</b> normal equation $\mathbf 1^\top\mathbf 1\,c=\mathbf 1^\top b$ <br> <i>Hint:</i> compute $\mathbf 1^\top\mathbf 1$ (the count) and $\mathbf 1^\top b$ (the sum) first, then divide. :: A: $\mathbf 1^\top\mathbf 1 = 3$, $\mathbf 1^\top b = 2+4+6=12$ <br> $3c=12$ <br> $\boxed{c=\operatorname{mean}(b)=4}$ (least squares with a constant model is the average)
- Q: Worked example (3D, metres, stone-scan plane fit). 4 mesh points on a sawn face fit $z=ax+by+c$, design rows $[x_i,y_i,1]$. Give the solve and the trust check. <br> <b>Formula:</b> $\hat x=(A^\top A)^{-1}A^\top z$ for $(a,b,c)$; trust gate is the smallest singular value $\sigma_{\min}(A)$ <br> <i>Hint:</i> a pass is $\sigma_{\min}$ well above tolerance; collinear points send $\sigma_{\min}\to0$. :: A: Solve $\hat x=(A^\top A)^{-1}A^\top z$ to get $(a,b,c)$ <br> compute the SVD of $A$ and inspect $\sigma_{\min}$ <br> $\boxed{(a,b,c)\ \text{trustworthy} \iff \sigma_{\min}\gg0}$ (near-collinear points collapse $\sigma_{\min}$)
- Q: Worked example (mm vs m scaling pitfall). A constant-$1$ column sits next to an $x$-column in millimetres near $3000$. Explain the ill-conditioning and the fix. <br> <b>Formula:</b> $\kappa(A^\top A)=\kappa(A)^2$; column norms differ by $\sim10^3$ <br> <i>Hint:</i> estimate how the mismatched column scales inflate $\kappa(A)$, then square it for the normal equations. :: A: Column norms differ by $\sim10^3$, so $\kappa(A)$ inflates by $\sim10^3$ <br> normal equations square it: $\kappa(A^\top A)$ worse by $\sim10^6$ <br> center and scale each column (or work in metres), solve, then unscale <br> $\boxed{\text{scale columns} \Rightarrow \kappa(A)\downarrow}$
- Q: Closing test (orthogonality check). For $y=\tfrac23+\tfrac12x$ on $(1,1),(2,2),(3,2)$, compute the residual $r$ and confirm $A^\top r=0$. <br> <b>Formula:</b> $r=b-A\hat x$; pass condition $A^\top r=0$ (both $\mathbf 1^\top r=0$ and $x^\top r=0$) <br> <i>Hint:</i> a pass means both column dot-products vanish; compute $\hat b=A\hat x$ first. :: A: $\hat b=(7/6,5/3,13/6)$, so $r=(-1/6,1/3,-1/6)$ <br> $\mathbf 1^\top r=-\tfrac16+\tfrac13-\tfrac16=0$ <br> $x^\top r=-\tfrac16+\tfrac23-\tfrac12=0$ <br> $\boxed{A^\top r=0\ \checkmark}$ and $\lVert r\rVert^2=\tfrac16$
- Q: Closing test (projection fixes the column space). Check that $P$ leaves any vector already in $\operatorname{Col}(A)$ unchanged, e.g. $b=Av$. <br> <b>Formula:</b> $P=A(A^\top A)^{-1}A^\top$; pass condition $Pb=b$ when $b=Av$ <br> <i>Hint:</i> a pass means $Pb=b$ with residual $0$; substitute $b=Av$ and cancel $A^\top A(A^\top A)^{-1}=I$. :: A: $Pb=A(A^\top A)^{-1}A^\top (Av)$ <br> $A^\top A v$ cancels with $(A^\top A)^{-1}$ <br> $Pb=Av=b$ <br> $\boxed{Pb=b\ \text{for }b\in\operatorname{Col}(A)\ \checkmark}$ (nothing to project, residual $0$)
- Q: When-to-use: pick the solver for small well-conditioned vs ill-conditioned vs suspected rank-deficient $A$. <br> <i>Hint:</i> the deciding factor is whether you can afford to square $\kappa$ and whether degeneracy is possible. :: A: Small + well-conditioned + speed: normal equations. <br> General default: QR (keeps $\kappa(A)$, no squaring). <br> Suspected degeneracy: SVD/COD (exposes and truncates small $\sigma_i$).
- Q: Pitfall: when do the normal equations actively fail, and why? <br> <b>Formula:</b> $\kappa(A^\top A)=\kappa(A)^2$, double-precision $\epsilon\approx10^{-16}$ <br> <i>Hint:</i> ask what $\kappa(A)^2$ becomes once $\kappa(A)$ passes $10^8$. :: A: Forming $A^\top A$ squares $\kappa$ <br> $\kappa(A)>10^8\Rightarrow\kappa(A^\top A)>10^{16}$, beyond double precision <br> Cholesky loses all accuracy or fails <br> $\boxed{\text{use QR or SVD instead}}$
- Q: Pitfall: for rank-deficient $A$, distinguish minimum-residual from minimum-norm. <br> <b>Formula:</b> pseudoinverse solution $\hat x=A^+b$ <br> <i>Hint:</i> note many $x$ tie on residual; recall which single one $A^+$ selects. :: A: Many $x$ achieve the same minimal residual <br> $A^+b$ returns the unique smallest-$\lVert x\rVert$ one <br> if you needed a specific parameterization, add a constraint rather than trust the default <br> $\boxed{A^+b=\text{min-norm minimizer}}$
- Q: Pitfall: why is least squares fragile to outliers, and what is the fix? <br> <b>Formula:</b> objective $\sum_i (a_i^\top x-b_i)^2$ squares each residual <br> <i>Hint:</i> see what squaring does to one very large $|a_i^\top x-b_i|$. :: A: Squaring lets one large error dominate the sum, tilting the fit <br> down-weight via weighted least squares ($w_i\downarrow$ for bad rows) <br> or move to RANSAC / robust M-estimators (IRLS on top of WLS) <br> $\boxed{\text{down-weight or go robust}}$
- CLOZE: The normal equations express that the residual is {{c1::orthogonal}} to the column space of $A$, i.e. $A^\top(b-A\hat x)=0$.
- CLOZE: Forming $A^\top A$ squares the conditioning: $\kappa(A^\top A)={{c1::\kappa(A)^2}}$. <br> <i>hint:</i> singular values of $A^\top A$ are $\sigma_i^2$.
- CLOZE: Under the Gauss-Markov assumptions, ordinary least squares is the {{c1::best linear unbiased estimator}} (BLUE).
- CLOZE: For rank-deficient $A$, the pseudoinverse solution $A^+b$ is the unique least squares minimizer of smallest {{c1::norm}} $\lVert x\rVert$.
- CLOZE: The Kabsch-Umeyama rigid alignment forms $H=\sum\tilde p_i\tilde q_i^\top$, takes its SVD $H=U\Sigma V^\top$, and sets $R=V\operatorname{diag}(1,1,{{c1::\det(VU^\top)}})U^\top$ to guarantee $\det R=+1$.
- Q: Eigen idiom for the recommended (QR) least squares solve? <br> <i>Hint:</i> the rank-revealing Householder QR method on $A$, then `.solve(b)`. :: A: `A.colPivHouseholderQr().solve(b)`.
- Q: Eigen idiom for the SVD-based least squares solve? <br> <i>Hint:</i> the divide-and-conquer SVD with thin $U$ and $V$, then `.solve(b)`. :: A: `A.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(b)`.
