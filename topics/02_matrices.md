# Matrices and linear systems

Prereqs: [[01_vectors]] => this => [[03_eig-svd-pca]], [[04_transforms-se3]], [[05_linear-least-squares]].

## Why it matters for the stone scan->CAD->IFC / CADFit workflow

Every coordinate you touch in the stone pipeline lives or dies by a matrix. A 3x3 rotation aligns a raw LiDAR scan of an irregular block to the cut template; a 4x4 homogeneous transform places that block in world space and again in the robot base frame and again in the IFC `IfcLocalPlacement`. PCA on a stone patch is an eigenproblem of a 3x3 covariance matrix, and the local frame it returns is literally three matrix columns. Fitting a plane, a cut line, or a B-spline control net to scattered scan points reduces to solving a linear system $Ax=b$ (or its least-squares cousin), so Gaussian elimination and the conditioning of $A$ decide whether the fit is stable or garbage. A CADFit-style program-synthesis layer that proposes an extrusion or revolution must compose and invert these transforms to test each candidate against the mesh, so matrix fluency is the floor you build the whole stack on.

## Intuition

A matrix is a machine that eats a vector and spits out a vector, and the only thing it is allowed to do is stretch, rotate, shear, reflect, and squash space uniformly. The whole machine is described by where it sends the basis arrows. In 2D, ask only "where does $\hat{x}=(1,0)$ go, and where does $\hat{y}=(0,1)$ go?" Those two answer-arrows, written as columns, *are* the matrix. Every other point rides along: if a point was $3\hat{x}+2\hat{y}$ before, it is $3(\text{new }\hat{x})+2(\text{new }\hat{y})$ after.

Analogy: a matrix is a stretchy rubber grid stamped on a sheet. You only need to know where two grid corners land to know where the entire sheet lands, because the grid stays a grid (parallel lines stay parallel, evenly spaced). The determinant is then just "how much did one grid square's area change, and did the sheet get flipped over?"

## The math (first principles)

### Definitions and conventions

Conventions used throughout this document:

- **Column vectors.** $x \in \mathbb{R}^n$ is an $n\times 1$ column. A linear map acts by left-multiplication: $y = Ax$.
- **Row-major vs column-major is a storage detail, not math.** Eigen defaults to column-major storage; OpenCV `cv::Mat` is row-major. The math below is storage-agnostic.
- **Right-handed frames.** Rotation matrices in $SO(3)$ have $\det = +1$ and preserve handedness. A $\det = -1$ orthogonal matrix is a reflection and must not be used as a rotation (see Pitfalls).
- **Tolerances are relative.** "Zero" for a pivot or singular value means $|v| < \tau$ with $\tau = \varepsilon \cdot \text{scale}$, never an exact compare. Default to a scale-relative threshold.

A matrix $A \in \mathbb{R}^{m\times n}$ represents a **linear map** $T:\mathbb{R}^n \to \mathbb{R}^m$ satisfying $T(\alpha u + \beta v) = \alpha T(u) + \beta T(v)$.

### Columns are images of basis vectors

This is the single most important geometric fact. Let $e_j$ be the $j$-th standard basis vector. Then

$$A e_j = a_j \quad \text{(the $j$-th column of $A$).}$$

So the columns of $A$ are exactly the images of the basis vectors. Any input is a weighted sum of columns:

$$A x = \sum_{j=1}^{n} x_j\, a_j .$$

The set of all reachable outputs is the **column space** (the span of the columns), which is the image of the map. Verified against standard linear-algebra references that define the matrix of a linear map as the coordinates of the images of the basis vectors.

### Matrix-vector and matrix-matrix product

Matrix-vector, row-dot-column form:

$$(Ax)_i = \sum_{j=1}^{n} A_{ij}\, x_j .$$

Matrix-matrix product $C = AB$ with $A\in\mathbb{R}^{m\times k}$, $B\in\mathbb{R}^{k\times n}$:

$$C_{ij} = \sum_{l=1}^{k} A_{il} B_{lj}, \qquad C \in \mathbb{R}^{m\times n}.$$

The geometric reading: $AB$ means "do $B$ first, then $A$." Composition reads right to left, like function composition $A(B(x))$. This is why **matrix multiplication does not commute**: $AB \neq BA$ in general. Rotating-then-translating is not translating-then-rotating.

### Identity, transpose, inverse

- **Identity** $I$: $Ix = x$. Columns are the basis vectors themselves.
- **Transpose** $A^\top$: $(A^\top)_{ij} = A_{ji}$. Rules: $(AB)^\top = B^\top A^\top$ and $(A^\top)^\top = A$.
- **Inverse** $A^{-1}$ (square $A$ only): the unique matrix with $A A^{-1} = A^{-1} A = I$. It undoes the map. Rule: $(AB)^{-1} = B^{-1} A^{-1}$ (order reverses).

A square matrix is **invertible** iff its columns are linearly independent iff $\det A \neq 0$ iff it has full rank.

The $2\times 2$ inverse is worth memorizing:

$$A = \begin{bmatrix} a & b \\ c & d \end{bmatrix}, \qquad A^{-1} = \frac{1}{ad-bc}\begin{bmatrix} d & -b \\ -c & a \end{bmatrix}, \qquad \det A = ad-bc.$$

### Determinant as signed area / volume scale

The determinant is the factor by which the map scales signed area (2D) or signed volume (3D), and its sign encodes orientation:

- $|\det A|$ = area/volume scaling factor.
- $\det A > 0$ => orientation preserved (handedness kept).
- $\det A < 0$ => orientation reversed (a reflection is included).
- $\det A = 0$ => the map collapses space onto a lower dimension (singular, non-invertible).

In 2D, $\det\begin{bmatrix} a & b \\ c & d\end{bmatrix} = ad - bc$ is the signed area of the parallelogram spanned by the two columns. In 3D it is the signed volume of the parallelepiped (the scalar triple product of the columns). Key identities: $\det(AB) = \det(A)\det(B)$ and $\det(A^{-1}) = 1/\det(A)$. Verified against multiple geometric-determinant references.

### Rank and singularity

The **rank** of $A$ is the dimension of its column space, equivalently the number of linearly independent columns (= number of independent rows). For a square $n\times n$ matrix:

- $\operatorname{rank}(A) = n$ => full rank => invertible => $\det A \neq 0$.
- $\operatorname{rank}(A) < n$ => rank-deficient => singular => $\det A = 0$.

Rank is a *numerical* notion in floating point: you count singular values above a scale-relative threshold, not exact nonzeros.

### Conditioning

Even when $A$ is invertible, it can be *nearly* singular. The 2-norm **condition number** measures how much $Ax=b$ amplifies input error:

$$\kappa(A) = \|A\|_2 \, \|A^{-1}\|_2 = \frac{\sigma_{\max}(A)}{\sigma_{\min}(A)},$$

the ratio of largest to smallest singular value. $\kappa \approx 1$ is well-conditioned; $\kappa$ huge means the solve loses digits and the problem is ill-posed; $\kappa = \infty$ means singular. Confirmed against standard condition-number references.

### Solving $Ax=b$ via Gaussian elimination and LU

To solve $Ax = b$ for square nonsingular $A$, factor with **LU decomposition with partial pivoting**:

$$P A = L U,$$

where $P$ is a permutation (row swaps), $L$ is unit lower-triangular, $U$ is upper-triangular. Then solve in two cheap triangular sweeps:

$$L y = P b \quad (\text{forward substitution}), \qquad U x = y \quad (\text{back substitution}).$$

**Partial pivoting** swaps rows so the largest-magnitude available entry becomes the pivot. This keeps all multipliers $|\ell_{ij}| \le 1$ and prevents tiny pivots from blowing up roundoff. Gaussian elimination without pivoting is numerically unstable unless $A$ is diagonally dominant or symmetric positive definite. Cost is $O(n^3)$ to factor, $O(n^2)$ per right-hand side, so factor once and reuse for many $b$. Verified against Gaussian-elimination / GEPP numerical-stability references.

Rule of thumb: **never form $A^{-1}$ to solve a system.** Use a factorization and `solve`. The inverse is less accurate and slower.

## Worked example

Solve

$$A = \begin{bmatrix} 2 & 1 \\ 4 & 3 \end{bmatrix}, \quad b = \begin{bmatrix} 1 \\ 7 \end{bmatrix}, \qquad A x = b.$$

**Determinant.** $\det A = (2)(3) - (1)(4) = 6 - 4 = 2$. Nonzero, so invertible; the map scales area by a factor of 2 and preserves orientation.

**LU with partial pivoting.** Column-1 candidates are $2$ and $4$; the larger magnitude is $4$, so swap rows. $P$ swaps rows 1 and 2:

$$PA = \begin{bmatrix} 4 & 3 \\ 2 & 1 \end{bmatrix}, \qquad Pb = \begin{bmatrix} 7 \\ 1 \end{bmatrix}.$$

Eliminate: multiplier $\ell_{21} = 2/4 = 0.5$. Row2 $\leftarrow$ Row2 $- 0.5\cdot$Row1:

$$U = \begin{bmatrix} 4 & 3 \\ 0 & 1 - 0.5\cdot 3 \end{bmatrix} = \begin{bmatrix} 4 & 3 \\ 0 & -0.5 \end{bmatrix}, \qquad L = \begin{bmatrix} 1 & 0 \\ 0.5 & 1 \end{bmatrix}.$$

**Forward solve** $Ly = Pb$: $y_1 = 7$; $0.5\,y_1 + y_2 = 1 \Rightarrow y_2 = 1 - 3.5 = -2.5$.

**Back solve** $Ux = y$: $-0.5\,x_2 = -2.5 \Rightarrow x_2 = 5$; $4 x_1 + 3 x_2 = 7 \Rightarrow 4 x_1 = 7 - 15 = -8 \Rightarrow x_1 = -2$.

So $x = (-2,\ 5)$.

**Check.** $Ax = \begin{bmatrix} 2(-2)+1(5)\\ 4(-2)+3(5)\end{bmatrix} = \begin{bmatrix} 1 \\ 7\end{bmatrix} = b$. Correct.

**Cross-check via the 2x2 inverse.** $A^{-1} = \frac{1}{2}\begin{bmatrix} 3 & -1 \\ -4 & 2\end{bmatrix}$, so $x = A^{-1}b = \frac{1}{2}\begin{bmatrix} 3(1)-1(7)\\ -4(1)+2(7)\end{bmatrix} = \frac{1}{2}\begin{bmatrix} -4 \\ 10\end{bmatrix} = \begin{bmatrix} -2 \\ 5\end{bmatrix}$. Matches.

## Code

Idiomatic C++ with Eigen first. The comments map each call to the math symbol.

```cpp
// matrices_demo.cpp
// Build: g++ -O2 -I /path/to/eigen matrices_demo.cpp -o matrices_demo
#include <Eigen/Dense>
#include <iostream>

int main() {
    using Eigen::Matrix2d;
    using Eigen::Vector2d;

    // A is the linear map; its COLUMNS are the images of the basis vectors.
    // col(0) = A*e1, col(1) = A*e2.
    Matrix2d A;
    A << 2, 1,
         4, 3;
    Vector2d b(1, 7);

    // Columns-as-images check: A * e1 must equal A.col(0).
    Vector2d e1(1, 0);
    std::cout << "A*e1   = " << (A * e1).transpose() << "  (== col 0)\n";

    // det(A): signed area scale. >0 keeps orientation, ==0 means singular.
    double detA = A.determinant();           // det A = ad - bc = 2
    std::cout << "det(A) = " << detA << "\n";

    // Solve A x = b. For a general square invertible matrix, partialPivLu()
    // (LU with partial pivoting) is the fast, reliable default.
    // It does P A = L U internally, then forward+back substitution.
    Vector2d x = A.partialPivLu().solve(b);  // x = (-2, 5)
    std::cout << "x      = " << x.transpose() << "\n";

    // Residual: how well A x reproduces b. Should be ~machine epsilon.
    double residual = (A * x - b).norm();
    std::cout << "||Ax-b|| = " << residual << "\n";

    // Rank / singularity check on a POSSIBLY rank-deficient matrix.
    // colPivHouseholderQr is rank-revealing and works for any shape.
    Matrix2d S;
    S << 1, 2,
         2, 4;                                // row2 = 2*row1 => rank 1, singular
    auto qr = S.colPivHouseholderQr();
    std::cout << "rank(S) = " << qr.rank()
              << ", det(S) = " << S.determinant() << "\n";  // rank 1, det 0

    // Inverse exists for fixed-size <=4x4; Eigen uses a closed-form formula.
    // Prefer solve() over inverse() for actually solving systems.
    std::cout << "A^-1 =\n" << A.inverse() << "\n";

    // Non-commutativity: A*B != B*A in general.
    Matrix2d B;
    B << 0, -1,
         1,  0;                               // 90-degree rotation
    std::cout << "AB - BA =\n" << (A * B - B * A) << "\n";  // not zero
    return 0;
}
```

Solver-choice cheat sheet (Eigen, all expose `.solve(b)`):

- `A.partialPivLu().solve(b)` — general square, full rank. Fast default.
- `A.fullPivLu().solve(b)` — slowest, most robust; also gives `.rank()`, `.kernel()`.
- `A.colPivHouseholderQr().solve(b)` — any shape, rank-revealing, good compromise.
- `A.llt().solve(b)` — symmetric positive definite (e.g. a normal-equations matrix $A^\top A$). Fastest when it applies.

A short Python mirror, only to show the same three operations in the language you prototype in:

```python
import numpy as np
A = np.array([[2., 1.], [4., 3.]])
b = np.array([1., 7.])
x = np.linalg.solve(A, b)          # LU under the hood; never use inv() to solve
print(x, np.linalg.det(A), np.linalg.matrix_rank(A))  # [-2. 5.] 2.0 2
```

## Connections

- [[01_vectors]] — matrices act on vectors; a matrix is built from vectors as its columns. Prerequisite.
- [[03_eig-svd-pca]] — eigenvalues, SVD, and PCA all decompose a matrix; rank, condition number, and the rigid-fit covariance step are read off the SVD. Direct unlock.
- [[04_transforms-se3]] — rigid transforms are $4\times 4$ matrices $\begin{bmatrix} R & t \\ 0 & 1\end{bmatrix}$; composition is matrix multiplication and inversion undoes a pose. Direct unlock.
- [[05_linear-least-squares]] — overdetermined fits (plane, line, B-spline net) solve $A^\top A\,x = A^\top b$, a linear system whose conditioning this topic governs. Direct unlock.
- [[06_dot-cross-norm-projection]] — the dot product is the row-times-column atom inside every matrix product; projection is a matrix.

## Pitfalls and failure modes

- **Frame inconsistency.** Multiplying matrices expressed in different frames (scan frame vs robot base vs IFC placement) silently produces nonsense. Track which frame each matrix maps from and to; the product $A_{c\leftarrow w}\,A_{w\leftarrow s}$ must have matching inner subscripts.
- **Non-commutativity.** $AB \neq BA$. Rotate-then-translate differs from translate-then-rotate. Composition is right-to-left.
- **Ill-conditioning.** A large $\kappa(A)$ means the solve amplifies noise; a plane fit through nearly-collinear stone points has a near-singular normal matrix. Check $\kappa$ or the smallest singular value, do not trust the answer blindly.
- **Forming the inverse to solve.** `x = A.inverse() * b` is slower and less accurate than `A.partialPivLu().solve(b)`. Only form $A^{-1}$ when you genuinely need the matrix itself.
- **Exact zero tests.** Testing `det(A) == 0` or a pivot `== 0` fails in floating point. Use a scale-relative threshold; rely on a rank-revealing decomposition.
- **Reflection masquerading as rotation.** An orthogonal matrix with $\det = -1$ is a reflection, not a rotation. After an SVD-based rigid fit you must flip the last singular vector if $\det(VU^\top) < 0$, or the stone gets mirrored.
- **No pivoting.** Hand-rolled Gaussian elimination without pivoting divides by tiny pivots and explodes. Always pivot (or use a library that does).
- **Storage-order confusion.** Eigen is column-major, OpenCV `cv::Mat` is row-major. A raw `memcpy` between them transposes your matrix. Convert explicitly.
- **Dimension mismatch.** $AB$ needs `A.cols() == B.rows()`. Eigen catches this at compile time for fixed sizes, at runtime (assert) for dynamic; do not disable asserts in debug.

## Practice

1. **By hand.** Compute $\det$, then $A^{-1}$ for $A = \begin{bmatrix} 3 & 2 \\ 1 & 4\end{bmatrix}$. Verify $A A^{-1} = I$. Then solve $Ax = (5, 6)^\top$ two ways: via $A^{-1}b$ and via LU-by-hand with partial pivoting. They must agree.

2. **Columns as images.** In Eigen, build a $2\times 2$ matrix that rotates by $30^\circ$. Print `A * e1` and `A * e2` and confirm they equal `A.col(0)` and `A.col(1)`. Then print `A.determinant()` and explain why it is exactly $1$.

3. **Orientation flip.** Construct a $2\times 2$ reflection (e.g. $\operatorname{diag}(1,-1)$). Show its determinant is $-1$, apply it to the unit-square corners, and confirm the corners' winding order reverses.

4. **Conditioning experiment.** Build $A = \begin{bmatrix} 1 & 1 \\ 1 & 1+\epsilon\end{bmatrix}$ for $\epsilon = 10^{-1}, 10^{-4}, 10^{-8}$. Solve $Ax=b$ for a fixed $b$, perturb $b$ by $10^{-6}$, and tabulate how much $x$ changes. Relate the swing to $\kappa(A) = \sigma_{\max}/\sigma_{\min}$ from `A.jacobiSvd()`.

5. **Rank detection.** Generate 50 nearly-coplanar 3D points (a thin slab of stone). Form the $3\times 3$ covariance matrix, print `colPivHouseholderQr(...).rank()` with a sensible threshold, and confirm it drops to 2 when the slab is flattened. This is the seed of plane fitting.

6. **Solve reuse.** Factor one $A$ with `partialPivLu()` and solve it against 1000 different right-hand sides in a loop. Confirm refactoring per solve is the slow way and the factor-once pattern is $O(n^2)$ per solve.

## References

- Lynch & Park, *Modern Robotics: Mechanics, Planning, and Control* (book + course site). Assumes matrix and $SE(3)$ fluency early. https://hades.mech.northwestern.edu/index.php/Modern_Robotics
- OpenCV `calib3d` camera model documentation (pinhole model $s\,p = K[R\mid t]P_w$, where the matrices of this topic appear). https://docs.opencv.org/4.x/d9/d0c/group__calib3d.html
- Eigen, *Linear algebra and decompositions* tutorial (`partialPivLu`, `fullPivLu`, `colPivHouseholderQr`, `.rank()`, `.determinant()`, `.inverse()`). https://eigen.tuxfamily.org/dox/group__TutorialLinearAlgebra.html
- Stanford CS205A, *Linear Systems and the LU Decomposition* (GEPP, stability). https://graphics.stanford.edu/courses/cs205a-13-fall/assets/notes/chapter2.pdf
- 3.3 Partial Pivoting, Chasnov, *Numerical Methods* (LibreTexts). https://math.libretexts.org/Bookshelves/Applied_Mathematics/Numerical_Methods_(Chasnov)/03%3A_System_of_Equations/3.03%3A_Partial_Pivoting
- Math Insight, *Determinants and linear transformations* (signed area/volume, orientation). https://mathinsight.org/determinant_linear_transformation
- Moler, *What is the Condition Number of a Matrix?* (Cleve's Corner, MathWorks). https://blogs.mathworks.com/cleve/2017/07/17/what-is-the-condition-number-of-a-matrix/
- OpenCASCADE Reference Manual (B-spline fitting and transforms downstream of this topic). https://dev.opencascade.org/doc/overview/html/
- IfcOpenShell geometry-processing docs (4x4 placements in IFC). https://docs.ifcopenshell.org/
- Nehme, Whalen & Ahmed, *CADFit: Precise Mesh-to-CAD Program Generation with Hybrid Optimization*, arXiv:2605.01171. https://arxiv.org/abs/2605.01171 ; code https://github.com/ghadinehme/CADFit

## Flashcard seeds

- Q: Derive what column $j$ of a matrix $A$ must equal, the image of a basis vector. <br> <b>Formula:</b> start from $Ax=\sum_j x_j (Ae_j)$, target $a_j = Ae_j$. <br> <i>Hint:</i> substitute $x=e_j$ so every weight is $0$ except $x_j=1$. :: A: Set $x=e_j$, so all weights are $0$ except $x_j=1$ <br> Then $Ax=Ae_j$ picks out a single column <br> By the definition of the matrix-vector product, $Ae_j$ is the $j$-th column $a_j$ <br> $\boxed{a_j = A e_j}$ (each column is the image of a basis vector).
- Q: Derive the expansion of $Ax$ over the columns of $A$ (step 1 of 2). <br> <b>Formula:</b> start from $x=\sum_j x_j e_j$, use linearity $A(\alpha u+\beta v)=\alpha Au+\beta Av$. <br> <i>Hint:</i> apply $A$ to both sides, then move $A$ inside the sum and past each scalar. :: A: Apply $A$ to both sides: $Ax = A\big(\sum_j x_j e_j\big)$ <br> Linearity moves $A$ inside the sum and past the scalars: $Ax=\sum_j x_j (A e_j)$ <br> $\boxed{Ax = \sum_{j} x_j\,(A e_j)}$.
- Q: Continue: derive $Ax$ as a combination of columns (step 2 of 2). <br> <b>Formula:</b> start from $Ax=\sum_j x_j (Ae_j)$ and the column identity $Ae_j=a_j$. <br> <i>Hint:</i> substitute $Ae_j=a_j$ into the sum. :: A: Substitute the column identity $Ae_j=a_j$ <br> $Ax=\sum_j x_j a_j$ <br> So the output is a weighted sum of the columns, weights $=$ entries of $x$ <br> $\boxed{Ax = \sum_{j=1}^{n} x_j\,a_j}$ (the column space is the span of the $a_j$).
- Q: Derive the row-dot-column formula for the $i$-th entry $(Ax)_i$. <br> <b>Formula:</b> start from $Ax=\sum_j x_j a_j$, target $(Ax)_i=\sum_j A_{ij}x_j$. <br> <i>Hint:</i> take component $i$ of each column, using $(a_j)_i=A_{ij}$. :: A: The $i$-th entry of column $a_j$ is $A_{ij}$ <br> Take component $i$ of both sides: $(Ax)_i=\sum_j x_j (a_j)_i=\sum_j x_j A_{ij}$ <br> $\boxed{(Ax)_i = \sum_{j=1}^{n} A_{ij}\,x_j}$.
- Q: Derive the product entry $C_{ij}$ for $C=AB$. <br> <b>Formula:</b> column $j$ of $C$ is $Ab_j$; entry $i$ is $(Ab_j)_i=\sum_l A_{il}(b_j)_l$. <br> <i>Hint:</i> write $b_j$ as column $j$ of $B$, then use $(b_j)_l=B_{lj}$. :: A: Column $j$ of $C$ is $A b_j$ where $b_j$ is column $j$ of $B$ <br> Its $i$-th entry is $(Ab_j)_i=\sum_l A_{il}(b_j)_l=\sum_l A_{il}B_{lj}$ <br> $\boxed{C_{ij}=\sum_{l=1}^{k} A_{il}B_{lj}}$.
- Q: Derive why composition $C=AB$ reads right-to-left (do $B$ first). <br> <b>Formula:</b> $(AB)x = A(Bx)$ by associativity. <br> <i>Hint:</i> read the right side from the inside out, like the function composition $A(B(x))$. :: A: For any $x$, $(AB)x = A(Bx)$ by associativity <br> So $B$ acts first, then $A$ acts on the result <br> This mirrors function composition $A(B(x))$ <br> $\boxed{(AB)x=A(Bx)\ \Rightarrow\ \text{do }B\text{ first, then }A}$ (hence $AB\neq BA$ in general).
- Q: Derive the transpose-of-a-product rule $(AB)^\top$. <br> <b>Formula:</b> entry definition $(M^\top)_{ij}=M_{ji}$ and product $(AB)_{ji}=\sum_l A_{jl}B_{li}$. <br> <i>Hint:</i> expand $((AB)^\top)_{ij}=(AB)_{ji}$, then rewrite each factor as a transpose entry. :: A: $((AB)^\top)_{ij}=(AB)_{ji}=\sum_l A_{jl}B_{li}$ <br> Rewrite each factor as a transpose entry: $A_{jl}=(A^\top)_{lj}$, $B_{li}=(B^\top)_{il}$ <br> $=\sum_l (B^\top)_{il}(A^\top)_{lj}=(B^\top A^\top)_{ij}$ <br> $\boxed{(AB)^\top = B^\top A^\top}$.
- Q: Derive the inverse-of-a-product rule $(AB)^{-1}$. <br> <b>Formula:</b> the inverse satisfies $(AB)(AB)^{-1}=I$; target $(AB)^{-1}=B^{-1}A^{-1}$. <br> <i>Hint:</i> guess the undo-order $B^{-1}A^{-1}$ and check it cancels from the inside out. :: A: Guess the undo-order: apply $B^{-1}$ then $A^{-1}$ <br> Check: $(AB)(B^{-1}A^{-1})=A(BB^{-1})A^{-1}=A I A^{-1}=I$ <br> Uniqueness of the inverse then forces the guess <br> $\boxed{(AB)^{-1}=B^{-1}A^{-1}}$ (order reverses).
- Q: Derive $\det A$ of $A=\begin{bmatrix}a&b\\c&d\end{bmatrix}$ as the signed area of the column parallelogram (step 1 of 2). <br> <b>Formula:</b> signed area $=$ cross-product $z$-component of columns $(a,c)$ and $(b,d)$. <br> <i>Hint:</i> the 2D cross-product of $(a,c)$ and $(b,d)$ is $a\cdot d - c\cdot b$. :: A: Columns $(a,c)$ and $(b,d)$ span a parallelogram <br> Signed area $=$ cross-product $z$-component $=ad-cb$ <br> $\boxed{\det A = ad-bc}$.
- Q: Derive the $2\times2$ inverse $A^{-1}$ (step 2 of 2). <br> <b>Formula:</b> use $\det A=ad-bc$ and the requirement $AA^{-1}=I$; target $A^{-1}=\frac{1}{ad-bc}\begin{bmatrix}d&-b\\-c&a\end{bmatrix}$. <br> <i>Hint:</i> try the adjugate (swap the diagonal, negate the off-diagonal), multiply by $A$, then divide by $\det A$. :: A: Try $A^{-1}=\frac{1}{\det A}\begin{bmatrix}d&-b\\-c&a\end{bmatrix}$ (swap diagonal, negate off-diagonal) <br> Check: $\begin{bmatrix}a&b\\c&d\end{bmatrix}\begin{bmatrix}d&-b\\-c&a\end{bmatrix}=\begin{bmatrix}ad-bc&0\\0&ad-bc\end{bmatrix}=(\det A) I$ <br> Divide by $\det A$ to get $I$ <br> $\boxed{A^{-1}=\frac{1}{ad-bc}\begin{bmatrix}d&-b\\-c&a\end{bmatrix}}$ (needs $\det A\neq0$).
- Q: Derive what $|\det A|$ and $\operatorname{sign}(\det A)$ mean geometrically, from the unit square. <br> <b>Formula:</b> $A$ maps the unit square (area $1$) to the parallelogram on columns $a_1,a_2$ with signed area $\det A$. <br> <i>Hint:</i> compare output area to input area $1$; a negative sign means the square was flipped over. :: A: $A$ sends the unit square (area $1$) to the parallelogram on columns $a_1,a_2$ <br> Its signed area is $\det A$, so area scales by $|\det A|$ <br> A negative sign means the square was flipped over <br> $\boxed{|\det A|=\text{area/volume scale},\ \ \det A<0\Rightarrow\text{orientation reversed}}$.
- Q: Derive why $\det A=0$ forces $A$ to be singular. <br> <b>Formula:</b> $|\det A|$ is the area/volume scale of the column parallelogram. <br> <i>Hint:</i> set that area to zero and ask what it says about the columns' independence. :: A: $\det A=0$ means the output parallelogram has zero area <br> So the columns are linearly dependent and span $<n$ dimensions <br> The map collapses space onto a lower dimension, so no inverse can recover it <br> $\boxed{\det A=0\iff \text{rank-deficient}\iff \text{singular}}$.
- Q: Compute $\det A$ for the $2\times2$ matrix $A=\begin{bmatrix}2&1\\4&3\end{bmatrix}$ and state the area/orientation effect. <br> <b>Formula:</b> $\det\begin{bmatrix}a&b\\c&d\end{bmatrix}=ad-bc$. <br> <i>Hint:</i> compute the main-diagonal product $2\cdot3$ first, then subtract the off-diagonal product $1\cdot4$. :: A: $\det A=(2)(3)-(1)(4)=6-4=2$ <br> Nonzero, so invertible <br> $\boxed{\det A=2}$: area scaled $\times 2$, orientation preserved.
- Q: Build the $2\times2$ rotation by $30^\circ$ and compute its determinant. <br> <b>Formula:</b> $R=\begin{bmatrix}\cos\theta&-\sin\theta\\ \sin\theta&\cos\theta\end{bmatrix}$, $\det R=\cos^2\theta+\sin^2\theta$. <br> <i>Hint:</i> plug in $\cos30^\circ=\tfrac{\sqrt3}{2},\ \sin30^\circ=\tfrac12$, then add the squares. :: A: $R=\begin{bmatrix}\cos30^\circ&-\sin30^\circ\\ \sin30^\circ&\cos30^\circ\end{bmatrix}=\begin{bmatrix}\tfrac{\sqrt3}{2}&-\tfrac12\\[2pt]\tfrac12&\tfrac{\sqrt3}{2}\end{bmatrix}$ <br> $\det R=\cos^2\theta+\sin^2\theta=\tfrac34+\tfrac14$ <br> $\boxed{\det R=1}$: rotations preserve area and orientation.
- Q: For the reflection $A=\operatorname{diag}(1,-1)$ compute $\det A$ and say what it does to a stone scan. <br> <b>Formula:</b> for a diagonal matrix, $\det\operatorname{diag}(d_1,d_2)=d_1 d_2$. <br> <i>Hint:</i> multiply the two diagonal entries; check the sign for the orientation flip. :: A: $\det A=(1)(-1)=-1$ <br> $|\det A|=1$ so area is preserved <br> Sign $<0$ so winding/handedness flips <br> $\boxed{\det A=-1}$: it mirrors the block (must not be used as a rotation).
- Q: A uniform scale $A=\operatorname{diag}(s,s,s)$ maps a $2\,\text{m}$ stone block; find the volume factor for $s=1.5$. <br> <b>Formula:</b> $\det\operatorname{diag}(s,s,s)=s^3$. <br> <i>Hint:</i> cube $s=1.5$; the block's size in metres does not enter. :: A: $\det A=s\cdot s\cdot s=s^3$ <br> For $s=1.5$: $\det A=1.5^3=3.375$ <br> $\boxed{\det A=3.375}$: volume scales $\times3.375$ regardless of the block's metres.
- Q: Factor $A=\begin{bmatrix}2&1\\4&3\end{bmatrix}$ by LU with partial pivoting (step 1 of 2): pivot and eliminate. <br> <b>Formula:</b> $PA=LU$; multiplier $\ell_{21}=A_{21}/A_{11}$ after pivoting; Row2 $\leftarrow$ Row2 $-\ell_{21}\cdot$Row1. <br> <i>Hint:</i> compare the column-1 magnitudes $2$ and $4$ first, swap so $4$ is the pivot, then compute the multiplier. :: A: Column-1 candidates $2,4$; larger is $4$, so swap rows ($P$ swaps rows 1,2): $PA=\begin{bmatrix}4&3\\2&1\end{bmatrix}$ <br> Multiplier $\ell_{21}=2/4=0.5$; Row2 $\leftarrow$ Row2 $-0.5\cdot$Row1 <br> $\boxed{U=\begin{bmatrix}4&3\\0&-0.5\end{bmatrix},\ L=\begin{bmatrix}1&0\\0.5&1\end{bmatrix}}$.
- Q: Solve $Ax=b$ with $b=(1,7)$ using the $PA=LU$ factors of the previous card (step 2 of 2). <br> <b>Formula:</b> $Ly=Pb$ (forward substitution), then $Ux=y$ (back substitution). <br> <i>Hint:</i> apply $P$ to $b$ first to get $Pb=(7,1)$, then solve $L$ top-down before $U$ bottom-up. :: A: $Pb=(7,1)$. Forward $Ly=Pb$: $y_1=7$, $0.5(7)+y_2=1\Rightarrow y_2=-2.5$ <br> Back $Ux=y$: $-0.5\,x_2=-2.5\Rightarrow x_2=5$; $4x_1+3(5)=7\Rightarrow x_1=-2$ <br> $\boxed{x=(-2,\ 5)}$.
- Q: Closing test: verify the solution $x=(-2,5)$ for $A=\begin{bmatrix}2&1\\4&3\end{bmatrix},\,b=(1,7)$. <br> <b>Formula:</b> substitute $x$ back and check the residual $\|Ax-b\|=0$. <br> <i>Hint:</i> a pass means each component of $Ax$ matches $b$ exactly; compute $Ax$ row by row. :: A: $Ax=\begin{bmatrix}2(-2)+1(5)\\4(-2)+3(5)\end{bmatrix}=\begin{bmatrix}-4+5\\-8+15\end{bmatrix}=\begin{bmatrix}1\\7\end{bmatrix}$ <br> Equals $b$, residual $\|Ax-b\|=0$ <br> $\boxed{Ax=b\ \checkmark}$ solution verified.
- Q: Closing test: re-solve $A=\begin{bmatrix}2&1\\4&3\end{bmatrix},\,b=(1,7)$ via $A^{-1}b$ and confirm it matches the LU answer. <br> <b>Formula:</b> $A^{-1}=\frac{1}{ad-bc}\begin{bmatrix}d&-b\\-c&a\end{bmatrix}$, then $x=A^{-1}b$. <br> <i>Hint:</i> here $\det A=2$, so build $A^{-1}=\tfrac12\begin{bmatrix}3&-1\\-4&2\end{bmatrix}$ first, then multiply by $b$. :: A: $A^{-1}=\tfrac12\begin{bmatrix}3&-1\\-4&2\end{bmatrix}$ <br> $x=A^{-1}b=\tfrac12\begin{bmatrix}3(1)-1(7)\\-4(1)+2(7)\end{bmatrix}=\tfrac12\begin{bmatrix}-4\\10\end{bmatrix}=\begin{bmatrix}-2\\5\end{bmatrix}$ <br> $\boxed{x=(-2,5)}$ matches the LU answer $\checkmark$.
- Q: Show the near-singular matrix $A=\begin{bmatrix}1&1\\1&1+\epsilon\end{bmatrix}$ is ill-conditioned as $\epsilon\to0$. <br> <b>Formula:</b> $\det A=ad-bc$ and $\kappa(A)=\sigma_{\max}/\sigma_{\min}$. <br> <i>Hint:</i> compute $\det A$ first; as it $\to0$ the columns turn parallel so $\sigma_{\min}\to0$ while $\sigma_{\max}\approx2$. :: A: $\det A=(1)(1+\epsilon)-(1)(1)=\epsilon\to0$ <br> Columns become nearly parallel, so $\sigma_{\min}\to0$ while $\sigma_{\max}\approx2$ <br> $\kappa=\sigma_{\max}/\sigma_{\min}\to\infty$ <br> $\boxed{\epsilon\to0\Rightarrow\kappa(A)\to\infty}$: the plane/line fit through near-collinear scan points is ill-posed.
- Q: Closing test: check that $\kappa(A)\ge1$ always, with equality for a pure rotation. <br> <b>Formula:</b> $\kappa(A)=\sigma_{\max}/\sigma_{\min}$, and an orthogonal $R$ has all $\sigma_i=1$. <br> <i>Hint:</i> a pass means the ratio bottoms out at $1$; use $\sigma_{\max}\ge\sigma_{\min}$, then set every $\sigma_i=1$. :: A: Singular values obey $\sigma_{\max}\ge\sigma_{\min}>0$, so the ratio is $\ge1$ <br> A rotation is orthogonal: all $\sigma_i=1$, so $\sigma_{\max}=\sigma_{\min}=1$ <br> $\boxed{\kappa(R)=1/1=1}$: rotations are perfectly conditioned (lower bound hit).
- CLOZE: Each column $j$ of $A$ is the image of the basis vector: $a_j = {{c1::A e_j}}$.
- CLOZE: $Ax$ expands as $\sum_j x_j a_j$, a {{c1::linear combination of the columns of $A$}} weighted by the entries of $x$.
- CLOZE: A square matrix is singular exactly when $\det A = {{c1::0}}$, equivalently when it is rank-deficient.
- CLOZE: To solve a general square system in Eigen, call {{c1::A.partialPivLu().solve(b)}} rather than forming the inverse. <br> <i>hint:</i> LU with partial pivoting, the fast default.
- CLOZE: For a symmetric positive definite system such as the normal matrix $A^\top A$, the fastest Eigen solver is {{c1::llt()}} (Cholesky).
- Q: When-to-use: why factor once with `partialPivLu()` and reuse it for many right-hand sides $b$? <br> <b>Formula:</b> factor cost is $O(n^3)$; each triangular solve is $O(n^2)$ per $b$. <br> <i>Hint:</i> compare the one-time cubic cost against the per-$b$ quadratic cost. :: A: Factoring costs $O(n^3)$ but each back/forward solve is only $O(n^2)$ <br> Reuse amortizes the cubic cost across all the $b$ <br> $\boxed{\text{factor once }O(n^3)\text{, then }O(n^2)\text{ per }b}$.
- Q: Pitfall: a rigid-fit rotation comes out with $\det=-1$. What is wrong and the fix? <br> <b>Formula:</b> a valid rotation in $SO(3)$ has $\det=+1$; after SVD check $\det(VU^\top)$. <br> <i>Hint:</i> $\det=-1$ flags a reflection; flip the last singular vector so the product determinant becomes $+1$. :: A: It is a reflection, not a rotation <br> Flip the last singular vector so $\det(VU^\top)=+1$, else the stone is mirrored <br> $\boxed{\det=-1\Rightarrow\text{reflection; flip last column to force }\det=+1}$.
- Q: Pitfall: why is `det(A) == 0` a bad singularity test in floating point, and what to use instead? <br> <b>Formula:</b> use a scale-relative threshold $|v|<\tau=\varepsilon\cdot\text{scale}$, or a rank-revealing decomposition. <br> <i>Hint:</i> roundoff never lands exactly on $0$; reach for `colPivHouseholderQr().rank()`. :: A: Roundoff never lands exactly on zero, so an exact compare misfires <br> Use a scale-relative threshold $|v|<\tau=\varepsilon\cdot\text{scale}$ or a rank-revealing decomposition (`colPivHouseholderQr().rank()`) <br> $\boxed{\text{test }|v|<\varepsilon\cdot\text{scale, not }==0}$.
