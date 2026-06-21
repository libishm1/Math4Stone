# Vectors: dot, cross, norm, projection (2D & 3D)

Prereqs: none => **this** => unlocks: [[02_matrices]], [[03_transforms-se3]], [[04_discrete-curvature]].

## Why it matters for the stone scan->CAD->IFC / CADFit workflow

A raw stone scan is just a bag of points $\{p_i \in \mathbb{R}^3\}$, and almost every downstream step is vector arithmetic on those points. When you estimate a surface normal at a scan point with PCA, you read it off an eigenvector and orient it with a dot-product sign test; that normal is what feeds Poisson reconstruction and the point-to-plane ICP that registers scans together. When you project a messy point onto a fitting template to build the ordered UV grid that `GeomAPI_PointsToBSplineSurface` needs, you are computing a vector projection. When CADFit scores a candidate CAD program by Chamfer distance, each term is a squared norm $\lVert p_i - q_i\rVert^2$. When you write the B-Rep face normal into IFC, its sign and unit length are exactly the cross-product and normalization rules below. Get these four operations right and the rest of the pipeline has solid ground to stand on; get the frame or the sign wrong and every later transform inherits the bug.

## Intuition

A vector is an arrow: it has a length and a direction, and you are free to slide it around without changing it. Points are places; vectors are displacements between places. Two questions dominate geometry work. First, "how much do these two arrows agree?" That is the dot product: large and positive when they point the same way, zero when perpendicular, negative when opposed. Second, "what is the flat patch these two arrows span, and which way does it face?" That is the cross product: its length is the area of the parallelogram they make, and its direction is the normal sticking out of that patch.

Analogy: think of the dot product as a shadow. Stand a flashlight directly behind vector $b$ and shine it down $b$'s direction. The length of $a$'s shadow on the line of $b$ is the scalar projection. The dot product is that shadow length scaled by how long $b$ is.

## The math (first principles)

**Conventions.** Work in a right-handed Cartesian frame. Coordinates are column vectors. Units are whatever the document declares; in this pipeline that is metres at site/quarry scale and millimetres in the shop. Geometric tolerances are relative: a length is "zero" when it falls below $\varepsilon = \max(\varepsilon_{\text{floor}}, 10^{-3} L)$ for a characteristic length $L$, never a bare hard-coded constant. The cross product is defined in $\mathbb{R}^3$ only.

**Vector, addition, scaling.** A vector in $\mathbb{R}^n$ is $a = (a_1, \dots, a_n)$. Addition is componentwise, $a + b = (a_1+b_1, \dots)$. Scaling by $s \in \mathbb{R}$ is $sa = (sa_1, \dots)$. These two operations define a vector space.

**Norm (length).** The Euclidean norm is

$$\lVert a \rVert = \sqrt{a_1^2 + a_2^2 + \cdots + a_n^2} = \sqrt{a \cdot a}.$$

A **unit vector** has $\lVert a \rVert = 1$. Normalizing is $\hat{a} = a / \lVert a \rVert$, defined only when $\lVert a \rVert > 0$.

**Dot product.** Algebraically $a \cdot b = \sum_i a_i b_i$. Geometrically

$$a \cdot b = \lVert a \rVert\, \lVert b \rVert \cos\theta,$$

where $\theta \in [0, \pi]$ is the angle between $a$ and $b$. The two definitions agree by the law of cosines applied to the triangle with sides $a$, $b$, $a-b$. Consequences: $a \cdot b = 0$ iff the vectors are orthogonal (for nonzero vectors); $a \cdot a = \lVert a \rVert^2$; and the Cauchy-Schwarz bound $\lvert a \cdot b\rvert \le \lVert a\rVert\,\lVert b\rVert$ holds because $\lvert\cos\theta\rvert \le 1$. To recover the angle safely,

$$\theta = \operatorname{atan2}\big(\lVert a \times b \rVert,\; a \cdot b\big),$$

which stays well conditioned near $0$ and $\pi$ where $\arccos$ of a near-$\pm 1$ argument loses precision.

**Scalar and vector projection.** The scalar projection (signed shadow length) of $a$ onto $b$ is

$$\operatorname{comp}_b a = \frac{a \cdot b}{\lVert b \rVert}.$$

The vector projection is that shadow placed back on $b$'s direction:

$$\operatorname{proj}_b a = \frac{a \cdot b}{\lVert b \rVert^2}\, b = (a \cdot \hat b)\,\hat b.$$

*Derivation.* The projection is parallel to $b$, so write it as $C\,b$. The residual $a - C b$ must be orthogonal to $b$, so $(a - Cb)\cdot b = 0$, giving $a\cdot b = C\,(b\cdot b)$ and $C = (a\cdot b)/(b\cdot b)$. The orthogonal remainder

$$\operatorname{perp}_b a = a - \operatorname{proj}_b a$$

is the component of $a$ orthogonal to $b$. This split (parallel plus orthogonal) is the engine behind Gram-Schmidt and behind point-to-plane residuals.

**Cross product (3D).** For $a, b \in \mathbb{R}^3$,

$$a \times b = \begin{pmatrix} a_2 b_3 - a_3 b_2 \\ a_3 b_1 - a_1 b_3 \\ a_1 b_2 - a_2 b_1 \end{pmatrix}.$$

Its magnitude is

$$\lVert a \times b \rVert = \lVert a \rVert\,\lVert b \rVert \sin\theta,$$

which equals the **area of the parallelogram** spanned by $a$ and $b$ (so the triangle area is half that). The direction is orthogonal to both $a$ and $b$, fixed by the right-hand rule. Key identities: $a \times b = -(b \times a)$ (anticommutative), $a \times a = 0$, and the Lagrange identity $\lVert a \times b\rVert^2 = \lVert a\rVert^2\lVert b\rVert^2 - (a\cdot b)^2$, which is just $\sin^2 + \cos^2 = 1$ in disguise. A **surface normal** from three non-collinear points $p_0, p_1, p_2$ is $n = (p_1 - p_0) \times (p_2 - p_0)$, then $\hat n = n/\lVert n\rVert$. The winding order of the points sets the sign of $n$; this is the IFC face-orientation convention you must keep consistent.

**Linear independence, basis, orthonormal frame.** Vectors are linearly independent if no nonzero combination sums to zero. A **basis** of $\mathbb{R}^n$ is $n$ independent vectors; every vector has unique coordinates in it. An **orthonormal basis** $\{e_1, e_2, e_3\}$ satisfies $e_i \cdot e_j = \delta_{ij}$. In such a frame, coordinates are dot products: $a = \sum_i (a \cdot e_i)\, e_i$. **Gram-Schmidt** builds an orthonormal frame from independent vectors by subtracting projections:

$$u_1 = v_1,\quad u_k = v_k - \sum_{j<k} \operatorname{proj}_{u_j} v_k,\quad e_k = \frac{u_k}{\lVert u_k\rVert}.$$

In 3D the fast route is: $e_1 = \hat v_1$, $e_3 = \widehat{e_1 \times v_2}$, $e_2 = e_3 \times e_1$. This is exactly how you turn a noisy scan normal plus an arbitrary tangent hint into a clean local frame for UV resampling. A right-handed orthonormal frame stacked as columns is a rotation matrix $R \in SO(3)$, which is the bridge to [[03_transforms-se3]].

## Worked example

Let $a = (3, 4, 0)$ and $b = (1, 0, 0)$, all in metres.

**Norms.** $\lVert a\rVert = \sqrt{9 + 16 + 0} = 5$. $\lVert b\rVert = 1$.

**Dot product.** $a \cdot b = 3\cdot 1 + 4\cdot 0 + 0 = 3$.

**Angle.** $\cos\theta = (a\cdot b)/(\lVert a\rVert\lVert b\rVert) = 3/5 = 0.6$, so $\theta = \arccos 0.6 \approx 53.13^\circ$.

**Scalar projection** of $a$ onto $b$: $(a\cdot b)/\lVert b\rVert = 3/1 = 3$.

**Vector projection**: $\operatorname{proj}_b a = (3/1^2)(1,0,0) = (3, 0, 0)$. **Orthogonal part**: $a - (3,0,0) = (0, 4, 0)$. Check: $(0,4,0)\cdot(1,0,0) = 0$. Correct, the residual is perpendicular to $b$.

**Cross product.** $a \times b = (4\cdot 0 - 0\cdot 0,\; 0\cdot 1 - 3\cdot 0,\; 3\cdot 0 - 4\cdot 1) = (0, 0, -4)$. Magnitude $4$. Cross-check via $\lVert a\rVert\lVert b\rVert\sin\theta = 5 \cdot 1 \cdot \sin(53.13^\circ) = 5 \cdot 0.8 = 4$. Correct. The parallelogram on $a$ and $b$ has area $4\,\text{m}^2$; the triangle has area $2\,\text{m}^2$. The normal points in $-z$, consistent with the right-hand rule for $a$ swept toward $b$.

**Lagrange check.** $\lVert a\rVert^2\lVert b\rVert^2 - (a\cdot b)^2 = 25\cdot 1 - 9 = 16 = \lVert a\times b\rVert^2$. Correct.

## Code

Idiomatic C++ with Eigen. Eigen gives `dot`, `cross`, `norm`, `squaredNorm`, `normalized`, and `stableNorm` directly, so the math maps one-to-one onto method calls.

```cpp
// vectors.cpp  --  build: g++ -I/path/to/eigen -O2 vectors.cpp -o vectors
#include <Eigen/Dense>
#include <iostream>
#include <cmath>

using Eigen::Vector3d;

// scalar projection comp_b(a) = (a . b) / ||b||
double scalarProjection(const Vector3d& a, const Vector3d& b) {
    return a.dot(b) / b.norm();                  // a . b  over  ||b||
}

// vector projection proj_b(a) = (a . b / ||b||^2) b
Vector3d vectorProjection(const Vector3d& a, const Vector3d& b) {
    return (a.dot(b) / b.squaredNorm()) * b;     // squaredNorm == b . b
}

// robust angle theta = atan2(||a x b||, a . b), stable near 0 and pi
double angleBetween(const Vector3d& a, const Vector3d& b) {
    return std::atan2(a.cross(b).norm(), a.dot(b));
}

// unit normal of triangle (p0,p1,p2); winding order sets the sign
Vector3d triangleNormal(const Vector3d& p0, const Vector3d& p1, const Vector3d& p2) {
    Vector3d n = (p1 - p0).cross(p2 - p0);       // a x b  =>  normal
    double len = n.norm();
    if (len < 1e-12) return Vector3d::Zero();    // degenerate: collinear points
    return n / len;                              // n_hat = n / ||n||
}

// right-handed orthonormal frame from a normal and a tangent hint (Gram-Schmidt, 3D)
void orthonormalFrame(const Vector3d& normal, const Vector3d& tangentHint,
                      Vector3d& e1, Vector3d& e2, Vector3d& e3) {
    e3 = normal.normalized();                    // e3 = n_hat
    e1 = (tangentHint - e3 * tangentHint.dot(e3)).normalized();  // remove component along e3
    e2 = e3.cross(e1);                           // completes the right-handed frame
}

int main() {
    Vector3d a(3, 4, 0), b(1, 0, 0);
    std::cout << "norm(a)        = " << a.norm() << "\n";          // 5
    std::cout << "a . b          = " << a.dot(b) << "\n";          // 3
    std::cout << "angle (deg)    = " << angleBetween(a, b) * 180.0 / M_PI << "\n"; // 53.13
    std::cout << "scalar proj    = " << scalarProjection(a, b) << "\n";            // 3
    std::cout << "vector proj    = " << vectorProjection(a, b).transpose() << "\n";// 3 0 0
    std::cout << "a x b          = " << a.cross(b).transpose() << "\n";            // 0 0 -4
    std::cout << "parallelogram  = " << a.cross(b).norm() << "\n";                 // 4
    return 0;
}
```

Notes mapping math to code. `squaredNorm()` is preferred over `norm()*norm()` inside comparisons and projections: it avoids a square root and the associated rounding. For inputs with extreme magnitudes use `stableNorm()`/`stableNormalize()`, which rescale to avoid overflow/underflow. Eigen's `cross()` is defined only for size-3 vectors at compile time; calling it on a `Vector2d` is a build error, which is the type system protecting you.

A short Python mirror for prototyping (NumPy):

```python
import numpy as np
a = np.array([3.0, 4.0, 0.0]); b = np.array([1.0, 0.0, 0.0])
proj = (a @ b) / (b @ b) * b              # vector projection
angle = np.arctan2(np.linalg.norm(np.cross(a, b)), a @ b)  # robust angle
n = np.cross(b, a); n = n / np.linalg.norm(n)              # unit normal
```

## Connections

- [[02_matrices]]: stacking vectors as rows or columns gives matrices; the dot product is a $1\times n$ by $n\times 1$ matrix product, and a matrix-vector product is a vector of dot products.
- [[03_transforms-se3]]: an orthonormal frame built here by Gram-Schmidt is exactly the rotation block $R \in SO(3)$ inside a rigid transform $T \in SE(3)$; translation is vector addition.
- [[04_discrete-curvature]]: vertex normals, angle defect, and the cotangent weights in the Laplacian are all built from dot and cross products of edge vectors.
- ICP and registration (later topic): the rigid-fit error is a sum of squared norms, and point-to-plane residuals are dot products of displacement with the surface normal.
- NURBS / B-spline fitting (later topic): control-point arithmetic, tangents, and surface normals are vector operations; ordered UV resampling uses projection.

## Pitfalls and failure modes

- **Frame inconsistency.** A vector is meaningless without its frame. Mixing a normal expressed in the scanner frame with points in the world frame silently corrupts every dot product downstream. Label frames explicitly; most "math bugs" in this pipeline are frame bugs.
- **Handedness / winding flips the normal.** $a \times b = -(b \times a)$. A reversed triangle winding flips the IFC face normal and can invert inside/outside tests for booleans. Pick a winding convention once and enforce it.
- **`arccos` near $0$ and $\pi$.** Computing angle via `acos((a·b)/(|a||b|))` loses precision and can return NaN when the argument drifts just past $\pm 1$ due to rounding (a documented Eigen/PCL gotcha). Use `atan2(||a×b||, a·b)` instead.
- **Normalizing a near-zero vector.** Dividing by $\lVert a\rVert$ when the vector is tiny amplifies noise into a garbage direction. Guard with a relative tolerance and treat below-threshold lengths as degenerate.
- **Degenerate inputs.** Collinear points give a zero cross product and no normal; coincident points give a zero displacement. Detect these before they propagate.
- **Hard-coded absolute tolerances.** A $10^{-6}$ floor is wrong at quarry scale (metres) and wrong at shop scale (millimetres). Scale tolerances by a characteristic length.
- **Projection onto a non-unit vector.** Forgetting the $1/\lVert b\rVert^2$ factor (or using $\hat b$ where you needed $b$) is the single most common projection bug. The form $(a\cdot\hat b)\,\hat b$ avoids it.

## Practice

1. By hand, take $a=(1,2,2)$ and $b=(0,0,3)$. Compute $\lVert a\rVert$, $a\cdot b$, the angle in degrees, $\operatorname{proj}_b a$, and the orthogonal remainder. Verify the remainder is perpendicular to $b$.
2. By hand, compute $a\times b$ for the same vectors, its magnitude, and confirm the Lagrange identity numerically.
3. Write a C++ function `bool isRightHanded(e1,e2,e3)` that returns true when $(e_1\times e_2)\cdot e_3 > 0$. Feed it the output of `orthonormalFrame` and print the result.
4. Given three scan points of your choice, compute the unit triangle normal, then reverse the winding order and confirm the normal flips sign. Print both.
5. Implement scalar projection two ways: with `norm()` and with `stableNorm()`. Feed a vector scaled by $10^{200}$ and compare. Observe where the naive version overflows.
6. Take a small synthetic point cloud (say 8 points on a plane), fit a normal by averaging the cross products of consecutive edge vectors about the centroid, normalize it, and orient it so it points toward a reference viewpoint using a dot-product sign test.

## References

- Lynch, K. and Park, F. *Modern Robotics: Mechanics, Planning, and Control*. Cambridge University Press. Foundations of $SE(3)$, frames, and rigid motion. http://hades.mech.northwestern.edu/index.php/Modern_Robotics
- OpenCV. *Camera Calibration and 3D Reconstruction (calib3d) documentation*. Pinhole projection and homogeneous coordinates. https://docs.opencv.org/4.x/d9/d0c/group__calib3d.html
- *Vector Projection - Formula, Definition, Derivation* (Cuemath); *Scalar projection* (Wikipedia). Projection formulas and derivation. https://en.wikipedia.org/wiki/Scalar_projection
- *Cross product* (Wikipedia). Magnitude as parallelogram area, right-hand rule, Lagrange identity. https://en.wikipedia.org/wiki/Cross_product
- Eigen documentation. `dot`, `cross`, `norm`, `squaredNorm`, `normalized`, `stableNorm`. https://eigen.tuxfamily.org/dox/group__TutorialMatrixArithmetic.html
- *Gram-Schmidt process* (Purdue MA262 notes; Wikipedia). Orthonormal basis construction by projection subtraction. https://en.wikipedia.org/wiki/Gram%E2%80%93Schmidt_process
- Crane, K. *Discrete Differential Geometry: An Applied Introduction* (CMU course notes). Vertex normals and angle-based geometry on meshes. https://www.cs.cmu.edu/~kmcrane/Projects/DDG/
- OpenCASCADE Technology reference manual. `GeomAPI_PointsToBSplineSurface` ordered-point fitting consumes vector data. https://dev.opencascade.org/doc/overview/html/
- Nehme, G., Whalen, E., Ahmed, F. *CADFit: Precise Mesh-to-CAD Program Generation with Hybrid Optimization*. arXiv:2605.01171. IoU/Chamfer objectives built on norms. https://arxiv.org/abs/2605.01171 ; code https://github.com/ghadinehme/CADFit

## Flashcard seeds

- Q: Derive the geometric dot product (step 1 of 2): get an expression for $\lVert a-b\rVert^2$ by expanding it algebraically. <br> <b>Formula:</b> $\lVert a-b\rVert^2 = (a-b)\cdot(a-b)$ and the law of cosines $\lVert a-b\rVert^2 = \lVert a\rVert^2 + \lVert b\rVert^2 - 2\lVert a\rVert\lVert b\rVert\cos\theta$ <br> <i>Hint:</i> Expand $(a-b)\cdot(a-b)$ using distributivity, then set it equal to the law-of-cosines form. :: A: Algebraic expansion: $\lVert a-b\rVert^2 = (a-b)\cdot(a-b) = \lVert a\rVert^2 - 2(a\cdot b) + \lVert b\rVert^2$ <br> Law of cosines on the triangle with sides $a$, $b$, $a-b$: $\lVert a-b\rVert^2 = \lVert a\rVert^2 + \lVert b\rVert^2 - 2\lVert a\rVert\lVert b\rVert\cos\theta$ <br> Set the two equal <br> $\boxed{\lVert a\rVert^2 - 2(a\cdot b) + \lVert b\rVert^2 = \lVert a\rVert^2 + \lVert b\rVert^2 - 2\lVert a\rVert\lVert b\rVert\cos\theta}$
- Q: Derive the geometric dot product (step 2 of 2): solve the matched expression for $a\cdot b$. <br> <b>Formula:</b> $\lVert a\rVert^2 - 2(a\cdot b) + \lVert b\rVert^2 = \lVert a\rVert^2 + \lVert b\rVert^2 - 2\lVert a\rVert\lVert b\rVert\cos\theta$; target $a\cdot b = \lVert a\rVert\,\lVert b\rVert\cos\theta$ <br> <i>Hint:</i> Cancel the $\lVert a\rVert^2+\lVert b\rVert^2$ that appears on both sides, then divide by $-2$. :: A: Cancel $\lVert a\rVert^2 + \lVert b\rVert^2$ from both sides: $-2(a\cdot b) = -2\lVert a\rVert\lVert b\rVert\cos\theta$ <br> Divide by $-2$ <br> $\boxed{a\cdot b = \lVert a\rVert\,\lVert b\rVert\cos\theta}$
- Q: Derive the robust angle formula and say why it beats $\arccos$. <br> <b>Formula:</b> $a\cdot b=\lVert a\rVert\lVert b\rVert\cos\theta$, $\lVert a\times b\rVert=\lVert a\rVert\lVert b\rVert\sin\theta$; target $\theta = \operatorname{atan2}(\lVert a\times b\rVert,\, a\cdot b)$ <br> <i>Hint:</i> Divide the cross magnitude by the dot product so the $\lVert a\rVert\lVert b\rVert$ factors cancel into $\tan\theta$. :: A: Divide the two: $\dfrac{\lVert a\times b\rVert}{a\cdot b} = \dfrac{\sin\theta}{\cos\theta} = \tan\theta$ <br> Invert with the two-argument arctangent so the quadrant and sign are kept <br> $\boxed{\theta = \operatorname{atan2}\big(\lVert a\times b\rVert,\; a\cdot b\big)}$ <br> Stable near $0$ and $\pi$, where $\arccos$ of a near-$\pm1$ argument loses precision or returns NaN.
- Q: Derive the vector projection (step 1 of 2): the projection of $a$ onto $b$ is parallel to $b$, so write it as $C\,b$ and solve for the scalar $C$. <br> <b>Formula:</b> residual orthogonality $(a - Cb)\cdot b = 0$; target $C = \dfrac{a\cdot b}{\lVert b\rVert^2}$ <br> <i>Hint:</i> Expand the dot product $(a-Cb)\cdot b=0$ and isolate $C$, using $b\cdot b=\lVert b\rVert^2$. :: A: The residual $a - Cb$ must be orthogonal to $b$, so $(a - Cb)\cdot b = 0$ <br> Expand: $a\cdot b - C\,(b\cdot b) = 0$ <br> $\boxed{C = \dfrac{a\cdot b}{b\cdot b} = \dfrac{a\cdot b}{\lVert b\rVert^2}}$
- Q: Derive the vector projection (step 2 of 2): write the final projection and its unit-vector form. <br> <b>Formula:</b> $C=\dfrac{a\cdot b}{\lVert b\rVert^2}$, $\operatorname{proj}_b a = C\,b$; target $\operatorname{proj}_b a = \dfrac{a\cdot b}{\lVert b\rVert^2}\,b = (a\cdot\hat b)\,\hat b$ <br> <i>Hint:</i> Substitute $C$ into $C\,b$, then split one $\lVert b\rVert$ off each side to form $\hat b = b/\lVert b\rVert$. :: A: Substitute $C$ back into $C\,b$: $\operatorname{proj}_b a = \dfrac{a\cdot b}{\lVert b\rVert^2}\,b$ <br> Group $b/\lVert b\rVert = \hat b$ to get the safe form $\operatorname{proj}_b a = \left(a\cdot\dfrac{b}{\lVert b\rVert}\right)\dfrac{b}{\lVert b\rVert}$ <br> $\boxed{\operatorname{proj}_b a = \dfrac{a\cdot b}{\lVert b\rVert^2}\,b = (a\cdot\hat b)\,\hat b}$
- Q: Derive the scalar projection (signed shadow length) of $a$ onto $b$. <br> <b>Formula:</b> start from $\operatorname{proj}_b a = \dfrac{a\cdot b}{\lVert b\rVert^2}\,b$; target $\operatorname{comp}_b a = \dfrac{a\cdot b}{\lVert b\rVert}$ <br> <i>Hint:</i> Take the signed length by multiplying the scalar coefficient by $\lVert b\rVert$, then cancel one factor of $\lVert b\rVert$. :: A: Take the signed length: $\operatorname{comp}_b a = \dfrac{a\cdot b}{\lVert b\rVert^2}\,\lVert b\rVert$ <br> Cancel one factor of $\lVert b\rVert$ <br> $\boxed{\operatorname{comp}_b a = \dfrac{a\cdot b}{\lVert b\rVert}}$
- Q: Derive the cross-product magnitude (step 1 of 2): express $\lVert a\times b\rVert$ via $\theta$. <br> <b>Formula:</b> Lagrange identity $\lVert a\times b\rVert^2 = \lVert a\rVert^2\lVert b\rVert^2 - (a\cdot b)^2$ and $a\cdot b = \lVert a\rVert\lVert b\rVert\cos\theta$; target $\lVert a\times b\rVert = \lVert a\rVert\lVert b\rVert\sin\theta$ <br> <i>Hint:</i> Substitute the geometric dot product into the Lagrange identity, then use $1-\cos^2\theta=\sin^2\theta$. :: A: Substitute $a\cdot b = \lVert a\rVert\lVert b\rVert\cos\theta$: $\lVert a\times b\rVert^2 = \lVert a\rVert^2\lVert b\rVert^2 - \lVert a\rVert^2\lVert b\rVert^2\cos^2\theta$ <br> Factor and use $1-\cos^2\theta=\sin^2\theta$: $= \lVert a\rVert^2\lVert b\rVert^2\sin^2\theta$ <br> Take the root ($\sin\theta\ge0$ for $\theta\in[0,\pi]$) <br> $\boxed{\lVert a\times b\rVert = \lVert a\rVert\lVert b\rVert\sin\theta}$
- Q: Derive the parallelogram area (step 2 of 2): show $\lVert a\times b\rVert$ equals the area of the parallelogram on $a$ and $b$, and give the triangle area. <br> <b>Formula:</b> $\lVert a\times b\rVert = \lVert a\rVert\lVert b\rVert\sin\theta$; area $=$ base $\times$ height; target $A_{\triangle} = \tfrac12\lVert a\times b\rVert$ <br> <i>Hint:</i> Take base $=\lVert a\rVert$ and height $=\lVert b\rVert\sin\theta$ (the part of $b$ perpendicular to $a$); the triangle is half. :: A: The parallelogram has base $\lVert a\rVert$ and height $\lVert b\rVert\sin\theta$ (the component of $b$ perpendicular to $a$) <br> Area $=$ base $\times$ height $= \lVert a\rVert\,(\lVert b\rVert\sin\theta) = \lVert a\times b\rVert$ <br> The triangle is half the parallelogram <br> $\boxed{A_{\triangle} = \tfrac12\lVert a\times b\rVert}$
- Q: Derive a unit surface normal from three points (step 1 of 2): given non-collinear $p_0,p_1,p_2$, build a vector orthogonal to both edges and make it unit length. <br> <b>Formula:</b> $\hat n = \dfrac{(p_1-p_0)\times(p_2-p_0)}{\lVert(p_1-p_0)\times(p_2-p_0)\rVert}$ <br> <i>Hint:</i> Form the two edge vectors from $p_0$ first, cross them, then divide by the magnitude. :: A: Edge vectors from $p_0$: $e_1 = p_1-p_0$, $e_2 = p_2-p_0$ <br> Their cross product is orthogonal to both: $n = e_1\times e_2$ <br> Normalize <br> $\boxed{\hat n = \dfrac{(p_1-p_0)\times(p_2-p_0)}{\lVert(p_1-p_0)\times(p_2-p_0)\rVert}}$
- Q: Orient a normal toward a viewpoint (step 2 of 2): derive the sign test that flips a raw unit normal $\hat n$ to face viewpoint $v$ (the IFC/scan orientation step). <br> <b>Formula:</b> view direction $d = v - p_0$; flip rule $\hat n \leftarrow \operatorname{sign}(\hat n\cdot d)\,\hat n$ <br> <i>Hint:</i> Compute the dot product of $\hat n$ with the direction to the viewer; positive means it already faces the viewer. :: A: Form the view direction from the surface point $p_0$: $d = v - p_0$ <br> $\hat n$ faces the viewer iff their dot product is positive: $\hat n\cdot d > 0$ <br> Flip when it is not <br> $\boxed{\hat n \leftarrow \operatorname{sign}(\hat n\cdot d)\,\hat n}$
- Q: Derive how a single coordinate is recovered in an orthonormal basis. <br> <b>Formula:</b> $a=\sum_j c_j e_j$ with $e_i\cdot e_j=\delta_{ij}$; target $c_i = a\cdot e_i$ <br> <i>Hint:</i> Dot both sides with $e_i$; the Kronecker delta kills every term except $j=i$. :: A: Dot both sides with $e_i$: $a\cdot e_i = \sum_j c_j (e_j\cdot e_i)$ <br> Only the $j=i$ term survives because $e_j\cdot e_i=\delta_{ij}$ <br> $\boxed{c_i = a\cdot e_i,\quad\text{so } a=\sum_i (a\cdot e_i)\,e_i}$
- Q: (worked example, mm scale) In the shop, $a=(6,8,0)\,\text{mm}$ and $b=(2,0,0)\,\text{mm}$. Compute $\lVert a\rVert$ and $a\cdot b$. <br> <b>Formula:</b> $\lVert a\rVert=\sqrt{a_1^2+a_2^2+a_3^2}$, $a\cdot b=\sum_i a_ib_i$ <br> <i>Hint:</i> Square and sum the components of $a$ for the norm; multiply matched components and add for the dot product. :: A: $\lVert a\rVert=\sqrt{6^2+8^2+0}=\sqrt{36+64}=\sqrt{100}=10\,\text{mm}$ <br> $a\cdot b = 6\cdot2+8\cdot0+0=12\,\text{mm}^2$ <br> $\boxed{\lVert a\rVert=10\,\text{mm},\ a\cdot b=12\,\text{mm}^2}$
- Q: (worked example, mm scale) For $a=(6,8,0)\,\text{mm}$, $b=(2,0,0)\,\text{mm}$, compute $\operatorname{proj}_b a$ and the orthogonal remainder $\operatorname{perp}_b a$. <br> <b>Formula:</b> $\operatorname{proj}_b a = \dfrac{a\cdot b}{\lVert b\rVert^2}\,b$, $\operatorname{perp}_b a = a - \operatorname{proj}_b a$ <br> <i>Hint:</i> Compute $\lVert b\rVert^2$ and $a\cdot b$ first, get the scalar $C=(a\cdot b)/\lVert b\rVert^2$, then scale $b$ and subtract. :: A: $\lVert b\rVert^2 = 4$, $a\cdot b=12$, so $C=12/4=3$ <br> $\operatorname{proj}_b a = 3\,(2,0,0)=(6,0,0)\,\text{mm}$ <br> $\operatorname{perp}_b a = (6,8,0)-(6,0,0)=(0,8,0)\,\text{mm}$ <br> $\boxed{\operatorname{proj}_b a=(6,0,0),\ \operatorname{perp}_b a=(0,8,0)}$
- Q: (worked example, site/quarry scale, m) A scan displacement $a=(300,400,0)\,\text{m}$ projected onto the easting axis $b=(1000,0,0)\,\text{m}$. Compute $\operatorname{proj}_b a$. <br> <b>Formula:</b> $\operatorname{proj}_b a = \dfrac{a\cdot b}{\lVert b\rVert^2}\,b$ <br> <i>Hint:</i> Compute $a\cdot b$ and $\lVert b\rVert^2$ first; notice the answer should pick off the $x$-component just like the mm case. :: A: $a\cdot b = 300\cdot1000 = 3\times10^5$, $\lVert b\rVert^2 = 10^6$ <br> $\operatorname{proj}_b a = \dfrac{3\times10^5}{10^6}(1000,0,0) = (300,0,0)\,\text{m}$ <br> Same shape as the mm case (projection picks off the $x$-component) confirming the formula is scale-free <br> $\boxed{\operatorname{proj}_b a=(300,0,0)\,\text{m}}$
- Q: (worked example, 3D normal) Three scan points $p_0=(0,0,0)$, $p_1=(2,0,0)$, $p_2=(0,3,0)$ (m). Compute the unit normal $\hat n$. <br> <b>Formula:</b> $n=(p_1-p_0)\times(p_2-p_0)$, $\hat n = n/\lVert n\rVert$, with $a\times b=(a_2b_3-a_3b_2,\ a_3b_1-a_1b_3,\ a_1b_2-a_2b_1)$ <br> <i>Hint:</i> Compute the two edge vectors first, then plug into the cross-product component formula. :: A: $e_1=p_1-p_0=(2,0,0)$, $e_2=p_2-p_0=(0,3,0)$ <br> $n=e_1\times e_2=(0\cdot0-0\cdot3,\ 0\cdot0-2\cdot0,\ 2\cdot3-0\cdot0)=(0,0,6)$ <br> $\lVert n\rVert=6$, so $\hat n=(0,0,1)$ <br> $\boxed{\hat n=(0,0,1)}$
- Q: (closing test, close the loop) For $a=(6,8,0)$, $b=(2,0,0)$ you found $\operatorname{perp}_b a=(0,8,0)$. Verify it is orthogonal to $b$ and that $\operatorname{proj}+\operatorname{perp}$ rebuilds $a$. <br> <b>Formula:</b> orthogonality test $\operatorname{perp}_b a\cdot b = 0$; reconstruction $\operatorname{proj}_b a + \operatorname{perp}_b a = a$ <br> <i>Hint:</i> A pass means the dot product is exactly $0$ and the two parts sum back to $a$ componentwise. :: A: Orthogonality: $(0,8,0)\cdot(2,0,0)=0+0+0=0$, so $\operatorname{perp}_b a\perp b$ <br> Reconstruction: $(6,0,0)+(0,8,0)=(6,8,0)=a$ <br> $\boxed{\operatorname{perp}_b a\cdot b=0\ \text{and}\ \operatorname{proj}+\operatorname{perp}=a\ \checkmark}$
- Q: (closing test, close the loop) For $a=(1,2,2)$, $b=(0,0,3)$, verify the Lagrange identity numerically. <br> <b>Formula:</b> $\lVert a\times b\rVert^2=\lVert a\rVert^2\lVert b\rVert^2-(a\cdot b)^2$, with $a\times b=(a_2b_3-a_3b_2,\ a_3b_1-a_1b_3,\ a_1b_2-a_2b_1)$ <br> <i>Hint:</i> Compute $a\times b$ and its squared norm (left side) and the right side separately; a pass means both equal the same number. :: A: $a\times b=(2\cdot3-2\cdot0,\ 2\cdot0-1\cdot3,\ 1\cdot0-2\cdot0)=(6,-3,0)$, so $\lVert a\times b\rVert^2=36+9+0=45$ <br> RHS: $\lVert a\rVert^2\lVert b\rVert^2-(a\cdot b)^2 = 9\cdot9 - 6^2 = 81-36=45$ <br> Both sides equal $45$ <br> $\boxed{45=45\ \checkmark}$
- Q: (closing test, limiting case) Check the angle of perpendicular vectors $a=(1,0,0)$, $b=(0,1,0)$ two ways. <br> <b>Formula:</b> $\theta=\operatorname{atan2}(\lVert a\times b\rVert,\,a\cdot b)$; cross-check via $a\cdot b=0\Rightarrow\theta=90^\circ$ <br> <i>Hint:</i> Compute $a\times b$ and $a\cdot b$ first; a pass means both methods give $90^\circ$. :: A: $a\times b=(0,0,1)$ so $\lVert a\times b\rVert=1$; $a\cdot b=0$ <br> $\theta=\operatorname{atan2}(1,0)=\pi/2=90^\circ$ <br> Cross-check: $a\cdot b=0\Rightarrow\cos\theta=0\Rightarrow\theta=90^\circ$ <br> $\boxed{\theta=90^\circ\ \checkmark}$
- Q: Definition: when is $\hat a=a/\lVert a\rVert$ well defined, and what is the pitfall otherwise? <br> <b>Formula:</b> $\hat a = a/\lVert a\rVert$ with tolerance $\varepsilon=\max(\varepsilon_{\text{floor}},10^{-3}L)$ <br> <i>Hint:</i> The danger is dividing by a length that is below $\varepsilon$. :: A: Defined only when $\lVert a\rVert>0$; guard with a relative tolerance $\varepsilon=\max(\varepsilon_{\text{floor}},10^{-3}L)$ <br> $\boxed{\text{Normalizing a near-zero vector amplifies noise into a garbage direction}}$
- Q: When (and why) do you prefer `squaredNorm()` over `norm()` in projection and comparison code? <br> <b>Formula:</b> $\lVert b\rVert^2 = b\cdot b$; projection scalar $C=(a\cdot b)/(b\cdot b)$ <br> <i>Hint:</i> Think about what computing a square root costs and whether the projection even needs it. :: A: Inside projections and length comparisons it avoids a square root and its rounding error; $\lVert b\rVert^2=b\cdot b$ is exactly what the projection scalar $C=(a\cdot b)/(b\cdot b)$ needs.
- Q: Pitfall: what is the single most common projection bug, and which form avoids it? <br> <b>Formula:</b> $\operatorname{proj}_b a = \dfrac{a\cdot b}{\lVert b\rVert^2}\,b = (a\cdot\hat b)\,\hat b$ <br> <i>Hint:</i> Compare the $1/\lVert b\rVert^2$ factor against using $\hat b$, and which form makes the unit vectors explicit. :: A: Forgetting the $1/\lVert b\rVert^2$ factor (or using $\hat b$ where you needed $b$). The form $(a\cdot\hat b)\,\hat b$ avoids it because both unit factors are explicit.
- Q: When should you prefer `stableNorm()` over `norm()` in Eigen? <br> <i>Hint:</i> Think about inputs scaled by something like $10^{200}$ or $10^{-200}$. :: A: When vector components may be extremely large or small, to avoid overflow/underflow.
- CLOZE: The cross product is {{c1::anticommutative}}, so $a\times b = -(b\times a)$.
- CLOZE: In an orthonormal basis the coordinates of $a$ are recovered as {{c1::dot products}}: $a=\sum_i (a\cdot e_i)e_i$.
- CLOZE: The cross product is defined only in {{c1::three}} dimensions.
- CLOZE: Reversing a triangle's {{c1::winding order}} flips the sign of its computed normal, which can invert IFC face orientation.
- CLOZE: A length counts as "zero" when it falls below the relative tolerance {{c1::$\varepsilon=\max(\varepsilon_{\text{floor}},10^{-3}L)$}}, never a bare hard-coded constant.
