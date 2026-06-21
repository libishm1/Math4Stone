# Homogeneous coordinates and rigid transforms SE(3)

Prereqs: [[03_matrices]], [[04_trig-rotations]] => this => [[06_camera-pinhole]], [[07_forward-kinematics]], [[18_brep-boolean-occt]]

## Why it matters for the stone scan->CAD->IFC / CADFit workflow

Every stage of the stone pipeline lives in a different coordinate frame, and the rigid transform is what bolts them together. A laser scanner returns points in the scanner frame; you register them into a site frame, then into a block frame whose axes line up with the saw bed, then into the IFC project frame for export. Each of those moves is one $4\times4$ matrix $T \in SE(3)$, and getting one of them backwards or in the wrong units is the single most common bug in a scan-to-CAD system. When you fit a B-spline patch in OpenCASCADE and place it as an `IfcLocalPlacement`, that placement is literally an $SE(3)$ element. When CADFit synthesizes a CAD program that re-creates the block, every `Translate`/`Rotate` operation it emits composes into the same chain of rigid transforms. Master this object once and the camera model, the robot arm, and the CAD kernel all become the same arithmetic wearing different clothes.

## Intuition

A frame is a little set of three perpendicular arrows (x, y, z) glued to an origin, sitting somewhere in space. A rigid transform is the instruction "pick up that frame, carry it to a new spot without bending or stretching it, and set it down rotated." Nothing about the object changes shape; only where it sits and which way it points.

The trick of homogeneous coordinates is to glue a fourth number, $w$, onto every 3D vector so that "rotate" and "move" both become a single matrix multiply. Rotation alone is linear and a $3\times3$ matrix handles it. Translation is not linear (it does not fix the origin), so a $3\times3$ matrix cannot do it. By lifting to 4D we smuggle translation into the last column of a $4\times4$ matrix, and now one multiplication does both.

Analogy: think of a point as a person standing on a turntable that is itself on a moving cart. The turntable spin is the rotation $R$; the cart's displacement is the translation $t$. The $w=1$ flag says "you are a location, you ride the cart." A direction (like which way the person is facing) carries a $w=0$ flag: it spins with the turntable but ignores the cart, because a heading has no position to move.

## The math (first principles)

### Homogeneous coordinates

A 3D point $p = (x, y, z)$ is written as a 4-vector by appending $w$:

$$\tilde{p} = \begin{bmatrix} x \\ y \\ z \\ w \end{bmatrix}.$$

Two conventions, and the distinction is load-bearing:

- A **point** (a location) uses $w = 1$: $\tilde{p} = [x, y, z, 1]^\top$. It is affected by translation.
- A **direction / free vector** (a displacement, axis, or velocity) uses $w = 0$: $\tilde{v} = [x, y, z, 0]^\top$. Translation does not move it.

For a general projective point with $w \neq 0$, the Euclidean coordinates are recovered by the perspective divide $(x/w, y/w, z/w)$. In pure $SE(3)$ work $w$ stays at $1$ or $0$; the divide matters once you reach [[06_camera-pinhole]]. The whole reason homogeneous coordinates exist: they make affine maps (rotation plus translation) and projective maps (perspective) into a single linear matrix multiply, so composition is just matrix product. (OpenCV's camera model states this explicitly.)

### The rigid transform $T \in SE(3)$

$SE(3)$ is the special Euclidean group: all rigid motions of 3D space (rotation plus translation, no scaling, no reflection, no shear). An element is the block matrix

$$T = \begin{bmatrix} R & t \\ 0^\top & 1 \end{bmatrix}, \qquad R \in SO(3),\; t \in \mathbb{R}^3,$$

where $R$ is a $3\times3$ rotation matrix and $t$ is the $3\times1$ translation. The membership conditions are exact:

$$R^\top R = I, \qquad \det(R) = +1.$$

The first says columns are orthonormal (a true rotation preserves lengths and angles). The second, $\det = +1$, rules out reflections (handedness flips). $SE(3)$ has six degrees of freedom: three for rotation, three for translation.

**Frame-naming convention.** Write $T_{sb}$ for "the configuration of frame {b} expressed in frame {s}." Its columns are frame {b}'s axes and origin written in {s} coordinates. This subscript discipline is the convention in *Modern Robotics* and it makes composition read like cancelling units.

### Transforming a point

To move a point $p_b$ given in frame {b} into frame {s}:

$$\tilde{p}_s = T_{sb}\, \tilde{p}_b, \qquad \begin{bmatrix} p_s \\ 1 \end{bmatrix} = \begin{bmatrix} R & t \\ 0^\top & 1 \end{bmatrix}\begin{bmatrix} p_b \\ 1 \end{bmatrix} = \begin{bmatrix} R\,p_b + t \\ 1 \end{bmatrix}.$$

The bottom row $[0,0,0,1]$ regenerates the $w=1$ flag, so the output is again a valid point. In Euclidean terms this is the affine map $p \mapsto Rp + t$.

### Transforming a direction

A free vector $v$ uses $w=0$, so the translation column is multiplied by $0$ and drops out:

$$\begin{bmatrix} v_s \\ 0 \end{bmatrix} = \begin{bmatrix} R & t \\ 0^\top & 1 \end{bmatrix}\begin{bmatrix} v_b \\ 0 \end{bmatrix} = \begin{bmatrix} R\,v_b \\ 0 \end{bmatrix}.$$

A direction only rotates: $v \mapsto Rv$. This is exactly why the $w$ tag exists. The vector from point $a$ to point $b$ is $\tilde{b} - \tilde{a}$, whose $w$ component is $1 - 1 = 0$ automatically. The arithmetic enforces the semantics for you.

### Composition and inverse

Transforms chain by matrix multiplication, and the subscripts cancel:

$$T_{sa} = T_{sb}\, T_{ba}.$$

Read right to left: take something in {a}, express it in {b}, then in {s}. **Order matters**: $SE(3)$ is non-commutative, $T_1 T_2 \neq T_2 T_1$ in general. Rotating then translating is not the same as translating then rotating.

The inverse undoes the motion. You do **not** compute a generic $4\times4$ inverse; exploit the structure. Since $R^{-1} = R^\top$,

$$T^{-1} = \begin{bmatrix} R^\top & -R^\top t \\ 0^\top & 1 \end{bmatrix}.$$

Derivation: we need $T^{-1}T = I$. Write the inverse as $\begin{bmatrix} A & b \\ 0 & 1\end{bmatrix}$. Then $\begin{bmatrix} A & b \\ 0 & 1\end{bmatrix}\begin{bmatrix} R & t \\ 0 & 1\end{bmatrix} = \begin{bmatrix} AR & At + b \\ 0 & 1\end{bmatrix}$. For this to be $I$ we need $AR = I$ so $A = R^\top$, and $At + b = 0$ so $b = -R^\top t$. This closed form is faster and numerically cleaner than a general inverse, and in subscript notation $T_{bs} = (T_{sb})^{-1}$.

### Transforming normals: the inverse-transpose

A surface normal is **not** a free vector and must not be transformed with $R$ in the general affine case. A normal is defined by a perpendicularity condition: $n^\top v = 0$ for every tangent $v$ on the surface. We want the transformed normal $n'$ to stay perpendicular to the transformed tangent $v' = Mv$, where $M$ is the linear part. Require $n'^\top v' = 0$:

$$n'^\top (M v) = 0 \quad\text{for all } v\text{ with } n^\top v = 0.$$

Choose $n' = (M^{-1})^\top n$. Then $n'^\top M v = n^\top M^{-1} M v = n^\top v = 0$. So the correct rule is

$$n' = (M^{-1})^\top\, n, \qquad \text{then renormalize } n' \leftarrow n'/\lVert n' \rVert.$$

This is the **inverse-transpose** (a.k.a. the cofactor / normal matrix). The key special case: for a **pure rotation** $M = R$, we have $(R^{-1})^\top = (R^\top)^\top = R$, so normals rotate with the same $R$ as everything else. The inverse-transpose only differs from $M$ when there is non-uniform scale or shear. Stone scan pipelines are rigid-only most of the time, so $n' = R n$ holds; but the moment you introduce anisotropic scaling (unit conversion on one axis, a stretch correction) you must switch to the inverse-transpose or your lighting, offsets, and boolean orientations silently corrupt.

### Conventions to pin down before any code

- **Handedness.** Right-handed throughout ($z = x \times y$). Rhino, OpenCASCADE, and IFC are right-handed. OpenCV's camera frame is also right-handed but with +z forward and +y down; mismatches here cause mirrored scans. $\det(R) = +1$ guards against accidental flips.
- **Units.** Pick one length unit per frame and document it. Mixing meters (site) and millimeters (shop) is a classic $1000\times$ blow-up. The rotation block is unitless; only $t$ carries length.
- **Active vs passive.** $T_{sb}\tilde p_b$ here is a passive change of coordinates (same point, new frame). The identical matrix can be read actively as "move the point." Decide which story a given $T$ tells and stay consistent.
- **Storage order.** Eigen is column-major by default; OpenGL expects column-major; raw row-major buffers from some I/O need transposing. A "transposed rotation" bug usually traces here.

## Worked example

Frame {b} is rotated $+90°$ about the world z-axis and its origin sits at $(2, 0, 0)$ in world frame {s}. Build $T_{sb}$, transform a point and a direction, then invert.

Rotation by $\theta = 90°$ about z: $\cos 90° = 0$, $\sin 90° = 1$, so

$$R = \begin{bmatrix} 0 & -1 & 0 \\ 1 & 0 & 0 \\ 0 & 0 & 1 \end{bmatrix}, \qquad t = \begin{bmatrix} 2 \\ 0 \\ 0 \end{bmatrix}, \qquad T_{sb} = \begin{bmatrix} 0 & -1 & 0 & 2 \\ 1 & 0 & 0 & 0 \\ 0 & 0 & 1 & 0 \\ 0 & 0 & 0 & 1 \end{bmatrix}.$$

**Transform a point.** Take $p_b = (1, 0, 0)$ (one unit along {b}'s x-axis). Then

$$p_s = R\,p_b + t = \begin{bmatrix} 0\cdot1 - 1\cdot0 + 0 \\ 1\cdot1 + 0 + 0 \\ 0 \end{bmatrix} + \begin{bmatrix} 2 \\ 0 \\ 0 \end{bmatrix} = \begin{bmatrix} 0 \\ 1 \\ 0 \end{bmatrix} + \begin{bmatrix} 2 \\ 0 \\ 0 \end{bmatrix} = \begin{bmatrix} 2 \\ 1 \\ 0 \end{bmatrix}.$$

Sanity check: {b}'s x-axis points along world +y after the $90°$ spin, and the origin offset adds $(2,0,0)$. Result $(2,1,0)$. Correct.

**Transform a direction.** The same arrow $v_b = (1,0,0)$ as a free vector ($w=0$):

$$v_s = R\,v_b = (0, 1, 0).$$

The translation is gone. A heading points along world +y, with no $(2,0,0)$ shift. This is the point-vs-direction distinction in one line.

**Invert.** Using $T^{-1} = [R^\top, -R^\top t; 0, 1]$:

$$R^\top = \begin{bmatrix} 0 & 1 & 0 \\ -1 & 0 & 0 \\ 0 & 0 & 1 \end{bmatrix}, \qquad -R^\top t = -\begin{bmatrix} 0 & 1 & 0 \\ -1 & 0 & 0 \\ 0 & 0 & 1 \end{bmatrix}\begin{bmatrix} 2 \\ 0 \\ 0 \end{bmatrix} = -\begin{bmatrix} 0 \\ -2 \\ 0 \end{bmatrix} = \begin{bmatrix} 0 \\ 2 \\ 0 \end{bmatrix}.$$

$$T_{bs} = T_{sb}^{-1} = \begin{bmatrix} 0 & 1 & 0 & 0 \\ -1 & 0 & 0 & 2 \\ 0 & 0 & 1 & 0 \\ 0 & 0 & 0 & 1 \end{bmatrix}.$$

Verify by mapping $p_s = (2,1,0)$ back: $R^\top p_s - R^\top t = (1, -2, 0) + (0, 2, 0) = (1, 0, 0) = p_b$. Round-trip recovered. The arithmetic checks.

## Code

Idiomatic C++ with Eigen first. `Eigen::Isometry3d` is the natural type for a rigid transform: it stores a $4\times4$ but knows the bottom row is $[0,0,0,1]$, and its `.inverse()` uses the cheap $R^\top$ shortcut automatically.

```cpp
// se3_demo.cpp
// Build: g++ -std=c++17 -I /path/to/eigen se3_demo.cpp -o se3_demo
#include <Eigen/Geometry>
#include <iostream>

int main() {
    using Eigen::Isometry3d;   // 4x4 rigid transform, last row fixed [0 0 0 1]
    using Eigen::Vector3d;
    using Eigen::AngleAxisd;

    // T_sb: frame {b} = rotate +90 deg about world z, origin at (2,0,0) in {s}.
    Isometry3d T_sb = Isometry3d::Identity();
    T_sb.linear()      = AngleAxisd(M_PI / 2.0, Vector3d::UnitZ()).toRotationMatrix(); // R
    T_sb.translation() = Vector3d(2.0, 0.0, 0.0);                                      // t

    // Transform a POINT (w = 1 is implicit for Isometry * Vector3d):
    Vector3d p_b(1.0, 0.0, 0.0);
    Vector3d p_s = T_sb * p_b;                 // p_s = R*p_b + t -> (2,1,0)

    // Transform a DIRECTION (w = 0): apply only the linear/rotation block.
    Vector3d v_b(1.0, 0.0, 0.0);
    Vector3d v_s = T_sb.linear() * v_b;        // v_s = R*v_b      -> (0,1,0)

    // Transform a NORMAL: inverse-transpose of the linear part, then renormalize.
    // For a pure rotation this equals R, but write it the general way once and reuse.
    Eigen::Matrix3d N = T_sb.linear().inverse().transpose(); // (M^-1)^T
    Vector3d n_b(1.0, 0.0, 0.0);
    Vector3d n_s = (N * n_b).normalized();     // -> (0,1,0)

    // Inverse via the SE(3) shortcut (Eigen uses R^T internally, not a full 4x4 inverse):
    Isometry3d T_bs = T_sb.inverse();          // [R^T  -R^T t]
    Vector3d back = T_bs * p_s;                // -> (1,0,0), the original p_b

    // Composition: T_sa = T_sb * T_ba (subscripts cancel, right-to-left).
    Isometry3d T_ba = Isometry3d::Identity();
    T_ba.translation() = Vector3d(0.0, 0.0, 5.0);
    Isometry3d T_sa = T_sb * T_ba;             // chain frames

    std::cout << "p_s   = " << p_s.transpose()   << "\n";   // 2 1 0
    std::cout << "v_s   = " << v_s.transpose()   << "\n";   // 0 1 0
    std::cout << "n_s   = " << n_s.transpose()   << "\n";   // 0 1 0
    std::cout << "back  = " << back.transpose()  << "\n";   // 1 0 0
    std::cout << "T_sa.t= " << T_sa.translation().transpose() << "\n"; // 2 0 5
    return 0;
}
```

Mapping math to code: `T_sb.linear()` is $R$, `T_sb.translation()` is $t$, `T_sb * p_b` is $Rp_b + t$ (point, $w=1$), `T_sb.linear() * v_b` is $Rv_b$ (direction, $w=0$), `.inverse().transpose()` is $(M^{-1})^\top$ (normal), and `T_sb.inverse()` is the structured $SE(3)$ inverse.

Bridge to OCCT, where the same $T$ becomes a `gp_Trsf` and is the math behind every `IfcLocalPlacement`:

```cpp
// occt_placement.cpp  -- gp_Trsf is OCCT's rigid (+uniform-scale) transform.
#include <gp_Trsf.hxx>
#include <gp_Ax3.hxx>
#include <gp_Pnt.hxx>
#include <gp_Dir.hxx>

gp_Trsf placeBlockFrame() {
    // Define block frame {b} by origin + local Z (normal) + local X (reference dir).
    gp_Ax3 blockFrame(gp_Pnt(2.0, 0.0, 0.0),  // origin t
                      gp_Dir(0.0, 0.0, 1.0),  // local Z
                      gp_Dir(0.0, 1.0, 0.0)); // local X -> 90 deg about Z
    gp_Trsf T;                 // identity == world frame {s}
    T.SetTransformation(blockFrame);          // builds T_bs (world -> block)
    // T.Inverted() gives T_sb (block -> world); pick the direction your data needs.
    return T;
}
```

A short Python check (NumPy) to inspect the numbers interactively:

```python
import numpy as np
c, s = 0.0, 1.0  # cos90, sin90
R = np.array([[c, -s, 0], [s, c, 0], [0, 0, 1.0]])
t = np.array([2.0, 0, 0])
T = np.eye(4); T[:3, :3] = R; T[:3, 3] = t
p = T @ np.array([1, 0, 0, 1.0])   # point, w=1 -> [2,1,0,1]
v = T @ np.array([1, 0, 0, 0.0])   # direction, w=0 -> [0,1,0,0]
Tinv = np.eye(4); Tinv[:3, :3] = R.T; Tinv[:3, 3] = -R.T @ t
print(p[:3], v[:3], np.allclose(Tinv @ T, np.eye(4)))  # [2 1 0] [0 1 0] True
```

## Connections

- [[03_matrices]]: $T$ is a block matrix; composition is matrix product and the inverse uses orthogonality $R^{-1}=R^\top$. Without matrix fluency none of this reads.
- [[04_trig-rotations]]: the $R$ block is built from sines and cosines; $\det R = +1$ and $R^\top R = I$ come straight from rotation trig.
- [[06_camera-pinhole]]: projection is $s\,\tilde p = K[R\mid t]\tilde P_w$; the $[R\mid t]$ is exactly this rigid transform, and the perspective divide is the $w \neq 1$ generalization of homogeneous coordinates.
- [[07_forward-kinematics]]: a robot pose is a product $T = T_1 T_2 \cdots T_n$ of $SE(3)$ links; FK is repeated composition of this object.
- [[18_brep-boolean-occt]]: OCCT places every face with a `gp_Trsf` (an $SE(3)$ element); boolean operations depend on consistent face normals, which is where the inverse-transpose rule bites.
- [[08_least-squares-svd]] / ICP: rigid alignment solves for the $R, t$ of an $SE(3)$ element that best matches two scans, so $SE(3)$ is the unknown being optimized.

## Pitfalls and failure modes

- **Frame inconsistency.** Mixing $T_{sb}$ and $T_{bs}$ is the top bug. Use subscript naming so $T_{sb}T_{bc}$ visibly cancels; if subscripts do not chain, the product is wrong. Most "vision/robotics/CAD math bugs" are actually frame bugs.
- **Point treated as direction (or vice versa).** Forgetting $w=1$ on a point drops the translation; forgetting $w=0$ on a normal adds a spurious offset. Always tag explicitly.
- **Transforming normals with $R$ under non-uniform scale.** Fine for pure rotation, silently wrong once a non-uniform scale or shear enters. Use $(M^{-1})^\top$ and renormalize. Symptom: lighting goes dark on stretched geometry, or boolean cuts flip inside/out.
- **Computing the $4\times4$ inverse generically.** Slower and less accurate. Use $[R^\top, -R^\top t; 0, 1]$. A general solver can also drift $R$ off $SO(3)$.
- **Non-orthonormal $R$ drift.** After many compositions, floating-point error makes $R^\top R \neq I$. Re-orthonormalize periodically (Gram-Schmidt or the SVD projection $R \leftarrow UV^\top$). Skipping this slowly skews long transform chains.
- **Unit mismatch.** Meters vs millimeters on $t$ produces $1000\times$ errors that pass type-checks. Fix one unit per frame and record it next to the matrix.
- **Non-commutativity ignored.** $T_1 T_2 \neq T_2 T_1$. Translate-then-rotate differs from rotate-then-translate. Build the chain in the physical order of operations.
- **Handedness flip.** A reflection sneaks in ($\det R = -1$) from a bad axis cross or a mirrored import; geometry comes out mirrored. Assert $\det R \approx +1$.
- **Storage-order transpose.** Row-major vs column-major buffers across an API boundary give a transposed (hence inverse) rotation. Verify on a known $90°$ case.

## Practice

1. By hand, build $T_{sb}$ for a $30°$ rotation about the y-axis with origin at $(0, 1, 0)$. Transform the point $(1, 0, 0)$ and the direction $(1, 0, 0)$. Confirm the direction ignores the translation. (Inspectable: two 3-vectors.)
2. Prove algebraically that $(T_{sb})^{-1} = T_{bs} = [R^\top,\, -R^\top t;\, 0, 1]$ by multiplying out $T_{bs}T_{sb}$ and showing it equals $I_4$.
3. In Eigen, compose three transforms $T_{sa} = T_{sb} T_{bc} T_{cd}$ with one $90°$ rotation in each about a different axis. Print the final rotation's determinant and confirm it is $+1$ (no reflection crept in).
4. Take the unit normal $(0, 0, 1)$ on a plane and apply a non-uniform scale $\mathrm{diag}(2, 1, 1)$ first with $R$ (wrong) and then with $(M^{-1})^\top$ (right). Show the wrong version is no longer perpendicular to a tangent like $(2, 0, 0)$ after scaling.
5. Write a function that takes two paired point sets and returns the $SE(3)$ that best aligns them (centroid subtraction + SVD of the covariance, fixing $\det(VU^\top)$). Verify it recovers a known $T$ you applied to a scan. (Inspectable: recovered $R, t$ vs ground truth.) This is the seed of ICP and connects to [[08_least-squares-svd]].
6. Re-orthonormalize a slightly corrupted rotation matrix (add $10^{-3}$ noise to each entry) using the SVD projection $R \leftarrow U V^\top$ and report $\lVert R^\top R - I\rVert$ before and after.

## References

- Kevin M. Lynch and Frank C. Park, *Modern Robotics: Mechanics, Planning, and Control*, Cambridge University Press, 2017. Chapter 3 (rigid-body motions, $SE(3)$, homogeneous transforms, the inverse $[R^\top, -R^\top t; 0,1]$, the three uses of a transform). Free PDF and course site: https://modernrobotics.northwestern.edu/ and https://hades.mech.northwestern.edu/index.php/Modern_Robotics
- Modern Robotics resource page, 3.3.1 Homogeneous Transformation Matrices: https://modernrobotics.northwestern.edu/nu-gm-book-resource/3-3-1-homogeneous-transformation-matrices/
- OpenCV documentation, Camera Calibration and 3D Reconstruction (calib3d), pinhole model $s\,p = K[R\mid t]P_w$ and the role of homogeneous coordinates: https://docs.opencv.org/4.x/d9/d0c/group__calib3d.html
- Eigen documentation, Space transformations tutorial (`Isometry3d`, `Affine3d`, `.linear()`, `.translation()`, transforming points vs vectors): https://eigen.tuxfamily.org/dox/group__TutorialGeometry.html
- Scratchapixel, "Transforming Normals" (derivation of the inverse-transpose rule): https://www.scratchapixel.com/lessons/mathematics-physics-for-computer-graphics/geometry/transforming-normals.html
- Nathan Reed, "Normals and the Inverse Transpose" (geometric justification): https://www.reedbeta.com/blog/normals-inverse-transpose-part-1/
- OpenCASCADE Technology reference manual, `gp_Trsf` / `gp_Ax3` and Foundation Classes (transforms in the CAD kernel): https://dev.opencascade.org/doc/refman/html/classgp___trsf.html
- IfcOpenShell documentation, geometry and `IfcLocalPlacement` (placements as rigid transforms): https://docs.ifcopenshell.org/
- Peter Corke and Jesse Haviland, *Robotics, Vision and Control* / Spatial Math toolbox (pose as $SE(3)$): https://petercorke.com/

## Flashcard seeds

- Q: Derive the action of $T$ on a point. Start from the block form and multiply it by the homogeneous point. <br> <b>Formula:</b> $T=\begin{bmatrix} R & t \\ 0^\top & 1\end{bmatrix}$, $\tilde p_b=\begin{bmatrix} p_b \\ 1\end{bmatrix}$, target $p_s = R\,p_b + t$ <br> <i>Hint:</i> multiply block-by-block; the top block row is $R\,p_b + t\cdot 1$, the bottom row regenerates $w$. :: A: Top row: $R\,p_b + t\cdot 1 = R\,p_b + t$ <br> bottom row: $0^\top p_b + 1\cdot 1 = 1$, which regenerates the $w=1$ flag <br> $\boxed{p_s = R\,p_b + t}$
- Q: Derive the action of $T$ on a free direction and show why translation drops out. <br> <b>Formula:</b> $\tilde v_b=\begin{bmatrix} v_b \\ 0\end{bmatrix}$, $T=\begin{bmatrix} R & t \\ 0^\top & 1\end{bmatrix}$, target $v_s = R\,v_b$ <br> <i>Hint:</i> the $w=0$ tag multiplies the translation column $t$ by $0$; compute the top row first. :: A: Top row: $R\,v_b + t\cdot 0 = R\,v_b$ <br> the $w=0$ tag zeroes the translation column $t$ <br> bottom row: $0^\top v_b + 1\cdot 0 = 0$, preserving the direction tag <br> $\boxed{v_s = R\,v_b}$
- Q: Derive the $SE(3)$ inverse, step 1 of 2: set up the block equation for the unknown inverse. <br> <b>Formula:</b> require $\begin{bmatrix} A & b \\ 0^\top & 1\end{bmatrix}\begin{bmatrix} R & t \\ 0^\top & 1\end{bmatrix}=I$ <br> <i>Hint:</i> multiply the two blocks out, then match each block of the product to the identity. :: A: $\begin{bmatrix} A & b \\ 0^\top & 1\end{bmatrix}\begin{bmatrix} R & t \\ 0^\top & 1\end{bmatrix}=\begin{bmatrix} AR & At+b \\ 0^\top & 1\end{bmatrix}$ <br> match to $I$ (top-left $=I$, top-right $=0$) <br> $\boxed{AR=I \text{ and } At+b=0}$
- Q: Derive the $SE(3)$ inverse, step 2 of 2: from $AR=I$ and $At+b=0$, solve for $A$ and $b$. <br> <b>Formula:</b> orthogonality gives $R^{-1}=R^\top$; target $T^{-1}=\begin{bmatrix} R^\top & -R^\top t \\ 0^\top & 1\end{bmatrix}$ <br> <i>Hint:</i> get $A$ from the first equation, then substitute it into the second to isolate $b$. :: A: From $AR=I$: $A=R^{-1}=R^\top$ (orthogonality) <br> From $At+b=0$: $b=-At=-R^\top t$ <br> $\boxed{T^{-1}=\begin{bmatrix} R^\top & -R^\top t \\ 0^\top & 1\end{bmatrix}}$
- Q: Derive the normal-transform rule, step 1 of 2: set up the condition that perpendicularity must survive. <br> <b>Formula:</b> a tangent maps as $v'=Mv$; require $n'^\top v'=0$ for every $v$ with $n^\top v=0$ <br> <i>Hint:</i> substitute $v'=Mv$ into $n'^\top v'=0$, then move $M$ onto $n'$ using $(n'^\top M)=(M^\top n')^\top$. :: A: Demand perpendicularity survives: $n'^\top (Mv)=0$ for all $v$ with $n^\top v=0$ <br> rewrite as $(M^\top n')^\top v = 0$ <br> $\boxed{M^\top n' \parallel n \text{ (must be proportional to } n)}$
- Q: Derive the normal-transform rule, step 2 of 2: solve for $n'$ so the constraint holds for all tangents. <br> <b>Formula:</b> need $n'^\top M v = n^\top v$ for all $v$; target $n'=(M^{-1})^\top n$ <br> <i>Hint:</i> guess $n'=(M^{-1})^\top n$ and check it; note $n'^\top M = n^\top M^{-1} M = n^\top$. :: A: Try $n'=(M^{-1})^\top n$ <br> then $n'^\top M v = n^\top M^{-1} M v = n^\top v = 0$, satisfying the constraint <br> $\boxed{n'=(M^{-1})^\top n,\ \text{then } n'\leftarrow n'/\lVert n'\rVert}$
- Q: Derive what the normal rule becomes for a pure rotation $M=R$. <br> <b>Formula:</b> $n'=(M^{-1})^\top n$ with $M=R$ and $R^{-1}=R^\top$ <br> <i>Hint:</i> substitute $M=R$, then apply $R^{-1}=R^\top$ and simplify the double transpose $(R^\top)^\top$. :: A: Substitute $M=R$: $n'=(R^{-1})^\top n$ <br> orthogonality gives $R^{-1}=R^\top$, so $(R^{-1})^\top=(R^\top)^\top=R$ <br> $\boxed{n'=R\,n}$ — normals rotate with the same $R$
- Q: Derive the composition rule $T_{sa}=T_{sb}T_{ba}$ by tracking one point through the frames. <br> <b>Formula:</b> $\tilde p_b = T_{ba}\tilde p_a$, then $\tilde p_s = T_{sb}\tilde p_b$; also $\tilde p_s = T_{sa}\tilde p_a$ <br> <i>Hint:</i> substitute the first map into the second, then equate with the direct map $T_{sa}\tilde p_a$. :: A: Start with $\tilde p_a$ in {a}; map to {b}: $\tilde p_b = T_{ba}\tilde p_a$ <br> map to {s}: $\tilde p_s = T_{sb}\tilde p_b = T_{sb}T_{ba}\tilde p_a$ <br> but also $\tilde p_s = T_{sa}\tilde p_a$, so $\boxed{T_{sa}=T_{sb}T_{ba}}$ (subscripts cancel)
- Q: Build $T_{sb}$ for frame {b} = a $+90°$ rotation about world z with origin at $(2,0,0)$. <br> <b>Formula:</b> $R=\begin{bmatrix} \cos\theta & -\sin\theta & 0 \\ \sin\theta & \cos\theta & 0 \\ 0 & 0 & 1\end{bmatrix}$, $T_{sb}=\begin{bmatrix} R & t \\ 0^\top & 1\end{bmatrix}$ <br> <i>Hint:</i> plug $\theta=90°$ ($\cos=0,\ \sin=1$) into $R$ first, then drop $t=(2,0,0)^\top$ into the last column. :: A: $\cos 90°=0,\ \sin 90°=1$ so $R=\begin{bmatrix} 0 & -1 & 0 \\ 1 & 0 & 0 \\ 0 & 0 & 1\end{bmatrix}$ <br> stack $t=(2,0,0)^\top$ in the last column <br> $\boxed{T_{sb}=\begin{bmatrix} 0 & -1 & 0 & 2 \\ 1 & 0 & 0 & 0 \\ 0 & 0 & 1 & 0 \\ 0 & 0 & 0 & 1\end{bmatrix}}$
- Q: Worked example (point, 3D): with the $90°$-about-z rotation $R=\begin{bmatrix} 0 & -1 & 0 \\ 1 & 0 & 0 \\ 0 & 0 & 1\end{bmatrix}$ and $t=(2,0,0)$, transform $p_b=(1,0,0)$. <br> <b>Formula:</b> $p_s = R\,p_b + t$ <br> <i>Hint:</i> compute $R\,p_b$ first (row by row), then add $t$ component-wise. :: A: $R\,p_b=(0\cdot1-1\cdot0,\ 1\cdot1+0,\ 0)=(0,1,0)$ <br> add $t$: $(0,1,0)+(2,0,0)$ <br> $\boxed{p_s=(2,1,0)}$
- Q: Worked example (direction, 3D): same $R=\begin{bmatrix} 0 & -1 & 0 \\ 1 & 0 & 0 \\ 0 & 0 & 1\end{bmatrix}$ and $t=(2,0,0)$, transform the free vector $v_b=(1,0,0)$ ($w=0$). <br> <b>Formula:</b> $v_s = R\,v_b$ (translation skipped because $w=0$) <br> <i>Hint:</i> apply only $R$; do not add $t$. :: A: Only $R$ acts: $v_s=R\,v_b=(0,1,0)$ <br> translation $(2,0,0)$ is skipped because $w=0$ <br> $\boxed{v_s=(0,1,0)}$
- Q: Worked example (small scale, mm): a scanner frame is offset $t=(5,0,0)\,\text{mm}$ with no rotation ($R=I$). Transform the point $p_b=(2,3,0)\,\text{mm}$. <br> <b>Formula:</b> $p_s = R\,p_b + t = I\,p_b + t$ <br> <i>Hint:</i> with $R=I$ this is pure translation — just add $t$ component-wise. :: A: $p_s = I\,p_b + t = (2,3,0)+(5,0,0)$ <br> $\boxed{p_s=(7,3,0)\,\text{mm}}$ — pure translation, all axes carry mm
- Q: Worked example (large scale, m): a block frame is $-90°$ about z with origin $t=(10,4,0)\,\text{m}$ in the site frame. Transform $p_b=(0,2,0)\,\text{m}$. <br> <b>Formula:</b> for $-90°$, $R=\begin{bmatrix} 0 & 1 & 0 \\ -1 & 0 & 0 \\ 0 & 0 & 1\end{bmatrix}$; $p_s = R\,p_b + t$ <br> <i>Hint:</i> compute $R\,p_b$ first, then add $t$; remember $\sin(-90°)=-1$ flips the sign pattern. :: A: $R\,p_b=(0\cdot0+1\cdot2,\ -1\cdot0+0,\ 0)=(2,0,0)$ <br> add $t$: $(2,0,0)+(10,4,0)$ <br> $\boxed{p_s=(12,4,0)\,\text{m}}$
- Q: Worked example (inverse): invert the $90°$-about-z transform with $R=\begin{bmatrix} 0 & -1 & 0 \\ 1 & 0 & 0 \\ 0 & 0 & 1\end{bmatrix}$, $t=(2,0,0)$. <br> <b>Formula:</b> $T^{-1}=\begin{bmatrix} R^\top & -R^\top t \\ 0^\top & 1\end{bmatrix}$ <br> <i>Hint:</i> transpose $R$ first (swap off-diagonals), then compute $R^\top t$ and negate it. :: A: $R^\top=\begin{bmatrix} 0 & 1 & 0 \\ -1 & 0 & 0 \\ 0 & 0 & 1\end{bmatrix}$ <br> $R^\top t = (0,-2,0)$, so $-R^\top t = (0,2,0)$ <br> $\boxed{T_{bs}=\begin{bmatrix} 0 & 1 & 0 & 0 \\ -1 & 0 & 0 & 2 \\ 0 & 0 & 1 & 0 \\ 0 & 0 & 0 & 1\end{bmatrix}}$
- Q: Closing test (round-trip): verify the inverse by mapping $p_s=(2,1,0)$ back to {b} and checking you recover $p_b=(1,0,0)$. Use $R=\begin{bmatrix} 0 & -1 & 0 \\ 1 & 0 & 0 \\ 0 & 0 & 1\end{bmatrix}$, $t=(2,0,0)$. <br> <b>Formula:</b> $p_b = R^\top p_s - R^\top t$ <br> <i>Hint:</i> a pass means the result equals the original $p_b=(1,0,0)$ exactly. :: A: $R^\top p_s = (1,-2,0)$ and $-R^\top t = (0,2,0)$ <br> $p_b = (1,-2,0)+(0,2,0) = (1,0,0)$, the original $p_b$ <br> $\boxed{\text{round-trip recovered} \Rightarrow T_{bs}T_{sb}=I}$
- Q: Closing test (normal under scale): a normal $(0,0,1)$ undergoes the scale $M=\mathrm{diag}(2,1,1)$. Check the inverse-transpose normal stays perpendicular to the transformed tangent $(2,0,0)$. <br> <b>Formula:</b> $n'=(M^{-1})^\top n$, $v'=Mv$, pass when $n'^\top v'=0$ <br> <i>Hint:</i> a pass is $n'^\top v'=0$; compute $(M^{-1})^\top$ (reciprocate the diagonal) first. :: A: $(M^{-1})^\top=\mathrm{diag}(1/2,1,1)$, so $n'=(0,0,1)$ (unchanged) <br> transform tangent: $v'=Mv=(4,0,0)$ <br> check $n'^\top v' = 0\cdot4+0\cdot0+1\cdot0 = 0$ <br> $\boxed{n'\perp v' \Rightarrow \text{rule verified; naive } Rn \text{ would tilt } n}$
- Q: Closing test (determinant): after composing several $SE(3)$ transforms you suspect a reflection crept in. State the check that confirms $R$ is still a proper rotation. <br> <b>Formula:</b> proper rotation needs $\det(R)=+1$ (and $R^\top R=I$) <br> <i>Hint:</i> a pass is $\det(R)\approx +1$; $\det(R)=-1$ means a handedness flip. :: A: Compute $\det(R)$ of the composed rotation block <br> proper rotation $\Rightarrow \det(R)=+1$; reflection $\Rightarrow \det(R)=-1$ <br> $\boxed{\text{pass when } \det(R)\approx +1}$
- Q: What is the block-matrix form of a rigid transform $T \in SE(3)$, and what are the membership conditions? <br> <i>Hint:</i> rotation block on top-left, translation column on the right, fixed bottom row. :: A: $T = \begin{bmatrix} R & t \\ 0^\top & 1\end{bmatrix}$ with $R^\top R=I$, $\det R=+1$, $t\in\mathbb{R}^3$.
- Q: When must you use the inverse-transpose $(M^{-1})^\top$ instead of $M$ for normals? <br> <i>Hint:</i> think about what breaks orthogonality — scale that differs per axis, or shear. :: A: Whenever $M$ has non-uniform scale or shear; for a pure rotation (or uniform scale) $R$ suffices because $(R^{-1})^\top=R$.
- Q: Why use the structured inverse $[R^\top, -R^\top t]$ instead of a generic $4\times4$ inverse? <br> <i>Hint:</i> think speed, numerical accuracy, and staying on $SO(3)$. :: A: Faster, more accurate, and it cannot drift $R$ off $SO(3)$.
- Q: What does $T_{sb}$ mean in Modern Robotics subscript notation? <br> <i>Hint:</i> read it as "frame {b} expressed in frame {s}"; its columns are {b}'s axes and origin. :: A: The configuration of frame {b} expressed in frame {s} (its columns are {b}'s axes and origin in {s}).
- CLOZE: A point in homogeneous coordinates uses $w = {{c1::1}}$ and a free direction uses $w = {{c2::0}}$.
- CLOZE: $SE(3)$ composition is {{c1::non-commutative}}, so $T_1 T_2 \neq T_2 T_1$ in general.
- CLOZE: After many compositions, re-orthonormalize $R$ via the SVD projection {{c1::$R \leftarrow UV^\top$}} to undo floating-point drift.
- CLOZE: The most common class of transform bug is a {{c1::frame-consistency}} bug — using $T_{sb}$ where $T_{bs}$ is needed.
- CLOZE: A reflection has sneaked in when {{c1::$\det(R)=-1$}}; assert $\det(R)\approx +1$ to catch a mirrored import.
