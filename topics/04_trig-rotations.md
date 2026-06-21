# Trigonometry and rotations (SO(3), axis-angle, quaternions)

Prereqs: [[03_vectors]] => this topic => unlocks [[05_transforms-se3]], [[06_twists-screws-exp]].

## Why it matters for the stone scan->CAD->IFC / CADFit workflow

Every stage of the stone pipeline moves geometry between frames, and a rotation sits inside each move. When you register a raw scan to a site frame or template, ICP solves for a rotation $R$ plus a translation $t$; the rotation is the hard part, and it must stay a valid member of $SO(3)$ or the fitted CAD surface drifts. When you place a fitted B-Rep block into an IFC model, the `IfcLocalPlacement` is a rotation plus an origin, so a wrong handedness or a non-orthonormal matrix produces a silently mirrored or skewed stone. When a CADFit-style synthesis layer emits a CAD program (extrude this profile, rotate it by this axis-angle, boolean-subtract), the rotation parameters it predicts must be expressible in the exact convention OCCT consumes. Getting rotations right, and knowing which representation to use where, is the difference between a pipeline that closes the loop and one that accumulates frame bugs you cannot see until the stone is cut.

## Intuition

A rotation is a rigid turn about a fixed point. It keeps every length and every angle, and it does not flip the object into its mirror image. In 3D, Euler's rotation theorem says any such turn, no matter how complicated, is equivalent to a single spin by some angle $\theta$ about some axis $\hat\omega$. That is the whole story: pick an axis, pick an angle, turn. Everything else (matrices, quaternions, Euler angles) is just a different way to write down "this axis, this angle."

Analogy: think of spinning a globe. You can describe where it ends up by three separate dial turns (Euler angles), but the globe itself only ever did one clean spin about the pole through the rotation axis. The single-spin view (axis-angle and its cousin the quaternion) is the honest one; the three-dial view is convenient but can jam (gimbal lock). Order matters too: spin the globe 90 degrees about the vertical, then 90 about the horizontal, and you land somewhere different than if you swap the two. Rotation composition does not commute.

## The math (first principles)

### Conventions used here

- Right-handed coordinate frames. Positive angle is counterclockwise when looking down the axis toward the origin (right-hand rule).
- Vectors are columns. A rotation acts as $v' = R v$. This is the active convention: the matrix moves the point, the frame stays fixed.
- Angles in radians. $\hat\omega$ denotes a unit vector ($\lVert\hat\omega\rVert = 1$).
- Quaternion order is $q = (w, x, y, z)$ with $w$ the scalar part. Hamilton convention.
- Tolerances are relative. Treat a vector as zero-length below $\varepsilon = 10^{-12}$; treat $\theta$ near $0$ or $\pi$ as a degeneracy needing special handling.

### Trig refresher

For a right triangle with angle $\theta$, $\sin\theta = \text{opp}/\text{hyp}$, $\cos\theta = \text{adj}/\text{hyp}$, $\tan\theta = \sin\theta/\cos\theta$. On the unit circle a point at angle $\theta$ is $(\cos\theta, \sin\theta)$. The identities you actually reuse:

$$\sin^2\theta + \cos^2\theta = 1,$$
$$\cos(\alpha+\beta) = \cos\alpha\cos\beta - \sin\alpha\sin\beta, \qquad \sin(\alpha+\beta) = \sin\alpha\cos\beta + \cos\alpha\sin\beta.$$

The angle-sum identities are not trivia. They are exactly why composing two 2D rotations adds their angles, which you can verify by multiplying the matrices below. Similar triangles (equal angles imply proportional sides) underlie every projection and scaling step downstream.

### 2D rotation matrix

A rotation by $\theta$ in the plane maps the basis vectors $\hat e_1 \mapsto (\cos\theta, \sin\theta)$ and $\hat e_2 \mapsto (-\sin\theta, \cos\theta)$. Stacking those images as columns:

$$R(\theta) = \begin{bmatrix} \cos\theta & -\sin\theta \\ \sin\theta & \cos\theta \end{bmatrix}.$$

Properties: $R(\theta)^{-1} = R(-\theta) = R(\theta)^{\mathsf T}$, $\det R = 1$, and $R(\alpha)R(\beta) = R(\alpha+\beta)$. The last identity is the angle-sum formula in disguise. In 2D rotations commute because they all share the same (out-of-plane) axis.

### The rotation group SO(3)

A 3D rotation is a $3\times3$ real matrix $R$ satisfying two conditions:

$$R^{\mathsf T} R = I \quad\text{(orthonormal)}, \qquad \det R = +1 \quad\text{(no reflection)}.$$

The set of all such matrices is the special orthogonal group $SO(3)$. Orthonormality means the columns are mutually perpendicular unit vectors, so $R$ maps an orthonormal frame to another orthonormal frame. The $\det = +1$ condition rules out reflections, which would have $\det = -1$. $SO(3)$ is a group: the product of two rotations is a rotation, the inverse of a rotation is a rotation ($R^{-1} = R^{\mathsf T}$), and the identity $I$ is the no-op rotation. It is a 3-dimensional manifold, not a flat vector space, which is why you cannot just add two rotation matrices and expect a rotation.

### Axis-angle and the Rodrigues formula

By Euler's theorem, $R$ corresponds to a unit axis $\hat\omega$ and angle $\theta$. Let $[\hat\omega]$ be the skew-symmetric cross-product matrix:

$$[\hat\omega] = \begin{bmatrix} 0 & -\omega_z & \omega_y \\ \omega_z & 0 & -\omega_x \\ -\omega_y & \omega_x & 0 \end{bmatrix}, \qquad [\hat\omega]\,v = \hat\omega \times v.$$

Rodrigues' rotation formula gives the matrix:

$$R = I + \sin\theta\,[\hat\omega] + (1-\cos\theta)\,[\hat\omega]^2.$$

This is verified against Modern Robotics and Wikipedia. It is the closed form of the matrix exponential $R = e^{[\hat\omega]\theta} = \exp(\theta[\hat\omega])$, which is why axis-angle is also called exponential coordinates. The equivalent vector form, rotating a single vector $v$ about $\hat\omega$, is:

$$v_{\text{rot}} = v\cos\theta + (\hat\omega \times v)\sin\theta + \hat\omega\,(\hat\omega \cdot v)(1-\cos\theta).$$

Derivation sketch: $[\hat\omega]^2 = \hat\omega\hat\omega^{\mathsf T} - I$ for a unit axis. Decompose $v$ into a part along $\hat\omega$ (unchanged) and a part perpendicular (which rotates by $\theta$ in its plane). The perpendicular part uses $\cos\theta$ along itself and $\sin\theta$ along $\hat\omega\times v$. Collecting terms gives the formula above.

The inverse direction (matrix logarithm) recovers axis-angle from $R$:

$$\theta = \arccos\!\left(\frac{\operatorname{tr}(R)-1}{2}\right), \qquad [\hat\omega] = \frac{1}{2\sin\theta}\,(R - R^{\mathsf T}).$$

Two degeneracies, both confirmed by Modern Robotics. When $\theta = 0$, $R = I$ and the axis is undefined (any axis works). When $\theta = \pi$, $\sin\theta = 0$ so the formula above divides by zero; recover the axis from the diagonal of $R$ instead, since $R = I + 2[\hat\omega]^2$ there and $\hat\omega_i^2 = (R_{ii}+1)/2$.

### Euler angles and gimbal lock

Euler angles factor $R$ into three rotations about coordinate axes, for example the $ZYX$ (yaw-pitch-roll) convention:

$$R = R_z(\alpha)\,R_y(\beta)\,R_x(\gamma).$$

There are twelve valid conventions (six Tait-Bryan like $ZYX$, six proper-Euler like $ZXZ$), so "Euler angles" alone is ambiguous; always state the axis order and whether it is intrinsic or extrinsic. The fatal flaw is gimbal lock: when the middle angle reaches its degenerate value (here $\beta = \pm90^\circ$), the first and third axes align, the rotation loses a degree of freedom, and $\alpha$ and $\gamma$ become indistinguishable. Near that configuration the conversion from $R$ back to angles is ill-conditioned. Use Euler angles for human-readable input and output, never as the internal state of a solver.

### Unit quaternions

A quaternion is $q = w + x i + y j + z k$ with $i^2 = j^2 = k^2 = ijk = -1$. Write it as scalar-vector $q = (w, \mathbf v)$ with $\mathbf v = (x,y,z)$. A rotation by angle $\theta$ about unit axis $\hat\omega$ is the unit quaternion:

$$q = \left(\cos\tfrac{\theta}{2},\ \hat\omega\,\sin\tfrac{\theta}{2}\right), \qquad \lVert q\rVert = 1.$$

Note the half-angle: this is verified against the standard references. A consequence is the double cover, $q$ and $-q$ represent the same rotation, because adding $2\pi$ to $\theta$ flips the sign of both $\cos(\theta/2)$ and $\sin(\theta/2)$.

Composition uses the Hamilton product. For $q_1=(w_1,\mathbf v_1)$ and $q_2=(w_2,\mathbf v_2)$:

$$q_1 q_2 = \big(w_1 w_2 - \mathbf v_1\cdot\mathbf v_2,\ \ w_1\mathbf v_2 + w_2\mathbf v_1 + \mathbf v_1\times\mathbf v_2\big).$$

This is non-commutative (the cross product term changes sign when you swap operands), exactly mirroring matrix multiplication. Applying $q$ to a vector $v$ uses the sandwich $v' = q\,(0,v)\,q^{-1}$, and for a unit quaternion $q^{-1} = (w, -\mathbf v)$ (the conjugate). The equivalent rotation matrix, verified entry-by-entry against euclideanspace.com, is:

$$R(q) = \begin{bmatrix} 1-2(y^2+z^2) & 2(xy-wz) & 2(xz+wy) \\ 2(xy+wz) & 1-2(x^2+z^2) & 2(yz-wx) \\ 2(xz-wy) & 2(yz+wx) & 1-2(x^2+y^2) \end{bmatrix}.$$

### Slerp

To interpolate between two orientations $q_0$ and $q_1$ with constant angular velocity, use spherical linear interpolation. With $\cos\Omega = q_0 \cdot q_1$ (the 4D dot product):

$$\operatorname{slerp}(q_0, q_1; t) = \frac{\sin((1-t)\Omega)}{\sin\Omega}\,q_0 + \frac{\sin(t\Omega)}{\sin\Omega}\,q_1, \qquad t\in[0,1].$$

This traces the shortest great-circle arc on the unit sphere $S^3$. Two practical rules, both confirmed by the references. First, if $q_0\cdot q_1 < 0$, negate one quaternion so you take the short way around (the double cover gives you a choice). Second, when $\Omega$ is tiny, $\sin\Omega \to 0$ and the formula is numerically unstable; fall back to normalized linear interpolation (nlerp).

## Worked example

Rotate the vector $v = (1, 0, 0)$ by $\theta = 90^\circ = \pi/2$ about the $z$-axis $\hat\omega = (0,0,1)$. We expect $(1,0,0)$ to land on $(0,1,0)$.

Rodrigues vector form. $\cos\theta = 0$, $\sin\theta = 1$, $\hat\omega\cdot v = 0$, and $\hat\omega \times v = (0,0,1)\times(1,0,0) = (0,1,0)$.

$$v_{\text{rot}} = v\cdot 0 + (0,1,0)\cdot 1 + (0,0,1)\cdot 0 \cdot 1 = (0,1,0). \checkmark$$

Matrix form. $[\hat\omega]^2$ for $\hat\omega=(0,0,1)$ is $\operatorname{diag}(-1,-1,0)$, so

$$R = I + 1\cdot[\hat\omega] + (1-0)[\hat\omega]^2 = \begin{bmatrix} 0 & -1 & 0 \\ 1 & 0 & 0 \\ 0 & 0 & 1 \end{bmatrix},$$

and $R(1,0,0)^{\mathsf T} = (0,1,0)^{\mathsf T}$. This is the 2D rotation by $90^\circ$ embedded in the $xy$-plane.

Quaternion form. $q = (\cos45^\circ, 0, 0, \sin45^\circ) = (\tfrac{\sqrt2}{2}, 0, 0, \tfrac{\sqrt2}{2})$. Plug $w=z=\tfrac{\sqrt2}{2}$, $x=y=0$ into $R(q)$. Top-left $= 1 - 2(0+\tfrac12) = 0$; entry $(0,1) = 2(0 - \tfrac{\sqrt2}{2}\cdot\tfrac{\sqrt2}{2}) = -1$; entry $(1,0) = 2(0 + \tfrac12) = 1$. Same matrix. $\checkmark$

Non-commutativity check. Let $A = R_z(90^\circ)$ and $B = R_x(90^\circ)$. Apply to $(1,0,0)$: $A$ then $B$ gives $B(A v) = B(0,1,0) = (0,0,1)$. Reverse order $A(B v) = A(1,0,0) = (0,1,0)$. Different results, so $AB \neq BA$.

## Code

Idiomatic C++ with Eigen. Eigen's `Geometry` module provides `AngleAxisd`, `Quaterniond`, and `Matrix3d`, and the conversions and slerp shown here match the official Eigen tutorial.

```cpp
#include <Eigen/Geometry>   // AngleAxisd, Quaterniond, Matrix3d
#include <iostream>
#include <cmath>

using Eigen::Vector3d;
using Eigen::Matrix3d;
using Eigen::Quaterniond;
using Eigen::AngleAxisd;

// Rodrigues by hand: R = I + sin(t)[w] + (1-cos t)[w]^2, w a UNIT axis.
// Maps math symbol omega-hat -> axis, theta -> angle.
Matrix3d rodrigues(const Vector3d& axis, double theta) {
    Vector3d w = axis.normalized();          // enforce ||w|| = 1
    Matrix3d K;                              // K = [w], the skew matrix
    K <<      0, -w.z(),  w.y(),
         w.z(),      0, -w.x(),
        -w.y(),  w.x(),      0;
    return Matrix3d::Identity()
         + std::sin(theta) * K
         + (1.0 - std::cos(theta)) * (K * K);
}

// Matrix log: recover (axis, angle) from R. Handles the two degeneracies.
void rotationLog(const Matrix3d& R, Vector3d& axis, double& theta) {
    double c = (R.trace() - 1.0) * 0.5;
    c = std::max(-1.0, std::min(1.0, c));    // clamp for acos safety
    theta = std::acos(c);
    if (theta < 1e-9) {                       // theta ~ 0: axis undefined
        axis = Vector3d::UnitZ();             // pick anything
    } else if (M_PI - theta < 1e-6) {         // theta ~ pi: skew part vanishes
        // recover from diagonal: w_i^2 = (R_ii + 1)/2
        Vector3d d((R(0,0)+1)*0.5, (R(1,1)+1)*0.5, (R(2,2)+1)*0.5);
        axis = d.cwiseMax(0.0).cwiseSqrt().normalized();
    } else {
        axis = Vector3d(R(2,1)-R(1,2),
                        R(0,2)-R(2,0),
                        R(1,0)-R(0,1)) / (2.0 * std::sin(theta));
    }
}

int main() {
    // Eigen idioms: axis-angle -> quaternion -> matrix, all consistent.
    AngleAxisd aa(M_PI / 2.0, Vector3d::UnitZ());  // 90 deg about z
    Quaterniond q(aa);                              // q = (cos45, 0,0, sin45)
    Matrix3d R = q.toRotationMatrix();              // same as rodrigues(z, pi/2)

    std::cout << "q.w=" << q.w() << "  R*x =\n"
              << (R * Vector3d::UnitX()).transpose() << "\n"; // -> (0,1,0)

    // Composition is matrix/quaternion product, left-to-right = applied last-first.
    Quaterniond qz(AngleAxisd(M_PI/2, Vector3d::UnitZ()));
    Quaterniond qx(AngleAxisd(M_PI/2, Vector3d::UnitX()));
    Quaterniond zx = qz * qx;   // apply x-rotation first, then z
    Quaterniond xz = qx * qz;   // different result: non-commutative
    std::cout << "qz*qx != qx*qz : "
              << (zx.angularDistance(xz) > 1e-9) << "\n";

    // Slerp: shortest-arc constant-speed interpolation. Eigen negates internally.
    Quaterniond mid = qz.slerp(0.5, qx);   // halfway orientation
    std::cout << "slerp midpoint w=" << mid.normalized().w() << "\n";

    Vector3d axis; double theta;
    rotationLog(R, axis, theta);
    std::cout << "log: theta=" << theta << " axis=" << axis.transpose() << "\n";
    return 0;
}
```

Short Python check with SciPy, useful for prototyping before the C++ port:

```python
import numpy as np
from scipy.spatial.transform import Rotation as Rot

r = Rot.from_rotvec(np.pi/2 * np.array([0, 0, 1]))   # axis-angle = rotvec
print(r.apply([1, 0, 0]))            # -> [0, 1, 0]
print(r.as_quat())                   # SciPy order is (x, y, z, w)!  -> (0,0,.707,.707)
# Non-commutativity:
rz = Rot.from_euler('z', 90, degrees=True)
rx = Rot.from_euler('x', 90, degrees=True)
print((rz * rx).apply([1,0,0]), (rx * rz).apply([1,0,0]))   # differ
```

Note the convention trap: SciPy stores quaternions as $(x,y,z,w)$ while Eigen stores $(w,x,y,z)$. Mixing them silently corrupts orientations.

## Connections

- [[03_vectors]]: rotations act on vectors; the cross product builds the skew matrix $[\hat\omega]$ and the quaternion product, and the dot product gives $\cos\Omega$ in slerp. Without vector fluency the formulas are opaque.
- [[05_transforms-se3]]: a rigid transform $T \in SE(3)$ bolts a translation $t$ onto a rotation $R$. Everything here is the $R$ block; $SE(3)$ is the next layer up and is the actual state object for ICP, IFC placement, and robot pose.
- [[06_twists-screws-exp]]: axis-angle is the exponential map on $SO(3)$. Twists and screws generalize it to $SE(3)$, where a single screw axis captures simultaneous rotation and translation. The Rodrigues derivation is the warm-up for the screw exponential.
- ICP and registration (downstream): the rotation solve inside ICP returns an $SO(3)$ matrix via SVD, which you must keep orthonormal using the tools here.

## Pitfalls and failure modes

- Frame and convention mismatch. Active vs passive, intrinsic vs extrinsic Euler, $(w,x,y,z)$ vs $(x,y,z,w)$, row-vector vs column-vector. Most "rotation bugs" are silent convention swaps between two libraries. Write the convention down at every boundary.
- Non-orthonormal matrices. Repeated float products drift $R$ off $SO(3)$, so $R^{\mathsf T}R \neq I$ and stone gets sheared. Re-orthonormalize periodically (Gram-Schmidt or, better, project via SVD: $R \leftarrow UV^{\mathsf T}$).
- Gimbal lock. Euler angles lose a degree of freedom at the singular middle angle and the back-conversion is ill-conditioned there. Never store orientation as Euler angles in a solver; convert only at I/O.
- Matrix-log degeneracies. $\theta=0$ leaves the axis undefined; $\theta=\pi$ makes $R-R^{\mathsf T}=0$ so the skew formula divides by zero. Branch on both and clamp the $\arccos$ argument to $[-1,1]$ before calling it.
- Slerp at small angles. $\sin\Omega \to 0$ blows up the weights; fall back to nlerp. And always check the dot-product sign to avoid the long way around.
- Forgetting non-commutativity. $R_1 R_2 \neq R_2 R_1$. Compose in the order operations are physically applied, and remember left-multiply applies last.
- Unnormalized quaternions. A quaternion that drifts off the unit sphere encodes a rotation plus a scale. Renormalize after every few products.

## Practice

1. By hand, build $R_x(90^\circ)$ and $R_y(90^\circ)$, compute both $R_x R_y$ and $R_y R_x$, and apply each to $(1,1,1)$. Confirm the two results differ, demonstrating non-commutativity. Inspect the numbers.
2. Implement `rodrigues(axis, theta)` and `rotationLog(R)` in C++/Eigen. Round-trip a random axis-angle through matrix and back; assert the recovered angle matches to $10^{-9}$. Test $\theta$ near $0$ and near $\pi$ explicitly.
3. Write a converter `quatToMatrix` from the 9-entry formula and compare against `Quaterniond::toRotationMatrix()` on 1000 random unit quaternions. Report the max Frobenius-norm error.
4. Implement slerp from the formula. Interpolate from identity to a $170^\circ$ rotation in 10 steps and print the angle of each intermediate (should increase linearly by $17^\circ$). Then flip the sign of $q_1$ and confirm your dot-product check still yields the short path.
5. Take a small synthetic point cloud, apply a known rotation, then recover it with the SVD rigid-alignment trick (Kabsch). Verify the recovered $R$ equals the planted one and that $\det R = +1$.
6. Construct a $ZYX$ Euler triple with $\beta = 90^\circ$ (gimbal lock), convert to a matrix, convert back, and show that many $(\alpha,\gamma)$ pairs produce the same matrix. Inspect the ambiguity numerically.

## References

- Kevin M. Lynch and Frank C. Park, *Modern Robotics: Mechanics, Planning, and Control*, Cambridge University Press, 2017. Free PDF: https://hades.mech.northwestern.edu/images/7/7f/MR.pdf . Chapter 3 (rigid-body motions): $SO(3)$, Rodrigues' formula, exponential coordinates, matrix log.
- Modern Robotics course resource, "Exponential Coordinates of Rotation (Part 2)": https://modernrobotics.northwestern.edu/nu-gm-book-resource/3-2-3-exponential-coordinates-of-rotation-part-2-of-2/
- Wikipedia, "Rodrigues' rotation formula": https://en.wikipedia.org/wiki/Rodrigues%27_rotation_formula
- Wikipedia, "Slerp": https://en.wikipedia.org/wiki/Slerp
- Ken Shoemake, "Animating Rotation with Quaternion Curves", SIGGRAPH 1985; tutorial PDF: https://www.ljll.fr/~frey/papers/scientific%20visualisation/Shoemake%20K.,%20Quaternions.pdf . Origin of quaternion slerp.
- EuclideanSpace, "Quaternion to Matrix": https://www.euclideanspace.com/maths/geometry/rotations/conversions/quaternionToMatrix/ . Verified 9-entry conversion.
- Eigen Geometry tutorial (AngleAxis, Quaternion, slerp, toRotationMatrix): https://libeigen.gitlab.io/eigen/docs-nightly/group__TutorialGeometry.html
- OpenCV calib3d documentation (Rodrigues, projection conventions): https://docs.opencv.org/4.x/d9/d0c/group__calib3d.html
- OpenCASCADE reference manual (gp_Trsf, gp_Quaternion for CAD placement): https://dev.opencascade.org/doc/refman/html/
- CADFit (CAD-program synthesis), arXiv:2605.01171 and https://github.com/ghadinehme/CADFit . Predicted CAD ops include axis-angle rotations that must match the kernel convention.

## Flashcard seeds

- Q: Derive the 2D rotation matrix $R(\theta)$ by tracking where the basis vectors land. <br> <b>Formula:</b> a point at angle $\theta$ on the unit circle is $(\cos\theta,\sin\theta)$; target $R(\theta)=\begin{bmatrix}\cos\theta & -\sin\theta\\ \sin\theta & \cos\theta\end{bmatrix}$ <br> <i>Hint:</i> first write down the image of $\hat e_1$, then note $\hat e_2$ leads it by $90^\circ$. :: A: A point at angle $\theta$ on the unit circle is $(\cos\theta,\sin\theta)$, so $\hat e_1\mapsto(\cos\theta,\sin\theta)$ <br> $\hat e_2$ leads $\hat e_1$ by $90^\circ$, so $\hat e_2\mapsto(-\sin\theta,\cos\theta)$ <br> Stack the images as columns: $\boxed{R(\theta)=\begin{bmatrix}\cos\theta & -\sin\theta\\ \sin\theta & \cos\theta\end{bmatrix}}$.
- Q: Prove that composing two 2D rotations adds their angles: $R(\alpha)R(\beta)=R(\alpha+\beta)$. <br> <b>Formula:</b> $\cos(\alpha+\beta)=\cos\alpha\cos\beta-\sin\alpha\sin\beta$ and $\sin(\alpha+\beta)=\sin\alpha\cos\beta+\cos\alpha\sin\beta$ <br> <i>Hint:</i> multiply the matrices; compute the top-left entry first and match the cosine identity, then the bottom-left and match the sine identity. :: A: Top-left of the product $=\cos\alpha\cos\beta-\sin\alpha\sin\beta=\cos(\alpha+\beta)$ <br> Bottom-left $=\sin\alpha\cos\beta+\cos\alpha\sin\beta=\sin(\alpha+\beta)$ <br> The remaining entries follow by the same identities, giving the matrix $R(\alpha+\beta)$ <br> $\boxed{R(\alpha)R(\beta)=R(\alpha+\beta)}$, so composing 2D rotations adds angles.
- Q: From the $SO(3)$ orthonormality condition, derive that inverting a rotation is just transposing it. <br> <b>Formula:</b> $R^{\mathsf T}R=I$; target $R^{-1}=R^{\mathsf T}$ <br> <i>Hint:</i> read $R^{\mathsf T}R=I$ as "$R^{\mathsf T}$ is a left inverse" and use that for square matrices a left inverse is the inverse. :: A: $R^{\mathsf T}R=I$ says $R^{\mathsf T}$ is a left inverse of $R$ <br> For square matrices a left inverse is the inverse <br> $\boxed{R^{-1}=R^{\mathsf T}}$, so inverting a rotation is just transposing it.
- Q: Prove the skew-matrix identity for a unit axis, $[\hat\omega]^2=\hat\omega\hat\omega^{\mathsf T}-I$. <br> <b>Formula:</b> $[\hat\omega]v=\hat\omega\times v$ and BAC-CAB: $a\times(b\times c)=b(a\cdot c)-c(a\cdot b)$ <br> <i>Hint:</i> apply $[\hat\omega]^2$ to a test vector $v$ so it becomes $\hat\omega\times(\hat\omega\times v)$, then expand with BAC-CAB and use $\hat\omega\cdot\hat\omega=1$. :: A: $[\hat\omega]^2 v=\hat\omega\times(\hat\omega\times v)$ <br> BAC-CAB: $\hat\omega\times(\hat\omega\times v)=\hat\omega(\hat\omega\cdot v)-v(\hat\omega\cdot\hat\omega)$ <br> With $\hat\omega\cdot\hat\omega=1$ this is $(\hat\omega\hat\omega^{\mathsf T}-I)v$ for all $v$ <br> $\boxed{[\hat\omega]^2=\hat\omega\hat\omega^{\mathsf T}-I}$.
- Q: Derive Rodrigues' vector form by splitting $v$ into parallel and perpendicular parts. <br> <b>Formula:</b> $v_\parallel=\hat\omega(\hat\omega\cdot v)$, $v_\perp=v-v_\parallel$; target $v_{\text{rot}}=v\cos\theta+(\hat\omega\times v)\sin\theta+\hat\omega(\hat\omega\cdot v)(1-\cos\theta)$ <br> <i>Hint:</i> the parallel part is unchanged; rotate $v_\perp$ in its plane using $\cos\theta$ along itself and $\sin\theta$ along $\hat\omega\times v$, then substitute $v_\perp=v-\hat\omega(\hat\omega\cdot v)$. :: A: $v_\parallel$ is on the axis and is unchanged <br> $v_\perp$ rotates in its plane: $\cos\theta\,v_\perp$ along itself plus $\sin\theta\,(\hat\omega\times v)$ perpendicular <br> Add and use $v_\perp=v-\hat\omega(\hat\omega\cdot v)$: $\boxed{v_{\text{rot}}=v\cos\theta+(\hat\omega\times v)\sin\theta+\hat\omega(\hat\omega\cdot v)(1-\cos\theta)}$.
- Q: Convert Rodrigues' vector form into the matrix form $R=I+\sin\theta[\hat\omega]+(1-\cos\theta)[\hat\omega]^2$. <br> <b>Formula:</b> $v_{\text{rot}}=v\cos\theta+(\hat\omega\times v)\sin\theta+\hat\omega(\hat\omega\cdot v)(1-\cos\theta)$ with $\hat\omega\times v=[\hat\omega]v$ and $\hat\omega\hat\omega^{\mathsf T}=[\hat\omega]^2+I$ <br> <i>Hint:</i> rewrite each term as (matrix)$\,v$, then collect the $I$ terms since $\cos\theta+(1-\cos\theta)=1$. :: A: Write each term as a matrix times $v$: $v\cos\theta=\cos\theta\,I\,v$; $(\hat\omega\times v)=[\hat\omega]v$; $\hat\omega(\hat\omega\cdot v)=\hat\omega\hat\omega^{\mathsf T}v=([\hat\omega]^2+I)v$ <br> So $R=\cos\theta\,I+\sin\theta[\hat\omega]+(1-\cos\theta)([\hat\omega]^2+I)$ <br> Collect the $I$ terms ($\cos\theta+1-\cos\theta=1$): $\boxed{R=I+\sin\theta[\hat\omega]+(1-\cos\theta)[\hat\omega]^2}$.
- Q: Derive the angle-recovery formula $\theta=\arccos((\operatorname{tr}R-1)/2)$ from Rodrigues' $R$. <br> <b>Formula:</b> $R=I+\sin\theta[\hat\omega]+(1-\cos\theta)[\hat\omega]^2$ with $\operatorname{tr}I=3$, $\operatorname{tr}[\hat\omega]=0$, $\operatorname{tr}[\hat\omega]^2=-2$ <br> <i>Hint:</i> take the trace term by term, simplify to $1+2\cos\theta$, then solve for $\theta$. :: A: $\operatorname{tr}R=\operatorname{tr}I+\sin\theta\,\operatorname{tr}[\hat\omega]+(1-\cos\theta)\operatorname{tr}[\hat\omega]^2$ <br> $=3+0-2(1-\cos\theta)=1+2\cos\theta$ <br> Solve for $\theta$: $\boxed{\theta=\arccos\!\big((\operatorname{tr}R-1)/2\big)}$.
- Q: Derive the skew (axis) part of the matrix log, $[\hat\omega]=\frac{1}{2\sin\theta}(R-R^{\mathsf T})$. <br> <b>Formula:</b> $R=I+\sin\theta[\hat\omega]+(1-\cos\theta)[\hat\omega]^2$ with $[\hat\omega]^{\mathsf T}=-[\hat\omega]$ and $[\hat\omega]^2$ symmetric <br> <i>Hint:</i> form $R-R^{\mathsf T}$; the $I$ and symmetric $[\hat\omega]^2$ terms cancel, leaving only the antisymmetric term. :: A: $R-R^{\mathsf T}=\sin\theta([\hat\omega]-[\hat\omega]^{\mathsf T})+(1-\cos\theta)([\hat\omega]^2-([\hat\omega]^2)^{\mathsf T})$ <br> Symmetric parts cancel; $[\hat\omega]^{\mathsf T}=-[\hat\omega]$ gives $R-R^{\mathsf T}=2\sin\theta[\hat\omega]$ <br> $\boxed{[\hat\omega]=\dfrac{1}{2\sin\theta}(R-R^{\mathsf T})}$.
- Q: At $\theta=\pi$ the skew log divides by zero; derive the diagonal axis-recovery formula $\hat\omega_i^2=(R_{ii}+1)/2$. <br> <b>Formula:</b> at $\theta=\pi$, $R=I+2[\hat\omega]^2$, and $[\hat\omega]^2=\hat\omega\hat\omega^{\mathsf T}-I$ <br> <i>Hint:</i> substitute the second into the first to get $R=2\hat\omega\hat\omega^{\mathsf T}-I$, then read off a diagonal entry $R_{ii}$. :: A: At $\theta=\pi$: $\sin\theta=0$, $1-\cos\theta=2$, so $R=I+2[\hat\omega]^2$ <br> Substitute $[\hat\omega]^2=\hat\omega\hat\omega^{\mathsf T}-I$: $R=2\hat\omega\hat\omega^{\mathsf T}-I$ <br> Diagonal entry $R_{ii}=2\hat\omega_i^2-1$, so $\boxed{\hat\omega_i^2=(R_{ii}+1)/2}$.
- Q: Derive the half-angle in the quaternion: for axis $\hat z$ and angle $\theta$, find $\phi$ in $q=(\cos\phi,0,0,\sin\phi)$ so $R(q)=R_z(\theta)$. <br> <b>Formula:</b> $R(q)_{11}=1-2z^2$ with $z=\sin\phi$, double-angle $1-2\sin^2\phi=\cos2\phi$, and $R_z(\theta)_{11}=\cos\theta$ <br> <i>Hint:</i> set the two top-left entries equal so $\cos2\phi=\cos\theta$, then solve for $\phi$. :: A: With $w=\cos\phi$, $z=\sin\phi$, $R(q)_{11}=1-2z^2=1-2\sin^2\phi=\cos 2\phi$ <br> This must equal $\cos\theta$, so $2\phi=\theta$ <br> $\boxed{\phi=\theta/2,\quad q=(\cos\tfrac{\theta}{2},\hat\omega\sin\tfrac{\theta}{2})}$.
- Q: From the half-angle form, derive why $q$ and $-q$ represent the same rotation (the double cover). <br> <b>Formula:</b> $q=(\cos\tfrac{\theta}{2},\hat\omega\sin\tfrac{\theta}{2})$, and $\cos\tfrac{\theta+2\pi}{2}=-\cos\tfrac{\theta}{2}$, $\sin\tfrac{\theta+2\pi}{2}=-\sin\tfrac{\theta}{2}$ <br> <i>Hint:</i> add $2\pi$ to $\theta$ (same physical rotation) and watch what happens to the sign of both components. :: A: Replacing $\theta\to\theta+2\pi$ leaves the physical rotation identical <br> But $\cos\tfrac{\theta+2\pi}{2}=-\cos\tfrac{\theta}{2}$ and $\sin\tfrac{\theta+2\pi}{2}=-\sin\tfrac{\theta}{2}$ <br> So the quaternion flips sign while the rotation is unchanged: $\boxed{q\text{ and }-q\text{ map to the same }R\text{ (double cover)}}$.
- Q: A scanned 3 m slab is mis-registered by $\theta=0.01$ rad about $\hat z$. By how much does its far corner at $(3,0,0)$ m move sideways? <br> <b>Formula:</b> $R_z(\theta)(x,0,0)=(x\cos\theta,\,x\sin\theta,\,0)$ <br> <i>Hint:</i> the sideways shift is the $y$-component $3\sin\theta$; compute $\sin 0.01$ first (radians, so $\approx0.01$). :: A: $R_z(0.01)(3,0,0)=(3\cos 0.01,\,3\sin 0.01,\,0)$ <br> $3\sin 0.01\approx3(0.0099998)\approx0.0300$ m lateral; $x$ drops by $3(1-\cos0.01)\approx1.5\times10^{-4}$ m <br> $\boxed{\approx30\text{ mm sideways}}$, showing tiny angle errors ruin a 3 m placement.
- Q: Rotate $v=(1,0,0)$ by $\theta=120^\circ$ about $\hat\omega=(1,1,1)/\sqrt3$ using Rodrigues' vector form. <br> <b>Formula:</b> $v_{\text{rot}}=v\cos\theta+(\hat\omega\times v)\sin\theta+\hat\omega(\hat\omega\cdot v)(1-\cos\theta)$ <br> <i>Hint:</i> compute the three scalars first ($\cos120^\circ=-\tfrac12$, $\sin120^\circ=\tfrac{\sqrt3}{2}$, $\hat\omega\cdot v=1/\sqrt3$) and the vector $\hat\omega\times v$, then plug in. :: A: $\cos120^\circ=-\tfrac12$, $\sin120^\circ=\tfrac{\sqrt3}{2}$, $\hat\omega\cdot v=1/\sqrt3$, $\hat\omega\times v=(0,\tfrac{1}{\sqrt3},-\tfrac{1}{\sqrt3})$ <br> $v_{\text{rot}}=-\tfrac12(1,0,0)+\tfrac{\sqrt3}{2}(0,\tfrac1{\sqrt3},-\tfrac1{\sqrt3})+\tfrac{1}{\sqrt3}\tfrac{3}{2}(\tfrac1{\sqrt3},\tfrac1{\sqrt3},\tfrac1{\sqrt3})$ <br> $\boxed{v_{\text{rot}}=(0,1,0)}$, the axis-cyclic permutation $x\to y$.
- Q: Build the unit quaternion for a $60^\circ$ rotation about $\hat\omega=(0,0,1)$ and state its scalar part $w$. <br> <b>Formula:</b> $q=(\cos\tfrac{\theta}{2},\,\hat\omega\sin\tfrac{\theta}{2})$ <br> <i>Hint:</i> use the half-angle $\theta/2=30^\circ$ first, then read off $w=\cos30^\circ$ (not $\cos60^\circ$). :: A: $q=(\cos30^\circ,\,0,0,\sin30^\circ)$ <br> $=(\tfrac{\sqrt3}{2},0,0,\tfrac12)$ <br> $\boxed{w=\cos30^\circ\approx0.866}$ (half-angle, not $\cos60^\circ$).
- Q: Closing test: you got $\theta=\arccos((\operatorname{tr}R-1)/2)$ from an $R$ with $\operatorname{tr}R=2$. Find $\theta$ and confirm it against a $60^\circ$ rotation's trace. <br> <b>Formula:</b> forward check $\operatorname{tr}R=1+2\cos\theta$ <br> <i>Hint:</i> a pass means the recovered $\theta$ plugged back into $1+2\cos\theta$ returns the original trace $2$. :: A: $\theta=\arccos((2-1)/2)=\arccos(0.5)=60^\circ$ <br> Check forward: $\operatorname{tr}R=1+2\cos60^\circ=1+2(0.5)=2$ <br> $\boxed{\theta=60^\circ,\ \text{trace round-trips to }2}$, loop closed.
- Q: Closing test: slerp from identity to a $90^\circ$ $z$-rotation; predict the rotation angle at $t=1/3$. <br> <b>Formula:</b> $\cos\Omega=q_0\cdot q_1$, and constant angular speed means the angle scales linearly with $t$ <br> <i>Hint:</i> a pass means the angle is exactly $t$ times the total; the quaternions sit $\Omega=45^\circ$ apart (half-angle of $90^\circ$). :: A: Quaternions $90^\circ$ apart means $\Omega=\arccos(q_0\cdot q_1)=45^\circ$ (half-angle) <br> Constant speed: at $t=1/3$ the rotation angle is $\tfrac13(90^\circ)=30^\circ$ <br> Check: $\cos(\theta/2)=\cos15^\circ$ matches the slerp scalar part, so $\boxed{\theta=30^\circ}$.
- Q: State the two conditions that define a matrix as a member of $SO(3)$. <br> <b>Formula:</b> orthonormality and a determinant condition <br> <i>Hint:</i> one condition makes the columns orthonormal, the other rules out reflections. :: A: $\boxed{R^{\mathsf T}R = I\ \text{and}\ \det R = +1}$ (orthonormal columns, no reflection).
- Q: Which rotation representation should be a solver's internal state, and why? <br> <i>Hint:</i> name the one with no gimbal lock that is cheap to compose and renormalize. :: A: $\boxed{\text{Unit quaternions}}$: no gimbal lock, cheap to compose and renormalize; convert to Euler only at I/O.
- Q: When is the matrix-log skew formula $\frac{1}{2\sin\theta}(R-R^{\mathsf T})$ invalid, and what is the fix at $\theta=\pi$? <br> <b>Formula:</b> diagonal fallback $\hat\omega_i^2=(R_{ii}+1)/2$ <br> <i>Hint:</i> spot the two angles where $\sin\theta=0$; at $\theta=\pi$ recover the axis from the diagonal instead. :: A: At $\theta=0$ (axis undefined, pick any) and at $\theta=\pi$ ($\sin\theta=0$, divide-by-zero) <br> $\boxed{\theta=\pi:\ \text{recover axis from the diagonal },\hat\omega_i^2=(R_{ii}+1)/2}$.
- Q: Why must you check the sign of $q_0\cdot q_1$ before slerp, and what do you do if it is negative? <br> <b>Formula:</b> the 4D dot $q_0\cdot q_1=\cos\Omega$; double cover means $q_1\equiv-q_1$ <br> <i>Hint:</i> a negative dot means $\Omega>90^\circ$, i.e. the long arc; negating one quaternion fixes it. :: A: The double cover means $q_1$ and $-q_1$ are the same orientation <br> If $q_0\cdot q_1<0$ the formula traces the long way around <br> $\boxed{\text{Negate one quaternion when }q_0\cdot q_1<0\text{ to take the short arc}}$.
- Q: How do you re-orthonormalize a drifted rotation matrix via SVD, and how do you build/slerp rotations in Eigen? <br> <b>Formula:</b> SVD $R=U\Sigma V^{\mathsf T}$; nearest rotation is $R\leftarrow UV^{\mathsf T}$ <br> <i>Hint:</i> drop the singular-value matrix $\Sigma$ and keep $UV^{\mathsf T}$; in Eigen build with `AngleAxisd(angle, axis)` and interpolate with `q1.slerp(t, q2)`. :: A: SVD: compute $R=U\Sigma V^{\mathsf T}$ and drop $\Sigma$, giving $\boxed{R\leftarrow UV^{\mathsf T}}$ (nearest $SO(3)$) <br> Eigen: build with `Eigen::AngleAxisd(angle, axis)` (axis normalized), interpolate with `q1.slerp(t, q2)`.
- CLOZE: Rotation composition is {{c1::non-commutative}}, so $R_1 R_2 \neq R_2 R_1$ in general.
- CLOZE: The skew-symmetric matrix satisfies $[\hat\omega]\,v = {{c1::\hat\omega \times v}}$, and for a unit axis $[\hat\omega]^2 = {{c2::\hat\omega\hat\omega^{\mathsf T}-I}}$. <br> <i>hint:</i> the first is the cross product; expand the second with BAC-CAB.
- CLOZE: Eigen stores quaternions in order $(w,x,y,z)$ while SciPy uses {{c1::(x,y,z,w)}}.
- CLOZE: Axis-angle is also called {{c1::exponential coordinates}} because $R = e^{[\hat\omega]\theta}$.
- CLOZE: Taking the trace of Rodrigues' formula gives $\operatorname{tr}R = {{c1::1+2\cos\theta}}$, which inverts to the angle-recovery formula.
