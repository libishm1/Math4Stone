# Twists, screw axes, and the matrix exponential se(3)->SE(3)

Prereqs: [[03_trig-rotations]], [[11_transforms-se3]] => this => [[13_forward-kinematics]], [[14_manipulator-jacobian]]

## Why it matters for the stone scan->CAD->IFC / CADFit workflow

Every rigid motion in your pipeline is one screw motion: a rotation about an axis plus a translation along it. When you register a scan to a CAD template, the ICP output is a transform $T \in SE(3)$, and the cleanest way to iterate, interpolate, or average those transforms is in the $se(3)$ tangent space via the log map. When CADFit fits an extrusion or revolution, the revolution itself is a one-parameter screw and the placement of each fitted primitive is a $T$ you must compose and invert without drift. When you serialize to IFC, every `IfcLocalPlacement` is a rigid transform, and re-expressing a feature frame inside a parent frame is exactly an adjoint operation $\mathrm{Ad}_T$. The exponential map gives you a minimal, singularity-free 6-number coordinate (a twist) for any pose, which is the right state to optimize over when you align messy stone geometry. This same machinery is the foundation for the Product of Exponentials forward kinematics in the next two topics.

## Intuition

Chasles' theorem says any rigid-body displacement is a screw motion: turn by some angle about a fixed axis in space while sliding along that same axis. Picture tightening a bolt. The bolt rotates and advances together, locked by the thread pitch. That single axis plus the rotation angle and slide distance fully describe where the body went.

A twist is the velocity version of this. It packs the angular velocity $\omega$ and a linear velocity $v$ into one 6-vector. The matrix exponential is the "integrator": feed it a twist and an amount $\theta$, and it produces the finite pose $T$ you reach by riding that screw for parameter $\theta$. The log map is the inverse: hand it a pose, and it tells you which screw axis and how far to travel to get there. Rotations and translations do not commute, so you cannot just add them; the exponential handles the coupling correctly.

## The math (first principles)

### Conventions

- Right-handed frames, column vectors, active transforms (the matrix moves points, not axes).
- $SO(3)$ = $3\times3$ rotation matrices, $R^\top R = I$, $\det R = +1$.
- $SE(3)$ = $4\times4$ homogeneous transforms $T = \begin{bmatrix} R & p \\ 0 & 1 \end{bmatrix}$, $R \in SO(3)$, $p \in \mathbb{R}^3$.
- Angles in radians. Lengths in your working unit (meters at site scale, millimeters at shop scale; keep one consistent unit per solve).
- $\omega$ angular velocity (rad/s or unit axis), $v$ linear velocity of the body-fixed point currently at the origin.

### Skew / hat operator

The hat operator turns a 3-vector into a skew-symmetric matrix so the cross product becomes matrix multiply:

$$[\omega] = \hat\omega = \begin{bmatrix} 0 & -\omega_3 & \omega_2 \\ \omega_3 & 0 & -\omega_1 \\ -\omega_2 & \omega_1 & 0 \end{bmatrix}, \qquad [\omega]\,x = \omega \times x .$$

This $3\times3$ matrix is an element of $so(3)$, the Lie algebra of $SO(3)$.

### Twist and its matrix form in se(3)

A twist (spatial velocity) is the 6-vector $\mathcal{V} = \begin{bmatrix} \omega \\ v \end{bmatrix}$. Its matrix form lives in $se(3)$:

$$[\mathcal{V}] = \begin{bmatrix} [\omega] & v \\ 0 & 0 \end{bmatrix} \in se(3).$$

The relation $\dot T\, T^{-1} = [\mathcal{V}_s]$ defines the **spatial** twist (fixed-frame), and $T^{-1}\dot T = [\mathcal{V}_b]$ defines the **body** twist. Same motion, two frames.

### Screw axis

A normalized screw axis is a twist scaled so the "speed" is unit. If there is rotation, scale so $\|\omega\| = 1$:

$$\mathcal{S} = \begin{bmatrix} \omega \\ v \end{bmatrix}, \quad \|\omega\| = 1, \qquad v = -\omega \times q + h\,\omega,$$

where $q$ is any point on the axis and $h$ is the pitch (translation per radian). Pure rotation has $h = 0$, so $v = -\omega \times q$. Pure translation has $\omega = 0$ and $\|v\| = 1$. With $\mathcal{S}$ normalized, the scalar $\theta$ is the angle of rotation (or the distance, for the pure-translation case). The pair $\mathcal{S}\theta$ are the **exponential coordinates** of the displacement.

### Rodrigues formula (so(3) -> SO(3))

For a unit axis $\omega$ and angle $\theta$, the matrix exponential of $[\omega]\theta$ is the rotation

$$R = e^{[\omega]\theta} = I + \sin\theta\,[\omega] + (1-\cos\theta)\,[\omega]^2 .$$

*Derivation sketch:* $[\omega]^3 = -[\omega]$ for a unit axis, so the power series $\sum \frac{([\omega]\theta)^k}{k!}$ collapses. Odd terms sum to $\sin\theta\,[\omega]$, even terms (past identity) sum to $(1-\cos\theta)[\omega]^2$. This is the closed form used everywhere for axis-angle to matrix.

### Matrix exponential (se(3) -> SE(3))

For a screw axis $\mathcal{S} = (\omega, v)$ with $\|\omega\| = 1$, riding the screw for $\theta$ gives

$$e^{[\mathcal{S}]\theta} = \begin{bmatrix} e^{[\omega]\theta} & G(\theta)\,v \\ 0 & 1 \end{bmatrix}, \qquad G(\theta) = I\theta + (1-\cos\theta)\,[\omega] + (\theta - \sin\theta)\,[\omega]^2 .$$

$G(\theta)$ is the integral $\int_0^\theta e^{[\omega]s}\,ds$. It couples rotation into the translation, which is why you cannot simply use $p = v\theta$.

For the pure-translation case $\omega = 0$, $\|v\| = 1$:

$$e^{[\mathcal{S}]\theta} = \begin{bmatrix} I & v\theta \\ 0 & 1 \end{bmatrix}.$$

### Log map (SE(3) -> se(3))

Given $T = (R, p)$, recover the screw axis and $\theta$.

1. If $R = I$: then $\omega = 0$, $\theta = \|p\|$, $v = p/\|p\|$ (pure translation). If also $p = 0$, the motion is identity ($\theta = 0$).
2. Otherwise recover the rotation first:

$$\theta = \cos^{-1}\!\left(\tfrac{\mathrm{tr}\,R - 1}{2}\right), \qquad [\omega] = \frac{1}{2\sin\theta}\,(R - R^\top).$$

3. Then recover $v = G^{-1}(\theta)\,p$ with the closed-form inverse

$$G^{-1}(\theta) = \tfrac{1}{\theta} I - \tfrac{1}{2}[\omega] + \left(\tfrac{1}{\theta} - \tfrac{1}{2}\cot\tfrac{\theta}{2}\right)[\omega]^2 .$$

The $\theta = \pi$ case needs the standard special handling (pick $\omega$ from a positive-diagonal column of $R$) because $\sin\theta \to 0$.

### Adjoint: moving twists between frames

A twist is a $6$-vector, so a $4\times4$ transform cannot multiply it directly. The $6\times6$ adjoint of $T = (R, p)$ does:

$$\mathrm{Ad}_T = \begin{bmatrix} R & 0 \\ [p]\,R & R \end{bmatrix}, \qquad \mathcal{V}_a = \mathrm{Ad}_{T_{ab}}\,\mathcal{V}_b .$$

It re-expresses the same physical twist in another frame. It also relates the two screw forms: $[\mathcal{V}_s] = T[\mathcal{V}_b]T^{-1}$ is the matrix-conjugation version, and $\mathcal{V}_s = \mathrm{Ad}_T \mathcal{V}_b$ is the 6-vector version. The adjoint is what builds each column of the Product-of-Exponentials Jacobian in the next topic.

## Worked example

Take a **pure rotation** by $\theta = \pi/2$ about the $z$-axis, where the axis passes through the point $q = (2, 0, 0)$. Pitch $h = 0$.

**Build the screw axis.** $\omega = (0,0,1)$, already unit. Then

$$v = -\omega \times q = -\,(0,0,1)\times(2,0,0) = -(0,2,0) = (0,-2,0).$$

So $\mathcal{S} = (0,0,1,\;0,-2,0)$.

**Rotation part (Rodrigues).** With $\theta = \pi/2$, $\sin\theta = 1$, $\cos\theta = 0$:

$$R = e^{[\omega]\theta} = \begin{bmatrix} 0 & -1 & 0 \\ 1 & 0 & 0 \\ 0 & 0 & 1 \end{bmatrix}.$$

**Translation part.** Compute $G(\theta) = I\theta + (1-\cos\theta)[\omega] + (\theta - \sin\theta)[\omega]^2$. With $\theta = \pi/2 \approx 1.5708$, $(1-\cos\theta) = 1$, $(\theta - \sin\theta) \approx 0.5708$. For $\omega = (0,0,1)$, $[\omega]^2 = \mathrm{diag}(-1,-1,0)$. Then

$$p = G(\theta)\,v = \begin{bmatrix} 2 \\ -2 \\ 0 \end{bmatrix}.$$

**Check by hand.** The point at the origin sits at distance $2$ from the axis. Rotating it $+90^\circ$ about the center $(2,0)$ in the $xy$-plane: displacement from center is $(-2, 0)$, rotated $90^\circ$ becomes $(0, -2)$, added back to center gives $(2, -2)$. So the origin maps to $(2,-2,0)$. The full transform is

$$T = e^{[\mathcal{S}]\theta} = \begin{bmatrix} 0 & -1 & 0 & 2 \\ 1 & 0 & 0 & -2 \\ 0 & 0 & 1 & 0 \\ 0 & 0 & 0 & 1 \end{bmatrix}, \qquad T\,(0,0,0,1)^\top = (2,-2,0,1)^\top. \checkmark$$

**Log map round-trip.** From $T$: $\mathrm{tr}\,R = 1$, so $\theta = \cos^{-1}(0) = \pi/2$. $\frac{1}{2\sin\theta}(R - R^\top)$ gives $[\omega]$ with $\omega = (0,0,1)$. Then $v = G^{-1}(\theta)p = (0,-2,0)$. We recover the original screw exactly.

## Code

C++ with Eigen. These functions map one-to-one to the formulas above and are the kernel you will reuse in forward kinematics.

```cpp
// twist_exp.hpp  -- se(3) <-> SE(3), verified against Modern Robotics ch.3
#pragma once
#include <Eigen/Dense>
#include <cmath>

namespace twist {

using Mat3 = Eigen::Matrix3d;
using Mat4 = Eigen::Matrix4d;
using Mat6 = Eigen::Matrix<double, 6, 6>;
using Vec3 = Eigen::Vector3d;
using Vec6 = Eigen::Matrix<double, 6, 1>;  // S = [omega; v]

constexpr double kEps = 1e-9;  // scale this to your model size if working in mm vs m

// [omega] : the hat / skew operator, so that hat(w)*x == w.cross(x)
inline Mat3 hat(const Vec3& w) {
    Mat3 W;
    W <<    0.0, -w.z(),  w.y(),
         w.z(),    0.0, -w.x(),
        -w.y(),  w.x(),    0.0;
    return W;
}

// Rodrigues: e^{[omega] theta}, omega need not be unit (theta absorbed into ||omega*theta||)
inline Mat3 expSO3(const Vec3& w, double theta) {
    const double n = w.norm();
    if (n < kEps) return Mat3::Identity();          // no rotation
    const Vec3 a = w / n;                            // unit axis
    const Mat3 W = hat(a);
    const double t = n * theta;                      // total angle
    return Mat3::Identity() + std::sin(t) * W + (1.0 - std::cos(t)) * (W * W);
}

// Matrix exponential e^{[S] theta}, S = [omega; v]. Handles both screw and pure-translation.
inline Mat4 expSE3(const Vec6& S, double theta) {
    const Vec3 w = S.head<3>();
    const Vec3 v = S.tail<3>();
    Mat4 T = Mat4::Identity();

    if (w.norm() < kEps) {                           // pure translation: omega = 0
        T.block<3,1>(0,3) = v * theta;
        return T;
    }
    // screw with rotation (assumes ||omega|| = 1 for the canonical screw form)
    const Mat3 W = hat(w);
    // G(theta) = I*theta + (1-cos)[w] + (theta - sin)[w]^2
    const Mat3 G = Mat3::Identity() * theta
                 + (1.0 - std::cos(theta)) * W
                 + (theta - std::sin(theta)) * (W * W);
    T.block<3,3>(0,0) = expSO3(w, theta);
    T.block<3,1>(0,3) = G * v;
    return T;
}

// Log map SE(3) -> (S, theta). Returns S = [omega; v] with theta packed out.
inline void logSE3(const Mat4& T, Vec6& S, double& theta) {
    const Mat3 R = T.block<3,3>(0,0);
    const Vec3 p = T.block<3,1>(0,3);

    if ((R - Mat3::Identity()).norm() < kEps) {      // pure translation
        S.head<3>().setZero();
        theta = p.norm();
        S.tail<3>() = (theta < kEps) ? Vec3::Zero() : (p / theta).eval();
        if (theta < kEps) theta = 0.0;
        return;
    }
    const double c = std::clamp((R.trace() - 1.0) * 0.5, -1.0, 1.0);
    theta = std::acos(c);                            // tr R = 1 + 2 cos theta
    const Mat3 W = (R - R.transpose()) / (2.0 * std::sin(theta));  // [omega]
    const Vec3 w(W(2,1), W(0,2), W(1,0));
    // G^{-1}(theta) = (1/theta) I - (1/2)[w] + (1/theta - (1/2)cot(theta/2)) [w]^2
    const Mat3 Ginv = Mat3::Identity() / theta
                    - 0.5 * W
                    + (1.0 / theta - 0.5 / std::tan(theta * 0.5)) * (W * W);
    S.head<3>() = w;
    S.tail<3>() = Ginv * p;                           // v
}

// Adjoint Ad_T : 6x6 that re-expresses a twist in another frame. V_a = Ad_{T_ab} V_b
inline Mat6 adjoint(const Mat4& T) {
    const Mat3 R = T.block<3,3>(0,0);
    const Vec3 p = T.block<3,1>(0,3);
    Mat6 Ad = Mat6::Zero();
    Ad.block<3,3>(0,0) = R;
    Ad.block<3,3>(3,3) = R;
    Ad.block<3,3>(3,0) = hat(p) * R;                 // [p] R
    return Ad;
}

} // namespace twist
```

Minimal Python check (NumPy), useful for prototyping the same numbers before porting:

```python
import numpy as np
def hat(w): return np.array([[0,-w[2],w[1]],[w[2],0,-w[0]],[-w[1],w[0],0]])
def exp_se3(S, th):
    w, v = S[:3], S[3:]
    if np.linalg.norm(w) < 1e-9:
        T = np.eye(4); T[:3,3] = v*th; return T
    W = hat(w)
    R = np.eye(3) + np.sin(th)*W + (1-np.cos(th))*(W@W)
    G = np.eye(3)*th + (1-np.cos(th))*W + (th-np.sin(th))*(W@W)
    T = np.eye(4); T[:3,:3] = R; T[:3,3] = G@v; return T
# pure rotation pi/2 about z through (2,0,0): origin -> (2,-2,0)
print(exp_se3(np.array([0,0,1,0,-2,0.0]), np.pi/2) @ np.array([0,0,0,1.0]))  # [2,-2,0,1]
```

## Connections

- [[03_trig-rotations]]: Rodrigues' formula is the bridge from axis-angle (sin/cos of $\theta$) to a rotation matrix; the whole exp map is built on it.
- [[11_transforms-se3]]: $SE(3)$ supplies the $4\times4$ group elements; this topic gives the minimal tangent-space coordinates (twists) and the maps between them.
- [[13_forward-kinematics]]: Product of Exponentials writes $T(\theta) = e^{[\mathcal{S}_1]\theta_1}\cdots e^{[\mathcal{S}_n]\theta_n} M$, one exponential per joint, directly from this kernel.
- [[14_manipulator-jacobian]]: each Jacobian column is a screw axis pushed forward by $\mathrm{Ad}$ of the preceding exponentials, defined here.
- ICP / registration topics: the log map gives the residual and update step when you optimize a pose over $se(3)$, and lets you average and interpolate transforms without leaving the manifold.

## Pitfalls and failure modes

- **Frame consistency.** Spatial twist $\mathcal{V}_s$ and body twist $\mathcal{V}_b$ are different 6-vectors for the same motion; mixing them silently corrupts every downstream pose. Subscript-cancel with $\mathrm{Ad}$ and label frames in code.
- **Small-angle division.** $G^{-1}$ divides by $\theta$ and by $\sin\theta$; near $\theta = 0$ use the pure-translation branch (or a Taylor expansion), never the general formula.
- **$\theta = \pi$ singularity.** At $\theta = \pi$, $\sin\theta = 0$ so $(R - R^\top)/(2\sin\theta)$ blows up. Use the special-case axis extraction from $R$'s diagonal.
- **Unnormalized screw axis.** The closed-form $G(\theta)$ assumes $\|\omega\| = 1$. If $\|\omega\| \neq 1$, either normalize and fold the magnitude into $\theta$, or use the general $\|\omega\|$-aware form. Do not silently pass a non-unit axis into the unit-axis formula.
- **Non-commutativity.** $e^{[\mathcal{S}_1]\theta_1}e^{[\mathcal{S}_2]\theta_2} \neq e^{[\mathcal{S}_1]\theta_1 + [\mathcal{S}_2]\theta_2}$ in general. You cannot add twists to compose motions; you compose by multiplying the resulting $T$ matrices.
- **Rotation drift.** Repeated matrix products accumulate non-orthogonality in $R$. Re-orthonormalize periodically (SVD or Gram-Schmidt) or re-derive $R$ from a clean exp.
- **Units and tolerance.** A fixed $10^{-9}$ epsilon that is fine in meters is far too tight at quarry scale or too loose in millimeters. Scale `kEps` to model size, consistent with the per-application tolerance budget.
- **Outliers in registration.** A few bad correspondences can dominate the log-space residual; pair the exp/log machinery with robust loss or correspondence rejection upstream.

## Practice

1. **Hat round-trip.** Implement `hat` and its inverse `vee`. Verify $\mathrm{vee}(\mathrm{hat}(w)) = w$ and $\mathrm{hat}(w)\,x = w \times x$ for ten random $w, x$. Print max error; it should be near machine epsilon.
2. **Exp/log inverse.** For 100 random screw axes and $\theta \in (0, \pi)$, compute $T = e^{[\mathcal{S}]\theta}$, then `logSE3(T)`, and confirm you recover $(\mathcal{S}, \theta)$ to $10^{-10}$. Add a $\theta$ near $0$ and near $\pi$ to expose the singular branches.
3. **Screw from two poses.** Given two poses $T_0, T_1$, compute the relative motion $T = T_1 T_0^{-1}$ and its screw axis via the log map. Verify the axis passes through the expected point for a known pure rotation.
4. **Adjoint check.** Pick a twist $\mathcal{V}_b$ and a transform $T$. Confirm $\mathrm{Ad}_T \mathcal{V}_b$ equals the 6-vector read from $T[\mathcal{V}_b]T^{-1}$. They must match.
5. **Pose interpolation (ScLERP).** Interpolate between $T_0$ and $T_1$ by $T(s) = T_0\,\exp\!\big(s\,\log(T_0^{-1}T_1)\big)$ for $s \in [0,1]$. Render the intermediate frames; the path should be a single smooth screw, not a separate lerp of position and slerp of rotation.
6. **Stone registration step.** Given an ICP rotation+translation as $T$, compute the $se(3)$ log, scale it by a damping factor, re-exponentiate, and apply. Confirm the damped step reduces alignment error monotonically on a synthetic block.

## References

- Kevin M. Lynch and Frank C. Park, *Modern Robotics: Mechanics, Planning, and Control*, Cambridge University Press, 2017. Chapter 3 (rigid-body motions, twists, exponential coordinates, adjoint). http://modernrobotics.northwestern.edu/
- Modern Robotics online resource, 3.3.3 Exponential Coordinates of Rigid-Body Motion. https://modernrobotics.northwestern.edu/nu-gm-book-resource/3-3-3-exponential-coordinates-of-rigid-body-motion/
- Modern Robotics online resource, 3.3.2 Twists (Parts 1 and 2). https://modernrobotics.northwestern.edu/nu-gm-book-resource/3-3-2-twists-part-1-of-2/
- Modern Robotics online resource, 3.3.1 Homogeneous Transformation Matrices (adjoint definition). https://modernrobotics.northwestern.edu/nu-gm-book-resource/3-3-1-homogeneous-transformation-matrices/
- Robotics Unveiled, Representations of Rigid Body Motions: Exponential Coordinates (closed-form $G(\theta)$ and $G^{-1}(\theta)$). https://www.roboticsunveiled.com/robotics-representations-of-rigid-body-motions-exponential-coordinates/
- Eigen documentation (Geometry and Dense modules). https://eigen.tuxfamily.org/dox/
- Ghadi Nehme, Eamon Whalen, Faez Ahmed, *CADFit: Precise Mesh-to-CAD Program Generation with Hybrid Optimization*, arXiv:2605.01171. https://arxiv.org/abs/2605.01171 and code https://github.com/ghadinehme/CADFit

## Flashcard seeds

- Q: Derive the linear part $v$ of a pure-rotation screw axis (axis through point $q$, no pitch). <br> <b>Formula:</b> a body-fixed point at $r$ relative to axis point $q$ has velocity $\dot r = \omega \times r$; target result $v = -\omega \times q$. <br> <i>Hint:</i> apply the velocity rule to the point currently at the origin $O$, so $r = O - q$, then set $O = 0$. :: A: offset of origin from axis point: $r = O - q$ <br> its velocity: $v = \omega \times (O - q)$ <br> set $O = 0$: $v = \omega \times (-q) = -\omega \times q$ <br> $\boxed{v = -\omega \times q}$ (pure rotation, $h=0$).
- Q: Extend the screw linear part to include pitch (slide along the axis). <br> <b>Formula:</b> start from $v = -\omega \times q$; a pitch $h$ adds translation $h\,\omega$ per radian; target $v = -\omega \times q + h\,\omega$. <br> <i>Hint:</i> the pitch term is parallel to $\omega$, so just superpose it onto the rotation term. :: A: rotation contributes $-\omega \times q$ <br> pitch slide contributes $+h\,\omega$ (parallel to the axis) <br> superpose: $v = -\omega \times q + h\,\omega$ <br> $\boxed{v = -\omega \times q + h\,\omega}$.
- Q: Derive the power reduction identity that collapses Rodrigues' series (step 1 of 2). <br> <b>Formula:</b> for a unit axis, $[\omega]^2 = \omega\omega^\top - I$; target $[\omega]^3 = -[\omega]$. <br> <i>Hint:</i> multiply $[\omega]\cdot[\omega]^2$ and use $[\omega]\,\omega = \omega\times\omega = 0$. :: A: for a unit axis, $[\omega]^2 = \omega\omega^\top - I$ <br> $[\omega]^3 = [\omega](\omega\omega^\top - I) = ([\omega]\omega)\omega^\top - [\omega] = 0 - [\omega]$ <br> so higher powers cycle: $[\omega]^4 = -[\omega]^2$, $[\omega]^5 = [\omega]$, ... <br> $\boxed{[\omega]^3 = -[\omega]}$.
- Q: Collapse the matrix-exponential series into Rodrigues' formula (step 2 of 2). <br> <b>Formula:</b> $e^{[\omega]\theta} = I + \theta[\omega] + \tfrac{\theta^2}{2}[\omega]^2 + \tfrac{\theta^3}{6}[\omega]^3 + \cdots$ with $[\omega]^3 = -[\omega]$; target $R = I + \sin\theta\,[\omega] + (1-\cos\theta)[\omega]^2$. <br> <i>Hint:</i> group the odd powers (all multiples of $[\omega]$) and the even powers past $I$ (all multiples of $[\omega]^2$) and recognize the two Taylor series. :: A: odd powers reduce to $\pm[\omega]$: coefficient $\theta - \tfrac{\theta^3}{6} + \cdots = \sin\theta$ <br> even powers past $I$ reduce to $\pm[\omega]^2$: coefficient $\tfrac{\theta^2}{2} - \tfrac{\theta^4}{24} + \cdots = 1 - \cos\theta$ <br> regroup: $I + \sin\theta\,[\omega] + (1-\cos\theta)[\omega]^2$ <br> $\boxed{R = e^{[\omega]\theta} = I + \sin\theta\,[\omega] + (1-\cos\theta)[\omega]^2}$.
- Q: Derive the translation factor $G(\theta)$ in the SE(3) exponential. <br> <b>Formula:</b> $p = \big(\int_0^\theta e^{[\omega]s}\,ds\big)v$ with $e^{[\omega]s} = I + \sin s\,[\omega] + (1-\cos s)[\omega]^2$; target $G(\theta) = I\theta + (1-\cos\theta)[\omega] + (\theta-\sin\theta)[\omega]^2$. <br> <i>Hint:</i> integrate the Rodrigues form term by term over $s$ from $0$ to $\theta$. :: A: $\int_0^\theta I\,ds = I\theta$ <br> $\int_0^\theta \sin s\,[\omega]\,ds = (1-\cos\theta)[\omega]$ <br> $\int_0^\theta (1-\cos s)[\omega]^2\,ds = (\theta - \sin\theta)[\omega]^2$ <br> sum: $G(\theta) = I\theta + (1-\cos\theta)[\omega] + (\theta-\sin\theta)[\omega]^2$ <br> $\boxed{e^{[\mathcal{S}]\theta} = \begin{bmatrix} e^{[\omega]\theta} & G(\theta)v \\ 0 & 1 \end{bmatrix}}$.
- Q: Derive how to recover $\theta$ from $R$ in the log map. <br> <b>Formula:</b> $R = I + \sin\theta[\omega] + (1-\cos\theta)[\omega]^2$; traces $\mathrm{tr}\,I=3$, $\mathrm{tr}\,[\omega]=0$, $\mathrm{tr}\,[\omega]^2=-2$; target $\theta = \cos^{-1}\!\big(\tfrac{\mathrm{tr}\,R-1}{2}\big)$. <br> <i>Hint:</i> take the trace of both sides of Rodrigues, then solve for $\cos\theta$. :: A: $\mathrm{tr}\,R = 3 + \sin\theta\cdot 0 + (1-\cos\theta)(-2)$ <br> $= 3 - 2 + 2\cos\theta = 1 + 2\cos\theta$ <br> solve: $\cos\theta = \tfrac{\mathrm{tr}\,R - 1}{2}$ <br> $\boxed{\theta = \cos^{-1}\!\big(\tfrac{\mathrm{tr}\,R - 1}{2}\big)}$.
- Q: Derive how to recover the axis $[\omega]$ from $R$ in the log map. <br> <b>Formula:</b> $R = I + \sin\theta[\omega] + (1-\cos\theta)[\omega]^2$ with $[\omega]^\top = -[\omega]$ and $([\omega]^2)^\top = [\omega]^2$; target $[\omega] = \tfrac{1}{2\sin\theta}(R - R^\top)$. <br> <i>Hint:</i> write $R^\top$ from Rodrigues, then subtract $R - R^\top$ so the symmetric parts cancel. :: A: $R^\top = I - \sin\theta[\omega] + (1-\cos\theta)[\omega]^2$ <br> subtract: $R - R^\top = 2\sin\theta\,[\omega]$ (symmetric parts cancel) <br> solve: $[\omega] = \tfrac{1}{2\sin\theta}(R - R^\top)$ <br> $\boxed{[\omega] = \tfrac{1}{2\sin\theta}(R - R^\top)}$ (valid for $\sin\theta \neq 0$).
- Q: Build the $6\times6$ adjoint $\mathrm{Ad}_T$ from the conjugation identity. <br> <b>Formula:</b> $[\mathcal{V}_s] = T[\mathcal{V}_b]T^{-1}$ with $T=(R,p)$; target $\mathrm{Ad}_T = \begin{bmatrix} R & 0 \\ [p]R & R \end{bmatrix}$. <br> <i>Hint:</i> expand the $4\times4$ block product and read off the $\omega$-row and $v$-row separately; recall $p\times(R\omega_b) = [p]R\,\omega_b$. :: A: expanding gives $\omega_s = R\,\omega_b$ <br> and $v_s = R\,v_b + p \times (R\,\omega_b) = [p]R\,\omega_b + R\,v_b$ <br> stack as a $6\times6$ acting on $[\omega_b;v_b]$ <br> $\boxed{\mathrm{Ad}_T = \begin{bmatrix} R & 0 \\ [p]R & R \end{bmatrix}}$.
- Q: Build the screw axis for a $\pi/2$ rotation about $z$ through $q=(2,0,0)$, pitch $h=0$. <br> <b>Formula:</b> $\omega=(0,0,1)$, $v = -\omega \times q$, $\mathcal{S} = [\omega;\,v]$. <br> <i>Hint:</i> compute the cross product $\omega \times q$ first, then negate it. :: A: $\omega \times q = (0,0,1)\times(2,0,0) = (0,2,0)$ <br> $v = -(0,2,0) = (0,-2,0)$ <br> stack: $\mathcal{S} = [\omega;v] = (0,0,1,\,0,-2,0)$ <br> $\boxed{\mathcal{S} = (0,0,1,\;0,-2,0)}$.
- Q: Find where the origin lands under that $\pi/2$-about-$z$ screw. <br> <b>Formula:</b> $p = G(\theta)v$, $G(\theta) = I\theta + (1-\cos\theta)[\omega] + (\theta-\sin\theta)[\omega]^2$, with $\mathcal{S}=(0,0,1,0,-2,0)$, $\theta=\pi/2$. <br> <i>Hint:</i> evaluate the scalars first: $(1-\cos\theta)=1$, $(\theta-\sin\theta)\approx0.5708$, and $[\omega]^2=\mathrm{diag}(-1,-1,0)$. :: A: scalars: $(1-\cos\theta)=1$, $(\theta-\sin\theta)\approx 0.5708$, $[\omega]^2=\mathrm{diag}(-1,-1,0)$ <br> $p = G(\theta)v = (2,-2,0)$ <br> origin maps to $T(0,0,0,1)^\top = (2,-2,0,1)^\top$ <br> $\boxed{O \mapsto (2,-2,0)}$.
- Q: Worked example (large scale, $\theta=\pi$). A quarry block rotates $\theta=\pi$ about $z$ through $q=(5,0,0)$ m, $h=0$. Build $\mathcal{S}$ and find where the origin lands. <br> <b>Formula:</b> $v = -\omega\times q$; for $\theta=\pi$ the rotation flips the offset-from-center sign. <br> <i>Hint:</i> compute $v$ first; then the origin's offset from the center $(5,0)$ is $(-5,0)$, which a $180^\circ$ turn sends to $(5,0)$. :: A: $v = -\omega \times q = -(0,0,1)\times(5,0,0) = (0,-5,0)$ <br> offset of $O$ from center $(5,0)$ is $(-5,0)$; $180^\circ$ sends it to $(5,0)$ <br> $O \mapsto (5,0)+(5,0) = (10,0,0)$ m <br> $\boxed{\mathcal{S}=(0,0,1,\,0,-5,0),\ O\mapsto(10,0,0)\text{ m}}$.
- Q: Worked example (pure translation, shop/mm scale). A slide of $\theta = 30$ mm along $x$: $\omega=0$, $v=(1,0,0)$. Build $T = e^{[\mathcal{S}]\theta}$. <br> <b>Formula:</b> pure-translation branch $e^{[\mathcal{S}]\theta} = \begin{bmatrix} I & v\theta \\ 0 & 1 \end{bmatrix}$. <br> <i>Hint:</i> $\omega=0$ means no rotation block; just compute $v\theta$ for the translation column. :: A: $\omega = 0$ triggers the pure-translation branch <br> $v\theta = (1,0,0)\cdot 30 = (30,0,0)$ mm <br> $e^{[\mathcal{S}]\theta} = \begin{bmatrix} I & v\theta \\ 0 & 1 \end{bmatrix}$ <br> $\boxed{T = \begin{bmatrix} I & (30,0,0)^\top \\ 0 & 1 \end{bmatrix}}$.
- Q: Worked example (small-angle, mm-scale ICP). A fine correction: $\theta = 0.01$ rad about $z$ through $q=(0,0,0)$, $h=0$. Approximate the resulting $R$ to first order. <br> <b>Formula:</b> $R = I + \sin\theta[\omega] + (1-\cos\theta)[\omega]^2$ with $\sin\theta\approx\theta$, $1-\cos\theta\approx\theta^2/2$. <br> <i>Hint:</i> first note $v=-\omega\times q=0$ since the axis passes through the origin; then drop the $\theta^2$ term and keep $R \approx I + \theta[\omega]$. :: A: $v = -\omega \times q = 0$ (axis through origin) <br> small-angle: $\sin\theta\approx\theta$, $1-\cos\theta\approx\theta^2/2\approx 0$ <br> $R \approx I + \theta[\omega] = I + 0.01[\omega]$ <br> $\boxed{R \approx \begin{bmatrix} 1 & -0.01 & 0 \\ 0.01 & 1 & 0 \\ 0 & 0 & 1 \end{bmatrix}}$.
- Q: Closing test (close the loop). For $\theta=\pi/2$ about $z$ through $(2,0,0)$, verify the origin image $(2,-2,0)$ by geometry, then round-trip through the log map. <br> <b>Formula:</b> geometry — rotate offset-from-center by $90^\circ$; log map — $\theta=\cos^{-1}(\tfrac{\mathrm{tr}\,R-1}{2})$, $[\omega]=\tfrac{1}{2\sin\theta}(R-R^\top)$, $v=G^{-1}(\theta)p$. <br> <i>Hint:</i> a pass means both the geometric image and the recovered $(\mathcal{S},\theta)$ match the originals exactly. :: A: geometry: offset of $O$ from center $(2,0)$ is $(-2,0)$; rotate $+90^\circ\Rightarrow(0,-2)$; add center $\Rightarrow(2,-2)$ ✓ <br> log map: $\mathrm{tr}\,R=1\Rightarrow\theta=\cos^{-1}(0)=\pi/2$ ✓ <br> $\tfrac{1}{2\sin\theta}(R-R^\top)\Rightarrow\omega=(0,0,1)$, then $v=G^{-1}(\theta)p=(0,-2,0)$ ✓ <br> $\boxed{\text{forward and inverse agree; screw recovered exactly}}$.
- Q: Closing test (units + limit check). For the quarry block ($\theta=\pi$ about $z$ through $q=(5,0,0)$ m), confirm $O\mapsto(10,0,0)$ m has consistent units and survives the $\theta=0$ limit. <br> <b>Formula:</b> $v=-\omega\times q$ ($[\,1\,]\times[\text{m}]=[\text{m}]$); limit $G(\theta)\to 0$ as $\theta\to 0$ so $p\to 0$. <br> <i>Hint:</i> a pass means $v$ and $p$ carry meters, a zero rotation moves nothing, and the result is the mirror of $O$ across the axis at $x=5$. :: A: $v = -\omega\times q$: $[1]\times[\text{m}] = [\text{m}]$, so $v$ and $p$ are in meters ✓ <br> limit $\theta\to 0$: $G(\theta)\to 0$ so $p\to 0$, origin stays put ✓ <br> at $\theta=\pi$ the result $(10,0,0)$ is the reflection of $O$ through the axis at $x=5$ ✓ <br> $\boxed{\text{units in m, limit consistent: result valid}}$.
- Q: What does Chasles' theorem say about any rigid-body displacement? <br> <i>Hint:</i> think of tightening a bolt — one axis, turn plus slide together. :: A: It equals a single screw motion: a rotation about a fixed axis plus a translation along that same axis.
- Q: When should you use the matrix exponential / log map instead of just storing $R$ and $p$? <br> <i>Hint:</i> think about what coordinate you want to optimize or interpolate over without leaving the manifold. :: A: When you need a minimal singularity-free 6-number coordinate (a twist) to optimize, interpolate, or average poses on the $SE(3)$ manifold without drift, e.g. ICP scan-to-CAD registration.
- Q: Pitfall — why can't you compose two screw motions by adding their twists? <br> <b>Formula:</b> $e^{[\mathcal{S}_1]\theta_1}e^{[\mathcal{S}_2]\theta_2} \neq e^{[\mathcal{S}_1]\theta_1+[\mathcal{S}_2]\theta_2}$. <br> <i>Hint:</i> recall what property of the matrix exponential fails when the matrices do not commute. :: A: The matrix exponential is non-commutative, so adding the exponents is wrong; compose instead by multiplying the resulting $T$ matrices.
- Q: Pitfall — which $\theta$ values trigger singular branches in the SE(3) log map, and what do you do? <br> <b>Formula:</b> $[\omega]=\tfrac{1}{2\sin\theta}(R-R^\top)$ and $G^{-1}$ divides by $\theta$ and $\sin\theta$. <br> <i>Hint:</i> look at where the denominators $\theta$ and $\sin\theta$ vanish. :: A: $\theta=0$: use the pure-translation branch ($G^{-1}$ divides by $\theta$ and $\sin\theta$). <br> $\theta=\pi$: $\sin\theta=0$ blows up $(R-R^\top)/(2\sin\theta)$, so extract $\omega$ from a positive-diagonal column of $R$.
- Q: Pitfall — what guard runs in code before the unit-axis exp formula? <br> <b>Formula:</b> test `w.norm() < kEps`; pure-translation form $\begin{bmatrix} I & v\theta \\ 0 & 1 \end{bmatrix}$. <br> <i>Hint:</i> the closed-form $G(\theta)$ assumes $\|\omega\|=1$, so catch the near-zero-axis case before dividing by the angle. :: A: Check `w.norm() < kEps`; if true, take the pure-translation branch instead of dividing by the angle. <br> Scale `kEps` to model size (mm vs m) so the tolerance is neither too tight nor too loose.
- CLOZE: A twist is the 6-vector {{c1::$\mathcal{V} = [\omega;\,v]$}} packing angular velocity and linear velocity, with matrix form $[\mathcal{V}] = [[\omega],v;\,0,0] \in se(3)$.
- CLOZE: For a normalized screw with rotation, scale the axis so {{c1::$\|\omega\| = 1$}}, and then the scalar $\theta$ is exactly the rotation angle.
- CLOZE: For a pure-translation screw, $\omega = 0$ and {{c1::$\|v\| = 1$}}, giving $e^{[\mathcal{S}]\theta} = [I,\ v\theta;\ 0,\ 1]$.
- CLOZE: The hat operator satisfies {{c1::$[\omega]\,x = \omega \times x$}}, turning a cross product into a matrix multiply.
- CLOZE: The maps $\dot T T^{-1} = [\mathcal{V}_s]$ and $T^{-1}\dot T = [\mathcal{V}_b]$ define the {{c1::spatial}} and body twists respectively. <br> <i>hint:</i> $\mathcal{V}_s$ uses $\dot T$ on the left ("fixed/space" frame).
- CLOZE: The exponential coordinates of a displacement are the 6-vector {{c1::$\mathcal{S}\theta$}}, the screw axis times the distance traveled.
