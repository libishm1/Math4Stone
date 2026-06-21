# Forward kinematics: DH and Product of Exponentials

Prereqs: [[11_twists-screws-exp]], [[12_transforms-se3]] => this => [[14_manipulator-jacobian]], [[15_inverse-kinematics]]

## Why it matters for the stone scan->CAD->IFC / CADFit workflow

A scan-to-CAD-to-IFC pipeline for irregular stone runs on a robot or a multi-axis CNC, and that machine is a serial chain of joints. Forward kinematics (FK) is the function that turns the joint angles your motion planner outputs into the actual pose of the cutting tool or the depth camera flange in world coordinates. You need FK twice: once to know where a hand-held or arm-mounted scanner sees the block from (so you can register partial scans into one cloud), and once to confirm the toolpath you generated against the reconstructed CAD model lands the bit where the geometry says it should. In a CADFit-style synthesis loop, the candidate CAD program is evaluated against the as-built stone, and the as-built measurement only exists because FK placed the sensor and the probe in a shared frame. Get the chain of transforms wrong by one frame convention and every downstream IoU score, IFC placement, and cut offset inherits that error silently.

## Intuition

A robot arm is a stack of turntables and sliders bolted end to end. Each joint adds one rigid motion to the part of the arm beyond it. To find where the tool ends up, you ride the chain from the base outward: stand at the base, apply joint 1's motion, then from that new vantage apply joint 2's motion, and so on, multiplying the motions together. Because rigid motions are 4x4 matrices, "apply then apply" is matrix multiplication, and order matters: turning then sliding is not the same as sliding then turning.

Analogy: think of nested Russian dolls where each shell can twist independently. The painted face on the innermost doll (the tool tip) ends up pointing somewhere that depends on every twist applied to every shell around it, composed from the outside in. FK just bookkeeps those twists exactly.

Two bookkeeping schemes exist. DH assigns four numbers per joint and chains per-link frames. Product of Exponentials (PoE) writes each joint as a screw motion about a fixed axis and multiplies one exponential per joint times a single home pose. PoE needs fewer arbitrary frame choices and treats sliding and twisting joints with the same formula, which is why production code prefers it.

## The math (first principles)

### Setup, frames, conventions

We work in $SE(3)$. A pose is a homogeneous transform
$$T = \begin{bmatrix} R & p \\ 0 & 1 \end{bmatrix} \in SE(3), \qquad R \in SO(3),\; p \in \mathbb{R}^3.$$
Conventions used throughout: right-handed frames, column vectors, transforms act on the left ($p' = Tp$ with $p$ in homogeneous coordinates), and angles in radians. Rotations are active (they move the body, not the frame description). Lengths carry whatever unit the model uses; for stone fabrication keep the chain in metres at the cell scale and convert at the IFC boundary. A serial chain has $n$ one-degree-of-freedom joints, indexed $1$ at the base to $n$ at the tool. The fixed base frame is $\{s\}$ (space); the end-effector frame is $\{b\}$ (body).

### Forward kinematics as an ordered product

Each joint-link pair contributes a transform that depends on its joint variable. The end-effector pose is the ordered product
$$T(\theta) = A_1(\theta_1)\, A_2(\theta_2) \cdots A_n(\theta_n).$$
This is exact, not an approximation. The only content of a kinematics scheme is how you parameterize each $A_i$.

### Classical Denavit-Hartenberg

DH assigns four parameters per joint: $a_i$ (link length), $\alpha_i$ (link twist), $d_i$ (link offset), and $\theta_i$ (joint angle). The standard (distal) homogeneous transform is
$$A_i = \operatorname{Rot}_z(\theta_i)\, \operatorname{Trans}_z(d_i)\, \operatorname{Trans}_x(a_i)\, \operatorname{Rot}_x(\alpha_i).$$
Because a translation and a rotation about the *same* axis commute, $\operatorname{Trans}_z(d_i)$ and $\operatorname{Rot}_z(\theta_i)$ may be written in either order; both give the same $A_i$. Written out:
$$A_i = \begin{bmatrix}
\cos\theta_i & -\sin\theta_i\cos\alpha_i & \sin\theta_i\sin\alpha_i & a_i\cos\theta_i \\
\sin\theta_i & \cos\theta_i\cos\alpha_i & -\cos\theta_i\sin\alpha_i & a_i\sin\theta_i \\
0 & \sin\alpha_i & \cos\alpha_i & d_i \\
0 & 0 & 0 & 1
\end{bmatrix}.$$
For a revolute joint $\theta_i$ is the variable and $d_i$ is fixed. For a prismatic joint $d_i$ is the variable and $\theta_i$ is fixed. Then $T(\theta) = A_1 A_2 \cdots A_n$. Warning: a competing convention (Craig, "modified" or proximal DH) reverses the factor order, $A_i = \operatorname{Rot}_x(\alpha_{i-1})\operatorname{Trans}_x(a_{i-1})\operatorname{Rot}_z(\theta_i)\operatorname{Trans}_z(d_i)$, and shifts the parameter indices. A DH table is meaningless without stating which convention produced it.

### Product of Exponentials

PoE writes each joint as a screw motion. A screw axis is the twist
$$\mathcal{S} = \begin{bmatrix} \omega \\ v \end{bmatrix} \in \mathbb{R}^6, \qquad \omega, v \in \mathbb{R}^3,$$
normalized so that for a revolute joint $\|\omega\| = 1$ and for a prismatic joint $\omega = 0$ with $\|v\| = 1$. The screw axes are expressed in the fixed space frame $\{s\}$ at the home (zero) configuration.

- Revolute joint with unit axis direction $\omega$ through a point $q$ on the axis: $v = -\omega \times q$. This $v$ is the linear velocity of the space-frame origin point under unit angular velocity about the axis.
- Prismatic joint sliding along unit direction $\hat{v}$: $\omega = 0$, $v = \hat{v}$.

The $4\times4$ matrix form of a twist is
$$[\mathcal{S}] = \begin{bmatrix} [\omega] & v \\ 0 & 0 \end{bmatrix}, \qquad
[\omega] = \begin{bmatrix} 0 & -\omega_3 & \omega_2 \\ \omega_3 & 0 & -\omega_1 \\ -\omega_2 & \omega_1 & 0 \end{bmatrix},$$
where $[\omega]$ is the skew-symmetric ("hat") matrix with $[\omega]x = \omega \times x$.

The space-frame Product of Exponentials formula is
$$\boxed{\,T(\theta) = e^{[\mathcal{S}_1]\theta_1}\, e^{[\mathcal{S}_2]\theta_2} \cdots e^{[\mathcal{S}_n]\theta_n}\, M\,}$$
where $M = T(0) \in SE(3)$ is the home configuration: the pose of $\{b\}$ relative to $\{s\}$ when every joint variable is zero. Reading right to left, $M$ places the tool at home, then each exponential, applied from the base side, screws the whole outboard chain about its space-frame axis.

### The matrix exponential of a twist (closed form)

For a unit rotation axis $\omega$, Rodrigues' formula gives the rotation part:
$$e^{[\omega]\theta} = I + \sin\theta\,[\omega] + (1-\cos\theta)\,[\omega]^2.$$
For the full twist exponential with $\|\omega\| = 1$,
$$e^{[\mathcal{S}]\theta} = \begin{bmatrix} e^{[\omega]\theta} & G(\theta)\,v \\ 0 & 1 \end{bmatrix}, \qquad
G(\theta) = I\theta + (1-\cos\theta)[\omega] + (\theta - \sin\theta)[\omega]^2.$$
For a prismatic joint ($\omega = 0$) the exponential degenerates to a pure translation:
$$e^{[\mathcal{S}]\theta} = \begin{bmatrix} I & v\,\theta \\ 0 & 1 \end{bmatrix}.$$
This single pair of cases handles both joint types with no special bookkeeping, the practical reason PoE wins in code.

### DH vs PoE tradeoff

DH is compact (four numbers per joint) and pervasive in industrial datasheets, so you must read it. It is fragile: it forces per-link frame placement rules, has singular parameterizations when consecutive axes are parallel, and splits revolute and prismatic into different rows. PoE needs only $\{s\}$, $\{b\}$, and the home pose $M$; the screw axes have a direct geometric meaning; the same exponential serves both joint types; and the Jacobian (next topic) falls out of the same exponentials via the adjoint. Prefer PoE in code. Keep a DH-to-PoE converter for when you only have a datasheet.

## Worked example

Take a planar 2R arm (two revolute joints, both axes along $+z$, link lengths $L_1 = 2$, $L_2 = 1$), and evaluate at $\theta_1 = 90^\circ$, $\theta_2 = 0^\circ$. We expect the tip to swing from the $+x$ axis up to point along $+y$, landing at $(0, 3)$ since the arm is straight.

### PoE setup

Home configuration: arm stretched along $+x$, tool at $(L_1 + L_2, 0) = (3, 0)$, tool frame aligned with $\{s\}$:
$$M = \begin{bmatrix} 1 & 0 & 0 & 3 \\ 0 & 1 & 0 & 0 \\ 0 & 0 & 1 & 0 \\ 0 & 0 & 0 & 1 \end{bmatrix}.$$
Screw axes, both with $\omega = (0,0,1)$. Axis 1 passes through the origin $q_1 = (0,0,0)$, so $v_1 = -\omega \times q_1 = 0$. Axis 2 passes through $q_2 = (2,0,0)$:
$$v_2 = -\omega \times q_2 = -\begin{bmatrix}0\\0\\1\end{bmatrix}\times\begin{bmatrix}2\\0\\0\end{bmatrix} = -\begin{bmatrix}0\\2\\0\end{bmatrix} = \begin{bmatrix}0\\-2\\0\end{bmatrix}.$$
So $\mathcal{S}_1 = (0,0,1,\;0,0,0)$ and $\mathcal{S}_2 = (0,0,1,\;0,-2,0)$.

### Exponentials

With $\theta_2 = 0$, $e^{[\mathcal{S}_2]\theta_2} = I$. For $\theta_1 = 90^\circ = \pi/2$, $\omega = (0,0,1)$ so $\sin\theta_1 = 1$, $\cos\theta_1 = 0$:
$$e^{[\omega]\theta_1} = \begin{bmatrix} 0 & -1 & 0 \\ 1 & 0 & 0 \\ 0 & 0 & 1 \end{bmatrix}.$$
Since $v_1 = 0$, the translation part $G(\theta_1)v_1 = 0$, so
$$e^{[\mathcal{S}_1]\theta_1} = \begin{bmatrix} 0 & -1 & 0 & 0 \\ 1 & 0 & 0 & 0 \\ 0 & 0 & 1 & 0 \\ 0 & 0 & 0 & 1 \end{bmatrix}.$$

### Compose

$$T = e^{[\mathcal{S}_1]\theta_1}\, I \, M
= \begin{bmatrix} 0 & -1 & 0 & 0 \\ 1 & 0 & 0 & 0 \\ 0 & 0 & 1 & 0 \\ 0 & 0 & 0 & 1 \end{bmatrix}
\begin{bmatrix} 1 & 0 & 0 & 3 \\ 0 & 1 & 0 & 0 \\ 0 & 0 & 1 & 0 \\ 0 & 0 & 0 & 1 \end{bmatrix}
= \begin{bmatrix} 0 & -1 & 0 & 0 \\ 1 & 0 & 0 & 3 \\ 0 & 0 & 1 & 0 \\ 0 & 0 & 0 & 1 \end{bmatrix}.$$
Tip position $p = (0, 3, 0)$. Matches the hand prediction.

### Cross-check with DH

DH for this planar arm: $\theta_i$ variable, $d_i = 0$, $\alpha_i = 0$, $a_1 = 2$, $a_2 = 1$. Tip position by the planar formula $x = a_1\cos\theta_1 + a_2\cos(\theta_1+\theta_2)$, $y = a_1\sin\theta_1 + a_2\sin(\theta_1+\theta_2)$. With $\theta_1 = 90^\circ, \theta_2 = 0$: $x = 2\cdot 0 + 1\cdot 0 = 0$, $y = 2\cdot 1 + 1\cdot 1 = 3$. Same answer, as it must be.

## Code

C++ with Eigen first. The FK kernel is three functions: hat of a 3-vector, the SE(3) twist exponential, and the left-to-right product.

```cpp
// fk_poe.hpp  -- Product-of-Exponentials forward kinematics (space frame).
// Build: g++ -std=c++17 -I /path/to/eigen fk_poe.cpp
#pragma once
#include <Eigen/Dense>
#include <vector>
#include <stdexcept>
#include <cmath>

namespace fk {

using Vec3 = Eigen::Vector3d;
using Vec6 = Eigen::Matrix<double, 6, 1>;   // screw axis S = [omega; v]
using Mat3 = Eigen::Matrix3d;
using Mat4 = Eigen::Matrix4d;               // a pose T in SE(3)

// [omega] : skew matrix so that hat(w) * x == w.cross(x)
inline Mat3 hat(const Vec3& w) {
    Mat3 W;
    W <<    0.0, -w.z(),  w.y(),
         w.z(),    0.0, -w.x(),
        -w.y(),  w.x(),    0.0;
    return W;
}

// e^{[omega] theta} : Rodrigues, omega assumed unit.
inline Mat3 expSO3(const Vec3& w, double theta) {
    const Mat3 W = hat(w);
    // R = I + sin(theta)[w] + (1 - cos(theta))[w]^2
    return Mat3::Identity()
         + std::sin(theta) * W
         + (1.0 - std::cos(theta)) * (W * W);
}

// e^{[S] theta} : closed-form twist exponential. S = [omega; v].
// Revolute: ||omega|| == 1. Prismatic: omega == 0 -> pure translation v*theta.
inline Mat4 expSE3(const Vec6& S, double theta) {
    const Vec3 w = S.head<3>();
    const Vec3 v = S.tail<3>();
    Mat4 T = Mat4::Identity();

    if (w.norm() < 1e-12) {                 // prismatic joint
        T.block<3, 1>(0, 3) = v * theta;
        return T;
    }

    const Mat3 W = hat(w);                  // omega is unit for revolute
    // G(theta) = I*theta + (1 - cos)[w] + (theta - sin)[w]^2
    const Mat3 G = Mat3::Identity() * theta
                 + (1.0 - std::cos(theta)) * W
                 + (theta - std::sin(theta)) * (W * W);
    T.block<3, 3>(0, 0) = expSO3(w, theta);
    T.block<3, 1>(0, 3) = G * v;
    return T;
}

// T(theta) = e^{[S1]t1} e^{[S2]t2} ... e^{[Sn]tn} M   (space frame PoE)
inline Mat4 fkSpace(const std::vector<Vec6>& Slist,
                    const std::vector<double>& theta,
                    const Mat4& M) {
    if (Slist.size() != theta.size())
        throw std::runtime_error("Slist/theta size mismatch");
    Mat4 T = Mat4::Identity();
    for (std::size_t i = 0; i < Slist.size(); ++i)
        T = T * expSE3(Slist[i], theta[i]);   // left-multiply, base side first
    return T * M;                             // home configuration on the right
}

} // namespace fk
```

Driver reproducing the worked example:

```cpp
// fk_poe.cpp
#include "fk_poe.hpp"
#include <iostream>

int main() {
    using namespace fk;
    Vec6 S1; S1 << 0,0,1,  0, 0,0;   // axis +z through origin
    Vec6 S2; S2 << 0,0,1,  0,-2,0;   // axis +z through (2,0,0): v = -w x q

    Mat4 M = Mat4::Identity();
    M(0,3) = 3.0;                    // home tip at (3,0,0), L1+L2 = 3

    std::vector<Vec6> Slist{S1, S2};
    std::vector<double> theta{M_PI/2.0, 0.0};

    Mat4 T = fkSpace(Slist, theta, M);
    std::cout << "tip = " << T.block<3,1>(0,3).transpose() << "\n"; // 0 3 0
    return 0;
}
```

A DH builder for the legacy path, so you can convert a datasheet table into the same $T$:

```cpp
// Standard (distal) DH link: Rot_z(theta) Trans_z(d) Trans_x(a) Rot_x(alpha).
inline fk::Mat4 dhLink(double theta, double d, double a, double alpha) {
    using std::sin; using std::cos;
    fk::Mat4 A;
    A << cos(theta), -sin(theta)*cos(alpha),  sin(theta)*sin(alpha), a*cos(theta),
         sin(theta),  cos(theta)*cos(alpha), -cos(theta)*sin(alpha), a*sin(theta),
                0.0,             sin(alpha),             cos(alpha),            d,
                0.0,                    0.0,                    0.0,          1.0;
    return A;
}
// T = A_1 * A_2 * ... * A_n  (chain the links in order).
```

Short Python check with `numpy`, useful for a quick scan-rig calibration script:

```python
import numpy as np

def hat(w):
    return np.array([[0, -w[2], w[1]],
                     [w[2], 0, -w[0]],
                     [-w[1], w[0], 0.0]])

def exp_se3(S, theta):                       # S = [w(3); v(3)]
    w, v = S[:3], S[3:]
    T = np.eye(4)
    if np.linalg.norm(w) < 1e-12:            # prismatic
        T[:3, 3] = v * theta
        return T
    W = hat(w)
    R = np.eye(3) + np.sin(theta)*W + (1-np.cos(theta))*(W@W)
    G = np.eye(3)*theta + (1-np.cos(theta))*W + (theta-np.sin(theta))*(W@W)
    T[:3, :3], T[:3, 3] = R, G @ v
    return T

def fk_space(Slist, theta, M):
    T = np.eye(4)
    for S, t in zip(Slist, theta):
        T = T @ exp_se3(S, t)
    return T @ M

S = [np.array([0,0,1, 0,0,0.]), np.array([0,0,1, 0,-2,0.])]
M = np.eye(4); M[0,3] = 3.0
print(fk_space(S, [np.pi/2, 0.0], M)[:3, 3])  # [0. 3. 0.]
```

## Connections

- [[11_twists-screws-exp]]: a joint motion is a screw, and $e^{[\mathcal{S}]\theta}$ is the matrix exponential of its twist; PoE is literally a product of those exponentials, so this is the algebraic backbone.
- [[12_transforms-se3]]: every $A_i$, every exponential, and the home pose $M$ live in $SE(3)$; FK is just disciplined $SE(3)$ composition, and the left-multiply convention comes from how $SE(3)$ acts on points.
- [[14_manipulator-jacobian]]: differentiating $T(\theta)$ with respect to the joint variables gives the space Jacobian, whose $i$-th column is $\operatorname{Ad}_{(e^{[\mathcal{S}_1]\theta_1}\cdots e^{[\mathcal{S}_{i-1}]\theta_{i-1}})}\mathcal{S}_i$, built from the same exponentials.
- [[15_inverse-kinematics]]: IK inverts FK numerically by Newton or damped least squares on the Jacobian, so a correct FK is the precondition for solving joint angles that hit a target tool pose.
- [[09_rigid-registration-icp]]: the same $SE(3)$ poses FK produces are what you feed into multi-view scan registration, aligning each camera shot of the stone into one cloud.

## Pitfalls and failure modes

- Frame inconsistency. Mixing space-frame and body-frame screw axes, or chaining DH links from a table built in the other convention, silently corrupts $T$. State the convention next to every parameter table. A space-frame screw left-multiplies; a body-frame screw right-multiplies. Do not interleave.
- Order and non-commutativity. $A_i A_j \neq A_j A_i$. Reversing two joints, or applying $M$ on the left instead of the right, gives a plausible-looking but wrong pose. The worked example anchors the correct order.
- Unnormalized screw axes. The closed-form exponential assumes $\|\omega\| = 1$ for revolute joints. If $\|\omega\| \neq 1$, $\theta$ no longer equals the physical angle and $G(\theta)$ is wrong. Normalize $\omega$ at construction and absorb its magnitude into $\theta$.
- Prismatic branch. Forgetting the $\omega = 0$ special case makes Rodrigues divide by a zero rotation angle or return identity translation. Branch on $\|\omega\| < \varepsilon$ explicitly.
- Numerical drift in $R$. Repeated products accumulate floating-point error and $T$'s rotation block stops being orthonormal. Re-project $R$ onto $SO(3)$ (nearest orthonormal via SVD) before exporting a pose to IFC placement, where a non-orthonormal axis triple is invalid.
- Degrees vs radians, units of length. The trig functions take radians; a stray degree input is off by $57\times$. Length units must match the model. For stone work keep the chain in metres and convert once at the IFC/CAD boundary.
- Parameterization singularities. DH is ill-posed when consecutive joint axes are parallel (the common-normal direction is undefined). PoE has no such parameterization singularity; only the Jacobian can be singular, which is a kinematic property, not a bookkeeping artifact.
- Wrong home pose. $M$ must be the actual zero-configuration pose measured on the machine. A guessed $M$ offsets every output pose by a constant rigid error that calibration, not formula, must remove.

## Practice

1. Hand-evaluate the worked 2R arm at $\theta_1 = 0$, $\theta_2 = 90^\circ$. Predict the tip first, then confirm via both the PoE product and the planar DH formula. Verify the tool frame orientation, not only position.
2. Add a third revolute joint coaxial with joint 2 but a prismatic joint instead, sliding along $+x$ with range $[0, 0.5]$. Write its screw ($\omega = 0$, $v = (1,0,0)$), extend `Slist`, and print the tip as the slide goes $0 \to 0.5$. Check the tip moves a straight 0.5 along the rotated $x$ direction.
3. Build the standard DH table for the 2R arm, compute $T$ with `dhLink`, and assert it equals the PoE result to $10^{-9}$ across 20 random joint pairs. This is your regression test that the two paths agree.
4. Take a published UR5 or Panda screw-axis list (Modern Robotics provides UR5e $\mathcal{S}_i$ and $M$). Implement `fkSpace`, evaluate at the all-zero and a random configuration, and compare against the vendor or Robotics Toolbox FK to $10^{-6}$.
5. Corrupt the chain on purpose: feed $\theta$ in degrees, then with an unnormalized $\omega = (0,0,2)$. Record how far the tip moves in each case. This builds intuition for which mistakes are catastrophic versus subtle.
6. Write the $SO(3)$ re-projection (SVD of $R$, set singular values to 1) and apply it after 10000 random product compositions. Measure $\|R^\top R - I\|$ before and after to see drift and its cure.

## References

- Kevin M. Lynch and Frank C. Park, *Modern Robotics: Mechanics, Planning, and Control*, Cambridge University Press, 2017. Chapter 4 (forward kinematics, PoE space and body frames). Full text PDF: https://hades.mech.northwestern.edu/images/2/2a/Park-lynch.pdf
- Modern Robotics course resource, "4.1.1 Product of Exponentials Formula in the Space Frame": https://modernrobotics.northwestern.edu/nu-gm-book-resource/4-1-1-product-of-exponentials-formula-in-the-space-frame/
- Jesse Haviland and Peter Corke, "Manipulator Differential Kinematics: Part 1: Kinematics, Velocity, and Applications," IEEE Robotics & Automation Magazine, 2023. arXiv:2207.01796. https://arxiv.org/abs/2207.01796 (elementary transform sequence, a DH-free FK modelling approach).
- "Product of exponentials formula," Wikipedia (Rodrigues and twist-exponential closed forms): https://en.wikipedia.org/wiki/Product_of_exponentials_formula
- Robotics Toolbox for Python (Corke), ETS models documentation: https://petercorke.github.io/robotics-toolbox-python/arm_ets.html
- Ghadi Nehme et al., "CADFit: Precise Mesh-to-CAD Program Generation with Hybrid Optimization," arXiv:2605.01171. https://arxiv.org/abs/2605.01171 and code https://github.com/ghadinehme/CADFit (the CAD-program-synthesis target the FK-placed measurements feed into).
- Open Cascade Technology reference documentation (TopLoc/gp transforms for placing CAD in world coordinates): https://dev.opencascade.org/doc/overview/html/
- IfcOpenShell documentation (IfcLocalPlacement / object placement in IFC): https://docs.ifcopenshell.org/

## Flashcard seeds

- Q: Each joint $i$ contributes a transform $A_i(\theta_i)$ acting on everything outboard of it. Ride the chain from the base to write the end-effector pose $T(\theta)$. <br> <b>Formula:</b> $T(\theta)=A_1(\theta_1)\,A_2(\theta_2)\cdots A_n(\theta_n)$ <br> <i>Hint:</i> start at the base with $A_1$, then right-multiply each next joint's motion as you move outward. :: A: stand at the base and apply joint 1's motion $A_1$ <br> from that new vantage apply joint 2's motion, composing on the right: $A_1 A_2$ <br> continue outward to joint $n$, each new motion right-multiplied <br> $\boxed{T(\theta) = A_1(\theta_1)\,A_2(\theta_2)\cdots A_n(\theta_n)}$ exact, not an approximation.
- Q: Build the standard DH link matrix, step 1 of 2: compose the first three screw moves and give the translation column so far. <br> <b>Formula:</b> $A_i'=\operatorname{Rot}_z(\theta_i)\operatorname{Trans}_z(d_i)\operatorname{Trans}_x(a_i)$ <br> <i>Hint:</i> $\operatorname{Rot}_z$ fills the upper-left $2\times2$ with $\cos\theta_i,\sin\theta_i$; then read off where each translation lands $p$. :: A: $\operatorname{Rot}_z(\theta_i)$ puts $\cos\theta_i,\sin\theta_i$ in the upper-left $2\times2$ <br> $\operatorname{Trans}_z(d_i)$ adds $d_i$ in the $z$ row of $p$ (commutes with $\operatorname{Rot}_z$) <br> $\operatorname{Trans}_x(a_i)$ adds an $x$-offset rotated by $\theta_i$: $p_x=a_i\cos\theta_i,\;p_y=a_i\sin\theta_i$ <br> $\boxed{A_i' = \operatorname{Rot}_z(\theta_i)\operatorname{Trans}_z(d_i)\operatorname{Trans}_x(a_i),\quad p=(a_i\cos\theta_i,\,a_i\sin\theta_i,\,d_i)}$.
- Q: Finish the DH link matrix, step 2 of 2: right-multiply the step-1 result by $\operatorname{Rot}_x(\alpha_i)$ and write the full $A_i$. <br> <b>Formula:</b> $A_i=\begin{bmatrix}\cos\theta_i & -\sin\theta_i\cos\alpha_i & \sin\theta_i\sin\alpha_i & a_i\cos\theta_i\\ \sin\theta_i & \cos\theta_i\cos\alpha_i & -\cos\theta_i\sin\alpha_i & a_i\sin\theta_i\\ 0 & \sin\alpha_i & \cos\alpha_i & d_i\\ 0&0&0&1\end{bmatrix}$ <br> <i>Hint:</i> $\operatorname{Rot}_x(\alpha_i)$ only mixes columns 2 and 3 by $\cos\alpha_i,\sin\alpha_i$; column 1 and the translation $p$ are untouched. :: A: $\operatorname{Rot}_x(\alpha_i)$ mixes the $y,z$ columns by $\cos\alpha_i,\sin\alpha_i$, leaving $p$ unchanged <br> column 1 stays $(\cos\theta_i,\sin\theta_i,0)$; columns 2,3 pick up the $\alpha_i$ factors <br> $\boxed{A_i=\begin{bmatrix}\cos\theta_i & -\sin\theta_i\cos\alpha_i & \sin\theta_i\sin\alpha_i & a_i\cos\theta_i\\ \sin\theta_i & \cos\theta_i\cos\alpha_i & -\cos\theta_i\sin\alpha_i & a_i\sin\theta_i\\ 0 & \sin\alpha_i & \cos\alpha_i & d_i\\ 0&0&0&1\end{bmatrix}}$.
- Q: A point coincident with the space origin moves under unit angular velocity $\omega$ about an axis through point $q$. Derive the screw's linear part $v$. <br> <b>Formula:</b> $\dot r=\omega\times(r-q)$, target $v=-\omega\times q$ <br> <i>Hint:</i> substitute $r=0$ (the space origin) into the velocity formula. :: A: velocity of a point $r$ on a body rotating about an axis through $q$ is $\dot r = \omega\times(r-q)$ <br> evaluate at the space-frame origin $r=0$: $\dot r = \omega\times(-q) = -\omega\times q$ <br> the screw's linear part is exactly this origin velocity <br> $\boxed{v = -\omega\times q}$.
- Q: Build screw axis 2 of the planar 2R arm: revolute axis $+z$ through the point $q_2=(2,0,0)$. Compute $\mathcal{S}_2=(\omega,v)$. <br> <b>Formula:</b> $\omega=(0,0,1)$, $v=-\omega\times q_2$ <br> <i>Hint:</i> compute the cross product $(0,0,1)\times(2,0,0)$ first, then negate it. :: A: revolute so $\omega=(0,0,1)$ with $\|\omega\|=1$ <br> $(0,0,1)\times(2,0,0)=(0,2,0)$ <br> $v=-\omega\times q_2 = -(0,2,0)=(0,-2,0)$ <br> $\boxed{\mathcal{S}_2=(0,0,1,\;0,-2,0)}$.
- Q: Derive the space-frame PoE formula, step 1 of 2: from home pose $M=T(0)$, move only joint $n$ and write $T$. <br> <b>Formula:</b> $T=e^{[\mathcal{S}_n]\theta_n}\,M$ <br> <i>Hint:</i> a space-frame screw acts on its base side, so place the exponential to the left of $M$. :: A: at home the tool sits at $M$ <br> moving only joint $n$ screws the outboard chain about $\mathcal{S}_n$ by $\theta_n$, applied on the base side of $M$ <br> a space-frame screw left-multiplies what it acts on <br> $\boxed{T = e^{[\mathcal{S}_n]\theta_n}\,M}$.
- Q: Finish the PoE formula, step 2 of 2: add joints $n-1,\dots,1$, each a space-frame screw on the base side. Write the full $T(\theta)$. <br> <b>Formula:</b> $T(\theta)=e^{[\mathcal{S}_1]\theta_1}e^{[\mathcal{S}_2]\theta_2}\cdots e^{[\mathcal{S}_n]\theta_n}\,M$ <br> <i>Hint:</i> each earlier joint prepends its exponential (left-multiplies), so chain from $i=n$ down to $i=1$. :: A: each earlier joint $i$ screws everything outboard (which already includes $M$ and later exponentials) about $\mathcal{S}_i$ <br> a space-frame screw left-multiplies, so joint $i$ prepends $e^{[\mathcal{S}_i]\theta_i}$ <br> chaining from $i=n$ down to $1$ <br> $\boxed{T(\theta)=e^{[\mathcal{S}_1]\theta_1}e^{[\mathcal{S}_2]\theta_2}\cdots e^{[\mathcal{S}_n]\theta_n}\,M}$.
- Q: Derive Rodrigues' rotation formula by summing the matrix-exponential series for unit $\omega$. <br> <b>Formula:</b> $e^{[\omega]\theta}=\sum_k \tfrac{(\theta[\omega])^k}{k!}$, identity $[\omega]^3=-[\omega]$, target $e^{[\omega]\theta}=I+\sin\theta\,[\omega]+(1-\cos\theta)[\omega]^2$ <br> <i>Hint:</i> use $[\omega]^3=-[\omega]$ to collapse every power onto $I$, $[\omega]$, or $[\omega]^2$, then recognize the two trig series. :: A: $[\omega]^3=-[\omega]\Rightarrow$ even powers collapse to $\pm[\omega]^2$, odd to $\pm[\omega]$ <br> group the series: $\big(\theta-\tfrac{\theta^3}{3!}+\cdots\big)[\omega]+\big(\tfrac{\theta^2}{2!}-\tfrac{\theta^4}{4!}+\cdots\big)[\omega]^2$, plus $I$ <br> the brackets are $\sin\theta$ and $(1-\cos\theta)$ <br> $\boxed{e^{[\omega]\theta}=I+\sin\theta\,[\omega]+(1-\cos\theta)[\omega]^2}$.
- Q: Derive the prismatic twist exponential by exponentiating the twist with $\omega=0$. <br> <b>Formula:</b> $[\mathcal{S}]=\begin{bmatrix}[\omega]&v\\0&0\end{bmatrix}$ with $\omega=0$, target $e^{[\mathcal{S}]\theta}=\begin{bmatrix}I&v\theta\\0&1\end{bmatrix}$ <br> <i>Hint:</i> set $\omega=0$, show $[\mathcal{S}]^2=0$, so the series truncates after the linear term. :: A: with $\omega=0$, $[\mathcal{S}]=\begin{bmatrix}0&v\\0&0\end{bmatrix}$ is nilpotent: $[\mathcal{S}]^2=0$ <br> series truncates: $e^{[\mathcal{S}]\theta}=I+\theta[\mathcal{S}]$ <br> $\boxed{e^{[\mathcal{S}]\theta}=\begin{bmatrix}I & v\,\theta\\0&1\end{bmatrix}}$ a pure translation.
- Q: Worked example (2D, near reach). Planar 2R arm $L_1=2,L_2=1$ at $\theta_1=0,\theta_2=90^\circ$. Find the tip. <br> <b>Formula:</b> $x=L_1\cos\theta_1+L_2\cos(\theta_1+\theta_2)$, $y=L_1\sin\theta_1+L_2\sin(\theta_1+\theta_2)$ <br> <i>Hint:</i> compute the combined angle $\theta_1+\theta_2=90^\circ$ first, then plug each term in. :: A: combined angle $\theta_1+\theta_2=0+90^\circ=90^\circ$ <br> $x=2\cos0+1\cos90^\circ=2+0=2$ <br> $y=2\sin0+1\sin90^\circ=0+1=1$ <br> $\boxed{p=(2,1,0)}$.
- Q: Worked example (cell scale, metres). Planar 2R arm $L_1=0.3\,\text{m},L_2=0.2\,\text{m}$ at $\theta_1=\theta_2=90^\circ$. Find the tip. <br> <b>Formula:</b> $x=L_1\cos\theta_1+L_2\cos(\theta_1+\theta_2)$, $y=L_1\sin\theta_1+L_2\sin(\theta_1+\theta_2)$ <br> <i>Hint:</i> the combined angle is $180^\circ$, so $\cos180^\circ=-1$ and $\sin180^\circ=0$ in the second term. :: A: combined angle $\theta_1+\theta_2=180^\circ$ <br> $x=0.3\cos90^\circ+0.2\cos180^\circ=0-0.2=-0.2$ <br> $y=0.3\sin90^\circ+0.2\sin180^\circ=0.3+0=0.3$ <br> $\boxed{p=(-0.2,\,0.3,\,0)\ \text{m}}$ (folded back over the base).
- Q: Worked example (revolute + prismatic). Base revolute $\mathcal{S}_1=(0,0,1,\,0,0,0)$ then prismatic $\mathcal{S}_2=(0,0,0,\,1,0,0)$, $M=I$, at $\theta_1=90^\circ$, $d=0.5$. Find the tip. <br> <b>Formula:</b> $T=e^{[\mathcal{S}_1]\theta_1}e^{[\mathcal{S}_2]d}M$, prismatic block $e^{[\mathcal{S}_2]d}=\begin{bmatrix}I & v\,d\\0&1\end{bmatrix}$ with $v=(1,0,0)$ <br> <i>Hint:</i> the slide contributes $v\,d=(0.5,0,0)$, then the $90^\circ$ rotation about $+z$ sends $+x\to+y$. :: A: prismatic block contributes the slide $v\,d=(0.5,0,0)$ <br> $e^{[\mathcal{S}_1]90^\circ}$ rotates $+x\to+y$, sending the slide to $(0,0.5,0)$ <br> $T=e^{[\mathcal{S}_1]\theta_1}e^{[\mathcal{S}_2]d}M$ has this tip column <br> $\boxed{p=(0,\,0.5,\,0)}$ slide runs along the rotated $x$ axis.
- Q: Worked example (DH path, mid reach). DH planar 2R $a_1=2,a_2=1,d=\alpha=0$ at $\theta_1=30^\circ,\theta_2=40^\circ$. Find the tip from $A_1A_2$. <br> <b>Formula:</b> $x=a_1\cos\theta_1+a_2\cos(\theta_1+\theta_2)$, $y=a_1\sin\theta_1+a_2\sin(\theta_1+\theta_2)$ <br> <i>Hint:</i> the second term uses $\theta_1+\theta_2=70^\circ$; evaluate $\cos70^\circ\approx0.3420$, $\sin70^\circ\approx0.9397$. :: A: combined angle $\theta_1+\theta_2=70^\circ$ <br> $x=2\cos30^\circ+1\cos70^\circ=1.7321+0.3420=2.0741$ <br> $y=2\sin30^\circ+1\sin70^\circ=1.0000+0.9397=1.9397$ <br> $\boxed{p\approx(2.0741,\,1.9397,\,0)}$.
- Q: Test + verify (straight arm closes the loop). 2R arm $L_1=2,L_2=1$ at $\theta_1=45^\circ,\theta_2=0$. Find the tip, then check the reach. <br> <b>Formula:</b> $x=L_1\cos\theta_1+L_2\cos(\theta_1+\theta_2)$, $y=L_1\sin\theta_1+L_2\sin(\theta_1+\theta_2)$, check $r=\sqrt{x^2+y^2}$ <br> <i>Hint:</i> with $\theta_2=0$ the arm is straight, so a pass means $r=L_1+L_2=3$. :: A: $\theta_2=0$ so both terms share angle $45^\circ$ <br> $x=2\cos45^\circ+1\cos45^\circ=3\cos45^\circ=2.1213$, $y=3\sin45^\circ=2.1213$ <br> verify radius: $\sqrt{x^2+y^2}=\sqrt{2\cdot2.1213^2}=3.000=L_1+L_2$ <br> $\boxed{p=(2.1213,\,2.1213,\,0),\ r=3=L_1+L_2}$.
- Q: Test + verify (PoE and DH must agree). Worked 2R arm at $\theta_1=90^\circ,\theta_2=0$. Show the PoE tip, then cross-check with the DH planar formula. <br> <b>Formula:</b> PoE $T=e^{[\mathcal{S}_1]\pi/2}\,I\,M$ with $M$ tip at $(3,0,0)$; DH $x=a_1\cos\theta_1+a_2\cos(\theta_1+\theta_2)$, $y=a_1\sin\theta_1+a_2\sin(\theta_1+\theta_2)$ <br> <i>Hint:</i> a pass means both tip columns match exactly; compute each independently then compare. :: A: PoE: $T=e^{[\mathcal{S}_1]\pi/2}\,I\,M$ gives tip column $(0,3,0)$ <br> DH: $x=2\cos90^\circ+1\cos90^\circ=0$, $y=2\sin90^\circ+1\sin90^\circ=3$ <br> the two independent paths land identically <br> $\boxed{p_{\text{PoE}}=(0,3,0)=p_{\text{DH}}}$.
- Q: Test + verify (zero-angle limit). Show that the revolute twist exponential $\to I$ as $\theta\to0$, confirming "no motion at home". <br> <b>Formula:</b> $e^{[\mathcal{S}]\theta}=\begin{bmatrix}e^{[\omega]\theta} & G(\theta)v\\0&1\end{bmatrix}$, $G(\theta)=I\theta+(1-\cos\theta)[\omega]+(\theta-\sin\theta)[\omega]^2$ <br> <i>Hint:</i> a pass means both the rotation block $\to I$ and $G(\theta)v\to0$; substitute $\theta=0$ into each. :: A: $e^{[\omega]\theta}\to I$ since $\sin0=0,\,1-\cos0=0$ <br> $G(0)=I\cdot0+0\cdot[\omega]+0\cdot[\omega]^2=0$, so translation $G(0)v=0$ <br> $\boxed{e^{[\mathcal{S}]0}=I}$ no motion at zero angle, as required.
- Q: Explain why screw axis 1 of the 2R arm has zero linear part, $\mathcal{S}_1=(0,0,1,0,0,0)$. <br> <b>Formula:</b> $v_1=-\omega\times q_1$ with axis $+z$ through $q_1=(0,0,0)$ <br> <i>Hint:</i> substitute $q_1=0$ into the cross product. :: A: axis 1 is $+z$ through the origin $q_1=(0,0,0)$ <br> $v_1=-\omega\times q_1=-\omega\times 0=0$ <br> $\boxed{\mathcal{S}_1=(0,0,1,\;0,0,0)}$ the space origin lies on the axis, so it has no velocity.
- Q: Define the forward kinematics function of a serial chain. <br> <i>Hint:</i> it is a map from one set of numbers to a pose in $SE(3)$. :: A: the map from joint variables to end-effector pose <br> $\boxed{T(\theta)\in SE(3)}$.
- Q: In the PoE formula, what is $M$ and where does it sit? <br> <b>Formula:</b> $T(\theta)=e^{[\mathcal{S}_1]\theta_1}\cdots e^{[\mathcal{S}_n]\theta_n}\,M$ <br> <i>Hint:</i> set every $\theta_i=0$ and see which factor survives. :: A: $M=T(0)$, the pose of $\{b\}$ relative to $\{s\}$ at the zero configuration <br> it sits on the far right, and must be the measured machine home, not a guess <br> $\boxed{M=T(0)}$.
- Q: State when classical DH is ill-posed but PoE is not. <br> <i>Hint:</i> think about what the DH common normal needs in order to be defined. :: A: DH is ill-posed when consecutive joint axes are parallel (the common normal is undefined) <br> PoE has no parameterization singularity, only kinematic Jacobian singularities <br> $\boxed{\text{DH fails on parallel consecutive axes; PoE does not}}$.
- Q: How do you branch the twist exponential for a prismatic joint in code? <br> <b>Formula:</b> prismatic $e^{[\mathcal{S}]\theta}=\begin{bmatrix}I&v\theta\\0&1\end{bmatrix}$ <br> <i>Hint:</i> test the magnitude of $\omega$ against a small $\varepsilon$ before calling Rodrigues. :: A: test $\|\omega\|<\varepsilon$ <br> if true, set rotation to identity and translation to $v\theta$, skipping Rodrigues (which would divide by a zero angle) <br> $\boxed{\text{if }\|\omega\|<\varepsilon:\ R=I,\ p=v\theta}$.
- Q: After many product compositions, $T$'s rotation block drifts. What invariant breaks, and how do you restore it before IFC export? <br> <b>Formula:</b> orthonormality $R^\top R=I$ <br> <i>Hint:</i> a valid rotation has all singular values equal to 1; use the SVD to force that. :: A: the invariant that breaks is orthonormality: $R$ must stay in $SO(3)$ with $R^\top R=I$ <br> re-project via SVD, setting all singular values to 1 ($R\leftarrow UV^\top$) <br> else IFC rejects the non-orthonormal axis triple <br> $\boxed{\text{re-orthonormalize }R\text{ via SVD before export}}$.
- CLOZE: Forward kinematics is the ordered product {{c1::$T(\theta)=A_1A_2\cdots A_n$}}, exact and not an approximation.
- CLOZE: In the space-frame PoE product each exponential {{c1::left-multiplies}} and the home pose {{c2::$M$}} sits on the far right. <br> <i>hint:</i> set all $\theta_i=0$ to see which factor remains.
- CLOZE: The revolute closed-form exponential requires {{c1::$\|\omega\|=1$}} so that $\theta$ equals the physical rotation angle.
- CLOZE: For a revolute joint with unit axis $\omega$ through point $q$, the screw linear part is {{c1::$v=-\omega\times q$}}, the velocity of the space origin under unit rotation.
- CLOZE: A space-frame screw {{c1::left}}-multiplies and a body-frame screw {{c2::right}}-multiplies; never interleave the two conventions.
