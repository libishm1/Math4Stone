# The manipulator Jacobian, adjoints, and singularities

Prereqs: [[12_forward-kinematics]], [[09_calculus-jacobians]] => this => [[15_inverse-kinematics]]

## Why it matters for the stone scan->CAD->IFC / CADFit workflow

A scanned stone block does not move itself. A robot arm carries the scanner over it, then later carries a cutting or polishing tool back to the same surface. The Jacobian is the matrix that converts joint speeds into end-effector twist, so it is the exact object you differentiate to keep a tool moving at constant speed along a fitted B-spline edge of the stone. When you later solve inverse kinematics to reach a target pose on the CAD model (the next topic), the same Jacobian is the local linear model that the IK iteration steps along. The manipulability ellipsoid tells you which scan poses are dexterous and which are near-singular, so you can plan a scan path that never freezes mid-acquisition. The force-velocity duality of the Jacobian is what tells you whether the arm can actually push hard enough to drag a probe across a rough granite face in a given configuration.

## Intuition

Forward kinematics is a function: feed it joint angles, get back the end-effector pose. The Jacobian is the derivative of that function. It answers a velocity question: if I turn each joint a little bit, how does the gripper move and spin right now?

Picture a desk lamp with several hinges. Hold the bulb and wiggle one hinge. The bulb traces a little arc. The direction and speed of that arc depend on where all the other hinges currently sit. Each joint contributes one such "instantaneous motion direction" to the bulb. The Jacobian stacks those contributions as columns. Multiply the column vector of joint speeds by the Jacobian and you get the bulb's linear and angular velocity, combined into a single 6-vector called a twist.

The key twist of the story: column $i$ is the screw axis of joint $i$, but not in its home pose. It is that axis carried into its current location by all the joints that come before it. The tool that "carries a screw axis to a new frame" is the adjoint. So the Jacobian is just a list of adjoint-transported screw axes.

When two columns line up, the arm loses the ability to move in some direction. That is a singularity. The math sees it as the matrix dropping rank.

## The math (first principles)

### Conventions

All formulas use the Product of Exponentials (PoE) convention from Modern Robotics (Lynch and Park). Spatial twists and screw axes are stacked angular-part-first:

$$\mathcal{V} = \begin{bmatrix} \omega \\ v \end{bmatrix} \in \mathbb{R}^6, \qquad \mathcal{S}_i = \begin{bmatrix} \omega_i \\ v_i \end{bmatrix}.$$

Here $\omega \in \mathbb{R}^3$ is angular velocity in rad/s and $v \in \mathbb{R}^3$ is linear velocity in m/s. Frames are right-handed. A pose is $T = \begin{bmatrix} R & p \\ 0 & 1 \end{bmatrix} \in SE(3)$ with $R \in SO(3)$ and $p \in \mathbb{R}^3$ in meters. Units of $J$ columns are mixed: the angular rows are dimensionless-per-radian, the linear rows are meters-per-radian. This mixed unit matters for manipulability (see pitfalls).

### Screw axis and the skew operator

The skew (hat) of a 3-vector $\omega = (\omega_x, \omega_y, \omega_z)$ is

$$[\omega] = \begin{bmatrix} 0 & -\omega_z & \omega_y \\ \omega_z & 0 & -\omega_x \\ -\omega_y & \omega_x & 0 \end{bmatrix}, \qquad [\omega]\,a = \omega \times a.$$

A revolute joint with unit axis $\hat{s}$ through a point $q$ has screw axis $\mathcal{S} = (\hat{s},\, -\hat{s}\times q)$. A prismatic joint along $\hat{s}$ has $\mathcal{S} = (0,\, \hat{s})$.

### The adjoint

Given $T = (R, p) \in SE(3)$, the adjoint is the $6\times 6$ matrix that transports a twist from one frame to another:

$$[\mathrm{Ad}_T] = \begin{bmatrix} R & 0 \\ [p]\,R & R \end{bmatrix} \in \mathbb{R}^{6\times 6}.$$

Its defining property is $[\mathrm{Ad}_T]\,\mathcal{V} = (T [\mathcal{V}] T^{-1})^{\vee}$, that is, conjugation of a twist by a rigid motion. Useful identity: $[\mathrm{Ad}_T]^{-1} = [\mathrm{Ad}_{T^{-1}}]$.

### The space Jacobian

Forward kinematics in PoE form is

$$T(\theta) = e^{[\mathcal{S}_1]\theta_1}\, e^{[\mathcal{S}_2]\theta_2}\cdots e^{[\mathcal{S}_n]\theta_n}\, M,$$

where $M$ is the home pose. Differentiate and rearrange. The end-effector twist in the space frame is linear in the joint rates:

$$\boxed{\;\mathcal{V}_s = J_s(\theta)\,\dot{\theta}\;}$$

Column $i$ is the $i$-th screw axis pushed forward by the product of the preceding exponentials, via the adjoint:

$$\boxed{\;J_{s,i}(\theta) = \left[\mathrm{Ad}_{e^{[\mathcal{S}_1]\theta_1}\cdots\, e^{[\mathcal{S}_{i-1}]\theta_{i-1}}}\right]\mathcal{S}_i\;}$$

with $J_{s,1} = \mathcal{S}_1$ (empty product is identity). Each column is the instantaneous screw the end-effector inherits from joint $i$ in the current configuration.

Short derivation of column $i$: write $P_{i-1} = e^{[\mathcal{S}_1]\theta_1}\cdots e^{[\mathcal{S}_{i-1}]\theta_{i-1}}$. Only $\theta_i$ multiplies $\mathcal{S}_i$ inside $P_{i-1}\,e^{[\mathcal{S}_i]\theta_i}$, so $\partial T/\partial\theta_i$ contributes the twist $P_{i-1}[\mathcal{S}_i]P_{i-1}^{-1}$ in space coordinates. Taking the vee of that conjugation gives exactly $[\mathrm{Ad}_{P_{i-1}}]\mathcal{S}_i$.

### The body Jacobian

The same end-effector twist, expressed in the body (end-effector) frame, is

$$\mathcal{V}_b = J_b(\theta)\,\dot{\theta}, \qquad J_{b,i}(\theta) = \left[\mathrm{Ad}_{e^{-[\mathcal{B}_n]\theta_n}\cdots\, e^{-[\mathcal{B}_{i+1}]\theta_{i+1}}}\right]\mathcal{B}_i,$$

where $\mathcal{B}_i$ are screw axes expressed in the body frame at the home pose. Here column $i$ depends on the joints after $i$, with $J_{b,n} = \mathcal{B}_n$. The two Jacobians are related by the adjoint of the current pose $T_{sb} = T(\theta)$:

$$\boxed{\;J_s = [\mathrm{Ad}_{T_{sb}}]\,J_b, \qquad J_b = [\mathrm{Ad}_{T_{bs}}]\,J_s = [\mathrm{Ad}_{T_{sb}^{-1}}]\,J_s\;}$$

Use the space Jacobian when reasoning in the fixed/world frame (scan rig coordinates). Use the body Jacobian when the task is naturally tool-relative (feed rate along a tool axis).

### Velocity, force, and statics

The Jacobian is the velocity map $\mathcal{V} = J\dot\theta$. By the principle of virtual work, the transpose is the force map. A wrench $\mathcal{F}$ (a stacked moment-and-force 6-vector) applied at the end-effector requires joint torques

$$\boxed{\;\tau = J^{\top}(\theta)\,\mathcal{F}\;}$$

So one matrix governs both differential motion and static force. This is why the Jacobian is the central object of manipulator control.

### Manipulability and singularities

Take a unit ball of joint rates, $\|\dot\theta\| = 1$, and map it through $J$. The image is the manipulability ellipsoid. Its shape is governed by

$$A = J J^{\top}.$$

The principal axes of the ellipsoid are the eigenvectors of $A$, and the axis half-lengths are $\sqrt{\lambda_i(A)} = \sigma_i(J)$, the singular values of $J$. Three standard scalar measures (Yoshikawa) are:

- Isotropy ratio: $\mu_1 = \dfrac{\sqrt{\lambda_{\max}}}{\sqrt{\lambda_{\min}}} = \dfrac{\sigma_{\max}}{\sigma_{\min}} \ge 1$. This is the condition number of $J$. Equal to 1 means isotropic (equally easy to move in every direction).
- Condition number of $A$: $\mu_2 = \dfrac{\lambda_{\max}}{\lambda_{\min}} = \mu_1^2$.
- Volume measure: $\mu_3 = \sqrt{\det A} = \prod_i \sigma_i$, proportional to the ellipsoid volume.

The force ellipsoid is the dual. It is governed by $A^{-1} = (JJ^{\top})^{-1}$, so its axes are the reciprocals of the velocity axes. Directions where the arm moves fast are exactly the directions where it can push only weakly, and vice versa.

A singularity is a configuration where $J$ loses rank. For a non-redundant arm ($n = 6$), $\mathrm{rank}\,J < 6 \iff \det J = 0$. At a singularity $\sigma_{\min} \to 0$, so $\mu_3 \to 0$ (zero-volume, flattened ellipsoid), $\mu_1 \to \infty$ (infinite condition number), and the arm cannot generate end-effector velocity in at least one direction. IK and resolved-rate control blow up there because $J^{-1}$ or $J^{+}$ has unbounded entries.

## Worked example

Take a planar 2R arm (two revolute joints, link lengths $L_1 = L_2 = 1$). Drop the angular row redundancy and use the familiar planar position Jacobian (this is the $2\times 2$ task-space block; it is the right object to inspect singularities by hand).

End-effector position:

$$x = L_1 c_1 + L_2 c_{12}, \qquad y = L_1 s_1 + L_2 s_{12},$$

with $c_1 = \cos\theta_1$, $c_{12} = \cos(\theta_1+\theta_2)$, etc. Differentiate:

$$J = \begin{bmatrix} -L_1 s_1 - L_2 s_{12} & -L_2 s_{12} \\ \;\;L_1 c_1 + L_2 c_{12} & \;\;L_2 c_{12} \end{bmatrix}.$$

Evaluate at $\theta_1 = 0$, $\theta_2 = \pi/2$. Then $s_1=0,\,c_1=1$, and $\theta_1+\theta_2 = \pi/2$ so $s_{12}=1,\,c_{12}=0$:

$$J = \begin{bmatrix} 0 - 1 & -1 \\ 1 + 0 & 0 \end{bmatrix} = \begin{bmatrix} -1 & -1 \\ 1 & 0 \end{bmatrix}.$$

Determinant: $\det J = (-1)(0) - (-1)(1) = 1 \ne 0$. Non-singular. Now compute manipulability:

$$A = JJ^{\top} = \begin{bmatrix} -1 & -1 \\ 1 & 0 \end{bmatrix}\begin{bmatrix} -1 & 1 \\ -1 & 0 \end{bmatrix} = \begin{bmatrix} 2 & -1 \\ -1 & 1 \end{bmatrix}.$$

Eigenvalues solve $\lambda^2 - 3\lambda + 1 = 0$, giving $\lambda = \tfrac{3\pm\sqrt5}{2} \approx 2.618,\, 0.382$. So $\sigma_{\max} \approx 1.618$, $\sigma_{\min} \approx 0.618$. Measures: $\mu_3 = \sqrt{\det A} = \sqrt{2\cdot1 - 1} = 1$ (equals $|\det J|$, as it must in the square case), and $\mu_1 = 1.618/0.618 \approx 2.618$.

Now the singular case $\theta_2 = 0$ (arm fully stretched), $\theta_1 = 0$: $s_{12}=0,\,c_{12}=1$:

$$J = \begin{bmatrix} 0 & 0 \\ 2 & 1 \end{bmatrix}, \qquad \det J = 0.$$

Rank 1. The arm can only move along $y$ at this instant; it cannot produce any $x$-velocity no matter how you spin the joints. Both joints push the tip in the same direction, so the columns are parallel. That is the boundary singularity.

## Code

Idiomatic C++ with Eigen. The twist convention is $[\omega; v]$ (angular first), matching Modern Robotics. Symbols map directly: `Slist[i]` is $\mathcal{S}_i$, `theta[i]` is $\theta_i$, `expSE3` is $e^{[\mathcal{S}]\theta}$, `adjoint(T)` is $[\mathrm{Ad}_T]$.

```cpp
// jacobian.hpp  -- space/body Jacobian, manipulability, singularity check
#pragma once
#include <Eigen/Dense>
#include <vector>
#include <cmath>

namespace kin {

using Vec6 = Eigen::Matrix<double, 6, 1>;
using Mat4 = Eigen::Matrix4d;
using Mat6 = Eigen::Matrix<double, 6, 6>;

// [w] : skew-symmetric matrix so that [w]*a == w x a
inline Eigen::Matrix3d hat3(const Eigen::Vector3d& w) {
    Eigen::Matrix3d W;
    W <<    0.0, -w.z(),  w.y(),
          w.z(),    0.0, -w.x(),
         -w.y(),  w.x(),    0.0;
    return W;
}

// e^{[w]theta} via Rodrigues; w need not be unit (theta absorbs |w|)
inline Eigen::Matrix3d expSO3(const Eigen::Vector3d& w, double theta) {
    const double n = w.norm();
    if (n < 1e-12) return Eigen::Matrix3d::Identity();
    const Eigen::Matrix3d W = hat3(w / n);
    const double th = n * theta;
    return Eigen::Matrix3d::Identity()
         + std::sin(th) * W
         + (1.0 - std::cos(th)) * (W * W);
}

// e^{[S]theta} for a screw S = [w; v], handling the prismatic case w == 0
inline Mat4 expSE3(const Vec6& S, double theta) {
    const Eigen::Vector3d w = S.head<3>();
    const Eigen::Vector3d v = S.tail<3>();
    Mat4 T = Mat4::Identity();
    if (w.norm() < 1e-12) {                 // prismatic joint
        T.block<3,1>(0,3) = v * theta;
        return T;
    }
    const double n  = w.norm();
    const double th = n * theta;
    const Eigen::Matrix3d W = hat3(w / n);
    const Eigen::Matrix3d G = Eigen::Matrix3d::Identity() * th
                            + (1.0 - std::cos(th)) * W
                            + (th - std::sin(th)) * (W * W);
    T.block<3,3>(0,0) = expSO3(w, theta);
    T.block<3,1>(0,3) = (G * (v / n));      // scale v by 1/|w| to keep theta the path param
    return T;
}

// [Ad_T] = [[R, 0], [ [p]R, R ]]
inline Mat6 adjoint(const Mat4& T) {
    const Eigen::Matrix3d R = T.block<3,3>(0,0);
    const Eigen::Vector3d p = T.block<3,1>(0,3);
    Mat6 Ad = Mat6::Zero();
    Ad.block<3,3>(0,0) = R;
    Ad.block<3,3>(3,3) = R;
    Ad.block<3,3>(3,0) = hat3(p) * R;
    return Ad;
}

// Space Jacobian: column i = [Ad_{prod of preceding exp}] * S_i
// Mirrors Modern Robotics JacobianSpace.
inline Eigen::MatrixXd spaceJacobian(const std::vector<Vec6>& Slist,
                                     const std::vector<double>& theta) {
    const int n = static_cast<int>(Slist.size());
    Eigen::MatrixXd J(6, n);
    Mat4 T = Mat4::Identity();               // running product P_{i-1}
    for (int i = 0; i < n; ++i) {
        J.col(i) = (i == 0) ? Slist[i] : (adjoint(T) * Slist[i]);
        T = T * expSE3(Slist[i], theta[i]);  // fold in joint i for the next column
    }
    return J;
}

// Yoshikawa volume measure mu3 = sqrt(det(J J^T)) = product of singular values
inline double manipulabilityVolume(const Eigen::MatrixXd& J) {
    Eigen::MatrixXd A = J * J.transpose();   // 6x6 (or task-block) SPD-ish matrix
    return std::sqrt(std::max(0.0, A.determinant()));
}

// Condition number mu1 = sigma_max / sigma_min ; large => near singular
inline double jacobianConditionNumber(const Eigen::MatrixXd& J) {
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(J);   // singular values only
    const auto& s = svd.singularValues();
    const double smin = s(s.size() - 1);
    if (smin < 1e-12) return std::numeric_limits<double>::infinity();
    return s(0) / smin;                      // sigma_max / sigma_min
}

} // namespace kin
```

A 10-line check that ties the math together: at a singularity `manipulabilityVolume` returns ~0 and `jacobianConditionNumber` returns a huge number. Compute both along a scan path and refuse waypoints where the condition number exceeds a threshold (a common choice is 20 to 50 for the dimensionless angular block).

Python sketch (NumPy), only to show the same three lines compactly:

```python
import numpy as np
# Js : 6xn space Jacobian already built from adjoints
A = Js @ Js.T
mu3 = np.sqrt(max(0.0, np.linalg.det(A)))      # volume measure
sv  = np.linalg.svd(Js, compute_uv=False)      # singular values
cond = sv[0] / sv[-1] if sv[-1] > 1e-12 else np.inf
```

## Connections

- [[12_forward-kinematics]]: the Jacobian is the configuration-derivative of the PoE forward-kinematics map $T(\theta)$; you cannot build $J$ without the same exponentials.
- [[09_calculus-jacobians]]: this is the calculus Jacobian specialized to a map from joint space into $SE(3)$ velocities; the chain rule and partial-derivative-as-column idea carry over directly.
- [[15_inverse-kinematics]]: IK iterates $\Delta\theta = J^{+}e$ (or damped least squares near singularities); the Jacobian is the local linear model it steps along, and the condition number warns when the step is unreliable.
- [[03_rigid-transforms-se3]]: the adjoint $[\mathrm{Ad}_T]$ is built from $R$ and $p$ of an $SE(3)$ element; understanding twists and conjugation there is prerequisite to the column formula here.
- [[11_svd-eigendecomposition]]: the manipulability ellipsoid, singular values, and rank-loss test are pure SVD/eigendecomposition applied to $JJ^{\top}$.

## Pitfalls and failure modes

- Twist ordering. Modern Robotics stacks $[\omega; v]$ (angular first). Many other texts and ROS use $[v; \omega]$ (linear first). Mixing the two silently corrupts every adjoint and column. Pick one and assert it at every boundary.
- Space versus body confusion. $J_s$ and $J_b$ are different matrices for the same configuration, related by $[\mathrm{Ad}_{T_{sb}}]$. A velocity command written in the wrong frame sends the tool the wrong way. Tag every Jacobian with its frame.
- Mixed units in manipulability. The angular rows are rad-based and the linear rows are meter-based, so $JJ^{\top}$ mixes units and the ellipsoid volume is not frame-invariant. Compare manipulability only at a fixed length scale, or weight the rows by a characteristic length before taking $\det$.
- Near-singular conditioning. Do not test $\det J = 0$ exactly. Floating point never lands on zero. Test $\sigma_{\min}$ against a tolerance, or test the condition number against a ceiling. Near singularities, switch IK to damped least squares; a raw $J^{-1}$ produces enormous joint velocities.
- Wrist/representation singularities. Singularities are real geometric events (aligned axes), not just numerical noise. They depend on the kinematics, not the solver. Plan paths that route around them rather than patching the solver afterward.
- Redundant arms. For $n > 6$, $J$ is wide and $\det J$ is undefined; use $\sqrt{\det(JJ^{\top})}$ and the pseudoinverse, not the inverse.
- Non-commutativity. The column formula uses the product of preceding exponentials in order. Reordering joints changes the answer. $e^{[\mathcal{S}_1]\theta_1}e^{[\mathcal{S}_2]\theta_2} \ne e^{[\mathcal{S}_2]\theta_2}e^{[\mathcal{S}_1]\theta_1}$ in general.
- Forgetting $1/\|\omega\|$ in $\exp$. When scaling the translation part of $e^{[\mathcal{S}]\theta}$ you must divide $v$ by $\|\omega\|$ if $\mathcal{S}$ is not unit-normalized; otherwise the screw pitch is wrong and the column is silently off.

## Practice

1. Implement `bodyJacobian` from the formula $J_{b,i} = [\mathrm{Ad}_{\exp(-[\mathcal{B}_n]\theta_n)\cdots\exp(-[\mathcal{B}_{i+1}]\theta_{i+1})}]\mathcal{B}_i$. Verify numerically that $J_s = [\mathrm{Ad}_{T(\theta)}]\,J_b$ for a random 6R arm and a random $\theta$. Print the max element-wise difference; it should be near machine epsilon.
2. For the planar 2R arm above, sweep $\theta_2$ from $0$ to $\pi$ in 50 steps with $\theta_1 = 0$. Plot $\mu_3 = |\det J|$ versus $\theta_2$. Identify the two singularities (fully extended and fully folded) as the zeros of the curve.
3. Verify the Jacobian by finite differences. For a 6R arm, compute $J_s$ analytically, then compute column $i$ as $(\text{vec of } T(\theta + \epsilon e_i)T(\theta)^{-1})/\epsilon$ for small $\epsilon$. Confirm they agree to $O(\epsilon)$.
4. Draw the manipulability ellipse. At three configurations of the 2R arm, compute the eigenvectors and eigenvalues of $JJ^{\top}$ and render the ellipse axes. Show that the ellipse flattens to a line as you approach the singularity.
5. (Stretch) Add a condition-number gate to a resolved-rate controller: command a straight-line tool motion and refuse/clamp the step whenever $\sigma_{\max}/\sigma_{\min}$ exceeds 30. Compare joint-velocity peaks with and without the gate.
6. (Stretch) Show the force-velocity duality: at one configuration, compute the velocity ellipsoid from $JJ^{\top}$ and the force ellipsoid from $(JJ^{\top})^{-1}$ and confirm their principal axis lengths are reciprocal.

## References

- Kevin M. Lynch and Frank C. Park, *Modern Robotics: Mechanics, Planning, and Control*, Cambridge University Press, 2017. Chapter 5 (Velocity Kinematics and Statics). Full PDF: https://hades.mech.northwestern.edu/images/2/2a/Park-lynch.pdf
- Modern Robotics course resource, 5.1.2 Body Jacobian: https://modernrobotics.northwestern.edu/nu-gm-book-resource/5-1-2-body-jacobian/
- Modern Robotics course resource, 5.4 Manipulability: https://modernrobotics.northwestern.edu/nu-gm-book-resource/5-4-manipulability/
- Modern Robotics course resource, 3.3.2 Twists and the adjoint: https://modernrobotics.northwestern.edu/nu-gm-book-resource/3-3-2-twists-part-2-of-2/
- ModernRobotics reference code (NxRLab), `JacobianSpace`, `JacobianBody`, `Adjoint`: https://github.com/NxRLab/ModernRobotics
- T. Yoshikawa, "Manipulability of Robotic Mechanisms," *International Journal of Robotics Research*, 1985 (origin of the volume manipulability measure).
- Wikipedia, Manipulability ellipsoid (SVD axes, condition number, isotropy): https://en.wikipedia.org/wiki/Manipulability_ellipsoid
- Eigen documentation, `JacobiSVD` and dense decompositions: https://eigen.tuxfamily.org/dox/group__TutorialLinearAlgebra.html
- G. Nehme, E. Whalen, F. Ahmed, "CADFit: Precise Mesh-to-CAD Program Generation with Hybrid Optimization," arXiv:2605.01171, 2026: https://arxiv.org/abs/2605.01171 (motivates the downstream CAD-program-synthesis layer this pipeline feeds).

## Flashcard seeds

- Q: A planar 2R arm has tip position $x=L_1c_1+L_2c_{12}$, $y=L_1s_1+L_2s_{12}$. Find the position Jacobian $J=\partial(x,y)/\partial(\theta_1,\theta_2)$. <br> <b>Formula:</b> $J_{rc}=\partial(\text{row }r)/\partial(\text{col }\theta_c)$, with $\tfrac{d}{d\theta}\cos\theta=-\sin\theta$, and $c_{12}=\cos(\theta_1+\theta_2)$ depends on BOTH $\theta_1,\theta_2$. <br> <i>Hint:</i> differentiate $x$ first; remember $c_{12}$ contributes to both columns since $\theta_1+\theta_2$ contains each angle. :: A: $\partial x/\partial\theta_1=-L_1s_1-L_2s_{12}$, $\partial x/\partial\theta_2=-L_2s_{12}$ <br> $\partial y/\partial\theta_1=L_1c_1+L_2c_{12}$, $\partial y/\partial\theta_2=L_2c_{12}$ <br> $\boxed{J=\begin{bmatrix}-L_1s_1-L_2s_{12} & -L_2s_{12}\\ L_1c_1+L_2c_{12} & L_2c_{12}\end{bmatrix}}$.
- Q: Evaluate the general 2R Jacobian numerically at $L_1=L_2=1,\ \theta_1=0,\ \theta_2=\pi/2$. <br> <b>Formula:</b> $J=\begin{bmatrix}-L_1s_1-L_2s_{12} & -L_2s_{12}\\ L_1c_1+L_2c_{12} & L_2c_{12}\end{bmatrix}$. <br> <i>Hint:</i> first compute $s_1,c_1$ at $\theta_1=0$ and $s_{12},c_{12}$ at $\theta_1+\theta_2=\pi/2$, then substitute. :: A: $s_1=0,c_1=1$; $\theta_1+\theta_2=\pi/2$ so $s_{12}=1,c_{12}=0$ <br> top row $(0-1,\,-1)$, bottom row $(1+0,\,0)$ <br> $\boxed{J=\begin{bmatrix}-1 & -1\\ 1 & 0\end{bmatrix}}$.
- Q: Differentiate the PoE map $T=P_{i-1}e^{[\mathcal{S}_i]\theta_i}\cdots M$ w.r.t. $\theta_i$ to get the space-frame twist of column $i$. <br> <b>Formula:</b> $\frac{\partial}{\partial\theta_i}e^{[\mathcal{S}_i]\theta_i}=[\mathcal{S}_i]e^{[\mathcal{S}_i]\theta_i}$; with $P_{i-1}=e^{[\mathcal{S}_1]\theta_1}\cdots e^{[\mathcal{S}_{i-1}]\theta_{i-1}}$; target twist $=P_{i-1}[\mathcal{S}_i]P_{i-1}^{-1}$. <br> <i>Hint:</i> only the $i$-th factor depends on $\theta_i$; differentiate it, then form $(\partial T/\partial\theta_i)\,T^{-1}$ so the trailing factors cancel. :: A: $\partial T/\partial\theta_i=P_{i-1}[\mathcal{S}_i]e^{[\mathcal{S}_i]\theta_i}\cdots M$ <br> right-multiply by $T^{-1}$; everything from $e^{[\mathcal{S}_i]\theta_i}$ rightward cancels, leaving $P_{i-1}[\mathcal{S}_i]P_{i-1}^{-1}$ <br> $\boxed{(\partial T/\partial\theta_i)T^{-1}=P_{i-1}[\mathcal{S}_i]P_{i-1}^{-1}}$.
- Q: From the conjugation $P_{i-1}[\mathcal{S}_i]P_{i-1}^{-1}$, derive the space-Jacobian column formula. <br> <b>Formula:</b> adjoint identity $(T[\mathcal{V}]T^{-1})^\vee=[\mathrm{Ad}_T]\mathcal{V}$; target $J_{s,i}=[\mathrm{Ad}_{P_{i-1}}]\mathcal{S}_i$. <br> <i>Hint:</i> take the vee ($^\vee$) of the conjugation and read off the adjoint of $P_{i-1}$. :: A: conjugating a twist by a rigid motion is the adjoint action, $(T[\mathcal{V}]T^{-1})^\vee=[\mathrm{Ad}_T]\mathcal{V}$ <br> apply with $T=P_{i-1}$, $\mathcal{V}=\mathcal{S}_i$ <br> $\boxed{J_{s,i}=[\mathrm{Ad}_{P_{i-1}}]\,\mathcal{S}_i}$, with $J_{s,1}=\mathcal{S}_1$.
- Q: Map a unit ball of joint rates $\|\dot\theta\|=1$ through $J$ and find the matrix whose eigenvalues set the velocity-ellipsoid axes. <br> <b>Formula:</b> $\mathcal{V}=J\dot\theta$; on the row space $\dot\theta=J^+\mathcal{V}$; target quadric $\mathcal{V}^\top(JJ^\top)^{-1}\mathcal{V}=1$. <br> <i>Hint:</i> substitute $\dot\theta=J^+\mathcal{V}$ into $\dot\theta^\top\dot\theta=1$ and simplify $ (J^+)^\top J^+$. :: A: $\dot\theta^\top\dot\theta=1$ becomes $\mathcal{V}^\top(J^+)^\top J^+\mathcal{V}=\mathcal{V}^\top(JJ^\top)^{-1}\mathcal{V}=1$, a quadric <br> the shape matrix is the inverse of $JJ^\top$ <br> $\boxed{A=JJ^\top}$ governs the velocity ellipsoid.
- Q: From $A=JJ^\top$, express the ellipsoid axis half-lengths using the SVD $J=U\Sigma V^\top$. <br> <b>Formula:</b> $A=JJ^\top=U\Sigma^2U^\top$; eigenvalues $\lambda_i(A)=\sigma_i(J)^2$; half-length$_i=\sqrt{\lambda_i}$. <br> <i>Hint:</i> substitute $J=U\Sigma V^\top$ into $JJ^\top$ and use $V^\top V=I$ to collapse the middle. :: A: $JJ^\top=U\Sigma V^\top V\Sigma U^\top=U\Sigma^2U^\top$, so $\lambda_i=\sigma_i^2$ <br> axis half-length is $\sqrt{\lambda_i}$ <br> $\boxed{\text{half-length}_i=\sqrt{\lambda_i(JJ^\top)}=\sigma_i(J)}$.
- Q: Derive the joint torques needed to hold an end-effector wrench $\mathcal{F}$, starting from virtual work. <br> <b>Formula:</b> virtual-work balance $\tau^\top\delta\theta=\mathcal{F}^\top\delta\mathcal{X}$ for all $\delta\theta$; kinematics $\delta\mathcal{X}=J\,\delta\theta$; target $\tau=J^\top\mathcal{F}$. <br> <i>Hint:</i> substitute $\delta\mathcal{X}=J\,\delta\theta$, factor out $\delta\theta$, then equate the coefficients since it holds for every $\delta\theta$. :: A: $\tau^\top\delta\theta=\mathcal{F}^\top J\,\delta\theta$ for all $\delta\theta$ <br> so $\tau^\top=\mathcal{F}^\top J$; transpose both sides <br> $\boxed{\tau=J^\top(\theta)\,\mathcal{F}}$.
- Q: For the worked 2R $J=\begin{bmatrix}-1&-1\\1&0\end{bmatrix}$, compute $A=JJ^\top$ and its eigenvalues. <br> <b>Formula:</b> $A=JJ^\top$; eigenvalues solve $\lambda^2-(\mathrm{tr}\,A)\lambda+\det A=0$. <br> <i>Hint:</i> form $JJ^\top$ first (multiply $J$ by its transpose), then read off $\mathrm{tr}\,A$ and $\det A$ for the quadratic. :: A: $A=\begin{bmatrix}-1&-1\\1&0\end{bmatrix}\begin{bmatrix}-1&1\\-1&0\end{bmatrix}=\begin{bmatrix}2&-1\\-1&1\end{bmatrix}$ <br> $\mathrm{tr}=3,\ \det=2-1=1$, so $\lambda^2-3\lambda+1=0$ <br> $\boxed{\lambda=\tfrac{3\pm\sqrt5}{2}\approx 2.618,\ 0.382}$.
- Q: From eigenvalues $\lambda\approx 2.618,\ 0.382$, compute the singular values and the condition number $\mu_1$. <br> <b>Formula:</b> $\sigma_i=\sqrt{\lambda_i}$; $\mu_1=\sigma_{\max}/\sigma_{\min}$. <br> <i>Hint:</i> take square roots of the eigenvalues first, then divide the larger by the smaller. :: A: $\sigma_{\max}=\sqrt{2.618}\approx1.618,\ \sigma_{\min}=\sqrt{0.382}\approx0.618$ <br> $\mu_1=1.618/0.618$ <br> $\boxed{\mu_1\approx 2.618}$ (finite, so non-singular).
- Q: From the general 2R $J$, derive its determinant and the condition on $\theta_2$ that makes it singular. <br> <b>Formula:</b> $J=\begin{bmatrix}-L_1s_1-L_2s_{12} & -L_2s_{12}\\ L_1c_1+L_2c_{12} & L_2c_{12}\end{bmatrix}$; $\det\begin{bmatrix}a&b\\c&d\end{bmatrix}=ad-bc$; identity $\sin(B-A)=s_Bc_A-c_Bs_A$. <br> <i>Hint:</i> after $ad-bc$ the $L_2^2$ terms cancel; group the rest as $L_1L_2(s_{12}c_1-c_{12}s_1)$ and apply the sine-difference identity with $B=\theta_1+\theta_2,\ A=\theta_1$. :: A: $\det J=(-L_1s_1-L_2s_{12})(L_2c_{12})-(-L_2s_{12})(L_1c_1+L_2c_{12})$ <br> $L_2^2$ terms cancel: $=L_1L_2(s_{12}c_1-c_{12}s_1)=L_1L_2\sin\theta_2$ <br> $\boxed{\det J=L_1L_2\sin\theta_2=0\iff\theta_2=0\text{ or }\pi}$ (arm straight or folded).
- Q: At the fully stretched config $\theta_2=0,\ \theta_1=0,\ L_1=L_2=1$, build $J$ and show it is singular. <br> <b>Formula:</b> $J=\begin{bmatrix}-s_1-s_{12} & -s_{12}\\ c_1+c_{12} & c_{12}\end{bmatrix}$ (with $L_1=L_2=1$); singular $\iff\det J=0$. <br> <i>Hint:</i> evaluate $s_1,c_1$ at $\theta_1=0$ and $s_{12},c_{12}$ at $\theta_1+\theta_2=0$, then check the determinant. :: A: $s_1=0,c_1=1$; $\theta_1+\theta_2=0$ so $s_{12}=0,c_{12}=1$ <br> $J=\begin{bmatrix}0 & 0\\ 2 & 1\end{bmatrix}$, $\det J=0\cdot1-0\cdot2=0$ <br> $\boxed{\text{rank }1:\text{ tip moves only along }y,\text{ no }x\text{-velocity}}$.
- Q: Worked example (2D, near-folded). For the 2R arm at $\theta_2=0.1$ rad, $L_1=L_2=1$, compute $\mu_3=|\det J|$. <br> <b>Formula:</b> $\mu_3=|\det J|=L_1L_2|\sin\theta_2|$. <br> <i>Hint:</i> plug the lengths and angle straight in; use $\sin(0.1)\approx0.0998$. :: A: $\mu_3=1\cdot1\cdot|\sin(0.1)|$ <br> $\sin(0.1)\approx0.0998$ <br> $\boxed{\mu_3\approx 0.0998}$ (small, nearly singular near full extension).
- Q: Worked example (2D, dexterous). For the 2R arm at $\theta_2=\pi/2$, $L_1=L_2=1$, compute $\mu_3=|\det J|$ and compare to the near-folded $0.0998$ case. <br> <b>Formula:</b> $\mu_3=L_1L_2|\sin\theta_2|$. <br> <i>Hint:</i> $|\sin\theta_2|$ peaks at $\theta_2=\pi/2$; evaluate, then divide by $0.0998$ for the ratio. :: A: $\mu_3=1\cdot1\cdot|\sin(\pi/2)|=1$ <br> the right-angle elbow maximizes $|\sin\theta_2|$; ratio $1/0.0998\approx10$ <br> $\boxed{\mu_3=1,\ \approx10\times\text{ the near-folded value}}$.
- Q: Worked example (3D, mm scale). A revolute joint has unit axis $\hat s=(0,0,1)$ through point $q=(50,0,0)$ mm. Compute its screw axis. <br> <b>Formula:</b> $\mathcal{S}=(\hat s,\ -\hat s\times q)$, with $\hat s\times q=(s_yq_z-s_zq_y,\ s_zq_x-s_xq_z,\ s_xq_y-s_yq_x)$. <br> <i>Hint:</i> compute the cross product $\hat s\times q$ first, then negate it for the linear part. :: A: $\hat s\times q=(0\cdot0-1\cdot0,\ 1\cdot50-0\cdot0,\ 0\cdot0-0\cdot0)=(0,50,0)$ <br> $-\hat s\times q=(0,-50,0)$ <br> $\boxed{\mathcal{S}=(0,0,1,\ 0,-50,0)}$ (linear part in mm).
- Q: Worked example (statics, m scale). At a config with $J^\top=\begin{bmatrix}1&0\\0&2\end{bmatrix}$ (m/rad), an end-effector wrench $\mathcal{F}=(10,5)$ N is applied. Compute the joint torques. <br> <b>Formula:</b> $\tau=J^\top\mathcal{F}$. <br> <i>Hint:</i> it is a $2\times2$ times $2\times1$ matrix-vector product; multiply each row of $J^\top$ by $\mathcal{F}$. :: A: $\tau=\begin{bmatrix}1&0\\0&2\end{bmatrix}\begin{bmatrix}10\\5\end{bmatrix}$ <br> $=(1\cdot10+0\cdot5,\ 0\cdot10+2\cdot5)$ <br> $\boxed{\tau=(10,\ 10)\text{ N·m}}$.
- Q: Closing test (square-case identity). For the 2R arm at $\theta_2=\pi/2$ you found $J=\begin{bmatrix}-1&-1\\1&0\end{bmatrix}$ and eigenvalues $2.618,0.382$ of $A=JJ^\top$. Verify $\sqrt{\det A}=|\det J|$. <br> <b>Formula:</b> $\det A=\lambda_{\max}\lambda_{\min}$; $\mu_3=\sqrt{\det A}$; for square $J$ this must equal $|\det J|$. <br> <i>Hint:</i> a pass means both numbers land on $1$; compute $\det A$ from the eigenvalue product and $\det J$ from $ad-bc$. :: A: $\det A=2.618\times0.382\approx1.000$, so $\sqrt{\det A}=1$ <br> $|\det J|=|(-1)(0)-(-1)(1)|=1$ <br> $\boxed{\sqrt{\det(JJ^\top)}=|\det J|=1\ \checkmark}$.
- Q: Closing test (force-velocity duality). The velocity-ellipsoid axes are $\sigma_{\max}\approx1.618,\ \sigma_{\min}\approx0.618$. Verify the force-ellipsoid axes are their reciprocals. <br> <b>Formula:</b> velocity ellipsoid $\sim A=JJ^\top$; force ellipsoid $\sim A^{-1}$; force axis half-lengths $=1/\sigma_i$. <br> <i>Hint:</i> a pass means the fast-motion axis matches the weak-force axis; take $1/\sigma_{\max}$ and $1/\sigma_{\min}$ and see they swap. :: A: force ellipsoid governed by $A^{-1}$, eigenvalues $1/\lambda_i$, axes $1/\sigma_i$ <br> $1/1.618\approx0.618$ and $1/0.618\approx1.618$ <br> $\boxed{\text{fast-motion axis}=\text{weak-force axis},\ 1.618\leftrightarrow0.618\ \checkmark}$.
- Q: Closing test (singularity limit). Take $\det J=L_1L_2\sin\theta_2$ and let $\theta_2\to0$. Verify $\mu_3$, $\mu_1$, and rank all flag a singularity. <br> <b>Formula:</b> $\mu_3=\sqrt{\det A}=|\det J|$; $\mu_1=\sigma_{\max}/\sigma_{\min}$; singular $\iff\sigma_{\min}\to0$. <br> <i>Hint:</i> a pass means $\mu_3\to0$, $\mu_1\to\infty$, rank drops; track what $\det J\to0$ does to $\sigma_{\min}$. :: A: $\det J\to0\Rightarrow\sigma_{\min}\to0$ <br> $\mu_3=|\det J|\to0$, $\mu_1=\sigma_{\max}/\sigma_{\min}\to\infty$, rank drops $2\to1$ <br> $\boxed{\theta_2=0\text{ is a boundary singularity; all three measures confirm it}}$.
- Q: What does the manipulator Jacobian map, and in which direction? <br> <i>Hint:</i> think of it as the derivative of forward kinematics. :: A: it maps joint velocities to the end-effector twist, $\mathcal{V}=J(\theta)\dot\theta$ (a velocity / derivative map).
- Q: When do you use the space Jacobian vs the body Jacobian? <br> <b>Formula:</b> $J_s=[\mathrm{Ad}_{T_{sb}}]\,J_b$. <br> <i>Hint:</i> match the frame to where the task is described — fixed world vs the tool itself. :: A: space $J_s$ for fixed/world-frame reasoning (scan-rig coords); body $J_b$ for tool-relative tasks (feed rate along the tool axis) <br> they convert via the adjoint of the current pose, $J_s=[\mathrm{Ad}_{T_{sb}}]J_b$.
- Q: Pitfall: why should you never test $\det J=0$ exactly to detect a singularity? <br> <i>Hint:</i> think about floating-point arithmetic and what to test instead. :: A: floating point never lands exactly on zero <br> instead test $\sigma_{\min}$ against a tolerance or cap the condition number $\mu_1$, and switch IK to damped least squares near singularities.
- Q: Pitfall: why is the manipulability-ellipsoid volume not frame-invariant? <br> <i>Hint:</i> look at the units of the angular vs linear rows of $J$. :: A: $J$'s angular rows are rad-based and its linear rows are meter-based, so $JJ^\top$ mixes units <br> compare manipulability only at a fixed length scale, or weight the rows by a characteristic length before $\det$.
- CLOZE: In the space Jacobian, column $i$ is the screw axis $\mathcal{S}_i$ transported by the {{c1::adjoint}} of the product of the preceding exponentials.
- CLOZE: The space and body Jacobians satisfy $J_s=[\mathrm{Ad}_{{{c1::T_{sb}}}}]\,J_b$, the adjoint of the current pose. <br> <i>hint:</i> it is the pose carrying body-frame into space-frame.
- CLOZE: A configuration is singular when the Jacobian {{c1::loses rank}}; for a square non-redundant arm this is $\det J=0$.
- CLOZE: The Modern Robotics twist convention stacks {{c1::angular velocity}} first, $\mathcal{V}=[\omega;v]$, while ROS uses $[v;\omega]$.
- CLOZE: For a redundant arm ($n>6$) the manipulability measure is {{c1::$\sqrt{\det(JJ^\top)}$}} because $\det J$ is undefined for a non-square $J$.
- CLOZE: Velocity and force ellipsoids are {{c1::dual}}: directions of fast motion are directions of weak force, and vice versa.
