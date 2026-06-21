# Inverse kinematics: pseudoinverse and damped least squares

Prereqs: [[14_manipulator-jacobian]], [[08_nonlinear-least-squares]] => this => (none)

## Why it matters for the stone scan->CAD->IFC / CADFit workflow

In the stone pipeline a robot arm carries the scanner over a block, then later carries a tool back to a target pose on the fitted CAD surface. Forward kinematics gives you the tool pose from joint angles. Inverse kinematics is the reverse and the one you actually need at runtime: given a target frame on the B-spline face (a polishing start point, a drill entry, a contact-probe pose), find joint angles that put the tool there. That map is nonlinear and often has several valid arm postures, so there is no clean formula. You iterate with the Jacobian. Near a stretched-out or folded arm configuration the Jacobian goes singular and a naive solver explodes, which on real granite means a violent unplanned motion. Damped least squares is the one-line fix that keeps the step bounded and the arm safe, and it is the same Levenberg-Marquardt damping you already met fitting noisy scan data.

## Intuition

Forward kinematics is a known function from joints to pose. Inverse kinematics asks for its inverse, and the inverse has no formula. So you guess, measure how wrong you are, and step.

Think of reaching for a coffee cup in the dark. You move your hand, feel the gap between hand and cup, and move again to close that gap. Each correction is roughly "move in the direction that shrinks the error." The Jacobian is what converts a desired hand motion into joint motions. So one IK step is: compute the pose error, ask the Jacobian which joint moves shrink it, take that step, repeat.

The trouble starts when the arm is nearly straight or nearly folded. Then some hand directions cost almost infinite joint motion. The plain Jacobian-inverse step asks for those infinite motions and blows up. Damping says "do not trust the cheap directions too much," trading a little accuracy for a step that never explodes. The damping knob $\lambda$ is exactly that distrust dial.

## The math (first principles)

### Conventions

Same Product of Exponentials convention as [[14_manipulator-jacobian]] (Modern Robotics, Lynch and Park). Twists stack angular-part-first:

$$\mathcal{V} = \begin{bmatrix} \omega \\ v \end{bmatrix} \in \mathbb{R}^6, \qquad \omega \text{ in rad/s}, \quad v \text{ in m/s}.$$

A pose is $T = \begin{bmatrix} R & p \\ 0 & 1\end{bmatrix} \in SE(3)$, $R \in SO(3)$, $p \in \mathbb{R}^3$ in meters, frames right-handed. The Jacobian $J(\theta)$ maps joint rates to a twist, $\mathcal{V} = J(\theta)\dot\theta$. We use the body Jacobian $J_b$ (twist in the end-effector frame) for the full-pose solver, matching Modern Robotics. Joint angles $\theta \in \mathbb{R}^n$. Task dimension $m$ (here $m=6$ for a full pose, or $m=3$ for position only).

### The IK problem and its linearization

We want $\theta$ with $T(\theta) = T_{sd}$ (desired pose). Define a task-space error and linearize forward kinematics about the current guess $\theta_i$:

$$e \approx J(\theta_i)\,\Delta\theta, \qquad \Delta\theta = \theta_{i+1} - \theta_i.$$

This is one Gauss-Newton-style step on the nonlinear least-squares problem $\min_\theta \tfrac12\|e(\theta)\|^2$ from [[08_nonlinear-least-squares]]. The whole game is: choose how to solve $e = J\,\Delta\theta$ for $\Delta\theta$ when $J$ is not square or not full rank.

### Defining the pose error correctly

For a full $SE(3)$ target you cannot subtract poses. Use the matrix logarithm to turn the relative pose into a body twist. In the body frame Modern Robotics writes

$$[\mathcal{V}_b] = \log\!\big(T_{sb}^{-1}(\theta_i)\,T_{sd}\big),$$

so $\mathcal{V}_b = (\omega_b, v_b)$ is the constant-twist motion that would carry the current frame onto the goal frame in unit time. This $\mathcal{V}_b$ plays the role of $e$. Using $\log$ (not raw matrix subtraction) is what makes the rotational part correct and frame-consistent.

### Pseudoinverse step (Moore-Penrose)

Replace the inverse with the pseudoinverse $J^\dagger$:

$$\Delta\theta = J^\dagger e, \qquad \theta_{i+1} = \theta_i + J_b^\dagger(\theta_i)\,\mathcal{V}_b.$$

For a **fat** Jacobian (redundant robot, $n > m$, $J$ full row rank), Modern Robotics gives the right inverse

$$J^\dagger = J^T (J J^T)^{-1}, \qquad J J^\dagger = I.$$

Among all exact solutions of $J\,\Delta\theta = e$, this picks the **minimum-norm** $\Delta\theta$. For a **tall** $J$ ($n < m$) or a rank-deficient $J$ it instead returns the least-squares solution minimizing $\|J\,\Delta\theta - e\|$. The catch: as $J$ approaches a singularity, $JJ^T$ becomes nearly singular and $\|\Delta\theta\|$ blows up.

### Why singularities break the pseudoinverse (SVD view)

Write $J = U\Sigma V^T$ with singular values $\sigma_1 \ge \dots \ge \sigma_m \ge 0$ (see [[03_eig-svd-pca]]). Then

$$J^\dagger = V\,\Sigma^\dagger\,U^T, \qquad \Sigma^\dagger \text{ inverts each nonzero } \sigma_i \text{ as } 1/\sigma_i.$$

The step in SVD coordinates is $\Delta\theta = \sum_i \frac{1}{\sigma_i}\,(u_i^T e)\,v_i$. When $\sigma_i \to 0$, the factor $1/\sigma_i \to \infty$. A tiny error component along a nearly-lost direction demands a huge joint motion. That is the explosion.

### Damped least squares (the fix)

Instead of solving $e = J\,\Delta\theta$ exactly, minimize a regularized objective that also penalizes large steps:

$$\Delta\theta = \arg\min_{\Delta\theta} \; \|J\,\Delta\theta - e\|^2 + \lambda^2 \|\Delta\theta\|^2.$$

Setting the gradient to zero gives the normal equations $(J^T J + \lambda^2 I)\,\Delta\theta = J^T e$. Two algebraically equal closed forms result, and the second is the one quoted in the report and almost always cheaper because $m \le n$:

$$\Delta\theta = (J^T J + \lambda^2 I)^{-1} J^T e = J^T (J J^T + \lambda^2 I)^{-1} e.$$

The identity $(J^T J + \lambda^2 I)^{-1} J^T = J^T (J J^T + \lambda^2 I)^{-1}$ holds for any $\lambda > 0$. The right form inverts an $m \times m$ matrix; for a 6-DoF arm that is a fixed $6\times 6$ or $3\times 3$ solve regardless of joint count.

In SVD coordinates the damping replaces the reciprocal with a bounded filter:

$$\frac{1}{\sigma_i} \;\longrightarrow\; \frac{\sigma_i}{\sigma_i^2 + \lambda^2}.$$

For $\sigma_i \gg \lambda$ this is $\approx 1/\sigma_i$ (full accuracy in well-conditioned directions). For $\sigma_i \to 0$ this is $\approx \sigma_i/\lambda^2 \to 0$ (the dangerous direction is smoothly switched off). The maximum of $\sigma/(\sigma^2+\lambda^2)$ is $1/(2\lambda)$ at $\sigma = \lambda$, so **no** gain ever exceeds $1/(2\lambda)$. That single bound is why the step stays finite at and through singularities.

### Relation to Levenberg-Marquardt

This is exactly Levenberg-Marquardt applied to the IK residual. LM solves $\min_x \|r(x)\|^2$ with the step $(J^T J + \lambda^2 I)\,\Delta x = -J^T r$, interpolating between Gauss-Newton ($\lambda \to 0$) and gradient descent ($\lambda \to \infty$, where $\Delta\theta \to \tfrac{1}{\lambda^2}J^T e$, the Jacobian-transpose method up to a scale). DLS is the constant-$\lambda$, fixed-Jacobian special case used in robotics. Adaptive-$\lambda$ LM (shrink $\lambda$ on success, grow it on failure) is the production upgrade.

### Joint limits

The bare step ignores joint limits. Two cheap, real-world remedies:

1. **Clamp and re-solve.** After the step, clamp $\theta$ to $[\theta_{\min}, \theta_{\max}]$. If a joint is clamped, optionally drop its column from $J$ and re-solve so other joints take up the slack.
2. **Null-space projection (redundant arms only).** For $n > m$, any $\Delta\theta_0 = J^\dagger e$ can absorb a secondary motion in the null space without changing the task: $\Delta\theta = J^\dagger e + (I - J^\dagger J)\,z$. Choose $z = -k\,\nabla H(\theta)$ where $H$ penalizes nearness to limits. This steers away from limits while still hitting the target.

## Worked example

Planar 2R robot, both links 1 m, position only ($m=2$, $n=2$). This mirrors Example 6.1 in Modern Robotics. Forward kinematics of the tip:

$$x = \cos\theta_1 + \cos(\theta_1+\theta_2), \qquad y = \sin\theta_1 + \sin(\theta_1+\theta_2).$$

Position Jacobian (derivative of $(x,y)$ w.r.t. $(\theta_1,\theta_2)$), writing $s_1=\sin\theta_1$, $s_{12}=\sin(\theta_1+\theta_2)$, etc.:

$$J = \begin{bmatrix} -s_1 - s_{12} & -s_{12} \\ \;\;c_1 + c_{12} & \;\;c_{12} \end{bmatrix}.$$

Target tip $(x,y) = (0.366, 1.366)$, which is reached at $\theta_d = (30^\circ, 90^\circ)$. Start at $\theta_0 = (0^\circ, 30^\circ)$.

**Step at $\theta_0$.** With $\theta_1=0, \theta_2=30^\circ$: $\theta_1+\theta_2=30^\circ$, so $s_1=0,\ c_1=1,\ s_{12}=0.5,\ c_{12}=0.866$.

Current tip: $x = 1 + 0.866 = 1.866$, $y = 0 + 0.5 = 0.500$. Error $e = (0.366-1.866,\ 1.366-0.500) = (-1.500,\ 0.866)$.

$$J = \begin{bmatrix} -0.5 & -0.5 \\ 1.866 & 0.866 \end{bmatrix}, \qquad \det J = (-0.5)(0.866) - (-0.5)(1.866) = 0.500.$$

This $J$ is square and well away from singular, so the plain inverse already works here. $J^{-1} = \frac{1}{0.5}\begin{bmatrix} 0.866 & 0.5 \\ -1.866 & -0.5\end{bmatrix} = \begin{bmatrix} 1.732 & 1.0 \\ -3.732 & -1.0 \end{bmatrix}$.

$$\Delta\theta = J^{-1} e = \begin{bmatrix} 1.732(-1.5) + 1.0(0.866) \\ -3.732(-1.5) + (-1.0)(0.866) \end{bmatrix} = \begin{bmatrix} -1.732 \\ 4.732 \end{bmatrix}\ \text{rad}.$$

That is a $-99.2^\circ,\ +271^\circ$ step: a huge overshoot, because the linear model is only local. New guess (before any clamping) $\theta_1 \approx -99^\circ$. This is exactly why you want damping or a step-length limit. Now redo the step with **DLS**, $\lambda = 0.5$:

$$J J^T = \begin{bmatrix} 0.5 & -1.366 \\ -1.366 & 4.234 \end{bmatrix}, \quad J J^T + \lambda^2 I = \begin{bmatrix} 0.75 & -1.366 \\ -1.366 & 4.484 \end{bmatrix}.$$

$\det = 0.75(4.484) - (-1.366)^2 = 3.363 - 1.866 = 1.497$. Solve $(JJ^T+\lambda^2 I)\,y = e$:

$$y = \frac{1}{1.497}\begin{bmatrix} 4.484 & 1.366 \\ 1.366 & 0.75 \end{bmatrix}\begin{bmatrix} -1.5 \\ 0.866 \end{bmatrix} = \frac{1}{1.497}\begin{bmatrix} -5.543 \\ -1.399 \end{bmatrix} = \begin{bmatrix} -3.703 \\ -0.935 \end{bmatrix}.$$

$$\Delta\theta = J^T y = \begin{bmatrix} -0.5 & 1.866 \\ -0.5 & 0.866 \end{bmatrix}\begin{bmatrix} -3.703 \\ -0.935 \end{bmatrix} = \begin{bmatrix} 0.107 \\ 1.041 \end{bmatrix}\ \text{rad}.$$

The damped step is $+6.1^\circ,\ +59.6^\circ$: a much gentler, safer move in the right direction. Modern Robotics, using the full-pose body-twist solver (undamped), reports the first iterate $\theta_1 = (34.23^\circ, 79.18^\circ)$ and convergence to $(30.00^\circ, 90.00^\circ)$ within tolerance after 3 iterations. DLS trades a few extra iterations for never producing the wild $271^\circ$ swing the raw inverse asked for.

## Code

Idiomatic C++ with Eigen. This is the `ikStepPositionDLS` idiom from the report, extended to a full damped step plus a loop. Symbols map to the math: `J` is $J$, `lambda` is $\lambda$, `e` is the error $e$, `dtheta` is $\Delta\theta$.

```cpp
// ik_dls.hpp  -- damped least squares IK step and loop (Eigen).
#pragma once
#include <Eigen/Dense>
#include <vector>
#include <cmath>

namespace ik {

// One damped-least-squares step: dtheta = J^T (J J^T + lambda^2 I)^-1 e.
// J is m x n (m = task dim, n = joints). e is the m-vector task error.
// Uses the right form so the solve is always m x m (cheap for m <= 6).
inline Eigen::VectorXd dlsStep(const Eigen::MatrixXd& J,
                               const Eigen::VectorXd& e,
                               double lambda)
{
    const int m = static_cast<int>(J.rows());
    // A = J J^T + lambda^2 I   (the m x m damped Gram matrix)
    Eigen::MatrixXd A = J * J.transpose()
                      + lambda * lambda * Eigen::MatrixXd::Identity(m, m);
    // y = A^-1 e via Cholesky (A is SPD for lambda > 0), then dtheta = J^T y.
    Eigen::VectorXd y = A.ldlt().solve(e);
    return J.transpose() * y;             // dtheta = J^T y
}

// A problem supplies forward kinematics, the Jacobian, and the error.
struct IkProblem {
    // pose error e(theta) in task space (e.g. body twist from log(Tsb^-1 Tsd)).
    std::function<Eigen::VectorXd(const Eigen::VectorXd&)> error;
    // task Jacobian J(theta), m x n.
    std::function<Eigen::MatrixXd(const Eigen::VectorXd&)> jacobian;
    Eigen::VectorXd lower, upper;         // joint limits (empty => ignore).
};

struct IkResult { Eigen::VectorXd theta; int iters; bool converged; };

// DLS IK loop. tol on ||e||; max_iters cap; lambda is the damping.
inline IkResult solveDLS(const IkProblem& p, Eigen::VectorXd theta,
                         double lambda = 0.05, double tol = 1e-4,
                         int max_iters = 100)
{
    for (int i = 0; i < max_iters; ++i) {
        Eigen::VectorXd e = p.error(theta);
        if (e.norm() < tol) return {theta, i, true};
        Eigen::MatrixXd J = p.jacobian(theta);
        theta += dlsStep(J, e, lambda);   // theta_{i+1} = theta_i + dtheta
        // Joint-limit clamp: keep theta feasible after each step.
        if (p.lower.size() == theta.size())
            theta = theta.cwiseMax(p.lower).cwiseMin(p.upper);
    }
    return {theta, max_iters, p.error(theta).norm() < tol};
}

} // namespace ik
```

Plain Moore-Penrose, for comparison and for when you are safely away from singularities. Prefer the SVD form for numerical robustness:

```cpp
// Minimum-norm / least-squares pseudoinverse step via SVD: dtheta = J^+ e.
// Small singular values are truncated below eps to avoid blow-up.
inline Eigen::VectorXd pinvStep(const Eigen::MatrixXd& J,
                                const Eigen::VectorXd& e, double eps = 1e-6)
{
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(J, Eigen::ComputeThinU | Eigen::ComputeThinV);
    const auto& s = svd.singularValues();
    Eigen::VectorXd inv = s;
    for (int i = 0; i < s.size(); ++i)
        inv(i) = (s(i) > eps) ? 1.0 / s(i) : 0.0;   // 1/sigma, truncated
    // J^+ = V diag(inv) U^T ; dtheta = J^+ e.
    return svd.matrixV() * inv.asDiagonal() * (svd.matrixU().transpose() * e);
}
```

For the SE(3) error itself, build the body twist from the matrix log: `Vb = vee(MatrixLog6(Tsb.inverse() * Tsd))`, where `MatrixLog6` is the $SE(3)$ log from [[12_twists-screws-exp]]. That `Vb` is the `e` you pass to `dlsStep`.

## Connections

- [[14_manipulator-jacobian]]: supplies $J(\theta)$, the local linear model every IK step inverts; its singularities are exactly where DLS earns its keep.
- [[08_nonlinear-least-squares]]: IK is nonlinear least squares; the pseudoinverse step is Gauss-Newton and the damped step is Levenberg-Marquardt.
- [[13_forward-kinematics]]: provides $T(\theta)$ used to compute the pose error $T_{sb}^{-1}T_{sd}$ each iteration.
- [[12_twists-screws-exp]]: the matrix log $\log(T_{sb}^{-1}T_{sd})$ that turns a relative pose into the body twist error $\mathcal{V}_b$.
- [[05_transforms-se3]]: the $SE(3)$ algebra for poses, inverses, and frame composition the error depends on.
- [[03_eig-svd-pca]]: the SVD that explains why $1/\sigma$ blows up and why $\sigma/(\sigma^2+\lambda^2)$ does not.
- [[06_linear-least-squares]]: the normal equations $(J^TJ + \lambda^2 I)\Delta\theta = J^Te$ are damped (Tikhonov-regularized) least squares.

## Pitfalls and failure modes

- **Frame inconsistency.** Mixing a space-frame error with a body Jacobian (or vice versa) silently diverges. Match the error frame to the Jacobian: body twist $\mathcal{V}_b$ with $J_b$, spatial twist $\mathcal{V}_s$ with $J_s$.
- **Subtracting poses.** Never form rotational error by subtracting rotation matrices or naive Euler angles. Use $\log(T_{sb}^{-1}T_{sd})$. Euler subtraction breaks near wrap-around and is non-commutative.
- **Twist ordering.** Modern Robotics stacks $(\omega, v)$ angular-first. Many libraries stack $(v, \omega)$. A swapped convention scrambles $J$ and the error. Fix one convention repo-wide.
- **Undamped near singularities.** Plain $J^\dagger$ produces near-infinite $\Delta\theta$ when $\sigma_i \to 0$. Always damp, or cap the step length, when the arm may pass through stretched or folded poses.
- **Picking $\lambda$.** Too small gives no protection; too large gives sluggish, biased convergence (you never quite reach the target because every direction is attenuated). Make $\lambda$ adaptive: small when well-conditioned, larger as $\sigma_{\min}$ drops. A common rule is $\lambda^2 = \lambda_{\max}^2(1 - (\sigma_{\min}/\epsilon)^2)$ inside a singular region, $0$ outside.
- **Ignoring joint limits.** A mathematically valid step can drive a joint past its mechanical stop. Clamp every iteration; for redundant arms steer with the null-space term.
- **Local minima and multiplicity.** IK is non-unique (elbow-up vs elbow-down) and the iteration only finds the basin of the initial guess. Seed from the previous solved pose along a trajectory, or try several seeds for a cold start.
- **Linearization overshoot.** The step assumes $J$ is constant over the step. Far from the target the model is wrong and you overshoot (see the $271^\circ$ in the worked example). A line search or step cap on $\|\Delta\theta\|$ tames this.
- **Wrong pseudoinverse branch.** $J^\dagger = J^T(JJ^T)^{-1}$ is only valid for full-row-rank fat $J$. At a singularity $JJ^T$ is itself singular; that is precisely when you must switch to the damped form (which never inverts a singular matrix).

## Practice

1. **Pseudoinverse vs DLS step size.** For the worked 2R robot at $\theta = (0^\circ, 30^\circ)$, compute $\|\Delta\theta\|$ for the plain inverse and for DLS with $\lambda = 0.1, 0.3, 0.5, 1.0$. Plot $\|\Delta\theta\|$ vs $\lambda$ and confirm it decreases monotonically.

2. **Singularity sweep.** Move the 2R target so the arm must pass through full extension ($\theta_2 = 0$). Print $\det(J)$, $\sigma_{\min}(J)$, and $\|\Delta\theta\|$ each iteration for the plain inverse, then for DLS. Show the plain inverse spikes and DLS does not.

3. **Convergence table.** Implement the full 6-DoF body-twist solver (`Vb` from the matrix log) and reproduce a Modern Robotics-style table of $(\theta, \|\omega_b\|, \|v_b\|)$ per iteration on any 6R chain. Check both thresholds $\varepsilon_\omega, \varepsilon_v$ trip together.

4. **Adaptive $\lambda$.** Add the singular-region damping rule from the pitfalls. Compare iteration count to reach $10^{-4}$ against a fixed $\lambda$, on a path that grazes a singularity.

5. **Joint limits via null space.** On a redundant 7-DoF arm, add the null-space term $(I - J^\dagger J)(-k\nabla H)$ with $H$ penalizing limit proximity. Verify the task error stays at machine zero while joints retreat from their stops.

6. **Multiplicity.** From 8 random seeds, solve IK to one fixed reachable pose. Cluster the converged $\theta$ and report how many distinct postures (elbow-up/down, etc.) you recover.

## References

- Kevin M. Lynch and Frank C. Park, *Modern Robotics: Mechanics, Planning, and Control*, Cambridge University Press, 2017. Chapter 6 (Inverse Kinematics), Section 6.2 numerical IK, Newton-Raphson body-twist algorithm and the 2R Example 6.1. Free PDF: https://hades.mech.northwestern.edu/images/7/7f/MR.pdf
- Modern Robotics course resource, 6.2 Numerical Inverse Kinematics: https://modernrobotics.northwestern.edu/nu-gm-book-resource/6-2-numerical-inverse-kinematics-part-1-of-2/
- Samuel R. Buss, "Introduction to Inverse Kinematics with Jacobian Transpose, Pseudoinverse and Damped Least Squares Methods," 2004. The canonical DLS-vs-pseudoinverse-vs-transpose reference with the SVD analysis. https://mathweb.ucsd.edu/~sbuss/ResearchWeb/ikmethods/iksurvey.pdf
- Levenberg-Marquardt algorithm (damped least squares), Wikipedia: https://en.wikipedia.org/wiki/Levenberg%E2%80%93Marquardt_algorithm
- Y. Nakamura and H. Hanafusa, "Inverse Kinematic Solutions with Singularity Robustness for Robot Manipulator Control," ASME J. Dyn. Sys. Meas. Control, 1986. Original singularity-robust DLS for robotics.
- Eigen `JacobiSVD` and `LDLT` documentation: https://eigen.tuxfamily.org/dox/group__SVD__Module.html
- Ghadi Nehme, Eamon Whalen, Faez Ahmed, "CADFit: Precise Mesh-to-CAD Program Generation with Hybrid Optimization," arXiv:2605.01171. CAD-program synthesis target of the stone pipeline. https://arxiv.org/abs/2605.01171 ; code: https://github.com/ghadinehme/CADFit

## Flashcard seeds

- Q: Derive the linearized IK step: starting from forward kinematics $T(\theta)$ and a desired pose, get the equation one Gauss-Newton step must solve for the joint update $\Delta\theta$. <br> <b>Formula:</b> Taylor expand $e(\theta_{i+1}) \approx e(\theta_i) + \frac{\partial e}{\partial\theta}\Delta\theta$ with $\frac{\partial e}{\partial\theta}=J(\theta_i)$, target $e \approx J(\theta_i)\,\Delta\theta$. <br> <i>Hint:</i> set the error after the step to zero, so $0 \approx e(\theta_i) - J\,\Delta\theta$; identify $e \equiv e(\theta_i)$. :: A: Define the task error $e$ between current and goal frame <br> Linearize FK about the current guess: $e(\theta_i+\Delta\theta)\approx e(\theta_i)-J(\theta_i)\Delta\theta$ <br> Demand the linear prediction hit the goal, set it to $0$ <br> $\boxed{e = J(\theta_i)\,\Delta\theta},\quad \Delta\theta=\theta_{i+1}-\theta_i$.
- Q: Derive step 1 of the DLS step: write the regularized objective in $\Delta\theta$ whose minimizer stays bounded near a singularity. <br> <b>Formula:</b> residual $\|J\Delta\theta-e\|^2$ plus a step penalty $\lambda^2\|\Delta\theta\|^2$. <br> <i>Hint:</i> just add the two terms; the penalty is the only new piece beyond plain least squares. :: A: Plain solve $e=J\Delta\theta$ explodes near singularities, so penalize step size <br> Add a Tikhonov term $\lambda^2\|\Delta\theta\|^2$ to the residual <br> $\boxed{\Delta\theta = \arg\min_{\Delta\theta}\ \|J\Delta\theta - e\|^2 + \lambda^2\|\Delta\theta\|^2}$.
- Q: Derive step 2 of the DLS step: from $\min_{\Delta\theta}\|J\Delta\theta-e\|^2+\lambda^2\|\Delta\theta\|^2$ obtain the closed form. <br> <b>Formula:</b> normal equations $(J^TJ+\lambda^2 I)\,\Delta\theta = J^T e$; target $\Delta\theta = J^T(JJ^T + \lambda^2 I)^{-1} e$. <br> <i>Hint:</i> take the gradient in $\Delta\theta$, set it to zero, then apply the push-through identity to move to the cheap $m\times m$ form. :: A: Set the gradient to zero: $2J^T(J\Delta\theta-e)+2\lambda^2\Delta\theta=0$ <br> Collect terms: $(J^TJ+\lambda^2 I)\,\Delta\theta = J^T e$ <br> Invert the SPD matrix: $\Delta\theta=(J^TJ+\lambda^2 I)^{-1}J^T e$ <br> $\boxed{\Delta\theta = J^T(JJ^T + \lambda^2 I)^{-1} e}$.
- Q: Prove the push-through identity that swaps the order of the damped inverse. <br> <b>Formula:</b> target $(J^TJ+\lambda^2 I)^{-1}J^T = J^T(JJ^T+\lambda^2 I)^{-1}$; start from $J^T(JJ^T+\lambda^2 I)=(J^TJ+\lambda^2 I)J^T$. <br> <i>Hint:</i> verify both sides of the start equation expand to $J^TJJ^T+\lambda^2 J^T$, then left- and right-multiply by the two inverses. :: A: Both sides expand to $J^TJJ^T+\lambda^2 J^T$, so $J^T(JJ^T+\lambda^2 I)=(J^TJ+\lambda^2 I)J^T$ <br> Left-multiply by $(J^TJ+\lambda^2 I)^{-1}$, right-multiply by $(JJ^T+\lambda^2 I)^{-1}$ <br> $\boxed{(J^TJ+\lambda^2 I)^{-1}J^T = J^T(JJ^T+\lambda^2 I)^{-1}}$, valid for any $\lambda>0$.
- Q: Derive step 1 of the pseudoinverse step in SVD coordinates from $J=U\Sigma V^T$. <br> <b>Formula:</b> $J^\dagger = V\Sigma^\dagger U^T$ ($\Sigma^\dagger$ sends each $\sigma_i\to 1/\sigma_i$); target $\Delta\theta = \sum_i \frac{1}{\sigma_i}(u_i^T e)\,v_i$. <br> <i>Hint:</i> apply $J^\dagger$ to $e$, then read $U^T e$ as the coefficients $u_i^T e$. :: A: $J^\dagger = V\Sigma^\dagger U^T$, with $\Sigma^\dagger$ inverting each $\sigma_i\to 1/\sigma_i$ <br> Apply to the error: $\Delta\theta = J^\dagger e = V\Sigma^\dagger U^T e$ <br> Expand the sum over columns <br> $\boxed{\Delta\theta = \sum_i \frac{1}{\sigma_i}\,(u_i^T e)\,v_i}$.
- Q: Derive step 2: from $\Delta\theta=\sum_i \frac{1}{\sigma_i}(u_i^Te)v_i$, show why a singularity blows the step up. <br> <b>Formula:</b> the coefficient on $v_i$ is $\frac{1}{\sigma_i}(u_i^Te)$; near a singularity some $\sigma_i\to 0$. <br> <i>Hint:</i> hold $u_i^Te$ fixed and let $\sigma_i\to0$; look at the limit of $1/\sigma_i$. :: A: As the arm nears a singularity, some $\sigma_i\to 0$ <br> The coefficient $1/\sigma_i\to\infty$ even for a tiny $u_i^Te$ <br> $\boxed{\sigma_i\to 0 \Rightarrow \|\Delta\theta\|\to\infty}$: one nearly-lost direction demands near-infinite joint motion.
- Q: Derive step 1 of the SVD form of the DLS step from $\Delta\theta=J^T(JJ^T+\lambda^2 I)^{-1}e$. <br> <b>Formula:</b> $J=U\Sigma V^T$, so $JJ^T+\lambda^2 I=U(\Sigma^2+\lambda^2 I)U^T$; target $\Delta\theta = V\,\Sigma(\Sigma^2+\lambda^2 I)^{-1}U^T e$. <br> <i>Hint:</i> substitute the SVD; the orthogonal $U^TU=I$ cancels in the middle. :: A: Substitute $J=U\Sigma V^T$, so $JJ^T=U\Sigma^2U^T$ and $JJ^T+\lambda^2 I=U(\Sigma^2+\lambda^2 I)U^T$ <br> Then $J^T(JJ^T+\lambda^2 I)^{-1}=V\Sigma U^T\,U(\Sigma^2+\lambda^2 I)^{-1}U^T$ <br> Cancel $U^TU=I$ <br> $\boxed{\Delta\theta = V\,\Sigma(\Sigma^2+\lambda^2 I)^{-1}U^T e}$.
- Q: Derive step 2: from $\Delta\theta=V\,\Sigma(\Sigma^2+\lambda^2 I)^{-1}U^T e$, read off the per-direction filter that replaces $1/\sigma_i$. <br> <b>Formula:</b> the $i$-th diagonal entry is $\frac{\sigma_i}{\sigma_i^2+\lambda^2}$; target $\frac{1}{\sigma_i}\to\frac{\sigma_i}{\sigma_i^2+\lambda^2}$. <br> <i>Hint:</i> compare the entry against $1/\sigma_i$ in the two regimes $\sigma_i\gg\lambda$ and $\sigma_i\to0$. :: A: Each diagonal entry maps $\sigma_i$ through $\sigma_i/(\sigma_i^2+\lambda^2)$ <br> For $\sigma_i\gg\lambda$ this $\approx 1/\sigma_i$; for $\sigma_i\to 0$ this $\approx \sigma_i/\lambda^2\to 0$ <br> $\boxed{\frac{1}{\sigma_i}\ \longrightarrow\ \frac{\sigma_i}{\sigma_i^2+\lambda^2}}$: the reciprocal becomes a bounded filter.
- Q: Derive the maximum gain of the DLS filter and where it occurs. <br> <b>Formula:</b> $g(\sigma)=\frac{\sigma}{\sigma^2+\lambda^2}$, with $g'(\sigma)=\frac{\lambda^2-\sigma^2}{(\sigma^2+\lambda^2)^2}$; target $g_{\max}=\frac{1}{2\lambda}$ at $\sigma=\lambda$. <br> <i>Hint:</i> set the numerator $\lambda^2-\sigma^2$ of $g'$ to zero to find the maximizer, then plug $\sigma=\lambda$ back into $g$. :: A: $g'(\sigma)=\frac{(\sigma^2+\lambda^2)-\sigma(2\sigma)}{(\sigma^2+\lambda^2)^2}=\frac{\lambda^2-\sigma^2}{(\sigma^2+\lambda^2)^2}$ <br> Set $g'=0\Rightarrow \sigma=\lambda$ <br> $g(\lambda)=\frac{\lambda}{2\lambda^2}$ <br> $\boxed{g_{\max}=\frac{1}{2\lambda}\ \text{at}\ \sigma=\lambda}$: no gain ever exceeds $1/(2\lambda)$.
- Q: Derive the large-$\lambda$ limit of the DLS step and name the method it becomes. <br> <b>Formula:</b> $\Delta\theta = J^T(JJ^T+\lambda^2 I)^{-1}e$; target $\Delta\theta \approx \frac{1}{\lambda^2}J^T e$. <br> <i>Hint:</i> when $\lambda^2$ dominates $\|JJ^T\|$, approximate $JJ^T+\lambda^2 I\approx\lambda^2 I$. :: A: $\Delta\theta = J^T(JJ^T+\lambda^2 I)^{-1}e$ <br> For $\lambda^2\gg\|JJ^T\|$, $JJ^T+\lambda^2 I\approx\lambda^2 I$ <br> $\boxed{\Delta\theta \approx \tfrac{1}{\lambda^2}J^T e}$: the Jacobian-transpose / gradient-descent method up to scale.
- Q: Derive the position Jacobian of the planar 2R tip from its forward kinematics. <br> <b>Formula:</b> $x=c_1+c_{12},\ y=s_1+s_{12}$ (with $c_1=\cos\theta_1,\ s_{12}=\sin(\theta_1+\theta_2)$, etc.); target $J=\frac{\partial(x,y)}{\partial(\theta_1,\theta_2)}$. <br> <i>Hint:</i> differentiate each of $x,y$ w.r.t. $\theta_1$ then $\theta_2$; note $\frac{\partial}{\partial\theta_1}c_{12}=-s_{12}$ and $\frac{\partial}{\partial\theta_2}c_{12}=-s_{12}$. :: A: $\partial x/\partial\theta_1=-s_1-s_{12}$, $\partial x/\partial\theta_2=-s_{12}$ <br> $\partial y/\partial\theta_1=c_1+c_{12}$, $\partial y/\partial\theta_2=c_{12}$ <br> $\boxed{J=\begin{bmatrix}-s_1-s_{12} & -s_{12}\\ c_1+c_{12} & c_{12}\end{bmatrix}}$.
- Q: Worked example (2R, well-conditioned, links in m). Compute the DLS step at $\theta=(0^\circ,30^\circ)$ with $J=\begin{bmatrix}-0.5&-0.5\\1.866&0.866\end{bmatrix}$, $e=(-1.5,0.866)$, $\lambda=0.5$. <br> <b>Formula:</b> $\Delta\theta = J^T(JJ^T+\lambda^2 I)^{-1} e$. <br> <i>Hint:</i> first build the $2\times2$ matrix $JJ^T+\lambda^2 I$, solve $(JJ^T+\lambda^2 I)\,y=e$ for $y$, then set $\Delta\theta=J^Ty$. :: A: $JJ^T+\lambda^2 I=\begin{bmatrix}0.75&-1.366\\-1.366&4.484\end{bmatrix}$, $\det=1.497$ <br> $y=(JJ^T+\lambda^2I)^{-1}e=(-3.703,-0.936)$ <br> $\Delta\theta=J^Ty=(0.107,1.041)$ rad <br> $\boxed{\Delta\theta\approx(+6.1^\circ,+59.7^\circ)}$, a gentle safe move.
- Q: Worked example (same 2R, plain inverse for comparison). Compute the undamped step at $\theta=(0^\circ,30^\circ)$ with $e=(-1.5,0.866)$ and $\det J=0.5$. <br> <b>Formula:</b> $\Delta\theta=J^{-1}e$ with $J^{-1}=\frac{1}{\det J}\begin{bmatrix}J_{22}&-J_{12}\\-J_{21}&J_{11}\end{bmatrix}$, $J=\begin{bmatrix}-0.5&-0.5\\1.866&0.866\end{bmatrix}$. <br> <i>Hint:</i> form $J^{-1}$ from the $2\times2$ adjugate over $\det J=0.5$ first, then multiply by $e$. :: A: $J^{-1}=\frac{1}{0.5}\begin{bmatrix}0.866&0.5\\-1.866&-0.5\end{bmatrix}=\begin{bmatrix}1.732&1.0\\-3.732&-1.0\end{bmatrix}$ <br> $\Delta\theta=J^{-1}e=(-1.732,4.732)$ rad <br> $\boxed{\Delta\theta\approx(-99^\circ,+271^\circ)}$: a wild overshoot the damped step tames.
- Q: Worked example (near-singular, arm near full extension). At $\theta=(0^\circ,10^\circ)$ the 2R has $\sigma_{\min}=0.078$; with $e=(0,0.1)$ compare the plain step size to DLS with $\lambda=0.5$. <br> <b>Formula:</b> plain gain $\approx 1/\sigma_{\min}$ so $\|\Delta\theta\|\approx\|e\|/\sigma_{\min}$; DLS gain capped at $\frac{1}{2\lambda}$ so $\|\Delta\theta\|\lesssim\|e\|/(2\lambda)$. <br> <i>Hint:</i> compute the plain magnitude $\|e\|/\sigma_{\min}$ first, then the DLS cap $\|e\|/(2\lambda)$, and compare. :: A: Plain step rides the $1/\sigma_{\min}$ gain ($\sim 1/0.078\approx 13$ on the lost direction): $\|\Delta\theta\|=\|J^{-1}e\|\approx 0.141$ rad <br> DLS gain capped at $1/(2\lambda)=1$, so $\|\Delta\theta\|\approx 0.043$ rad <br> $\boxed{\text{DLS}\approx 0.043\ \text{rad} \ll 0.141\ \text{rad}}$: damping shrinks the near-singular spike $\sim 3\times$.
- Q: Worked example (SVD filter, two regimes). With $\lambda=0.1$, evaluate the reciprocal vs the DLS filter for a small $\sigma=0.05$ and a large $\sigma=2.0$. <br> <b>Formula:</b> reciprocal $1/\sigma$ vs filter $\frac{\sigma}{\sigma^2+\lambda^2}$ (use $\lambda^2=0.01$). <br> <i>Hint:</i> plug each $\sigma$ into both expressions; watch the small $\sigma$ get attenuated while the large $\sigma$ is essentially unchanged. :: A: Small: $1/\sigma=20$ vs $0.05/(0.0025+0.01)=0.05/0.0125=4.0$ (attenuated) <br> Large: $1/\sigma=0.5$ vs $2/(4+0.01)=0.499$ (untouched) <br> $\boxed{\sigma=0.05\!:\,20\to 4;\quad \sigma=2\!:\,0.5\to 0.499}$.
- Q: Test yourself (verify the gain bound). With $\lambda=0.5$, check the DLS filter never exceeds $1/(2\lambda)$. <br> <b>Formula:</b> $g(\sigma)=\frac{\sigma}{\sigma^2+\lambda^2}$, bound $\frac{1}{2\lambda}$, peak at $\sigma=\lambda$. <br> <i>Hint:</i> a pass means $g(\lambda)$ equals the bound exactly and every other $\sigma$ gives a smaller $g$; test $\sigma=0.5,\,2,\,0.05$. :: A: Bound $1/(2\lambda)=1.0$ <br> At $\sigma=\lambda=0.5$: $g=0.5/(0.25+0.25)=1.0$ (hits the bound) <br> At $\sigma=2$: $g=2/4.25=0.47<1$; at $\sigma=0.05$: $g=0.05/0.2525=0.198<1$ <br> $\boxed{\max_\sigma g=1.0=1/(2\lambda)\ \checkmark}$: bound confirmed.
- Q: Test yourself (units check). In the body-twist error $\mathcal{V}_b=(\omega_b,v_b)$, justify why the solver uses two separate stopping thresholds. <br> <b>Formula:</b> stop when $\|\omega_b\|\le\varepsilon_\omega$ (rad) AND $\|v_b\|\le\varepsilon_v$ (m). <br> <i>Hint:</i> a pass means you can name the units of each part and explain why one scalar tolerance cannot compare rad against m. :: A: $\omega_b$ has units rad, $v_b$ has units m, so one scalar tolerance cannot gate both <br> Loop while $\|\omega_b\|>\varepsilon_\omega$ (rad) OR $\|v_b\|>\varepsilon_v$ (m) <br> $\boxed{\varepsilon_\omega,\varepsilon_v\ \text{distinct because rad and m are unlike units}}$.
- Q: Derive the body-frame pose error: explain why poses cannot be subtracted and build the correct error. <br> <b>Formula:</b> relative pose $T_{sb}^{-1}(\theta_i)\,T_{sd}$; target $[\mathcal{V}_b]=\log\!\big(T_{sb}^{-1}(\theta_i)\,T_{sd}\big)$. <br> <i>Hint:</i> first form the relative pose that carries current frame to goal, then take the matrix log to turn it into a body twist. :: A: $SE(3)$ matrices do not subtract to a meaningful pose; rotations live on a manifold <br> Form the relative pose $T_{sb}^{-1}(\theta_i)\,T_{sd}$, then take the matrix log <br> $\boxed{[\mathcal{V}_b]=\log\!\big(T_{sb}^{-1}(\theta_i)\,T_{sd}\big)}$: the constant body twist carrying current frame onto goal in unit time.
- Q: When-to-use: when must you switch from the plain right pseudoinverse to the damped form? <br> <b>Formula:</b> plain $J^\dagger=J^T(JJ^T)^{-1}$ vs damped $J^T(JJ^T+\lambda^2 I)^{-1}$. <br> <i>Hint:</i> ask what happens to $JJ^T$ at a singularity, and which of the two forms still inverts a non-singular matrix there. :: A: $J^\dagger=J^T(JJ^T)^{-1}$ needs full-row-rank fat $J$; at a singularity $JJ^T$ is singular and uninvertible <br> $\boxed{\text{Switch to } J^T(JJ^T+\lambda^2 I)^{-1}\ \text{near/through singularities}}$: $\lambda^2 I$ keeps the matrix invertible.
- Q: When-to-use: how do joint limits enter a redundant-arm ($n>m$) IK step without disturbing the task? <br> <b>Formula:</b> null-space projector $(I-J^\dagger J)$ and secondary motion $z=-k\nabla H(\theta)$; target $\Delta\theta = J^\dagger e + (I-J^\dagger J)\,z$. <br> <i>Hint:</i> add a term living in the null space of $J$, so $J(I-J^\dagger J)=0$ leaves the task error unchanged. :: A: Any $J^\dagger e$ admits a null-space addition that leaves the task error unchanged <br> Pick $z=-k\nabla H(\theta)$ with $H$ penalizing limit proximity <br> $\boxed{\Delta\theta = J^\dagger e + (I-J^\dagger J)\,z}$.
- Q: Pitfall: what goes wrong if you pair a body Jacobian with a spatial-frame error, and what is the rule? <br> <b>Formula:</b> consistent pairs $\mathcal{V}_b$ with $J_b$, $\mathcal{V}_s$ with $J_s$. <br> <i>Hint:</i> the fix is to match the error's frame to the Jacobian's frame; name what a mismatch does to the iteration. :: A: Frame mismatch makes the linear model inconsistent and the iteration silently diverges <br> $\boxed{\text{Match frames: }\mathcal{V}_b\ \text{with}\ J_b,\quad \mathcal{V}_s\ \text{with}\ J_s}$.
- CLOZE: The damped least squares method is the robotics name for {{c1::Levenberg-Marquardt}} regularization applied to the IK residual.
- CLOZE: The DLS filter $\sigma/(\sigma^2+\lambda^2)$ peaks at $\sigma={{c1::\lambda}}$ with value {{c2::1/(2\lambda)}}, which bounds every joint-rate gain. <br> <i>hint:</i> set the derivative's numerator $\lambda^2-\sigma^2$ to zero.
- CLOZE: For the DLS solve to use Cholesky, the matrix $JJ^T + \lambda^2 I$ must be {{c1::symmetric positive definite}}, which holds for any $\lambda > 0$.
- CLOZE: In Eigen, the damped step is computed as `J.transpose() * A.ldlt().solve(e)` where `A` equals {{c1::J * J.transpose() + lambda*lambda*I}}.
- CLOZE: Inverse kinematics is hard because the pose-to-joints map is {{c1::nonlinear and non-unique}}.
