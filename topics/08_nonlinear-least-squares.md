# Nonlinear least squares: Gauss-Newton and Levenberg-Marquardt

Prereqs: [[calculus-jacobians]], [[linear-least-squares]] => this => [[calibration-pnp]], [[inverse-kinematics]], [[cadfit-program-synthesis]]

## Why it matters for the stone scan->CAD->IFC / CADFit workflow

Almost every "make the model match the scan" step in the pipeline is a nonlinear least squares (NLLS) problem. When you register a raw point cloud of an irregular granite block to a reference frame (ICP), you minimize the sum of squared point-to-plane distances over a 6-DoF rigid pose. When you calibrate the camera or projector that captured the scan (PnP / bundle adjustment), you minimize squared reprojection error over intrinsics and poses. When you fit a primitive (plane, cylinder, extruded profile) to a noisy mesh patch, you minimize squared geometric residuals over the primitive parameters. And in a CADFit-style program-synthesis layer, once the discrete operation sequence is chosen, the continuous sketch and extrusion parameters that maximize mesh overlap are tuned by a local NLLS or trust-region refinement on a geometric objective. Gauss-Newton and Levenberg-Marquardt are the two workhorse solvers that make all of these converge from a rough initial guess to a tight fit.

## Intuition

You have a knob vector $x$ and a pile of error numbers $r_i(x)$ (residuals) that you want all near zero. You cannot solve them exactly because there are more errors than knobs and the errors bend nonlinearly. So you do the next best thing: at your current guess, pretend each residual is a straight line in $x$ (its tangent plane), solve the resulting *linear* least squares exactly, step there, and repeat.

That linearize-then-solve loop is Gauss-Newton. It is fast when the tangent-plane picture is accurate, which happens when you are close to the answer and the residuals are small. It misbehaves when the picture is bad: the linear model says "jump far that way" and the jump overshoots into a worse spot.

Levenberg-Marquardt adds a leash. A damping term $\lambda$ shrinks the step and bends its direction toward plain downhill (gradient descent). Trust the linear model (small $\lambda$, Gauss-Newton-like) when it keeps predicting real improvement. Distrust it (large $\lambda$, tiny cautious steps) when it lies. The analogy: Gauss-Newton is a confident driver taking the racing line; LM is the same driver who eases off the gas exactly when the road gets foggy.

## The math (first principles)

### Problem statement

We minimize a sum of squared residuals. Let $x \in \mathbb{R}^n$ be the parameters and $r(x) = (r_1(x), \dots, r_m(x))^\top \in \mathbb{R}^m$ the residual vector, with $m \ge n$. Define the cost

$$F(x) = \tfrac{1}{2} \sum_{i=1}^{m} r_i(x)^2 = \tfrac{1}{2}\, r(x)^\top r(x) = \tfrac{1}{2}\,\|r(x)\|^2.$$

The factor $\tfrac{1}{2}$ is a convention that cancels the 2 from differentiation. It does not change the minimizer. Some texts drop it; track which one your code uses so gradients match.

### Jacobian, gradient, Hessian

The Jacobian collects all first partials, $J(x) \in \mathbb{R}^{m \times n}$:

$$J_{ij}(x) = \frac{\partial r_i}{\partial x_j}(x).$$

The gradient of $F$ is exact:

$$\nabla F(x) = J(x)^\top r(x) \equiv g.$$

The exact Hessian is

$$\nabla^2 F(x) = J^\top J + \sum_{i=1}^{m} r_i\, \nabla^2 r_i.$$

The second term needs per-residual second derivatives and is weighted by the residuals $r_i$ themselves. The Gauss-Newton idea is to drop it:

$$\nabla^2 F(x) \approx J^\top J.$$

This is exact when the residuals are affine in $x$, and a good approximation when residuals are small or nearly linear near the solution. It is cheap (no second derivatives) and $J^\top J$ is symmetric positive semidefinite by construction.

### Gauss-Newton update

Linearize the residual at $x_k$: $r(x_k + \Delta x) \approx r(x_k) + J\,\Delta x$. Substituting into $F$ and minimizing the quadratic over $\Delta x$ gives the **normal equations**

$$(J^\top J)\,\Delta x = -\,J^\top r,$$

so the update is

$$\boxed{\,x_{k+1} = x_k - (J^\top J)^{-1} J^\top r\,}$$

with $J$, $r$ evaluated at $x_k$. In practice you never invert: you solve the linear system, ideally via the QR factorization of $J$ (more numerically stable than forming $J^\top J$, which squares the condition number).

Convergence is locally superlinear, and quadratic in the zero-residual case. It can diverge when started far from the solution or when $J$ is rank-deficient (then $J^\top J$ is singular).

### Levenberg-Marquardt damping

LM replaces the possibly-singular $J^\top J$ with a damped, always-invertible matrix. Two damping forms are standard:

$$(J^\top J + \lambda I)\,\Delta x = -\,J^\top r \qquad\text{(Levenberg)}$$

$$(J^\top J + \lambda\,\operatorname{diag}(J^\top J))\,\Delta x = -\,J^\top r \qquad\text{(Marquardt)}$$

The diagonal-scaled form makes the method invariant to the units/scale of each parameter, so it is the usual default (Marquardt's contribution). Interpolation behavior:

- $\lambda \to 0$: recover the Gauss-Newton step (fast, trusts the linear model).
- $\lambda \to \infty$: $\Delta x \to -\tfrac{1}{\lambda} J^\top r$, a short step along the negative gradient (steepest descent).

So $\lambda$ slides continuously between the two. Large $\lambda$ also shrinks $\|\Delta x\|$, which is the trust-region connection: limiting step length to a ball of radius $\Delta$ is equivalent, by Lagrange duality, to adding some $\lambda I$. LM is "implicit trust region" because it controls $\lambda$ directly instead of the radius.

### The gain ratio (step accept/reject and $\lambda$ control)

After a trial step, compare the actual cost drop to the drop the linear model predicted:

$$\rho = \frac{F(x_k) - F(x_k + \Delta x)}{\tfrac{1}{2}\,\Delta x^\top(\lambda\,\Delta x - g)}$$

where $g = J^\top r$ and the denominator is the predicted reduction (always positive). Decision rule (Madsen-Nielsen-Tingleff):

- If $\rho > 0$: accept the step. Decrease damping, e.g. $\lambda \gets \lambda \cdot \max\!\big(\tfrac{1}{3},\, 1 - (2\rho - 1)^3\big)$, and reset a multiplier $\nu \gets 2$.
- If $\rho \le 0$: reject the step (keep $x_k$). Increase damping, $\lambda \gets \lambda \cdot \nu$, then $\nu \gets 2\nu$.

A common initialization is $\lambda_0 = \tau \cdot \max_i (J^\top J)_{ii}$ with $\tau \in [10^{-6}, 10^{-3}]$ ($10^{-3}$ for a poor initial guess).

### Stopping criteria

Stop when the problem is solved to tolerance, using *relative* thresholds so the test is scale-aware:

- gradient small: $\|J^\top r\|_\infty \le \varepsilon_1$,
- step small: $\|\Delta x\| \le \varepsilon_2 (\|x\| + \varepsilon_2)$,
- iteration cap reached.

### Conventions to fix up front

- **Residual sign**: define $r = \text{model} - \text{measurement}$ consistently. Sign flips do not change $F$ but do flip $g$ and break any hand-checked update.
- **Weighting**: for heteroscedastic data use $r_i \gets r_i / \sigma_i$ (or a full $W = \Sigma^{-1/2}$). This is just a rescaling of residuals; the same solver applies.
- **Frame and parameterization** (geometry): for rigid pose, do not optimize a 9-entry rotation matrix directly. Use a minimal 3-parameter local update (an $\mathfrak{so}(3)$ tangent vector) composed onto the current rotation, so $J$ is $m \times 6$, not $m \times 12$. Optimizing over-parameterized rotations makes $J^\top J$ singular.

## Worked example

Fit the one-parameter model $y = e^{a t}$ to two data points, by hand. Data: $(t_1,y_1)=(0,\,2)$, $(t_2,y_2)=(1,\,3)$. There is no exact fit (at $t=0$ the model is forced to $1 \ne 2$), so this is a genuine least squares problem in the single parameter $a$.

Residuals (model minus data):

$$r_1(a) = e^{a\cdot 0} - 2 = 1 - 2 = -1 \quad(\text{constant}),\qquad r_2(a) = e^{a} - 3.$$

Jacobian (a $2\times 1$ column), $\partial r_i / \partial a$:

$$J = \begin{bmatrix} 0 \\ e^{a} \end{bmatrix}.$$

Start at $a_0 = 1$. Then $e^{1} = 2.7183$.

- $r = (-1,\; 2.7183 - 3) = (-1,\; -0.2817)$.
- $J = (0,\; 2.7183)^\top$.
- $J^\top J = 0^2 + 2.7183^2 = 7.3891$.
- $J^\top r = 0\cdot(-1) + 2.7183\cdot(-0.2817) = -0.7658$.

Gauss-Newton step:

$$\Delta a = -(J^\top J)^{-1} J^\top r = -\frac{-0.7658}{7.3891} = +0.1036.$$

So $a_1 = 1 + 0.1036 = 1.1036$.

Check the cost dropped. $F(a_0) = \tfrac12(1 + 0.2817^2) = \tfrac12(1.0794) = 0.5397$. At $a_1$: $e^{1.1036} = 3.0150$, so $r = (-1,\,0.0150)$, $F(a_1) = \tfrac12(1 + 0.000225) = 0.5001$. Cost fell from $0.5397$ to $0.5001$. The second residual went from $-0.282$ to $+0.015$, nearly zeroed in one step because near the solution the model is almost linear in $a$.

Now the LM version of the same first step with $\lambda = 1$ (Levenberg form): solve $(7.3891 + 1)\,\Delta a = 0.7658$, giving $\Delta a = 0.7658 / 8.3891 = 0.0913$, a slightly shorter, more cautious step in the same direction. With $\lambda = 100$ it becomes $0.7658/107.39 = 0.00713$, almost pure gradient descent. That is the damping knob in action.

## Code

Idiomatic, self-contained C++ with Eigen. This hand-rolled Levenberg-Marquardt (Madsen-Nielsen-Tingleff variant) is short, dependency-light, and shows every symbol. It fits an exponential decay $y = A\,e^{-b t}$ to data, a stand-in for any geometric residual you would drop in later.

```cpp
// g++ -O2 -I /path/to/eigen lm.cpp -o lm
#include <Eigen/Dense>
#include <cmath>
#include <iostream>

using Eigen::VectorXd;
using Eigen::MatrixXd;

// Model y = x(0) * exp(-x(1) * t). Replace residual()/jacobian() with your geometry.
struct ExpModel {
    VectorXd t, y;                       // measured data

    // r_i = model(t_i) - y_i   (sign convention: model minus measurement)
    VectorXd residual(const VectorXd& x) const {
        return (x(0) * (-x(1) * t.array()).exp()).matrix() - y;
    }
    // J_ij = d r_i / d x_j
    MatrixXd jacobian(const VectorXd& x) const {
        const auto e = (-x(1) * t.array()).exp();           // exp(-b t_i)
        MatrixXd J(t.size(), 2);
        J.col(0) = e.matrix();                              // d r / d A
        J.col(1) = (-x(0) * t.array() * e).matrix();        // d r / d b
        return J;
    }
};

VectorXd levenberg_marquardt(const ExpModel& m, VectorXd x, int max_iter = 100) {
    const double eps_grad = 1e-12, eps_step = 1e-12, tau = 1e-3;

    VectorXd r = m.residual(x);
    MatrixXd J = m.jacobian(x);
    MatrixXd A = J.transpose() * J;                         // approx Hessian J^T J
    VectorXd g = J.transpose() * r;                         // gradient J^T r
    double lambda = tau * A.diagonal().maxCoeff();          // initial damping
    double nu = 2.0;

    for (int k = 0; k < max_iter; ++k) {
        if (g.lpNorm<Eigen::Infinity>() <= eps_grad) break; // gradient ~ 0

        // (J^T J + lambda I) dx = -J^T r   (Levenberg form; swap I for A.diagonal() for Marquardt)
        MatrixXd H = A;
        H.diagonal().array() += lambda;
        VectorXd dx = H.ldlt().solve(-g);                   // SPD solve, no explicit inverse

        if (dx.norm() <= eps_step * (x.norm() + eps_step)) break;

        VectorXd x_new = x + dx;
        VectorXd r_new = m.residual(x_new);
        double F_old = 0.5 * r.squaredNorm();
        double F_new = 0.5 * r_new.squaredNorm();
        // predicted reduction = 0.5 * dx^T (lambda*dx - g)
        double pred = 0.5 * dx.dot(lambda * dx - g);
        double rho = (F_old - F_new) / (pred > 0 ? pred : 1e-30);

        if (rho > 0) {                                      // accept step
            x = x_new;
            r = r_new;
            J = m.jacobian(x);
            A = J.transpose() * J;
            g = J.transpose() * r;
            double f = 1.0 - std::pow(2.0 * rho - 1.0, 3);
            lambda *= std::max(1.0 / 3.0, f);               // shrink damping
            nu = 2.0;
        } else {                                            // reject: stay put, damp harder
            lambda *= nu;
            nu *= 2.0;
        }
    }
    return x;
}

int main() {
    // Synthetic data from A=5, b=0.5 plus the implied fit.
    VectorXd t(7); t << 0, 1, 2, 3, 4, 5, 6;
    VectorXd y(7);
    for (int i = 0; i < t.size(); ++i) y(i) = 5.0 * std::exp(-0.5 * t(i));

    ExpModel m{t, y};
    VectorXd x0(2); x0 << 1.0, 1.0;                          // deliberately rough guess
    VectorXd x = levenberg_marquardt(m, x0);
    std::cout << "A = " << x(0) << ", b = " << x(1) << "\n"; // -> ~5, ~0.5
    return 0;
}
```

To switch to Gauss-Newton, set `lambda = 0` and always accept; to switch to Marquardt scaling, replace `H.diagonal().array() += lambda;` with `H.diagonal() += lambda * A.diagonal();`.

For production, prefer a battle-tested library. [Ceres Solver](http://ceres-solver.org/) is the standard for vision/geometry: you write a residual functor, Ceres builds the Jacobian by automatic differentiation and runs LM with a proper trust region. Eigen also ships `unsupported/Eigen/LevenbergMarquardt`.

A tiny Python check with SciPy, useful to validate your C++ output:

```python
import numpy as np
from scipy.optimize import least_squares
t = np.arange(7.0)
y = 5.0 * np.exp(-0.5 * t)
resid = lambda x: x[0] * np.exp(-x[1] * t) - y     # model minus data
sol = least_squares(resid, x0=[1.0, 1.0], method="lm")  # Levenberg-Marquardt
print(sol.x)   # -> [5. 0.5]
```

## Connections

- [[calculus-jacobians]]: the Jacobian $J = \partial r / \partial x$ is the one object every step needs; getting its analytic form (or trusting autodiff) is prerequisite.
- [[linear-least-squares]]: each Gauss-Newton iteration *is* a linear least squares solve of $\min_{\Delta x}\|J\,\Delta x + r\|^2$; NLLS is LLS in a loop.
- [[calibration-pnp]]: PnP and camera calibration minimize squared reprojection error; LM is the solver behind OpenCV `solvePnP` refinement and `calibrateCamera`.
- [[inverse-kinematics]]: damped least squares IK, $\Delta\theta = J^\top(JJ^\top + \lambda^2 I)^{-1} e$, is the same damping idea applied to the kinematic Jacobian to stay bounded near singularities.
- [[cadfit-program-synthesis]]: once the discrete CAD operation graph is fixed, the continuous sketch/extrude parameters are refined against a geometric (IoU or chamfer) objective; LM-style local refinement tightens the fit even though the outer search is discrete.
- [[icp-registration]]: point-to-plane ICP is NLLS over a 6-DoF pose; LM gives the per-iteration pose increment.

## Pitfalls and failure modes

- **Rank-deficient or ill-conditioned $J^\top J$**: if parameters are redundant or data is degenerate (e.g. fitting a cylinder axis to a near-flat patch), $J^\top J$ is singular and Gauss-Newton blows up. LM damping hides this but the result is unconstrained along the null direction. Inspect $\operatorname{cond}(J)$ or the singular values; add priors or remove parameters.
- **Forming $J^\top J$ explicitly squares the condition number**. For poorly conditioned problems solve the linear least squares via QR or SVD of $J$ instead of factoring $J^\top J$.
- **Big residuals break the Gauss-Newton Hessian approximation**. The dropped term $\sum_i r_i \nabla^2 r_i$ matters when $r_i$ are large (outliers, bad model). Convergence slows or stalls; use LM, robust losses (Huber, Cauchy), or true Newton.
- **Over-parameterized rotations/poses**. Optimizing a full rotation matrix or unnormalized quaternion creates exact null directions in $J$. Parameterize on the tangent space ($\mathfrak{so}(3)$ / $\mathfrak{se}(3)$) and compose increments. Rotation increments do not commute, so apply the local update on the correct side.
- **Sign and frame inconsistency**. A flipped residual sign or mixing world/camera frames produces a gradient that points the wrong way; the solver then climbs the cost. Unit-test the gradient against a finite-difference of $F$.
- **Scale-dependent damping and tolerances**. A fixed $\lambda$ or absolute stopping tolerance fails when parameters span very different magnitudes (millimeters of position vs radians of angle). Use Marquardt's diagonal scaling and relative tolerances.
- **Local minima and initialization**. NLLS is local. ICP, PnP, and CAD fitting all need a decent starting pose/guess; LM converges to the nearest basin, not the global optimum. Do a coarse global step (RANSAC, feature match, parameter sweep) first.
- **Finite-difference Jacobians**. Numerical $J$ with a bad step size injects noise; too large a step biases, too small amplifies roundoff. Prefer analytic or automatic differentiation.

## Practice

1. **By hand, one LM step.** Reuse the worked example ($y = e^{a t}$, data $(0,2),(1,3)$, start $a_0 = 1$). Compute the LM step with $\lambda = 10$ in the Levenberg form and confirm it is shorter than and same-sign as the Gauss-Newton step.
2. **Gauss-Newton vs LM divergence.** Modify the C++ to a pure Gauss-Newton loop (no damping, always accept) and run it on the exponential fit from a very bad guess like `x0 = (10, 10)`. Observe divergence or oscillation, then re-enable damping and watch it converge. Print cost per iteration.
3. **Gradient check.** Add a finite-difference gradient $\tilde g_j = (F(x + h e_j) - F(x - h e_j)) / (2h)$ with $h = 10^{-6}$ and assert it matches `J.transpose() * r` to within $10^{-4}$. This is the single most useful debugging tool for new residuals.
4. **Robust loss.** Replace the squared cost with a Huber loss (linear past a threshold) by re-weighting residuals each iteration (iteratively reweighted least squares) and verify it ignores an injected outlier in the data.
5. **Plane fit to a scan patch (geometry).** Write residuals $r_i = n^\top p_i + d$ for points $p_i$, with $n$ a unit normal (parameterize by two angles to avoid the redundant scale) and offset $d$. Fit a plane to a noisy point set and report the RMS residual. This is the kernel inside point-to-plane ICP.
6. **Marquardt vs Levenberg scaling.** Run the plane fit with parameters in mixed units (meters and radians). Compare convergence of the $\lambda I$ form vs the $\lambda\operatorname{diag}(J^\top J)$ form and explain the difference.

## References

- K. Madsen, H. B. Nielsen, O. Tingleff, *Methods for Non-Linear Least Squares Problems*, 2nd ed., Technical University of Denmark, 2004. https://www2.imm.dtu.dk/pubdb/edoc/imm3215.pdf (gain ratio, damping update rule, stopping criteria used above).
- J. Nocedal, S. Wright, *Numerical Optimization*, 2nd ed., Springer, 2006, Ch. 10 (Least-Squares Problems: Gauss-Newton, Levenberg-Marquardt, trust region equivalence).
- D. W. Marquardt, "An Algorithm for Least-Squares Estimation of Nonlinear Parameters," *J. SIAM*, 11(2):431-441, 1963. https://doi.org/10.1137/0111030 (diagonal-scaled damping).
- Ethan Eade, "Gauss-Newton / Levenberg-Marquardt Optimization." https://www.ethaneade.com/optimization.pdf (concise derivation, manifold parameterization).
- Eigen unsupported, Non-linear optimization module (LevenbergMarquardt, NumericalDiff). https://eigen.tuxfamily.org/dox/unsupported/group__NonLinearOptimization__Module.html
- Ceres Solver documentation (trust-region NLLS, autodiff). http://ceres-solver.org/nnls_tutorial.html
- SciPy `scipy.optimize.least_squares` (method="lm" / "trf"). https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.least_squares.html
- OpenCV calib3d, `solvePnP` / `calibrateCamera` (LM-based pose and intrinsics refinement). https://docs.opencv.org/4.x/d9/d0c/group__calib3d.html
- K. Lynch, F. Park, *Modern Robotics*, Cambridge, 2017, Ch. 6 (numerical IK, damped least squares). https://hades.mech.northwestern.edu/index.php/Modern_Robotics
- G. Nehme, E. Whalen, F. Ahmed, "CADFit: Precise Mesh-to-CAD Program Generation with Hybrid Optimization," arXiv:2605.01171, 2026. https://arxiv.org/abs/2605.01171 and https://github.com/ghadinehme/CADFit (geometric IoU objective for CAD parameter fitting).

## Flashcard seeds

- Q: Derive the gradient $\nabla F$ of the cost (derive, step 1 of 2). <br> <b>Formula:</b> start from $F(x)=\tfrac12 r(x)^\top r(x)$, target $\nabla F = J^\top r$ with $J=\partial r/\partial x$. <br> <i>Hint:</i> use the chain rule on the quadratic in $r$; the $\tfrac12$ cancels the factor 2. :: A: $\nabla F = \big(\tfrac{\partial r}{\partial x}\big)^\top r$ <br> name the Jacobian $J=\partial r/\partial x$ <br> $\boxed{\nabla F = J^\top r \equiv g}$
- Q: From the gradient, derive the exact Hessian $\nabla^2 F$ (derive, step 2 of 2). <br> <b>Formula:</b> start from $g_j=\sum_i J_{ij} r_i$, target $\nabla^2 F = J^\top J + \sum_i r_i\,\nabla^2 r_i$. <br> <i>Hint:</i> differentiate $g_j$ w.r.t. $x_k$ with the product rule; one factor gives $J_{ik}$, the other gives $\partial^2 r_i$. :: A: $\partial g_j/\partial x_k=\sum_i J_{ij}J_{ik} + \sum_i r_i\,\partial^2 r_i/\partial x_j\partial x_k$ <br> the first sum is $(J^\top J)_{jk}$ <br> $\boxed{\nabla^2 F = J^\top J + \sum_i r_i\,\nabla^2 r_i}$
- Q: From the exact Hessian, derive the Gauss-Newton approximation and say why it is cheap. <br> <b>Formula:</b> exact $\nabla^2 F = J^\top J + \sum_i r_i\,\nabla^2 r_i$, target $\nabla^2 F \approx J^\top J$. <br> <i>Hint:</i> decide which term needs second derivatives and drop it; note when it is small. :: A: Drop the second-derivative term $\sum_i r_i\,\nabla^2 r_i$ (small when residuals are small or the model is nearly linear) <br> no per-residual second derivatives are needed, and what remains is symmetric PSD <br> $\boxed{\nabla^2 F \approx J^\top J}$
- Q: Derive the Gauss-Newton normal equations from the linearized residual (derive, step 1 of 2). <br> <b>Formula:</b> tangent model $r(x_k+\Delta x)\approx r+J\,\Delta x$; minimize $F=\tfrac12\|r+J\Delta x\|^2$; target $(J^\top J)\Delta x=-J^\top r$. <br> <i>Hint:</i> expand the squared norm in $\Delta x$, then set $\nabla_{\Delta x}F=0$. :: A: $F=\tfrac12(r^\top r + 2\,\Delta x^\top J^\top r + \Delta x^\top J^\top J\,\Delta x)$ <br> $\nabla_{\Delta x}F = J^\top r + J^\top J\,\Delta x = 0$ <br> $\boxed{(J^\top J)\,\Delta x = -J^\top r}$
- Q: From the normal equations, derive the Gauss-Newton update $x_{k+1}$ (derive, step 2 of 2). <br> <b>Formula:</b> $(J^\top J)\Delta x=-J^\top r$, target $x_{k+1}=x_k-(J^\top J)^{-1}J^\top r$. <br> <i>Hint:</i> solve for $\Delta x$, then add it to $x_k$. :: A: $\Delta x = -(J^\top J)^{-1}J^\top r$ <br> $x_{k+1}=x_k+\Delta x$ <br> $\boxed{x_{k+1}=x_k-(J^\top J)^{-1}J^\top r}$
- Q: Construct the Levenberg-damped form that is always invertible. <br> <b>Formula:</b> from $(J^\top J)\Delta x=-J^\top r$, target $(J^\top J + \lambda I)\,\Delta x = -J^\top r$. <br> <i>Hint:</i> $J^\top J$ can be singular; add a positive multiple of $I$ to make it positive definite. :: A: $J^\top J$ may be singular/rank-deficient <br> add $\lambda I$ ($\lambda>0$) so the matrix is symmetric positive definite <br> $\boxed{(J^\top J + \lambda I)\,\Delta x = -J^\top r}$
- Q: From the Levenberg form, write the Marquardt scale-invariant variant and say what changes. <br> <b>Formula:</b> from $(J^\top J+\lambda I)\Delta x=-J^\top r$, target $(J^\top J + \lambda\,\operatorname{diag}(J^\top J))\Delta x = -J^\top r$. <br> <i>Hint:</i> replace the unit $I$ by something carrying each parameter's own curvature scale. :: A: Replace $I$ by $\operatorname{diag}(J^\top J)$ so each parameter is damped in its own units (unit-invariant) <br> $\boxed{(J^\top J + \lambda\,\operatorname{diag}(J^\top J))\,\Delta x = -J^\top r}$
- Q: Take the limit $\lambda\to\infty$ of the Levenberg system and identify the method. <br> <b>Formula:</b> $(J^\top J+\lambda I)\Delta x=-g$ with $g=J^\top r$; target $\Delta x\to-\tfrac1\lambda J^\top r$. <br> <i>Hint:</i> for huge $\lambda$, drop $J^\top J$ next to $\lambda I$, then solve. :: A: $\lambda I\,\Delta x \approx -g$ <br> $\Delta x \approx -\tfrac1\lambda g = -\tfrac1\lambda J^\top r$ <br> $\boxed{\Delta x \to -\tfrac1\lambda\,J^\top r}$, a short steepest-descent step
- CLOZE: As $\lambda \to 0$ Levenberg-Marquardt reduces to {{c1::Gauss-Newton}}, and as $\lambda \to \infty$ it reduces to {{c1::gradient descent}}.
- Q: Construct the LM gain ratio $\rho$ used to accept/reject a trial step. <br> <b>Formula:</b> $\rho=\dfrac{F(x_k)-F(x_k+\Delta x)}{\tfrac12\Delta x^\top(\lambda\Delta x - g)}$, with $g=J^\top r$. <br> <i>Hint:</i> numerator = actual cost drop; denominator = the drop the linear model predicted (always positive). :: A: actual drop $=F(x_k)-F(x_k+\Delta x)$ <br> predicted drop $=\tfrac12\Delta x^\top(\lambda\Delta x - g)>0$ <br> $\boxed{\rho=\dfrac{F(x_k)-F(x_k+\Delta x)}{\tfrac12\Delta x^\top(\lambda\Delta x - g)}}$
- Q: Construct the damped-least-squares IK step (the LM idea applied to inverse kinematics). <br> <b>Formula:</b> kinematic error $e$, Jacobian $J$; target $\Delta\theta = J^\top(JJ^\top + \lambda^2 I)^{-1}e$. <br> <i>Hint:</i> start from $\Delta\theta=J^\dagger e$, then damp the Gram matrix $JJ^\top$ by $\lambda^2 I$ before inverting. :: A: bare $\Delta\theta=J^\dagger e$ blows up near singularities <br> damp the right-inverse Gram matrix by $\lambda^2 I$ <br> $\boxed{\Delta\theta = J^\top(JJ^\top + \lambda^2 I)^{-1}e}$
- Q: Worked example, GN step (small/near regime). Fit $y=e^{at}$ to $(0,2),(1,3)$ at $a_0=1$. <br> <b>Formula:</b> $\Delta a = -(J^\top J)^{-1}J^\top r$, given $J^\top J=7.3891$, $J^\top r=-0.7658$; then $a_1=a_0+\Delta a$. <br> <i>Hint:</i> compute $\Delta a$ first (watch the double negative), then add to $a_0=1$. :: A: $\Delta a = -(-0.7658)/7.3891$ <br> $\Delta a = +0.1036$ <br> $\boxed{a_1 = 1+0.1036 = 1.1036}$
- Q: Worked example, LM step (mild damping). Same $e^{at}$ setup. <br> <b>Formula:</b> $(J^\top J+\lambda)\Delta a = -J^\top r$ with $J^\top J=7.3891$, $-J^\top r=0.7658$, $\lambda=1$. <br> <i>Hint:</i> add $\lambda$ to the $7.3891$ in the denominator, then divide $0.7658$ by it. :: A: $(7.3891+1)\Delta a = 0.7658$ <br> $\Delta a = 0.7658/8.3891$ <br> $\boxed{\Delta a = 0.0913}$, shorter than the GN $0.1036$, same sign
- Q: Worked example, LM step (heavy damping / steepest-descent regime). Same $e^{at}$ setup. <br> <b>Formula:</b> $(J^\top J+\lambda)\Delta a = -J^\top r$ with $J^\top J=7.3891$, $-J^\top r=0.7658$, $\lambda=100$. <br> <i>Hint:</i> $\lambda$ now dwarfs $7.3891$, so the denominator is $\approx 107.39$; name the regime after dividing. :: A: $(7.3891+100)\Delta a = 0.7658$ <br> $\Delta a = 0.7658/107.39$ <br> $\boxed{\Delta a = 0.00713}$, almost pure gradient descent
- Q: Worked example, large-residual 2-point linear fit (scalar parameter). Residuals $r=(2,-4)$, constant Jacobian $J=(1,1)^\top$. <br> <b>Formula:</b> $\Delta x = -(J^\top J)^{-1}J^\top r$. <br> <i>Hint:</i> compute the scalars $J^\top J$ and $J^\top r$ first by summing the two components, then divide. :: A: $J^\top J = 1^2+1^2 = 2$, $J^\top r = 2+(-4)=-2$ <br> $\Delta x = -(-2)/2$ <br> $\boxed{\Delta x = +1}$ (exact in one step: residuals are affine)
- Q: Test yourself (close the loop). For the $e^{at}$ example, $a_1=1.1036$ gives $e^{1.1036}=3.0150$. <br> <b>Formula:</b> $F(a)=\tfrac12\sum_i r_i^2$ with $r=(e^{a\cdot0}-2,\,e^{a}-3)$; compare to $F(a_0)=0.5397$. <br> <i>Hint:</i> a pass is $F(a_1)<F(a_0)$; build the new residuals first, then square-and-halve. :: A: new residuals $r=(1-2,\,3.0150-3)=(-1,\,0.0150)$ <br> $F(a_1)=\tfrac12(1^2+0.0150^2)=\tfrac12(1.000225)=0.5001$ <br> $\boxed{0.5001 < 0.5397}$: cost fell, step accepted ($\rho>0$)
- Q: Test yourself (verify the minimizer). Confirm the linear-fit GN step $\Delta x=+1$ from $r=(2,-4)$, $J=(1,1)^\top$ actually minimizes $F$. <br> <b>Formula:</b> updated residual $r_{new}=r+J\Delta x$; minimum iff gradient $J^\top r_{new}=0$. <br> <i>Hint:</i> a pass is a zero gradient; form $r_{new}$ first, then dot with $J$. :: A: $r_{new}=(2+1,\,-4+1)=(3,-3)$ <br> $J^\top r_{new}=3+(-3)=0$ <br> $\boxed{\nabla F=0\Rightarrow \text{exact minimizer in one step}}$
- Q: Test yourself (verify damping shrinks the step). For $-J^\top r=0.7658$, $J^\top J=7.3891$, compare the $\lambda=1$ and $\lambda=100$ LM steps. <br> <b>Formula:</b> $|\Delta a|=0.7658/(7.3891+\lambda)$ for each $\lambda$. <br> <i>Hint:</i> a pass is ratio $<1$; compute both magnitudes, then divide the larger-$\lambda$ one by the smaller-$\lambda$ one. :: A: $|\Delta a|_{\lambda=1}=0.7658/8.3891=0.0913$; $|\Delta a|_{\lambda=100}=0.7658/107.39=0.00713$ <br> ratio $0.00713/0.0913\approx 0.078<1$ <br> $\boxed{\text{larger }\lambda\Rightarrow\text{shorter step, same direction}}$
- Q: When is the Gauss-Newton Hessian approximation $\nabla^2 F\approx J^\top J$ exact? <br> <b>Formula:</b> exact Hessian $\nabla^2 F = J^\top J + \sum_i r_i\,\nabla^2 r_i$. <br> <i>Hint:</i> ask when the dropped sum is exactly zero. :: A: When the residuals are affine (linear) in the parameters, so every $\nabla^2 r_i=0$ and the dropped term vanishes <br> $\boxed{r \text{ affine} \Rightarrow \nabla^2 F = J^\top J}$
- Q: Definition: what does nonlinear least squares minimize, and what is a residual? <br> <b>Formula:</b> $F(x)=\tfrac12\sum_i r_i(x)^2$ with $r_i=\text{model}-\text{measurement}$. <br> <i>Hint:</i> name what the cost sums and what a single residual measures. :: A: It minimizes $F(x)=\tfrac12\sum_i r_i(x)^2$ <br> each residual $r_i=\text{model}-\text{measurement}$ is an error number you want near zero <br> $\boxed{\min_x \tfrac12\|r(x)\|^2}$
- Q: When to use LM over plain Gauss-Newton? <br> <i>Hint:</i> think about the three cases where the undamped step $-(J^\top J)^{-1}J^\top r$ overshoots or fails. :: A: When started far from the solution, with large residuals/outliers, or a rank-deficient $J^\top J$ <br> there the undamped GN step overshoots or diverges, and $\lambda$ leashes it <br> $\boxed{\text{far start / big residuals / singular } J^\top J}$
- Q: Pitfall: why solve $\min_{\Delta x}\|J\Delta x+r\|$ via QR/SVD of $J$ rather than forming $J^\top J$? <br> <b>Formula:</b> $\operatorname{cond}(J^\top J)=\operatorname{cond}(J)^2$. <br> <i>Hint:</i> compare conditioning of working with $J$ directly vs its product with itself. :: A: Forming $J^\top J$ squares the condition number <br> on ill-conditioned problems that doubles the lost digits, so factor $J$ directly via QR/SVD <br> $\boxed{\operatorname{cond}(J^\top J)=\operatorname{cond}(J)^2}$
- Q: Pitfall: how must a rigid rotation be parameterized inside NLLS, and why? <br> <i>Hint:</i> count the free parameters of a rotation vs entries of a $3\times3$ matrix; the gap creates null directions in $J$. :: A: Use a 3-parameter $\mathfrak{so}(3)$ tangent increment composed onto the current rotation <br> a full 9-entry matrix or unnormalized quaternion is over-parameterized and makes $J^\top J$ singular <br> $\boxed{\text{optimize on the }\mathfrak{so}(3)\text{ tangent, } J \text{ is } m\times6}$
- CLOZE: The LM gain ratio is {{c1::actual cost reduction divided by predicted reduction}}; a value above zero means {{c2::accept the step and decrease $\lambda$}}.
- CLOZE: When a trial step is rejected ($\rho\le 0$) LM keeps $x_k$ and {{c1::increases $\lambda$}}, e.g. $\lambda\gets\lambda\nu$ with $\nu$ doubling. <br> <i>hint:</i> a bad step means distrust the linear model, so damp harder.
- CLOZE: Large LM damping both shrinks the step length and bends its direction toward {{c1::the negative gradient}}, the trust-region connection.
- CLOZE: The single most useful debugging check for a new residual is comparing `J.transpose()*r` against a {{c1::finite-difference gradient}} of the cost. <br> <i>hint:</i> $\tilde g_j=(F(x+he_j)-F(x-he_j))/(2h)$ with $h\approx10^{-6}$.
- Q: In Eigen, how do you solve the SPD system $(J^\top J+\lambda I)\Delta x=-g$ without inverting? <br> <b>Formula:</b> the system is symmetric positive definite, so use an $LDL^\top$ (Cholesky-type) solve. <br> <i>Hint:</i> form $H=J^\top J+\lambda I$, then call the factorization that takes $-g$ as the right-hand side. :: A: build $H = J^\top J + \lambda I$ <br> solve via the SPD factorization, no explicit inverse <br> $\boxed{\texttt{(JtJ + lambda*I).ldlt().solve(-g)}}$
