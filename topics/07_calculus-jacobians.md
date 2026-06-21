# Calculus for optimization: gradients, Jacobians, linearization

Prereqs: [[05_vectors]], [[06_matrices]] => this: **07 calculus-jacobians** => unlocks: [[08_nonlinear-least-squares]], [[09_manipulator-jacobian]], [[10_inverse-kinematics]].

## Why it matters for the stone scan->CAD->IFC / CADFit workflow

Almost nothing in a real stone pipeline has a closed-form answer. When you align a fresh scan to a previous one, fit a B-spline surface to a wall of irregular ashlar, calibrate the scanner camera, or recover a CAD program with CADFit, you are minimizing an error that depends nonlinearly on the unknowns (a pose, control points, intrinsics, extrusion lengths). The single tool that turns each of these nonlinear problems into something solvable is **linearization**: replace the curved error landscape near your current guess with a flat plane whose slope is the Jacobian, take one good step, and repeat. CADFit literally does this, fitting and validating parametric operations through an IoU-driven optimization that nudges parameters by repeated local steps. Once you understand gradients, Jacobians, and first-order Taylor expansion, every iterative solver in the pipeline (ICP, Gauss-Newton fitting, damped-least-squares inverse kinematics for a placing robot) stops being a black box and becomes the same loop with different residuals.

## Intuition

Imagine you are blindfolded on a hillside and want the lowest point. You cannot see the whole valley, but you can feel the slope under your feet. The slope tells you which way is downhill and how steep it is. You step a little, feel again, step again. That local slope is the **gradient**. The valley is the error function you want to minimize.

Now suppose instead of one altitude you are tracking several outputs at once: the x, y, and z error of every scan point. Each output has its own slope with respect to each input you can change. Stack all those slopes into a grid and you have the **Jacobian**: a matrix that says "if I wiggle input j a tiny bit, output i changes by this much." The Jacobian is the best flat approximation of a curved, multi-input multi-output map near your current point. Solving a hard nonlinear geometry problem is just: flatten with the Jacobian, solve the easy linear problem, move, and flatten again.

Analogy: a paper map is flat, but the Earth is curved. The map is a good linearization of the globe near your city and useless for the whole planet. The Jacobian is the local map of your problem, accurate close to where you stand and only there.

## The math (first principles)

### Conventions

Vectors are columns. $\mathbf{x} \in \mathbb{R}^n$ is the input (the unknowns), $\mathbf{f}: \mathbb{R}^n \to \mathbb{R}^m$ is a vector-valued function. Scalars are lowercase, vectors bold lowercase, matrices uppercase. The gradient is a column vector by convention here. Frames are right-handed unless stated. Lengths in this topic are unitless math; in the pipeline keep millimetres in shop space and metres in site space and never mix them inside one Jacobian, because the columns carry units and mismatched units wreck the conditioning (see Pitfalls).

### Derivative and partial derivative

For a scalar function of one variable, the derivative is the limit of the slope:

$$f'(x) = \lim_{h \to 0} \frac{f(x+h) - f(x)}{h}.$$

For a scalar function of several variables $f(x_1,\dots,x_n)$, the **partial derivative** $\partial f / \partial x_j$ is the ordinary derivative with all other variables held fixed:

$$\frac{\partial f}{\partial x_j} = \lim_{h \to 0} \frac{f(\dots, x_j + h, \dots) - f(\dots, x_j, \dots)}{h}.$$

### Gradient

The **gradient** collects all first partials of a scalar field into one column vector:

$$\nabla f(\mathbf{x}) = \left[ \frac{\partial f}{\partial x_1}, \frac{\partial f}{\partial x_2}, \dots, \frac{\partial f}{\partial x_n} \right]^{\top} \in \mathbb{R}^n.$$

It points in the direction of steepest increase, and its negative is the steepest-descent direction used by gradient descent. The first-order change in $f$ for a small step $\Delta\mathbf{x}$ is $\nabla f(\mathbf{x})^{\top} \Delta\mathbf{x}$.

### Jacobian of a vector-valued function

When the output is a vector $\mathbf{f}(\mathbf{x}) = [f_1(\mathbf{x}), \dots, f_m(\mathbf{x})]^{\top}$, the **Jacobian** $J$ is the $m \times n$ matrix of all first partials, one output per row, one input per column:

$$J(\mathbf{x}) = \frac{\partial \mathbf{f}}{\partial \mathbf{x}} = \begin{bmatrix}
\dfrac{\partial f_1}{\partial x_1} & \cdots & \dfrac{\partial f_1}{\partial x_n} \\
\vdots & \ddots & \vdots \\
\dfrac{\partial f_m}{\partial x_1} & \cdots & \dfrac{\partial f_m}{\partial x_n}
\end{bmatrix} \in \mathbb{R}^{m \times n}.$$

Row $i$ is $\nabla f_i^{\top}$. The Jacobian is the matrix that best linearly approximates $\mathbf{f}$ near $\mathbf{x}$. If $m = 1$ the Jacobian is a single row, the transpose of the gradient.

### Multivariate chain rule

If $\mathbf{f} = \mathbf{g}(\mathbf{h}(\mathbf{x}))$ with $\mathbf{h}: \mathbb{R}^n \to \mathbb{R}^k$ and $\mathbf{g}: \mathbb{R}^k \to \mathbb{R}^m$, the Jacobians multiply, outer first:

$$J_{\mathbf{f}}(\mathbf{x}) = J_{\mathbf{g}}\big(\mathbf{h}(\mathbf{x})\big)\, J_{\mathbf{h}}(\mathbf{x}).$$

This is matrix multiplication, and order matters because matrix multiplication does not commute. This rule is exactly how a camera projection Jacobian is built: world-to-camera transform, then perspective divide, then distortion, then pixel scaling, with each stage contributing one factor.

### First-order Taylor expansion (linearization)

For a scalar field, expand about $\mathbf{x}_0$:

$$f(\mathbf{x}_0 + \Delta\mathbf{x}) \approx f(\mathbf{x}_0) + \nabla f(\mathbf{x}_0)^{\top} \Delta\mathbf{x} + \tfrac{1}{2}\,\Delta\mathbf{x}^{\top} H(\mathbf{x}_0)\, \Delta\mathbf{x},$$

where $H$ is the **Hessian**, the symmetric $n \times n$ matrix of second partials $H_{ij} = \partial^2 f / (\partial x_i \partial x_j)$. Dropping the quadratic term leaves the **first-order (linear) approximation**, a tangent hyperplane.

For a vector field the first-order expansion uses the Jacobian:

$$\mathbf{f}(\mathbf{x}_0 + \Delta\mathbf{x}) \approx \mathbf{f}(\mathbf{x}_0) + J(\mathbf{x}_0)\, \Delta\mathbf{x}.$$

This is the workhorse equation of the whole topic. It says: near $\mathbf{x}_0$, the curved map $\mathbf{f}$ behaves like the affine map $\mathbf{x} \mapsto \mathbf{f}(\mathbf{x}_0) + J\,\Delta\mathbf{x}$. The approximation is good only for small $\Delta\mathbf{x}$, which is why solvers take small steps and re-linearize.

### From linearization to an iterative solver

A geometry problem is usually posed as nonlinear least squares: find $\mathbf{x}$ minimizing the sum of squared residuals

$$E(\mathbf{x}) = \tfrac{1}{2} \sum_{i} r_i(\mathbf{x})^2 = \tfrac{1}{2}\, \mathbf{r}(\mathbf{x})^{\top} \mathbf{r}(\mathbf{x}),$$

where each residual $r_i$ measures one mismatch (a point-to-point distance, a reprojection error, a pose error). Linearize the residual vector about the current guess $\mathbf{x}_k$ with $J = J_{\mathbf{r}}(\mathbf{x}_k)$:

$$\mathbf{r}(\mathbf{x}_k + \Delta\mathbf{x}) \approx \mathbf{r}(\mathbf{x}_k) + J\,\Delta\mathbf{x}.$$

Minimizing the squared length of that linear model over $\Delta\mathbf{x}$ is an ordinary linear least-squares problem. Setting its gradient to zero gives the **normal equations** and the **Gauss-Newton** step:

$$\big(J^{\top} J\big)\, \Delta\mathbf{x} = - J^{\top} \mathbf{r}(\mathbf{x}_k), \qquad \mathbf{x}_{k+1} = \mathbf{x}_k - \big(J^{\top} J\big)^{-1} J^{\top} \mathbf{r}(\mathbf{x}_k).$$

Gauss-Newton approximates the Hessian of $E$ by $J^{\top} J$ and so needs no second derivatives. When $J^{\top} J$ is ill-conditioned (near a singularity, or with too few constraints), add damping to get the **Levenberg-Marquardt** / **damped least squares** step:

$$\big(J^{\top} J + \lambda I\big)\, \Delta\mathbf{x} = - J^{\top} \mathbf{r}.$$

Small $\lambda$ behaves like Gauss-Newton (fast); large $\lambda$ behaves like gradient descent (safe). This single damped update is the same math used by camera calibration, B-spline fitting, ICP refinement, and the inverse-kinematics step that drives a placing robot. That is the conceptual bridge: closed-form solutions exist only for special linear cases, so the general nonlinear case is solved by **repeated linearization**.

## Worked example

Take a tiny two-residual fit so the arithmetic is checkable by hand. Model: predict a scalar $y$ from input $t$ with $\hat{y}(t) = a\,t + b$. Unknowns $\mathbf{x} = [a, b]^{\top}$. Two data points: $(t_1, y_1) = (1, 3)$ and $(t_2, y_2) = (2, 4)$. Residuals:

$$r_1 = a(1) + b - 3, \qquad r_2 = a(2) + b - 4.$$

Because the model is linear in $\mathbf{x}$, the Jacobian is constant:

$$J = \begin{bmatrix} \partial r_1/\partial a & \partial r_1/\partial b \\ \partial r_2/\partial a & \partial r_2/\partial b \end{bmatrix} = \begin{bmatrix} 1 & 1 \\ 2 & 1 \end{bmatrix}.$$

Start at the guess $\mathbf{x}_0 = [a, b]^{\top} = [0, 0]^{\top}$. Then $\mathbf{r}(\mathbf{x}_0) = [0 + 0 - 3,\; 0 + 0 - 4]^{\top} = [-3, -4]^{\top}$.

Build the normal-equation pieces:

$$J^{\top} J = \begin{bmatrix} 1 & 2 \\ 1 & 1 \end{bmatrix}\begin{bmatrix} 1 & 1 \\ 2 & 1 \end{bmatrix} = \begin{bmatrix} 5 & 3 \\ 3 & 2 \end{bmatrix}, \qquad J^{\top} \mathbf{r} = \begin{bmatrix} 1 & 2 \\ 1 & 1 \end{bmatrix}\begin{bmatrix} -3 \\ -4 \end{bmatrix} = \begin{bmatrix} -11 \\ -7 \end{bmatrix}.$$

Solve $\big(J^{\top} J\big)\,\Delta\mathbf{x} = -J^{\top}\mathbf{r} = [11, 7]^{\top}$. The determinant of $J^{\top}J$ is $5\cdot 2 - 3\cdot 3 = 1$, so the inverse is $\begin{bmatrix} 2 & -3 \\ -3 & 5 \end{bmatrix}$. Then

$$\Delta\mathbf{x} = \begin{bmatrix} 2 & -3 \\ -3 & 5 \end{bmatrix}\begin{bmatrix} 11 \\ 7 \end{bmatrix} = \begin{bmatrix} 22 - 21 \\ -33 + 35 \end{bmatrix} = \begin{bmatrix} 1 \\ 2 \end{bmatrix}.$$

So $\mathbf{x}_1 = \mathbf{x}_0 + \Delta\mathbf{x} = [1, 2]^{\top}$, that is $a = 1, b = 2$. Check: $\hat{y}(1) = 3$ (exact), $\hat{y}(2) = 4$ (exact). Both residuals are zero. Because this model is linear in the unknowns, Gauss-Newton lands on the exact answer in **one step**. For a nonlinear model the Jacobian would change each iteration and you would need several steps, but each step is this same computation.

## Code

Idiomatic C++ with Eigen first. The example builds a residual and its Jacobian, then takes Gauss-Newton steps. It also shows the numerical (finite-difference) Jacobian you reach for when an analytic one is painful, which is the everyday debugging tool.

```cpp
// gauss_newton.cpp  -- build with: g++ -I/path/to/eigen -O2 gauss_newton.cpp
#include <Eigen/Dense>
#include <functional>
#include <iostream>

using Eigen::VectorXd;
using Eigen::MatrixXd;

// Residual model: fit y = a*exp(b*t) to data. x = [a, b].
// This is nonlinear in b, so it needs several Gauss-Newton steps.
struct ExpFit {
    VectorXd t, y;                         // data
    // r_i(x) = a*exp(b*t_i) - y_i           (the residual vector r : R^2 -> R^m)
    VectorXd residual(const VectorXd& x) const {
        const double a = x(0), b = x(1);
        return (a * (b * t.array()).exp() - y.array()).matrix();
    }
    // Analytic Jacobian J = d r / d x   (m rows, 2 cols).
    // dr_i/da = exp(b t_i);   dr_i/db = a t_i exp(b t_i)
    MatrixXd jacobian(const VectorXd& x) const {
        const double a = x(0), b = x(1);
        const int m = static_cast<int>(t.size());
        MatrixXd J(m, 2);
        J.col(0) = (b * t.array()).exp().matrix();              // partial wrt a
        J.col(1) = (a * t.array() * (b * t.array()).exp()).matrix(); // partial wrt b
        return J;
    }
};

// Generic finite-difference Jacobian: column j = (r(x + h e_j) - r(x)) / h.
// Use this to CHECK an analytic Jacobian. If they disagree, the analytic one is wrong.
MatrixXd numericalJacobian(const std::function<VectorXd(const VectorXd&)>& r,
                           const VectorXd& x, double h = 1e-6) {
    VectorXd r0 = r(x);
    MatrixXd J(r0.size(), x.size());
    for (int j = 0; j < x.size(); ++j) {
        VectorXd xp = x; xp(j) += h;                 // perturb input j
        J.col(j) = (r(xp) - r0) / h;                 // slope of every output wrt input j
    }
    return J;
}

int main() {
    ExpFit p;
    p.t = VectorXd::LinSpaced(5, 0.0, 2.0);
    p.y = (2.0 * (0.5 * p.t.array()).exp()).matrix();   // truth: a=2, b=0.5

    VectorXd x(2); x << 1.0, 1.0;                       // initial guess (deliberately off)
    for (int it = 0; it < 20; ++it) {
        VectorXd r = p.residual(x);                     // r(x_k)
        MatrixXd J = p.jacobian(x);                     // J(x_k)
        // Solve the normal equations (J^T J) dx = -J^T r  for the Gauss-Newton step.
        VectorXd dx = (J.transpose() * J).ldlt().solve(-J.transpose() * r);
        x += dx;                                        // x_{k+1} = x_k + dx
        if (dx.norm() < 1e-10) break;                   // converged
    }
    std::cout << "fit a=" << x(0) << " b=" << x(1) << "\n"; // ~ a=2 b=0.5

    // Sanity check the analytic Jacobian against finite differences once.
    VectorXd x0(2); x0 << 1.5, 0.3;
    MatrixXd Ja = p.jacobian(x0);
    MatrixXd Jn = numericalJacobian([&](const VectorXd& z){ return p.residual(z); }, x0);
    std::cout << "max Jacobian mismatch: " << (Ja - Jn).cwiseAbs().maxCoeff() << "\n";
    return 0;
}
```

Domain note for the pipeline: when you call `cv::projectPoints` in OpenCV, the optional `jacobian` output is exactly the Jacobian of the projected pixel coordinates with respect to the rotation vector, translation, focal lengths, principal point, and distortion. Calibration feeds that Jacobian straight into a Levenberg-Marquardt loop. In OpenCASCADE, `GeomAPI_PointsToBSplineSurface` hides the linearized fitting loop, but the underlying minimization is the same normal-equation step. You rarely hand-write the camera Jacobian, but you must recognise it as the chain-rule product of the projection stages.

A short Python check clarifies the finite-difference idiom and nothing more:

```python
import numpy as np
def residual(x, t, y): a, b = x; return a*np.exp(b*t) - y
def num_jac(f, x, h=1e-6):
    r0 = f(x); J = np.zeros((r0.size, x.size))
    for j in range(x.size):
        xp = x.copy(); xp[j] += h
        J[:, j] = (f(xp) - r0) / h          # output slopes wrt input j
    return J
```

## Connections

- [[05_vectors]]: the gradient is a vector, and the dot product $\nabla f^{\top} \Delta\mathbf{x}$ is the first-order change; you need vector fluency before partials make sense.
- [[06_matrices]]: the Jacobian is a matrix, and the Gauss-Newton step is a linear solve $(J^{\top}J)\Delta\mathbf{x} = -J^{\top}\mathbf{r}$; conditioning and rank live in the matrix layer.
- [[08_nonlinear-least-squares]]: this topic supplies the residual, Jacobian, and linearized step that nonlinear least squares iterates; ICP, calibration, and surface fitting are instances.
- [[09_manipulator-jacobian]]: the manipulator Jacobian is the Jacobian of the forward-kinematics map from joint angles to end-effector pose, the same object applied to robot motion.
- [[10_inverse-kinematics]]: numerical IK is repeated linearization with a damped-least-squares step $\Delta\theta = J^{\top}(JJ^{\top} + \lambda^2 I)^{-1} e$, the robotics twin of Levenberg-Marquardt.
- [[03_svd-eigen]] (if present): SVD of the Jacobian reveals conditioning and singular directions, which is why damping is needed near singularities.

## Pitfalls and failure modes

- **Frame and units inconsistency.** A Jacobian column inherits the units of its input. Mixing metres and radians, or site-metres and shop-millimetres, in one parameter vector makes $J^{\top}J$ wildly anisotropic and slows or breaks convergence. Scale parameters to comparable magnitudes, or weight residuals, before solving.
- **Ill-conditioning and rank deficiency.** If two unknowns affect the residuals almost identically, $J^{\top}J$ is nearly singular, its inverse blows up, and the Gauss-Newton step explodes. This is the same disease as a robot near a kinematic singularity. The fix is damping ($\lambda I$) or removing the redundant parameter.
- **Stepping too far for the linear model.** First-order Taylor is local. If the residual is strongly curved and you take a full Gauss-Newton step, you can overshoot and diverge. Levenberg-Marquardt's adaptive $\lambda$, or a line search, keeps steps inside the region where the linearization is trustworthy.
- **Wrong or stale analytic Jacobian.** A sign error or a forgotten chain-rule factor produces a Jacobian that points the solver in a slightly wrong direction; it may still converge slowly and hide the bug. Always check the analytic Jacobian against a finite-difference Jacobian once, as in the code above.
- **Finite-difference step size.** Too large an $h$ adds truncation error; too small an $h$ drowns in floating-point cancellation. A step near $\sqrt{\epsilon_{\text{machine}}}$ ($\approx 10^{-6}$ to $10^{-8}$ for doubles) is the usual compromise. Central differences are more accurate but cost two evaluations per column.
- **Non-commutativity of the chain rule.** $J_{\mathbf{g}} J_{\mathbf{h}}$ is not $J_{\mathbf{h}} J_{\mathbf{g}}$. Build composite Jacobians outer-factor first, in the same order the functions apply.
- **Outliers.** Plain least squares squares the residual, so one bad scan point dominates the Jacobian-weighted normal equations. Use a robust loss (Huber, Tukey) or reject correspondences before the linear solve.
- **Parameterization of rotations.** Optimizing a $3\times3$ rotation matrix directly has nine entries but three degrees of freedom; the extra constraints corrupt the Jacobian. Parameterize rotations minimally (axis-angle / Lie-algebra increment) so the Jacobian has full, correct rank.

## Practice

1. By hand, compute the gradient of $f(x, y) = x^2 + 3xy + y^2$ at $(1, 2)$. Then confirm it numerically with a finite-difference snippet and print both. Expected gradient: $[2x + 3y,\; 3x + 2y] = [8, 7]$.
2. Write the Jacobian of $\mathbf{f}(x, y) = [x^2 - y,\; 2xy]^{\top}$ symbolically, evaluate it at $(1, 1)$, and verify against `numericalJacobian`. Inspect where the determinant of $J$ vanishes (a degeneracy locus).
3. Modify the C++ Gauss-Newton example to fit $y = a\sin(b t + c)$ to noisy samples. Print $\|\mathbf{r}\|$ each iteration and watch it decrease. Try a bad initial guess and observe divergence, then add Levenberg-Marquardt damping and watch it recover.
4. Take a $2\times2$ Jacobian whose columns are nearly parallel (e.g. $[[1, 1.001],[2, 2.002]]$). Print the condition number of $J^{\top}J$, attempt the undamped solve, then the damped solve, and compare the step magnitudes.
5. Build a tiny point-to-point registration: two paired 2D point sets related by a small rotation and translation, parameterize the pose as $[\theta, t_x, t_y]$, derive the $2N\times3$ residual Jacobian, and recover the pose by Gauss-Newton. Print the recovered pose against the truth.
6. (Stretch) For the projection chain world->camera->normalize->pixel, derive on paper the Jacobian of one pixel with respect to the rotation increment and translation, as the product of three stage Jacobians, then confirm the shape ($2\times6$) against `cv::projectPoints`.

## References

- Kevin Lynch and Frank Park, *Modern Robotics: Mechanics, Planning, and Control* (Cambridge University Press, 2017) and the companion site. Chapters on rigid-body motion, the manipulator Jacobian, and numerical inverse kinematics. https://hades.mech.northwestern.edu/index.php/Modern_Robotics
- OpenCV documentation, Camera Calibration and 3D Reconstruction (`calib3d`), including `projectPoints` and its Jacobian output, `calibrateCamera`, and `solvePnP`. https://docs.opencv.org/4.x/d9/d0c/group__calib3d.html
- Jorge Nocedal and Stephen Wright, *Numerical Optimization*, 2nd ed. (Springer, 2006). Gauss-Newton, Levenberg-Marquardt, Taylor expansion, and conditioning.
- K. Madsen, H. B. Nielsen, O. Tingleff, *Methods for Non-Linear Least Squares Problems*, 2nd ed. (Technical University of Denmark, 2004). Compact derivation of Gauss-Newton and Levenberg-Marquardt with Jacobians. https://www2.compute.dtu.dk/pubdb/pubs/3215-full.html
- Samuel Buss, *Introduction to Inverse Kinematics with Jacobian Transpose, Pseudoinverse and Damped Least Squares Methods* (2009). The damped-least-squares IK step $J^{\top}(JJ^{\top}+\lambda^2 I)^{-1}$.
- Eigen documentation, dense linear algebra (`ldlt`, `JacobiSVD`, array operations). https://eigen.tuxfamily.org/dox/
- Keenan Crane, CMU *Discrete Differential Geometry* course notes (gradients and Laplacians on meshes). https://www.cs.cmu.edu/~kmcrane/Projects/DDG/
- OpenCASCADE reference manual, `GeomAPI_PointsToBSplineSurface` (linearized surface fitting). https://dev.opencascade.org/doc/overview/html/
- Ghadi Nehme, Eamon Whalen, Faez Ahmed, *CADFit: Precise Mesh-to-CAD Program Generation with Hybrid Optimization*, arXiv:2605.01171 (2026). IoU-driven optimization over CAD programs. https://arxiv.org/abs/2605.01171 and code https://github.com/ghadinehme/CADFit

## Flashcard seeds

- Q: Build the Jacobian of $\mathbf{f}:\mathbb{R}^n\to\mathbb{R}^m$ column-by-column, starting from the limit definition of a partial derivative. <br> <b>Formula:</b> $\partial f_i/\partial x_j=\lim_{h\to 0}\frac{f_i(\dots,x_j+h,\dots)-f_i(\dots,x_j,\dots)}{h}$; target $J\in\mathbb{R}^{m\times n}$ <br> <i>Hint:</i> compute one partial $\partial f_i/\partial x_j$, then decide which index labels the row and which the column. :: A: for output $f_i$ the slope in direction $x_j$ is $\partial f_i/\partial x_j$ from the limit above <br> place outputs along rows ($i$) and inputs along columns ($j$), so entry $(i,j)=\partial f_i/\partial x_j$ <br> row $i$ is therefore $\nabla f_i^{\top}$ <br> $\boxed{J=\dfrac{\partial\mathbf{f}}{\partial\mathbf{x}}\in\mathbb{R}^{m\times n},\;J_{ij}=\partial f_i/\partial x_j}$.
- Q: Derive the first-order Taylor prediction of $\mathbf{f}$ at $\mathbf{x}_0+\Delta\mathbf{x}$ from the tangent picture (rise = slope $\times$ run). <br> <b>Formula:</b> rise $=J(\mathbf{x}_0)\,\Delta\mathbf{x}$; target $\mathbf{f}(\mathbf{x}_0+\Delta\mathbf{x})\approx\mathbf{f}(\mathbf{x}_0)+J\,\Delta\mathbf{x}$ <br> <i>Hint:</i> start from the known height $\mathbf{f}(\mathbf{x}_0)$ and add the rise along the tangent. :: A: the tangent through $(\mathbf{x}_0,\mathbf{f}(\mathbf{x}_0))$ has slope $J$, so the rise over run $\Delta\mathbf{x}$ is $J\,\Delta\mathbf{x}$ <br> add that rise to the starting height $\mathbf{f}(\mathbf{x}_0)$ <br> the dropped piece is the curvature term $O(\|\Delta\mathbf{x}\|^2)$ <br> $\boxed{\mathbf{f}(\mathbf{x}_0+\Delta\mathbf{x})\approx\mathbf{f}(\mathbf{x}_0)+J(\mathbf{x}_0)\,\Delta\mathbf{x}}$.
- Q: Reduce the full second-order Taylor expansion of a scalar field to its first-order (linear) approximation (step 1 of 2). <br> <b>Formula:</b> $f(\mathbf{x}_0+\Delta\mathbf{x})=f(\mathbf{x}_0)+\nabla f^{\top}\Delta\mathbf{x}+\tfrac12\Delta\mathbf{x}^{\top}H\,\Delta\mathbf{x}+\dots$ <br> <i>Hint:</i> compare how each term scales in $\|\Delta\mathbf{x}\|$ and drop the one that vanishes faster. :: A: the linear term scales as $\|\Delta\mathbf{x}\|$, the quadratic term as $\|\Delta\mathbf{x}\|^2$ <br> for small steps the quadratic term is negligible, so drop it <br> $\boxed{f(\mathbf{x}_0+\Delta\mathbf{x})\approx f(\mathbf{x}_0)+\nabla f(\mathbf{x}_0)^{\top}\Delta\mathbf{x}}$.
- Q: From the scalar linearization, derive the steepest-descent direction (step 2 of 2). <br> <b>Formula:</b> $\Delta f\approx\nabla f^{\top}\Delta\mathbf{x}=\|\nabla f\|\,\|\Delta\mathbf{x}\|\cos\theta$ <br> <i>Hint:</i> hold $\|\Delta\mathbf{x}\|$ fixed and ask which $\cos\theta$ makes $\Delta f$ most negative. :: A: write $\Delta f\approx\|\nabla f\|\,\|\Delta\mathbf{x}\|\cos\theta$ where $\theta$ is the angle between $\Delta\mathbf{x}$ and $\nabla f$ <br> for fixed step length this is most negative at $\cos\theta=-1$, i.e. $\Delta\mathbf{x}$ antiparallel to $\nabla f$ <br> $\boxed{\Delta\mathbf{x}\propto-\nabla f(\mathbf{x}_0)}$.
- Q: Set up the linear least-squares model whose minimizer is the Gauss-Newton step (step 1 of 2). <br> <b>Formula:</b> $E(\mathbf{x})=\tfrac12\,\mathbf{r}(\mathbf{x})^{\top}\mathbf{r}(\mathbf{x})$ and $\mathbf{r}(\mathbf{x}_k+\Delta\mathbf{x})\approx\mathbf{r}(\mathbf{x}_k)+J\,\Delta\mathbf{x}$ <br> <i>Hint:</i> substitute the linearized residual into the cost to get a quadratic in $\Delta\mathbf{x}$. :: A: cost is $E=\tfrac12\,\mathbf{r}^{\top}\mathbf{r}$ <br> linearize the residual about $\mathbf{x}_k$ with $J=J_{\mathbf{r}}(\mathbf{x}_k)$ <br> substitute the linear model in place of $\mathbf{r}$ <br> $\boxed{E\approx\tfrac12\|\mathbf{r}(\mathbf{x}_k)+J\,\Delta\mathbf{x}\|^2}$.
- Q: From the quadratic model $E\approx\tfrac12\|\mathbf{r}+J\,\Delta\mathbf{x}\|^2$, derive the normal equations and the Gauss-Newton step (step 2 of 2). <br> <b>Formula:</b> $E=\tfrac12(\mathbf{r}^{\top}\mathbf{r}+2\,\mathbf{r}^{\top}J\Delta\mathbf{x}+\Delta\mathbf{x}^{\top}J^{\top}J\Delta\mathbf{x})$ <br> <i>Hint:</i> set the gradient $\partial E/\partial\Delta\mathbf{x}$ to zero and solve for $\Delta\mathbf{x}$. :: A: expand the squared norm as shown <br> set $\partial E/\partial\Delta\mathbf{x}=J^{\top}\mathbf{r}+J^{\top}J\,\Delta\mathbf{x}=\mathbf{0}$ <br> so $(J^{\top}J)\,\Delta\mathbf{x}=-J^{\top}\mathbf{r}$ and $\Delta\mathbf{x}=-(J^{\top}J)^{-1}J^{\top}\mathbf{r}$ <br> $\boxed{\mathbf{x}_{k+1}=\mathbf{x}_k-(J^{\top}J)^{-1}J^{\top}\mathbf{r}(\mathbf{x}_k)}$.
- Q: From the Gauss-Newton normal equations, derive the Levenberg-Marquardt (damped) step and state its two limits. <br> <b>Formula:</b> start from $(J^{\top}J)\,\Delta\mathbf{x}=-J^{\top}\mathbf{r}$; target $(J^{\top}J+\lambda I)\,\Delta\mathbf{x}=-J^{\top}\mathbf{r}$ <br> <i>Hint:</i> add $\lambda I$ to the curvature matrix, then take $\lambda\to0$ and $\lambda\to\infty$. :: A: add damping $\lambda I$ to $J^{\top}J$ to keep it invertible: $(J^{\top}J+\lambda I)\,\Delta\mathbf{x}=-J^{\top}\mathbf{r}$ <br> $\lambda\to0$ recovers Gauss-Newton (fast) <br> $\lambda\to\infty$ gives $\Delta\mathbf{x}\approx-\tfrac1\lambda J^{\top}\mathbf{r}$, a scaled gradient-descent step (safe) <br> $\boxed{(J^{\top}J+\lambda I)\,\Delta\mathbf{x}=-J^{\top}\mathbf{r}}$.
- Q: Derive the damped-least-squares inverse-kinematics step and say why $\lambda$ is essential near a singularity. <br> <b>Formula:</b> minimize $\|J\,\Delta\theta-e\|^2+\lambda^2\|\Delta\theta\|^2$; target $\Delta\theta=J^{\top}(JJ^{\top}+\lambda^2 I)^{-1}e$ <br> <i>Hint:</i> set the gradient of the damped objective to zero, then rewrite the left-pseudoinverse form as the right-pseudoinverse form. :: A: task error $e\approx J\,\Delta\theta$; minimizing the damped objective gives $(J^{\top}J+\lambda^2 I)\Delta\theta=J^{\top}e$ <br> the algebra rearranges to the equivalent right-pseudoinverse form $\Delta\theta=J^{\top}(JJ^{\top}+\lambda^2 I)^{-1}e$ <br> near a singularity small singular values of $J$ would blow up the pure pseudoinverse; $\lambda^2$ floors them <br> $\boxed{\Delta\theta=J^{\top}(JJ^{\top}+\lambda^2 I)^{-1}e}$.
- Q: Derive the constant Jacobian of the worked line-fit. <br> <b>Formula:</b> residuals $r_1=a\cdot1+b-3$, $r_2=a\cdot2+b-4$, unknowns $\mathbf{x}=[a,b]^{\top}$; $J_{ij}=\partial r_i/\partial x_j$ <br> <i>Hint:</i> differentiate each residual w.r.t. $a$ then $b$; the $t$-coefficient becomes the $a$-column. :: A: $\partial r_1/\partial a=1,\ \partial r_1/\partial b=1$; $\partial r_2/\partial a=2,\ \partial r_2/\partial b=1$ <br> the model is linear in $\mathbf{x}$, so these are constants (no dependence on $a,b$) <br> $\boxed{J=\begin{bmatrix}1&1\\2&1\end{bmatrix}}$.
- Q: Worked example (line fit, one Gauss-Newton step). Compute $\Delta\mathbf{x}$ and $\mathbf{x}_1$ given $J=\begin{bmatrix}1&1\\2&1\end{bmatrix}$, $\mathbf{x}_0=[0,0]^{\top}$, $\mathbf{r}(\mathbf{x}_0)=[-3,-4]^{\top}$. <br> <b>Formula:</b> $(J^{\top}J)\,\Delta\mathbf{x}=-J^{\top}\mathbf{r}$, $\mathbf{x}_1=\mathbf{x}_0+\Delta\mathbf{x}$ <br> <i>Hint:</i> form $J^{\top}J$ and $J^{\top}\mathbf{r}$ first, then invert the $2\times2$ via $\det$. :: A: $J^{\top}J=\begin{bmatrix}5&3\\3&2\end{bmatrix}$, $J^{\top}\mathbf{r}=[-11,-7]^{\top}$ <br> $\det(J^{\top}J)=10-9=1$, so $(J^{\top}J)^{-1}=\begin{bmatrix}2&-3\\-3&5\end{bmatrix}$ <br> $\Delta\mathbf{x}=\begin{bmatrix}2&-3\\-3&5\end{bmatrix}\begin{bmatrix}11\\7\end{bmatrix}=[1,2]^{\top}$ <br> $\boxed{\mathbf{x}_1=[1,2]^{\top}\ (a=1,b=2)}$.
- Q: Closing test (line fit). You found $a=1,b=2$ in one step; verify the residuals are zero and explain why one step sufficed. <br> <b>Formula:</b> $\hat y(t)=at+b$, residuals $r_i=\hat y(t_i)-y_i$ at $(t_1,y_1)=(1,3),(t_2,y_2)=(2,4)$ <br> <i>Hint:</i> a pass means both residuals evaluate to $0$; then recall $J$ is constant for a linear model. :: A: $\hat y(1)=1\cdot1+2=3=y_1$ and $\hat y(2)=1\cdot2+2=4=y_2$, so $\mathbf{r}(\mathbf{x}_1)=\mathbf{0}$ <br> the model is linear in $\mathbf{x}$, so $J$ is constant and the linearization is exact <br> $\boxed{\text{Gauss-Newton is exact in one step for linear-in-parameter models}}$.
- Q: Worked example (small scale, mm shop space). Take one Gauss-Newton step on $r(a)=a^2-2$ from $a_0=2$. <br> <b>Formula:</b> $\Delta a=-(J^{\top}J)^{-1}J^{\top}r$ with $J=dr/da=2a$ <br> <i>Hint:</i> compute $r(2)$ and $J=2a$ first, then the scalar step is just $-r/J$. :: A: $r(2)=4-2=2$, $J=2a=4$ <br> $\Delta a=-\tfrac{1}{16}\cdot4\cdot2=-\tfrac{r}{J}=-0.5$ <br> $\boxed{a_1=2-0.5=1.5}$ (true root $\sqrt2\approx1.414$, one step closes most of the gap).
- Q: Closing test (root fit). Continue the $r(a)=a^2-2$ fit one more step from $a_1=1.5$ and check it nears $\sqrt2$. <br> <b>Formula:</b> $\Delta a=-(J^{\top}J)^{-1}J^{\top}r$, $J=2a$, $\sqrt2=1.41421$ <br> <i>Hint:</i> a pass means $|a_2-\sqrt2|$ shrinks well below the previous error $\approx0.086$. :: A: $r(1.5)=2.25-2=0.25$, $J=2(1.5)=3$ <br> $\Delta a=-\tfrac{1}{9}\cdot3\cdot0.25=-0.0833$ <br> $a_2=1.5-0.0833=1.4167$ vs $\sqrt2=1.41421$ <br> $\boxed{|a_2-\sqrt2|\approx2.5\times10^{-3}\ \text{— error shrank}\sim34\times,\text{ converging}}$.
- Q: Worked example (finite-difference step at two scales). Pick the FD step $h$ for a large input $\|x_j\|\!\sim\!10^3$ and a small one $\|x_j\|\!\sim\!1$ in doubles. <br> <b>Formula:</b> $h\approx\sqrt{\epsilon_{\text{mach}}}\,(1+\|x_j\|)\approx10^{-8}(1+\|x_j\|)$ <br> <i>Hint:</i> plug each magnitude into the formula; the $(1+\|x_j\|)$ factor scales $h$ to the input. :: A: large: $h\approx10^{-8}(1+10^3)\approx10^{-5}$ <br> small: $h\approx10^{-8}(1+1)\approx2\times10^{-8}$ <br> $\boxed{\text{scale }h\text{ relative to the input magnitude, never a fixed }10^{-6}}$.
- Q: Worked example (robot IK, task space, metres/radians). Compute the joint step $\Delta\theta$ given $e=[0.02,\,0]^{\top}$ m, $JJ^{\top}+\lambda^2 I=\begin{bmatrix}0.04&0\\0&0.04\end{bmatrix}$, $J^{\top}=\begin{bmatrix}0.1&0\\0&0.1\end{bmatrix}$. <br> <b>Formula:</b> $\Delta\theta=J^{\top}(JJ^{\top}+\lambda^2 I)^{-1}e$ <br> <i>Hint:</i> invert the diagonal $2\times2$ first (reciprocate the diagonal), then apply it to $e$, then multiply by $J^{\top}$. :: A: $(JJ^{\top}+\lambda^2 I)^{-1}=\begin{bmatrix}25&0\\0&25\end{bmatrix}$ <br> $(JJ^{\top}+\lambda^2 I)^{-1}e=[0.5,\,0]^{\top}$ <br> $\Delta\theta=J^{\top}[0.5,\,0]^{\top}=[0.05,\,0]^{\top}$ <br> $\boxed{\Delta\theta=[0.05,\,0]^{\top}\ \text{rad}}$.
- Q: Worked example (gradient, 2D, by hand). Compute $\nabla f$ of $f(x,y)=x^2+3xy+y^2$ at $(1,2)$. <br> <b>Formula:</b> $\nabla f=[\partial f/\partial x,\ \partial f/\partial y]^{\top}=[2x+3y,\ 3x+2y]^{\top}$ <br> <i>Hint:</i> compute each partial symbolically first, then substitute $x=1,y=2$. :: A: $\partial f/\partial x=2x+3y$, $\partial f/\partial y=3x+2y$ <br> at $(1,2)$: $2(1)+3(2)=8$ and $3(1)+2(2)=7$ <br> $\boxed{\nabla f(1,2)=[8,\,7]^{\top}}$.
- Q: Worked example (Jacobian of a 2D map). Write the Jacobian of $\mathbf{f}(x,y)=[x^2-y,\ 2xy]^{\top}$ symbolically and evaluate at $(1,1)$. <br> <b>Formula:</b> $J=\begin{bmatrix}\partial f_1/\partial x&\partial f_1/\partial y\\\partial f_2/\partial x&\partial f_2/\partial y\end{bmatrix}$ <br> <i>Hint:</i> fill the four partials of $f_1=x^2-y$ and $f_2=2xy$ first, then substitute $x=1,y=1$. :: A: $\partial f_1/\partial x=2x,\ \partial f_1/\partial y=-1$; $\partial f_2/\partial x=2y,\ \partial f_2/\partial y=2x$ <br> so $J=\begin{bmatrix}2x&-1\\2y&2x\end{bmatrix}$ <br> at $(1,1)$: $\boxed{J=\begin{bmatrix}2&-1\\2&2\end{bmatrix}}$.
- Q: Closing test (verify an analytic Jacobian against finite differences). For $\mathbf{f}(x,y)=[x^2-y,\ 2xy]^{\top}$ at $(1,1)$ with $h=10^{-3}$, check the first column. <br> <b>Formula:</b> FD column $j=(\mathbf{f}(\mathbf{x}+h\mathbf{e}_j)-\mathbf{f}(\mathbf{x}))/h$; analytic first column $=[2x,\ 2y]^{\top}=[2,2]^{\top}$ <br> <i>Hint:</i> a pass means the FD column matches $[2,2]^{\top}$ to about $O(h)$; perturb only $x$. :: A: $\mathbf{f}(1,1)=[0,\,2]^{\top}$; $\mathbf{f}(1+h,1)=[(1.001)^2-1,\ 2(1.001)]=[0.002001,\,2.002]^{\top}$ <br> FD col $=([0.002001,2.002]^{\top}-[0,2]^{\top})/10^{-3}=[2.001,\,2.000]^{\top}$ <br> analytic col $=[2,\,2]^{\top}$, mismatch $\approx10^{-3}=O(h)$ <br> $\boxed{\text{match within }O(h)\Rightarrow\text{analytic Jacobian verified}}$.
- Q: Closing test (conditioning / when to damp). For $J=\begin{bmatrix}1&1.001\\2&2.002\end{bmatrix}$, decide whether to use the plain or damped solve. <br> <b>Formula:</b> step explodes when $\det(J^{\top}J)\approx0$; damped solve uses $(J^{\top}J+\lambda I)\,\Delta\mathbf{x}=-J^{\top}\mathbf{r}$ <br> <i>Hint:</i> a "near-singular" verdict means $\det(J)\approx0$; compute $\det(J)=1\cdot2.002-1.001\cdot2$. :: A: $\det(J)=2.002-2.002=0$ exactly (columns parallel), so $J^{\top}J$ is singular <br> the undamped step needs $(J^{\top}J)^{-1}$, which does not exist, so it blows up <br> add $\lambda I$ to floor the small eigenvalue and get a finite step <br> $\boxed{\det\approx0\Rightarrow\text{use the damped (Levenberg-Marquardt) solve}}$.
- Q: What is the gradient of a scalar field $f:\mathbb{R}^n\to\mathbb{R}$, and how does it relate to a row of the Jacobian? <br> <b>Formula:</b> $\nabla f=[\partial f/\partial x_1,\dots,\partial f/\partial x_n]^{\top}$ <br> <i>Hint:</i> think about the special case $m=1$. :: A: $\nabla f$ is a column vector of all first partials <br> for $m=1$ the Jacobian has a single row, which is exactly $\nabla f^{\top}$ <br> $\boxed{J=\nabla f^{\top}\ \text{when }m=1}$.
- Q: When should you reach for Levenberg-Marquardt instead of plain Gauss-Newton? <br> <b>Formula:</b> $(J^{\top}J+\lambda I)\,\Delta\mathbf{x}=-J^{\top}\mathbf{r}$ <br> <i>Hint:</i> think about what makes $(J^{\top}J)^{-1}$ blow up or the linear model untrustworthy. :: A: when $J^{\top}J$ is ill-conditioned or rank-deficient (near a singularity, too few constraints), or when a full step overshoots a strongly curved residual <br> the $\lambda I$ damping keeps the step inside the trustworthy region <br> $\boxed{\text{damp when }J^{\top}J\text{ is near-singular or the step overshoots}}$.
- Q: Pitfall: why must you never mix metres and radians (or site-m and shop-mm) inside one parameter vector? <br> <b>Formula:</b> each Jacobian column inherits the units of its input $x_j$ <br> <i>Hint:</i> think about what mixed-unit columns do to the scaling (anisotropy) of $J^{\top}J$. :: A: a column with metres and a column with radians have wildly different magnitudes <br> this makes $J^{\top}J$ strongly anisotropic, wrecking conditioning and convergence <br> $\boxed{\text{scale parameters to comparable magnitudes (or weight residuals) before solving}}$.
- Q: Pitfall: what is the everyday check for a hand-derived analytic Jacobian, and what does failure look like? <br> <b>Formula:</b> FD Jacobian column $=(\mathbf{r}(\mathbf{x}+h\mathbf{e}_j)-\mathbf{r}(\mathbf{x}))/h$ <br> <i>Hint:</i> compare analytic vs FD once; a pass is a tiny max entrywise difference. :: A: build the finite-difference Jacobian and compare it once to the analytic one <br> a sign error or missing chain-rule factor shows as a large max entrywise mismatch, even if the solver still limps to an answer <br> $\boxed{\text{large }|J_{\text{analytic}}-J_{\text{FD}}|_{\max}\Rightarrow\text{the analytic Jacobian is wrong}}$.
- CLOZE: The Jacobian is the best {{c1::linear (first-order)}} approximation of a vector-valued map near a point, accurate only for {{c2::small}} steps.
- CLOZE: Gauss-Newton approximates the cost Hessian by {{c1::$J^{\top}J$}}, which is why it needs no {{c2::second}} derivatives.
- CLOZE: Setting $\partial E/\partial\Delta\mathbf{x}=0$ in the linear model yields the normal equations {{c1::$(J^{\top}J)\Delta\mathbf{x}=-J^{\top}\mathbf{r}$}}. <br> <i>hint:</i> curvature matrix times step equals minus $J^{\top}\mathbf{r}$.
- CLOZE: The matrix chain rule $J_{\mathbf{f}}=J_{\mathbf{g}}J_{\mathbf{h}}$ does not {{c1::commute}}, so build composite Jacobians {{c2::outer-factor}} first.
- CLOZE: A finite-difference Jacobian column is $(\mathbf{r}(\mathbf{x}+h\mathbf{e}_j)-\mathbf{r}(\mathbf{x}))/h$, with $h\approx${{c1::$\sqrt{\epsilon_{\text{machine}}}$}} for doubles.
- CLOZE: Optimizing a full $3\times3$ rotation matrix directly is wrong because it has nine entries but only {{c1::three}} degrees of freedom; parameterize with an {{c2::axis-angle / Lie-algebra}} increment instead.
- Q: In Eigen, how do you solve the Gauss-Newton normal equations for the step, and why `ldlt`? <br> <b>Formula:</b> $(J^{\top}J)\,\Delta\mathbf{x}=-J^{\top}\mathbf{r}$ <br> <i>Hint:</i> the call factors the left-hand matrix; recall what kind of matrix $J^{\top}J$ always is. :: A: `(J.transpose()*J).ldlt().solve(-J.transpose()*r)` <br> $J^{\top}J$ is symmetric positive-semidefinite, so the $LDL^{\top}$ factorization is the right, cheap solver <br> $\boxed{\text{use }LDL^{\top}\text{ because }J^{\top}J\text{ is symmetric PSD}}$.
- Q: What does `cv::projectPoints` return in its optional `jacobian` argument, and which chain-rule stages build it? <br> <b>Formula:</b> $J_{\mathbf{f}}=J_{\mathbf{g}}\,J_{\mathbf{h}}$ (Jacobians multiply, outer factor first) <br> <i>Hint:</i> list the projection stages in order and assign one Jacobian factor to each. :: A: it returns partials of the projected pixels w.r.t. rotation, translation, focal lengths, principal point, and distortion <br> it is the product of the world-to-camera, perspective-divide, distortion, and pixel-scaling stage Jacobians <br> $\boxed{J_{\text{pixel}}=J_{\text{scale}}J_{\text{distort}}J_{\text{divide}}J_{\text{w2c}}}$.
