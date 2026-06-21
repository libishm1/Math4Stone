# Optimization (foundations): gradient descent, convexity, and constraints

Placement: [[31_calculus-foundations]], [[06_linear-least-squares]] => **Optimization foundations** => [[08_nonlinear-least-squares]], [[15_inverse-kinematics]], [[27_supervised-learning-optimization]], [[25_cadfit-program-synthesis]]

## Why it matters for the stone scan->CAD->IFC / CADFit / robotic-fabrication workflow

Almost everything in this pipeline is an optimization problem wearing a domain costume. Camera calibration minimizes reprojection error. ICP minimizes alignment error. Inverse kinematics minimizes pose error. Surface fitting minimizes geometric deviation. CADFit searches CAD programs that minimize $-\mathrm{IoU}$ plus fabrication and assembly penalties. Once you see the shared structure "find $x$ that minimizes $f(x)$, maybe subject to constraints," you stop learning each algorithm as a separate trick and start recognizing one method specialized many ways. This node is the trunk that [[08_nonlinear-least-squares]], [[15_inverse-kinematics]], and [[27_supervised-learning-optimization]] branch from.

## Intuition

Optimization is descending a landscape to its lowest valley. The gradient is the uphill arrow; you step against it. If the landscape is a single smooth bowl (convex), any downhill path reaches the global bottom. If it is a mountain range (non-convex), where you start decides which valley you land in. Constraints are fences: you must find the lowest point that stays inside the fenced region, and the solution often presses right against a fence.

## The math (first principles)

**The problem.** Unconstrained: $\min_{x\in\mathbb{R}^n} f(x)$. The first-order optimality (stationarity) condition is

$$ \nabla f(x^\star) = 0. $$

A stationary point is a minimum if the Hessian $H=\nabla^2 f$ is positive semidefinite there.

**Gradient descent.** Step opposite the gradient with step size (learning rate) $\alpha$:

$$ x_{k+1} = x_k - \alpha\,\nabla f(x_k). $$

Too large $\alpha$ diverges; too small crawls. A backtracking **line search** picks $\alpha$ by shrinking it until $f$ decreases enough (Armijo condition). **Momentum** averages past gradients to speed up narrow valleys: $v_{k+1}=\beta v_k - \alpha\nabla f$, $x_{k+1}=x_k+v_{k+1}$.

**Newton's method.** Use curvature. From the second-order Taylor model, the step that zeros the model gradient is

$$ x_{k+1} = x_k - H(x_k)^{-1}\nabla f(x_k). $$

Newton converges quadratically near a minimum but needs the Hessian and a positive-definite $H$. **Gauss-Newton** approximates $H \approx J^\top J$ for least-squares objectives, and **Levenberg-Marquardt** damps it; both live in [[08_nonlinear-least-squares]].

**Convexity.** A set is convex if the segment between any two of its points stays inside it. A function is convex if $f(\lambda a + (1-\lambda)b) \le \lambda f(a)+(1-\lambda)f(b)$, equivalently if $H \succeq 0$ everywhere. For convex problems every local minimum is global, so gradient descent reaches the true optimum. Linear least squares is convex; neural-network training and CADFit search are not.

**Constrained optimization.** With equality constraints $g(x)=0$, form the Lagrangian

$$ \mathcal{L}(x,\lambda) = f(x) + \lambda^\top g(x), $$

and solve $\nabla_x\mathcal{L}=0,\ \nabla_\lambda\mathcal{L}=0$. Geometrically, at the optimum the cost gradient is parallel to the constraint gradient: $\nabla f = -\lambda^\top \nabla g$. With inequality constraints $h(x)\le 0$ the **KKT conditions** hold: stationarity $\nabla f + \lambda^\top\nabla g + \mu^\top\nabla h = 0$, primal feasibility, dual feasibility $\mu\ge 0$, and complementary slackness $\mu_i h_i(x)=0$ (a constraint is either inactive or its multiplier presses on it).

**The unifying view.** Calibration, ICP, IK, nonlinear least squares, and ML training all instantiate $\min_x f(x)$ with different residuals and constraints. Pick the solver by the structure: convex and linear -> closed form; smooth and non-convex -> gradient/Newton from a good start; least-squares residuals -> Gauss-Newton/LM; discrete/combinatorial (CADFit op sequences) -> search plus a continuous refine.

## Worked example

Minimize $f(x) = (x-3)^2$. The gradient is $f'(x) = 2(x-3)$; setting it to zero gives $x^\star = 3$. One gradient-descent step from $x_0 = 0$ with $\alpha = 0.1$: $x_1 = 0 - 0.1\cdot 2(0-3) = 0 + 0.6 = 0.6$. Next: $x_2 = 0.6 - 0.1\cdot 2(0.6-3) = 0.6 + 0.48 = 1.08$. The iterates march toward 3. Newton's method does it in one step: $H = f'' = 2$, so $x_1 = 0 - 2^{-1}\cdot(2)(0-3) = 0 + 3 = 3$. A quadratic is solved exactly by one Newton step.

Constrained: minimize $f(x,y)=x^2+y^2$ subject to $x+y=1$. Lagrangian $\mathcal{L}=x^2+y^2+\lambda(x+y-1)$. Stationarity: $2x+\lambda=0,\ 2y+\lambda=0 \Rightarrow x=y$. Feasibility $x+y=1 \Rightarrow x=y=1/2$. The closest point on the line to the origin is $(0.5,0.5)$, as expected.

## Code

Idiomatic C++ with Eigen: gradient descent with backtracking line search, and a Newton step. Both take an objective and its gradient (gradient-check it against [[31_calculus-foundations]]'s finite difference first).

```cpp
#include <Eigen/Dense>
#include <functional>

using Vec = Eigen::VectorXd;
using Mat = Eigen::MatrixXd;

// Gradient descent with Armijo backtracking line search.
Vec gradientDescent(const std::function<double(const Vec&)>& f,
                    const std::function<Vec(const Vec&)>& grad,
                    Vec x, int iters = 500, double tol = 1e-8) {
    const double c = 1e-4, shrink = 0.5;       // Armijo params
    for (int k = 0; k < iters; ++k) {
        Vec g = grad(x);
        if (g.norm() < tol) break;             // stationarity: ||grad f|| ~ 0
        double alpha = 1.0, fx = f(x);
        Vec dir = -g;                          // steepest-descent direction
        while (f(x + alpha * dir) > fx + c * alpha * g.dot(dir))
            alpha *= shrink;                   // shrink until sufficient decrease
        x += alpha * dir;
    }
    return x;
}

// One Newton step: x <- x - H^{-1} grad.  Use ldlt for an SPD Hessian.
Vec newtonStep(const Vec& x, const Vec& g, const Mat& H) {
    return x - H.ldlt().solve(g);
}
```

## Connections

- [[31_calculus-foundations]] gives the gradient and Hessian; "set $\nabla f=0$" is the optimality condition.
- [[06_linear-least-squares]] is the convex, closed-form special case (normal equations).
- [[08_nonlinear-least-squares]] specializes Newton to least-squares residuals via Gauss-Newton and Levenberg-Marquardt.
- [[15_inverse-kinematics]] solves a nonlinear optimization for joint angles; damped least squares is regularized Gauss-Newton.
- [[27_supervised-learning-optimization]] trains models by gradient descent on a loss; same update rule, many parameters.
- [[25_cadfit-program-synthesis]] is a discrete search (op sequences) wrapped around continuous parameter optimization.

## Pitfalls and failure modes

- **Step size.** Fixed $\alpha$ is fragile; prefer a line search. Divergence usually means $\alpha$ too big.
- **Non-convexity and initialization.** Gradient descent finds *a* local minimum, not the global one. ICP and IK fail from bad starts for this reason; warm-start them.
- **Ill-conditioning.** A long narrow valley (large Hessian condition number) makes plain gradient descent zig-zag; use momentum, preconditioning, or Newton.
- **Hessian cost and definiteness.** Newton needs $H$ and a positive-definite $H$; far from the optimum $H$ can be indefinite. Damping (trust region / Levenberg-Marquardt) fixes this.
- **Scaling.** Variables in different units (metres vs radians) wreck conditioning; normalize them.
- **Stopping criteria.** Stop on small gradient norm and/or small step, not just iteration count.

## Practice

1. Run `gradientDescent` on $f(x)=(x-3)^2$ from $x_0=0$; log $\alpha$ and confirm convergence to 3. Then remove the line search, set $\alpha=1.1$, and watch it diverge.
2. Minimize the Rosenbrock function $f(x,y)=(1-x)^2+100(y-x^2)^2$ with gradient descent and with Newton; compare iteration counts.
3. Solve $\min x^2+y^2$ s.t. $x+y=1$ with a Lagrange multiplier by hand, then verify numerically with a penalty method ($+\rho(x+y-1)^2$, increasing $\rho$).
4. Show that one Newton step solves any quadratic exactly, and explain why.
5. Take an ICP-like residual with two basins; run gradient descent from several starts and tabulate which basin each reaches (demonstrates non-convexity).

## References

- Stephen Boyd and Lieven Vandenberghe, *Convex Optimization* (free): https://web.stanford.edu/~boyd/cvxbook/
- Jorge Nocedal and Stephen Wright, *Numerical Optimization* (line search, Newton, trust region, KKT).
- Sebastian Ruder, "An overview of gradient descent optimization algorithms": https://www.ruder.io/optimizing-gradient-descent/
- Lynch and Park, *Modern Robotics* (numerical IK as optimization).

## Flashcard seeds

- Q: Derive the first-order optimality (stationarity) condition for $\min_x f(x)$, starting from the first-order Taylor expansion. <br> <b>Formula:</b> start from $f(x^\star+\delta)\approx f(x^\star)+\nabla f(x^\star)^\top\delta$; target $\nabla f(x^\star)=0$. <br> <i>Hint:</i> at a minimum no direction $\delta$ may lower $f$, so the linear term $\nabla f(x^\star)^\top\delta$ must vanish for every $\delta$. :: A: Near $x^\star$, $f(x^\star+\delta)\approx f(x^\star)+\nabla f(x^\star)^\top\delta$ <br> For $x^\star$ to be a min, no direction $\delta$ may lower $f$, so the linear term must vanish for all $\delta$ <br> $\nabla f(x^\star)^\top\delta=0\ \forall\delta$ forces the vector itself to be zero <br> $\boxed{\nabla f(x^\star)=0}$
- Q: Derive the gradient-descent direction (step 1 of 2): which direction $d$ gives the steepest local decrease of $f$ at $x_k$? <br> <b>Formula:</b> $f(x_k+d)\approx f(x_k)+\nabla f(x_k)^\top d$; target $d=-\nabla f(x_k)$. <br> <i>Hint:</i> minimize $\nabla f(x_k)^\top d$ over unit-length $d$; by Cauchy-Schwarz it is most negative when $d$ points opposite the gradient. :: A: $f(x_k+d)\approx f(x_k)+\nabla f(x_k)^\top d$ for small $d$ <br> minimize the change $\nabla f(x_k)^\top d$ over unit $d$; by Cauchy-Schwarz it is most negative when $d\parallel-\nabla f(x_k)$ <br> $\boxed{d=-\nabla f(x_k)}$ (steepest-descent direction)
- Q: Continue from the steepest-descent direction to the gradient-descent update (step 2 of 2). <br> <b>Formula:</b> $x_{k+1}=x_k+\alpha d$ with $d=-\nabla f(x_k)$; target $x_{k+1}=x_k-\alpha\nabla f(x_k)$. <br> <i>Hint:</i> just substitute $d=-\nabla f(x_k)$ into the step $x_{k+1}=x_k+\alpha d$. :: A: take a step of size $\alpha>0$ along $d$: $x_{k+1}=x_k+\alpha d$ <br> substitute $d=-\nabla f(x_k)$ <br> $\boxed{x_{k+1}=x_k-\alpha\,\nabla f(x_k)}$
- Q: Derive Newton's update from the second-order Taylor model of $f$ at $x_k$ (step 1 of 2): write the model and set its gradient to zero. <br> <b>Formula:</b> $m(p)=f(x_k)+\nabla f(x_k)^\top p+\tfrac12 p^\top H(x_k)p$; target $H(x_k)p=-\nabla f(x_k)$. <br> <i>Hint:</i> take $\nabla_p m(p)=\nabla f(x_k)+H(x_k)p$ and set it to $0$. :: A: $m(p)=f(x_k)+\nabla f(x_k)^\top p+\tfrac12 p^\top H(x_k)\,p$ <br> $\nabla_p m(p)=\nabla f(x_k)+H(x_k)\,p$ <br> set $=0$: $\boxed{H(x_k)\,p=-\nabla f(x_k)}$ (the Newton system)
- Q: From the Newton system finish the Newton update (step 2 of 2). <br> <b>Formula:</b> $H(x_k)p=-\nabla f(x_k)$ and $x_{k+1}=x_k+p$; target $x_{k+1}=x_k-H(x_k)^{-1}\nabla f(x_k)$. <br> <i>Hint:</i> solve the system for $p$ by left-multiplying by $H^{-1}$, then add $p$ to $x_k$. :: A: solve for the step $p=-H(x_k)^{-1}\nabla f(x_k)$ <br> apply it: $x_{k+1}=x_k+p$ <br> $\boxed{x_{k+1}=x_k-H(x_k)^{-1}\nabla f(x_k)}$
- Q: Derive the Lagrange tangency condition for $\min f(x)$ s.t. $g(x)=0$ (step 1 of 2). <br> <b>Formula:</b> feasible steps satisfy $\nabla g^\top\delta=0$ and at the optimum $\nabla f^\top\delta=0$; target $\nabla f\parallel\nabla g$. <br> <i>Hint:</i> a feasible move stays on $g=0$ (so it is tangent, $\nabla g\perp\delta$); at the optimum no feasible move lowers $f$, so $\nabla f\perp\delta$ too — both are orthogonal to the same tangent. :: A: A feasible step $\delta$ must keep $g=0$, so $\nabla g^\top\delta=0$ (tangent to the curve) <br> at the optimum no feasible step lowers $f$, so $\nabla f^\top\delta=0$ for every such $\delta$ <br> both $\nabla f$ and $\nabla g$ are orthogonal to the same tangent direction <br> $\boxed{\nabla f \parallel \nabla g}$ (collinear gradients at $x^\star$)
- Q: From the collinearity $\nabla f\parallel\nabla g$ derive the Lagrangian stationarity condition (step 2 of 2). <br> <b>Formula:</b> $\mathcal{L}(x,\lambda)=f(x)+\lambda^\top g(x)$; target $\nabla_x\mathcal{L}=0$. <br> <i>Hint:</i> collinear means $\nabla f=-\lambda\nabla g$ for some scalar $\lambda$; move everything to one side. :: A: collinear means $\nabla f=-\lambda\,\nabla g$ for some scalar $\lambda$ <br> rearrange: $\nabla f+\lambda\,\nabla g=0$ <br> this is exactly $\nabla_x\mathcal{L}=0$ for $\mathcal{L}(x,\lambda)=f(x)+\lambda^\top g(x)$ <br> $\boxed{\nabla_x\big(f+\lambda^\top g\big)=0}$
- Q: Solve $\min x^2+y^2$ s.t. $x+y=1$ in closed form via the Lagrangian. <br> <b>Formula:</b> $\mathcal{L}=x^2+y^2+\lambda(x+y-1)$; stationarity $\partial_x\mathcal{L}=\partial_y\mathcal{L}=0$ plus feasibility $x+y=1$. <br> <i>Hint:</i> set $\partial_x\mathcal{L}=2x+\lambda$ and $\partial_y\mathcal{L}=2y+\lambda$ to zero first; they force $x=y$, then use $x+y=1$. :: A: $\mathcal{L}=x^2+y^2+\lambda(x+y-1)$ <br> $\partial_x\mathcal{L}=2x+\lambda=0,\ \partial_y\mathcal{L}=2y+\lambda=0\Rightarrow x=y$ <br> feasibility $x+y=1\Rightarrow 2x=1$ <br> $\boxed{(x,y)=(\tfrac12,\tfrac12)}$
- Q: Show why one Newton step solves a quadratic $f(x)=\tfrac12 x^\top A x - b^\top x$ exactly, and state the step. <br> <b>Formula:</b> $\nabla f=Ax-b$, $H=A$, Newton step $x_1=x_0-H^{-1}\nabla f(x_0)$; target $x_1=A^{-1}b$. <br> <i>Hint:</i> substitute $\nabla f(x_0)=Ax_0-b$ and $H=A$ into the Newton step; the $x_0$ terms cancel. :: A: $\nabla f=Ax-b$, $H=A$ (constant), so the Taylor model equals $f$ <br> Newton step from any $x_0$: $x_1=x_0-A^{-1}(Ax_0-b)$ <br> $=x_0-x_0+A^{-1}b$ <br> $\boxed{x_1=A^{-1}b=x^\star}$ (exact in one step)
- Q: Worked example (1D, small scale): minimize $f(x)=(x-3)^2$ by gradient descent from $x_0=0$ with $\alpha=0.1$; compute $x_1$ and $x_2$. <br> <b>Formula:</b> $f'(x)=2(x-3)$ and $x_{k+1}=x_k-\alpha f'(x_k)$. <br> <i>Hint:</i> compute $f'(0)$ first, take one step to get $x_1$, then plug $x_1$ back into $f'$ for $x_2$. :: A: $f'(x)=2(x-3)$ <br> $x_1=0-0.1\cdot2(0-3)=0+0.6=0.6$ <br> $x_2=0.6-0.1\cdot2(0.6-3)=0.6+0.48=1.08$ <br> $\boxed{x_1=0.6,\ x_2=1.08\to 3}$
- Q: Worked example (same problem, curvature): minimize $f(x)=(x-3)^2$ with one Newton step from $x_0=0$. <br> <b>Formula:</b> $f'(x)=2(x-3)$, $f''=2$, and $x_1=x_0-f'(x_0)/f''$. <br> <i>Hint:</i> evaluate $f'(0)$ and divide by $f''=2$; for a quadratic this lands exactly on the minimum. :: A: $f'(x)=2(x-3),\ f''=2$ <br> $x_1=0-\tfrac{1}{2}\cdot2(0-3)=0+3$ <br> $\boxed{x_1=3}$ (quadratic solved in one step)
- Q: Worked example (large step, divergence regime): for $f(x)=(x-3)^2$ a fixed step $\alpha$ gives the error recursion below. For which $\alpha$ does it diverge, and check $\alpha=1.1$. <br> <b>Formula:</b> $x_{k+1}-3=(1-2\alpha)(x_k-3)$; it converges iff $|1-2\alpha|<1$. <br> <i>Hint:</i> solve $|1-2\alpha|<1$ for the convergent range, then evaluate the factor $1-2\alpha$ at $\alpha=1.1$ and compare its magnitude to $1$. :: A: $|1-2\alpha|<1\Rightarrow 0<\alpha<1$ converges; $\alpha>1$ diverges <br> at $\alpha=1.1$: factor $1-2.2=-1.2$, $|{-1.2}|>1$ <br> $\boxed{\alpha=1.1\Rightarrow|x_k-3|\ \text{grows},\ \text{diverges}}$
- Q: Worked example (2D, robotics scale, units in metres): minimize $f(x,y)=x^2+4y^2$ from $(1,1)$ with $\alpha=0.1$; compute one GD step. <br> <b>Formula:</b> $\nabla f=(2x,\,8y)$ and $(x,y)\leftarrow(x,y)-\alpha\nabla f$. <br> <i>Hint:</i> evaluate the gradient $(2x,8y)$ at $(1,1)$ first, then subtract $\alpha$ times each component. :: A: $\nabla f=(2x,8y)=(2,8)$ at $(1,1)$ <br> $(x,y)\leftarrow(1,1)-0.1(2,8)$ <br> $=(1-0.2,\,1-0.8)$ <br> $\boxed{(0.8,\,0.2)}$ (the $y$-axis, with curvature 4, moves faster)
- Q: Test yourself + verify (limiting case): plug the candidate $(\tfrac12,\tfrac12)$ back into $\min x^2+y^2$ s.t. $x+y=1$ and confirm it is optimal. <br> <b>Formula:</b> check feasibility $x+y=1$ and collinearity $\nabla f=-\lambda\nabla g$ with $\nabla f=(2x,2y)$, $\nabla g=(1,1)$; cross-check with point-to-line distance $\tfrac{|x+y-1|}{\sqrt2}$ evaluated at the origin. <br> <i>Hint:</i> a pass means feasibility holds AND the two gradients are parallel; the distance $\tfrac1{\sqrt2}$ should match the known closest-point distance. :: A: feasibility: $\tfrac12+\tfrac12=1$ holds <br> cost $=\tfrac14+\tfrac14=\tfrac12$; the multiplier $\lambda=-2x=-1$ so $\nabla f=(1,1)=-\lambda\nabla g=(1,1)$, gradients collinear <br> distance to origin $=\tfrac{1}{\sqrt2}$ equals the known point-to-line distance $\tfrac{|0+0-1|}{\sqrt2}$ <br> $\boxed{(\tfrac12,\tfrac12)\ \text{is the global optimum}}$
- Q: Test yourself + verify (units / limit check): for $x_{k+1}=x_k-\alpha\nabla f$ with $x$ in m and $f$ in m$^2$, check the units of $\alpha\nabla f$ and the $\alpha\to0$ limit. <br> <b>Formula:</b> $[\nabla f]=[f]/[x]=\text{m}^2/\text{m}=\text{m}$, and $\alpha\nabla f$ must have units of m. <br> <i>Hint:</i> a pass means $\alpha$ comes out dimensionless and $\alpha\to0$ leaves $x_{k+1}=x_k$ (no move). :: A: $\nabla f$ has units m$^2$/m $=$ m, so $\alpha$ must be dimensionless for $\alpha\nabla f$ to be in m <br> as $\alpha\to0$, $x_{k+1}\to x_k$ (no move, no divergence) <br> $\boxed{[\alpha]=1\ \text{(dimensionless)},\ \alpha\to0\Rightarrow\text{stationary}}$ — confirms the update is consistent
- Q: Define a convex function two equivalent ways (segment form and Hessian form). <br> <b>Formula:</b> $f(\lambda a+(1-\lambda)b)\le\lambda f(a)+(1-\lambda)f(b)$; equivalently $H=\nabla^2 f\succeq 0$ everywhere. <br> <i>Hint:</i> the first says the chord lies above the graph; the second says curvature is non-negative. :: A: function form: $f(\lambda a+(1-\lambda)b)\le\lambda f(a)+(1-\lambda)f(b)$ <br> curvature form: $\boxed{H=\nabla^2 f\succeq0\ \text{everywhere}}$
- CLOZE: For a {{c1::convex}} function every local minimum is also a {{c2::global}} minimum, so gradient descent reaches the true optimum.
- CLOZE: At a constrained optimum the cost gradient is {{c1::parallel}} to the constraint gradient: $\nabla f=-\lambda^\top\nabla g$.
- CLOZE: Gauss-Newton approximates the Hessian as {{c1::$J^\top J$}}, and Levenberg-Marquardt {{c2::damps}} it. <br> <i>hint:</i> the approximation drops the second-derivative term of the residuals.
- Q: Name the four KKT conditions for $\min f(x)$ s.t. $g(x)=0,\ h(x)\le0$. <br> <b>Formula:</b> stationarity $\nabla f+\lambda^\top\nabla g+\mu^\top\nabla h=0$; primal feasibility; dual feasibility $\mu\ge0$; complementary slackness $\mu_i h_i(x)=0$. <br> <i>Hint:</i> two are feasibility conditions, one signs the inequality multipliers, one ties each $\mu_i$ to its constraint being active. :: A: Stationarity $\nabla f+\lambda^\top\nabla g+\mu^\top\nabla h=0$; primal feasibility; dual feasibility $\mu\ge0$; complementary slackness $\mu_i h_i=0$.
- Q: When-to-use: pick the solver by problem structure. <br> <b>Formula:</b> match structure $\to$ method (convex+linear, smooth non-convex, least-squares residuals, discrete/combinatorial). <br> <i>Hint:</i> closed form when the optimum is unique and linear; gradient/Newton when smooth; Gauss-Newton/LM when residuals are squared; search+refine when discrete. :: A: convex+linear $\to$ closed form; smooth non-convex $\to$ gradient/Newton from a good start; least-squares residuals $\to$ Gauss-Newton/LM; discrete/combinatorial $\to$ search + continuous refine.
- Q: Pitfall: why do ICP and IK fail from a bad initial guess, and what is the fix? <br> <b>Formula:</b> gradient methods satisfy only the local condition $\nabla f(x^\star)=0$. <br> <i>Hint:</i> think about what "local minimum" means on a non-convex landscape. :: A: Their objectives are non-convex; gradient methods only find a local minimum near the start, so warm-start them.
- Q: Pitfall: name the symptom and fix for an ill-conditioned (zig-zagging) descent. <br> <b>Formula:</b> conditioning is governed by the Hessian condition number $\kappa(H)$; large $\kappa$ means a long narrow valley. <br> <i>Hint:</i> the cure either averages past steps or rescales by curvature. :: A: A long narrow valley (large Hessian condition number); use momentum, preconditioning, or Newton's method.
- Q: Pitfall: what is the right stopping criterion for an iterative optimizer? <br> <b>Formula:</b> stop when $\|\nabla f(x_k)\|<\text{tol}$ and/or $\|x_{k+1}-x_k\|<\text{tol}$. <br> <i>Hint:</i> tie the test to the optimality condition $\nabla f=0$, not to a fixed iteration budget. :: A: Small gradient norm and/or small step size, not iteration count alone.
- Q: Code idiom: which Eigen solver applies the Newton step $p=-H^{-1}g$ for an SPD Hessian? <br> <b>Formula:</b> solve $Hp=-g$ via an $LDL^\top$ factorization. <br> <i>Hint:</i> it is the Cholesky-style decomposition Eigen exposes for symmetric positive-definite matrices. :: A: `H.ldlt().solve(g)`.
