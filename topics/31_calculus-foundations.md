# Calculus (foundations): derivatives, integrals, and multivariable calculus

Placement: [[01_vectors]] => **Calculus foundations** => [[33_optimization-foundations]], [[07_calculus-jacobians]], [[34_differential-equations]], [[19_laplacian-poisson]], [[32_probability-statistics]]

## Why it matters for the stone scan->CAD->IFC / CADFit / robotic-fabrication workflow

Every step that *fits*, *smooths*, or *moves* runs on calculus. Surface fitting and registration minimize an error by following its derivative downhill, so you must be able to differentiate the error. Laplacian smoothing and the Signed Heat method are built from the differential operators gradient, divergence, and Laplacian acting on a mesh. The volume and mass of a carved block are integrals over its solid. A robot arm integrates joint velocity into joint position. Calculus is the shared language of all of these, so this file is the floor that [[33_optimization-foundations]], [[34_differential-equations]], and [[19_laplacian-poisson]] stand on.

## Intuition

A derivative is how fast something changes *right now*: the slope of the curve at a point. An integral is the *total accumulation*: the area swept under the curve. In many variables, the gradient is the "uphill" arrow on a landscape, pointing in the steepest-ascent direction with a length equal to how steep it is. To minimize error you walk opposite the gradient.

Picture a marble resting on a curved metal sheet. The gradient at the marble's position points uphill; gravity pulls it the other way, toward the lowest point. Optimization is just rolling that marble. Geometry processing is asking how the sheet bends (its second derivatives).

## The math (first principles)

**Derivative.** For a scalar function of one variable,

$$ f'(x) = \frac{df}{dx} = \lim_{h \to 0} \frac{f(x+h) - f(x)}{h}. $$

It is the slope of the tangent line and has units of (output units) per (input unit).

**Rules.** Power $\frac{d}{dx}x^n = n x^{n-1}$; product $(fg)' = f'g + fg'$; quotient $(f/g)' = (f'g - fg')/g^2$; and the **chain rule**

$$ \frac{d}{dx} f(g(x)) = f'(g(x))\, g'(x). $$

The chain rule is the engine behind backpropagation and behind the Jacobians in [[07_calculus-jacobians]].

**Integral and the Fundamental Theorem.** The definite integral is the limit of Riemann sums, $\int_a^b f(x)\,dx = \lim_{n\to\infty}\sum_i f(x_i)\,\Delta x$. The Fundamental Theorem of Calculus links the two operations:

$$ \int_a^b f'(x)\,dx = f(b) - f(a), \qquad \frac{d}{dx}\int_a^x f(t)\,dt = f(x). $$

Differentiation and integration are inverses.

**Partial derivatives and the gradient.** For $f:\mathbb{R}^n \to \mathbb{R}$, the partial derivative $\partial f/\partial x_i$ differentiates with respect to one coordinate and holds the rest fixed. The **gradient** stacks them:

$$ \nabla f = \left[\frac{\partial f}{\partial x_1}, \dots, \frac{\partial f}{\partial x_n}\right]^\top. $$

The **directional derivative** along a unit vector $u$ is $D_u f = \nabla f \cdot u$. Because $\nabla f \cdot u = \lVert\nabla f\rVert \cos\theta$, the gradient is the direction of steepest ascent and its norm is the rate of climb.

**Vector calculus operators.** For a vector field $F = (F_1,F_2,F_3)$ and scalar field $f$:

$$ \nabla \cdot F = \sum_i \frac{\partial F_i}{\partial x_i}\ \text{(divergence)},\qquad \nabla \times F\ \text{(curl)},\qquad \Delta f = \nabla\cdot\nabla f = \sum_i \frac{\partial^2 f}{\partial x_i^2}\ \text{(Laplacian)}. $$

The Laplacian measures how much a value differs from the average of its neighbours; its discrete cousin is the cotangent Laplacian in [[19_laplacian-poisson]].

**Taylor expansion.** Near a point $x$, with step $\delta$,

$$ f(x+\delta) \approx f(x) + \nabla f(x)^\top \delta + \tfrac{1}{2}\,\delta^\top H(x)\,\delta, $$

where $H = \nabla^2 f$ is the Hessian. The first-order term drives gradient descent; the second-order term drives Newton's method. Every iterative solver in this corpus is a Taylor expansion truncated somewhere.

## Worked example

Let $f(x) = x^2$. Then $f'(x) = 2x$, so the slope at $x=3$ is $6$. A central finite difference with $h=10^{-3}$ gives $\frac{(3.001)^2 - (2.999)^2}{0.002} = \frac{9.006001 - 8.994001}{0.002} = \frac{0.012}{0.002} = 6.000$. The definite integral $\int_0^2 x^2\,dx = \big[\tfrac{x^3}{3}\big]_0^2 = \tfrac{8}{3} \approx 2.6667$.

Now a two-variable case: $f(x,y) = x^2 + xy$. The partials are $\partial f/\partial x = 2x + y$ and $\partial f/\partial y = x$, so $\nabla f = [2x+y,\ x]^\top$. At $(1,2)$, $\nabla f = [4,\ 1]^\top$. The directional derivative along $u=(1,0)$ is $\nabla f \cdot u = 4$: moving in $+x$ raises $f$ at rate 4. Steepest ascent is along $[4,1]/\sqrt{17}$.

## Code

Idiomatic C++ with Eigen: a central-difference gradient for any scalar objective, plus Simpson's rule for a 1-D integral. The finite-difference gradient is also how you *gradient-check* an analytic derivative.

```cpp
#include <Eigen/Dense>
#include <functional>
#include <cmath>

using Vec = Eigen::VectorXd;

// Central-difference gradient of f at x.  grad_i ~ (f(x + h e_i) - f(x - h e_i)) / (2h)
// Central difference is O(h^2) accurate; pick h around sqrt(machine_eps)*|x_i|.
Vec numericalGradient(const std::function<double(const Vec&)>& f, const Vec& x) {
    Vec g(x.size());
    for (int i = 0; i < x.size(); ++i) {
        double xi = x(i);
        double h  = 1e-6 * std::max(1.0, std::abs(xi));   // scale step to the variable
        Vec xp = x; xp(i) = xi + h;
        Vec xm = x; xm(i) = xi - h;
        g(i) = (f(xp) - f(xm)) / (2.0 * h);               // partial df/dx_i
    }
    return g;
}

// Simpson's rule: integral of g on [a,b] with n (even) subintervals. Error O(h^4).
double integrateSimpson(const std::function<double(double)>& g, double a, double b, int n) {
    if (n % 2 != 0) ++n;                  // Simpson needs an even count
    const double h = (b - a) / n;
    double s = g(a) + g(b);
    for (int k = 1; k < n; ++k)
        s += (k % 2 ? 4.0 : 2.0) * g(a + k * h);
    return s * h / 3.0;
}

// Example: f(x,y) = x^2 + x*y ; check gradient at (1,2) equals [4,1].
//   auto f = [](const Vec& v){ return v(0)*v(0) + v(0)*v(1); };
//   Vec x(2); x << 1, 2;  numericalGradient(f, x)  ->  ~[4, 1]
//   integrateSimpson([](double t){ return t*t; }, 0, 2, 100)  ->  ~2.6667
```

## Connections

- [[01_vectors]] supplies the dot product used in the directional derivative and gradient.
- [[07_calculus-jacobians]] generalizes the derivative to vector-valued functions (the Jacobian) and is the applied, solver-facing layer above this foundation.
- [[33_optimization-foundations]] uses the gradient (first-order) and Hessian (second-order) directly; "set $\nabla f = 0$" is the optimality condition.
- [[34_differential-equations]] is equations *about* derivatives; the Laplacian here is the operator in the heat and Poisson equations.
- [[19_laplacian-poisson]] is the discrete Laplacian on meshes, the workhorse of geometry processing.
- [[32_probability-statistics]] needs integration to define continuous probability and expectation.

## Pitfalls and failure modes

- **Step size in finite differences.** Forward difference is only $O(h)$ accurate; central difference is $O(h^2)$. Too large $h$ gives truncation error, too small $h$ gives catastrophic cancellation. A good default is $h \approx \sqrt{\varepsilon}\,|x|$ for central differences.
- **Non-differentiable points.** $\lvert x\rvert$, ReLU, and rigid contact have kinks; the derivative does not exist there. Use subgradients or smooth approximations.
- **Gradient vs Jacobian.** Gradient is for scalar output; a vector-valued map has a Jacobian matrix. Mixing them up is a common bug.
- **Units and scaling.** A derivative carries units. If $x$ is in metres and $f$ is energy in joules, $\nabla f$ is in joules per metre. Badly scaled variables wreck conditioning downstream.
- **Always gradient-check.** Before trusting a hand-derived gradient, compare it to `numericalGradient`. Disagreement beyond a few $\times h^2$ means a bug in the analytic derivative.

## Practice

1. Differentiate $f(x) = \sin(x^2)$ by hand with the chain rule, then confirm with the central-difference code at $x = 1$.
2. Integrate $\int_0^\pi \sin x\,dx$ analytically (answer 2) and with `integrateSimpson` for $n = 10, 100, 1000$; tabulate the error and confirm the $O(h^4)$ rate.
3. For $f(x,y) = x^2 + 3xy + y^2$, compute $\nabla f$ and the Hessian by hand, then gradient-check $\nabla f$ at three random points.
4. Compute the steepest-ascent unit direction of $f(x,y)=e^{-(x^2+y^2)}$ at $(0.5, 0.5)$ and the directional derivative along $(1,1)/\sqrt2$.
5. Write a 2-D trapezoid integrator and estimate the area of a unit disk by integrating the indicator $x^2+y^2 \le 1$; compare to $\pi$.

## References

- James Stewart, *Calculus: Early Transcendentals* (single and multivariable chapters).
- Gilbert Strang, *Calculus*, MIT OpenCourseWare (free): https://ocw.mit.edu/courses/res-18-001-calculus-online-textbook-spring-2005/
- Keenan Crane, *Discrete Differential Geometry* course notes (vector calculus: grad, div, curl, Laplacian): https://www.cs.cmu.edu/~kmcrane/Projects/DDG/
- "The Essential Mathematics for Computational Design," Robert McNeel & Associates (bridge text for designers).

## Flashcard seeds

- Q: From first principles, derive $f'(x)$ for $f(x)=x^2$. <br> <b>Formula:</b> $f'(x)=\lim_{h\to0}\frac{f(x+h)-f(x)}{h}$, target $\boxed{2x}$ <br> <i>Hint:</i> substitute $f(x+h)=(x+h)^2$, expand the square, then cancel the $x^2$ terms before dividing by $h$. :: A: $f'(x)=\lim_{h\to0}\frac{(x+h)^2-x^2}{h}$ <br> $=\lim_{h\to0}\frac{x^2+2xh+h^2-x^2}{h}$ <br> $=\lim_{h\to0}\frac{2xh+h^2}{h}=\lim_{h\to0}(2x+h)$ <br> $\boxed{f'(x)=2x}$.
- Q: Derive the central finite-difference identity from two Taylor expansions (step 1 of 2). <br> <b>Formula:</b> $f(x\pm h)=f(x)\pm f'(x)h+\tfrac12 f''(x)h^2\pm\tfrac16 f'''(x)h^3+\dots$ <br> <i>Hint:</i> write both expansions, then SUBTRACT $f(x-h)$ from $f(x+h)$ so the even-power terms cancel. :: A: $f(x+h)=f(x)+f'(x)h+\tfrac12 f''(x)h^2+\tfrac16 f'''(x)h^3+\dots$ <br> $f(x-h)=f(x)-f'(x)h+\tfrac12 f''(x)h^2-\tfrac16 f'''(x)h^3+\dots$ <br> subtract: $f(x+h)-f(x-h)=2f'(x)h+\tfrac13 f'''(x)h^3+\dots$ <br> $\boxed{f(x+h)-f(x-h)=2f'(x)h+O(h^3)}$.
- Q: From $f(x+h)-f(x-h)=2f'(x)h+O(h^3)$, derive the central-difference estimator and its error order (step 2 of 2). <br> <b>Formula:</b> start from $f(x+h)-f(x-h)=2f'(x)h+O(h^3)$, target $\boxed{f'(x)\approx\frac{f(x+h)-f(x-h)}{2h}}$ <br> <i>Hint:</i> divide both sides by $2h$; note $O(h^3)/h=O(h^2)$. :: A: divide by $2h$: $\frac{f(x+h)-f(x-h)}{2h}=f'(x)+O(h^2)$ <br> the leading error term is $\tfrac16 f'''(x)h^2$ <br> $\boxed{f'(x)\approx\frac{f(x+h)-f(x-h)}{2h},\ \text{error }O(h^2)}$.
- Q: Derive why a FORWARD difference is only $O(h)$ accurate. <br> <b>Formula:</b> $f(x+h)=f(x)+f'(x)h+\tfrac12 f''(x)h^2+\dots$, target $\boxed{f'(x)+O(h)}$ <br> <i>Hint:</i> rearrange the single expansion to isolate $\frac{f(x+h)-f(x)}{h}$; the $\tfrac12 f''(x)h$ term has no symmetric partner to cancel it. :: A: $f(x+h)=f(x)+f'(x)h+\tfrac12 f''(x)h^2+\dots$ <br> $\frac{f(x+h)-f(x)}{h}=f'(x)+\tfrac12 f''(x)h+\dots$ <br> leading error $\propto h$ does not cancel (no symmetric term) <br> $\boxed{\frac{f(x+h)-f(x)}{h}=f'(x)+O(h)}$.
- Q: Differentiate $\sin(x^2)$ using the chain rule. <br> <b>Formula:</b> $\frac{d}{dx}f(g(x))=f'(g(x))\,g'(x)$, with $\frac{d}{dx}\sin u=\cos u$ and $\frac{d}{dx}x^2=2x$ <br> <i>Hint:</i> take the inner function $g(x)=x^2$ first, differentiate the outer $\sin$ at $g(x)$, then multiply by $g'(x)$. :: A: let $g(x)=x^2$, outer $f(u)=\sin u$ <br> $\frac{d}{dx}f(g(x))=f'(g(x))\,g'(x)$ <br> $f'(u)=\cos u$, $g'(x)=2x$ <br> $\boxed{\frac{d}{dx}\sin(x^2)=2x\cos(x^2)}$.
- Q: Show the directional derivative $D_u f$ is largest when $u$ points along $\nabla f$. <br> <b>Formula:</b> $D_u f=\nabla f\cdot u=\lVert\nabla f\rVert\,\lVert u\rVert\cos\theta$, $u$ a unit vector <br> <i>Hint:</i> set $\lVert u\rVert=1$, then ask which value of $\cos\theta$ maximizes the expression. :: A: $D_u f=\nabla f\cdot u=\lVert\nabla f\rVert\,\lVert u\rVert\cos\theta$ <br> with unit $u$, $\lVert u\rVert=1$, so $D_u f=\lVert\nabla f\rVert\cos\theta$ <br> maximized when $\cos\theta=1$, i.e. $u\parallel\nabla f$ <br> $\boxed{\max_u D_u f=\lVert\nabla f\rVert\ \text{at}\ u=\nabla f/\lVert\nabla f\rVert}$.
- Q: Find the gradient of $f(x,y)=x^2+xy$ at $(1,2)$ (marble-on-a-sheet, step 1 of 2). <br> <b>Formula:</b> $\nabla f=\left[\frac{\partial f}{\partial x},\frac{\partial f}{\partial y}\right]^\top$ <br> <i>Hint:</i> compute $\partial f/\partial x$ (treat $y$ as constant) and $\partial f/\partial y$ (treat $x$ as constant), then plug in $x=1,\ y=2$. :: A: $\partial f/\partial x=2x+y$, $\partial f/\partial y=x$ <br> $\nabla f=[2x+y,\ x]^\top$ <br> at $(1,2)$: $\nabla f=[2(1)+2,\ 1]^\top$ <br> $\boxed{\nabla f(1,2)=[4,\ 1]^\top}$.
- Q: From $\nabla f(1,2)=[4,1]^\top$, find the unit steepest-DESCENT direction the marble rolls (step 2 of 2). <br> <b>Formula:</b> descent direction $=-\nabla f$; unit vector $=\dfrac{-\nabla f}{\lVert\nabla f\rVert}$, $\lVert v\rVert=\sqrt{v_1^2+v_2^2}$ <br> <i>Hint:</i> negate the gradient first, then divide by its length $\sqrt{4^2+1^2}$. :: A: descent is $-\nabla f=[-4,-1]^\top$ <br> $\lVert\nabla f\rVert=\sqrt{4^2+1^2}=\sqrt{17}$ <br> normalize: $\boxed{\hat d=\frac{1}{\sqrt{17}}[-4,\ -1]^\top}$.
- Q: Derive the Newton step by minimizing the second-order Taylor model. <br> <b>Formula:</b> $m(\delta)=f(x)+\nabla f(x)^\top\delta+\tfrac12\delta^\top H\delta$, target $\boxed{\delta=-H^{-1}\nabla f(x)}$ <br> <i>Hint:</i> differentiate $m$ with respect to $\delta$, set the result to zero, then solve the linear system for $\delta$. :: A: $m(\delta)=f(x)+\nabla f(x)^\top\delta+\tfrac12\delta^\top H\delta$ <br> $\nabla_\delta m=\nabla f(x)+H\delta$ <br> set $=0$: $H\delta=-\nabla f(x)$ <br> $\boxed{\delta=-H^{-1}\nabla f(x)}$ (the Newton step).
- Q: Evaluate $\int_0^2 x^2\,dx$ using the Fundamental Theorem. <br> <b>Formula:</b> $\int_a^b f'(x)\,dx=F(b)-F(a)$ with antiderivative $F(x)=x^3/3$ (since $F'(x)=x^2$) <br> <i>Hint:</i> find the antiderivative first, then evaluate at $2$ minus at $0$. :: A: antiderivative of $x^2$ is $F(x)=x^3/3$ since $F'(x)=x^2$ <br> $\int_0^2 x^2\,dx=F(2)-F(0)$ <br> $=8/3-0$ <br> $\boxed{\int_0^2 x^2\,dx=8/3\approx2.667}$.
- Q: Worked example (small scale): estimate $f'(3)$ for $f(x)=x^2$ with a central difference, $h=10^{-3}$. <br> <b>Formula:</b> $f'(x)\approx\frac{f(x+h)-f(x-h)}{2h}$ <br> <i>Hint:</i> compute $(3.001)^2$ and $(2.999)^2$ first, then divide their difference by $2h=0.002$. :: A: $\frac{(3.001)^2-(2.999)^2}{0.002}=\frac{9.006001-8.994001}{0.002}$ <br> $=\frac{0.012}{0.002}=6.000$ <br> matches exact $f'(3)=2\cdot3=6$ <br> $\boxed{f'(3)\approx6.000}$.
- Q: Worked example (large scale, m units): a quarry block has cross-section $A(z)=4-z^2$ (m$^2$); find its volume over $z\in[0,2]$ m. <br> <b>Formula:</b> $V=\int_a^b A(z)\,dz$, antiderivative of $4-z^2$ is $4z-z^3/3$ <br> <i>Hint:</i> integrate term by term, then evaluate $[4z-z^3/3]$ at $z=2$ minus $z=0$. :: A: $V=\int_0^2(4-z^2)\,dz=[4z-z^3/3]_0^2$ <br> $=(8-8/3)-(0)$ <br> $=24/3-8/3=16/3$ <br> $\boxed{V=16/3\approx5.333\ \text{m}^3}$.
- Q: Worked example (mm shop scale): a fillet profile $y(x)=x^2$ (mm); find the area under it over $x\in[0,5]$ mm. <br> <b>Formula:</b> $A=\int_a^b y(x)\,dx$, antiderivative of $x^2$ is $x^3/3$ <br> <i>Hint:</i> evaluate $[x^3/3]$ at $x=5$ minus at $x=0$. :: A: $A=\int_0^5 x^2\,dx=[x^3/3]_0^5$ <br> $=125/3-0$ <br> $\boxed{A=125/3\approx41.67\ \text{mm}^2}$.
- Q: Worked example (3D gradient, registration error): for $E(x,y,z)=x^2+y^2+z^2$, find $\nabla E$ and the unit steepest-descent direction at $(1,-2,2)$. <br> <b>Formula:</b> $\nabla E=[2x,2y,2z]^\top$; unit descent $=\dfrac{-\nabla E}{\lVert\nabla E\rVert}$ <br> <i>Hint:</i> plug the point into $\nabla E$ first, then compute $\lVert\nabla E\rVert=\sqrt{2^2+4^2+4^2}$ before normalizing. :: A: $\nabla E=[2x,2y,2z]^\top$ <br> at $(1,-2,2)$: $\nabla E=[2,-4,4]^\top$ <br> descent $-\nabla E=[-2,4,-4]^\top$, $\lVert\nabla E\rVert=\sqrt{4+16+16}=6$ <br> $\boxed{\hat d=\tfrac16[-2,4,-4]^\top=[-\tfrac13,\tfrac23,-\tfrac23]^\top}$.
- Q: Closing test: gradient-check the analytic $\nabla f=[2x+y,x]^\top$ for $f=x^2+xy$ at $(1,2)$ against central differences ($h=10^{-3}$). <br> <b>Formula:</b> $\frac{\partial f}{\partial x_i}\approx\frac{f(x+he_i)-f(x-he_i)}{2h}$; pass if the two agree within $O(h^2)$ <br> <i>Hint:</i> compute the analytic vector $[4,1]$ first, then perturb $x$ and $y$ one at a time; a pass means the FD values round to $4$ and $1$. :: A: analytic $\nabla f=[4,1]^\top$ <br> FD in $x$: $\frac{f(1.001,2)-f(0.999,2)}{0.002}=\frac{(1.001^2+1.001\cdot2)-(0.999^2+0.999\cdot2)}{0.002}\approx4.000$ <br> FD in $y$: $\frac{f(1,2.001)-f(1,1.999)}{0.002}=\frac{(1+2.001)-(1+1.999)}{0.002}=1.000$ <br> $\boxed{[4.000,1.000]\approx[4,1]\ \checkmark}$ agreement within $O(h^2)$.
- Q: Closing test (units + limiting case): for energy $f$ in joules and position $x$ in metres, give the units of $\nabla f$ and check the FTC zero-width case $\int_a^a f'\,dx$. <br> <b>Formula:</b> $[\nabla f]=[f]/[x]$; FTC $\int_a^b f'\,dx=f(b)-f(a)$ <br> <i>Hint:</i> divide the unit of $f$ by the unit of $x$; for the integral, set both limits equal so $b=a$. A pass gives $0$. :: A: $\nabla f$ has units J/m (force-like, joules per metre) <br> FTC: $\int_a^a f'(x)\,dx=f(a)-f(a)=0$, consistent with zero-width interval <br> $\boxed{[\nabla f]=\text{J/m},\quad \int_a^a f'\,dx=0}$.
- Q: When do you use a gradient versus a Jacobian? <br> <i>Hint:</i> think about the SHAPE of the function's output: one number, or a vector? :: A: Gradient for a scalar-valued $f:\mathbb{R}^n\to\mathbb{R}$ (a vector); Jacobian (a matrix) for a vector-valued map $F:\mathbb{R}^n\to\mathbb{R}^m$.
- Q: Define the derivative $f'(x)$ as a limit (recall). :: A: $\boxed{f'(x)=\lim_{h\to0}\frac{f(x+h)-f(x)}{h}}$, the slope of the tangent line.
- Q: Pitfall: why can a finite-difference step $h$ be too SMALL? <br> <i>Hint:</i> think about what happens when you subtract two nearly equal floating-point numbers. :: A: Catastrophic cancellation: subtracting nearly equal floating-point numbers loses significant digits; default $h\approx\sqrt{\varepsilon}\,|x|$ balances truncation vs round-off.
- Q: When-to-use: which Taylor term drives gradient descent, which drives Newton's method? <br> <b>Formula:</b> $f(x+\delta)\approx f(x)+\nabla f(x)^\top\delta+\tfrac12\delta^\top H\delta$ <br> <i>Hint:</i> match the first-order term to one method, the second-order (Hessian) term to the other. :: A: The first-order (gradient) term drives gradient descent; the second-order (Hessian) term drives Newton's method.
- CLOZE: The {{c1::gradient}} points in the direction of steepest {{c2::ascent}}, and its norm equals the rate of climb in that direction.
- CLOZE: The Laplacian $\Delta f=\nabla\cdot\nabla f$ measures how much a value differs from the {{c1::average of its neighbours}}.
- CLOZE: A central finite difference is {{c1::$O(h^2)$}} accurate, while a forward difference is only {{c2::$O(h)$}} accurate. <br> <i>hint:</i> the symmetric scheme cancels one more error term.
- CLOZE: Divergence of a field $F$ is {{c1::$\nabla\cdot F=\sum_i \partial F_i/\partial x_i$}}; the Laplacian is the divergence of the {{c2::gradient}}.
- CLOZE: The Fundamental Theorem of Calculus states $\int_a^b f'(x)\,dx={{c1::f(b)-f(a)}}$, so differentiation and integration are {{c2::inverses}}.
