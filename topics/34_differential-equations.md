# Differential equations (ODEs, PDEs) and numerical integration

Placement: [[31_calculus-foundations]] => **Differential equations** => [[35_dynamics-control]], [[20_sdf-signed-heat]]

## Why it matters for the stone scan->CAD->IFC / CADFit / robotic-fabrication workflow

Two pillars of this corpus are differential equations in disguise. Geometry processing solves *partial* differential equations on meshes: Laplacian smoothing, the Poisson reconstruction, and the Signed Heat method are all diffusion or steady-state PDEs. Robot motion is governed by *ordinary* differential equations: the equations of motion in [[35_dynamics-control]] are integrated forward in time to simulate or control an arm carving stone. Knowing how to set up and numerically integrate these equations is what turns a static formula into a running simulation or controller.

## Intuition

A differential equation describes a rule for change rather than a value. "The cup cools at a rate proportional to how much hotter it is than the room" is a differential equation; its solution is the cooling curve over time. An ODE involves derivatives in one variable (usually time); a PDE involves partials in several (space and time). Solving means rolling the rule forward from a starting condition. Numerical integration is taking small time steps and updating by the local rule, like walking a trail by repeatedly looking at the slope under your feet.

## The math (first principles)

**ODE, initial value problem.** A first-order ODE is

$$ \dot y = \frac{dy}{dt} = f(t, y), \qquad y(t_0) = y_0. $$

Higher-order ODEs reduce to first-order systems by stacking derivatives into the state vector (e.g. position and velocity).

**Analytic vs numerical.** A few ODEs have closed forms (e.g. $\dot y = -ky$ gives $y=y_0 e^{-kt}$). Most do not, so we integrate numerically.

**Numerical integrators.**
- Forward (explicit) Euler: $y_{n+1} = y_n + h\,f(t_n, y_n)$. Simple, first-order accurate ($O(h)$), can be unstable for stiff problems.
- Backward (implicit) Euler: $y_{n+1} = y_n + h\,f(t_{n+1}, y_{n+1})$. Stable for stiff problems but requires solving an equation each step.
- Runge-Kutta 4 (RK4): the workhorse, fourth-order accurate ($O(h^4)$):

$$ \begin{aligned} k_1 &= f(t_n, y_n), & k_2 &= f(t_n + \tfrac{h}{2}, y_n + \tfrac{h}{2}k_1),\\ k_3 &= f(t_n + \tfrac{h}{2}, y_n + \tfrac{h}{2}k_2), & k_4 &= f(t_n + h, y_n + h k_3),\\ y_{n+1} &= y_n + \tfrac{h}{6}(k_1 + 2k_2 + 2k_3 + k_4). \end{aligned} $$

**Stability and stiffness.** A method is stable if errors do not blow up. Stiff systems mix very fast and very slow dynamics; explicit methods need tiny steps to stay stable, so implicit methods are preferred there.

**PDEs by example.** The **heat / diffusion equation** spreads a quantity over space and time:

$$ \frac{\partial u}{\partial t} = \alpha\,\Delta u, $$

where $\Delta$ is the Laplacian from [[31_calculus-foundations]]. Diffusing a field for a short time smooths it; this is exactly the first step of the Signed Heat method in [[20_sdf-signed-heat]] and of the heat method for geodesic distance. Its steady state ($\partial u/\partial t = 0$) is the **Poisson equation**

$$ \Delta u = f, $$

the central solve of [[19_laplacian-poisson]]. On a mesh, $\Delta$ becomes the cotangent Laplacian matrix and the PDE becomes a sparse linear system.

**Boundary and initial conditions.** A PDE needs conditions on the domain boundary (Dirichlet fixes values, Neumann fixes normal derivatives); an ODE needs an initial state. Without them the solution is not unique.

## Worked example

Cooling: $\dot y = -2y$, $y(0)=1$, true solution $y(t)=e^{-2t}$. At $t=1$, $y = e^{-2}\approx 0.1353$. One forward-Euler step with $h=0.5$: $y_1 = 1 + 0.5\cdot(-2\cdot 1) = 0$ then $y_2 = 0 + 0.5\cdot(-2\cdot 0)=0$; Euler with this big step is wildly off (it even crosses zero). With $h=0.1$, ten steps give $y_{10}=(1-0.2)^{10}=0.8^{10}\approx 0.107$, closer. RK4 with $h=0.5$ in two steps lands near $0.135$. This shows order of accuracy in action: RK4 with a coarse step beats Euler with a fine one.

## Code

Idiomatic C++ with Eigen: a generic RK4 integrator for a vector ODE $\dot y = f(t,y)$. The state can be a robot's joint positions and velocities (feeds [[35_dynamics-control]]).

```cpp
#include <Eigen/Dense>
#include <functional>

using Vec = Eigen::VectorXd;

// One RK4 step of dy/dt = f(t, y) with step size h.
Vec rk4Step(const std::function<Vec(double, const Vec&)>& f,
            double t, const Vec& y, double h) {
    Vec k1 = f(t,           y);
    Vec k2 = f(t + 0.5*h,   y + 0.5*h*k1);
    Vec k3 = f(t + 0.5*h,   y + 0.5*h*k2);
    Vec k4 = f(t + h,       y + h*k3);
    return y + (h/6.0) * (k1 + 2.0*k2 + 2.0*k3 + k4);   // 4th-order update
}

// Integrate from t0 to t1 in n steps.
Vec integrate(const std::function<Vec(double, const Vec&)>& f,
              double t0, Vec y, double t1, int n) {
    const double h = (t1 - t0) / n;
    double t = t0;
    for (int i = 0; i < n; ++i) { y = rk4Step(f, t, y, h); t += h; }
    return y;
}
// Example: spring-mass  y=[pos,vel],  f = [vel, -k*pos]  oscillates.
```

## Connections

- [[31_calculus-foundations]] defines the derivatives and the Laplacian operator these equations are written in.
- [[19_laplacian-poisson]] is the discrete Poisson solve, the steady state of the heat equation.
- [[20_sdf-signed-heat]] runs short-time heat diffusion as its first step.
- [[35_dynamics-control]] integrates the robot equations of motion with exactly these integrators.
- [[08_nonlinear-least-squares]] appears inside implicit (backward Euler) steps, which solve a nonlinear equation per step.

## Pitfalls and failure modes

- **Step size vs stability.** Explicit Euler diverges if $h$ is too large for the dynamics; reduce $h$ or switch to an implicit method.
- **Stiffness.** Mixed fast/slow dynamics force explicit methods to crawl; use implicit (backward Euler, BDF) integrators.
- **Order vs cost.** RK4 costs four function evaluations per step but its $O(h^4)$ accuracy usually beats many cheap Euler steps.
- **Missing or wrong boundary/initial conditions.** A PDE without proper boundary conditions has no unique solution; double-check Dirichlet vs Neumann.
- **Energy drift in long simulations.** Generic integrators leak or gain energy over time; for physical systems consider symplectic integrators.

## Practice

1. Integrate $\dot y=-2y,\ y(0)=1$ to $t=1$ with forward Euler ($h=0.1, 0.01$) and RK4 ($h=0.1$); tabulate error against $e^{-2}$.
2. Write the spring-mass system as a first-order ODE in state $[x,\dot x]$ and simulate one period; check the period against $2\pi\sqrt{m/k}$.
3. Make the spring stiff ($k$ large) and show explicit Euler needs a tiny step while backward Euler stays stable.
4. Discretize the 1-D heat equation on a line of nodes with the Laplacian stencil and watch an initial spike diffuse and flatten.
5. Solve the 1-D Poisson equation $u'' = f$ with Dirichlet ends as a tridiagonal linear system and compare to the analytic solution.

## References

- Lloyd N. Trefethen, *Finite Difference and Spectral Methods* / *Spectral Methods in MATLAB*.
- Ascher and Petzold, *Computer Methods for Ordinary Differential Equations and Differential-Algebraic Equations*.
- Crane, Weischedel, Wardetzky, "Geodesics in Heat" (the heat method): https://www.cs.cmu.edu/~kmcrane/Projects/HeatMethod/
- Strang, *Computational Science and Engineering* (PDEs, Laplacians, sparse solves).

## Flashcard seeds

- Q: Derive the forward Euler update from the slope rule $\dot y=f(t_n,y_n)$. <br> <b>Formula:</b> $\dfrac{y_{n+1}-y_n}{h}\approx f(t_n,y_n)$ (target: $y_{n+1}=y_n+h\,f(t_n,y_n)$) <br> <i>Hint:</i> replace the derivative with the forward finite difference, then multiply both sides by $h$. :: A: The rule is $\dot y=f(t_n,y_n)$ <br> approximate it over step $h$: $\dfrac{y_{n+1}-y_n}{h}\approx f(t_n,y_n)$ <br> multiply through by $h$ and solve for the next value <br> $\boxed{y_{n+1}=y_n+h\,f(t_n,y_n)}$.
- Q: Derive forward Euler as a truncated Taylor series (alternative route). <br> <b>Formula:</b> $y(t_n+h)=y(t_n)+h\,\dot y(t_n)+\tfrac{h^2}{2}\ddot y(\xi)$, with $\dot y=f(t,y)$ <br> <i>Hint:</i> substitute the rule for $\dot y(t_n)$, then drop the $\tfrac{h^2}{2}\ddot y$ remainder. :: A: Expand $y(t_n+h)=y(t_n)+h\,\dot y(t_n)+\tfrac{h^2}{2}\ddot y(\xi)$ <br> substitute the rule $\dot y(t_n)=f(t_n,y_n)$ <br> drop the $O(h^2)$ remainder term <br> $\boxed{y_{n+1}=y_n+h\,f(t_n,y_n)},\ \text{local error }O(h^2)$.
- Q: Show forward Euler's GLOBAL error is $O(h)$, i.e. it is first-order. <br> <b>Formula:</b> local error $=\tfrac{h^2}{2}\ddot y=O(h^2)$; steps to time $T$ are $N=T/h$ <br> <i>Hint:</i> multiply the per-step error by the number of steps $N=T/h$ and simplify. :: A: The local (one-step) error is the dropped term $\tfrac{h^2}{2}\ddot y=O(h^2)$ <br> reaching a fixed time $T$ needs $N=T/h$ steps <br> errors accumulate: global $\sim N\cdot O(h^2)=\tfrac{T}{h}O(h^2)$ <br> $\boxed{\text{global error}=O(h)\ \Rightarrow\ \text{first-order method}}$.
- Q: Derive the BACKWARD (implicit) Euler update by evaluating the slope at the END of the step. <br> <b>Formula:</b> $\dfrac{y_{n+1}-y_n}{h}\approx f(t_{n+1},y_{n+1})$ (target: $y_{n+1}=y_n+h\,f(t_{n+1},y_{n+1})$) <br> <i>Hint:</i> take the slope at $t_{n+1}=t_n+h$, then multiply by $h$; note $y_{n+1}$ appears on both sides. :: A: Use a backward finite difference: $\dfrac{y_{n+1}-y_n}{h}\approx f(t_{n+1},y_{n+1})$ <br> the slope is taken at the new point $t_{n+1}=t_n+h$ <br> rearrange for the unknown $y_{n+1}$ (it appears on both sides) <br> $\boxed{y_{n+1}=y_n+h\,f(t_{n+1},y_{n+1})}$ (solve an equation each step).
- Q: Derive the RK4 final combination from the four slope samples $k_1,k_2,k_3,k_4$ (step 2 of 2). <br> <b>Formula:</b> $y_{n+1}=y_n+\int_{t_n}^{t_n+h} f\,dt$; Simpson weights $(1,4,1)/6$ <br> <i>Hint:</i> apply Simpson's rule with the two midpoint samples $k_2,k_3$ sharing the weight-4 term, giving weights $(1,2,2,1)/6$. :: A: Approximate $y_{n+1}=y_n+\int_{t_n}^{t_n+h} f\,dt$ <br> Simpson's rule weights endpoints by 1 and midpoint by 4, and RK4 uses two midpoint samples $k_2,k_3$ <br> weights become $(1,2,2,1)/6$ over the four slopes <br> $\boxed{y_{n+1}=y_n+\tfrac{h}{6}(k_1+2k_2+2k_3+k_4)}$.
- Q: Derive the RK4 stage samples for $\dot y=f(t,y)$ with step $h$ (step 1 of 2). <br> <b>Formula:</b> $k_1=f(t_n,y_n)$; $k_{i+1}=f(t_n+ch,\,y_n+ch\,k_i)$ <br> <i>Hint:</i> compute $k_1$ at the start, then half-step along $k_1$ for $k_2$, along $k_2$ for $k_3$, and full-step along $k_3$ for $k_4$. :: A: Sample the slope at the start: $k_1=f(t_n,y_n)$ <br> step a half-$h$ along $k_1$ and resample: $k_2=f(t_n+\tfrac{h}{2},y_n+\tfrac{h}{2}k_1)$ <br> refine the midpoint slope: $k_3=f(t_n+\tfrac{h}{2},y_n+\tfrac{h}{2}k_2)$ <br> step a full $h$ along $k_3$: $\boxed{k_4=f(t_n+h,\,y_n+h\,k_3)}$.
- Q: Derive the steady-state PDE (Poisson equation) from the heat equation (step 2 of 2). <br> <b>Formula:</b> $\partial u/\partial t=\alpha\,\Delta u + f$ (target: $\Delta u=f$) <br> <i>Hint:</i> set $\partial u/\partial t=0$ for steady state, then divide through by $\alpha\neq0$ and rearrange. :: A: Steady state means the field stops changing in time: $\partial u/\partial t=0$ <br> so $0=\alpha\,\Delta u + \text{(source }f)$ with $\alpha\neq 0$ <br> moving the source to the right-hand side <br> $\boxed{\Delta u=f}$ (the central Poisson solve).
- Q: Derive the heat equation from conservation plus Fick's law (step 1 of 2). <br> <b>Formula:</b> conservation $\partial u/\partial t=-\nabla\!\cdot\mathbf{q}$; Fick/Fourier $\mathbf{q}=-\alpha\nabla u$ <br> <i>Hint:</i> substitute the flux law into conservation and use $\nabla\!\cdot\nabla=\Delta$. :: A: Conservation: $\partial u/\partial t=-\nabla\!\cdot\mathbf{q}$ <br> Fick/Fourier law: flux down the gradient, $\mathbf{q}=-\alpha\nabla u$ <br> substitute: $\partial u/\partial t=\nabla\!\cdot(\alpha\nabla u)$ <br> $\boxed{\partial u/\partial t=\alpha\,\Delta u}$ since $\nabla\!\cdot\nabla=\Delta$.
- Q: Derive the first-order state-space form of the spring-mass ODE $m\ddot x=-kx$. <br> <b>Formula:</b> state $\mathbf y=[x,v]^\top$; rows $\dot x=v$ and $\dot v=-\tfrac{k}{m}x$ <br> <i>Hint:</i> let $v=\dot x$ be the new variable, then solve Newton's law for $\dot v=\ddot x$. :: A: Define the state $\mathbf y=[x,\dot x]^\top=[x,v]^\top$ <br> the first row is the definition $\dot x=v$ <br> the second row is the law solved for acceleration: $\dot v=-\tfrac{k}{m}x$ <br> $\boxed{\dot{\mathbf y}=\begin{bmatrix}v\\ -\tfrac{k}{m}x\end{bmatrix},\ \mathbf y(0)=[x_0,v_0]^\top}$.
- Q: Derive the closed form of $\dot y=-ky$ by separation of variables. <br> <b>Formula:</b> $\dfrac{dy}{y}=-k\,dt$ (target: $y(t)=y_0 e^{-kt}$) <br> <i>Hint:</i> integrate both sides, then exponentiate and fix the constant with $y(0)=y_0$. :: A: Separate: $\dfrac{dy}{y}=-k\,dt$ <br> integrate both sides: $\ln y=-kt+C$ <br> exponentiate and fix $C$ with $y(0)=y_0$ <br> $\boxed{y(t)=y_0 e^{-kt}}$.
- Q: Worked example (small step). For $\dot y=-2y,\ y(0)=1$, take ONE forward-Euler step with $h=0.1$ and compare to exact $e^{-0.2}\approx0.8187$. <br> <b>Formula:</b> $y_1=y_0+h\,f(t_0,y_0)$ <br> <i>Hint:</i> compute the slope $f(0,1)=-2y_0$ first, then plug into the update. :: A: $f(0,1)=-2\cdot 1=-2$ <br> $y_1=1+0.1\cdot(-2)=0.8$ <br> exact $e^{-0.2}\approx0.8187$ <br> $\boxed{y_1=0.8,\ \text{error}\approx0.019}$ (small step tracks the curve).
- Q: Worked example (ten steps). For $\dot y=-2y,\ y(0)=1$, integrate to $t=1$ with forward Euler at $h=0.1$. <br> <b>Formula:</b> $y_{n+1}=(1+h\lambda)\,y_n$ with $\lambda=-2$ <br> <i>Hint:</i> compute the per-step factor $1+h\lambda$ first, then raise it to the power of 10 steps. :: A: Each step multiplies by $(1+h\cdot(-2))=(1-0.2)=0.8$ <br> ten steps: $y_{10}=0.8^{10}$ <br> $0.8^{10}\approx0.107$ <br> $\boxed{y_{10}\approx0.107\ \text{vs exact }e^{-2}\approx0.135}$.
- Q: Worked example (large step, instability). For $\dot y=-2y,\ y(0)=1$, take forward-Euler steps with $h=0.5$ to $t=1$. <br> <b>Formula:</b> $y_{n+1}=(1+h\lambda)\,y_n$ with $\lambda=-2$; stable iff $|1+h\lambda|<1$ <br> <i>Hint:</i> compute $1+h\lambda$ first; if it lands at $0$ or beyond $\pm1$ the step is too large. :: A: Step factor $(1+0.5\cdot(-2))=(1-1)=0$ <br> $y_1=1\cdot 0=0$, then $y_2=0$ <br> the step is too large: $|1+h\lambda|=1\not<1$ at the stability edge <br> $\boxed{y_2=0\ \text{(collapsed; useless, exact is }0.135)}$.
- Q: Worked example (3D robot regime). A joint state $\mathbf y=[q,\dot q]$ with $q=0.20$ rad, $\dot q=0$, under $\ddot q=-9\,q$, takes ONE forward-Euler step at $h=0.01$ s. <br> <b>Formula:</b> $q_{new}=q+h\,\dot q$, $\quad\dot q_{new}=\dot q+h\,\ddot q$ <br> <i>Hint:</i> update the position row with the current $\dot q$ first, then the velocity row with $\ddot q=-9q$. :: A: $\dot q$-row: $q_{new}=0.20+0.01\cdot 0=0.20$ <br> $\ddot q$-row: $\dot q_{new}=0+0.01\cdot(-9\cdot0.20)=-0.018$ <br> $\boxed{\mathbf y_1=[0.20,\,-0.018]}$ rad, rad/s.
- Q: Test yourself (limiting case). Show forward Euler reduces to the exact solution as $h\to0$ for $\dot y=-2y$. <br> <b>Formula:</b> Euler $y(h)=y_0(1-2h)$; exact $y_0 e^{-2h}=y_0(1-2h+2h^2-\dots)$ <br> <i>Hint:</i> Taylor-expand the exact $e^{-2h}$ and compare term by term; a pass means they agree to first order. :: A: One step gives $y(h)=y_0(1-2h)$ <br> compare to exact $y_0 e^{-2h}=y_0(1-2h+2h^2-\dots)$ <br> the two agree to first order; difference is $O(h^2)$ <br> $\boxed{\lim_{h\to0}\text{Euler}=\text{exact}}$ (consistency check passes).
- Q: Test yourself (stability). Backward Euler on $\dot y=-2y$ gives $y_{n+1}=y_n/(1+2h)$. Verify it cannot blow up for any $h>0$. <br> <b>Formula:</b> $y_{n+1}=y_n+h(-2y_{n+1})$; amplification factor $r=1/(1+2h)$ <br> <i>Hint:</i> solve for $y_{n+1}$, then check $|r|<1$; a pass means the factor stays in $(0,1)$. :: A: Solve $y_{n+1}=y_n+h(-2y_{n+1})\Rightarrow y_{n+1}(1+2h)=y_n$ <br> amplification factor $|1/(1+2h)|<1$ for all $h>0$ <br> so iterates decay monotonically like the true solution <br> $\boxed{\text{unconditionally stable: }0<\tfrac{1}{1+2h}<1}$.
- Q: Test yourself (period check). For the spring-mass state $\dot{\mathbf y}=[v,-\tfrac{k}{m}x]$ with $m=2,k=8$, predict the period and check units. <br> <b>Formula:</b> $\omega=\sqrt{k/m}$, $\quad T=2\pi/\omega$ <br> <i>Hint:</i> compute $\omega$ first, then $T$; a pass means $\sqrt{(\text{N/m})/\text{kg}}$ has units $1/\text{s}$ so $T$ is in seconds. :: A: Natural frequency $\omega=\sqrt{k/m}=\sqrt{8/2}=2$ rad/s <br> period $T=2\pi/\omega=2\pi/2=\pi\approx3.14$ s <br> units: $\sqrt{(\text{N/m})/\text{kg}}=\sqrt{1/\text{s}^2}=1/\text{s}$, so $T$ is in seconds <br> $\boxed{T=\pi\ \text{s, units consistent}}$.
- Q: When do you reach for an implicit method (backward Euler / BDF) over RK4? <br> <i>Hint:</i> think about which systems force explicit methods into tiny steps. :: A: For stiff systems that mix very fast and slow dynamics, where explicit methods need impractically tiny steps to stay stable.
- Q: Why does RK4 with a coarse step often beat Euler with a fine one (cost vs order)? <br> <b>Formula:</b> RK4 error $O(h^4)$ at 4 evals/step; Euler error $O(h)$ <br> <i>Hint:</i> compare how much each method's error drops when $h$ is halved. :: A: RK4 is $O(h^4)$, so halving $h$ cuts error 16x for only 4 evals/step; Euler is $O(h)$ and needs many more steps for the same accuracy.
- Q: Pitfall: a PDE solve gives a non-unique answer. What is the likely cause? <br> <i>Hint:</i> a PDE needs data on the domain boundary to be well-posed. :: A: Missing or wrong boundary conditions; a PDE needs Dirichlet (fixed value) or Neumann (fixed normal derivative) data to be well-posed.
- Q: Pitfall: a long rigid-body simulation slowly gains energy. What integrator fixes it? <br> <i>Hint:</i> name the integrator class that conserves phase-space volume. :: A: A symplectic integrator, which conserves energy/phase-space volume and avoids drift over many steps.
- CLOZE: Forward Euler is {{c1::$O(h)$}} (first-order) accurate, while RK4 is {{c2::$O(h^4)$}} (fourth-order) accurate.
- CLOZE: A method is stable when errors {{c1::do not blow up}}; explicit methods are only {{c2::conditionally}} stable, while backward Euler is {{c3::unconditionally}} stable.
- CLOZE: On a mesh the Laplacian $\Delta$ becomes the {{c1::cotangent Laplacian}} matrix and the heat/Poisson PDE becomes a {{c2::sparse linear system}}.
- CLOZE: Short-time heat diffusion $\partial u/\partial t=\alpha\Delta u$ is the first step of the {{c1::Signed Heat}} method and of the heat method for {{c2::geodesic distance}}.
- CLOZE: A higher-order ODE is reduced to first order by stacking the variable and its derivatives into a {{c1::state vector}}, e.g. $[x,\dot x]$.
