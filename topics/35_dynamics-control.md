# Dynamics and control: robot equations of motion, trajectories, and PID

Placement: [[15_inverse-kinematics]], [[34_differential-equations]] => **Dynamics and control** => (robotic fabrication: place, scribe, carve)

## Why it matters for the stone scan->CAD->IFC / CADFit / robotic-fabrication workflow

Kinematics tells a robot *where* to put its tool; dynamics tells it *what forces* are needed to get there smoothly, and control keeps it on track despite gravity, payload, and disturbances. Lifting and placing a heavy stone block, scribing a contact face, or dragging a carving tool all involve real forces, accelerations, and contact, not just static poses. A pick-and-place that ignores dynamics will overshoot, vibrate, or stall under load. This node closes the robotics branch of the master graph: kinematics ([[13_forward-kinematics]] to [[15_inverse-kinematics]]) plus the integrators of [[34_differential-equations]] become a controller that drives an arm.

## Intuition

Think of steering a loaded wheelbarrow. Kinematics is the path you want; dynamics is the heft you feel because of mass and slope; control is your hands constantly correcting to stay on the path. A robot controller does the same: it predicts the torques needed (feedforward from the model) and corrects the leftover error (feedback from sensors). Trajectory generation is planning a smooth speed-up and slow-down so nothing jerks.

## The math (first principles)

**Equations of motion.** For an $n$-joint arm with joint positions $q$, the rigid-body dynamics are

$$ M(q)\,\ddot q + C(q,\dot q)\,\dot q + g(q) = \tau, $$

where $M(q)$ is the symmetric positive-definite mass / inertia matrix, $C(q,\dot q)\dot q$ collects Coriolis and centrifugal terms, $g(q)$ is the gravity torque, and $\tau$ is the vector of joint torques (the control input). Add a friction term and an external wrench when modeling contact.

**Forward vs inverse dynamics.**
- Forward dynamics: given $\tau$, solve $\ddot q = M(q)^{-1}\big(\tau - C\dot q - g\big)$, then integrate with RK4 from [[34_differential-equations]] to simulate motion.
- Inverse dynamics: given a desired $q,\dot q,\ddot q$, compute the required $\tau$ directly from the equation above. This is the feedforward term.

**Formulations.** The same equation comes from the **Lagrangian** $L = T - V$ (kinetic minus potential energy) via $\frac{d}{dt}\frac{\partial L}{\partial \dot q} - \frac{\partial L}{\partial q} = \tau$, or recursively from the **Newton-Euler** force/torque balance link by link (efficient, $O(n)$).

**Trajectory generation.** A motion is a time-parameterized path $q(t)$ with smooth velocity and acceleration. A cubic or quintic polynomial interpolates between waypoints with matched boundary velocities; a **trapezoidal** (or S-curve) velocity profile accelerates at a constant rate, cruises, then decelerates, bounding velocity and acceleration so the arm does not jerk.

**Feedback control: PID.** With error $e(t) = q_{\text{desired}}(t) - q(t)$, the PID law is

$$ \tau = K_p\,e + K_i \int e\,dt + K_d\,\dot e. $$

The proportional term pulls toward the target, the integral term removes steady-state offset, and the derivative term damps oscillation. **Computed-torque (feedforward + feedback)** control combines the model and the corrector:

$$ \tau = \underbrace{M(q)\,\ddot q_d + C(q,\dot q)\dot q + g(q)}_{\text{feedforward (inverse dynamics)}} + \underbrace{M(q)\big(K_p e + K_d \dot e\big)}_{\text{feedback}}, $$

which linearizes the closed-loop error dynamics so the gains behave predictably.

**Stability.** Gains must keep the closed loop stable; too much $K_p$ causes oscillation, too much $K_d$ amplifies sensor noise. The error dynamics $\ddot e + K_d \dot e + K_p e = 0$ should be critically or slightly over-damped.

## Worked example

A single joint modeled as inertia $J=2$ with gravity load $g(q)=5\sin q$. To hold $q_d = \pi/2$ at rest ($\dot q=\ddot q=0$), inverse dynamics gives the feedforward torque $\tau_{ff}=J\cdot 0 + 5\sin(\pi/2)=5$. Suppose the joint sits at $q=1.4$ rad with $\dot q=0$; error $e=\pi/2-1.4\approx 0.171$. With $K_p=20,\ K_d=8$, the feedback torque is $20(0.171)+8(0)=3.42$, so total commanded $\tau=5+3.42=8.42$. The proportional term drives the joint up toward $\pi/2$ while the gravity feedforward already balances the weight, so the joint settles without a steady offset.

## Code

Idiomatic C++ with Eigen: forward-dynamics simulation of an arm with a computed-torque tracking controller, integrated by RK4 (reuse `rk4Step` from [[34_differential-equations]]).

```cpp
#include <Eigen/Dense>
using Vec = Eigen::VectorXd;
using Mat = Eigen::MatrixXd;

struct Arm {
    // Model terms (stub these with your robot's M, C, g).
    Mat M(const Vec& q) const;                 // inertia matrix, SPD
    Vec Cqd(const Vec& q, const Vec& qd) const;// Coriolis/centrifugal * qd
    Vec g(const Vec& q) const;                 // gravity torque
};

// Computed-torque control: feedforward inverse dynamics + PD feedback.
Vec computedTorque(const Arm& a, const Vec& q, const Vec& qd,
                   const Vec& qd_des, const Vec& qdd_des, const Vec& q_des,
                   const Mat& Kp, const Mat& Kd) {
    Vec e   = q_des  - q;
    Vec ed  = qd_des - qd;
    Vec ff  = a.M(q) * qdd_des + a.Cqd(q, qd) + a.g(q);  // inverse dynamics
    Vec fb  = a.M(q) * (Kp * e + Kd * ed);               // PD, mass-shaped
    return ff + fb;                                      // joint torques tau
}

// Forward dynamics: qddot = M^{-1} (tau - C qd - g).  State y = [q; qd].
Vec dynamicsRHS(const Arm& a, const Vec& q, const Vec& qd, const Vec& tau) {
    return a.M(q).ldlt().solve(tau - a.Cqd(q, qd) - a.g(q));   // = qddot
}
```

## Connections

- [[15_inverse-kinematics]] provides the joint targets $q_d(t)$ that this controller tracks.
- [[34_differential-equations]] supplies the RK4 integrator that turns torques into motion (forward dynamics).
- [[14_manipulator-jacobian]] maps joint to task space; needed for task-space (operational-space) control and for contact forces.
- [[07_calculus-jacobians]] and [[33_optimization-foundations]] underlie trajectory optimization (minimize jerk/energy subject to limits).
- [[33_optimization-foundations]] also frames gain tuning and optimal control as optimization.

## Pitfalls and failure modes

- **Kinematics is not enough.** Commanding poses without torque/acceleration limits causes overshoot and vibration under load.
- **Gain tuning.** Too high $K_p$ oscillates; too high $K_d$ amplifies measurement noise; integral windup saturates actuators (clamp the integral term).
- **Model error.** Wrong mass or gravity makes feedforward fight the robot; identify the model or lean more on feedback.
- **Non-smooth trajectories.** Position waypoints without matched velocities create jerk; use cubic/quintic or trapezoidal profiles.
- **Singularities in task-space control.** Near a Jacobian singularity, required joint velocities explode; damp them (same damped least squares as IK).
- **Contact and stiffness.** Carving involves contact forces; pure position control is unsafe, prefer impedance/force control.

## Practice

1. Simulate a single pendulum joint ($M=J$, $g(q)=mgl\sin q$) with `dynamicsRHS` + RK4; verify it swings at the small-angle period $2\pi\sqrt{J/(mgl)}$.
2. Add PD control to hold $q_d=\pi/2$; tune $K_p,K_d$ for critical damping of $\ddot e + K_d\dot e + K_p e=0$.
3. Generate a trapezoidal velocity profile between two joint angles with given max velocity and acceleration; plot $q,\dot q,\ddot q$.
4. Show that computed-torque control reduces the closed-loop error dynamics to $\ddot e + K_d\dot e + K_p e=0$.
5. Introduce a 20% mass error in the feedforward and measure the resulting steady-state tracking error; then add an integral term and show it vanishes.

## References

- Kevin Lynch and Frank Park, *Modern Robotics* (Ch. 8 dynamics, Ch. 11 control): https://hades.mech.northwestern.edu/index.php/Modern_Robotics
- Siciliano, Sciavicco, Villani, Oriolo, *Robotics: Modelling, Planning and Control*.
- Haviland and Corke, "Manipulator Differential Kinematics" tutorial (Parts 1-2).
- Spong, Hutchinson, Vidyasagar, *Robot Modeling and Control*.

## Flashcard seeds

- Q: Derive forward dynamics: solve the equations of motion for the joint acceleration $\ddot q$. <br> <b>Formula:</b> $M(q)\ddot q + C(q,\dot q)\dot q + g(q)=\tau$ (target: $\ddot q=M^{-1}(\tau-C\dot q-g)$) <br> <i>Hint:</i> move the $C\dot q$ and $g$ terms to the right, then left-multiply by $M^{-1}$ (it is SPD, hence invertible). :: A: Move non-inertial terms right: $M(q)\ddot q=\tau-C(q,\dot q)\dot q-g(q)$ <br> $M$ is SPD hence invertible, so left-multiply by $M^{-1}$ <br> $\boxed{\ddot q=M(q)^{-1}\big(\tau-C(q,\dot q)\dot q-g(q)\big)}$.
- Q: Rewrite the second-order dynamics as a first-order state ODE $\dot y=f(y)$ ready for RK4. <br> <b>Formula:</b> $\ddot q=M^{-1}(\tau-C\dot q-g)$, state $y=[q;\dot q]$ <br> <i>Hint:</i> the top block of $\dot y$ is just $\dot q$ (already a state); the bottom block is $\ddot q$. :: A: Choose state $y=[q;\dot q]$ <br> $\dot y_{\text{top}}=\dot q$ (already a state) <br> $\dot y_{\text{bot}}=\ddot q=M^{-1}(\tau-C\dot q-g)$ <br> $\boxed{\dot y=\begin{bmatrix}\dot q\\ M(q)^{-1}(\tau-C\dot q-g)\end{bmatrix}}$ then integrate $y$ with RK4.
- Q: Computed torque, step 1 of 3: pick $\tau$ so the closed loop becomes $\ddot q=a$ for a chosen input $a$. <br> <b>Formula:</b> $M\ddot q+C\dot q+g=\tau$; try $\tau=M\,a+C\dot q+g$ <br> <i>Hint:</i> substitute that $\tau$ into the EOM and cancel the matching $C\dot q+g$ on both sides. :: A: Pick $\tau=M\,a+C\dot q+g$ and substitute <br> $M\ddot q+C\dot q+g=M\,a+C\dot q+g$ <br> cancel $C\dot q+g$ and left-multiply by $M^{-1}$ <br> $\boxed{\ddot q=a}$ (feedback-linearized double integrator).
- Q: Computed torque, step 2 of 3: choose the input $a$ so the tracking error decays. <br> <b>Formula:</b> $\ddot q=a$, error $e=q_d-q$; try $a=\ddot q_d+K_d\dot e+K_p e$ <br> <i>Hint:</i> substitute, then use $\ddot e=\ddot q_d-\ddot q$ to eliminate $\ddot q$. :: A: Set $a=\ddot q_d+K_d\dot e+K_p e$ <br> substitute: $\ddot q=\ddot q_d+K_d\dot e+K_p e$ <br> use $\ddot e=\ddot q_d-\ddot q$ so $\ddot q_d-\ddot q=-(K_d\dot e+K_p e)$ <br> $\boxed{\ddot e+K_d\dot e+K_p e=0}$.
- Q: Computed torque, step 3 of 3: assemble the full torque and split it into feedforward and feedback. <br> <b>Formula:</b> $\tau=M\,a+C\dot q+g$ with $a=\ddot q_d+K_d\dot e+K_p e$ <br> <i>Hint:</i> substitute $a$, then group the $\ddot q_d,C\dot q,g$ terms (model) apart from the $K_p,K_d$ terms (corrector). :: A: Substitute $a$: $\tau=M(\ddot q_d+K_d\dot e+K_p e)+C\dot q+g$ <br> group feedforward (model) and feedback (corrector) <br> $\boxed{\tau=\underbrace{M\ddot q_d+C\dot q+g}_{\text{inverse dynamics}}+\underbrace{M(K_p e+K_d\dot e)}_{\text{feedback}}}$.
- Q: From the closed-loop error ODE, read off the natural frequency $\omega_n$ and damping ratio $\zeta$. <br> <b>Formula:</b> $\ddot e+K_d\dot e+K_p e=0$ vs the standard $\ddot e+2\zeta\omega_n\dot e+\omega_n^2 e=0$ <br> <i>Hint:</i> match the constant ($e$) coefficient first to get $\omega_n$, then the $\dot e$ coefficient. :: A: Compare to $\ddot e+2\zeta\omega_n\dot e+\omega_n^2 e=0$ <br> match coefficients: $\omega_n^2=K_p$ and $2\zeta\omega_n=K_d$ <br> $\boxed{\omega_n=\sqrt{K_p},\quad \zeta=\dfrac{K_d}{2\sqrt{K_p}}}$.
- Q: Derive the $K_d$ that gives critical damping. <br> <b>Formula:</b> $\zeta=K_d/(2\sqrt{K_p})$, critical means $\zeta=1$ <br> <i>Hint:</i> set $\zeta=1$ and solve for $K_d$. :: A: Critical damping means $\zeta=1$ <br> $1=K_d/(2\sqrt{K_p})$ <br> solve for $K_d$ <br> $\boxed{K_d=2\sqrt{K_p}}$ (e.g. $K_p=100\Rightarrow K_d=20$).
- Q: Derive the inverse-dynamics feedforward torque needed to hold a pose at rest. <br> <b>Formula:</b> $M\ddot q+C\dot q+g=\tau$ with $\dot q=\ddot q=0$ at rest <br> <i>Hint:</i> set both $\dot q$ and $\ddot q$ to zero and see which term survives. :: A: Set $\dot q=0,\ \ddot q=0$ in $M\ddot q+C\dot q+g=\tau$ <br> $M\cdot 0+C\cdot 0+g(q)=\tau$ <br> $\boxed{\tau_{ff}=g(q)}$ — at rest the feedforward only fights gravity.
- Q: Derive the steady-state error of a P controller when the gravity model is off by a constant bias $\Delta$. <br> <b>Formula:</b> at rest the residual proportional torque balances the bias: $K_p e_{ss}=\Delta$ <br> <i>Hint:</i> at steady state $\ddot q=\dot q=0$, so only $K_p e$ is left to supply the missing $\Delta$; solve for $e_{ss}$. :: A: At steady state $\ddot q=\dot q=0$, so the residual torque $K_p e$ balances the unmodeled bias $\Delta$ <br> $K_p e_{ss}=\Delta$ <br> $\boxed{e_{ss}=\Delta/K_p}$ (nonzero offset).
- Q: Show that adding an integral term drives the constant-bias offset $e_{ss}=\Delta/K_p$ to zero. <br> <b>Formula:</b> $\tau=\dots+K_i\int e\,dt$ <br> <i>Hint:</i> note a nonzero constant $e_{ss}$ makes $\int e\,dt$ ramp forever; ask what happens once $K_i\int e\,dt$ has grown to supply $\Delta$. :: A: A constant $e_{ss}$ makes $\int e\,dt$ ramp without bound, so $K_i\int e\,dt$ grows until it supplies $\Delta$ on its own <br> then the P term no longer needs error, so $e$ can return to $0$ <br> $\boxed{e_{ss}\to 0\ \text{(integral cancels the constant bias)}}$.
- Q: For a trapezoidal velocity profile, derive the distance covered during one acceleration ramp. <br> <b>Formula:</b> $t_a=v_{\max}/a_{\max}$ and $d_a=\tfrac12 a_{\max}t_a^2$ <br> <i>Hint:</i> compute $t_a$ first, then substitute it into the constant-acceleration distance. :: A: Time to reach cruise: $t_a=v_{\max}/a_{\max}$ <br> distance under constant accel: $d_a=\tfrac12 a_{\max}t_a^2$ <br> $\boxed{d_a=\dfrac{v_{\max}^2}{2a_{\max}}}$.
- Q: Derive the total move time of a trapezoidal profile over distance $D$. <br> <b>Formula:</b> $d_a=v_{\max}^2/(2a_{\max})$, $t_a=v_{\max}/a_{\max}$; cruise distance $d_c=D-2d_a$ at $v_{\max}$, total $T=2t_a+t_c$ <br> <i>Hint:</i> find the cruise time $t_c=d_c/v_{\max}$, then add the two ramps; the $d_a$ terms cancel. :: A: Two ramps cover $2d_a$; cruise distance $d_c=D-2d_a$ at $v_{\max}$ <br> $t_c=d_c/v_{\max}=D/v_{\max}-v_{\max}/a_{\max}$ <br> total $T=2t_a+t_c$ <br> $\boxed{T=\dfrac{D}{v_{\max}}+\dfrac{v_{\max}}{a_{\max}}}$ (valid when $D\ge v_{\max}^2/a_{\max}$).
- Q: Worked example (hold a pose, small joint). Single joint $J=2$, gravity $g(q)=5\sin q$, target $q_d=\pi/2$ at rest; joint at $q=1.4$ rad, $\dot q=0$, $K_p=20,K_d=8$. Find the commanded torque. <br> <b>Formula:</b> $\tau=\underbrace{g(q_d)}_{\tau_{ff}}+\underbrace{K_p e+K_d\dot e}_{\text{feedback}},\ e=q_d-q$ <br> <i>Hint:</i> compute the gravity feedforward $g(\pi/2)$ first, then the error $e=\pi/2-1.4$. :: A: Feedforward at rest: $\tau_{ff}=g(\pi/2)=5\sin(\pi/2)=5$ <br> error $e=\pi/2-1.4\approx0.1708$, $\dot e=0$ <br> feedback $K_p e+K_d\dot e=20(0.1708)+0=3.42$ <br> $\boxed{\tau=5+3.42=8.42}$.
- Q: Worked example (different regime: light joint, $\pi/3$). $J=0.5$, $g(q)=2\sin q$, $q_d=\pi/3$ at rest; $q=0.9$ rad, $\dot q=0$, $K_p=40,K_d=12$. Find $\tau$. <br> <b>Formula:</b> $\tau=g(q_d)+K_p e+K_d\dot e,\ e=q_d-q$ <br> <i>Hint:</i> compute $\tau_{ff}=2\sin(\pi/3)$ first, then $e=\pi/3-0.9$ ($\dot e=0$ so the $K_d$ term drops). :: A: Feedforward $\tau_{ff}=2\sin(\pi/3)=2(0.8660)=1.732$ <br> $e=\pi/3-0.9\approx0.1472$ <br> feedback $=40(0.1472)=5.888$ <br> $\boxed{\tau=1.732+5.888\approx7.62}$.
- Q: Worked example (critical damping, two gain sets). Tune $K_d$ for critical damping of $\ddot e+K_d\dot e+K_p e=0$ when (a) $K_p=20$ and (b) $K_p=100$. <br> <b>Formula:</b> $K_d=2\sqrt{K_p}$ <br> <i>Hint:</i> plug each $K_p$ straight into the formula; $\sqrt{20}\approx4.47$. :: A: Use $K_d=2\sqrt{K_p}$ <br> (a) $K_d=2\sqrt{20}\approx8.94$ <br> (b) $K_d=2\sqrt{100}=20$ <br> $\boxed{K_d=2\sqrt{K_p}:\ 8.94\ \text{and}\ 20}$.
- Q: Worked example (trapezoidal, m-scale move). Joint travels $D=3$ rad with $v_{\max}=2$ rad/s, $a_{\max}=4$ rad/s$^2$. Find the total time and check it reaches cruise. <br> <b>Formula:</b> cruise exists iff $D\ge v_{\max}^2/a_{\max}$; then $T=D/v_{\max}+v_{\max}/a_{\max}$ <br> <i>Hint:</i> first test the cruise condition $v_{\max}^2/a_{\max}$ against $D$, then plug into $T$. :: A: Reach-cruise distance $v_{\max}^2/a_{\max}=4/4=1\le3$, so a cruise phase exists <br> $T=D/v_{\max}+v_{\max}/a_{\max}=3/2+2/4=1.5+0.5$ <br> $\boxed{T=2.0\ \text{s}}$.
- Q: Worked example (mm-scale linear move). A scribe tip travels $D=120$ mm with $v_{\max}=60$ mm/s and $a_{\max}=300$ mm/s$^2$. Find the total move time. <br> <b>Formula:</b> cruise exists iff $D\ge v_{\max}^2/a_{\max}$; then $T=D/v_{\max}+v_{\max}/a_{\max}$ <br> <i>Hint:</i> check $v_{\max}^2/a_{\max}=3600/300$ against $D=120$ first, then evaluate $T$. :: A: Cruise check: $v_{\max}^2/a_{\max}=3600/300=12\,\text{mm}\le120$, so a cruise phase exists <br> $T=120/60+60/300=2.0+0.2$ <br> $\boxed{T=2.2\ \text{s}}$.
- Q: Test yourself (close the loop). For the $K_p=20$ joint above, the P controller leaves $e_{ss}=\Delta/K_p$ under a gravity bias $\Delta=1$. Compute $e_{ss}$ and verify it shrinks if you raise $K_p$ to $100$. <br> <b>Formula:</b> $e_{ss}=\Delta/K_p$ <br> <i>Hint:</i> a pass means $e_{ss}$ drops by the same factor $K_p$ rises (here $5\times$), and stays nonzero. :: A: $e_{ss}=1/20=0.05$ rad <br> raise gain: $e_{ss}'=1/100=0.01$ rad <br> $\boxed{0.05\to0.01\ \text{as }K_p:20\to100}$ — higher $K_p$ shrinks (but does not kill) the offset; only the integral drives it to $0$. Limit check: $K_p\to\infty\Rightarrow e_{ss}\to0$.
- Q: Test yourself (units + limiting case). Verify $T=D/v_{\max}+v_{\max}/a_{\max}$ is dimensionally a time, then check the limit $a_{\max}\to\infty$. <br> <b>Formula:</b> $T=D/v_{\max}+v_{\max}/a_{\max}$ (use rad, rad/s, rad/s$^2$) <br> <i>Hint:</i> a pass means each term reduces to seconds, and the second term vanishes as $a_{\max}\to\infty$. :: A: Units: $\frac{\text{rad}}{\text{rad/s}}+\frac{\text{rad/s}}{\text{rad/s}^2}=\text{s}+\text{s}=\text{s}$ (ok) <br> as $a_{\max}\to\infty$ the ramps vanish: $T\to D/v_{\max}$ <br> $\boxed{T\to D/v_{\max}}$ = the pure constant-velocity time, as expected.
- Q: Test yourself (forward-dynamics sanity). With $\tau=g(q)$ at rest ($\dot q=0$), plug into the forward dynamics and confirm equilibrium. <br> <b>Formula:</b> $\ddot q=M^{-1}(\tau-C\dot q-g)$ <br> <i>Hint:</i> a pass means $\ddot q=0$; note $C\dot q=0$ when $\dot q=0$, so $\tau-g$ should cancel. :: A: $C\dot q=0$ since $\dot q=0$ <br> $\ddot q=M^{-1}(g(q)-0-g(q))=M^{-1}\,0$ <br> $\boxed{\ddot q=0}$ — gravity-balancing torque holds the pose, confirming equilibrium.
- Q: When do you use computed-torque (inverse-dynamics) control instead of plain PID? <br> <i>Hint:</i> think about what extra ingredient computed torque needs, and what it buys for the error dynamics. :: A: When you have a good dynamics model and want the closed-loop error to follow $\ddot e+K_d\dot e+K_p e=0$ with model-independent, predictable gains (model cancels the nonlinearities).
- Q: Lagrangian route to the EOM: which energy expression and which operator give $\tau$? <br> <b>Formula:</b> $L=T-V$, operator $\frac{d}{dt}\frac{\partial L}{\partial\dot q}-\frac{\partial L}{\partial q}=\tau$ <br> <i>Hint:</i> $T$ is kinetic, $V$ is potential energy. :: A: $L=T-V$ (kinetic minus potential), with $\frac{d}{dt}\frac{\partial L}{\partial\dot q}-\frac{\partial L}{\partial q}=\tau$.
- Q: Pitfall — integral windup: what is it and one fix? <br> <i>Hint:</i> think about what the integrator does while the actuator is saturated and cannot deliver more torque. :: A: The integral keeps accumulating while the actuator is saturated, overshooting badly; fix with clamping / anti-windup on the integral term.
- Q: Pitfall — why is pure position control unsafe when carving stone? <br> <i>Hint:</i> think about what happens to force the instant a stiff tool meets a stiff stone surface. :: A: Carving creates contact forces; rigid position control spikes force on contact, so use impedance or force control.
- CLOZE: In PID, the {{c1::integral}} term removes steady-state error and the {{c2::derivative}} term damps oscillation.
- CLOZE: Computed-torque control linearizes the closed loop to {{c1::$\ddot e + K_d\dot e + K_p e = 0$}}, so the gains set $\omega_n=\sqrt{K_p}$ and $\zeta=K_d/(2\sqrt{K_p})$.
- CLOZE: Forward dynamics solves $\ddot q = {{c1::M(q)^{-1}(\tau - C\dot q - g)}}$ and is integrated with {{c2::RK4}}. <br> <i>hint:</i> isolate $\ddot q$ from $M\ddot q+C\dot q+g=\tau$.
- CLOZE: At rest the inverse-dynamics feedforward reduces to {{c1::$\tau_{ff}=g(q)$}}, just balancing gravity.
- Q: Code idiom: which Eigen call solves $M\ddot q=b$ for an SPD inertia matrix, and why that one? <br> <i>Hint:</i> the factorization name reflects "symmetric positive-definite". :: A: `M.ldlt().solve(b)` — LDL$^T$ is the fast, stable factorization for a symmetric positive-definite $M$.
