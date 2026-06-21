# State estimation and filtering (Kalman / EKF, sensor fusion)

Placement: [[32_probability-statistics]], [[06_linear-least-squares]] => **State estimation & filtering** => (live tool pose, sensor fusion)

## Why it matters for the stone scan->CAD->IFC / CADFit / robotic-fabrication workflow

During carving the robot must know where its tool actually is, right now, despite noisy sensors and drift. No single source is trusted: joint encoders drift, a scanner is noisy, a turntable wobbles. A filter fuses them into one best estimate with a quantified uncertainty, updated every cycle. This is the live-pose backbone that closes the loop between the scan/CAD world and the moving machine, and it sits directly under [[37_trajectory-planning-collision]] and [[35_dynamics-control]].

## Intuition

You have a prediction of where the tool is (from the motion model) and a measurement (from a sensor). Neither is exact. The filter forms a weighted blend: trust whichever is more certain. The weight is the gain. After fusing, the result is more certain than either input alone. See the hero diagram: a wide prior Gaussian, a measurement Gaussian, and a narrower posterior that lands between them, pulled toward the more certain one.

## The math (first principles)

The **Bayes filter** alternates predict and update on a belief $p(x)$. For linear systems with Gaussian noise it has a closed form, the **Kalman filter**, with state $x$, covariance $P$, motion model $F$, control $u$ via $B$, process noise $Q$, measurement model $H$, measurement noise $R$.

**Predict** (push the belief through the motion model, uncertainty grows):

$$ x^- = F x + B u, \qquad P^- = F P F^\top + Q. $$

**Update** (fuse the measurement $z$, uncertainty shrinks):

$$ \tilde y = z - H x^-, \quad S = H P^- H^\top + R, \quad K = P^- H^\top S^{-1}, $$
$$ x^+ = x^- + K\tilde y, \qquad P^+ = (I - K H) P^-. $$

$K$ is the **Kalman gain**: it is the ratio of prediction uncertainty to total uncertainty, so a confident prediction ($P^-$ small) yields small $K$ (trust the model), a confident sensor ($R$ small) yields large $K$ (trust the measurement). In **1-D** ($H=1$) this collapses to the picture in the diagram:

$$ K = \frac{P^-}{P^- + R}, \qquad x^+ = x^- + K(z - x^-), \qquad P^+ = (1-K)P^-. $$

For nonlinear motion $f$ or measurement $h$, the **Extended Kalman Filter (EKF)** linearizes with Jacobians $F = \partial f/\partial x$ and $H = \partial h/\partial x$ (the link to [[07_calculus-jacobians]]) and uses the same update. The Kalman filter is also **recursive least squares**: it is the running MLE under Gaussian noise, tying back to [[06_linear-least-squares]] and [[32_probability-statistics]].

## Worked example

1-D position. Prediction $x^- = 10$ with variance $P^- = 4$. A sensor reads $z = 12$ with variance $R = 1$ (more certain than the prediction). Gain $K = \frac{4}{4+1} = 0.8$. Updated estimate $x^+ = 10 + 0.8(12-10) = 11.6$. Updated variance $P^+ = (1-0.8)\cdot 4 = 0.8$. The estimate moved most of the way to the measurement, and the variance dropped to $0.8$, below both $4$ and $1$: fusing two noisy numbers beats either alone.

## Code

Idiomatic C++ with Eigen: one Kalman predict and update step. The same shapes hold for a 6-DoF tool pose.

```cpp
#include <Eigen/Dense>
using Vec = Eigen::VectorXd;
using Mat = Eigen::MatrixXd;

struct Kalman {
    Vec x;   // state estimate
    Mat P;   // state covariance

    void predict(const Mat& F, const Mat& Q, const Vec& u, const Mat& B) {
        x = F * x + B * u;             // x- = Fx + Bu
        P = F * P * F.transpose() + Q; // P- = FPF^T + Q  (uncertainty grows)
    }
    void update(const Vec& z, const Mat& H, const Mat& R) {
        Vec y = z - H * x;                          // innovation
        Mat S = H * P * H.transpose() + R;          // innovation covariance
        Mat K = P * H.transpose() * S.inverse();    // Kalman gain
        x = x + K * y;                              // x+ = x- + K y
        Mat I = Mat::Identity(P.rows(), P.cols());
        P = (I - K * H) * P;                        // P+ = (I - KH) P-  (shrinks)
    }
};
```

## Connections

- [[32_probability-statistics]] supplies the Gaussian and the Bayes update the filter implements.
- [[06_linear-least-squares]] is the static version; the Kalman filter is recursive least squares over time.
- [[07_calculus-jacobians]] provides the $F$, $H$ Jacobians that the EKF needs for nonlinear models.
- [[35_dynamics-control]] gives the motion model $F$ used in predict; control consumes the filtered state.
- [[16_registration-icp]] and [[11_calibration-pnp]] produce pose measurements $z$ to fuse.

## Pitfalls and failure modes

- **Tuning $Q$ and $R$.** Too small $R$ makes the filter chase noise; too small $Q$ makes it ignore real motion and diverge. They encode how much you trust sensor vs model.
- **Divergence.** Linearization error (EKF) or an over-confident $P$ can make the filter confidently wrong. Watch the innovation $\tilde y$; persistent large innovations mean the model is off.
- **Non-PSD covariance.** Round-off can make $P$ lose positive-definiteness; use the Joseph form $P^+ = (I-KH)P^-(I-KH)^\top + KRK^\top$ for numerical stability.
- **Observability.** If the measurements do not constrain a state direction, its uncertainty grows without bound.
- **Synchronization.** Fusing measurements with wrong timestamps injects error; timestamp and align first.

## Practice

1. Implement the 1-D filter and reproduce $x^+ = 11.6$, $P^+ = 0.8$ from the worked example.
2. Simulate a constant-velocity tool ($x = [\text{pos}, \text{vel}]$) with noisy position measurements; show the estimate tracks truth and $P$ converges.
3. Vary $R$ from small to large and plot how the gain $K$ and tracking lag change.
4. Build an EKF for a nonlinear measurement (range to a fixed point) and compare to the linear KF.
5. Break it: feed an unmodeled acceleration and watch the innovations blow up; then raise $Q$ to recover.

## References

- Welch and Bishop, "An Introduction to the Kalman Filter": https://www.cs.unc.edu/~welch/kalman/
- Thrun, Burgard, Fox, *Probabilistic Robotics* (Bayes filter, KF, EKF chapters).
- Roger Labbe, *Kalman and Bayesian Filters in Python* (free): https://github.com/rlabbe/Kalman-and-Bayesian-Filters-in-Python
- Simon, *Optimal State Estimation*.

## Flashcard seeds

- Q: Derive the 1-D Kalman gain from the general gain by setting the measurement model to identity. <br> <b>Formula:</b> general $K = P^- H^\top (H P^- H^\top + R)^{-1}$; target $K = \dfrac{P^-}{P^- + R}$ <br> <i>Hint:</i> substitute $H = 1$ (scalar) so every $H$ and $H^\top$ becomes $1$, then simplify. :: A: Set $H=1$: $K = P^-\cdot 1\cdot(1\cdot P^-\cdot 1 + R)^{-1}$ <br> $= P^- (P^- + R)^{-1}$ <br> $\boxed{K = \dfrac{P^-}{P^- + R}}$
- Q: Derive the 1-D covariance update from the general one with $H=1$. <br> <b>Formula:</b> general $P^+ = (I - KH)P^-$; target $P^+ = (1-K)P^-$ <br> <i>Hint:</i> replace $I$ with $1$ and $H$ with $1$ in the scalar case. :: A: With $H=1$, $I=1$: $P^+ = (1 - K\cdot 1)P^-$ <br> $\boxed{P^+ = (1-K)P^-}$
- Q: Show that with $H=1$ the general state update becomes the 1-D form. <br> <b>Formula:</b> general $x^+ = x^- + K(z - Hx^-)$; target $x^+ = x^- + K(z - x^-)$ <br> <i>Hint:</i> set $H=1$ inside the innovation $z - Hx^-$. :: A: Innovation $\tilde y = z - Hx^- = z - 1\cdot x^- = z - x^-$ <br> $x^+ = x^- + K\tilde y$ <br> $\boxed{x^+ = x^- + K(z - x^-)}$
- Q: From the predict–update diagram, write the two Kalman predict equations. <br> <b>Formula:</b> $x^- = Fx + Bu$, $\;P^- = FPF^\top + Q$ <br> <i>Hint:</i> one line pushes the state through the motion model, the other grows the covariance by adding $Q$. :: A: State: $x^- = Fx + Bu$ <br> Covariance: $P^- = FPF^\top + Q$ <br> $\boxed{x^- = Fx + Bu,\quad P^- = FPF^\top + Q}$
- Q: Write the full Kalman update sequence (innovation, its covariance, gain, state, covariance). <br> <b>Formula:</b> $\tilde y = z - Hx^-$, $\;S = HP^-H^\top + R$, $\;K = P^-H^\top S^{-1}$, $\;x^+ = x^- + K\tilde y$, $\;P^+ = (I-KH)P^-$ <br> <i>Hint:</i> compute $\tilde y$ and $S$ first, then $K$ from them, then update $x$ and $P$. :: A: $\tilde y = z - Hx^-$ <br> $S = HP^-H^\top + R$ <br> $K = P^-H^\top S^{-1}$ <br> $x^+ = x^- + K\tilde y$ <br> $\boxed{P^+ = (I-KH)P^-}$
- Q: Worked 1-D fuse. Prediction $x^-=10$, $P^-=4$; sensor $z=12$, $R=1$. Compute $K$, $x^+$, $P^+$. <br> <b>Formula:</b> $K = \dfrac{P^-}{P^-+R}$, $\;x^+ = x^- + K(z-x^-)$, $\;P^+ = (1-K)P^-$ <br> <i>Hint:</i> compute the gain $K$ first, then plug it into the state and covariance updates. :: A: $K = \dfrac{4}{4+1} = 0.8$ <br> $x^+ = 10 + 0.8(12-10) = 10 + 1.6 = 11.6$ <br> $P^+ = (1-0.8)\cdot 4 = 0.8$ <br> $\boxed{K=0.8,\; x^+=11.6,\; P^+=0.8}$
- Q: Worked 1-D fuse at mm scale with a noisier sensor. Prediction $x^-=50$ mm, $P^-=2$; sensor $z=53$ mm, $R=6$. Compute $K$ and $x^+$. <br> <b>Formula:</b> $K = \dfrac{P^-}{P^-+R}$, $\;x^+ = x^- + K(z-x^-)$ <br> <i>Hint:</i> $R$ is larger than $P^-$ here, so expect $K<0.5$ — the estimate stays near the prediction. :: A: $K = \dfrac{2}{2+6} = 0.25$ <br> $x^+ = 50 + 0.25(53-50) = 50 + 0.75 = 50.75$ <br> $\boxed{K=0.25,\; x^+=50.75\text{ mm}}$
- Q: Worked 1-D fuse at metre scale, very confident sensor. Prediction $x^-=3.0$ m, $P^-=0.9$; sensor $z=3.4$ m, $R=0.1$. Compute $K$, $x^+$, $P^+$. <br> <b>Formula:</b> $K = \dfrac{P^-}{P^-+R}$, $\;x^+ = x^- + K(z-x^-)$, $\;P^+ = (1-K)P^-$ <br> <i>Hint:</i> small $R$ means large $K$ — the estimate jumps most of the way to $z$. :: A: $K = \dfrac{0.9}{0.9+0.1} = 0.9$ <br> $x^+ = 3.0 + 0.9(3.4-3.0) = 3.0 + 0.36 = 3.36$ <br> $P^+ = (1-0.9)\cdot 0.9 = 0.09$ <br> $\boxed{K=0.9,\; x^+=3.36\text{ m},\; P^+=0.09}$
- Q: Worked predict-only step (2-D covariance grows). $P = \begin{psmallmatrix}1&0\\0&1\end{psmallmatrix}$, $F = I$, process noise $Q = \begin{psmallmatrix}0.2&0\\0&0.2\end{psmallmatrix}$. Compute $P^-$. <br> <b>Formula:</b> $P^- = FPF^\top + Q$ <br> <i>Hint:</i> with $F=I$ the term $FPF^\top = P$, so just add $Q$ entrywise. :: A: $FPF^\top = I\,P\,I = P = \begin{psmallmatrix}1&0\\0&1\end{psmallmatrix}$ <br> $P^- = P + Q = \begin{psmallmatrix}1.2&0\\0&1.2\end{psmallmatrix}$ <br> $\boxed{P^- = \begin{psmallmatrix}1.2&0\\0&1.2\end{psmallmatrix}}$ (uncertainty grew)
- Q: Closing test — does the fused variance beat both inputs? Use the worked numbers $P^-=4$, $R=1$, $P^+=0.8$. <br> <b>Formula:</b> check $P^+ < \min(P^-, R)$ <br> <i>Hint:</i> a pass means the posterior variance is below BOTH $P^-$ and $R$; compute $\min(4,1)$ and compare. :: A: $\min(P^-, R) = \min(4, 1) = 1$ <br> Is $0.8 < 1$? Yes <br> $\boxed{\text{Pass: } P^+=0.8 < 1 \text{ — fusing beats either alone}}$
- Q: Closing test — is the Kalman gain a valid weight? Use $P^-=4$, $R=1$. <br> <b>Formula:</b> $K = \dfrac{P^-}{P^-+R}$, valid iff $0 \le K \le 1$ <br> <i>Hint:</i> compute $K$ and confirm it lies in $[0,1]$; $K\to1$ as $R\to0$, $K\to0$ as $P^-\to0$. :: A: $K = \dfrac{4}{4+1} = 0.8$ <br> Is $0 \le 0.8 \le 1$? Yes <br> $\boxed{\text{Pass: } K=0.8 \in [0,1]}$
- Q: Closing test — should the covariance stay symmetric positive-definite? Given $P^+ = (I-KH)P^-$ lost symmetry to round-off, name the fix and write it. <br> <b>Formula:</b> Joseph form $P^+ = (I-KH)P^-(I-KH)^\top + KRK^\top$ <br> <i>Hint:</i> a pass is a $P$ that stays symmetric PSD; the Joseph form is symmetric by construction. :: A: Replace the plain update with the Joseph form <br> $P^+ = (I-KH)P^-(I-KH)^\top + KRK^\top$ <br> Each term is symmetric, so $P^+$ stays symmetric PSD <br> $\boxed{P^+ = (I-KH)P^-(I-KH)^\top + KRK^\top}$
- Q: How does the EKF handle nonlinear motion $f$ or measurement $h$? <br> <b>Formula:</b> $F = \dfrac{\partial f}{\partial x}$, $\;H = \dfrac{\partial h}{\partial x}$ <br> <i>Hint:</i> linearize, then reuse the standard KF predict/update with these Jacobians. :: A: Linearize $f,h$ at the current estimate using Jacobians $F=\partial f/\partial x$, $H=\partial h/\partial x$ <br> Then apply the ordinary KF predict and update <br> $\boxed{\text{EKF = KF on the linearized model}}$
- Q: What does a filter do that no single sensor can? <br> <i>Hint:</i> think about the posterior variance versus each input variance. :: A: It fuses prediction and measurement into one estimate that is more certain than either input alone <br> $\boxed{P^+ < \min(P^-, R)}$
- Q: The Kalman filter is the recursive form of which static estimator? <br> <i>Hint:</i> it runs the same Gaussian MLE one step at a time. :: A: Least squares <br> $\boxed{\text{recursive least squares (running Gaussian MLE)}}$
- Q: In predict vs update, when does uncertainty grow and when does it shrink? <br> <b>Formula:</b> $P^- = FPF^\top + Q$ (predict), $\;P^+ = (I-KH)P^-$ (update) <br> <i>Hint:</i> look at which step adds $Q$ and which multiplies by $(I-KH)$. :: A: Predict adds $Q$, so uncertainty grows <br> Update multiplies by $(I-KH)$ with $0\le K\le1$, so it shrinks <br> $\boxed{\text{grows in predict, shrinks in update}}$
- Q: Symptom and meaning of persistently large innovations $\tilde y = z - Hx^-$? <br> <b>Formula:</b> $\tilde y = z - Hx^-$ <br> <i>Hint:</i> the innovation is what the measurement says minus what the model predicted. :: A: Large, persistent $\tilde y$ means the prediction keeps missing the measurement <br> $\boxed{\text{model is wrong / mistuned — the filter is diverging}}$
- CLOZE: The Kalman gain weights {{c1::prediction uncertainty}} against {{c2::total uncertainty}}; small $R$ gives large $K$ (trust the sensor).
- CLOZE: For numerical stability use the {{c1::Joseph form}} of the covariance update to keep $P$ positive-definite. <br> <i>hint:</i> $P^+ = (I-KH)P^-(I-KH)^\top + KRK^\top$.
- CLOZE: Too-small $R$ makes the filter {{c1::chase sensor noise}}; too-small $Q$ makes it {{c2::ignore real motion and diverge}}.
- CLOZE: In code the innovation is {{c1::`y = z - H * x`}} — measurement minus predicted measurement.
