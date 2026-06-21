# Supervised learning as optimization

Prereqs: [[19_nonlinear-least-squares]], [[12_calculus-jacobians]] => this => unlocks [[28_transformers-attention]], [[29_gnn-message-passing]], [[30_pointnet-neural-implicits]].

## Why it matters for the stone scan->CAD->IFC / CADFit workflow

In your pipeline the deterministic stages (camera calibration, ICP registration, NURBS fitting, IK) are already nonlinear-least-squares problems: pick parameters, define residuals, minimize. Supervised learning is the *same loop* with a different, much larger parameterization of the function. You will reach for it exactly where hand-written geometry stalls: classifying which scan points are stone vs. saw-kerf vs. background, predicting per-point normals on noisy charnockite, segmenting a boulder into facets, or proposing a coarse pose before ICP refines it. CADFit (arXiv:2605.01171) is the canonical example for your end goal: it treats mesh-to-CAD as an optimization over a structured CAD program, scored by a volumetric IoU loss, which is supervised learning fused with a geometric kernel. The lesson to internalize now: a neural network is not a new kind of math, it is `min_theta (1/N) sum loss(f_theta(x_i), y_i) + lambda R(theta)` with a flexible `f_theta`, trained by the gradient descent you already use for calibration and IK.

## Intuition

You already fit a circle to scan points by nudging center and radius until the residual stops shrinking. Supervised learning does the same nudging, but the "shape" being fit is a function from inputs to outputs, and it has thousands or millions of knobs instead of three.

Analogy: think of `f_theta` as a mixing board with `theta` as the sliders. The training set is a stack of (input, correct output) recordings. The loss is a meter showing how far the board's output is from the correct recording, averaged over the stack. Gradient descent reads the meter, computes which direction to push every slider to lower it, and pushes a small step. Backpropagation is just the efficient bookkeeping that tells you the meter's sensitivity to every slider at once, by applying the chain rule from the output backward through the board. You stop when a *held-out* stack of recordings (the validation set) stops improving, because past that point you are memorizing the recordings instead of learning the mix.

## The math (first principles)

**Setup and conventions.** Training data is `N` pairs `{(x_i, y_i)}`, with inputs `x_i` (a scan patch, a point neighborhood, a feature vector) and targets `y_i` (a label, a normal, a class). The model is a parameterized function `f_theta : X -> Y` with parameter vector `theta in R^p`. We minimize a scalar objective. All gradients are column vectors of the same shape as `theta`. We use the convention that `nabla_theta` denotes the gradient with respect to all parameters.

**Empirical risk.** The true goal is low *expected* loss over the data distribution, the *risk*. We cannot see the distribution, only samples, so we minimize the *empirical risk* plus a regularizer:

$$ J(\theta) = \frac{1}{N}\sum_{i=1}^{N} \ell\big(f_\theta(x_i), y_i\big) + \lambda\, R(\theta). $$

This is the whole topic in one line. The first term fits the data, the second term `lambda R(theta)` penalizes complexity, and `lambda >= 0` trades them off. Verified against the ERM principle in statistical learning theory.

**Loss choices.** The loss `ell` measures per-example error. Two carry your whole workflow:

- Regression (predict a continuous target like a normal coordinate): squared error `ell(\hat y, y) = \tfrac12 \lVert \hat y - y\rVert^2`. This is exactly the nonlinear-least-squares residual from [[19_nonlinear-least-squares]].
- Classification (predict a discrete label like stone/not-stone): softmax cross-entropy. Given logits `z = f_theta(x) in R^K`, softmax produces probabilities

$$ p_k = \frac{e^{z_k}}{\sum_{j=1}^{K} e^{z_j}}, \qquad \ell(z, y) = -\sum_{k=1}^{K} y_k \log p_k, $$

where `y` is a one-hot target. The convention `0 log 0 = 0` applies.

**Regularizer.** The standard choice is L2 (ridge / weight decay): `R(\theta) = \tfrac12 \lVert \theta \rVert^2`. It shrinks parameters toward zero, which biases the model toward smoother functions and improves generalization. L1, `R(\theta)=\lVert\theta\rVert_1`, induces sparsity instead.

**Gradient descent.** Update parameters opposite the gradient with step size (learning rate) `eta > 0`:

$$ \theta \leftarrow \theta - \eta\, \nabla_\theta J(\theta). $$

With L2 regularization the gradient splits cleanly:

$$ \nabla_\theta J = \frac{1}{N}\sum_{i=1}^{N} \nabla_\theta\, \ell\big(f_\theta(x_i), y_i\big) + \lambda\, \theta. $$

The `+ lambda theta` term is why L2 is called *weight decay*: each step also multiplies `theta` by roughly `(1 - eta lambda)`. *Stochastic* gradient descent (SGD) replaces the full sum by an average over a small random *minibatch* `B`, giving a cheap unbiased estimate of the gradient. This is the only version that scales.

**Why this is the same optimizer as calibration/IK/ICP.** For the squared loss with a linear-in-residual model, `nabla J = J_r^T r` where `r` is the residual vector and `J_r` is its Jacobian. That is the Gauss-Newton / Levenberg-Marquardt right-hand side from [[19_nonlinear-least-squares]]. The damped-least-squares IK update `Delta theta = J^T (J J^T + lambda^2 I)^{-1} e` is the same residual-Jacobian machinery with Tikhonov damping. Supervised learning differs only in that `f_theta` is a deep composition, so we use plain (stochastic) first-order descent instead of forming `J^T J`, because `p` is too large to invert.

**Backpropagation = chain rule at scale.** A network is a composition `f_\theta(x) = g_L \circ g_{L-1} \circ \cdots \circ g_1(x)`, layer `l` mapping activation `a^{(l-1)}` to `a^{(l)} = g_l(a^{(l-1)}; \theta_l)`. Let the scalar loss be `ell`. Define the *adjoint* (the gradient of the loss w.r.t. a layer's output) `\delta^{(l)} = \partial \ell / \partial a^{(l)}`. The chain rule gives a backward recursion:

$$ \delta^{(L)} = \frac{\partial \ell}{\partial a^{(L)}}, \qquad \delta^{(l)} = \Big(\frac{\partial a^{(l+1)}}{\partial a^{(l)}}\Big)^{\!\top} \delta^{(l+1)}, \qquad \nabla_{\theta_l}\ell = \Big(\frac{\partial a^{(l)}}{\partial \theta_l}\Big)^{\!\top}\delta^{(l)}. $$

`\partial a^{(l+1)} / \partial a^{(l)}` is the layer Jacobian. Backpropagation is reverse-mode automatic differentiation: one forward pass caches activations, one backward pass propagates `delta` from output to input, and every parameter gradient costs about the same as one forward pass total. This is what makes training networks with millions of parameters feasible.

**The softmax cross-entropy gradient (a result worth memorizing).** Differentiating `ell` w.r.t. the *logits* `z`, the softmax and log compose to a clean form:

$$ \frac{\partial \ell}{\partial z_k} = p_k - y_k. $$

So `nabla_z ell = p - y`, predicted minus true. This is why softmax + cross-entropy is the default classification head: the output-layer adjoint is just the prediction error. Verified by the standard derivation.

**Train / validation / test split and overfitting.** Partition data into three disjoint sets. *Train* fits `theta`. *Validation* selects hyperparameters (`lambda`, `eta`, architecture, when to stop). *Test* is touched once, to estimate generalization. Overfitting is the gap: train loss keeps falling while validation loss rises, because the model is fitting noise specific to the training sample. The bias-variance tradeoff frames it: too-simple `f_theta` has high bias (underfits), too-flexible `f_theta` has high variance (overfits). Regularization, more data, and *early stopping* (halt when validation loss stops improving) all reduce variance.

## Worked example

A two-class logistic regression by hand. Model: `f_theta(x) = sigma(w x + b)` with `sigma(z) = 1/(1+e^{-z})`, parameters `theta = (w, b)`. Loss: binary cross-entropy `ell = -[y log p + (1-y) log(1-p)]` where `p = sigma(z)`, `z = wx + b`. No regularization (`lambda = 0`). Learning rate `eta = 1`.

Data, three points: `(x, y) = (1, 1), (2, 0), (3, 1)`. Initialize `w = 0, b = 0`.

**Forward pass.** With `w = b = 0`, every `z = 0`, so `p = sigma(0) = 0.5` for all three points.

**Per-example gradient.** For logistic regression the logit gradient is the same clean form as softmax: `d ell / d z = p - y`. Then by the chain rule `d ell / d w = (p - y) x` and `d ell / d b = (p - y)`.

| i | x | y | p   | p - y | (p-y)x |
|---|---|---|-----|-------|--------|
| 1 | 1 | 1 | 0.5 | -0.5  | -0.5   |
| 2 | 2 | 0 | 0.5 | +0.5  | +1.0   |
| 3 | 3 | 1 | 0.5 | -0.5  | -1.5   |

**Batch gradient (average over N = 3).**

$$ \frac{\partial J}{\partial w} = \frac{1}{3}(-0.5 + 1.0 - 1.5) = \frac{-1.0}{3} = -0.3333, \qquad \frac{\partial J}{\partial b} = \frac{1}{3}(-0.5 + 0.5 - 0.5) = -0.1667. $$

**One gradient-descent step** (`theta <- theta - eta nabla J`):

$$ w \leftarrow 0 - 1\cdot(-0.3333) = 0.3333, \qquad b \leftarrow 0 - 1\cdot(-0.1667) = 0.1667. $$

**Check it moved the right way.** Recompute `p` for point 1: `z = 0.3333(1) + 0.1667 = 0.5`, `p = sigma(0.5) = 0.6225`. Was `0.5`, target is `1`, so `p` rose toward 1. Point 2: `z = 0.3333(2) + 0.1667 = 0.8333`, `p = 0.6970`. That *rose* away from its target `0`, because point 3 (`x=3, y=1`) pulls `w` positive and dominates. This is the genuine tension a single linear boundary cannot resolve perfectly: the arithmetic is honest about the model's limited capacity. Initial average loss is `-log(0.5) = 0.6931`; after the step the average loss is `(−log 0.6225 − log(1−0.6970) − log 0.7503)/3 = (0.4741 + 1.1939 + 0.2874)/3 = 0.6518`, a decrease, confirming the step lowered `J`.

## Code

Idiomatic C++ with Eigen first: a full-batch L2-regularized logistic-regression trainer. The math symbols map directly to variable names in the comments. This is the smallest program that exercises the whole loop: `f_theta`, loss, gradient, descent step, train/validation monitoring.

```cpp
// logreg.cpp  -- supervised learning as optimization, minimal & runnable.
// Build: g++ -std=c++17 -I /path/to/eigen logreg.cpp -o logreg
#include <Eigen/Dense>
#include <iostream>
#include <cmath>

using Eigen::MatrixXd;
using Eigen::VectorXd;

// sigma(z) = 1/(1+e^{-z}), elementwise. f_theta(x) = sigma(X w + b).
static VectorXd sigmoid(const VectorXd& z) {
    return (1.0 / (1.0 + (-z.array()).exp())).matrix();
}

// Mean binary cross-entropy ell + (lambda/2)||w||^2 ; the objective J(theta).
static double objective(const MatrixXd& X, const VectorXd& y,
                        const VectorXd& w, double b, double lambda) {
    const VectorXd p = sigmoid((X * w).array() + b);          // predictions p_i
    const double eps = 1e-12;                                 // guard log(0)
    const Eigen::ArrayXd pa = p.array().max(eps).min(1.0 - eps);
    const double ce = -(y.array() * pa.log()
                        + (1.0 - y.array()) * (1.0 - pa).log()).mean();
    return ce + 0.5 * lambda * w.squaredNorm();               // + (lambda/2)||w||^2
}

int main() {
    // Toy 2D set: two blobs, label by which side of x0 = 0.
    const int N = 8;
    MatrixXd X(N, 2);
    VectorXd y(N);
    X <<  1.0,  0.5,   1.5, -0.2,   0.8,  1.0,   2.0,  0.3,    // class 1 (x0 > 0)
         -1.0,  0.4,  -1.6, -0.3,  -0.7,  0.9,  -2.1,  0.1;    // class 0 (x0 < 0)
    y <<  1, 1, 1, 1,  0, 0, 0, 0;

    VectorXd w = VectorXd::Zero(2);   // theta = (w, b)
    double   b = 0.0;
    const double eta = 0.5;           // learning rate
    const double lambda = 1e-3;       // L2 strength  R(theta) = (1/2)||w||^2

    for (int epoch = 0; epoch <= 200; ++epoch) {
        // Forward: p = f_theta(X).
        const VectorXd p = sigmoid((X * w).array() + b);
        // Output adjoint of cross-entropy + sigmoid is (p - y): grad w.r.t. logits.
        const VectorXd err = p - y;                           // p - y
        // Batch gradients: dJ/dw = (1/N) X^T (p - y) + lambda w ; dJ/db = mean(p - y).
        const VectorXd grad_w = (X.transpose() * err) / N + lambda * w;
        const double   grad_b = err.mean();
        // Gradient-descent step: theta <- theta - eta * grad.
        w -= eta * grad_w;
        b -= eta * grad_b;
        if (epoch % 50 == 0)
            std::cout << "epoch " << epoch
                      << "  J = " << objective(X, y, w, b, lambda)
                      << "  w = [" << w.transpose() << "]  b = " << b << "\n";
    }
    std::cout << "Learned boundary: " << w(0) << "*x0 + " << w(1)
              << "*x1 + " << b << " = 0\n";
    return 0;
}
```

Expected behavior: `J` falls monotonically, `w(0)` grows positive (the discriminating axis), `w(1)` stays near zero, and the printed boundary approximates `x0 = 0`. The line `err = p - y` is the load-bearing identity from the math section.

A short Python version clarifies the train/validation split and early stopping, which the C++ snippet omits for brevity:

```python
import numpy as np

def sigmoid(z): return 1.0 / (1.0 + np.exp(-z))

# X: (N,2) features, y: (N,) labels. Split 70/30 into train / validation.
rng = np.random.default_rng(0)
X = np.vstack([rng.normal([1, 0], 0.5, (40, 2)),
               rng.normal([-1, 0], 0.5, (40, 2))])
y = np.r_[np.ones(40), np.zeros(40)]
idx = rng.permutation(80); tr, va = idx[:56], idx[56:]

w, b, eta, lam = np.zeros(2), 0.0, 0.5, 1e-3
best_val, best, patience, wait = np.inf, None, 10, 0
for epoch in range(500):
    p = sigmoid(X[tr] @ w + b)
    err = p - y[tr]                          # p - y : output adjoint
    w -= eta * (X[tr].T @ err / len(tr) + lam * w)   # weight decay
    b -= eta * err.mean()
    pv = np.clip(sigmoid(X[va] @ w + b), 1e-12, 1 - 1e-12)
    val = -(y[va]*np.log(pv) + (1-y[va])*np.log(1-pv)).mean()  # validation loss
    if val < best_val - 1e-5:                # early stopping on validation loss
        best_val, best, wait = val, (w.copy(), b), 0
    else:
        wait += 1
        if wait >= patience: break
print("stopped epoch", epoch, "val loss", round(best_val, 4))
```

## Connections

- [[19_nonlinear-least-squares]] -- supervised learning *is* this with a deep, high-dimensional `f_theta`; squared-error regression reduces to it exactly, and the residual-Jacobian intuition transfers directly.
- [[12_calculus-jacobians]] -- backpropagation is the chain rule applied to layer Jacobians; every `delta^{(l)}` is a Jacobian-transpose times a vector.
- [[28_transformers-attention]] -- attention is one more differentiable `f_theta` block, `softmax(QK^T/sqrt(d_k))V`, trained by the same objective and backprop.
- [[29_gnn-message-passing]] -- a GNN parameterizes `f_theta` over graph neighborhoods (meshes, point neighborhoods); the optimization loop is unchanged.
- [[30_pointnet-neural-implicits]] -- PointNet and DeepSDF/Occupancy Networks swap the parameterization for permutation-invariant point sets and implicit fields, but minimize the same regularized empirical risk; this is the direct on-ramp to learned scan-to-CAD and CADFit.

## Pitfalls and failure modes

- **Train/test leakage.** If any validation or test example influenced training (shared scans, normalized using global statistics, near-duplicate point patches), reported accuracy is fiction. Split *before* preprocessing; fit scalers on train only.
- **Overfitting from too little data.** Your domain is data-starved (featureless convex stone blocks overfit instantly). Watch the train-vs-validation gap, not train loss. Add `lambda`, augment, or reduce capacity.
- **Learning rate ill-conditioning.** `eta` too large diverges (loss to NaN); too small crawls. Unscaled inputs make the loss surface elongated and badly conditioned, the same conditioning problem as in least squares. Standardize features (zero mean, unit variance) before training.
- **Numeric overflow in softmax / log.** `exp(z)` overflows for large logits; `log(0)` is `-inf`. Use the log-sum-exp trick (subtract `max(z)` before exponentiating) and clamp probabilities away from 0 and 1, as the code does with `eps`.
- **Frame and unit inconsistency in geometric targets.** If `y_i` is a normal or pose, mixing frames or handedness, or training in millimeters and inferring in meters, silently corrupts every gradient. Pin one frame, one unit, one handedness across the whole dataset (the same rule that governs ICP and IK here).
- **Wrong loss for the task.** Squared error on a classification problem, or cross-entropy without softmax, gives a malformed gradient. Match the loss to the output type; pair softmax with cross-entropy so the adjoint stays `p - y`.
- **Outliers dominate squared loss.** A few bad scan points pull a least-squares fit hard. Use a robust loss (Huber) or outlier rejection, as in robust registration.
- **Non-convexity and local minima.** Deep `f_theta` objectives are non-convex; SGD finds a good basin, not the global optimum. Initialization, seeds, and learning-rate schedule matter, and results are not bitwise reproducible without fixing all of them.

## Practice

1. **Hand-trace one step (warm-up).** Redo the worked example with `eta = 0.5` instead of `1.0`. Compute the new `w`, `b`, and the average loss after one step. Confirm the loss still decreases. Inspectable: a table of `p` before and after.
2. **Run and read the C++ trainer.** Compile `logreg.cpp`. Print `J` every epoch and plot it (or eyeball the printed values). Then set `lambda = 1.0` and observe `||w||` shrink and the boundary blur. Inspectable: the loss curve and final boundary coefficients.
3. **Verify the gradient numerically.** For the logistic model, implement finite differences `(J(theta + e) - J(theta - e)) / (2e)` with `e = 1e-5` and compare elementwise to the analytic `nabla J`. They should agree to about 5 decimals. This is the single most valuable debugging habit for any learned `f_theta`.
4. **Overfit on purpose.** Make a tiny dataset (5 points) and a high-capacity model (add polynomial features up to degree 8). Train to near-zero train loss, then evaluate on 50 held-out points. Plot train vs. validation loss against epoch and locate the early-stopping point. Inspectable: the diverging curves.
5. **Geometric target (stretch).** Replace the classification target with a regression target: predict a 2D unit normal from a local point feature, using squared loss with the network normalized to unit length. Confirm the gradient reduces to the least-squares form `J_r^T r`. Inspectable: predicted vs. true normals as arrows.
6. **Softmax head (stretch).** Extend to `K = 3` classes with a softmax output and cross-entropy loss. Implement the `p - y` logit gradient and verify it against finite differences. Inspectable: the `3 x p` gradient matrix matching numerics.

## References

- Goodfellow, Bengio, Courville, *Deep Learning* (2016), Ch. 5 (machine learning basics: ERM, capacity, regularization, bias-variance) and Ch. 6 (backpropagation as reverse-mode autodiff). https://www.deeplearningbook.org/
- Hardt and Recht, *Patterns, Predictions, and Actions* (mlstory.org), Optimization chapter (empirical risk minimization and gradient methods). https://mlstory.org/optimization.html
- Bishop, *Pattern Recognition and Machine Learning* (2006), Ch. 4-5 (logistic regression, cross-entropy, gradient training). https://www.microsoft.com/en-us/research/publication/pattern-recognition-machine-learning/
- "Derivative of the Softmax Function and the Categorical Cross-Entropy Loss" (derivation of the `p - y` logit gradient), GeeksforGeeks. https://www.geeksforgeeks.org/machine-learning/derivative-of-the-softmax-function-and-the-categorical-cross-entropy-loss/
- Levenberg-Marquardt / damped least squares (shared optimizer backbone with calibration and IK): see [[19_nonlinear-least-squares]] references.
- Vaswani et al., "Attention Is All You Need" (2017), for scaled dot-product attention as a differentiable `f_theta` block. https://arxiv.org/abs/1706.03762
- Kipf and Welling, "Semi-Supervised Classification with Graph Convolutional Networks" (2017). https://arxiv.org/abs/1609.02907
- Qi et al., "PointNet: Deep Learning on Point Sets for 3D Classification and Segmentation" (2017). https://arxiv.org/abs/1612.00593
- Park et al., "DeepSDF: Learning Continuous Signed Distance Functions for Shape Representation" (2019). https://arxiv.org/abs/1901.05103
- Mescheder et al., "Occupancy Networks: Learning 3D Reconstruction in Function Space" (2019). https://arxiv.org/abs/1812.03828
- Nehme, Whalen, Ahmed, "CADFit: Precise Mesh-to-CAD Program Generation with Hybrid Optimization" (arXiv:2605.01171); code at https://github.com/ghadinehme/CADFit . https://arxiv.org/abs/2605.01171
- Eigen documentation (array/matrix ops used in the C++ trainer). https://eigen.tuxfamily.org/

## Flashcard seeds

- Q: The objective is $J(\theta)$ and you are at a point where the gradient $\nabla_\theta J$ points uphill. Write the gradient-descent update that steps $\theta$ downhill. <br> <b>Formula:</b> $\theta \leftarrow \theta - \eta\,\nabla_\theta J(\theta)$ <br> <i>Hint:</i> downhill means the direction opposite the gradient; scale that direction by a small step size $\eta>0$. :: A: Steepest descent moves opposite the gradient: direction $-\nabla_\theta J$ <br> scale by step size $\eta>0$ <br> $\boxed{\theta \leftarrow \theta - \eta\,\nabla_\theta J(\theta)}$.
- Q: Find the gradient of the regularized risk $J(\theta)=\frac1N\sum_i\ell(f_\theta(x_i),y_i)+\lambda R(\theta)$ with $R(\theta)=\tfrac12\lVert\theta\rVert^2$ (derive, step 1 of 2). <br> <b>Formula:</b> target $\nabla_\theta J=\frac1N\sum_i\nabla_\theta\ell+\lambda\theta$; use $\nabla_\theta\big(\tfrac12\lVert\theta\rVert^2\big)=\theta$ <br> <i>Hint:</i> the gradient is linear, so differentiate the data term and the penalty term separately. :: A: Differentiate term by term <br> data term: $\nabla_\theta\big[\frac1N\sum_i\ell\big]=\frac1N\sum_i\nabla_\theta\ell(f_\theta(x_i),y_i)$ <br> penalty term: $\nabla_\theta\big[\tfrac{\lambda}{2}\lVert\theta\rVert^2\big]=\lambda\theta$ <br> $\boxed{\nabla_\theta J=\frac1N\sum_i\nabla_\theta\ell+\lambda\theta}$.
- Q: From $\nabla_\theta J=\frac1N\sum_i\nabla_\theta\ell+\lambda\theta$, show that L2 regularization makes each GD step shrink $\theta$ — "weight decay" (derive, step 2 of 2). <br> <b>Formula:</b> $\theta\leftarrow\theta-\eta\nabla_\theta J$; target $\theta\leftarrow(1-\eta\lambda)\theta-\tfrac{\eta}{N}\sum_i\nabla_\theta\ell$ <br> <i>Hint:</i> substitute the gradient into the update, then collect the two $\theta$ terms. :: A: Substitute into $\theta\leftarrow\theta-\eta\nabla_\theta J$ <br> $\theta\leftarrow\theta-\eta\big(\frac1N\sum_i\nabla_\theta\ell\big)-\eta\lambda\theta$ <br> group the $\theta$ terms: $\theta\leftarrow(1-\eta\lambda)\theta-\frac{\eta}{N}\sum_i\nabla_\theta\ell$ <br> $\boxed{\theta\leftarrow(1-\eta\lambda)\theta-\tfrac{\eta}{N}\textstyle\sum_i\nabla_\theta\ell}$ — each step shrinks $\theta$ by factor $(1-\eta\lambda)$.
- Q: Given logits $z\in\mathbb{R}^K$, build probabilities $p_k$ that are all positive and sum to 1 (derive softmax, step 1 of 3). <br> <b>Formula:</b> target $p_k=\dfrac{e^{z_k}}{\sum_{j=1}^{K}e^{z_j}}$ <br> <i>Hint:</i> first apply a function that makes everything positive, then divide by the total so the entries sum to 1. :: A: Exponentiate to force positivity: $e^{z_k}>0$ <br> normalize so the vector sums to 1: divide by $\sum_j e^{z_j}$ <br> $\boxed{p_k=\dfrac{e^{z_k}}{\sum_{j=1}^{K}e^{z_j}}}$.
- Q: From $p_k=e^{z_k}/\sum_j e^{z_j}$, find the softmax Jacobian $\partial p_k/\partial z_m$ (derive, step 2 of 3). <br> <b>Formula:</b> quotient rule $\frac{d}{dz}\frac{u}{S}=\frac{u'S-uS'}{S^2}$; target $\partial p_k/\partial z_m=p_k(\delta_{km}-p_m)$ <br> <i>Hint:</i> let $S=\sum_j e^{z_j}$ and split into the cases $m=k$ and $m\neq k$. :: A: Quotient rule with $S=\sum_j e^{z_j}$ <br> if $m=k$: $\partial p_k/\partial z_k=\frac{e^{z_k}S-e^{z_k}e^{z_k}}{S^2}=p_k(1-p_k)$ <br> if $m\neq k$: $\partial p_k/\partial z_m=\frac{-e^{z_k}e^{z_m}}{S^2}=-p_kp_m$ <br> $\boxed{\partial p_k/\partial z_m=p_k(\delta_{km}-p_m)}$.
- Q: From $\ell=-\sum_k y_k\log p_k$ and $\partial p_k/\partial z_m=p_k(\delta_{km}-p_m)$ with one-hot $y$, derive the cross-entropy logit gradient $\partial\ell/\partial z_m$ (step 3 of 3). <br> <b>Formula:</b> chain rule $\frac{\partial\ell}{\partial z_m}=\sum_k\frac{\partial\ell}{\partial p_k}\frac{\partial p_k}{\partial z_m}$; target $\partial\ell/\partial z_m=p_m-y_m$ <br> <i>Hint:</i> use $\partial\ell/\partial p_k=-y_k/p_k$, cancel the $p_k$, then use $\sum_k y_k=1$. :: A: $\frac{\partial\ell}{\partial z_m}=-\sum_k\frac{y_k}{p_k}\,p_k(\delta_{km}-p_m)=-\sum_k y_k(\delta_{km}-p_m)$ <br> $=-y_m+p_m\sum_k y_k$, and $\sum_k y_k=1$ <br> $\boxed{\partial\ell/\partial z_m=p_m-y_m}$ (predicted minus true).
- Q: For a network $a^{(l+1)}=g_{l+1}(a^{(l)})$ with adjoint $\delta^{(l)}=\partial\ell/\partial a^{(l)}$, derive the backward recursion relating $\delta^{(l)}$ to $\delta^{(l+1)}$ (step 1 of 2). <br> <b>Formula:</b> chain rule; target $\delta^{(l)}=\big(\partial a^{(l+1)}/\partial a^{(l)}\big)^{\top}\delta^{(l+1)}$ <br> <i>Hint:</i> route $\partial\ell/\partial a^{(l)}$ through $a^{(l+1)}$, then recognize $\partial\ell/\partial a^{(l+1)}$ as $\delta^{(l+1)}$. :: A: Chain rule through the next layer's input <br> $\dfrac{\partial\ell}{\partial a^{(l)}}=\Big(\dfrac{\partial a^{(l+1)}}{\partial a^{(l)}}\Big)^{\!\top}\dfrac{\partial\ell}{\partial a^{(l+1)}}$ <br> recognize the right factor as $\delta^{(l+1)}$ <br> $\boxed{\delta^{(l)}=\big(\partial a^{(l+1)}/\partial a^{(l)}\big)^{\top}\delta^{(l+1)}}$.
- Q: With layer params entering via $a^{(l)}=g_l(a^{(l-1)};\theta_l)$ and cached adjoint $\delta^{(l)}$, derive the parameter gradient $\nabla_{\theta_l}\ell$ (step 2 of 2). <br> <b>Formula:</b> chain rule through $a^{(l)}$; target $\nabla_{\theta_l}\ell=\big(\partial a^{(l)}/\partial\theta_l\big)^{\top}\delta^{(l)}$ <br> <i>Hint:</i> the only path from $\theta_l$ to the loss runs through $a^{(l)}$; multiply the layer Jacobian transpose by $\delta^{(l)}$. :: A: Chain rule from loss to $\theta_l$ routes through $a^{(l)}$ <br> $\nabla_{\theta_l}\ell=\Big(\dfrac{\partial a^{(l)}}{\partial\theta_l}\Big)^{\!\top}\dfrac{\partial\ell}{\partial a^{(l)}}$ <br> $\boxed{\nabla_{\theta_l}\ell=\big(\partial a^{(l)}/\partial\theta_l\big)^{\top}\delta^{(l)}}$ — a layer-Jacobian-transpose times the cached adjoint.
- Q: For squared loss $\ell=\tfrac12\lVert r(\theta)\rVert^2$ with residual $r$ and Jacobian $J_r=\partial r/\partial\theta$, derive $\nabla_\theta\ell$ and show it equals the Gauss-Newton right-hand side. <br> <b>Formula:</b> chain rule $\nabla_\theta\ell=(\partial r/\partial\theta)^{\top}\,\partial\ell/\partial r$; target $\nabla_\theta\ell=J_r^{\top}r$ <br> <i>Hint:</i> first compute the inner derivative $\partial\ell/\partial r$ for $\ell=\tfrac12\lVert r\rVert^2$. :: A: Chain rule: $\nabla_\theta\ell=\big(\partial r/\partial\theta\big)^{\top}\,\partial\ell/\partial r$ <br> $\partial\ell/\partial r=r$ <br> $\boxed{\nabla_\theta\ell=J_r^{\top}r}$ — the Gauss-Newton right-hand side, identical to calibration/ICP/IK.
- Q: Regression at mm scale: model output $\hat y=(0.30,0.40)$, target unit normal $y=(0,1)$. Compute the squared loss $\ell$ and the output gradient $\partial\ell/\partial\hat y$. <br> <b>Formula:</b> $\ell=\tfrac12\lVert\hat y-y\rVert^2$, $\ \partial\ell/\partial\hat y=\hat y-y$ <br> <i>Hint:</i> compute the residual $\hat y-y$ first; the gradient is exactly that residual. :: A: $\hat y-y=(0.30,-0.60)$ <br> $\ell=\tfrac12(0.30^2+0.60^2)=\tfrac12(0.09+0.36)=0.225$ <br> $\partial\ell/\partial\hat y=\hat y-y=\boxed{(0.30,-0.60)}$.
- Q: 3-class softmax (classify a scan point stone/kerf/background): logits $z=(2,1,0)$, true class index 1 so $y=(0,1,0)$. Compute $p$ and the logit gradient $p-y$. <br> <b>Formula:</b> $p_k=e^{z_k}/\sum_j e^{z_j}$, $\ \partial\ell/\partial z=p-y$ (use $e^2{=}7.389,\,e^1{=}2.718,\,e^0{=}1$) <br> <i>Hint:</i> exponentiate the logits and divide by their sum to get $p$, then subtract the one-hot $y$. :: A: $e^z=(7.389,2.718,1)$, sum $=11.107$ <br> $p=(0.665,0.245,0.090)$ <br> $\ell=-\log p_1=-\log 0.245=1.406$ <br> $p-y=\boxed{(0.665,-0.755,0.090)}$.
- Q: Logistic regression, one GD step: data $(x,y)=(1,1),(2,0),(3,1)$, start $w=b=0$, $\eta=1$. Compute the batch gradients $\partial J/\partial w,\ \partial J/\partial b$. <br> <b>Formula:</b> $p=\sigma(wx+b)$, $\ \partial\ell/\partial z=p-y$, $\ \partial\ell/\partial w=(p-y)x$, $\ \partial\ell/\partial b=(p-y)$; average over $N=3$ <br> <i>Hint:</i> at $w=b=0$ every $z=0$ so $p=\sigma(0)=0.5$; build the $(p-y)$ and $(p-y)x$ columns, then average. :: A: At $w=b=0$, every $z=0\Rightarrow p=\sigma(0)=0.5$ <br> $(p-y)$: $-0.5,\,+0.5,\,-0.5$; $(p-y)x$: $-0.5,\,+1.0,\,-1.5$ <br> $\partial J/\partial w=\tfrac13(-0.5+1.0-1.5)=\boxed{-0.3333}$, $\partial J/\partial b=\tfrac13(-0.5)=-0.1667$.
- Q: Apply the GD step: gradients $\partial J/\partial w=-0.3333,\ \partial J/\partial b=-0.1667$, learning rate $\eta=1$, start $w=b=0$. Compute updated $w,b$. <br> <b>Formula:</b> $\theta\leftarrow\theta-\eta\,\partial J/\partial\theta$ <br> <i>Hint:</i> subtract $\eta$ times each gradient; subtracting a negative gradient increases the parameter. :: A: $w\leftarrow 0-1\cdot(-0.3333)=0.3333$ <br> $b\leftarrow 0-1\cdot(-0.1667)=0.1667$ <br> $\boxed{(w,b)=(0.3333,\,0.1667)}$.
- Q: Large-scale weight decay: with $\eta=0.1,\ \lambda=0.01$, by what fraction does a weight shrink per step, and after 100 steps (ignoring the data gradient)? <br> <b>Formula:</b> per-step factor $(1-\eta\lambda)$; over $n$ steps $(1-\eta\lambda)^n$ <br> <i>Hint:</i> compute the single-step factor $1-\eta\lambda$ first, then raise it to the 100th power. :: A: per step factor $=1-0.1\cdot0.01=0.999$ <br> $\lVert\theta\rVert$ keeps $0.999$ of its value each step <br> per 100 steps: $0.999^{100}\approx 0.905$, $\boxed{\approx 9.5\%\text{ decay}}$.
- Q: Closing test (the step moved $p$ the right way): after the logistic step $(w,b)=(0.3333,0.1667)$, recompute $p_1$ for point 1 ($x=1,y=1$) and check it moved toward the target. <br> <b>Formula:</b> $z=wx+b$, $\ p=\sigma(z)=1/(1+e^{-z})$ <br> <i>Hint:</i> a pass means the new $p_1$ is larger than the old $0.5$ (closer to target $1$). :: A: $z_1=0.3333(1)+0.1667=0.5$ <br> $p_1=\sigma(0.5)=1/(1+e^{-0.5})=0.6225$ <br> was $0.5$, target $1$: $p_1$ rose $\Rightarrow$ verified moved correctly <br> $\boxed{p_1=0.6225\uparrow}$.
- Q: Closing test (loss decreased): stepped logistic model gives $p_1=0.6225,\ p_2=0.6971,\ p_3=0.7625$ for targets $1,0,1$. Compute the mean loss and check it fell below the initial $0.6931$. <br> <b>Formula:</b> $\ell_i=-[y_i\log p_i+(1-y_i)\log(1-p_i)]$, then average over $N=3$ <br> <i>Hint:</i> a pass means the mean is $<0.6931$; for $y=0$ use the $-\log(1-p)$ branch. :: A: terms $-\log0.6225=0.474$, $-\log(1-0.6971)=1.194$, $-\log0.7625=0.271$ <br> mean $=(0.474+1.194+0.271)/3=0.6465$ <br> $0.6465<0.6931$ $\Rightarrow$ $\boxed{J\text{ decreased, step valid}}$.
- Q: Closing test (stability/limit check): in the weight-decay update $\theta\leftarrow(1-\eta\lambda)\theta$, what range of $\eta\lambda$ keeps it a stable shrink, and what happens as $\lambda\to0$? <br> <b>Formula:</b> stable shrink requires $|1-\eta\lambda|<1$ <br> <i>Hint:</i> solve $|1-\eta\lambda|<1$ for $\eta\lambda$; then set $\lambda=0$ in the factor. :: A: need $0\le\eta\lambda<1$ so $|1-\eta\lambda|<1$ (shrinks, not flips/grows) <br> limit $\lambda\to0$: factor $\to1\Rightarrow$ no decay, plain GD recovered <br> $\boxed{0\le\eta\lambda<1,\ \ \lambda\to0\Rightarrow\text{plain GD}}$.
- CLOZE: Stochastic gradient descent replaces the full sum by an average over a small random {{c1::minibatch}}, an unbiased estimate of the gradient.
- CLOZE: The softmax+cross-entropy logit gradient is the clean form {{c1::$p-y$}}, predicted minus true, which is why they are paired.
- CLOZE: Each backprop step is a Jacobian-transpose times a vector: $\delta^{(l)}=(\partial a^{(l+1)}/\partial a^{(l)})^{\top}\,{{c1::\delta^{(l+1)}}}$.
- CLOZE: The {{c1::bias-variance}} tradeoff: too-simple $f_\theta$ underfits (high bias), too-flexible overfits (high variance).
- CLOZE: To avoid softmax overflow, subtract {{c1::$\max(z)$}} from the logits before exponentiating (log-sum-exp trick).
- Q: When-to-use: which loss for predicting a continuous geometric target (a per-point normal) versus a discrete label (stone/not-stone)? <br> <b>Formula:</b> regression $\ell=\tfrac12\lVert\hat y-y\rVert^2$; classification softmax + cross-entropy with adjoint $p-y$ <br> <i>Hint:</i> match the loss to the output type — continuous output to squared error, discrete output to cross-entropy. :: A: Continuous: squared error $\tfrac12\lVert\hat y-y\rVert^2$ (the NLS residual) <br> Discrete: softmax + cross-entropy, so the logit adjoint stays $p-y$ <br> $\boxed{\text{regression}\to\text{MSE},\ \text{classification}\to\text{softmax-CE}}$.
- Q: Definition: distinguish the train, validation, and test sets by their job. <br> <i>Hint:</i> three disjoint roles — fit parameters, choose hyperparameters, estimate generalization. :: A: Train fits $\theta$ <br> validation selects hyperparameters ($\lambda,\eta$, when to stop) <br> test is touched once to estimate generalization <br> $\boxed{\text{fit}\,/\,\text{select}\,/\,\text{estimate}}$.
- Q: Pitfall: name the train/test leakage that makes reported accuracy fiction, and the fix. <br> <i>Hint:</i> any way a validation/test example influenced training counts; the fix is about ordering preprocessing relative to the split. :: A: Leakage = validation/test influenced training (shared scans, global-stat normalization, near-duplicate patches) <br> Fix: split before preprocessing; fit scalers on train only <br> $\boxed{\text{split first, fit scalers on train only}}$.
- Q: Pitfall: why use plain first-order SGD for deep nets instead of Gauss-Newton? <br> <b>Formula:</b> Gauss-Newton needs $J^{\top}J\in\mathbb{R}^{p\times p}$ formed and inverted <br> <i>Hint:</i> think about the size of $p$ (the parameter count) for a deep network. :: A: Gauss-Newton forms and inverts $J^{\top}J$, which is $p\times p$ <br> deep $f_\theta$ has $p$ in the millions, too large to store or invert <br> $\boxed{p\text{ too large}\Rightarrow\text{first-order descent}}$.
- Q: Closing test (gradient check): how do you verify any analytic gradient $\nabla J$ is correct? <br> <b>Formula:</b> central difference $\big(J(\theta+e)-J(\theta-e)\big)/(2e)$ with $e\approx10^{-5}$ <br> <i>Hint:</i> a pass means the finite-difference value matches the analytic gradient elementwise to about 5 decimals. :: A: Perturb each component by $e\approx10^{-5}$ and form $(J(\theta+e)-J(\theta-e))/(2e)$ <br> compare elementwise to analytic $\nabla J$ <br> agreement to ~5 decimals $\Rightarrow$ $\boxed{\text{gradient verified}}$.
