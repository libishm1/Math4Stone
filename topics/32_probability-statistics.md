# Probability and statistics

Placement: [[31_calculus-foundations]] => **Probability and statistics** => [[27_supervised-learning-optimization]], [[09_robust-ransac]]

## Why it matters for the stone scan->CAD->IFC / CADFit / robotic-fabrication workflow

Scan data is noisy: every LiDAR or photogrammetry point carries measurement error, and probability is how you model and tame it. The reason least squares is everywhere in this corpus is a probabilistic fact: minimizing squared residuals is the maximum-likelihood estimate under Gaussian noise. Covariance matrices describe the uncertainty of a fitted plane, a registration pose, or a normal estimate, and they feed RANSAC's inlier thresholds in [[09_robust-ransac]]. The machine-learning layer ([[27_supervised-learning-optimization]] and the neural-geometry topics) is statistics with many parameters. So this node sits under both the optimization story and the AI story.

## Intuition

Probability quantifies uncertainty before you see data; statistics infers the unknowns after you see data. A distribution is a budget of belief spread over possible outcomes. The Gaussian (bell curve) shows up so often because sums of many small independent errors tend toward it (the central limit theorem), which is exactly what sensor noise is.

Analogy: you measure a stone's edge ten times and get ten slightly different numbers. The mean is your best single guess; the spread (variance) is your honesty about how unsure you are. A good pipeline carries that spread forward instead of pretending each measurement is exact.

## The math (first principles)

**Random variable and distribution.** A random variable $X$ takes values with a probability law. Discrete variables have a probability mass function $p(x)$ with $\sum_x p(x)=1$. Continuous variables have a probability density $p(x)$ with $\int p(x)\,dx = 1$ (note the integral, hence the [[31_calculus-foundations]] prerequisite).

**Expectation, variance, standard deviation.**

$$ \mathbb{E}[X] = \mu = \int x\,p(x)\,dx, \qquad \mathrm{Var}(X) = \mathbb{E}[(X-\mu)^2] = \sigma^2, \qquad \sigma = \sqrt{\mathrm{Var}(X)}. $$

The mean is the center of mass of the distribution; the variance is its spread.

**Gaussian / normal.** The one-dimensional normal density is

$$ p(x) = \frac{1}{\sqrt{2\pi}\,\sigma}\exp\!\left(-\frac{(x-\mu)^2}{2\sigma^2}\right). $$

The **multivariate Gaussian** for $x \in \mathbb{R}^n$ with mean $\mu$ and covariance $\Sigma$ (symmetric positive definite) is

$$ p(x) = \frac{1}{(2\pi)^{n/2}\lvert\Sigma\rvert^{1/2}}\exp\!\left(-\tfrac12 (x-\mu)^\top \Sigma^{-1}(x-\mu)\right). $$

**Covariance matrix.** For a random vector, $\Sigma_{ij} = \mathbb{E}[(X_i-\mu_i)(X_j-\mu_j)]$. The diagonal holds per-axis variances; off-diagonals hold correlations. Its eigenvectors are the principal axes of the uncertainty ellipsoid, which is the link to PCA in [[03_eig-svd-pca]].

**Conditional probability and Bayes.**

$$ p(A\mid B) = \frac{p(A\cap B)}{p(B)}, \qquad p(\theta\mid D) = \frac{p(D\mid \theta)\,p(\theta)}{p(D)}. $$

Bayes turns a likelihood and a prior into a posterior belief about parameters $\theta$ given data $D$.

**Maximum likelihood (MLE) and MAP.** The likelihood is $\mathcal{L}(\theta) = p(D\mid\theta)$. MLE picks $\hat\theta = \arg\max_\theta \mathcal{L}(\theta)$; MAP adds a prior, $\arg\max_\theta p(D\mid\theta)p(\theta)$. Maximizing the log-likelihood is standard because logs turn products into sums.

**The least-squares bridge.** Suppose residuals $r_i = y_i - f(x_i;\theta)$ are independent zero-mean Gaussian with variance $\sigma^2$. The log-likelihood is

$$ \log\mathcal{L}(\theta) = -\frac{1}{2\sigma^2}\sum_i r_i^2 + \text{const}. $$

Maximizing it is the same as minimizing $\sum_i r_i^2$. **This is why least squares is the default**: it is the MLE under i.i.d. Gaussian noise. With per-point variances $\sigma_i^2$, the MLE becomes *weighted* least squares with weights $1/\sigma_i^2$. Connects to [[06_linear-least-squares]] and [[08_nonlinear-least-squares]].

**Law of large numbers and central limit theorem.** Sample averages converge to the true mean as $n\to\infty$ (LLN). The distribution of a sum of many independent finite-variance variables approaches a Gaussian (CLT), which justifies the Gaussian noise model in the first place.

## Worked example

Five edge measurements (mm): $100.2, 99.8, 100.1, 99.9, 100.0$. The sample mean is $\bar x = \frac{500.0}{5} = 100.0$. The deviations are $+0.2, -0.2, +0.1, -0.1, 0.0$; their squares sum to $0.04+0.04+0.01+0.01+0 = 0.10$. The unbiased sample variance is $s^2 = \frac{0.10}{5-1} = 0.025$, so $s \approx 0.158$ mm. Your best estimate of the edge is $100.0 \pm 0.16$ mm (one standard deviation). If you instead fit a line to noisy points and the noise is Gaussian, minimizing the squared vertical residuals gives the maximum-likelihood line.

## Code

Idiomatic C++ with Eigen: compute the sample mean and covariance of a set of 3-D points (the same covariance feeds PCA normal estimation), and evaluate a multivariate Gaussian density.

```cpp
#include <Eigen/Dense>
#include <vector>
#include <cmath>

using Vec3 = Eigen::Vector3d;
using Mat3 = Eigen::Matrix3d;

// Sample mean and unbiased (N-1) covariance of a point set.
void meanCovariance(const std::vector<Vec3>& pts, Vec3& mean, Mat3& cov) {
    const int n = static_cast<int>(pts.size());
    mean.setZero();
    for (const auto& p : pts) mean += p;
    mean /= n;                                   // mu = (1/n) sum p_i

    cov.setZero();
    for (const auto& p : pts) {
        Vec3 d = p - mean;
        cov += d * d.transpose();                // (p - mu)(p - mu)^T
    }
    cov /= (n - 1);                              // unbiased estimator divides by n-1
}

// Multivariate Gaussian density N(x; mu, Sigma).
double gaussianPdf(const Vec3& x, const Vec3& mu, const Mat3& Sigma) {
    Mat3 Sinv = Sigma.inverse();
    double det = Sigma.determinant();
    Vec3 d = x - mu;
    double mahal = d.transpose() * Sinv * d;     // squared Mahalanobis distance
    double norm = std::pow(2.0 * M_PI, -1.5) / std::sqrt(det);
    return norm * std::exp(-0.5 * mahal);
}
// The eigenvectors of cov are the PCA axes (see 03_eig-svd-pca); the smallest-
// eigenvalue axis estimates the surface normal of a locally planar patch.
```

## Connections

- [[31_calculus-foundations]] supplies the integration that defines continuous densities and expectation.
- [[03_eig-svd-pca]] diagonalizes the covariance matrix; PCA is the eigen-decomposition of $\Sigma$.
- [[06_linear-least-squares]] and [[08_nonlinear-least-squares]] are the MLE under Gaussian noise; weights are inverse variances.
- [[09_robust-ransac]] handles the case where the noise is *not* Gaussian (outliers); inlier thresholds come from the assumed noise standard deviation.
- [[16_registration-icp]] alignment quality and pose uncertainty are described by covariance.
- [[27_supervised-learning-optimization]] is statistics with many parameters; the loss is a negative log-likelihood.

## Pitfalls and failure modes

- **$n$ vs $n-1$.** Dividing by $n$ underestimates variance; use $n-1$ (Bessel's correction) for an unbiased sample covariance.
- **Gaussian assumption breaks on outliers.** A single gross outlier shifts the mean and inflates the covariance arbitrarily. That is exactly why robust methods exist ([[09_robust-ransac]]).
- **Singular covariance.** With fewer points than dimensions, or collinear points, $\Sigma$ is rank-deficient and not invertible; regularize ($\Sigma + \lambda I$) before inverting.
- **Correlation is not causation, and uncorrelated is not independent** (except for jointly Gaussian variables).
- **Probabilities vs densities.** A density can exceed 1; only its integral over a region is a probability.

## Practice

1. Implement `meanCovariance` and confirm on 1000 samples drawn from a known $\mathcal{N}(\mu,\Sigma)$ that the estimates approach the truth as $n$ grows (law of large numbers).
2. Show algebraically that minimizing $\sum_i (y_i - a x_i - b)^2$ is the MLE of a line under i.i.d. Gaussian residuals.
3. Take a locally planar point patch, compute its covariance, and use the smallest eigenvector as the normal; compare to the analytic plane normal.
4. Give three measurements with one gross outlier; compute the mean with and without the outlier and report how far the estimate moves.
5. Derive weighted least squares from the log-likelihood with per-point variances $\sigma_i^2$.

## References

- Christopher Bishop, *Pattern Recognition and Machine Learning* (Ch. 1-2: probability, Gaussians, MLE).
- Larry Wasserman, *All of Statistics*.
- Kevin Murphy, *Probabilistic Machine Learning: An Introduction* (free): https://probml.github.io/pml-book/
- Sam Roweis, "Gaussian identities" notes (multivariate Gaussian algebra).

## Flashcard seeds

- Q: Derive the i.i.d. Gaussian log-likelihood for residuals $r_i=y_i-f(x_i;\theta)$ (step 1 of 2): turn the product of densities into a log-sum. <br> <b>Formula:</b> $p(r_i)=\frac{1}{\sqrt{2\pi}\,\sigma}\exp\!\big(-\frac{r_i^2}{2\sigma^2}\big)$, target $\log\mathcal{L}(\theta)=-\frac{1}{2\sigma^2}\sum_i r_i^2+\text{const}$ <br> <i>Hint:</i> use independence to write $\mathcal{L}=\prod_i p(r_i)$ first, then take $\log$ so the product becomes a sum. :: A: Independence multiplies the densities: $\mathcal{L}(\theta)=\prod_i p(r_i)$ <br> take logs to turn the product into a sum: $\log\mathcal{L}=\sum_i\big[-\tfrac12\log(2\pi\sigma^2)-\frac{r_i^2}{2\sigma^2}\big]$ <br> collect the $\theta$-independent terms into a constant: $\boxed{\log\mathcal{L}(\theta)=-\frac{1}{2\sigma^2}\sum_i r_i^2+\text{const}}$.
- Q: From the log-likelihood, show the MLE equals least squares (step 2 of 2). <br> <b>Formula:</b> $\log\mathcal{L}=-\frac{1}{2\sigma^2}\sum_i r_i^2+\text{const}$, target $\hat\theta_{\text{MLE}}=\arg\min_\theta\sum_i r_i^2$ <br> <i>Hint:</i> $\sigma^2>0$ is a fixed positive constant — maximizing a negative multiple of $\sum r_i^2$ is the same as minimizing $\sum r_i^2$. :: A: $\sigma^2>0$ fixed, so the only $\theta$ dependence is in $-\frac{1}{2\sigma^2}\sum_i r_i^2$ <br> maximizing a negative multiple of $\sum r_i^2$ means minimizing $\sum r_i^2$ <br> $\hat\theta=\arg\max_\theta\log\mathcal{L}=\arg\min_\theta\sum_i r_i^2$ <br> $\boxed{\hat\theta_{\text{MLE}}=\arg\min_\theta\sum_i\big(y_i-f(x_i;\theta)\big)^2}$.
- Q: Derive weighted least squares from per-point variances $\sigma_i^2$ (step 1 of 2): get the log-likelihood. <br> <b>Formula:</b> $p(r_i)=\frac{1}{\sqrt{2\pi}\,\sigma_i}\exp\!\big(-\frac{r_i^2}{2\sigma_i^2}\big)$, target $\log\mathcal{L}(\theta)=-\tfrac12\sum_i\frac{r_i^2}{\sigma_i^2}+\text{const}$ <br> <i>Hint:</i> take $\log$ of $\prod_i p(r_i)$; this time $\sigma_i$ stays inside the sum, so only the $\theta$-free pieces become the constant. :: A: $\log\mathcal{L}=\sum_i\big[-\tfrac12\log(2\pi\sigma_i^2)-\frac{r_i^2}{2\sigma_i^2}\big]$ <br> the $\sigma_i$ now sit inside the sum, so drop only the additive constants in $\theta$ <br> $\boxed{\log\mathcal{L}(\theta)=-\tfrac12\sum_i\frac{r_i^2}{\sigma_i^2}+\text{const}}$.
- Q: Finish the weighted least-squares objective (step 2 of 2). <br> <b>Formula:</b> $\log\mathcal{L}=-\tfrac12\sum_i \frac{r_i^2}{\sigma_i^2}$, target $\hat\theta=\arg\min_\theta\sum_i w_i r_i^2$ <br> <i>Hint:</i> drop the $-\tfrac12$ to flip max to min, then name the per-point factor $w_i=1/\sigma_i^2$. :: A: Maximizing $-\tfrac12\sum_i r_i^2/\sigma_i^2$ is minimizing $\sum_i r_i^2/\sigma_i^2$ <br> define weights $w_i=1/\sigma_i^2$ (inverse variance) <br> $\boxed{\hat\theta=\arg\min_\theta\sum_i w_i\,r_i^2,\quad w_i=1/\sigma_i^2}$.
- Q: Derive the MAP estimate from Bayes' theorem. <br> <b>Formula:</b> $p(\theta\mid D)=\frac{p(D\mid\theta)\,p(\theta)}{p(D)}$, target $\hat\theta_{\text{MAP}}=\arg\max_\theta\big[\log\mathcal{L}(\theta)+\log p(\theta)\big]$ <br> <i>Hint:</i> $p(D)$ has no $\theta$ in it, so drop it from the argmax, then take $\log$ to split the product into a sum. :: A: $p(D)$ does not depend on $\theta$, so it does not move the argmax <br> $\hat\theta_{\text{MAP}}=\arg\max_\theta p(D\mid\theta)p(\theta)$ <br> take logs: $\arg\max_\theta[\log p(D\mid\theta)+\log p(\theta)]$ <br> $\boxed{\hat\theta_{\text{MAP}}=\arg\max_\theta\big[\log\mathcal{L}(\theta)+\log p(\theta)\big]}$, the MLE plus a prior penalty.
- Q: Derive the precision-weighted mean of two Gaussian measurements. <br> <b>Formula:</b> minimize $J(\hat\mu)=\frac{(y_1-\hat\mu)^2}{\sigma_1^2}+\frac{(y_2-\hat\mu)^2}{\sigma_2^2}$, target $\hat\mu=\frac{w_1y_1+w_2y_2}{w_1+w_2}$ <br> <i>Hint:</i> set $\frac{dJ}{d\hat\mu}=0$, substitute $w_i=1/\sigma_i^2$, then solve the linear equation for $\hat\mu$. :: A: $\frac{dJ}{d\hat\mu}=-2\frac{y_1-\hat\mu}{\sigma_1^2}-2\frac{y_2-\hat\mu}{\sigma_2^2}=0$ <br> let $w_i=1/\sigma_i^2$: $w_1(y_1-\hat\mu)+w_2(y_2-\hat\mu)=0$ <br> $\hat\mu(w_1+w_2)=w_1y_1+w_2y_2$ <br> $\boxed{\hat\mu=\frac{w_1y_1+w_2y_2}{w_1+w_2}}$, the inverse-variance weighted average.
- Q: From the multivariate Gaussian exponent, identify the squared Mahalanobis distance and show the density is constant on ellipsoids. <br> <b>Formula:</b> $p(x)\propto\exp\!\big(-\tfrac12(x-\mu)^\top\Sigma^{-1}(x-\mu)\big)$, target $d_M^2(x)=(x-\mu)^\top\Sigma^{-1}(x-\mu)$ <br> <i>Hint:</i> name the quadratic form in the exponent $d_M^2$; $p$ depends on $x$ only through it, so level sets of $p$ are level sets of $d_M^2$. :: A: The exponent is $-\tfrac12(x-\mu)^\top\Sigma^{-1}(x-\mu)$ <br> name the quadratic form $d_M^2=(x-\mu)^\top\Sigma^{-1}(x-\mu)$ <br> $p(x)$ depends on $x$ only through $d_M^2$, so $p=$ const $\iff d_M^2=$ const <br> $\boxed{d_M^2(x)=(x-\mu)^\top\Sigma^{-1}(x-\mu)=\text{const are the uncertainty ellipsoids}}$.
- Q: Show that the diagonal of the covariance matrix holds the per-axis variances. <br> <b>Formula:</b> $\Sigma_{ij}=\mathbb{E}[(X_i-\mu_i)(X_j-\mu_j)]$, target $\Sigma_{ii}=\mathrm{Var}(X_i)$ <br> <i>Hint:</i> set $i=j$ and notice the two factors become identical, giving a square. :: A: Set $i=j$: $\Sigma_{ii}=\mathbb{E}[(X_i-\mu_i)(X_i-\mu_i)]$ <br> the product is the square of the $i$-th centered component <br> $\boxed{\Sigma_{ii}=\mathbb{E}[(X_i-\mu_i)^2]=\mathrm{Var}(X_i)=\sigma_i^2}$.
- Q: Worked example (small, mm). Five edge measurements $100.2,99.8,100.1,99.9,100.0$ mm — compute the sample mean $\bar x$ and the unbiased sample SD $s$. <br> <b>Formula:</b> $\bar x=\frac{1}{n}\sum_i x_i$, $s^2=\frac{1}{n-1}\sum_i (x_i-\bar x)^2$ <br> <i>Hint:</i> compute $\bar x$ first, then the five deviations $x_i-\bar x$, square and sum them, and divide by $n-1=4$. :: A: $\bar x=\frac{500.0}{5}=100.0$ <br> deviations $+0.2,-0.2,+0.1,-0.1,0$; squares sum to $0.04+0.04+0.01+0.01+0=0.10$ <br> $s^2=\frac{0.10}{5-1}=0.025$ <br> $\boxed{s=\sqrt{0.025}\approx0.158\text{ mm},\ \text{edge}=100.0\pm0.16\text{ mm}}$.
- Q: Worked example (fusing two scans, mm). Coarse scan $y_1=10.0$ mm, $\sigma_1=0.1$; fine scan $y_2=10.6$ mm, $\sigma_2=0.3$ — fuse them. <br> <b>Formula:</b> $w_i=1/\sigma_i^2$, $\hat\mu=\frac{w_1y_1+w_2y_2}{w_1+w_2}$, $\mathrm{Var}(\hat\mu)=\frac{1}{w_1+w_2}$ <br> <i>Hint:</i> compute the two weights $w_i=1/\sigma_i^2$ first, then plug into the weighted average. :: A: $w_1=1/0.1^2=100,\ w_2=1/0.3^2=11.11$ <br> $\hat\mu=\frac{100(10.0)+11.11(10.6)}{100+11.11}=\frac{1117.8}{111.11}$ <br> $\boxed{\hat\mu\approx10.06\text{ mm},\ \mathrm{Var}=\tfrac{1}{w_1+w_2}=0.009\Rightarrow\sigma\approx0.095\text{ mm}}$ (pulled toward the precise scan).
- Q: Worked example (3D, large scale, m). A point lies $2$ m east and $3$ m north of a fitted block center, with $\Sigma=\mathrm{diag}(4,9)\ \text{m}^2$ — compute its Mahalanobis distance. <br> <b>Formula:</b> $d_M^2=(x-\mu)^\top\Sigma^{-1}(x-\mu)$, with $\Sigma^{-1}=\mathrm{diag}(1/4,1/9)$ for a diagonal $\Sigma$ <br> <i>Hint:</i> for a diagonal $\Sigma$ the quadratic form is just $\sum_k \frac{(x_k-\mu_k)^2}{\sigma_k^2}$; compute each term, then $\sqrt{}$. :: A: $d_M^2=(x-\mu)^\top\Sigma^{-1}(x-\mu)$ with $\Sigma^{-1}=\mathrm{diag}(1/4,1/9)$ <br> $d_M^2=\frac{2^2}{4}+\frac{3^2}{9}=1+1$ <br> $\boxed{d_M^2=2,\ d_M=\sqrt2\approx1.41}$ standard deviations from the center (both axes one sigma out).
- Q: Worked example (3D normalizer). Evaluate the Gaussian normalizing constant for $n=3$ and $\Sigma=I_3$. <br> <b>Formula:</b> $(2\pi)^{-n/2}\lvert\Sigma\rvert^{-1/2}$ <br> <i>Hint:</i> compute $\lvert\Sigma\rvert=\det I_3=1$ first, so the determinant factor drops, leaving $(2\pi)^{-3/2}$. :: A: $|\Sigma|=1$, so $|\Sigma|^{-1/2}=1$ <br> $(2\pi)^{-3/2}=1/(2\pi)^{1.5}$ <br> $\boxed{(2\pi)^{-3/2}\approx0.0635}$, the peak density at $x=\mu$ for a unit isotropic 3D Gaussian.
- Q: Closing test (outlier sensitivity). Three readings $9.9,10.1,10.0$ mm, then a gross outlier $20.0$ is appended — compute both means and check the corruption claim. <br> <b>Formula:</b> $\bar x=\frac{1}{n}\sum_i x_i$, then compare $\Delta\bar x=\bar x_{\text{outlier}}-\bar x_{\text{clean}}$ <br> <i>Hint:</i> a single bad point should move the mean noticeably; a pass means $\Delta\bar x$ is large relative to the clean spread. :: A: clean mean $=\frac{9.9+10.1+10.0}{3}=10.0$ <br> with outlier $=\frac{30.0+20.0}{4}=12.5$ <br> verify: the estimate moved $12.5-10.0=2.5$ mm on one bad point, confirming a single outlier shifts the Gaussian mean arbitrarily <br> $\boxed{\Delta\bar x=2.5\text{ mm}\Rightarrow\text{use a robust method (RANSAC)}}$.
- Q: Closing test (units + equal-variance limit). In weighted least squares, check that $w_i r_i^2$ is dimensionless and that WLS reduces to ordinary LS when all variances are equal. <br> <b>Formula:</b> $w_i=1/\sigma_i^2$, objective $\sum_i w_i r_i^2$ <br> <i>Hint:</i> a pass = $w_i r_i^2$ has no units (proper likelihood) AND setting all $\sigma_i=\sigma$ leaves an argmin that does not depend on $\sigma$. :: A: $\sigma_i$ has the units of $y_i$, so $w_i=1/\sigma_i^2$ has units $1/[y]^2$, making $w_i r_i^2$ dimensionless (a proper likelihood) <br> set all $\sigma_i=\sigma$: $\sum w_i r_i^2=\frac{1}{\sigma^2}\sum r_i^2$, whose argmin is independent of $\sigma$ <br> $\boxed{\text{WLS}\to\text{ordinary LS when variances are equal — limit check passes}}$.
- Q: When should you reach for weighted (not ordinary) least squares? <br> <i>Hint:</i> think about whether every point's noise has the same variance. :: A: When per-point noise variances $\sigma_i^2$ differ (heteroscedastic data) — weight each residual by $1/\sigma_i^2$ so precise points count more.
- Q: Pitfall: why divide the sample variance by $n-1$ instead of $n$? <br> <b>Formula:</b> $s^2=\frac{1}{n-1}\sum_i (x_i-\bar x)^2$ <br> <i>Hint:</i> count how many degrees of freedom the sample mean has already used up. :: A: The sample mean already absorbs one degree of freedom; dividing by $n$ underestimates spread. $n-1$ (Bessel's correction) makes the estimator unbiased.
- Q: Pitfall: when is a covariance matrix singular, and what is the fix? <br> <b>Formula:</b> regularize as $\Sigma+\lambda I$ before inverting <br> <i>Hint:</i> think about rank — how many independent points vs dimensions you have. :: A: When points are fewer than dimensions or collinear, $\Sigma$ is rank-deficient and not invertible; regularize as $\Sigma+\lambda I$ before inverting.
- Q: Code idiom: how do you estimate a surface normal from a point patch? <br> <b>Formula:</b> $\Sigma=\frac{1}{n-1}\sum_i (p_i-\mu)(p_i-\mu)^\top$, then eigendecompose $\Sigma$ <br> <i>Hint:</i> the normal is the direction of least spread — pick the eigenvector of the smallest eigenvalue. :: A: Build the $3\times3$ patch covariance, eigendecompose it, and take the eigenvector of the smallest eigenvalue (the direction of least spread).
- CLOZE: Minimizing the sum of squared residuals is the {{c1::maximum-likelihood}} estimate under i.i.d. zero-mean {{c2::Gaussian}} noise.
- CLOZE: With per-point variances, the MLE becomes {{c1::weighted}} least squares with weights {{c2::$1/\sigma_i^2$}}.
- CLOZE: A probability {{c1::density}} can exceed 1; only its {{c2::integral}} over a region is a probability.
- CLOZE: The eigenvectors of the covariance matrix are the {{c1::principal axes}} of the uncertainty ellipsoid, linking statistics to {{c2::PCA}}.
- CLOZE: Uncorrelated implies independent only for {{c1::jointly Gaussian}} variables; in general the two are {{c2::not}} equivalent.
- CLOZE: The central limit theorem justifies the Gaussian noise model because sensor error is a {{c1::sum}} of many small {{c2::independent}} errors.
