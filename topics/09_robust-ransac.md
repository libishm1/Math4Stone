# Robust estimation and RANSAC

Prereqs: [[03_linear-least-squares]], [[04_eig-svd-pca]] => this => unlocks [[15_remeshing-segmentation]], [[28_cadfit-program-synthesis]].

## Why it matters for the stone scan->CAD->IFC / CADFit workflow

A raw stone scan is a point cloud where most points lie on the real surface but a hard minority are garbage: laser speckle off wet granite, registration ghosts, dust, and the saw-blade kerf wall bleeding into the slab face. Plain least squares fits all of them at once, so a single bright outlier tilts your plane and your "sawn face" detector misses. RANSAC plus a robust loss lets you fit a plane, cylinder, or sphere to the *consensus* surface and label every other point an outlier, which is exactly how you separate a flat sawn face from rough split rock before meshing. That segmentation feeds two downstream stages: feature-aware remeshing needs to know which regions are planar primitives, and the CADFit program-synthesis layer needs clean primitive parameters (a plane normal, a cylinder axis and radius) as the tokens it assembles into a parametric CAD program. Get the robust fit wrong and every later stage inherits a biased surface.

## Intuition

Least squares is a democracy where every voter shouts with a volume equal to the *square* of how wrong it is. One liar who is ten units off contributes a hundred to the total error, so the fit bends toward the liar to shut it up. That is the whole problem: squared error gives outliers disproportionate power.

Two fixes exist and they compose. First, change the voting rule so far-off points stop getting louder past a cutoff. That is a robust loss (Huber, Tukey). Second, stop trusting the whole crowd and instead poll tiny random juries: pick the minimum number of points that define a model, fit it, then count how many of *all* points agree. Repeat many times and keep the model with the biggest agreeing majority. That is RANSAC.

Analogy: you want the true height of a room's floor, but a few people are standing on ladders. Averaging everyone gives a floor that floats. Instead, pick 3 random people, assume they define the floor plane, then count how many others are within ankle height of it. The 3 standing on the actual floor produce a plane that hundreds agree with. The 3 on ladders produce a plane almost nobody agrees with. Keep the popular plane.

## The math (first principles)

### Setup and conventions

Data are points $x_i \in \mathbb{R}^3$ (a scan in metres, right-handed world frame). A model has parameter vector $\theta$. The residual $r_i(\theta)$ is the signed or unsigned geometric distance from point $i$ to the model surface, in the same length units as the data. A point is an **inlier** if $|r_i| \le t$ for a threshold $t$; otherwise an **outlier**. Conventions: residuals carry data units; the threshold $t$ is set from sensor noise, not guessed in the abstract; all primitive normals are unit length.

### Why one outlier wrecks least squares

Ordinary least squares minimises
$$ \theta^\star = \arg\min_\theta \sum_i r_i(\theta)^2 . $$
The cost of residual $r$ is $\rho_{\text{LS}}(r) = \tfrac12 r^2$, with influence (derivative) $\psi(r) = \rho'(r) = r$. The influence is **unbounded**: a residual twice as large pulls twice as hard, with no ceiling. So a single point at residual $R$ contributes gradient proportional to $R$, and as $R \to \infty$ its pull on $\theta$ grows without bound. The breakdown point of least squares is $0\%$: one bad point arbitrarily far away can move the estimate arbitrarily far. This is the formal statement of "one outlier wrecks the fit."

### M-estimators: bound the influence

An M-estimator replaces the squared cost with a robust loss $\rho$:
$$ \theta^\star = \arg\min_\theta \sum_i \rho\!\left(\frac{r_i(\theta)}{\sigma}\right), $$
where $\sigma$ is a scale (noise estimate). Robustness comes from making the **influence function** $\psi = \rho'$ bounded, so distant points stop dominating.

**Huber loss** (with tuning constant $k$, standardised residual $u = r/\sigma$):
$$ \rho_k(u) = \begin{cases} \tfrac12 u^2, & |u| \le k, \\[4pt] k\,|u| - \tfrac12 k^2, & |u| > k. \end{cases} $$
It is quadratic near zero (efficient like least squares for small residuals) and linear in the tails (influence saturates at $\pm k$). Influence: $\psi_k(u) = \max(-k,\min(k,u))$, which is *bounded but never zero*. Outliers are downweighted, not discarded. Conventional choice $k = 1.345$ gives about $95\%$ efficiency relative to OLS under Gaussian noise.

**Tukey biweight (bisquare) loss** (tuning constant $c$):
$$ \rho_c(u) = \begin{cases} \dfrac{c^2}{6}\!\left[\,1 - \left(1 - (u/c)^2\right)^3\right], & |u| \le c, \\[8pt] \dfrac{c^2}{6}, & |u| > c. \end{cases} $$
Beyond $|u| = c$ the loss is *constant*, so its influence is exactly zero: gross outliers are fully rejected, not merely shrunk. This is "redescending." Common choice $c = 4.685$ gives about $95\%$ Gaussian efficiency. Tukey is non-convex, so it needs a good initial guess; this is exactly what RANSAC provides.

### Solving M-estimators by IRLS

M-estimators are solved by **Iteratively Reweighted Least Squares**. Define a weight $w(u) = \psi(u)/u$, then each iteration runs weighted least squares with those weights and re-evaluates. For Huber:
$$ w_k(u) = \begin{cases} 1, & |u| \le k, \\[4pt] \dfrac{k}{|u|}, & |u| > k. \end{cases} $$
For Tukey: $w_c(u) = \left(1 - (u/c)^2\right)^2$ for $|u| \le c$, else $0$. Loop: compute residuals, standardise with a robust scale $\sigma = 1.4826 \cdot \text{MAD}$, set weights, solve $\min_\theta \sum_i w_i\, r_i(\theta)^2$, repeat until $\theta$ stops moving. The constant $1.4826$ makes the median-absolute-deviation a consistent estimate of the standard deviation for Gaussian data.

### RANSAC: the core loop

RANSAC (Fischler and Bolles, 1981) does not downweight; it *searches* for the inlier set. Each iteration:

1. Draw a **minimal sample** of $s$ points uniformly at random. ($s$ is the smallest number that determines the model: $3$ for a plane, $4$ for a sphere, $2$ for a 2D line.)
2. **Fit** the model exactly to those $s$ points.
3. **Score**: count consensus, the points with $|r_i| \le t$.
4. **Keep** the model with the largest consensus set so far.

After the loop, **refit** the best model to all its inliers (least squares or an M-estimator) for a low-variance final estimate.

### How many iterations N

We want at least one iteration whose minimal sample is *all inliers*. Let $w$ be the inlier ratio (fraction of points that are inliers). The probability that one minimal sample of size $s$ is all inliers is $w^s$ (sampling with replacement, or $n \gg s$). The probability that it is *not* all-inlier is $1 - w^s$. Over $N$ independent iterations the probability of *never* drawing a clean sample is $(1 - w^s)^N$. Set the success probability to $p$:
$$ 1 - (1 - w^s)^N = p \quad\Longrightarrow\quad \boxed{\,N = \dfrac{\log(1 - p)}{\log\!\left(1 - w^s\right)}\,}. $$
The spread of $N$ is captured by $\mathrm{SD}(N) = \dfrac{\sqrt{1 - w^s}}{w^s}$; budget a margin above the mean. Three facts to internalise: $N$ grows fast as $w$ drops, $N$ grows fast as $s$ rises (prefer minimal $s$), and $N$ is independent of the total point count.

### Adaptive termination

You rarely know $w$ in advance. Start with a pessimistic $w$ and a large $N$, then after every iteration update $w$ from the best consensus set found so far ($w \leftarrow |\text{inliers}|/n$) and recompute $N$ with the same boxed formula. The loop usually stops far earlier than the worst-case bound.

### Choosing the inlier threshold t

If measurement noise on the residual is Gaussian with standard deviation $\sigma$, the squared normalised residual $(r/\sigma)^2$ is chi-squared with degrees of freedom equal to the residual's codimension. For a point-to-plane (scalar) residual that is $1$ DOF, and a common $95\%$ choice is
$$ t = 1.96\,\sigma \approx 2\sigma . $$
Set $t$ from the sensor, not from the data you are fitting.

### MSAC: a better score for free

Plain RANSAC scores by counting (a top-hat: $0$ for outliers, $1$ for inliers, flat inside). **MSAC** replaces the count with a truncated quadratic cost, scoring inliers by *how well* they fit:
$$ C = \sum_i \min\!\big(r_i^2,\; t^2\big), \qquad \text{keep the model with minimum } C. $$
Inliers contribute $r_i^2$; outliers contribute the constant $t^2$. This is the same minimal-sample search with a robust ($\rho$-style) scoring rule, and it costs nothing extra.

## Worked example

Fit a line $y = ax + b$ to five 2D points by RANSAC, $s = 2$, threshold $t = 0.5$ (residual = vertical distance).

Points: $P_1(0,0)$, $P_2(1,1)$, $P_3(2,2)$, $P_4(3,3)$, $P_5(2, 8)$. Four lie on $y = x$; $P_5$ is an outlier.

**Iteration A.** Sample $\{P_1, P_5\}$. Line through $(0,0)$ and $(2,8)$ is $y = 4x$.
Residuals $|y_i - 4x_i|$: $P_1{:}\,0$, $P_2{:}\,|1-4|=3$, $P_3{:}\,|2-8|=6$, $P_4{:}\,|3-12|=9$, $P_5{:}\,0$.
Inliers ($\le 0.5$): $\{P_1, P_5\}$. Consensus $= 2$.

**Iteration B.** Sample $\{P_2, P_4\}$. Line through $(1,1)$ and $(3,3)$ is $y = x$.
Residuals $|y_i - x_i|$: $P_1{:}\,0$, $P_2{:}\,0$, $P_3{:}\,0$, $P_4{:}\,0$, $P_5{:}\,|8-2|=6$.
Inliers: $\{P_1, P_2, P_3, P_4\}$. Consensus $= 4$.

Keep B. **Refit** by least squares on the 4 inliers: they lie exactly on $y = x$, so $a = 1$, $b = 0$. The outlier $P_5$ never touched the final fit.

**Compare to plain least squares on all 5 points.** With $\bar{x} = 1.6$, $\bar{y} = 2.8$:
$\sum (x_i-\bar x)(y_i-\bar y) = 5.2$, $\sum (x_i-\bar x)^2 = 5.2$, so $a = 5.2/5.2 = 1.0$... but recompute with the outlier's leverage: $\sum x_i = 8$, $\sum y_i = 14$, $\sum x_i y_i = 0{+}1{+}4{+}9{+}16 = 30$, $\sum x_i^2 = 0{+}1{+}4{+}9{+}4 = 18$. Then
$$ a = \frac{n\sum x_iy_i - \sum x_i \sum y_i}{n\sum x_i^2 - (\sum x_i)^2} = \frac{5(30) - 8(14)}{5(18) - 64} = \frac{150 - 112}{90 - 64} = \frac{38}{26} \approx 1.46, $$
$$ b = \bar y - a\bar x = 2.8 - 1.46(1.6) \approx 0.46. $$
Least squares reports slope $1.46$ instead of $1.0$: a $46\%$ error from one outlier. RANSAC recovers the truth.

**Iteration-count sanity check.** With $w = 4/5 = 0.8$, $s = 2$, $p = 0.99$:
$$ N = \frac{\log(1 - 0.99)}{\log(1 - 0.8^2)} = \frac{\log 0.01}{\log 0.36} = \frac{-2}{-0.4437} \approx 4.5 \to 5 . $$
Five iterations suffice for $99\%$ confidence. If the inlier ratio were $0.3$, the same formula gives $N \approx \log(0.01)/\log(1 - 0.09) \approx 48.9 \to 49$, an order of magnitude more.

## Code

Idiomatic C++ with Eigen first. A header-only RANSAC plane fitter for a stone face, with the minimal-sample fit, MSAC scoring, adaptive termination, and a least-squares refit. Symbol map is in comments.

```cpp
// ransac_plane.hpp  -- robust plane fitting for a stone scan
#pragma once
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <vector>
#include <random>
#include <cmath>
#include <limits>

namespace robust {

using Point = Eigen::Vector3d;

// A plane: unit normal n and offset d, so the surface is { x : n.x + d = 0 }.
// Signed point-plane distance r_i = n . x_i + d   (data units, e.g. metres).
struct Plane { Eigen::Vector3d n{0,0,1}; double d{0.0}; };

inline double residual(const Plane& P, const Point& x) {
    return P.n.dot(x) + P.d;                  // r_i = n . x_i + d
}

// Least-squares plane through a set: centroid + smallest-eigenvector normal.
// This is the [[04_eig-svd-pca]] step: normal = eigvec of min eigenvalue of cov.
inline Plane fitPlaneLeastSquares(const std::vector<Point>& pts,
                                  const std::vector<int>& idx) {
    Point c = Point::Zero();
    for (int i : idx) c += pts[i];
    c /= static_cast<double>(idx.size());     // centroid c = (1/m) sum x_i

    Eigen::Matrix3d C = Eigen::Matrix3d::Zero();
    for (int i : idx) {                        // scatter matrix C = sum (x-c)(x-c)^T
        Eigen::Vector3d e = pts[i] - c;
        C += e * e.transpose();
    }
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> es(C);  // ascending eigenvalues
    Plane P;
    P.n = es.eigenvectors().col(0).normalized();           // smallest eigenvalue -> normal
    P.d = -P.n.dot(c);                                      // plane passes through centroid
    return P;
}

// RANSAC/MSAC plane fit. t = inlier threshold (~ 2*sigma_sensor); p = confidence.
inline Plane fitPlaneRansac(const std::vector<Point>& pts,
                            double t, double p = 0.99,
                            int maxIters = 10000, unsigned seed = 0) {
    const int n = static_cast<int>(pts.size());
    std::mt19937 rng(seed);
    std::uniform_int_distribution<int> pick(0, n - 1);

    Plane best; double bestCost = std::numeric_limits<double>::infinity();
    std::vector<int> bestInliers;
    const int s = 3;                          // minimal sample size for a plane
    int N = maxIters;                          // adaptive cap, shrinks as w grows

    for (int it = 0; it < N && it < maxIters; ++it) {
        int a, b, cc;                          // draw 3 distinct points
        do { a = pick(rng); b = pick(rng); cc = pick(rng); }
        while (a == b || b == cc || a == cc);

        Plane cand = fitPlaneLeastSquares(pts, {a, b, cc});  // exact fit on minimal set
        if (!cand.n.allFinite()) continue;

        double cost = 0.0;                     // MSAC score: sum min(r^2, t^2)
        std::vector<int> inl;
        for (int i = 0; i < n; ++i) {
            double r = residual(cand, pts[i]);
            double r2 = r * r;
            cost += std::min(r2, t * t);
            if (r2 <= t * t) inl.push_back(i); // |r_i| <= t  => inlier
        }

        if (cost < bestCost) {                 // keep lowest-cost model
            bestCost = cost; best = cand; bestInliers = inl;

            // adaptive N: w = inlier ratio so far, recompute log(1-p)/log(1-w^s)
            double w = static_cast<double>(inl.size()) / n;
            if (w > 0.0) {
                double denom = std::log(1.0 - std::pow(w, s));
                if (denom < 0.0)
                    N = static_cast<int>(std::log(1.0 - p) / denom) + 1;
            }
        }
    }
    // Final refit on all inliers -> low-variance estimate.
    if (bestInliers.size() >= 3) best = fitPlaneLeastSquares(pts, bestInliers);
    return best;
}

} // namespace robust
```

In production you would not hand-roll this. PCL ships tested sample-consensus models. The same plane fit:

```cpp
#include <pcl/sample_consensus/ransac.h>
#include <pcl/sample_consensus/sac_model_plane.h>
// model carries the cloud; RANSAC drives the minimal-sample search
pcl::SampleConsensusModelPlane<pcl::PointXYZ>::Ptr
    model(new pcl::SampleConsensusModelPlane<pcl::PointXYZ>(cloud));
pcl::RandomSampleConsensus<pcl::PointXYZ> ransac(model);
ransac.setDistanceThreshold(0.01);   // t = inlier threshold in metres
ransac.computeModel();
std::vector<int> inliers; ransac.getInliers(inliers);
Eigen::VectorXf coeff; ransac.getModelCoefficients(coeff);  // [a b c d], plane ax+by+cz+d=0
```

Swap `SampleConsensusModelPlane` for `SampleConsensusModelSphere` (minimal sample $4$) or `SampleConsensusModelCylinder` (needs normals; axis + radius) to segment rough vs sawn faces by primitive type.

A short Python sketch of IRLS-Huber clarifies the M-estimator half, which RANSAC does not cover:

```python
import numpy as np
def irls_line_huber(x, y, k=1.345, iters=10):
    A = np.c_[x, np.ones_like(x)]          # design matrix [x 1]
    theta = np.linalg.lstsq(A, y, None)[0] # OLS init
    for _ in range(iters):
        r = y - A @ theta                  # residuals
        s = 1.4826 * np.median(np.abs(r - np.median(r)))  # robust scale (MAD)
        u = r / max(s, 1e-12)
        w = np.where(np.abs(u) <= k, 1.0, k / np.abs(u))  # Huber weights
        W = np.sqrt(w)                      # weighted least squares
        theta = np.linalg.lstsq(A * W[:, None], y * W, None)[0]
    return theta                            # [slope, intercept]
```

## Connections

- [[03_linear-least-squares]] is the engine inside both halves: RANSAC's minimal-sample fit and final refit are least-squares solves, and IRLS is weighted least squares iterated.
- [[04_eig-svd-pca]] gives the closed-form plane normal as the smallest-eigenvalue eigenvector of the scatter matrix, the exact step `fitPlaneLeastSquares` performs.
- [[15_remeshing-segmentation]] consumes RANSAC primitive labels: knowing which points form a plane or cylinder drives feature-aware remeshing and region boundaries.
- [[28_cadfit-program-synthesis]] takes the robust primitive parameters (normal, axis, radius) as tokens it composes into a parametric CAD program, so a biased fit corrupts the synthesised program.
- ICP and registration (a sibling topic) use the same robust-loss machinery to reject bad correspondences during alignment.

## Pitfalls and failure modes

- **Frame and units in the threshold.** $t$ is a distance in data units. Fitting a scan in millimetres with $t = 0.01$ assuming metres rejects everything. Set $t \approx 2\sigma_{\text{sensor}}$ in the actual units.
- **Degenerate minimal samples.** Three collinear points do not define a plane; four coplanar points do not define a sphere. The scatter matrix becomes rank-deficient and the normal is garbage. Reject samples whose smallest eigenvalue is near zero, or whose points are nearly collinear.
- **Inlier ratio too low.** $N = \log(1-p)/\log(1-w^s)$ explodes as $w \to 0$. At $w = 0.1$, $s = 4$, $p = 0.99$ you need about $46{,}000$ iterations. Reduce $s$, pre-filter, or use a guided sampler (PROSAC) if $w$ is small.
- **Threshold too large vs too small.** Too large and every hypothesis looks equally good, so the search cannot discriminate; too small and the estimate is unstable and the true model gets few inliers. This is the central RANSAC tuning tension.
- **Counting ties (plain RANSAC).** The top-hat score cannot rank two models with equal inlier counts. MSAC's truncated-quadratic cost breaks ties by fit quality and is a free upgrade.
- **No final refit.** Returning the minimal-sample model directly is high variance: it was fit to only $s$ points. Always refit on the full inlier set.
- **Non-convex robust losses need a seed.** Tukey biweight has zero influence in the tails but a non-convex cost. Start IRLS-Tukey from a RANSAC or Huber estimate, not from cold OLS, or it lands in a bad local minimum.
- **Outliers that form a structure.** RANSAC finds the *largest* consensus set. If a flat saw-kerf wall has more points than the slab face, RANSAC fits the wall. Segment iteratively (fit, remove inliers, refit) and check primitive plausibility.

## Practice

1. **By hand.** Repeat the worked-example line fit but move $P_5$ to $(2, 100)$. Confirm RANSAC still returns $y = x$ while OLS slope blows up further. Report both slopes.
2. **Iteration budget.** Tabulate $N$ from $N = \log(1-p)/\log(1-w^s)$ for $p = 0.99$, $s \in \{2,3,4\}$, and $w \in \{0.9, 0.7, 0.5, 0.3\}$. Mark which cells exceed $10{,}000$ iterations.
3. **Plane on synthetic stone.** Generate $1000$ points on $z = 0$ with Gaussian noise $\sigma = 0.005$, plus $150$ uniform outliers in a box. Fit with `fitPlaneRansac` and with plain least squares. Print both normals and the angle between them.
4. **MSAC vs counting.** Modify the scorer to also track the plain inlier count. On a cloud with two parallel planes of nearly equal size, report whether count-RANSAC and MSAC pick the same plane and why.
5. **IRLS sweep.** Run `irls_line_huber` with $k \in \{0.5, 1.345, 5.0\}$ on the outlier-contaminated line. Plot the fitted slope vs $k$ and explain the bias-robustness trade-off.
6. **Sphere fitting.** Swap to PCL `SampleConsensusModelSphere` on a scanned cobble. Report centre, radius, and inlier fraction, then visualise inliers vs outliers as two coloured clouds.

## References

- Fischler, M. A. and Bolles, R. C. "Random Sample Consensus: A Paradigm for Model Fitting." *Communications of the ACM*, 24(6), 1981. https://dl.acm.org/doi/10.1145/358669.358692
- Random sample consensus (RANSAC) overview, iteration-count and threshold discussion. https://en.wikipedia.org/wiki/Random_sample_consensus
- Torr, P. H. S. and Zisserman, A. "MLESAC: A New Robust Estimator with Application to Estimating Image Geometry" (MSAC/MLESAC scoring). *CVIU*, 2000. https://www.robots.ox.ac.uk/~vgg/publications/2000/Torr00/torr00.pdf
- Huber, P. J. *Robust Statistics*. Wiley, 1981 (Huber loss and M-estimators). https://onlinelibrary.wiley.com/doi/book/10.1002/0471725250
- "Robust Regression" R Data Analysis Examples, UCLA OARC (Huber/Tukey weights, IRLS, $k=1.345$, $c=4.685$). https://stats.oarc.ucla.edu/r/dae/robust-regression/
- Eberly, D. "Least Squares Fitting of Data by Linear or Quadratic Structures" (plane/sphere fits via eigenanalysis). https://www.geometrictools.com/Documentation/LeastSquaresFitting.pdf
- PCL: How to use Random Sample Consensus model (plane/sphere C++ API). https://pointclouds.org/documentation/tutorials/random_sample_consensus.html
- OpenCV calib3d: RANSAC-based estimators (`findHomography`, `solvePnPRansac`). https://docs.opencv.org/4.x/d9/d0c/group__calib3d.html
- CADFit: learning CAD program synthesis. arXiv:2605.01171; code https://github.com/ghadinehme/CADFit

## Flashcard seeds

- Q: Derive the RANSAC iteration count $N$ (step 1 of 3): the probability one minimal sample of $s$ points is all-inlier, given inlier ratio $w$. <br> <b>Formula:</b> $P(\text{one draw is inlier})=w$, and independent events multiply. <br> <i>Hint:</i> the $s$ draws are independent, so multiply $w$ by itself $s$ times. :: A: Sampling is independent (with replacement, or $n\gg s$), so each draw is an inlier with probability $w$ <br> the $s$ draws are independent, so probabilities multiply: $w\cdot w\cdots w$ ($s$ times) <br> $\boxed{P(\text{sample all-inlier}) = w^s}$
- Q: Derive the RANSAC iteration count $N$ (step 2 of 3): the probability that across $N$ independent iterations you NEVER draw a clean sample. <br> <b>Formula:</b> start from $P(\text{clean})=w^s$; complement and independence give $(1-w^s)^N$. <br> <i>Hint:</i> first write the chance ONE iteration is not clean, then raise it to the $N$th power. :: A: One iteration fails to be clean with probability $1-w^s$ <br> the $N$ iterations are independent, so failures multiply <br> $\boxed{P(\text{no clean sample in }N) = (1-w^s)^N}$
- Q: Derive the RANSAC iteration count $N$ (step 3 of 3): solve for $N$ given target success probability $p$. <br> <b>Formula:</b> start from $1-(1-w^s)^N = p$; target $N = \dfrac{\log(1-p)}{\log(1-w^s)}$. <br> <i>Hint:</i> isolate $(1-w^s)^N$, then take $\log$ of both sides to pull $N$ down. :: A: Success means at least one clean sample: $1-(1-w^s)^N = p$ <br> rearrange: $(1-w^s)^N = 1-p$ <br> take logs: $N\log(1-w^s) = \log(1-p)$ <br> $\boxed{N = \dfrac{\log(1-p)}{\log(1-w^s)}}$
- Q: From the OLS cost, derive the influence function $\psi$ and state the breakdown point. <br> <b>Formula:</b> $\rho_{\text{LS}}(r)=\tfrac12 r^2$, influence $\psi(r)=\rho'(r)$. <br> <i>Hint:</i> differentiate $\tfrac12 r^2$ with respect to $r$, then ask if the result is bounded as $r\to\infty$. :: A: Influence is the derivative $\psi(r)=\rho'(r)$ <br> $\psi(r)=\frac{d}{dr}\big(\tfrac12 r^2\big)=r$ <br> $\psi$ is unbounded, so as $r\to\infty$ the pull on $\theta$ grows without limit <br> $\boxed{\psi_{\text{LS}}(r)=r,\ \text{breakdown point}=0\%}$
- Q: Derive the Huber tail influence (step 1 of 2): find $\psi_k(u)$ for $|u|>k$. <br> <b>Formula:</b> Huber tail loss $\rho_k(u)=k|u|-\tfrac12 k^2$, influence $\psi_k=\rho_k'$. <br> <i>Hint:</i> differentiate term by term and use $\frac{d}{du}|u|=\operatorname{sgn}(u)$. :: A: $\psi_k(u)=\rho_k'(u)=\frac{d}{du}\big(k|u|-\tfrac12 k^2\big)$ <br> $\frac{d}{du}|u|=\operatorname{sgn}(u)$ for $u\ne 0$ <br> $\boxed{\psi_k(u)=k\,\operatorname{sgn}(u)}$ (saturates at $\pm k$)
- Q: Derive the Huber IRLS weight in the tail (step 2 of 2): find $w_k(u)$ for $|u|>k$. <br> <b>Formula:</b> $\psi_k(u)=k\,\operatorname{sgn}(u)$, weight $w(u)=\psi(u)/u$. <br> <i>Hint:</i> divide by $u$ and simplify $\operatorname{sgn}(u)/u=1/|u|$. :: A: $w_k(u)=\dfrac{\psi_k(u)}{u}=\dfrac{k\,\operatorname{sgn}(u)}{u}$ <br> $\dfrac{\operatorname{sgn}(u)}{u}=\dfrac{1}{|u|}$ <br> $\boxed{w_k(u)=\dfrac{k}{|u|}}$ for $|u|>k$ (and $w_k=1$ for $|u|\le k$)
- Q: From the Tukey loss, derive the IRLS weight for $|u|\le c$. <br> <b>Formula:</b> $\rho_c(u)=\tfrac{c^2}{6}\big[1-(1-(u/c)^2)^3\big]$, then $\psi_c=\rho_c'$ and $w_c=\psi_c/u$; target $w_c=(1-(u/c)^2)^2$. <br> <i>Hint:</i> differentiate the cube via the chain rule (outer $3(\cdot)^2$, inner $-2u/c^2$), then divide by $u$. :: A: $\psi_c(u)=\rho_c'(u)=\tfrac{c^2}{6}\cdot\big(-3(1-(u/c)^2)^2\big)\cdot(-2u/c^2)$ <br> $=u\,(1-(u/c)^2)^2$ <br> $w_c(u)=\psi_c(u)/u$ <br> $\boxed{w_c(u)=(1-(u/c)^2)^2}$ for $|u|\le c$, else $0$
- Q: Show the MSAC score caps an outlier's contribution at $t^2$. <br> <b>Formula:</b> $C=\sum_i\min(r_i^2,t^2)$. <br> <i>Hint:</i> split into the inlier case $|r_i|\le t$ and outlier case $|r_i|>t$, and ask which argument the $\min$ selects. :: A: Top-hat RANSAC counts inliers but gives outliers influence $0$ and inliers influence $1$, ignoring fit quality <br> truncate the squared cost at the threshold: contribution $=\min(r_i^2,t^2)$ <br> if $|r_i|\le t$ it pays $r_i^2$; if $|r_i|>t$ then $r_i^2>t^2$ so the min selects $t^2$ <br> $\boxed{\text{outlier contributes the constant }t^2}$
- Q: Derive the inlier threshold $t$ for a scalar (point-to-plane) residual with Gaussian noise $\sigma$ at 95%. <br> <b>Formula:</b> $(r/\sigma)^2\sim\chi^2_1$; the $95\%$ one-DOF cutoff is $1.96^2$. <br> <i>Hint:</i> set $(r/\sigma)^2\le 1.96^2$ and take the square root to bound $|r|$. :: A: $(r/\sigma)^2\sim\chi^2$ with $1$ DOF (residual codimension $1$) <br> 95% one-DOF cutoff: $(r/\sigma)^2\le 1.96^2$ <br> so $|r|\le 1.96\,\sigma$ <br> $\boxed{t=1.96\,\sigma\approx 2\sigma}$
- Q: Worked example (2D, small). Fit $y=ax+b$ by RANSAC to $P_1(0,0),P_2(1,1),P_3(2,2),P_4(3,3),P_5(2,8)$ with $t=0.5$; sample $\{P_2,P_4\}$ and report the model and consensus count. <br> <b>Formula:</b> slope $a=\dfrac{y_2-y_1}{x_2-x_1}$, residual $|y_i-(ax_i+b)|$, inlier if $\le t$. <br> <i>Hint:</i> compute the slope from the two sampled points first, then score all five residuals. :: A: Line through $(1,1),(3,3)$ has slope $1$, intercept $0$, so $y=x$ <br> residuals $|y_i-x_i|$: $P_1{=}0,P_2{=}0,P_3{=}0,P_4{=}0,P_5{=}|8-2|=6$ <br> inliers ($\le0.5$): $P_1,P_2,P_3,P_4$ <br> $\boxed{\text{consensus}=4,\ \text{model }y=x}$
- Q: Worked example (the outlier sample). Same 5 points and $t=0.5$; instead sample $\{P_1,P_5\}=\{(0,0),(2,8)\}$ and report the model and consensus. <br> <b>Formula:</b> slope $a=\dfrac{y_2-y_1}{x_2-x_1}$, residual $|y_i-ax_i|$, inlier if $\le t$. <br> <i>Hint:</i> get the slope through the two points first; note one is the outlier, so expect a tiny consensus. :: A: Slope $=(8-0)/(2-0)=4$, through origin, so $y=4x$ <br> residuals $|y_i-4x_i|$: $P_1{=}0,P_2{=}3,P_3{=}6,P_4{=}9,P_5{=}0$ <br> inliers ($\le0.5$): only $P_1,P_5$ <br> $\boxed{\text{consensus}=2}$, far worse than the clean sample
- Q: Worked example (OLS contrast). For the same 5 points compute the plain least-squares slope. Givens: $n=5,\ \sum x_i=8,\ \sum y_i=14,\ \sum x_iy_i=30,\ \sum x_i^2=18$. <br> <b>Formula:</b> $a=\dfrac{n\sum x_iy_i-\sum x_i\sum y_i}{n\sum x_i^2-(\sum x_i)^2}$. <br> <i>Hint:</i> plug the five given sums straight in; compute the numerator and denominator separately. :: A: $a=\dfrac{n\sum x_iy_i-\sum x_i\sum y_i}{n\sum x_i^2-(\sum x_i)^2}$ <br> $=\dfrac{5(30)-8(14)}{5(18)-64}=\dfrac{150-112}{90-64}=\dfrac{38}{26}$ <br> $\boxed{a\approx1.46}$ vs the true $1.0$: a $46\%$ error from one outlier
- Q: Worked example (iteration budget, easy regime). With $w=0.8$, $s=2$, $p=0.99$, compute $N$. <br> <b>Formula:</b> $N=\dfrac{\log(1-p)}{\log(1-w^s)}$. <br> <i>Hint:</i> compute $w^s=0.8^2$ first, then $1-w^s$, then the ratio of logs; round up. :: A: $w^s=0.8^2=0.64$, so $1-w^s=0.36$ <br> $N=\dfrac{\log 0.01}{\log 0.36}=\dfrac{-2}{-0.4437}$ <br> $\approx4.5$ <br> $\boxed{N=5\text{ iterations}}$
- Q: Worked example (iteration budget, hard regime). Same $p=0.99$, $s=2$, but the inlier ratio drops to $w=0.3$; compute $N$. <br> <b>Formula:</b> $N=\dfrac{\log(1-p)}{\log(1-w^s)}$. <br> <i>Hint:</i> recompute $w^s=0.3^2$ and $1-w^s$; the denominator log is now tiny, so $N$ balloons. :: A: $w^s=0.3^2=0.09$, so $1-w^s=0.91$ <br> $N=\dfrac{\log 0.01}{\log 0.91}=\dfrac{-2}{-0.04096}$ <br> $\approx48.9$ <br> $\boxed{N=49}$, about $10\times$ the $w=0.8$ case
- Q: Worked example (3D, larger $s$, low $w$). Sphere fit with $w=0.1$, $s=4$, $p=0.99$; estimate $N$. <br> <b>Formula:</b> $N=\dfrac{\log(1-p)}{\log(1-w^s)}$, with $\log(1-x)\approx -x$ for small $x$. <br> <i>Hint:</i> $w^s=0.1^4$ is tiny, so approximate $\log(1-w^s)\approx -w^s$ before dividing. :: A: $w^s=0.1^4=10^{-4}$, so $1-w^s\approx0.9999$ <br> $N=\dfrac{\log 0.01}{\log 0.9999}\approx\dfrac{-2}{-4.34\times10^{-5}}$ <br> $\boxed{N\approx46{,}000}$ — minimal $s$ and pre-filtering matter when $w$ is small
- Q: Closing test (verify $N$). You computed $N=5$ for $w=0.8,s=2$; check what success probability $5$ iterations actually deliver. <br> <b>Formula:</b> $P(\text{success})=1-(1-w^s)^N$. <br> <i>Hint:</i> a pass is $P_{\text{actual}}\ge p=0.99$; plug $1-w^s=0.36$ and $N=5$ in. :: A: $P(\text{success})=1-(1-w^s)^N=1-(0.36)^5$ <br> $0.36^5\approx0.0060$ <br> $1-0.0060=0.994$ <br> $\boxed{p_{\text{actual}}\approx0.994\ge0.99}$ — the rounded-up $N$ meets the target
- Q: Closing test (units check on $t$). A scan is in millimetres with sensor noise $\sigma=2\,\text{mm}$; a colleague sets $t=0.01$. Find the correct $t$ and show why $0.01$ fails. <br> <b>Formula:</b> $t\approx 2\sigma$, in the SAME units as the data. <br> <i>Hint:</i> a pass keeps the true surface points; compute $2\sigma$ in mm and compare it to $0.01$ mm. :: A: Correct $t\approx2\sigma=2(2\,\text{mm})=4\,\text{mm}$ <br> $t=0.01$ has no stated unit; if read as mm it admits only $|r|\le0.01\,\text{mm}\ll\sigma$ <br> almost every true point becomes an outlier <br> $\boxed{t=4\,\text{mm},\ \text{never a unitless }0.01}$
- Q: Closing test (Huber continuity). Verify the Huber weight $w_k(u)$ is continuous at the knee $|u|=k$. <br> <b>Formula:</b> $w_k=1$ for $|u|\le k$ and $w_k=k/|u|$ for $|u|>k$. <br> <i>Hint:</i> a pass is the two pieces agreeing; evaluate the tail branch at $|u|=k$. :: A: Quadratic side: $w_k(u)=1$ for $|u|\le k$ <br> tail side: $w_k(u)=k/|u|$, evaluate at $|u|=k$ <br> $k/k=1$ <br> $\boxed{w_k(k^-)=w_k(k^+)=1}$ — the weight is continuous, as a valid loss requires
- CLOZE: The probability a single minimal sample is all-inlier is {{c1::$w^s$}}, so the chance of no clean sample in $N$ tries is {{c2::$(1-w^s)^N$}}.
- CLOZE: In RANSAC, $p$ is the desired success probability, $w$ is the {{c1::inlier ratio}}, and $s$ is the {{c2::minimal sample size}}.
- CLOZE: M-estimators are solved by {{c1::iteratively reweighted least squares}}, where the weight is $w(u)={{c2::\psi(u)/u}}$. <br> <i>hint:</i> influence over residual.
- CLOZE: Adaptive RANSAC recomputes $N$ each iteration after updating $w$ to {{c1::inliers found so far / total points}}.
- CLOZE: A least-squares plane normal is the {{c1::smallest-eigenvalue}} eigenvector of the centred scatter matrix, and the plane passes through the {{c2::centroid}}.
- Q: When to use a robust loss (Huber/Tukey) versus RANSAC? <br> <i>Hint:</i> think contamination level and whether you already have a good initial estimate. :: A: Robust loss/IRLS handles moderate contamination from a good seed and gives a smooth downweighted estimate; RANSAC searches for the inlier set when outliers are gross or numerous; they compose (RANSAC seeds IRLS).
- Q: Minimal sample size $s$ for a 2D line, a plane, and a sphere? <br> <i>Hint:</i> $s$ is the smallest count of points that uniquely determines the model. :: A: $2$ (line), $3$ (plane), $4$ (sphere) — the smallest count that determines the model.
- Q: Key difference between Huber and Tukey influence in the tails? <br> <i>Hint:</i> ask whether the tail influence stays nonzero or returns to exactly zero. :: A: Huber influence saturates to a nonzero constant $\pm k$ (downweight, never reject); Tukey redescends to exactly zero past $|u|=c$ (full rejection).
- Q: Pitfall: why must you always refit the best RANSAC model on all its inliers? <br> <i>Hint:</i> count how many points the kept model was actually fit to. :: A: The kept model was fit to only $s$ points and is high variance; refitting on the full inlier set lowers variance and removes minimal-sample bias.
- Q: Pitfall: why can RANSAC return the saw-kerf wall instead of the slab face? <br> <i>Hint:</i> RANSAC optimises one quantity across hypotheses — which one? :: A: RANSAC keeps the LARGEST consensus set; if the wall has more points it wins. Segment iteratively (fit, remove inliers, refit) and check primitive plausibility.
- Q: Pitfall: why does Tukey IRLS need a RANSAC or Huber seed? <br> <i>Hint:</i> recall the shape of the redescending cost surface. :: A: Its redescending cost is non-convex, so a cold OLS start can land in a bad local minimum; a good seed places it in the correct basin.
- Q: How do you set a robust scale $\sigma$ from residuals for IRLS standardisation? <br> <b>Formula:</b> $\sigma=1.4826\cdot\text{MAD}$, $\text{MAD}=\operatorname{median}_i|r_i-\operatorname{median}_j r_j|$. <br> <i>Hint:</i> compute the median of residuals, then the median absolute deviation, then scale by $1.4826$. :: A: $\sigma=1.4826\cdot\text{MAD}$; the constant $1.4826$ makes the median absolute deviation a consistent estimate of the Gaussian standard deviation.
