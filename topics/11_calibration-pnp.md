# Camera calibration, distortion, and PnP

Prereqs: [[09_camera-pinhole]], [[10_nonlinear-least-squares]] => this ([[11_calibration-pnp]]) => unlocks [[12_registration-icp]].

## Why it matters for the stone scan->CAD->IFC / CADFit workflow

Every photo you take of an irregular stone is a 2D shadow of a 3D object, and the camera that took it bends straight edges into curves and shifts pixels off where a perfect pinhole would put them. Calibration recovers the camera's true internal geometry (focal length, principal point, lens distortion) so that pixels become metric rays you can trust. PnP then takes a handful of points you can recognize in both the image and the 3D model (for example marker corners on a fixture next to the stone, or feature points on a prior scan) and recovers the rigid pose of the camera relative to that object. That pose is the seed transform you hand to ICP to register a photogrammetry or structured-light cloud into the same frame as your CAD template, and it is the camera-to-world link a CADFit-style synthesis layer needs before it can reason about which primitives explain the observed silhouette. The central lesson is that vision never outputs the final answer; it outputs a transform plus a confidence (the reprojection error), and the rest of the pipeline consumes that estimate.

## Intuition

A pinhole camera is a perfect projector: a 3D point, a straight ray through the lens center, a dot on the sensor. Real lenses are not pinholes. They are stacks of curved glass that pull light slightly off the ideal ray, so straight lines near the image edge bow outward (barrel) or inward (pincushion). Calibration is the act of measuring that bending once, with a known object (a checkerboard), so you can undo it forever after.

Think of calibration like adjusting a bathroom scale before you weigh anything. You place known weights on it (the checkerboard squares are your known geometry), read what the scale says (the detected corner pixels), and turn the dial until reading matches truth. The "dial" here is a bundle of numbers: focal lengths, principal point, and distortion coefficients. PnP is the next step: now that the scale is honest, you put an unknown object on it (a new view) and read off its pose. The residual you minimize in both cases is the same: reprojection error, the pixel gap between where a 3D point lands when you project it and where you actually saw it.

## The math (first principles)

### Conventions

- World/object points $X = (X, Y, Z)$ in metric units (meters for stone, millimeters for shop fixtures; stay consistent).
- Camera frame: $+Z$ points forward out of the lens, $+X$ right, $+Y$ down. Right-handed.
- Image coordinates $(u, v)$ in pixels, origin top-left, $u$ right, $v$ down.
- $R \in SO(3)$, $t \in \mathbb{R}^3$ map world to camera: $X_c = R X_w + t$. So $[R \mid t]$ is the extrinsic, the pose of the world expressed in the camera frame.
- Rotation is parameterized minimally as a Rodrigues vector $\mathbf{r} \in \mathbb{R}^3$, where $\theta = \lVert \mathbf{r} \rVert$ is the angle and $\hat{\mathbf{r}} = \mathbf{r}/\theta$ the axis. This is 3 numbers, not 9, which keeps the optimization on the manifold.

### The pinhole projection chain

OpenCV's model is a fixed pipeline: world -> camera -> normalize -> distort -> pixels. The undistorted core is the projective equation

$$ s \begin{bmatrix} u \\ v \\ 1 \end{bmatrix} = K \, [R \mid t] \begin{bmatrix} X \\ Y \\ Z \\ 1 \end{bmatrix}, \qquad K = \begin{bmatrix} f_x & 0 & c_x \\ 0 & f_y & c_y \\ 0 & 0 & 1 \end{bmatrix}. $$

Here $s$ is an arbitrary scale (depth) absorbed by the homogeneous equality, $K$ is the intrinsic matrix, $f_x, f_y$ are focal lengths in pixels, and $(c_x, c_y)$ is the principal point. The five intrinsic numbers $(f_x, f_y, c_x, c_y, \text{skew})$ are fixed for a lens; skew is almost always assumed zero.

### Normalized camera coordinates

Project to the $z=1$ plane first. With $X_c = R X_w + t = (X', Y', Z')$,

$$ x' = \frac{X'}{Z'}, \qquad y' = \frac{Y'}{Z'}. $$

These are dimensionless. Distortion acts here, before $K$ scales them to pixels.

### Distortion model (radial + tangential, rational extension)

Let $r^2 = x'^2 + y'^2$. The distorted normalized coordinates are

$$ x'' = x' \, \frac{1 + k_1 r^2 + k_2 r^4 + k_3 r^6}{1 + k_4 r^2 + k_5 r^4 + k_6 r^6} + 2 p_1 x' y' + p_2 (r^2 + 2 x'^2), $$

$$ y'' = y' \, \frac{1 + k_1 r^2 + k_2 r^4 + k_3 r^6}{1 + k_4 r^2 + k_5 r^4 + k_6 r^6} + p_1 (r^2 + 2 y'^2) + 2 p_2 x' y'. $$

- $k_1, k_2, k_3$ (and the denominator $k_4, k_5, k_6$ in the rational model) are **radial** coefficients. They model the symmetric barrel/pincushion bending. Most lenses need only $k_1, k_2, k_3$; the denominator is opt-in via `CALIB_RATIONAL_MODEL`.
- $p_1, p_2$ are **tangential** coefficients, from lens elements not being perfectly parallel to the sensor (decentering).
- OpenCV packs the vector in the order $(k_1, k_2, p_1, p_2, k_3, [k_4, k_5, k_6, \dots])$. The transposed $p_1/p_2$ ordering is a classic bug source.

Finally apply $K$ to get pixels:

$$ u = f_x \, x'' + c_x, \qquad v = f_y \, y'' + c_y. $$

Call this whole composed map $\pi(K, \mathbf{d}, R, t; X)$, where $\mathbf{d}$ is the distortion vector. This is exactly `cv::projectPoints`.

### Calibration as nonlinear least squares

Given $M$ views of a planar target with $N$ known 3D corners $X_j$ and detected pixels $\hat{u}_{ij}$, calibration solves

$$ \min_{K, \, \mathbf{d}, \, \{R_i, t_i\}} \ \sum_{i=1}^{M} \sum_{j=1}^{N} \big\lVert \, \hat{u}_{ij} - \pi(K, \mathbf{d}, R_i, t_i; X_j) \, \big\rVert^2. $$

This is one big nonlinear least-squares problem ([[10_nonlinear-least-squares]]) over intrinsics, distortion, and one pose per view. The residual $r_{ij} = \hat{u}_{ij} - \pi(\cdot)$ is the **reprojection error** in pixels. OpenCV solves it with Levenberg-Marquardt, using the analytic Jacobians `projectPoints` returns. The optimizer is initialized in closed form: Zhang's method estimates a homography per planar view, factors out a linear guess for $K$, then `solvePnP` gives each initial $(R_i, t_i)$, with distortion started at zero.

The RMS reprojection error reported at the end is your confidence number. Sub-pixel (well under 1.0) is healthy. A large value means bad corner detections, too few or too-coplanar views, or an over-parameterized distortion model fitting noise.

### PnP: pose from known intrinsics

Once $K$ and $\mathbf{d}$ are known, **Perspective-n-Point** recovers a single pose from $n$ 2D-3D correspondences:

$$ \min_{R, \, t} \ \sum_{j=1}^{n} \big\lVert \, \hat{u}_j - \pi(K, \mathbf{d}, R, t; X_j) \, \big\rVert^2 . $$

Same residual, far fewer unknowns: just 6 DoF $(R, t)$. Counting DoF, each correspondence gives 2 equations and the pose has 6 unknowns, so you need at least 3 points (P3P, up to 4 discrete solutions disambiguated by a 4th) and in practice $\geq 4$. OpenCV's `solvePnP` flags:

- `SOLVEPNP_ITERATIVE`: LM refinement, needs an initial guess or a planar/6+ configuration. The general-purpose default.
- `SOLVEPNP_P3P` / `SOLVEPNP_AP3P`: closed-form minimal solvers, exactly 4 points.
- `SOLVEPNP_EPnP`: efficient $O(n)$ non-iterative solver, good initializer.
- `SOLVEPNP_IPPE` / `SOLVEPNP_IPPE_SQUARE`: for coplanar points / square markers.
- `SOLVEPNP_SQPNP`: globally optimal, works from 3 points.

The output `rvec` is a Rodrigues vector; convert with `cv::Rodrigues` to get $R$. `tvec` is $t$. Together $X_c = R X_w + t$.

### Rodrigues map (rotation vector to matrix)

For $\mathbf{r}$ with $\theta = \lVert \mathbf{r} \rVert$, $\hat{\mathbf{r}} = \mathbf{r}/\theta$, and $K_\times$ the skew-symmetric cross-product matrix of $\hat{\mathbf{r}}$,

$$ R = I + \sin\theta \, K_\times + (1 - \cos\theta) \, K_\times^2 . $$

This is the exponential map $\exp([\mathbf{r}]_\times)$ on $SO(3)$. Using 3 parameters instead of a 9-entry matrix keeps every LM step a valid rotation.

## Worked example

Take a trivial pinhole, no distortion ($\mathbf{d} = 0$), with

$$ K = \begin{bmatrix} 800 & 0 & 320 \\ 0 & 800 & 240 \\ 0 & 0 & 1 \end{bmatrix}. $$

Place the camera so the pose is identity ($R = I$, $t = 0$): the camera frame equals the world frame. Project the 3D point $X = (0.1, -0.05, 2.0)$ meters.

Normalize:

$$ x' = \frac{0.1}{2.0} = 0.05, \qquad y' = \frac{-0.05}{2.0} = -0.025. $$

No distortion, so $x'' = x' = 0.05$, $y'' = y' = -0.025$. Apply $K$:

$$ u = 800 \cdot 0.05 + 320 = 40 + 320 = 360, $$
$$ v = 800 \cdot (-0.025) + 240 = -20 + 240 = 220. $$

So $X$ lands at pixel $(360, 220)$. Now check reprojection error: suppose the corner detector actually found $(361.5, 219.0)$. The residual is

$$ r = (361.5 - 360, \ 219.0 - 220) = (1.5, -1.0), \qquad \lVert r \rVert = \sqrt{1.5^2 + 1.0^2} = \sqrt{3.25} \approx 1.80 \text{ px}. $$

That 1.8 px is exactly what calibration and PnP sum-of-squares minimize. Add a mild radial distortion $k_1 = 0.1$ to see its effect: $r^2 = 0.05^2 + 0.025^2 = 0.003125$, so the radial factor is $1 + 0.1 \cdot 0.003125 = 1.0003125$. Then $x'' = 0.05 \cdot 1.0003125 = 0.05001563$ and $u = 800 \cdot 0.05001563 + 320 = 360.0125$. The point moved $\approx 0.0125$ px outward. Distortion is small near the center and grows with $r^2$ toward the edges, which is why corner coverage across the whole frame matters during calibration.

## Code

C++ with OpenCV (the calib3d module is the production-standard implementation of this exact math).

```cpp
#include <opencv2/calib3d.hpp>
#include <opencv2/core.hpp>
#include <iostream>
#include <vector>

int main() {
    // --- Intrinsics K and distortion d (would come from calibrateCamera) ---
    cv::Matx33d K(800, 0, 320,
                    0, 800, 240,
                    0,   0,   1);
    // distCoeffs order: (k1, k2, p1, p2, k3) -- here zero for clarity.
    cv::Mat distCoeffs = cv::Mat::zeros(5, 1, CV_64F);

    // --- Known 3D object points X_w (a small planar fiducial, meters) ---
    std::vector<cv::Point3f> objectPoints = {
        {0.00f, 0.00f, 0.0f}, {0.10f, 0.00f, 0.0f},
        {0.10f, 0.10f, 0.0f}, {0.00f, 0.10f, 0.0f}
    };

    // --- Matching 2D detections u_hat (pixels) from the image ---
    std::vector<cv::Point2f> imagePoints = {
        {360.0f, 220.0f}, {410.0f, 222.0f},
        {408.0f, 272.0f}, {358.0f, 270.0f}
    };

    // --- PnP: solve min over (R,t) of sum || u_hat - pi(K,d,R,t;X) ||^2 ---
    cv::Mat rvec, tvec;                       // rvec = Rodrigues vector r, tvec = t
    bool ok = cv::solvePnP(objectPoints, imagePoints, K, distCoeffs,
                            rvec, tvec, false, cv::SOLVEPNP_ITERATIVE);
    if (!ok) { std::cerr << "PnP failed\n"; return 1; }

    cv::Matx33d R;                            // R = exp([r]_x), the rotation
    cv::Rodrigues(rvec, R);                   // rotation vector -> matrix
    std::cout << "R =\n" << cv::Mat(R) << "\nt = " << tvec.t() << "\n";

    // --- Confidence: recompute the reprojection error pi(...) ourselves ---
    std::vector<cv::Point2f> reproj;
    cv::projectPoints(objectPoints, rvec, tvec, K, distCoeffs, reproj);
    double sse = 0.0;
    for (size_t j = 0; j < imagePoints.size(); ++j) {
        cv::Point2f e = imagePoints[j] - reproj[j];   // residual r_j in pixels
        sse += e.x * e.x + e.y * e.y;
    }
    double rms = std::sqrt(sse / imagePoints.size());
    std::cout << "RMS reprojection error = " << rms << " px\n";
    // Vision returns a transform (R,t) AND a confidence (rms), not a final answer.
    return 0;
}
```

Calibration itself, abbreviated, is the same residual summed over many views:

```cpp
// objectPointsList[i] = 3D corners of view i ; imagePointsList[i] = detected pixels
cv::Mat K, distCoeffs;
std::vector<cv::Mat> rvecs, tvecs;            // one (R,t) per view, also solved for
double rms = cv::calibrateCamera(objectPointsList, imagePointsList, imageSize,
                                 K, distCoeffs, rvecs, tvecs);
// rms is the global RMS reprojection error in pixels. Want it well under 1.0.
```

Python equivalent, only to show the identical call shape:

```python
import numpy as np, cv2
ok, rvec, tvec = cv2.solvePnP(objectPoints, imagePoints, K, distCoeffs,
                              flags=cv2.SOLVEPNP_ITERATIVE)
R, _ = cv2.Rodrigues(rvec)                    # r -> R
reproj, _ = cv2.projectPoints(objectPoints, rvec, tvec, K, distCoeffs)
rms = np.sqrt(np.mean(np.sum((imagePoints - reproj.reshape(-1, 2))**2, axis=1)))
```

## Connections

- [[09_camera-pinhole]]: the noise-free projective core $s\,p = K[R\mid t]X$. Calibration adds distortion and fits the constants; this topic is the pinhole model made real.
- [[10_nonlinear-least-squares]]: both calibration and PnP are instances of $\min_x \sum r_i(x)^2$ solved by Gauss-Newton / Levenberg-Marquardt; the residual is reprojection error and the Jacobians come from `projectPoints`.
- [[12_registration-icp]]: PnP supplies the initial camera-to-object pose that ICP refines against dense geometry; ICP minimizes a 3D point/plane residual where PnP minimized a 2D pixel residual.
- [[06_se3-rigid-transforms]]: $[R\mid t]$ is an element of $SE(3)$; the Rodrigues vector is the $\mathfrak{so}(3)$ tangent used to parameterize $R$ on its manifold.
- [[02_svd-and-least-squares]]: Zhang's calibration initializer and EPnP both reduce to SVD/eigen problems before the nonlinear refine.

## Pitfalls and failure modes

- **Frame and convention mismatch.** OpenCV gives world-to-camera $(R, t)$ with $X_c = R X_w + t$. To get the camera's position in world coordinates you need $-R^\top t$, not $t$. Mixing these up flips your pipeline inside out. Also confirm the $+Z$-forward, $+Y$-down convention before feeding poses to a robot or CAD frame that uses $+Z$-up.
- **Distortion coefficient order.** The vector is $(k_1, k_2, p_1, p_2, k_3)$. Swapping $p_1/p_2$ or putting $k_3$ before the tangentials silently corrupts everything. Pass the array OpenCV produced, do not retype it.
- **Degenerate / coplanar geometry.** PnP from coplanar points has a near-mirror ambiguity (two poses with similar reprojection error). Use `SOLVEPNP_IPPE` for planar targets and check both returned solutions, or add non-coplanar points.
- **Ill-conditioning from poor view diversity.** If all calibration views are at the same distance, same angle, or fill only the image center, $f_x/f_y$ trade off against $t_z$ and distortion is unobservable at the edges. Vary depth, tilt, and cover the corners. A low RMS on a bad set is overfitting, not accuracy.
- **Over-parameterization.** Turning on $k_4, k_5, k_6$ (rational model) with too few points lets the optimizer absorb noise into distortion. Use the simplest model that drives RMS down; fix high-order $k$ to zero otherwise.
- **Outlier correspondences.** A single mismatched point pulls the least-squares pose badly because the residual is squared. Use `solvePnPRansac` to reject outliers, then `solvePnPRefineLM` on the inliers.
- **Non-commutativity of pose composition.** Chaining camera-to-world with world-to-robot is matrix multiplication in $SE(3)$; the order matters and rotations do not commute. Compose with explicit $4\times4$ matrices, not by adding Euler angles.
- **Units and scale.** PnP returns $t$ in the same units as your object points. Feed millimeter object points and you get millimeter translation. A 1000x error here is almost always a meters/mm mixup, the same bug class flagged across the stone fixtures.

## Practice

1. **Hand projection.** With the worked-example $K$ and identity pose, project $X = (-0.2, 0.1, 4.0)$. Then add $k_1 = -0.2$ (barrel) and report how many pixels the point shifts. Verify your radial factor by hand.
2. **Reprojection error sweep.** Take 4 known 3D points, project them with a chosen ground-truth pose to make synthetic detections, then perturb one detection by 5 px. Run `solvePnP` and print how the RMS and the recovered $t$ change. Confirm the squared residual makes one outlier dominate.
3. **Calibrate from images.** Print a checkerboard, capture 15+ varied views, run `findChessboardCorners` + `cornerSubPix` + `calibrateCamera`. Inspect the per-view reprojection error and identify which view is worst and why.
4. **Undistort and check straight lines.** Use your calibrated $K, \mathbf{d}$ with `cv::undistort` on an edge-rich image. Overlay a straight ruler in the scene; confirm it is curved before and straight after.
5. **Planar ambiguity.** Build a coplanar 4-point set, solve with `SOLVEPNP_IPPE`, and print both candidate poses with their reprojection errors. Explain which one is physically correct and how you would disambiguate.
6. **RANSAC vs plain LM.** Add 30% gross-outlier correspondences and compare `solvePnP` against `solvePnPRansac` followed by `solvePnPRefineLM`. Report inlier count and final pose error.

## References

- OpenCV, "Camera Calibration and 3D Reconstruction" (calib3d module), https://docs.opencv.org/4.x/d9/d0c/group__calib3d.html (projection chain, distortion model, `calibrateCamera`, `projectPoints`).
- OpenCV, "Perspective-n-Point (PnP) pose computation", https://docs.opencv.org/4.x/d5/d1f/calib3d_solvePnP.html (solvePnP flags, minimum points, RANSAC and LM refinement).
- OpenCV, "Camera Calibration" tutorial, https://docs.opencv.org/4.x/dc/dbb/tutorial_py_calibration.html (radial/tangential formulas, coefficient order, reprojection-error recipe).
- Z. Zhang, "A Flexible New Technique for Camera Calibration," IEEE TPAMI 22(11):1330-1334, 2000, https://www.microsoft.com/en-us/research/publication/a-flexible-new-technique-for-camera-calibration/ (planar-homography calibration initializer).
- V. Lepetit, F. Moreno-Noguer, P. Fua, "EPnP: An Accurate O(n) Solution to the PnP Problem," IJCV 81(2):155-166, 2009 (`SOLVEPNP_EPNP`).
- G. Terzakis, M. Lourakis, "A Consistently Fast and Globally Optimal Solution to the Perspective-n-Point Problem," ECCV 2020, https://www.ecva.net/papers/eccv_2020/papers_ECCV/papers/123460460.pdf (`SOLVEPNP_SQPNP`).
- T. Collins, A. Bartoli, "Infinitesimal Plane-Based Pose Estimation," IJCV 2014 (`SOLVEPNP_IPPE`, coplanar pose).
- K. Lynch, F. Park, *Modern Robotics: Mechanics, Planning, and Control*, Cambridge, 2017, https://hades.mech.northwestern.edu/index.php/Modern_Robotics (SE(3), exponential map, Rodrigues).
- CADFit, "CAD program synthesis," arXiv:2605.01171 and https://github.com/ghadinehme/CADFit (downstream consumer of calibrated camera poses for CAD-fitting).

## Flashcard seeds

- Q: Derive the normalized camera coordinates (step 1 of 3) by projecting the camera-frame point $X_c=(X',Y',Z')$ onto the $z=1$ plane. <br> <b>Formula:</b> a ray scales every coordinate by $1/Z'$, giving target $x'=\dfrac{X'}{Z'},\ y'=\dfrac{Y'}{Z'}$ <br> <i>Hint:</i> divide each of $X'$ and $Y'$ by the depth $Z'$. :: A: A ray from the center through $(X',Y',Z')$ meets $z=1$ where each coordinate is scaled by $1/Z'$ <br> divide out depth: $x'=X'\cdot\frac{1}{Z'}$, $y'=Y'\cdot\frac{1}{Z'}$ <br> $\boxed{x'=\dfrac{X'}{Z'},\quad y'=\dfrac{Y'}{Z'}}$ (dimensionless; distortion acts here).
- Q: Derive the distorted radial coords (step 2 of 3) from $(x',y')$ using only radial coefficients $k_1,k_2,k_3$ (denominator $=1$). <br> <b>Formula:</b> $r^2=x'^2+y'^2$, factor $L(r)=1+k_1r^2+k_2r^4+k_3r^6$, target $x''=x'\,L(r)$ <br> <i>Hint:</i> build $L(r)$ first, then multiply each normalized coord by it. :: A: Compute $r^2=x'^2+y'^2$ <br> build the even-power radial factor $L(r)=1+k_1 r^2+k_2 r^4+k_3 r^6$ <br> scale both coords: $x''=x'\,L(r)$, $y''=y'\,L(r)$ <br> $\boxed{x''=x'(1+k_1 r^2+k_2 r^4+k_3 r^6)}$ (tangential terms add separately).
- Q: Derive the pixel mapping (step 3 of 3) from the distorted coords $(x'',y'')$ by applying the intrinsic matrix $K$. <br> <b>Formula:</b> target $u=f_x x''+c_x,\ v=f_y y''+c_y$ <br> <i>Hint:</i> multiply each normalized coord by its focal length, then add the principal-point offset. :: A: $K$ scales normalized units to pixels by focal length and shifts by the principal point <br> $u=f_x x''+c_x$, $v=f_y y''+c_y$ <br> $\boxed{u=f_x x''+c_x,\quad v=f_y y''+c_y}$ — this composed map is $\pi(K,\mathbf{d},R,t;X)$.
- Q: Derive the full projective equation (step 1 of 2) from world$\to$camera $X_c=[R\mid t]\tilde X$ then pixels $=KX_c$, in homogeneous form. <br> <b>Formula:</b> target $\tilde p = K\,[R\mid t]\,[X\ Y\ Z\ 1]^\top$ <br> <i>Hint:</i> substitute $X_c=[R\mid t]\tilde X$ into $\tilde p=KX_c$. :: A: Stack the extrinsic: $X_c=[R\mid t]\,[X\ Y\ Z\ 1]^\top$ <br> the intrinsic maps that camera ray to homogeneous pixels: $\tilde p = K X_c$ <br> $\boxed{\tilde p = K\,[R\mid t]\,[X\ Y\ Z\ 1]^\top}$.
- Q: From $\tilde p = K[R\mid t]\tilde X$ derive why a free scale $s$ appears (step 2 of 2). <br> <b>Formula:</b> $\tilde p=(su,sv,s)$ with $s=Z'$; target $s\,[u\ v\ 1]^\top = K[R\mid t]\tilde X$ <br> <i>Hint:</i> a homogeneous pixel is defined only up to scale; pull that scale out as $s$. :: A: Pixels live in a projective plane; $\tilde p=(su,sv,s)$ is defined only up to the scale $s=Z'$ <br> write the equality with the scale explicit <br> $\boxed{s\,[u\ v\ 1]^\top = K\,[R\mid t]\,[X\ Y\ Z\ 1]^\top}$; recover $u,v$ by dividing the third row.
- Q: Derive the calibration objective (step 1 of 2): write the single residual for detected pixel $\hat u_{ij}$ minus its projection. <br> <b>Formula:</b> target $r_{ij}=\hat u_{ij}-\pi(K,\mathbf{d},R_i,t_i;X_j)$ <br> <i>Hint:</i> subtract the projected point $\pi(\cdots)$ from the detected pixel for corner $j$ in view $i$. :: A: For corner $j$ in view $i$ the error is the detected pixel minus the projection <br> $r_{ij}=\hat u_{ij}-\pi(K,\mathbf{d},R_i,t_i;X_j)$ — the reprojection error in pixels <br> $\boxed{r_{ij}=\hat u_{ij}-\pi(K,\mathbf{d},R_i,t_i;X_j)}$.
- Q: Derive the full calibration least-squares problem (step 2 of 2) from the per-corner residual $r_{ij}$. <br> <b>Formula:</b> target $\min_{K,\mathbf{d},\{R_i,t_i\}}\ \sum_{i=1}^{M}\sum_{j=1}^{N}\lVert r_{ij}\rVert^2$ <br> <i>Hint:</i> square each residual and double-sum over $N$ corners and $M$ views, minimizing over all unknowns. :: A: Square each residual and sum over all $N$ corners and $M$ views <br> minimize jointly over intrinsics, distortion, and one pose per view <br> $\boxed{\min_{K,\mathbf{d},\{R_i,t_i\}}\ \sum_{i=1}^{M}\sum_{j=1}^{N}\lVert \hat u_{ij}-\pi(K,\mathbf{d},R_i,t_i;X_j)\rVert^2}$.
- Q: Derive the PnP objective from the calibration objective when $K,\mathbf{d}$ are already known. <br> <b>Formula:</b> target $\min_{R,t}\ \sum_{j=1}^{n}\lVert \hat u_j-\pi(K,\mathbf{d},R,t;X_j)\rVert^2$ <br> <i>Hint:</i> treat $K,\mathbf{d}$ as fixed constants, use one view, and drop the view index. :: A: Fix $K,\mathbf{d}$ as constants and use a single view, so only one pose remains <br> drop the view index and minimize over $(R,t)$ alone <br> $\boxed{\min_{R,t}\ \sum_{j=1}^{n}\lVert \hat u_j-\pi(K,\mathbf{d},R,t;X_j)\rVert^2}$ — 6 DoF instead of the full intrinsic set.
- Q: Derive the minimum point count for PnP from a degrees-of-freedom count. <br> <b>Formula:</b> pose DoF $=6$, equations per correspondence $=2$; solve $2n\ge 6$ <br> <i>Hint:</i> set the equation count $2n$ at least equal to the 6 pose unknowns. :: A: A pose has 6 unknowns; each 2D-3D correspondence supplies 2 scalar equations <br> require $2n\ge 6\Rightarrow n\ge 3$ (P3P, up to 4 discrete solutions) <br> $\boxed{n\ge 3\ \text{minimal; a 4th point disambiguates}}$.
- Q: Derive the camera position in world coordinates from the world-to-camera relation $X_c=RX_w+t$. <br> <b>Formula:</b> set $X_c=0$ and use $R^{-1}=R^\top$; target $C=-R^\top t$ <br> <i>Hint:</i> the camera center is the world point that maps to $X_c=0$ — set it to zero and solve for $X_w$. :: A: The camera center is the world point that maps to $X_c=0$ <br> $0=R X_w+t\Rightarrow R X_w=-t\Rightarrow X_w=-R^{-1}t$ <br> since $R\in SO(3)$, $R^{-1}=R^\top$ <br> $\boxed{C = -R^\top t}$ (not $t$).
- Q: Derive the Rodrigues rotation matrix from the exponential map $R=\exp(\theta K_\times)$, with $K_\times$ the skew matrix of the unit axis $\hat{\mathbf r}$. <br> <b>Formula:</b> use $K_\times^3=-K_\times$; target $R=I+\sin\theta\,K_\times+(1-\cos\theta)K_\times^2$ <br> <i>Hint:</i> expand the series and collapse odd/even powers with the identity $K_\times^3=-K_\times$. :: A: Expand $\exp(\theta K_\times)$ and use $K_\times^3=-K_\times$ to collapse the series <br> odd powers sum to $\sin\theta\,K_\times$, even powers to $(1-\cos\theta)K_\times^2$ <br> $\boxed{R=I+\sin\theta\,K_\times+(1-\cos\theta)K_\times^2}$, with $\theta=\lVert\mathbf r\rVert$.
- Q: Worked example (near, no distortion). With $f_x=f_y=800$, $c_x=320$, $c_y=240$ and identity pose, project $X=(0.1,-0.05,2.0)$ m. <br> <b>Formula:</b> $x'=X/Z$, $y'=Y/Z$; $u=f_x x'+c_x$, $v=f_y y'+c_y$ <br> <i>Hint:</i> divide by $Z=2.0$ first, then plug into the pixel equations. :: A: Normalize: $x'=0.1/2=0.05$, $y'=-0.05/2=-0.025$ <br> no distortion so $x''=x'$, $y''=y'$ <br> $u=800(0.05)+320=360$, $v=800(-0.025)+240=220$ <br> $\boxed{(u,v)=(360,220)}$.
- Q: Worked example (far / large depth, same $K$, identity pose). Project $X=(-0.2,0.1,4.0)$ m, no distortion. <br> <b>Formula:</b> $x'=X/Z$, $y'=Y/Z$; $u=800x'+320$, $v=800y'+240$ <br> <i>Hint:</i> divide by $Z=4.0$ first; expect a smaller off-center offset than at $Z=2$. :: A: Normalize: $x'=-0.2/4=-0.05$, $y'=0.1/4=0.025$ <br> $u=800(-0.05)+320=-40+320=280$ <br> $v=800(0.025)+240=20+240=260$ <br> $\boxed{(u,v)=(280,260)}$ — doubling $Z$ halves the off-center pixel offset.
- Q: Worked example (mm shop fixture, different scale). Same $K$, identity pose, project $X=(10,-5,200)$ mm. <br> <b>Formula:</b> $x'=X/Z$, $y'=Y/Z$; $u=800x'+320$, $v=800y'+240$ <br> <i>Hint:</i> compute the ratios $X/Z$ and $Y/Z$ first — normalization is unit-free. :: A: Normalization is scale-free: $x'=10/200=0.05$, $y'=-5/200=-0.025$ <br> identical normalized coords to the 2.0 m case <br> $u=800(0.05)+320=360$, $v=800(-0.025)+240=220$ <br> $\boxed{(u,v)=(360,220)}$ — pixels depend on the ratio, not the unit; but $t$ returns in mm.
- Q: Worked example (distortion regime). With $k_1=-0.2$ (barrel), $x'=0.05$, $y'=-0.025$, recompute $x''$. <br> <b>Formula:</b> $r^2=x'^2+y'^2$, $L=1+k_1r^2$, $x''=x'\,L$ <br> <i>Hint:</i> compute $r^2$ first, then the scalar factor $L$, then multiply. :: A: $r^2=0.05^2+0.025^2=0.003125$ <br> radial factor $L=1+(-0.2)(0.003125)=1-0.000625=0.999375$ <br> $x''=0.05\cdot0.999375=0.04996875$ <br> $\boxed{x''\approx0.0499688}$ — barrel pulls the point inward (factor $<1$).
- Q: Test yourself (close the loop). You projected $X$ to $(360,220)$; the detector found $(361.5,219.0)$. Compute the reprojection residual and its norm. <br> <b>Formula:</b> $r=\hat u-u_{\text{proj}}$, $\lVert r\rVert=\sqrt{r_x^2+r_y^2}$ <br> <i>Hint:</i> subtract componentwise first, then take the norm; a healthy single-point error is sub-pixel-ish. :: A: $r=(361.5-360,\ 219.0-220)=(1.5,-1.0)$ px <br> $\lVert r\rVert=\sqrt{1.5^2+1.0^2}=\sqrt{3.25}\approx1.80$ <br> both terms are pixel differences, so the norm is in px <br> $\boxed{\lVert r\rVert\approx1.80\ \text{px}}$ — exactly the per-point term the sum-of-squares minimizes.
- Q: Test yourself (limiting-case check). Confirm the radial factor $L(r)=1+k_1r^2+k_2r^4+k_3r^6$ is identity at the principal point, then apply $k_1=0.1$ at $r^2=0.003125$. <br> <b>Formula:</b> evaluate $L$ at $r=0$, then at $r^2=0.003125$; $u=800x''+320$ <br> <i>Hint:</i> a pass means $L=1$ exactly when $r=0$ (no on-axis distortion). :: A: At the center $r=0\Rightarrow L=1$, so $x''=x'$ (no distortion on-axis) — pass; this is why corner coverage matters <br> with $k_1=0.1$, $r^2=0.003125$: $L=1.0003125$, $x''=0.05\cdot1.0003125=0.05001563$, $u=800(0.05001563)+320=360.0125$ <br> $\boxed{\Delta u\approx+0.0125\ \text{px outward}}$ — small near center, grows with $r^2$.
- Q: Test yourself (sanity-check a calibration run). After `calibrateCamera` you get RMS $=0.4$ px on 15 varied views. Decide pass/fail and state the criterion. <br> <b>Formula:</b> RMS $=\sqrt{\frac{1}{NM}\sum\lVert r_{ij}\rVert^2}$; healthy when well under $1.0$ px <br> <i>Hint:</i> a pass is sub-pixel RMS on diverse views; suspect overfitting if a tiny RMS comes from too few or too-coplanar views. :: A: Criterion: sub-pixel RMS (well under $1.0$ px) on depth/tilt/corner-varied views indicates healthy calibration <br> $0.4<1.0$ and the views are varied, so this is a pass <br> $\boxed{\text{PASS: }0.4\ \text{px}<1.0\ \text{px}}$ — but a low RMS on a non-diverse set would be overfitting, not accuracy.
- CLOZE: Normalized coords $x'=X'/Z'$, $y'={{c1::Y'/Z'}}$ are computed {{c2::before}} distortion is applied. <br> <i>hint:</i> divide each camera-frame coordinate by the depth $Z'$.
- CLOZE: The OpenCV distortion vector is packed in the order $(k_1,k_2,{{c1::p_1,p_2}},k_3[,k_4,k_5,k_6])$, with the tangentials between $k_2$ and $k_3$.
- CLOZE: Calibration and PnP are both {{c1::nonlinear least squares}} problems whose residual is the {{c2::reprojection error}}, solved with Levenberg-Marquardt.
- CLOZE: A vision stage outputs a {{c1::transform}} plus a {{c2::confidence}} (RMS reprojection error), not a final answer; healthy RMS is well under {{c3::1.0}} px.
- CLOZE: To get the camera center in world coordinates use $C={{c1::-R^\top t}}$, not $t$. <br> <i>hint:</i> invert $X_c=RX_w+t$ at $X_c=0$ and use $R^{-1}=R^\top$.
- Q: When do you reach for `SOLVEPNP_IPPE` instead of the iterative default, and why? <br> <b>Formula:</b> coplanar PnP has a near-mirror ambiguity (two poses with similar reprojection error) <br> <i>Hint:</i> think about planar targets / square markers. :: A: Use it for coplanar point sets / planar targets (`SOLVEPNP_IPPE_SQUARE` for square markers) <br> plain coplanar PnP has a near-mirror ambiguity, so check both returned poses <br> $\boxed{\text{coplanar / planar targets} \Rightarrow \texttt{IPPE}}$.
- Q: Give a robust PnP recipe for outlier-contaminated correspondences, and say why plain LM fails. <br> <b>Formula:</b> the residual is squared, so one bad match dominates the sum <br> <i>Hint:</i> reject first, refine second. :: A: A single mismatch dominates because residuals are squared <br> reject outliers with `solvePnPRansac`, then refine on the inliers with `solvePnPRefineLM` <br> $\boxed{\texttt{solvePnPRansac} \to \texttt{solvePnPRefineLM}}$.
- Q: Why does the same PnP code give a 1000x-wrong translation across stone fixtures, and how do you fix it? <br> <b>Formula:</b> $t$ returns in the same units as the 3D object points $X_j$ <br> <i>Hint:</i> look for a meters/millimeters mixup ($1\,\text{m}=1000\,\text{mm}$). :: A: $t$ comes out in the units of the 3D object points; mixing meters and millimeters gives a 1000x scale error <br> keep object points and reported $t$ in one consistent unit system <br> $\boxed{\text{units}(t)=\text{units}(X_j)}$.
- Q: Which `solvePnP` flag is the general-purpose default, and what does it require to run? <br> <b>Formula:</b> needs an initial guess or a planar / $\ge 6$-point configuration <br> <i>Hint:</i> it is the LM-refinement flag. :: A: `SOLVEPNP_ITERATIVE` (Levenberg-Marquardt refinement) <br> it needs an initial guess or a planar / 6+ point configuration <br> $\boxed{\texttt{SOLVEPNP\_ITERATIVE}}$.
