# The pinhole camera model and projection

Prereqs: [[09_transforms-se3]] => this topic => unlocks [[11_calibration-pnp]].

## Why it matters for the stone scan->CAD->IFC / CADFit workflow

Every photogrammetry scan of a stone starts as a set of 2D pixels, and the pinhole model is the exact rule that says which 3D ray each pixel came from. When you triangulate a dense point cloud from multiple photos, the reconstruction engine inverts this model: it knows $K$ and $[R\mid t]$ for each shot and intersects rays to recover $P_w$. If you want to overlay a CAD edge or an IFC face back onto a source photo for visual validation, you forward-project it through $s\,p = K[R\mid t]P_w$ and check that the projected outline lands on the real stone in the image. The same matrix lets you back-project a hand-clicked silhouette point in an image to a 3D ray that you can intersect with the OCCT B-Rep, which is the manual fitting fallback when a CADFit-style synthesizer needs a human-seeded constraint. Get the frame conventions wrong here and your reprojected model is shifted, mirrored, or scaled, and no amount of downstream fitting will fix it.

## Intuition

Hold a shoebox with a tiny needle hole in one end and a sheet of tracing paper at the other. Light from the world passes through the single hole and paints an upside-down image on the paper. That hole is the camera center. Everything the camera sees is a bundle of straight lines (rays) through that one point.

A pinhole camera throws away one number per ray. A 3D point and the whole ray behind it land on the same pixel, so depth along the ray is lost. Projection is the act of dividing out that depth. The mental picture is a similar-triangles diagram: a point at depth $Z$ and height $Y$ projects to image height $f\,Y/Z$, because the small triangle from the hole to the image plane is similar to the big triangle from the hole to the point. Calibration later turns those metric heights into integer pixel coordinates. The whole model is just that one similar-triangles fact, dressed in matrices so it composes cleanly with rigid motion.

## The math (first principles)

### Frames and conventions

We use the OpenCV / computer-vision convention. The camera frame has $+Z$ pointing forward out of the lens along the optical axis, $+X$ to the right, and $+Y$ down. This is a right-handed frame. Pixel coordinates $(u,v)$ have origin at the top-left, $u$ increasing right, $v$ increasing down, which is why $+Y$ points down: it keeps image and camera axes aligned. Points in front of the camera have $Z > 0$.

Verified against the OpenCV `calib3d` documentation. A different library (for example a graphics/OpenGL convention with $-Z$ forward and $+Y$ up) negates axes and silently flips your image. State your frame before writing any code.

Units: $X,Y,Z$ and $t$ are in scene units (meters for site-scale stone work). Focal lengths $f_x,f_y$ and principal point $c_x,c_y$ are in pixels. The world point is $P_w = [X_w, Y_w, Z_w, 1]^T$ in homogeneous coordinates.

### The two stages: extrinsics then intrinsics

Stage 1, extrinsics. Move the world point into the camera frame with a rigid transform from $SE(3)$:

$$
\begin{bmatrix} x \\ y \\ z \end{bmatrix} = R \begin{bmatrix} X_w \\ Y_w \\ Z_w \end{bmatrix} + t,
\qquad R \in SO(3),\; t \in \mathbb{R}^3.
$$

Here $R$ and $t$ are the extrinsic parameters. $R$ is the rotation of world axes into camera axes and $t$ is the world origin expressed in the camera frame. A common bug source: $t$ is **not** the camera position in the world. The camera center in world coordinates is $C = -R^{T} t$. Equivalently $[R\mid t]$ maps world to camera, and its inverse maps camera to world.

Stage 2, the pinhole division. Project onto the normalized image plane at $z = 1$ by dividing out depth:

$$
x' = \frac{x}{z}, \qquad y' = \frac{y}{z}.
$$

These are the **normalized image coordinates**: dimensionless, independent of any sensor, the pure geometry of the ray. The division is where depth is destroyed.

Stage 3, intrinsics. Convert normalized coordinates to pixels with focal lengths and principal point:

$$
u = f_x\,x' + c_x, \qquad v = f_y\,y' + c_y.
$$

$f_x = f \cdot m_x$ and $f_y = f \cdot m_y$, where $f$ is focal length in scene units and $m_x, m_y$ are pixels per unit along each axis. Two focal lengths allow non-square pixels. $(c_x, c_y)$ is the principal point, where the optical axis pierces the image, near but not exactly the image center.

### Packing it into matrices

Collect the intrinsics into the calibration matrix $K$ and write the whole chain in homogeneous form. This is the equation OpenCV documents:

$$
s \begin{bmatrix} u \\ v \\ 1 \end{bmatrix}
=
\underbrace{\begin{bmatrix} f_x & 0 & c_x \\ 0 & f_y & c_y \\ 0 & 0 & 1 \end{bmatrix}}_{K}
\underbrace{\begin{bmatrix} r_{11} & r_{12} & r_{13} & t_1 \\ r_{21} & r_{22} & r_{23} & t_2 \\ r_{31} & r_{32} & r_{33} & t_3 \end{bmatrix}}_{[R\,\mid\,t]}
\begin{bmatrix} X_w \\ Y_w \\ Z_w \\ 1 \end{bmatrix}.
$$

The product $P = K[R\mid t]$ is the $3\times 4$ **projection matrix**. The scalar $s$ on the left is the **homogeneous scale factor**, and it is exactly the camera-frame depth:

$$
s = z = r_{31}X_w + r_{32}Y_w + r_{33}Z_w + t_3.
$$

To recover actual pixels you divide the first two entries of the right-hand $3$-vector by the third entry. That division is the projective step; $s$ is what you divide by. You cannot recover $s$ from $(u,v)$ alone, which is the precise statement that monocular projection loses depth.

### The general form with skew

Hartley and Zisserman write the calibration matrix for a finite projective camera with a skew term:

$$
K = \begin{bmatrix} \alpha_x & s & x_0 \\ 0 & \alpha_y & y_0 \\ 0 & 0 & 1 \end{bmatrix},
$$

where $\alpha_x, \alpha_y$ are the focal lengths in pixels, $(x_0, y_0)$ is the principal point, and $s$ is the **skew** (note: a different $s$ from the scale factor above; skew is often written $\gamma$ to avoid the clash). Skew is nonzero only if the pixel axes are not perpendicular, which never happens on modern CCD/CMOS sensors. OpenCV's standard pinhole model fixes skew to zero, which is why its $K$ has a $0$ in the $(1,2)$ slot. Set $\gamma = 0$ unless you have a documented reason not to.

### What the model omits

This is the ideal pinhole. Real lenses bend rays, so radial and tangential distortion are applied to $(x', y')$ **before** the intrinsic step. Distortion lives between stage 2 and stage 3 and is covered with calibration in [[11_calibration-pnp]]. For now assume an undistorted (rectified) image.

## Worked example

Take a camera with $f_x = f_y = 800$ px, principal point $c_x = 320$, $c_y = 240$ (a $640\times480$ image), zero skew.

Place the camera at the world origin looking straight down $+Z$, so $R = I$ and $t = 0$. World point $P_w = (0.1,\ -0.2,\ 2.0)$ meters.

Stage 1, extrinsics ($R=I, t=0$):

$$
(x, y, z) = (0.1,\ -0.2,\ 2.0).
$$

Stage 2, normalized coordinates:

$$
x' = \frac{0.1}{2.0} = 0.05, \qquad y' = \frac{-0.2}{2.0} = -0.10.
$$

Stage 3, pixels:

$$
u = 800 \cdot 0.05 + 320 = 40 + 320 = 360,
$$
$$
v = 800 \cdot (-0.10) + 240 = -80 + 240 = 160.
$$

So $P_w$ images at pixel $(360, 160)$, and the scale factor is $s = z = 2.0$.

Sanity check via the full matrix. With $R=I, t=0$:

$$
K [R\mid t] P_w = K \begin{bmatrix} 0.1 \\ -0.2 \\ 2.0 \end{bmatrix}
= \begin{bmatrix} 800(0.1) + 320(2.0) \\ 800(-0.2) + 240(2.0) \\ 2.0 \end{bmatrix}
= \begin{bmatrix} 720 \\ 320 \\ 2.0 \end{bmatrix}.
$$

Divide by the third entry: $u = 720/2.0 = 360$, $v = 320/2.0 = 160$. Matches.

Now move the camera back by $1$ m along world $-Z$, i.e. the camera sits at $C = (0,0,-1)$. Then $t = -R\,C = -I\,(0,0,-1) = (0,0,1)$, so $z = 2.0 + 1 = 3.0$. Then $x' = 0.1/3.0 = 0.0333$, $u = 800(0.0333)+320 = 346.7$. The point moves toward the principal point as the camera retreats, exactly as a farther object shrinks toward center. Check $t \ne C$: here $C_z = -1$ but $t_3 = +1$. Confusing the two is the most common projection bug.

## Code

C++ with Eigen first. This is the minimal forward-projection kernel. Symbol map is in the comments.

```cpp
// pinhole.hpp  — ideal pinhole forward projection (no distortion).
// Build: g++ -I/path/to/eigen pinhole_demo.cpp -o demo
#pragma once
#include <Eigen/Dense>
#include <optional>

namespace pinhole {

// K: 3x3 calibration matrix [[fx,0,cx],[0,fy,cy],[0,0,1]].
inline Eigen::Matrix3d makeK(double fx, double fy, double cx, double cy,
                             double skew = 0.0) {
    Eigen::Matrix3d K = Eigen::Matrix3d::Zero();
    K(0,0) = fx;  K(0,1) = skew;  K(0,2) = cx;   // skew = gamma, normally 0
    K(1,1) = fy;  K(1,2) = cy;
    K(2,2) = 1.0;
    return K;
}

// Project one world point P_w through extrinsics (R,t) and intrinsics K.
// Returns pixel (u,v); std::nullopt if the point is at/behind the camera (z<=0).
//   stage 1:  Xc = R*Pw + t          (world -> camera, SE3)
//   stage 2:  x' = Xc.x/Xc.z, ...     (perspective divide, s = z)
//   stage 3:  [u,v] = K * [x',y',1]   (camera -> pixels)
inline std::optional<Eigen::Vector2d>
project(const Eigen::Matrix3d& K,
        const Eigen::Matrix3d& R,
        const Eigen::Vector3d& t,
        const Eigen::Vector3d& Pw) {
    Eigen::Vector3d Xc = R * Pw + t;        // camera-frame point (x,y,z)
    const double s = Xc.z();                 // homogeneous scale = depth
    if (s <= 1e-12) return std::nullopt;     // reject behind-camera points
    Eigen::Vector3d xn(Xc.x()/s, Xc.y()/s, 1.0);  // normalized coords (x',y',1)
    Eigen::Vector3d px = K * xn;             // px = [u,v,1]
    return Eigen::Vector2d(px.x(), px.y());  // already homogeneous-normalized
}

// Back-project a pixel to a unit ray direction in the camera frame.
// Depth is unknown, so this returns the ray, not a 3D point.
inline Eigen::Vector3d backprojectRay(const Eigen::Matrix3d& K,
                                      const Eigen::Vector2d& uv) {
    Eigen::Vector3d px(uv.x(), uv.y(), 1.0);
    Eigen::Vector3d dir = K.inverse() * px;  // x' = K^{-1} * [u,v,1]
    return dir.normalized();                  // ray through camera center
}

} // namespace pinhole
```

Driver that reproduces the worked example:

```cpp
#include "pinhole.hpp"
#include <iostream>

int main() {
    auto K = pinhole::makeK(800, 800, 320, 240);   // fx,fy,cx,cy
    Eigen::Matrix3d R = Eigen::Matrix3d::Identity();
    Eigen::Vector3d t(0, 0, 0);
    Eigen::Vector3d Pw(0.1, -0.2, 2.0);

    if (auto uv = pinhole::project(K, R, t, Pw))
        std::cout << "pixel = " << uv->transpose() << "\n";  // 360 160
    return 0;
}
```

The OpenCV idiom for the same operation, when you already work in that library, is `cv::projectPoints(objectPts, rvec, tvec, K, distCoeffs, imagePts)`, where `rvec` is the Rodrigues axis-angle vector of $R$. The `s` divide happens internally.

Short Python check (NumPy), useful to validate the C++ output:

```python
import numpy as np
K = np.array([[800, 0, 320], [0, 800, 240], [0, 0, 1.0]])
Rt = np.hstack([np.eye(3), np.zeros((3, 1))])   # [R|t]
Pw = np.array([0.1, -0.2, 2.0, 1.0])            # homogeneous world point
h = K @ Rt @ Pw                                 # h = [s*u, s*v, s]
print(h[:2] / h[2])                             # [360. 160.]
```

## Connections

- [[09_transforms-se3]]: the extrinsic block $[R\mid t]$ is exactly an $SE(3)$ rigid transform; everything about composition, inversion, and the camera center $C = -R^T t$ comes from there.
- [[11_calibration-pnp]]: this topic gives the forward map $P_w \to (u,v)$; calibration/PnP runs it backward to estimate $K$, $R$, $t$ from known correspondences, and adds lens distortion.
- [[12_epipolar-stereo]] (if present): two pinhole cameras define the epipolar geometry and the fundamental matrix, the basis of stereo triangulation that turns photos into the stone point cloud.
- [[13_least-squares-svd]] (if present): the projection matrix $P$ is recovered by a homogeneous least-squares / DLT solve, then decomposed into $K[R\mid t]$ by RQ decomposition.
- [[01_homogeneous-coordinates]] (if present): the perspective divide and the scale factor $s$ are pure projective-geometry mechanics.

## Pitfalls and failure modes

- **$t$ is not the camera position.** $t$ is the world origin in camera coordinates; the camera center is $C = -R^T t$. Plotting $t$ as the camera location flips and offsets everything.
- **Wrong axis convention.** OpenCV uses $+Z$ forward, $+Y$ down. OpenGL/graphics stacks use $-Z$ forward, $+Y$ up. Mixing them mirrors the image vertically and can put points behind the camera.
- **Points behind the camera.** When $z \le 0$ the perspective divide produces a finite but meaningless pixel (a point behind you appears as if in front). Always test $z > 0$ before dividing; the code above does.
- **Forgetting the divide.** $K[R\mid t]P_w$ gives $[su, sv, s]$, not $[u,v,1]$. You must divide the first two by the third. Reading off $su$ and $sv$ as pixels scales your image by the depth.
- **Skew confusion.** The scale factor $s$ and the skew $\gamma$ are different quantities that share a letter in some texts. Set skew to zero for real sensors.
- **Pixel non-squareness ignored.** Assuming $f_x = f_y$ when the sensor or resampling made pixels non-square stretches the projection. Keep both focal lengths.
- **Distortion applied in the wrong place.** Lens distortion acts on normalized coordinates $(x', y')$ before $K$, not on final pixels. Applying it after $K$ is wrong.
- **Ill-conditioning at depth.** As $z \to 0$ the divide blows up and small 3D errors become huge pixel errors; near-camera geometry is numerically fragile.
- **Non-commutativity of the stages.** You must rotate-then-translate-then-project in that order. $R$ and $t$ do not commute, and projection is not linear, so reordering breaks the result silently.

## Practice

1. By hand, project $P_w = (0.3, 0.0, 1.5)$ through $f_x=f_y=600$, $c_x=320$, $c_y=240$, $R=I$, $t=0$. Then move the camera to look from $C=(1,0,0)$ still facing $+Z$ (so $t = (-1,0,0)$) and recompute. Report both pixels and explain the horizontal shift.

2. Implement `backprojectRay` from the code, pick the pixel $(360,160)$ from the worked example, back-project it to a ray, then verify the ray passes through the original 3D point $(0.1,-0.2,2.0)$ up to scale. Print the ray and the normalized version of the point.

3. Write a loop that projects a $10\times10$ planar grid of 3D points lying on $z=2$ and writes the pixel coordinates to a CSV. Plot them. Confirm the grid stays a regular grid (an ideal pinhole maps planes to planes with no curvature).

4. Build $K$, $R$ (a $30^\circ$ rotation about the camera $Y$ axis), and $t$, then verify in code that $C = -R^T t$ reproduces the camera center you intended. Print $C$.

5. Given a known $K$ and a single image point with assumed depth $z=2.5$, recover the full 3D point $P_w$ in the camera frame, then express it in the world frame using $P_w^{world} = R^T(X_c - t)$. Confirm a round trip back through `project` returns the same pixel.

6. (Stretch) Take two synthetic cameras with known $K$, $R$, $t$, project one 3D point into both, then triangulate it back by intersecting the two back-projected rays in a least-squares sense. Report the residual. This is the core of turning stone photos into a point cloud.

## References

- OpenCV `calib3d` documentation, "Camera Calibration and 3D Reconstruction": the verified source for $s\,p = K[R\mid t]P_w$, the matrix $K$, and the $x'=x/z$, $u=f_x x'+c_x$ decomposition. https://docs.opencv.org/4.x/d9/d0c/group__calib3d.html (stable mirror: https://docs.opencv.org/2.4.13.7/modules/calib3d/doc/camera_calibration_and_3d_reconstruction.html)
- `cv::projectPoints` reference, OpenCV. https://docs.opencv.org/4.x/d9/d0c/group__calib3d.html#ga1019495a2c8d1743ed5cc23fa0daff8c
- Richard Hartley and Andrew Zisserman, *Multiple View Geometry in Computer Vision*, 2nd ed., Cambridge University Press, 2004. Chapter 6 (camera models), the finite projective camera and the calibration matrix $K$ with skew. https://www.robots.ox.ac.uk/~vgg/hzbook/
- Kevin Lynch and Frank Park, *Modern Robotics*, Cambridge University Press, 2017: $SE(3)$ background for the extrinsic transform. https://hades.mech.northwestern.edu/index.php/Modern_Robotics
- Eigen documentation (linear algebra used in the code). https://eigen.tuxfamily.org/dox/
- CADFit: Precise Mesh-to-CAD Program Generation with Hybrid Optimization, Ghadi Nehme et al., arXiv:2605.01171; code at https://github.com/ghadinehme/CADFit. Motivates projecting CAD candidates back into source images for validation. https://arxiv.org/abs/2605.01171

## Flashcard seeds

- Q: Derive the normalized image coordinate $x'$ from similar triangles (the small triangle from pinhole to the $z=1$ plane is similar to the big triangle from pinhole to the camera point $(x,y,z)$). <br> <b>Formula:</b> $\dfrac{x'}{1}=\dfrac{x}{z}\ \Rightarrow\ x'=\dfrac{x}{z}$ <br> <i>Hint:</i> equate (height across axis)/(depth along axis) for the two triangles, then solve for $x'$. :: A: Small triangle: height $x'$ at depth $1$; big triangle: height $x$ at depth $z$ <br> equate the ratios: $\dfrac{x'}{1}=\dfrac{x}{z}$ <br> $\boxed{x'=\dfrac{x}{z}}$ (and identically $y'=y/z$).
- Q: Starting from $x'=x/z$, derive the pixel column $u$ for focal length $f_x$ (px) and principal point $c_x$. <br> <b>Formula:</b> $u=f_x\,x'+c_x$, with $x'=x/z$ <br> <i>Hint:</i> first scale $x'$ by $f_x$ (pixels-per-unit), then add $c_x$ to shift the origin to top-left; finally substitute $x'=x/z$. :: A: Scale into pixels: $u=f_x\,x'+c_x$ <br> substitute $x'=x/z$ <br> $\boxed{u=f_x\dfrac{x}{z}+c_x}$ (and $v=f_y\dfrac{y}{z}+c_y$).
- Q: Derive why the homogeneous scale factor equals depth, starting from the pixel maps $u=f_x x/z+c_x$, $v=f_y y/z+c_y$. <br> <b>Formula:</b> $K[x\ y\ z]^T=[f_x x+c_x z,\ f_y y+c_y z,\ z]^T$, target $s=z$ <br> <i>Hint:</i> multiply $K$ by $[x\ y\ z]^T$ without dividing, then divide all three entries by the third and read off what cancels. :: A: $K[x\ y\ z]^T=[f_x x+c_x z,\ f_y y+c_y z,\ z]^T$ <br> divide all three by the third entry: $u=(f_x x+c_x z)/z=f_x x/z+c_x$, same for $v$ <br> the divisor is $z$, so $\boxed{s=z}$.
- Q: From the extrinsic map $x=RX_w+t$, write the world-to-camera step as one matrix (step 1 of 2). <br> <b>Formula:</b> $[x\ y\ z]^T=[R\mid t]\,[X_w\ Y_w\ Z_w\ 1]^T$ <br> <i>Hint:</i> stack $R$ (3x3) and $t$ (3x1) side by side into a 3x4 block, then apply it to the homogeneous point $P_w$. :: A: Concatenate so one matrix rotates-then-translates: $[R\mid t]$ is $3\times4$ <br> apply to $P_w=[X_w\ Y_w\ Z_w\ 1]^T$ <br> $\boxed{[x\ y\ z]^T=[R\mid t]\,P_w}$.
- Q: From $[x\ y\ z]^T=[R\mid t]P_w$ and $s=z$, derive the OpenCV projection equation (step 2 of 2). <br> <b>Formula:</b> $s[u\ v\ 1]^T=K[x\ y\ z]^T$, target $s[u\ v\ 1]^T=K[R\mid t]P_w$ <br> <i>Hint:</i> apply $K$ to the camera point, then substitute the extrinsic block in for $[x\ y\ z]^T$. :: A: Apply intrinsics with the scale on the left: $[su\ sv\ s]^T=K[x\ y\ z]^T$, $s=z$ <br> substitute $[x\ y\ z]^T=[R\mid t]P_w$ <br> $\boxed{s[u\ v\ 1]^T=K[R\mid t]P_w}$; recover $(u,v)$ by dividing the first two entries by the third.
- Q: Derive the camera center $C$ in world coordinates from $x=RX_w+t$. <br> <b>Formula:</b> $0=RC+t$, target $C=-R^T t$ <br> <i>Hint:</i> set the camera-frame point to $x=0$, solve for $C$, and use $R^{-1}=R^T$ since $R\in SO(3)$. :: A: The camera center maps to the camera origin, so set $0=R\,C+t$ <br> $R\,C=-t$ <br> left-multiply by $R^T$: $\boxed{C=-R^T t}$ (so $t\ne C$ in general).
- Q: Derive the back-projection ray direction $d$ for a pixel $[u\ v\ 1]^T$ (camera frame, $R=I,t=0$). <br> <b>Formula:</b> $s[u\ v\ 1]^T=K[x\ y\ z]^T$, target $d=K^{-1}[u\ v\ 1]^T$ <br> <i>Hint:</i> left-multiply both sides by $K^{-1}$; the unknown $s=z$ scales the result, so you get a direction, not a point. :: A: $K^{-1}[u\ v\ 1]^T=\tfrac{1}{s}[x\ y\ z]^T$ <br> the right side is the camera point up to the unknown scale $s=z$, i.e. a direction <br> $\boxed{d=K^{-1}[u\ v\ 1]^T}$ (then normalize); depth stays unknown.
- Q: Derive the image-plane offset $\Delta v$ of an object of true height $Y$ at depth $Z$, then take $Z\to\infty$. <br> <b>Formula:</b> $y'=Y/Z$, $\Delta v=f_y\,y'$ <br> <i>Hint:</i> get $y'$ from similar triangles first, multiply by $f_y$, then let $Z\to\infty$. :: A: Similar triangles: $y'=Y/Z$ on the $z=1$ plane <br> pixel offset from principal point: $\Delta v=f_y\,Y/Z$ <br> $\boxed{\Delta v=f_y\dfrac{Y}{Z}}$; as $Z\to\infty$, $\Delta v\to0$ — distant objects collapse toward the principal point.
- Q: Worked example (near, 3D, meters). $f_x=f_y=800$, $c_x=320$, $c_y=240$, $R=I$, $t=0$, $P_w=(0.1,-0.2,2.0)$ m. Find the pixel. <br> <b>Formula:</b> $u=f_x\dfrac{x}{z}+c_x$, $v=f_y\dfrac{y}{z}+c_y$ <br> <i>Hint:</i> with $R=I,t=0$ the camera point equals $P_w$; compute $x'=x/z$ and $y'=y/z$ first, then scale and shift. :: A: Extrinsics: $(x,y,z)=(0.1,-0.2,2.0)$ <br> normalize: $x'=0.1/2.0=0.05$, $y'=-0.2/2.0=-0.10$ <br> pixels: $u=800(0.05)+320=360$, $v=800(-0.10)+240=160$ <br> $\boxed{(u,v)=(360,160),\ s=z=2.0}$.
- Q: Worked example (camera retreat, regime change). Same $f=800$ camera, moved back $1$ m so $C=(0,0,-1)$, giving $t=(0,0,1)$. Re-project $P_w=(0.1,-0.2,2.0)$. <br> <b>Formula:</b> $z=Z_w+t_3$, $u=f_x\dfrac{x}{z}+c_x$ <br> <i>Hint:</i> recompute the depth $z=2.0+1$ first, then redo $x'=x/z$ and $u$; compare $u$ to $c_x=320$. :: A: New depth $z=2.0+1=3.0$ <br> $x'=0.1/3.0=0.0333$, so $u=800(0.0333)+320=346.7$ <br> $\boxed{u\approx346.7}$ — the point slides toward $c_x=320$ as the camera retreats, just as a farther object shrinks toward center.
- Q: Worked example (far / site-scale, large numbers). $f_x=f_y=1200$ px, $c_x=960$, $c_y=540$, $R=I$, $t=0$. Stone edge at $P_w=(2.0,-1.0,40.0)$ m. Find the pixel. <br> <b>Formula:</b> $u=f_x\dfrac{x}{z}+c_x$, $v=f_y\dfrac{y}{z}+c_y$ <br> <i>Hint:</i> divide by the large depth $z=40$ first to get small $x',y'$, then scale and shift. :: A: $x'=2.0/40=0.05$, $y'=-1.0/40=-0.025$ <br> $u=1200(0.05)+960=60+960=1020$ <br> $v=1200(-0.025)+540=-30+540=510$ <br> $\boxed{(u,v)=(1020,510),\ s=40.0}$.
- Q: Worked example (mm shop-scale, non-square pixels). $f_x=900$, $f_y=850$ px, $c_x=640$, $c_y=512$, $R=I$, $t=0$. Drill-hole edge at $P_w=(15,-8,300)$ mm. Find the pixel. <br> <b>Formula:</b> $u=f_x\dfrac{x}{z}+c_x$, $v=f_y\dfrac{y}{z}+c_y$ <br> <i>Hint:</i> compute $x'=x/z$ and $y'=y/z$ once, then use the two different focal lengths $f_x\ne f_y$ on each axis. :: A: $x'=15/300=0.05$, $y'=-8/300=-0.02\overline{6}$ <br> $u=900(0.05)+640=45+640=685$ <br> $v=850(-0.0267)+512=-22.7+512\approx489.3$ <br> $\boxed{(u,v)\approx(685,489.3)}$; $f_x\ne f_y$ stretches the axes differently.
- Q: Closing test. Project $P_w=(0.1,-0.2,2.0)$ with the $f=800$ camera by the staged formula, then by the full matrix, and confirm they match. <br> <b>Formula:</b> $K[R\mid t]P_w=[su,\,sv,\,s]^T$, then divide by $s$ <br> <i>Hint:</i> a pass means the matrix's $[su,sv,s]^T$ divided by its third entry equals the staged $(360,160)$. :: A: Staged: $(u,v)=(360,160)$ <br> matrix: $K[0.1,-0.2,2.0]^T=[800(0.1)+320(2.0),\,800(-0.2)+240(2.0),\,2.0]^T=[720,320,2.0]^T$ <br> divide by third: $720/2.0=360$, $320/2.0=160$ <br> $\boxed{(360,160)\ \text{matches}}$.
- Q: Closing test (round trip via back-projection). Back-project pixel $(360,160)$ with the $f=800$ camera and confirm the ray passes through $P_w=(0.1,-0.2,2.0)$ up to scale. <br> <b>Formula:</b> $d=K^{-1}[u\ v\ 1]^T$; compare to $P_w/z$ <br> <i>Hint:</i> a pass means $d$ equals $P_w$ divided by its own $z$ (both reduced to $z=1$). :: A: $d=K^{-1}[360,160,1]^T=[(360-320)/800,\,(160-240)/800,\,1]^T=[0.05,-0.10,1]^T$ <br> scale the point to $z=1$: $(0.1,-0.2,2.0)/2.0=[0.05,-0.10,1]^T$ <br> $\boxed{d\parallel P_w/z\ \checkmark}$ — same ray, depth ambiguous.
- Q: Closing test (units + limiting case) for $\Delta v=f_y\,Y/Z$. Check dimensions and the $Y=0$ case. <br> <b>Formula:</b> $\Delta v=f_y\dfrac{Y}{Z}$ <br> <i>Hint:</i> a pass means units reduce to pixels and $Y=0$ gives $\Delta v=0$ (on-axis lands on the principal point). :: A: Units: $[\text{px}]\cdot[\text{m}]/[\text{m}]=[\text{px}]$, a pure pixel offset $\checkmark$ <br> $Y=0$: $\Delta v=0$, on-axis point images at the principal point $\checkmark$ <br> $\boxed{\text{dimensionally consistent; on-axis case exact}}$.
- Q: Define the projection matrix $P$ and give its shape. <br> <b>Formula:</b> $P=K[R\mid t]$ <br> <i>Hint:</i> multiply the $3\times3$ intrinsics by the $3\times4$ extrinsics and read off the shape. :: A: $P=K[R\mid t]$ — $3\times3$ intrinsics times $3\times4$ extrinsics <br> $\boxed{P=K[R\mid t]\in\mathbb{R}^{3\times4}}$.
- Q: When should you use back-projection instead of forward projection? <br> <i>Hint:</i> think about what is known (pixel vs 3D point) and what the perspective divide destroyed. :: A: Use it when you have a pixel and want the 3D ray it came from <br> because depth $s=z$ was destroyed by the divide, you recover a ray, not a single point.
- Q: Pitfall: why is $t$ not the camera position in the world? <br> <b>Formula:</b> $C=-R^T t$ <br> <i>Hint:</i> recall what frame $t$ lives in versus where the camera center $C$ lives. :: A: $t$ is the world origin expressed in the camera frame <br> the camera center is $C=-R^T t$; plotting $t$ as the camera location flips and offsets everything.
- Q: Pitfall: what must you check before the perspective divide, and why? <br> <b>Formula:</b> require $z>0$ <br> <i>Hint:</i> consider what the divide $x/z$ does when $z\le0$ for a point behind the lens. :: A: Test $z>0$ (point in front of the camera) <br> for $z\le0$ the divide yields a finite but meaningless pixel, placing a behind-camera point as if it were in front.
- CLOZE: The pinhole model is one similar-triangles fact: a point at depth $Z$, height $Y$ images to height {{c1::$f\,Y/Z$}}, dressed in matrices to compose with rigid motion.
- CLOZE: The extrinsic transform maps a world point to the camera frame by {{c1::$x = R X_w + t$}}, and its inverse-derived camera center is {{c1::$C=-R^T t$}}.
- CLOZE: In OpenCV's camera convention, $+Z$ points {{c1::forward out of the lens}} and $+Y$ points {{c1::down}}; points in front satisfy {{c1::$z>0$}}.
- CLOZE: The {{c1::skew}} term in $K$, often written $\gamma$, is zero for real rectangular-pixel sensors and must not be confused with the scale factor $s=z$.
- CLOZE: Lens distortion is applied to the {{c1::normalized coordinates $(x',y')$}} before multiplying by $K$, never to the final pixels.
- CLOZE: The homogeneous product $K[R\mid t]P_w$ returns $[su,\,sv,\,s]$; you recover pixels by {{c1::dividing the first two entries by the third}}.
