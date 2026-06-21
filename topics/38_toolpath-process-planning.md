# Toolpath / process planning (CAM for stone)

Placement: [[23_nurbs-surfaces-fitting]], [[24_brep-boolean-occt]], [[18_discrete-curvature]] => **Toolpath / process planning** => (machine-ready motion)

## Why it matters for the stone scan->CAD->IFC / CADFit / robotic-fabrication workflow

This is the final translation: a validated CAD/scan surface becomes ordered, machine-ready tool motion that actually carves the stone. The choices here set the surface finish, the cycle time, and whether the tool gouges the work. It consumes the geometry from [[23_nurbs-surfaces-fitting]] and [[24_brep-boolean-occt]], emits tool poses (SE(3), [[05_transforms-se3]]) for [[37_trajectory-planning-collision]] and [[35_dynamics-control]], and is where "the model" finally becomes "the part."

## Intuition

Cover the target surface with a family of passes, like mowing a lawn in stripes. Between two adjacent stripes a thin ridge of uncut material is left: the **cusp** (scallop). Closer stripes (smaller stepover) leave a smaller cusp and a smoother surface, but take longer. The hero diagram shows two ball-end tool circles one stepover apart and the cusp between them, with the right triangle you solve to relate stepover, tool radius, and cusp height.

## The math (first principles)

**Tool positioning.** For a ball-end tool of radius $R$ touching a surface point $p$ with outward unit normal $n$, the **tool center** is offset along the normal:

$$ c = p + R\,n. $$

The toolpath is this offset surface; sweeping the contact point traces the path.

**Toolpath patterns.** Parallel / zig-zag (constant direction), contour / offset (constant-Z or boundary-parallel), and spiral. Pattern choice trades finish direction, retracts, and machining time.

**Stepover and cusp height.** Place two adjacent passes a stepover $s$ apart with tool radius $R$. From the right triangle in the diagram (half-stepover $s/2$, vertical leg $R-h$, hypotenuse $R$):

$$ R^2 = \left(\tfrac{s}{2}\right)^2 + (R-h)^2 \ \Longrightarrow\ h = R - \sqrt{R^2 - \left(\tfrac{s}{2}\right)^2}. $$

For $s \ll R$, a Taylor expansion gives the standard CAM approximation

$$ \boxed{\,h \approx \dfrac{s^2}{8R}\,}, \qquad\text{equivalently}\qquad s \approx \sqrt{8 R h}. $$

So to hit a finish target $h$, pick the stepover $s = \sqrt{8Rh}$.

**Material removal rate** (productivity): $\mathrm{MRR} = a_p\, a_e\, v_f$ with axial depth $a_p$, radial width $a_e$, feed $v_f$.

**Gouge avoidance.** A ball-end tool of radius $R$ cannot reach into a concave surface feature whose radius of curvature $\rho$ is smaller than $R$. Local gouging occurs when

$$ R > \rho_{\text{concave}}, $$

so the tool radius must be at most the smallest concave radius of curvature on the surface (the link to [[18_discrete-curvature]]). Global gouging (the shank or holder hitting a far part of the work) is checked by collision, link to [[37_trajectory-planning-collision]].

**Post-processing.** Each contact point plus the tool axis defines a tool pose in SE(3); the sequence of poses is post-processed to G-code or robot joint targets.

## Worked example

Botticino marble finish pass, ball-end tool radius $R = 6$ mm, target cusp height $h = 0.01$ mm. Stepover $s = \sqrt{8 R h} = \sqrt{8 \cdot 6 \cdot 0.01} = \sqrt{0.48} \approx 0.69$ mm. Check the exact formula: $h = 6 - \sqrt{36 - (0.345)^2} = 6 - \sqrt{36 - 0.119} = 6 - \sqrt{35.881} = 6 - 5.9901 = 0.0099$ mm, matching $0.01$ mm. The approximation $s^2/(8R)$ is excellent because $s \ll R$.

## Code

Idiomatic C++ with Eigen: choose a stepover for a target cusp, then generate parallel ball-end toolpath passes over a height-field surface, offsetting each contact point by the surface normal.

```cpp
#include <Eigen/Dense>
#include <vector>
#include <functional>
#include <cmath>
using Vec3 = Eigen::Vector3d;

// Stepover that yields a target cusp height for a ball-end tool of radius R.
double stepoverForCusp(double R, double h) { return std::sqrt(8.0 * R * h); }   // s = sqrt(8 R h)

// Parallel toolpath over a surface z = f(x,y) with normal n(x,y); ball-end radius R.
std::vector<std::vector<Vec3>> parallelToolpath(
    double x0, double x1, double y0, double y1, double R, double cusp,
    const std::function<double(double,double)>& f,
    const std::function<Vec3(double,double)>& n) {
    double s = stepoverForCusp(R, cusp);           // pass spacing
    double feed = s / 4.0;                          // along-pass sampling
    std::vector<std::vector<Vec3>> passes;
    for (double y = y0; y <= y1; y += s) {
        std::vector<Vec3> pass;
        for (double x = x0; x <= x1; x += feed) {
            Vec3 p(x, y, f(x, y));                  // contact point on surface
            Vec3 c = p + R * n(x, y);               // tool center c = p + R n
            pass.push_back(c);
        }
        passes.push_back(pass);
    }
    return passes;
}
```

## Connections

- [[23_nurbs-surfaces-fitting]] and [[24_brep-boolean-occt]] provide the surface the toolpath rides on.
- [[18_discrete-curvature]] gives the surface curvature that bounds the tool radius (gouge avoidance).
- [[05_transforms-se3]] is the tool pose (position + axis) emitted per contact point.
- [[37_trajectory-planning-collision]] checks global gouging and links passes with safe retracts.
- [[35_dynamics-control]] executes the path with feed/acceleration limits.

## Pitfalls and failure modes

- **Gouging on concave curvature.** If $R > \rho_{\text{concave}}$ the tool removes material it should not; pick a tool no larger than the smallest concave radius.
- **Finish vs time.** Cusp height scales as $s^2/(8R)$, so halving the stepover quarters the cusp but doubles the passes; choose deliberately.
- **Tool axis (3- vs 5-axis).** A fixed tool axis cannot reach undercuts or steep walls; stone often needs tilt (an extra orientation degree of freedom).
- **Units.** Shop units are millimetres; mixing with metre-scale site coordinates corrupts the stepover by 1000x (see [[26_ifc-bim]] units note).
- **Retracts and engagement.** Plunging at full feed or leaving the tool fully engaged at corners overloads it; ramp and control engagement.

## Practice

1. Derive $h \approx s^2/(8R)$ from $R^2=(s/2)^2+(R-h)^2$ via a Taylor expansion of the square root.
2. Implement `stepoverForCusp` and tabulate $s$ for $R=3,6,10$ mm at $h=0.005, 0.01, 0.02$ mm.
3. Generate a parallel toolpath over a synthetic bump surface and plot the passes; halve the stepover and compare pass count.
4. Flag gouge points: mark surface samples where the local concave radius is smaller than $R$.
5. Convert one pass of tool centers + a fixed tool axis into SE(3) poses (link to [[05_transforms-se3]]).

## References

- ISO 10303 / common CAM texts on scallop height and stepover (the $h \approx s^2/8R$ relation).
- Choi and Jerard, *Sculptured Surface Machining: Theory and Applications*.
- Marsh, *Applied Geometry for Computer Graphics and CAD* (offset surfaces).
- LinuxCNC / generic post-processor docs (toolpath to G-code).

## Flashcard seeds

- Q: Derive the exact cusp height $h$ in terms of stepover $s$ and tool radius $R$ from the diagram's right triangle. <br> <b>Formula:</b> start from $R^2 = (s/2)^2 + (R-h)^2$, target $h = R - \sqrt{R^2 - (s/2)^2}$. <br> <i>Hint:</i> isolate the $(R-h)^2$ term, then take the square root of both sides and solve for $h$. :: A: $R^2 = (s/2)^2 + (R-h)^2$ <br> $(R-h)^2 = R^2 - (s/2)^2$ <br> $R-h = \sqrt{R^2 - (s/2)^2}$ (positive root, since $h<R$) <br> $\boxed{h = R - \sqrt{R^2 - (s/2)^2}}$
- Q: Derive the small-stepover approximation $h \approx s^2/(8R)$ from the exact cusp height. <br> <b>Formula:</b> $h = R - \sqrt{R^2 - (s/2)^2}$, target $h \approx \dfrac{s^2}{8R}$. <br> <i>Hint:</i> factor $R$ out of the square root and Taylor-expand $\sqrt{1-x}\approx 1-\tfrac{x}{2}$ for small $x$. :: A: $h = R - R\sqrt{1 - \dfrac{(s/2)^2}{R^2}}$ <br> Let $x = \dfrac{(s/2)^2}{R^2} = \dfrac{s^2}{4R^2}$ (small) <br> $\sqrt{1-x} \approx 1 - \dfrac{x}{2}$, so $h \approx R - R\left(1 - \dfrac{x}{2}\right) = R\dfrac{x}{2}$ <br> $h \approx \dfrac{R}{2}\cdot\dfrac{s^2}{4R^2} = \boxed{\dfrac{s^2}{8R}}$
- Q: Derive the stepover $s$ that hits a target cusp height $h$. <br> <b>Formula:</b> start from $h \approx \dfrac{s^2}{8R}$, target $s \approx \sqrt{8Rh}$. <br> <i>Hint:</i> multiply both sides by $8R$ to clear the fraction, then take the square root. :: A: $h \approx \dfrac{s^2}{8R}$ <br> $8Rh \approx s^2$ <br> $\boxed{s \approx \sqrt{8 R h}}$
- Q: Worked example (mm scale): a finish pass uses ball-end radius $R=6$ mm with target cusp $h=0.01$ mm. Find the stepover $s$. <br> <b>Formula:</b> $s = \sqrt{8 R h}$. <br> <i>Hint:</i> substitute the numbers under the root and multiply $8\cdot 6\cdot 0.01$ first. :: A: $s = \sqrt{8\cdot 6\cdot 0.01}$ <br> $= \sqrt{0.48}$ <br> $\boxed{s \approx 0.69\text{ mm}}$
- Q: Worked example (coarser finish): ball-end radius $R=3$ mm, target cusp $h=0.02$ mm. Find the stepover $s$. <br> <b>Formula:</b> $s = \sqrt{8 R h}$. <br> <i>Hint:</i> compute $8\cdot 3\cdot 0.02$ inside the root, then take the square root. :: A: $s = \sqrt{8\cdot 3\cdot 0.02}$ <br> $= \sqrt{0.48}$ <br> $\boxed{s \approx 0.69\text{ mm}}$
- Q: Worked example (larger tool): ball-end radius $R=10$ mm, target cusp $h=0.005$ mm. Find the stepover $s$. <br> <b>Formula:</b> $s = \sqrt{8 R h}$. <br> <i>Hint:</i> compute $8\cdot 10\cdot 0.005$ inside the root first. :: A: $s = \sqrt{8\cdot 10\cdot 0.005}$ <br> $= \sqrt{0.40}$ <br> $\boxed{s \approx 0.63\text{ mm}}$
- Q: Worked example (predict the finish): with $R=6$ mm and stepover $s=0.69$ mm, what cusp height $h$ results? <br> <b>Formula:</b> $h \approx \dfrac{s^2}{8R}$. <br> <i>Hint:</i> square the stepover first, then divide by $8R = 48$. :: A: $s^2 = 0.69^2 = 0.476$ <br> $h \approx \dfrac{0.476}{8\cdot 6} = \dfrac{0.476}{48}$ <br> $\boxed{h \approx 0.0099\text{ mm}}$
- Q: Worked example (3-axis tool center): a contact point $p=(10,20,5)$ mm has outward unit normal $n=(0,0,1)$; the ball-end radius is $R=6$ mm. Find the tool-center point $c$. <br> <b>Formula:</b> $c = p + R\,n$. <br> <i>Hint:</i> scale the normal by $R$ first, then add it componentwise to $p$. :: A: $R\,n = 6\cdot(0,0,1) = (0,0,6)$ <br> $c = (10,20,5) + (0,0,6)$ <br> $\boxed{c = (10,20,11)\text{ mm}}$
- Q: Worked example (MRR): axial depth $a_p=2$ mm, radial width $a_e=3$ mm, feed $v_f=500$ mm/min. Find the material removal rate. <br> <b>Formula:</b> $\mathrm{MRR} = a_p\, a_e\, v_f$. <br> <i>Hint:</i> just multiply the three numbers together. :: A: $\mathrm{MRR} = 2\cdot 3\cdot 500$ <br> $\boxed{\mathrm{MRR} = 3000\text{ mm}^3/\text{min}}$
- Q: State the Pythagorean relation encoded by the cusp diagram's right triangle (legs $s/2$ and $R-h$, hypotenuse $R$). <br> <b>Formula:</b> $R^2 = (s/2)^2 + (R-h)^2$. <br> <i>Hint:</i> hypotenuse squared equals the sum of the two legs squared. :: A: The hypotenuse is the tool radius $R$; the legs are the half-stepover $s/2$ and the vertical $R-h$. <br> $\boxed{R^2 = (s/2)^2 + (R-h)^2}$
- Q: For a ball-end tool, where does the tool center sit relative to the contact point $p$ with unit normal $n$? <br> <b>Formula:</b> $c = p + R\,n$. <br> <i>Hint:</i> the offset is along the surface normal by exactly the tool radius. :: A: The toolpath is the surface offset outward along the normal by $R$. <br> $\boxed{c = p + R\,n}$
- Q: Closing test (finish target met?): for $R=6$ mm and chosen stepover $s=0.69$ mm, verify the cusp meets the $h=0.01$ mm target using the exact formula. <br> <b>Formula:</b> $h = R - \sqrt{R^2 - (s/2)^2}$. <br> <i>Hint:</i> a pass means the computed $h$ comes out at or below $0.01$ mm. :: A: $s/2 = 0.345$ <br> $h = 6 - \sqrt{36 - 0.345^2} = 6 - \sqrt{36 - 0.119} = 6 - \sqrt{35.881}$ <br> $= 6 - 5.9901 = 0.0099$ mm $\le 0.01$ mm <br> $\boxed{\text{PASS}}$
- Q: Closing test (gouge check): a ball-end tool of radius $R=6$ mm machines a concave fillet of radius of curvature $\rho=4$ mm. Does it gouge? <br> <b>Formula:</b> local gouge when $R > \rho_{\text{concave}}$. <br> <i>Hint:</i> a pass (no gouge) needs $R \le \rho_{\text{concave}}$. :: A: Compare: $R = 6 > \rho = 4$ <br> The tool is larger than the smallest concave radius, so it cannot reach in. <br> $\boxed{\text{GOUGES — fail; pick } R \le 4\text{ mm}}$
- Q: Write the local gouge condition for a ball-end tool. <br> <b>Formula:</b> $R > \rho_{\text{concave}}$. <br> <i>Hint:</i> name what each side is — tool radius vs smallest concave radius. :: A: Gouging occurs when the tool radius exceeds the smallest concave radius of curvature. <br> $\boxed{R > \rho_{\text{concave}}}$
- Q: Write the material removal rate in terms of axial depth, radial width, and feed. <br> <b>Formula:</b> $\mathrm{MRR} = a_p\, a_e\, v_f$. <br> <i>Hint:</i> it is just the product of the three cutting quantities. :: A: Axial depth $a_p$, radial width $a_e$, feed $v_f$. <br> $\boxed{\mathrm{MRR} = a_p\, a_e\, v_f}$
- Q: Which earlier topic supplies the curvature that bounds the tool radius for gouge avoidance? <br> <i>Hint:</i> it is the per-point bending of the surface. :: A: Discrete curvature ([[18_discrete-curvature]]) gives $\rho_{\text{concave}}$, the smallest concave radius the tool must respect.
- Q: Why might stone carving need a 5-axis (tilting) toolpath instead of fixed 3-axis? <br> <i>Hint:</i> think about features a straight-down tool axis cannot enter. :: A: A fixed tool axis cannot reach undercuts or steep walls; tilting adds an orientation degree of freedom to reach them.
- Q: What does each toolpath contact point plus the tool axis define, and what is it turned into? <br> <i>Hint:</i> position plus orientation is a rigid pose. :: A: A tool pose in SE(3) ([[05_transforms-se3]]); the sequence is post-processed to G-code or robot joint targets.
- Q: Units pitfall in toolpaths — what goes wrong mixing shop and site units? <br> <i>Hint:</i> shop is millimetres, site is metres. :: A: Mixing millimetre shop units with metre site coordinates scales the stepover by 1000×, corrupting the path.
- CLOZE: Halving the stepover {{c1::quarters}} the cusp height (since $h\propto s^2$) but {{c2::doubles}} the number of passes.
- CLOZE: The toolpath is the surface {{c1::offset}} by the tool radius along the {{c2::normal}}. <br> <i>hint:</i> $c = p + R\,n$.
- CLOZE: For a ball-end tool, local gouging happens when the tool radius $R$ {{c1::exceeds}} the smallest {{c2::concave}} radius of curvature.
