# B-splines and NURBS curves

Prereqs: [[01_matrices]], [[02_calculus-jacobians]] => this => unlocks [[23_nurbs-surfaces-fitting]], [[24_brep-boolean-occt]].

## Why it matters for the stone scan->CAD->IFC / CADFit workflow

A raw stone scan is a triangle soup or point cloud. To become CAD it must become smooth parametric geometry, and in every industrial kernel that geometry is a NURBS. When you call `GeomAPI_PointsToBSpline` in OpenCASCADE to fit a profile curve through scan samples, the output is exactly a B-spline curve defined by the basis functions and knot vector below. When you later promote that curve to an IFC edge (`IfcBSplineCurveWithKnots`), the file literally stores the degree, knot vector, multiplicities, control points, and weights that this topic defines, so you cannot serialize correct IFC without understanding them. A CADFit-style program-synthesis layer that reconstructs a CAD program from a scan must emit these same parameters, so the synthesizer's output space is the space of knot vectors and control polygons. Getting parameterization and knot placement right is the difference between a clean fabricable edge and a wobbly curve that fails the Brep validator downstream.

## Intuition

A B-spline curve is a flexible wire bent by a small number of magnets called control points. The curve does not pass through the magnets; it is pulled toward them and stays inside the polygon they form. Each magnet only influences a short stretch of the wire, so moving one magnet near the start does nothing to the end. This is local control, and it is why splines are the right tool for editing scanned shapes: you can fix one bump without disturbing the rest.

The analogy: think of a string of weighted averages sliding along the curve. At any parameter $u$, the point on the curve is a blend of nearby control points, and the blend weights are the basis functions $N_{i,p}(u)$. The knot vector is the schedule that says when each magnet switches on and off as $u$ advances. NURBS add one more knob: a weight $w_i$ per magnet that makes it pull harder, which is what lets the wire trace a perfect circle, something no plain polynomial can do.

## The math (first principles)

### Conventions

- Parameter $u$ is dimensionless, increasing along the curve. Control points $P_i$ carry the model's length units (meters at site scale, millimeters at shop scale per the project unit policy).
- Degree is $p$, order is $p+1$. A cubic has $p=3$.
- Indices are 0-based: control points $P_0,\dots,P_n$ (so $n+1$ control points).
- Right-handed coordinate frame. Curves are frame-agnostic but the control points inherit the document frame.
- Tolerance for "is this knot equal to that knot" should be relative to the knot-vector span, not absolute, to stay scale-invariant.

### Knot vector

A knot vector is a non-decreasing sequence

$$U = \{u_0, u_1, \dots, u_m\}, \qquad u_0 \le u_1 \le \dots \le u_m.$$

The counts are linked by the fundamental relation

$$m = n + p + 1,$$

so the number of knots equals (number of control points) + (degree) + 1. A knot vector is **clamped** (also called open) when the first and last knots are repeated $p+1$ times. Clamping forces the curve to start at $P_0$ and end at $P_n$, which is what you almost always want for a CAD edge.

A knot value repeated $k$ times has **multiplicity** $k$. Multiplicity controls smoothness: at a knot of multiplicity $k$ the curve is

$$C^{\,p-k}\text{ continuous}.$$

So for a cubic ($p=3$): a simple interior knot ($k=1$) gives $C^2$, a double knot gives $C^1$, a triple knot gives $C^0$ (a corner), and multiplicity $p+1=4$ breaks the curve apart. This is the single most useful fact for stone work: you insert a triple knot to put a deliberate sharp edge where a chisel line should be, and you keep multiplicity 1 elsewhere for a smooth fabricable surface. (Verified against The NURBS Book continuity rule.)

### Cox-de Boor recursion (the basis functions)

The B-spline basis functions are defined recursively. Degree 0 is a box function:

$$N_{i,0}(u) = \begin{cases} 1 & u_i \le u < u_{i+1} \\ 0 & \text{otherwise} \end{cases}$$

Degree $p \ge 1$ blends two lower-degree functions:

$$N_{i,p}(u) = \frac{u - u_i}{u_{i+p} - u_i}\, N_{i,p-1}(u) \;+\; \frac{u_{i+p+1} - u}{u_{i+p+1} - u_{i+1}}\, N_{i+1,p-1}(u).$$

A fraction with a zero denominator is taken to be 0 (this happens at repeated knots). This is the Cox-de Boor recurrence, the definition every kernel implements. (Verified against The NURBS Book and standard course notes.)

Four properties fall out of this definition and you should memorize them:

1. **Non-negativity**: $N_{i,p}(u) \ge 0$ for all $u$.
2. **Partition of unity**: $\sum_{i} N_{i,p}(u) = 1$ on the valid domain $[u_p, u_{n+1}]$.
3. **Local support**: $N_{i,p}(u) = 0$ outside $[u_i, u_{i+p+1})$. So at most $p+1$ basis functions are nonzero at any $u$.
4. **Smoothness**: $N_{i,p}$ is $C^\infty$ between knots and $C^{p-k}$ across a knot of multiplicity $k$.

### B-spline curve

A (non-rational) B-spline curve is the basis-weighted sum of control points:

$$C(u) = \sum_{i=0}^{n} N_{i,p}(u)\, P_i, \qquad u \in [u_p, u_{n+1}].$$

Partition of unity makes each curve point an affine combination of control points. Two consequences:

- **Convex hull property**: $C(u)$ lies inside the convex hull of the control points whose basis functions are nonzero there. Because of local support, it lies in the hull of just $p+1$ nearby control points, not all of them.
- **Affine invariance**: to transform the curve, transform the control points. $T(C(u)) = \sum_i N_{i,p}(u)\, T(P_i)$ for any affine $T$. You never transform sampled points.

### NURBS curve (rational form)

Attach a positive weight $w_i$ to each control point. The Non-Uniform Rational B-Spline curve is

$$C(u) = \frac{\sum_{i=0}^{n} w_i\, P_i\, N_{i,p}(u)}{\sum_{i=0}^{n} w_i\, N_{i,p}(u)}.$$

Define the rational basis function

$$R_{i,p}(u) = \frac{w_i\, N_{i,p}(u)}{\sum_{j=0}^{n} w_j\, N_{j,p}(u)},$$

so that $C(u) = \sum_i R_{i,p}(u)\, P_i$, with $\sum_i R_{i,p}(u) = 1$ still holding. The $R_{i,p}$ keep non-negativity, local support, and partition of unity, so the convex hull property survives. (Verified against Piegl NURBS survey and the seed report equation.)

**Why rational**: a polynomial curve cannot trace a circle, but a rational one can. The trick is to lift to homogeneous coordinates. Build 4D control points $P_i^w = (w_i x_i,\, w_i y_i,\, w_i z_i,\, w_i)$, run the ordinary non-rational B-spline algorithm in 4D, then project back by dividing by the last coordinate. NURBS is "a polynomial B-spline in 4D, seen in perspective." This gives **projective invariance**: a projective (perspective) transform of the curve equals the same transform applied to the homogeneous control points.

### Exact conics

A rational quadratic ($p=2$) Bezier (a one-segment NURBS) with control points $P_0, P_1, P_2$ where $P_0P_1 = P_1P_2$ traces a circular arc when the middle weight is

$$w_1 = \cos(\theta),$$

where $\theta$ is the half-angle of the arc at $P_1$. A 90-degree arc (quarter circle) uses half-angle 45 degrees, so $w_1 = \cos 45^\circ = \tfrac{\sqrt 2}{2} \approx 0.7071$, with $w_0 = w_2 = 1$. Chaining such arcs gives an exact full circle, which is why CAD kernels store circles as NURBS. (Verified against the rational-Bezier conic construction.)

### Knot insertion (Boehm)

Inserting a knot $\bar u$ into span $[u_k, u_{k+1})$ raises the control-point count by one without changing the curve's shape (it is a pure refinement). The new control points $Q_i$ are corner cuts of the old polygon:

$$Q_i = (1 - \alpha_i)\, P_{i-1} + \alpha_i\, P_i, \qquad \alpha_i = \frac{\bar u - u_i}{u_{i+p} - u_i},$$

for $k - p + 1 \le i \le k$. Control points outside that range are copied unchanged. For NURBS, run the same formula on the homogeneous points $P_i^w$. (Verified against the Boehm single-insertion formula.) Knot insertion is the workhorse: it underlies curve splitting, de Boor evaluation, and refining a fit before a Brep boolean.

### Degree elevation

Degree elevation rewrites a degree-$p$ curve as an exact degree-$(p+1)$ curve. The shape is identical; the control-point count grows and interior knot multiplicities each increase by one (to preserve the same continuity). You elevate when a kernel operation needs matched degrees, for example sewing two edges into one face. The closed-form coefficients are in The NURBS Book; in practice you call the kernel.

## Worked example

Take a **quadratic** B-spline ($p = 2$) with four control points

$$P_0 = (0,0),\quad P_1 = (1,2),\quad P_2 = (3,2),\quad P_3 = (4,0),$$

so $n = 3$. The number of knots is $m = n + p + 1 = 3 + 2 + 1 = 6$, giving 7 knot values $u_0\dots u_6$. Use the clamped uniform vector

$$U = \{0, 0, 0,\ \tfrac12,\ 1, 1, 1\}.$$

First and last knots have multiplicity $p+1 = 3$, so the curve is clamped. The single interior knot $u_3 = \tfrac12$ has multiplicity 1, so the curve is $C^{p-k} = C^{1}$ there.

**Evaluate $C(u)$ at $u = \tfrac12$.** At $u=\tfrac12$ the nonzero quadratic basis functions are $N_{1,2}, N_{2,2}, N_{3,2}$ (local support: $p+1 = 3$ of them). Compute bottom-up.

Degree 0 at $u=\tfrac12$: the only active box is the one whose span contains $\tfrac12$. With $u_3 = \tfrac12$ and the convention $N_{i,0}=1$ on $[u_i, u_{i+1})$, the active box is $N_{3,0}$ on $[u_3,u_4)=[\tfrac12,1)$. So $N_{3,0}(\tfrac12)=1$, all other degree-0 are 0.

Degree 1 (only those reaching $N_{3,0}$ are nonzero):

$$N_{2,1}(\tfrac12) = \frac{u_4 - u}{u_4 - u_3}N_{3,0} = \frac{1 - \tfrac12}{1 - \tfrac12}\cdot 1 = 1,$$
$$N_{3,1}(\tfrac12) = \frac{u - u_3}{u_4 - u_3}N_{3,0} = \frac{\tfrac12 - \tfrac12}{1 - \tfrac12}\cdot 1 = 0.$$

Degree 2:

$$N_{1,2}(\tfrac12) = \frac{u - u_1}{u_3 - u_1}N_{2,1} = \frac{\tfrac12 - 0}{\tfrac12 - 0}\cdot 1 = 1\ \cdot\ \text{(but check the second term)}.$$

Carefully, $N_{1,2}(u) = \frac{u-u_1}{u_3-u_1}N_{1,1} + \frac{u_4-u}{u_4-u_2}N_{2,1}$. Here $N_{1,1}(\tfrac12)=0$ (it feeds $N_{2,0}$ which is 0), and $u_1=0,u_2=0,u_4=1$:

$$N_{1,2}(\tfrac12) = 0 + \frac{1 - \tfrac12}{1 - 0}\cdot 1 = \tfrac12.$$

$$N_{2,2}(\tfrac12) = \frac{u-u_2}{u_4-u_2}N_{2,1} + \frac{u_5-u}{u_5-u_3}N_{3,1} = \frac{\tfrac12-0}{1-0}\cdot 1 + \frac{1-\tfrac12}{1-\tfrac12}\cdot 0 = \tfrac12.$$

$$N_{3,2}(\tfrac12) = \frac{u-u_3}{u_5-u_3}N_{3,1} + \frac{u_6-u}{u_6-u_4}N_{4,1}.$$

$N_{3,1}(\tfrac12)=0$ and there is no $P_4$, so the active contributions are $N_{1,2}=\tfrac12$, $N_{2,2}=\tfrac12$, and we confirm partition of unity: $\tfrac12 + \tfrac12 = 1$. Good.

Curve point:

$$C(\tfrac12) = \tfrac12 P_1 + \tfrac12 P_2 = \tfrac12(1,2) + \tfrac12(3,2) = (2,\ 2).$$

That is the midpoint of the top edge of the control polygon, which matches intuition: at the interior knot of a symmetric polygon the curve sits centered and pulled up toward $P_1,P_2$. It is strictly inside the convex hull, as the convex-hull property promises.

**Sanity check at the ends**: $C(0) = P_0 = (0,0)$ and $C(1) = P_3 = (4,0)$ because the vector is clamped.

## Code

C++ first, using OpenCASCADE, the kernel you will actually ship with. This builds the worked-example curve, evaluates it, and inserts a knot. Symbol map is in comments.

```cpp
// build_bspline.cpp
// Link against TKMath TKG3d TKGeomBase TKBRep (OCCT 7.x).
#include <Geom_BSplineCurve.hxx>
#include <TColgp_Array1OfPnt.hxx>   // poles P_i
#include <TColStd_Array1OfReal.hxx> // knots, weights
#include <TColStd_Array1OfInteger.hxx> // multiplicities
#include <gp_Pnt.hxx>
#include <iostream>

int main() {
    const Standard_Integer p = 2;          // degree p

    // Control points P_0..P_3  (n = 3, so n+1 = 4 poles)
    TColgp_Array1OfPnt poles(1, 4);        // OCCT arrays are 1-based
    poles.SetValue(1, gp_Pnt(0, 0, 0));    // P_0
    poles.SetValue(2, gp_Pnt(1, 2, 0));    // P_1
    poles.SetValue(3, gp_Pnt(3, 2, 0));    // P_2
    poles.SetValue(4, gp_Pnt(4, 0, 0));    // P_3

    // Knot vector stored as DISTINCT knots + multiplicities.
    // Full vector {0,0,0, 1/2, 1,1,1} == knots {0, 1/2, 1} with mult {3,1,3}.
    TColStd_Array1OfReal    knots(1, 3);
    TColStd_Array1OfInteger mults(1, 3);
    knots.SetValue(1, 0.0);   mults.SetValue(1, p + 1); // clamp: mult = p+1
    knots.SetValue(2, 0.5);   mults.SetValue(2, 1);     // interior, C^1
    knots.SetValue(3, 1.0);   mults.SetValue(3, p + 1); // clamp

    // Check m = n + p + 1: sum(mults) = 3+1+3 = 7 = (3) + (2) + 1 + 1 = m+1 knot values. OK.
    Handle(Geom_BSplineCurve) curve =
        new Geom_BSplineCurve(poles, knots, mults, p);

    // Evaluate C(u) at u = 0.5  (expect (2, 2, 0))
    gp_Pnt pt;
    curve->D0(0.5, pt);  // D0 = position; D1 adds first derivative, D2 second
    std::cout << "C(0.5) = (" << pt.X() << ", " << pt.Y() << ", " << pt.Z() << ")\n";

    // Boehm knot insertion: add a knot at u = 0.75 (shape unchanged, +1 pole)
    curve->InsertKnot(0.75);
    std::cout << "after insert: " << curve->NbPoles() << " poles\n";

    // For a NURBS (rational) curve, use the weighted constructor:
    //   Geom_BSplineCurve(poles, weights, knots, mults, degree)
    // and supply w_i in a TColStd_Array1OfReal of the same length as poles.
    return 0;
}
```

A compact Python reference for the raw Cox-de Boor recursion, useful for unit-testing your understanding before trusting the kernel:

```python
# coxdeboor.py  -- verified against the worked example: C(0.5) == (2, 2)
def N(i, p, u, U):
    """B-spline basis N_{i,p}(u) over knot vector U (Cox-de Boor)."""
    if p == 0:
        # half-open spans; special-case the very last knot
        return 1.0 if (U[i] <= u < U[i+1]) or (u == U[-1] and U[i] <= u <= U[i+1]) else 0.0
    a_den = U[i+p]   - U[i]
    b_den = U[i+p+1] - U[i+1]
    a = ((u - U[i])   / a_den) * N(i,   p-1, u, U) if a_den > 0 else 0.0
    b = ((U[i+p+1]-u) / b_den) * N(i+1, p-1, u, U) if b_den > 0 else 0.0
    return a + b

def curve_point(u, p, P, U):
    """C(u) = sum_i N_{i,p}(u) P_i for a non-rational B-spline."""
    n = len(P) - 1
    return tuple(sum(N(i, p, u, U) * P[i][d] for i in range(n+1)) for d in range(len(P[0])))

if __name__ == "__main__":
    P = [(0,0), (1,2), (3,2), (4,0)]
    U = [0,0,0, 0.5, 1,1,1]
    print(curve_point(0.5, 2, P, U))  # -> (2.0, 2.0)
    print(curve_point(0.0, 2, P, U))  # -> (0.0, 0.0) clamped start
```

## Connections

- [[01_matrices]]: control points are stored and transformed as matrices; affine invariance is a matrix multiply on the poles.
- [[02_calculus-jacobians]]: curve derivatives $C'(u)$ (tangents) and the Jacobian of a fit drive both fairing and the Gauss-Newton step in least-squares fitting.
- [[23_nurbs-surfaces-fitting]]: the tensor-product surface $S(u,v) = \sum_i \sum_j R_{i,j}(u,v) P_{ij}$ is this topic squared; fitting a scan to a surface reuses the basis, knots, and parameterization machinery here.
- [[24_brep-boolean-occt]]: a Brep face is a trimmed NURBS surface bounded by NURBS edges; degree elevation and knot insertion are what OCCT runs internally to make two faces compatible before a boolean.
- The scan-to-CAD pipeline upstream depends on clean parameterization, which connects to remeshing and mesh-conditioning topics that produce the ordered point arrays a fitter needs.

## Pitfalls and failure modes

- **Wrong knot count**. If $m \ne n + p + 1$ the constructor throws or produces garbage. Always check `len(knots_full) == n + p + 2` (number of values), or in OCCT check that $\sum \text{mults} = n + p + 2$.
- **Unclamped when you wanted clamped**. Forgetting multiplicity $p+1$ at the ends means the curve does not reach $P_0$ and $P_n$, leaving gaps that fail edge sewing.
- **Bad parameterization on irregular samples**. Uniform parameterization on unevenly spaced scan points produces loops and overshoots. Use chord-length or, for noisy spacing, **centripetal** parameterization. OCCT's `GeomAPI_PointsToBSpline` exposes this; pick it for scan data.
- **Over-fitting noise**. Interpolating every scan point reproduces sensor noise as wiggles. Prefer approximation with a tolerance over exact interpolation, and use the fewest control points that meet tolerance.
- **Negative or zero weights**. A weight $w_i \le 0$ can push the denominator $\sum w_j N_{j,p}$ to zero and create a pole (division by zero) or a curve that escapes the convex hull. Keep all $w_i > 0$.
- **Ill-conditioned fits**. Too many control points relative to data points makes the fitting matrix near-singular. Symptom: huge oscillating control points. Fix by reducing control count or adding smoothing regularization.
- **Non-commutativity of refinements**. Degree elevation then knot insertion is not the same control polygon as the reverse order, even though the curve is identical. Compare curves geometrically, never by comparing raw control-point arrays.
- **Absolute tolerance at the wrong scale**. A $10^{-6}$ knot-equality tolerance is fine for a unit domain but wrong if your parameter range is huge. Keep tolerances relative to span, consistent with the project's scale-relative-epsilon rule.
- **Floating-point at repeated knots**. The zero-denominator convention in Cox-de Boor must be coded as an exact `den > 0` guard, not a near-zero compare, or you get NaNs at multiplicity.

## Practice

1. **Knot bookkeeping.** A cubic ($p=3$) clamped B-spline has 8 control points. How many knot values are in $U$? Write the multiplicities of the two end knots. (Answer check: $m = n+p+1 = 7+3+1 = 11$, so 12 knot values; end multiplicities are 4 each.)
2. **Continuity by hand.** Given a cubic with interior knots at $0.3$ (mult 1), $0.6$ (mult 2), $0.8$ (mult 3), state the continuity $C^{?}$ at each. Identify which one is a corner.
3. **Evaluate a basis.** Using the Python `N` function, plot $N_{0,2}, N_{1,2}, N_{2,2}, N_{3,2}$ over $U=\{0,0,0,0.5,1,1,1\}$ on a grid of $u$ and confirm they sum to 1 at every sample. Inspect the support of each.
4. **Exact quarter circle.** Build a rational quadratic NURBS with $P_0=(1,0)$, $P_1=(1,1)$, $P_2=(0,1)$ and middle weight $w_1=\tfrac{\sqrt2}{2}$ (ends weight 1). Sample 20 points and verify each lies on the unit circle to $10^{-9}$. Then set $w_1=1$ and watch it fail.
5. **Knot insertion is shape-preserving.** Take the worked-example curve, insert a knot at $u=0.75$ with Boehm's formula by hand (compute the $\alpha_i$ and the new $Q_i$), then evaluate the new curve at several $u$ and confirm the points are unchanged versus the original.
6. **Scan-fit comparison (OCCT).** Sample 30 points off a known arc, add small Gaussian noise, then fit with `GeomAPI_PointsToBSpline` once with uniform and once with centripetal parameterization. Compare max deviation and control-point count. Report which parameterization wins on noisy data.

## References

- Les Piegl and Wayne Tiller, *The NURBS Book*, 2nd ed., Springer 1997. Definitive source for Cox-de Boor, knot insertion, degree elevation, and continuity rules.
- Les Piegl, "On NURBS: A Survey," *IEEE Computer Graphics and Applications*, 11(1), 1991. https://ieeexplore.ieee.org/document/67687 — affine/projective invariance, rational form, why NURBS unify CAD.
- C.-K. Shene, CS3621 Geometric Modeling notes (Michigan Tech): B-spline basis https://pages.mtu.edu/~shene/COURSES/cs3621/NOTES/spline/B-spline/bspline-basis.html , knot insertion https://pages.mtu.edu/~shene/COURSES/cs3621/NOTES/spline/NURBS-knot-insert.html , rational-Bezier conics https://pages.mtu.edu/~shene/COURSES/cs3621/NOTES/spline/NURBS/RB-conics.html
- OpenCASCADE Technology reference manual, `Geom_BSplineCurve`: https://dev.opencascade.org/doc/refman/html/class_geom___b_spline_curve.html and `GeomAPI_PointsToBSpline` for fitting.
- OpenCASCADE Modeling Data user guide (knots, multiplicities, poles): https://dev.opencascade.org/doc/overview/html/occt_user_guides__modeling_data.html
- IfcOpenShell documentation and the IFC `IfcBSplineCurveWithKnots` / `IfcRationalBSplineCurveWithKnots` entities: https://ifcopenshell.org/
- buildingSMART IFC4 spec, B-spline curve entities: https://standards.buildingsmart.org/IFC/RELEASE/IFC4/ADD2_TC1/HTML/

## Flashcard seeds

- Q: Derive the knot-count relation $m=n+p+1$ (step 1 of 2): find the highest-indexed knot that the basis functions reference. <br> <b>Formula:</b> $C(u)=\sum_{i=0}^{n}N_{i,p}(u)P_i$, and Cox-de Boor makes $N_{i,p}$ reference knots up to $u_{i+p+1}$. <br> <i>Hint:</i> there are $n+1$ basis functions; plug the largest index $i=n$ into $u_{i+p+1}$. :: A: There are $n+1$ basis functions, indexed $i=0\dots n$ <br> $N_{i,p}$ reaches up to knot $u_{i+p+1}$ <br> the highest, $N_{n,p}$, references $u_{n+p+1}$ <br> $\boxed{u_{n+p+1}\text{ is the last knot needed}}$.
- Q: Finish the knot-count derivation (step 2 of 2): solve for $m$ and the number of knot values. <br> <b>Formula:</b> knots are indexed $u_0\dots u_m$ and the last one needed is $u_{n+p+1}$. <br> <i>Hint:</i> set the last index equal to $m$, then count values from $0$ to $m$ inclusive. :: A: The last index must equal $m$, so $m=n+p+1$ <br> knot values run $u_0\dots u_m$, a count of $m+1=n+p+2$ <br> $\boxed{m=n+p+1,\quad \#\text{knots}=n+p+2}$.
- Q: Derive the support interval of $N_{i,p}$ (step 1 of 2): find the interval on which $N_{i,1}$ can be nonzero. <br> <b>Formula:</b> $N_{i,0}\ne0$ only on $[u_i,u_{i+1})$; the recursion blends $N_{i,p-1}$ with $N_{i+1,p-1}$. <br> <i>Hint:</i> $N_{i,1}$ mixes $N_{i,0}$ and $N_{i+1,0}$ — take the union of their two boxes. :: A: $N_{i,1}$ blends $N_{i,0}$ on $[u_i,u_{i+1})$ and $N_{i+1,0}$ on $[u_{i+1},u_{i+2})$ <br> the union spans $[u_i,u_{i+2})$ <br> $\boxed{\operatorname{supp}N_{i,1}=[u_i,u_{i+2})}$.
- Q: Finish the support derivation (step 2 of 2): generalize $\operatorname{supp}N_{i,p}$ and state how many basis functions are nonzero at a generic $u$. <br> <b>Formula:</b> each degree step extends the right end by one knot: $[u_i,u_{i+2})\to[u_i,u_{i+3})\to\dots$ <br> <i>Hint:</i> start at width $2$ for $p=1$ and add one knot of reach per degree, $p$ times. :: A: After $p$ steps the support is $[u_i,u_{i+p+1})$ <br> a generic $u$ lies inside $p+1$ overlapping supports <br> $\boxed{\operatorname{supp}N_{i,p}=[u_i,u_{i+p+1}),\ \text{at most }p+1\text{ nonzero}}$.
- Q: Derive the convex-hull property of a B-spline from two basis properties. <br> <b>Formula:</b> $C(u)=\sum_i N_{i,p}(u)P_i$, with $N_{i,p}\ge0$ and $\sum_i N_{i,p}(u)=1$. <br> <i>Hint:</i> nonnegative weights that sum to $1$ define a convex combination — name where such a combination must lie. :: A: $C(u)$ is a weighted sum of the $P_i$ with weights $N_{i,p}\ge0$ summing to $1$ <br> that is by definition a convex combination of the $P_i$ <br> a convex combination lies in the convex hull; local support restricts it to the $p+1$ active points <br> $\boxed{C(u)\in\operatorname{conv}\{P_i:N_{i,p}(u)>0\}}$.
- Q: Derive affine invariance: show $T(C(u))=\sum_i N_{i,p}(u)\,T(P_i)$ for affine $T$. <br> <b>Formula:</b> $T(x)=Ax+b$, $C(u)=\sum_i N_{i,p}(u)P_i$, and $\sum_i N_{i,p}(u)=1$. <br> <i>Hint:</i> apply $T$, then use partition of unity to rewrite the lone $b$ as $\big(\sum_i N_{i,p}\big)b$ so it folds into the sum. :: A: $T(C(u))=A\big(\sum_i N_{i,p}P_i\big)+b=\sum_i N_{i,p}(AP_i)+b$ <br> rewrite $b=\big(\sum_i N_{i,p}\big)b=\sum_i N_{i,p}\,b$ via partition of unity <br> $=\sum_i N_{i,p}(AP_i+b)=\sum_i N_{i,p}\,T(P_i)$ <br> $\boxed{T(C(u))=\sum_i N_{i,p}(u)\,T(P_i)}$ — transform only the poles.
- Q: Derive the rational basis $R_{i,p}$ and show it still sums to $1$. <br> <b>Formula:</b> $C(u)=\dfrac{\sum_i w_iP_iN_{i,p}(u)}{\sum_j w_jN_{j,p}(u)}$; target $R_{i,p}=\dfrac{w_iN_{i,p}}{\sum_j w_jN_{j,p}}$. <br> <i>Hint:</i> pull the common denominator inside the sum so each $P_i$ gets its own coefficient, then add those coefficients. :: A: Move the common denominator inside: $C(u)=\sum_i\dfrac{w_iN_{i,p}}{\sum_j w_jN_{j,p}}P_i$ <br> define $R_{i,p}=\dfrac{w_iN_{i,p}}{\sum_j w_jN_{j,p}}$ <br> then $\sum_i R_{i,p}=\dfrac{\sum_i w_iN_{i,p}}{\sum_j w_jN_{j,p}}=1$ <br> $\boxed{C(u)=\sum_i R_{i,p}(u)P_i,\ \sum_i R_{i,p}=1}$.
- Q: Derive the NURBS projection formula from the homogeneous (4D) lift. <br> <b>Formula:</b> 4D poles $P_i^w=(w_ix_i,w_iy_i,w_iz_i,w_i)$, run $C^w(u)=\sum_i N_{i,p}(u)P_i^w$, then divide by the last coordinate. <br> <i>Hint:</i> write out the four components of $C^w(u)$, then perspective-divide the first three by the fourth. :: A: $C^w(u)=\big(\sum_i N_{i,p}w_ix_i,\dots,\sum_i N_{i,p}w_i\big)$ <br> divide first three by the 4th: $x=\dfrac{\sum_i N_{i,p}w_ix_i}{\sum_i N_{i,p}w_i}$ (same for $y,z$) <br> $\boxed{C(u)=\dfrac{\sum_i w_iP_iN_{i,p}(u)}{\sum_i w_iN_{i,p}(u)}}$.
- Q: Derive the continuity at an interior knot of multiplicity $k$ in a degree-$p$ spline. <br> <b>Formula:</b> $N_{i,p}$ is $C^\infty$ between knots and loses one derivative per coincident knot; target $C^{\,p-k}$. <br> <i>Hint:</i> a simple knot ($k=1$) gives $C^{p-1}$; subtract one more derivative for each extra coincidence, $k-1$ of them. :: A: Between knots $N_{i,p}$ is $C^\infty$; across a simple knot it drops to $C^{p-1}$ <br> raising multiplicity by one removes one more derivative <br> after $k$ coincidences: $p-1-(k-1)=p-k$ <br> $\boxed{C^{\,p-k}\text{ continuous at a knot of multiplicity }k}$.
- Q: Derive the middle weight that makes a rational quadratic Bezier an exact circular arc. <br> <b>Formula:</b> symmetric control legs $|P_0P_1|=|P_1P_2|$, end weights $1$, arc half-angle $\theta$ at $P_1$; target $w_1=\cos\theta$. <br> <i>Hint:</i> the conic is a circle exactly when the apex weight equals the cosine of the half-angle subtended at $P_1$. :: A: A rational quadratic with symmetric legs traces a conic <br> it is a circular arc exactly when $w_1$ equals the cosine of the half-angle at the apex $P_1$ <br> $\boxed{w_1=\cos\theta}$.
- Q: Evaluate the worked-example basis (step 1 of 2): find the nonzero degree-2 basis values at $u=\tfrac12$. <br> <b>Formula:</b> $U=\{0,0,0,\tfrac12,1,1,1\}$, $p=2$; $N_{i,2}=\frac{u-u_i}{u_{i+2}-u_i}N_{i,1}+\frac{u_{i+3}-u}{u_{i+3}-u_{i+1}}N_{i+1,1}$, with active $N_{2,1}=1$. <br> <i>Hint:</i> only the box $N_{3,0}$ is active, so $N_{2,1}=1,N_{3,1}=0$; substitute the knot values into $N_{1,2}$ and $N_{2,2}$. :: A: At $u=\tfrac12$ only $N_{3,0}$ is the active box, so $N_{2,1}=1,\ N_{3,1}=0$ <br> $N_{1,2}=\frac{u_4-u}{u_4-u_2}N_{2,1}=\frac{1-1/2}{1-0}=\tfrac12$ <br> $N_{2,2}=\frac{u-u_2}{u_4-u_2}N_{2,1}=\frac{1/2}{1}=\tfrac12$ <br> $\boxed{N_{1,2}(\tfrac12)=N_{2,2}(\tfrac12)=\tfrac12}$.
- Q: Finish the worked-example evaluation (step 2 of 2): compute $C(\tfrac12)$. <br> <b>Formula:</b> $C(\tfrac12)=N_{1,2}P_1+N_{2,2}P_2$ with $N_{1,2}=N_{2,2}=\tfrac12$, $P_1=(1,2)$, $P_2=(3,2)$. <br> <i>Hint:</i> it is just the average of $P_1$ and $P_2$ — add the coordinates and halve. :: A: $C(\tfrac12)=\tfrac12P_1+\tfrac12P_2$ <br> $=\tfrac12(1,2)+\tfrac12(3,2)=(2,2)$ <br> $\boxed{C(\tfrac12)=(2,2)}$ — midpoint of the top edge, inside the hull.
- Q: State the Boehm knot-insertion corner-cut coefficient $\alpha_i$. <br> <b>Formula:</b> $Q_i=(1-\alpha_i)P_{i-1}+\alpha_iP_i$ with $\alpha_i=\dfrac{\bar u-u_i}{u_{i+p}-u_i}$, for $k-p+1\le i\le k$ when $\bar u\in[u_k,u_{k+1})$. <br> <i>Hint:</i> $\alpha_i$ is just how far $\bar u$ sits along the interval $[u_i,u_{i+p}]$. :: A: $\alpha_i$ measures the position of $\bar u$ within $[u_i,u_{i+p}]$ <br> $\boxed{\alpha_i=\dfrac{\bar u-u_i}{u_{i+p}-u_i},\quad k-p+1\le i\le k}$.
- Q: Worked example (mm shop scale, small): a cubic ($p=3$) clamped B-spline fits a chamfer with $n+1=6$ control points in millimeters. Find the number of knot values and the two end multiplicities. <br> <b>Formula:</b> $m=n+p+1$, knot values $=m+1$, clamp multiplicity $=p+1$. <br> <i>Hint:</i> first get $n$ from $n+1=6$, then plug into $m$. :: A: $n=5$, so $m=n+p+1=5+3+1=9$ <br> knot values $=m+1=10$ <br> clamped ends each have multiplicity $p+1=4$ <br> $\boxed{10\text{ knots},\ \text{end mult}=4,4}$ (units don't affect counts).
- Q: Worked example (site scale, large, meters): a quartic ($p=4$) NURBS edge of a quarry-block boundary has $14$ control points in meters. Find the knot count and the clamp multiplicity. <br> <b>Formula:</b> $m=n+p+1$, knot values $=m+1$, clamp multiplicity $=p+1$. <br> <i>Hint:</i> $n=14-1=13$; substitute into $m$, then add one for the value count. :: A: $n=13$, $m=n+p+1=13+4+1=18$ <br> knot values $=m+1=19$ <br> clamp mult $=p+1=5$ <br> $\boxed{19\text{ knots},\ \text{end mult}=5}$.
- Q: Worked example (continuity regime): a cubic has interior knots at $0.3$ (mult 1), $0.6$ (mult 2), $0.8$ (mult 3). Give $C^{?}$ at each and name the corner. <br> <b>Formula:</b> continuity $=C^{\,p-k}$ with $p=3$. <br> <i>Hint:</i> compute $3-k$ for each $k$; the corner is the knot where continuity hits $C^0$. :: A: $p=3$, so continuity $=C^{3-k}$ <br> at $0.3$: $C^{3-1}=C^2$; at $0.6$: $C^{3-2}=C^1$; at $0.8$: $C^{3-3}=C^0$ <br> $\boxed{C^2,\,C^1,\,C^0;\ u=0.8\text{ is the corner}}$.
- Q: Worked example (exact quarter circle, 2D): build a rational quadratic arc $P_0=(1,0),P_1=(1,1),P_2=(0,1)$ for a $90^\circ$ arc; find $w_1$ and the end weights. <br> <b>Formula:</b> $w_1=\cos\theta$ with $\theta$ the arc half-angle; end weights $w_0=w_2=1$. <br> <i>Hint:</i> a $90^\circ$ arc has half-angle $\theta=45^\circ$; evaluate $\cos45^\circ$. :: A: half-angle $\theta=45^\circ$, so $w_1=\cos45^\circ=\tfrac{\sqrt2}{2}\approx0.7071$ <br> ends weight $w_0=w_2=1$ <br> $\boxed{w_1=\tfrac{\sqrt2}{2},\ w_0=w_2=1}$.
- Q: Closing test (verify the arc lies on the circle): sample the quarter-circle NURBS at $u=\tfrac12$ and check the radius. <br> <b>Formula:</b> midpoint $=\dfrac{\tfrac14 P_0+\tfrac{w_1}{2} P_1+\tfrac14 P_2}{\tfrac14+\tfrac{w_1}{2}+\tfrac14}$ with $w_1=\tfrac{\sqrt2}{2}$; pass if $x^2+y^2=1$. <br> <i>Hint:</i> by symmetry $x=y$; compute one component, then check $x^2+y^2$ equals $1$. :: A: numerator $=\tfrac14(1,0)+\tfrac{\sqrt2}{4}(1,1)+\tfrac14(0,1)=(\tfrac14+\tfrac{\sqrt2}{4},\ \tfrac{\sqrt2}{4}+\tfrac14)$ <br> denominator $=\tfrac12+\tfrac{\sqrt2}{4}$ <br> by symmetry $x=y=\dfrac{1+\sqrt2}{2+\sqrt2}\approx0.7071$ <br> $x^2+y^2=2(0.7071)^2=1.000$ <br> $\boxed{\text{on the unit circle, radius }1\ \checkmark}$.
- Q: Closing test (knot insertion is shape-preserving): insert $\bar u=0.75$ into the worked-example quadratic via Boehm and verify the shape is unchanged. <br> <b>Formula:</b> new poles $Q_i=(1-\alpha_i)P_{i-1}+\alpha_iP_i$ (convex, $\alpha_i\in[0,1]$); clamp keeps end multiplicity $p+1$. Pass if endpoints still match. <br> <i>Hint:</i> a clamped curve interpolates its ends, so re-evaluate $C(0)$ and $C(1)$ — they should be unchanged. :: A: insertion adds one pole and re-blends affected interior poles by convex combos ($\alpha_i\in[0,1]$) <br> end knots keep multiplicity $p+1=3$, so clamping is untouched <br> re-evaluate: $C(0)=P_0=(0,0)$, $C(1)=P_3=(4,0)$ as before <br> $\boxed{\text{geometry identical, }+1\text{ pole}\ \checkmark}$.
- Q: Closing test (IFC knot-count validation): you serialize an `IfcBSplineCurveWithKnots` with degree $p=3$, $7$ control points, and $10$ knot values. Does it validate? <br> <b>Formula:</b> required knot values $=n+p+2$. <br> <i>Hint:</i> $n=7-1=6$; compute $n+p+2$ and compare to $10$ — a pass means they're equal. :: A: required count $=n+p+2=6+3+2=11$ <br> you wrote $10\ne11$ <br> $\boxed{\text{invalid: need }11\text{ knot values, not }10}$.
- CLOZE: At a knot of multiplicity $k$ in a degree-$p$ B-spline, the curve is {{c1::$C^{p-k}$}} continuous.
- CLOZE: A clamped (open) knot vector repeats the first and last knots {{c1::$p+1$}} times so the curve interpolates its end control points.
- CLOZE: The B-spline basis functions satisfy partition of unity, meaning {{c1::$\sum_i N_{i,p}(u) = 1$}} on the valid domain.
- CLOZE: For irregularly spaced scan samples, prefer {{c1::centripetal}} (or chord-length) parameterization over uniform to avoid loops and overshoot. <br> <i>hint:</i> uniform parameterization on uneven spacing causes loops and overshoot.
- CLOZE: A NURBS curve equals a non-rational B-spline run in {{c1::homogeneous (4D)}} coordinates and then projected by dividing by the weight coordinate.
- Q: When-to-use: why insert a triple knot in a cubic for stone work? <br> <b>Formula:</b> continuity $=C^{\,p-k}$; here $p=3$, $k=3$. <br> <i>Hint:</i> compute $p-k$ to see the continuity a triple knot forces. :: A: a triple knot ($k=p=3$) gives $C^{3-3}=C^0$, a deliberate corner for a chisel line <br> keep multiplicity 1 elsewhere ($C^2$) for a smooth fabricable edge <br> $\boxed{\text{triple knot}\Rightarrow C^0\text{ corner}}$.
- Q: Pitfall: what happens if a NURBS weight is zero or negative? <br> <b>Formula:</b> denominator $D(u)=\sum_j w_jN_{j,p}(u)$ in $C(u)=\dfrac{\sum_i w_iP_iN_{i,p}}{D(u)}$. <br> <i>Hint:</i> ask what a non-positive $w_i$ does to $D(u)$ and whether the convex-combination argument still holds. :: A: a weight $w_i\le0$ can drive $D(u)$ to zero (a pole / division by zero) <br> the weights are no longer all positive, so the convex-hull guarantee fails and the curve can escape the hull <br> $\boxed{\text{keep all }w_i>0}$.
- Q: Pitfall: in the OCCT `Geom_BSplineCurve` constructor, how is the knot vector supplied? :: A: as distinct knots plus a parallel array of integer multiplicities, not the expanded vector <br> e.g. $\{0,0,0,\tfrac12,1,1,1\}$ is knots $\{0,\tfrac12,1\}$ with mults $\{3,1,3\}$ <br> $\boxed{\text{knots}+\text{mults arrays}}$.
- Q: When-to-use: which OCCT class fits/approximates a B-spline through ordered scan points, and which parameterization for noisy data? :: A: `GeomAPI_PointsToBSpline` fits/approximates through ordered points <br> choose centripetal parameterization for irregular/noisy spacing <br> $\boxed{\texttt{GeomAPI\_PointsToBSpline},\ \text{centripetal}}$.
