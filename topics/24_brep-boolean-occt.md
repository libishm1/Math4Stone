# B-Rep construction, booleans, sewing, and validation (OpenCASCADE)

Prereqs: [[nurbs-surfaces-fitting]], [[transforms-se3]] => this => [[cadfit-program-synthesis]], [[ifc-bim]]

## Why it matters for the stone scan->CAD->IFC / CADFit workflow

Fitting a NURBS patch to a stone scan gives you geometry, but a geometry is not a solid. IFC wants a closed, watertight, orientable solid (`IfcAdvancedBrep` or a faceted fallback), and to get there you must glue the fitted patches into one consistent topological skin and then prove it is valid before you serialize. That gluing step is sewing, the proving step is validity analysis, and the building blocks are faces, prisms, and booleans. When you cut a saw kerf or subtract a drill hole from the block model, that is a boolean Cut on B-Reps. When CADFit synthesizes a CAD program and you execute it, every extrude is a `MakePrism` and every union/subtract is a `BRepAlgoAPI` call, so the same kernel that exports IFC also runs the synthesized program. If you skip the validity gate, the bug does not surface in OpenCASCADE; it surfaces three steps later as a corrupt IFC file or a self-intersecting solid that the fabrication software silently rejects.

## Intuition

A B-Rep (boundary representation) describes a solid by its skin, not its interior. You never store "the inside of the block." You store the surfaces that bound it, the curves where surfaces meet, the points where curves meet, and a small web of "who touches whom." Think of a paper model of a polyhedron: the flat faces are surfaces, the fold creases are edges, the corners are vertices, and the tape holding adjacent panels together is the shared topology. A loose pile of cut-out panels is geometry; once you tape every seam so no gap remains, it becomes a closed shell, and a closed shell with a consistent inside-vs-outside is a solid.

Two ideas make this work in practice. First, every join is fuzzy: two panels that should share an edge never share it exactly in floating point, so the kernel sews them if their boundaries fall within a tolerance band. Second, orientation is load-bearing: each face has an outward side, and the kernel must agree on which side is "out" across the whole skin, or booleans cut the wrong region away.

## The math (first principles)

### The topology-geometry split

A B-Rep is two coupled layers. The **topology** is a graph of cells with adjacency. The **geometry** attaches a real surface/curve/point to each cell. OpenCASCADE names the topological cells, from high to low dimension:

$$\texttt{Compound} \supset \texttt{CompSolid} \supset \texttt{Solid} \supset \texttt{Shell} \supset \texttt{Face} \supset \texttt{Wire} \supset \texttt{Edge} \supset \texttt{Vertex}.$$

The geometry carriers are: a `Face` carries a `Geom_Surface` (often the NURBS patch you fitted), an `Edge` carries a 3D `Geom_Curve` and one or two 2D `p-curves` (the curve in each adjacent face's $(u,v)$ parameter space), and a `Vertex` carries a `gp_Pnt`. A `Wire` is an ordered loop of edges bounding a face. A `Shell` is a set of faces joined along edges. A `Solid` is the region bounded by one or more shells (an outer shell plus optional inner shells for cavities).

The split matters because validity is checked on **both** layers. The graph must be well formed (every edge bounded by vertices, every face bounded by closed wires) and the geometry must be consistent with it (the 3D curve, the two p-curves, and the surfaces must agree to within tolerance).

### Tolerances

Every topological entity carries a tolerance. A vertex tolerance $t_V$ is the radius of a sphere inside which the vertex "is." An edge tolerance $t_E$ is the radius of a tube around its curve. A face tolerance $t_F$ is the half-thickness of a slab around its surface. The kernel requires nesting:

$$t_V \ge t_E \ge t_F,$$

so a vertex's sphere always contains the edge tube which always contains the face slab. Default kernel confusion is `Precision::Confusion()` $= 10^{-7}$ in model units. Convention: OCCT is dimensionless and **right-handed**; you pick the length unit (the stone pipeline uses meters at site scale, millimeters at shop scale). A face's outward normal is defined by the surface's $(u,v)$ parameterization right-hand rule, $\hat n \propto S_u \times S_v$, **combined with** the face's stored `Orientation` flag (`FORWARD` or `REVERSED`), which can flip it. The effective outward normal is

$$\hat n_{\text{out}} = \sigma \,\frac{S_u \times S_v}{\lVert S_u \times S_v \rVert}, \qquad \sigma = \begin{cases} +1 & \text{FORWARD} \\ -1 & \text{REVERSED.}\end{cases}$$

Getting $\sigma$ wrong is the classic "boolean cut removed the wrong half" bug, and it traces straight back to the normal-orientation rule from [[transforms-se3]].

### Making a face and a prism

`BRepBuilderAPI_MakeFace` builds a `Face` from a surface plus a UV box, or from a surface plus a bounding wire, or from a planar wire alone. Given the fitted surface $S(u,v)$ on $[u_0,u_1]\times[v_0,v_1]$ and a fit tolerance $\varepsilon$, the face is the trimmed patch with that tolerance assigned to its edges.

`BRepPrimAPI_MakePrism` sweeps a shape along a vector $\vec d \in \mathbb{R}^3$, raising dimension by one: a vertex sweeps to an edge, an edge to a face, a wire to a shell, and a **face to a solid**. The swept solid from face $F$ along $\vec d$ is

$$\text{Prism}(F, \vec d) = \{\, p + s\,\vec d : p \in F,\; s \in [0,1] \,\}.$$

The result is a closed solid with the original face on top, a translated copy on the bottom, and side faces generated by sweeping the boundary wire. This is how a fitted top surface becomes a slab of thickness $\lVert \vec d \rVert$.

### Sewing: gluing faces into a shell

After fitting several patches you have a `Compound` of disconnected faces. Their shared boundaries do not coincide exactly. `BRepBuilderAPI_Sewing` with sewing tolerance $\tau$ merges two free edges into one shared edge when the maximum gap between them is within $\tau$:

$$\max_{t}\, \min_{t'}\, \lVert C_1(t) - C_2(t') \rVert \;\le\; \tau \quad\Longrightarrow\quad \text{merge } C_1, C_2 \text{ into one edge.}$$

Each edge of the assembled shell is then classified by how many faces share it. Let $n_e$ be that count:

- $n_e = 1$: **free edge** (boundary or a gap the tolerance missed),
- $n_e = 2$: **contiguous edge** (a properly sewn manifold seam),
- $n_e \ge 3$: **multiple / non-manifold edge**.

A watertight 2-manifold solid has **zero free edges** and every edge contiguous. The sewing class exposes `NbFreeEdges()`, `NbContigousEdges()` (note the kernel's spelling), and `NbMultipleEdges()` so you can audit this directly. Choosing $\tau$ is a trade: too small and real seams stay open (free edges remain); too large and you collapse distinct features or weld geometry that should not touch. A robust default ties $\tau$ to the data scale, e.g. $\tau \approx 10^{-3}\,L$ for characteristic size $L$, never below `Precision::Confusion()`.

### Boolean operations

Booleans are set operations on the point sets occupied by two solids $A$ and $B$. The B-Rep boolean computes the boundary of the result by intersecting all faces of $A$ with all faces of $B$, splitting both along the intersection curves, then keeping the pieces that belong to the target region:

$$\partial(A \cup B), \qquad \partial(A \setminus B), \qquad \partial(A \cap B).$$

In code these are `BRepAlgoAPI_Fuse` ($A \cup B$, union), `BRepAlgoAPI_Cut` ($A \setminus B$, difference, order matters), and `BRepAlgoAPI_Common` ($A \cap B$, intersection). `BRepAlgoAPI_Section` returns only the intersection **curves** $\partial A \cap \partial B$. Which faces survive is decided by classifying each split face's center point against the other solid (inside / outside / on), using the **outward normals**, which is exactly why orientation $\sigma$ must be consistent before you start.

Booleans are **not commutative for Cut** ($A \setminus B \ne B \setminus A$) and are sensitive to near-tangency. When two faces nearly touch within machine precision, the intersection is ill-posed and the boolean fails or produces sliver faces. The fix is a **fuzzy** boolean: set an extra tolerance $f$ via `SetFuzzyValue(f)` so the algorithm treats entities within $f$ as coincident. The default fuzzy value is $0$ (no fuzzing); raising it to roughly the gap size $f \sim$ (max gap from import/translation) makes imperfect, imported, or scan-derived geometry survive the operation.

### Validity: Euler-Poincaré and the analyzer

A necessary topological condition for a valid closed orientable B-Rep is the **Euler-Poincaré formula**:

$$V - E + F = 2(S - G),$$

where $V, E, F$ are vertex/edge/face counts, $S$ is the number of shells (the solid itself counts, so $S \ge 1$), and $G$ is the genus (number of through-handles). For a single simply-connected solid shell ($S=1$, $G=0$) this is the familiar $V - E + F = 2$. A cube: $V=8$, $E=12$, $F=6$, and $8 - 12 + 6 = 2$. Checked.

Crucial caveat: this is a **one-way test**. If $V - E + F \ne 2(S - G)$, the B-Rep is definitely broken. If it equals, the solid **may still be invalid** (the formula counts cells but does not see self-intersections, wrong orientations, or geometry-topology mismatch). So Euler-Poincaré is a cheap necessary filter, not a proof.

The real gate is `BRepCheck_Analyzer`, which checks both layers: closed wires, valid p-curves, `SameParameter` flags, self-intersecting wires, face-on-surface consistency, and orientation. `IsValid()` returns true only if no defect is found on the shape or any subshape. When it returns false, `Result(subshape).Status()` yields a list of `BRepCheck_Status` codes (e.g. `BRepCheck_InvalidCurveOnSurface`, `BRepCheck_SelfIntersectingWire`) telling you exactly what is wrong, which you then repair with the `ShapeFix` package.

### Controlled tessellation

For preview, fallback export, or any consumer that wants triangles, `BRepMesh_IncrementalMesh` tessellates the B-Rep. It is governed by two deflection bounds. The **linear deflection** $\delta_{\text{lin}}$ caps the chordal distance between the true curve/surface and its facet:

$$\max_{\text{facet}}\; \operatorname{dist}\big(S(u,v),\ \text{triangle plane}\big) \;\le\; \delta_{\text{lin}}.$$

The **angular deflection** $\delta_{\text{ang}}$ (radians) caps the angle between consecutive segments of a discretized edge, so high-curvature regions get more segments:

$$\angle(\,\vec t_i,\ \vec t_{i+1}\,) \;\le\; \delta_{\text{ang}}.$$

Smaller deflections mean finer, heavier meshes. If the `isRelative` flag is true, $\delta_{\text{lin}}$ is multiplied by each edge's size, giving scale-adaptive density; if false it is an absolute model-unit distance. The OCCT default angular deflection is $0.5$ rad; there is no default linear deflection (you must supply it). The rule: tessellate after validity passes, never before, because meshing a broken B-Rep hides the defect behind plausible-looking triangles.

## Worked example

Build a unit cube as a prism, sew it, check Euler-Poincaré by hand, then cut a corner.

Start with the bottom square face $F$ on the plane $z=0$ with corners $(0,0,0),(1,0,0),(1,1,0),(0,1,0)$. Sweep it along $\vec d = (0,0,1)$. The prism is the unit cube $[0,1]^3$.

Count its topology:

- Vertices: the 8 corners, so $V = 8$.
- Edges: 12 cube edges, $E = 12$.
- Faces: 6 cube faces, $F = 6$.
- Shells: one outer shell, $S = 1$. Genus: no handles, $G = 0$.

Euler-Poincaré: $V - E + F = 8 - 12 + 6 = 2$, and $2(S - G) = 2(1 - 0) = 2$. Equal. The necessary condition holds.

Now sewing. Each of the 12 edges is shared by exactly 2 faces, so after sewing `NbFreeEdges() = 0` and `NbContigousEdges() = 12`. Zero free edges confirms watertightness.

Now cut a smaller cube $B = [0.5,1.5]^3$ from $A = [0,1]^3$ with `BRepAlgoAPI_Cut(A, B)`. The overlap region is $[0.5,1]^3$, a corner of $A$. The result $A \setminus B$ is the cube with that corner notched out. Recount:

The notch removes the corner vertex $(1,1,1)$ and introduces new vertices where $B$'s faces slice $A$. The three planes $x=0.5$, $y=0.5$, $z=0.5$ each cut into $A$. Working it through, the L-shaped (in 3D, a cube-minus-corner) result has $V = 10$, $E = 15$, $F = 7$. Check: $10 - 15 + 7 = 2 = 2(1 - 0)$. Still genus-0, single shell, Euler holds. The arithmetic confirms the boolean preserved a valid closed solid. If instead you had computed $8 - 12 + 6 = 2$ and the boolean had silently left a free edge, the analyzer (not Euler) would catch it.

A trap worth noting: if you had run `BRepAlgoAPI_Cut(B, A)` you would get $B \setminus A = [1,1.5]^3 \cup \ldots$, a different shape, because Cut is non-commutative. Order is part of the math.

## Code

Idiomatic C++ with OpenCASCADE. This is the full path: fit-surface-stand-in -> face -> prism -> sew -> validate -> boolean -> revalidate -> tessellate. Each call is annotated with the math symbol it implements.

```cpp
// brep_pipeline.cpp
// Build (Linux, OCCT 7.x):
//   g++ -std=c++17 brep_pipeline.cpp -o brep_pipeline \
//       -I/usr/include/opencascade \
//       -lTKernel -lTKMath -lTKG3d -lTKGeomBase -lTKBRep \
//       -lTKTopAlgo -lTKPrim -lTKBO -lTKShHealing -lTKMesh
#include <gp_Pnt.hxx>
#include <gp_Vec.hxx>
#include <BRepBuilderAPI_MakePolygon.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepPrimAPI_MakePrism.hxx>
#include <BRepPrimAPI_MakeBox.hxx>
#include <BRepBuilderAPI_Sewing.hxx>
#include <BRepAlgoAPI_Cut.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <BRepMesh_IncrementalMesh.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>
#include <TopAbs.hxx>
#include <iostream>

// Count V, E, F and report the Euler-Poincare value V - E + F.
static void eulerReport(const TopoDS_Shape& s, const char* tag) {
    auto countType = [&](TopAbs_ShapeEnum t) {
        int n = 0;
        // TopExp_Explorer walks unique subshapes of a given type.
        for (TopExp_Explorer e(s, t); e.More(); e.Next()) ++n;
        return n;
    };
    int V = countType(TopAbs_VERTEX);   // note: counts oriented uses; see pitfalls
    int E = countType(TopAbs_EDGE);
    int F = countType(TopAbs_FACE);
    std::cout << tag << ":  V=" << V << " E=" << E << " F=" << F
              << "  V-E+F=" << (V - E + F) << "\n";
}

int main() {
    const double tol = 1.0e-7;   // ~ Precision::Confusion(), model units

    // --- 1. Face: a unit square wire on z=0, closed into a planar face. -------
    // Stands in for BRepBuilderAPI_MakeFace(fittedSurface, tol) in the real
    // pipeline, where the surface is the NURBS patch from prereq fitting.
    BRepBuilderAPI_MakePolygon poly(
        gp_Pnt(0,0,0), gp_Pnt(1,0,0), gp_Pnt(1,1,0), gp_Pnt(0,1,0), true /*close*/);
    TopoDS_Wire   wire = poly.Wire();
    TopoDS_Face   face = BRepBuilderAPI_MakeFace(wire, true /*planar*/);

    // --- 2. Prism: sweep the face along d = (0,0,1) -> a unit-cube solid. -----
    // Prism(F, d) = { p + s*d : p in F, s in [0,1] }.
    TopoDS_Shape A = BRepPrimAPI_MakePrism(face, gp_Vec(0,0,1));
    eulerReport(A, "cube (prism)");                // expect V-E+F = 2

    // --- 3. Sew: glue any contiguous faces; audit free vs contiguous edges. ---
    // tolerance tau, then four option flags: sew, analyze-degenerate,
    // cut-free-edges, non-manifold.
    BRepBuilderAPI_Sewing sewer(1.0e-6, true, true, true, false);
    sewer.Add(A);
    sewer.Perform();
    std::cout << "free=" << sewer.NbFreeEdges()
              << " contiguous=" << sewer.NbContigousEdges()
              << " multiple=" << sewer.NbMultipleEdges() << "\n"; // free=0 wanted
    TopoDS_Shape Asewn = sewer.SewedShape();

    // --- 4. Validate BEFORE any boolean or export. ---------------------------
    if (!BRepCheck_Analyzer(Asewn).IsValid()) {
        std::cerr << "cube invalid after sewing\n"; return 1;
    }

    // --- 5. Boolean Cut: A \ B, removing the [0.5,1.5]^3 corner. --------------
    // Cut is non-commutative: Cut(A,B) != Cut(B,A).
    TopoDS_Shape B   = BRepPrimAPI_MakeBox(gp_Pnt(0.5,0.5,0.5), 1.0, 1.0, 1.0).Shape();
    BRepAlgoAPI_Cut cut(Asewn, B);
    // cut.SetFuzzyValue(1.0e-5);   // uncomment for imperfect/imported geometry
    cut.Build();
    if (!cut.IsDone()) { std::cerr << "boolean failed\n"; return 1; }
    TopoDS_Shape result = cut.Shape();
    eulerReport(result, "cube minus corner");      // expect V-E+F = 2 again

    // --- 6. Re-validate the boolean output. ----------------------------------
    if (!BRepCheck_Analyzer(result).IsValid()) {
        std::cerr << "result invalid after Cut\n"; return 1;
    }

    // --- 7. Tessellate LAST, only after validity passes. ---------------------
    // linDefl = 0.01 absolute, isRelative=false, angDefl=0.5 rad, parallel=true.
    BRepMesh_IncrementalMesh mesher(result, 0.01, false, 0.5, true);
    mesher.Perform();
    std::cout << "tessellation done: " << (mesher.IsDone() ? "ok" : "FAILED") << "\n";

    return 0;
}
```

Symbol map: `MakeFace` is the trimmed patch with tolerance $\varepsilon$; `MakePrism(face, gp_Vec)` is $\text{Prism}(F,\vec d)$; the `BRepBuilderAPI_Sewing` constructor args are $(\tau,\,\text{sew},\,\text{analyze},\,\text{cut-free},\,\text{non-manifold})$; `NbFreeEdges` is the count of $n_e = 1$ edges; `BRepAlgoAPI_Cut` is $A \setminus B$; `SetFuzzyValue` is $f$; `IsValid` is the analyzer gate; `BRepMesh_IncrementalMesh(s, linDefl, isRelative, angDefl, parallel)` are $(\delta_{\text{lin}},\,\text{rel},\,\delta_{\text{ang}})$.

A tiny Python (pythonOCC) version of the validity gate, useful for quick scripting against the same kernel:

```python
# pip install pythonocc-core  (conda-forge build recommended)
from OCC.Core.BRepPrimAPI import BRepPrimAPI_MakeBox
from OCC.Core.BRepAlgoAPI import BRepAlgoAPI_Cut
from OCC.Core.BRepCheck import BRepCheck_Analyzer

A = BRepPrimAPI_MakeBox(1, 1, 1).Shape()          # unit cube at origin
B = BRepPrimAPI_MakeBox(gp_Pnt(0.5, 0.5, 0.5), 1, 1, 1).Shape()  # offset cube
cut = BRepAlgoAPI_Cut(A, B)                        # A \ B
cut.SetFuzzyValue(1e-5)                            # f, for robustness
cut.Build()
assert cut.IsDone()
result = cut.Shape()
assert BRepCheck_Analyzer(result).IsValid()        # the gate, same as C++
```

## Connections

- [[nurbs-surfaces-fitting]]: the `Geom_Surface` carried by each `Face` is the fitted NURBS patch; `MakeFace(surface, tol)` is the bridge from a fitted patch to a topological face, and the fit tolerance seeds the face tolerance $t_F$.
- [[transforms-se3]]: every face is placed by a `gp_Trsf` ($SE(3)$), and the outward-normal/orientation flag $\sigma$ is the inverse-transpose normal rule applied at the topology level; consistent normals are what make booleans cut the right region.
- [[cadfit-program-synthesis]]: a synthesized CAD program is a sequence of `MakePrism` extrudes and `BRepAlgoAPI` booleans on this same kernel, so executing the program and validating it reuses everything here.
- [[ifc-bim]]: a validated, watertight `Solid` becomes an `IfcAdvancedBrep` (exact) or, if it fails the gate, a tessellated `IfcTriangulatedFaceSet` fallback produced by `BRepMesh_IncrementalMesh`.
- [[robust-ransac]] / segmentation: when one patch cannot span an overhang, you segment into multiple patches and sew them; the free-edge audit tells you whether the segmentation seam actually closed.

## Pitfalls and failure modes

- **Orientation inconsistency.** If a face's `Orientation` flag $\sigma$ disagrees with its neighbors, the boolean classifies the wrong side as inside and removes the wrong region, or produces an inside-out solid. Symptom: Cut "deletes the whole part" or leaves a shell with no volume. Fix orientations with `ShapeFix` or rebuild the shell with consistent face normals.
- **Sewing tolerance too small or too large.** Too small leaves real seams as free edges (not watertight, IFC export rejects it). Too large welds distinct features or collapses thin geometry. Always check `NbFreeEdges() == 0` and tie $\tau$ to data scale, never below `Precision::Confusion()`.
- **Tolerance non-nesting.** The kernel requires $t_V \ge t_E \ge t_F$. After heavy editing this can break and the analyzer reports `BRepCheck_InvalidTolerance...`. Running `ShapeFix_Shape` re-imposes the hierarchy.
- **Boolean near-tangency / coincident faces.** Two faces touching within machine precision make the intersection ill-conditioned; the boolean fails or emits sliver faces. Set a `SetFuzzyValue(f)` near the real gap size. This is the conditioning issue from least-squares carried into geometry.
- **Trusting Euler-Poincaré as proof.** $V - E + F = 2(S - G)$ is necessary, not sufficient. A self-intersecting or wrongly oriented solid can still satisfy it. Always run `BRepCheck_Analyzer` as the actual gate.
- **Counting subshapes naively.** `TopExp_Explorer` can return oriented uses, so a raw walk may double-count shared vertices/edges if you do not deduplicate (use a `TopTools_IndexedMapOfShape`). A wrong $V$ makes a valid solid look like it fails Euler.
- **Cut non-commutativity.** $A \setminus B \ne B \setminus A$. Always pass (tool-to-keep, tool-to-remove) in that order; swapping silently produces the complementary chunk.
- **Tessellating before validating.** Meshing a broken B-Rep hides the defect behind plausible triangles, and the corrupt solid reaches IFC anyway. Validity gate first, mesh last.
- **Degenerate faces from a bad NURBS fit.** A rippled or self-overlapping fitted patch produces a face whose p-curves do not match the 3D curve; the analyzer flags `BRepCheck_InvalidCurveOnSurface`. Re-fit with better parameterization upstream rather than patching the topology.

## Practice

1. By hand, verify Euler-Poincaré for a triangular prism (5 faces, 9 edges, 6 vertices) and for a torus-like solid with one through-hole ($G=1$). Confirm $V - E + F$ equals $2$ and $0$ respectively. (Inspectable: two integers.)
2. In OCCT, build a unit cube two ways: `BRepPrimAPI_MakeBox` and a `MakePrism` of a square face. Run `BRepCheck_Analyzer` on both, and print `NbFreeEdges()` after sewing each. Confirm both give zero free edges. (Inspectable: validity bool + edge counts.)
3. Deliberately break a shell: take the cube, remove one face (rebuild the shell without it), sew, and report `NbFreeEdges()`. Confirm it jumps to 4 and that `IsValid()` is false. Then add the face back and confirm recovery.
4. Cut a cylindrical hole through a slab with `BRepAlgoAPI_Cut`, then recompute $V - E + F$ and confirm $G$ increased to $1$ (so $V - E + F = 0$). Validate, then tessellate at $\delta_{\text{lin}} = 0.1$ and again at $0.01$ and compare triangle counts. (Inspectable: Euler value + two mesh sizes.)
5. Construct two boxes that share a face exactly, Fuse them, and observe whether you get one solid or a shell with a free internal edge. Then offset one box by $10^{-9}$ and repeat with and without `SetFuzzyValue(1e-7)`. Report which combinations produce a valid single solid. (Inspectable: validity table.)
6. Take a fitted NURBS patch (from the prereq topic), `MakeFace` it, `MakePrism` to a slab, sew, validate, and export both an exact and a tessellated fallback. Diff the two: measure the max chordal deviation and confirm it is below your chosen $\delta_{\text{lin}}$.

## References

- Open CASCADE Technology, Modeling Algorithms user guide (faces, prisms, Boolean operations, sewing, validity): https://dev.opencascade.org/doc/overview/html/occt_user_guides__modeling_algos.html
- Open CASCADE Technology, Boolean Operations specification (fuzzy operations, General Fuse, `SetFuzzyValue`): https://dev.opencascade.org/doc/overview/html/specification__boolean_operations.html
- OCCT reference manual, `BRepBuilderAPI_Sewing` (constructor `tolerance=1e-6` plus four option flags, `NbFreeEdges`, `NbContigousEdges`, `SewedShape`): https://dev.opencascade.org/doc/refman/html/class_b_rep_builder_a_p_i___sewing.html
- OCCT reference manual, `BRepMesh_IncrementalMesh` (constructor `(shape, linDefl, isRelative=false, angDefl=0.5, isInParallel=false)`): https://dev.opencascade.org/doc/refman/html/class_b_rep_mesh___incremental_mesh.html
- OCCT reference manual, `BRepCheck_Analyzer` (`IsValid`, `Result`, `BRepCheck_Status`): https://dev.opencascade.org/doc/refman/html/class_b_rep_check___analyzer.html
- OCCT reference manual, `BRepAlgoAPI_BooleanOperation` / `BRepAlgoAPI_Cut` / `_Fuse` / `_Common` / `_Section`: https://dev.opencascade.org/doc/refman/html/class_b_rep_algo_a_p_i___boolean_operation.html
- Open CASCADE Technology, Mesh user guide (linear vs angular deflection meaning): https://dev.opencascade.org/doc/overview/html/occt_user_guides__mesh.html
- Ching-Kuang Shene, CS3621 notes, "The Euler-Poincaré Formula" ($V - E + F = 2(S - G)$, one-way test): https://pages.mtu.edu/~shene/COURSES/cs3621/NOTES/model/euler.html
- pythonOCC documentation (Python bindings to the same kernel): https://github.com/tpaviot/pythonocc-core
- IfcOpenShell documentation (geometry processing, `IfcAdvancedBrep` and triangulated fallback paths): https://docs.ifcopenshell.org/
- G. Nehme, E. Whalen, F. Ahmed, "CADFit: Precise Mesh-to-CAD Program Generation with Hybrid Optimization," arXiv:2605.01171, and code https://github.com/ghadinehme/CADFit (the synthesized programs execute as prism/boolean sequences on this kernel).

## Flashcard seeds

- Q: Derive the un-flipped unit normal $\hat n_{uv}$ of a face from its surface parameterization $S(u,v)$ (step 1 of 2). <br> <b>Formula:</b> $S_u=\partial S/\partial u,\ S_v=\partial S/\partial v,\ \hat n_{uv}=\dfrac{S_u\times S_v}{\lVert S_u\times S_v\rVert}$ <br> <i>Hint:</i> the cross product of the two tangents is normal to the patch; divide by its length to make it a unit vector. :: A: The two tangent directions are $S_u=\partial S/\partial u$ and $S_v=\partial S/\partial v$ <br> their cross product is normal to the patch, with magnitude $\lVert S_u\times S_v\rVert$ (the area scale) <br> normalize to get the parameterization normal <br> $\boxed{\hat n_{uv}=\dfrac{S_u\times S_v}{\lVert S_u\times S_v\rVert}}$
- Q: Fold the stored `Orientation` flag into $\hat n_{uv}$ to get the effective outward normal $\hat n_{\text{out}}$ (step 2 of 2). <br> <b>Formula:</b> $\hat n_{\text{out}}=\sigma\,\hat n_{uv}$, with $\sigma=+1$ (FORWARD), $\sigma=-1$ (REVERSED) <br> <i>Hint:</i> just multiply the step-1 unit normal by the sign $\sigma$; REVERSED flips it. :: A: The face stores a sign $\sigma=+1$ for FORWARD and $\sigma=-1$ for REVERSED, which may flip the parameterization normal <br> multiply through so the result points out of the solid regardless of how $(u,v)$ was laid out <br> $\boxed{\hat n_{\text{out}}=\sigma\,\dfrac{S_u\times S_v}{\lVert S_u\times S_v\rVert}}$
- Q: Derive the point set of the solid produced by sweeping a face $F$ along $\vec d$ (the `MakePrism` set). <br> <b>Formula:</b> $\text{Prism}(F,\vec d)=\{\,p+s\vec d:\;p\in F,\;s\in[0,1]\,\}$ <br> <i>Hint:</i> let each base point $p$ ride the direction $\vec d$ as the travel fraction $s$ runs from 0 to 1, then take the union. :: A: Every base point $p\in F$ travels along the fixed direction $\vec d$ <br> a fractional travel $s\in[0,1]$ places it at $p+s\vec d$ ($s=0$ is the base face, $s=1$ the translated copy) <br> the solid is the union over all base points and all travel fractions <br> $\boxed{\text{Prism}(F,\vec d)=\{\,p+s\vec d:\;p\in F,\;s\in[0,1]\,\}}$
- Q: From the prism set, derive the slab thickness $h$ produced by sweeping a planar face along $\vec d$. <br> <b>Formula:</b> $\text{Prism}(F,\vec d)=\{p+s\vec d:p\in F,s\in[0,1]\}\Rightarrow h=\lVert\vec d\rVert$ <br> <i>Hint:</i> the base is at $s=0$ and its copy at $s=1$; the separation is the full length of $\vec d$. :: A: The base face sits at $s=0$, its copy at $s=1$, displaced by the full vector $\vec d$ <br> for a planar face the thickness is the distance between those parallel copies, i.e. the length of $\vec d$ <br> $\boxed{h=\lVert\vec d\rVert}$
- Q: Define the worst-case gap $g$ between two free edges carrying curves $C_1(t),C_2(t')$ (sewing, step 1 of 2). <br> <b>Formula:</b> $g=\max_{t}\,\min_{t'}\,\lVert C_1(t)-C_2(t')\rVert$ <br> <i>Hint:</i> for each point on $C_1$ first take the closest approach on $C_2$ (the inner $\min$), then the worst of those over $C_1$ (the outer $\max$). :: A: For each point $C_1(t)$ the closest approach on the other curve is $\min_{t'}\lVert C_1(t)-C_2(t')\rVert$ <br> the curves are at most as far apart as the worst such closest approach <br> $\boxed{g=\max_{t}\,\min_{t'}\,\lVert C_1(t)-C_2(t')\rVert}$ (a one-sided Hausdorff gap)
- Q: State the sewing decision rule that merges two free edges, given the gap $g$ and tolerance $\tau$ (step 2 of 2). <br> <b>Formula:</b> $g\le\tau\;\Longrightarrow\;\text{merge }C_1,C_2\text{ into one edge}$ <br> <i>Hint:</i> compare the worst gap $g$ to the tolerance band $\tau$; merge only when it fits inside. :: A: Two free edges become one shared edge only if their worst gap fits inside the tolerance band <br> compare the gap to $\tau$ <br> $\boxed{g\le\tau\;\Longrightarrow\;\text{merge }C_1,C_2\text{ into one edge}}$
- Q: Derive a scale-robust default sewing tolerance $\tau$ for data of characteristic size $L$. <br> <b>Formula:</b> size-relative term $10^{-3}L$, confusion floor $10^{-7}$, so $\tau=\max(10^{-3}L,\,10^{-7})$ <br> <i>Hint:</i> take the larger of the size-relative term and the kernel floor so $\tau$ is never below `Precision::Confusion()`. :: A: A fixed absolute $\tau$ is wrong across scales (too big in mm, too small in m), so tie it to size: $\tau\approx 10^{-3}L$ <br> but it must never drop below the kernel confusion floor $10^{-7}$ <br> $\boxed{\tau=\max\!\left(10^{-3}L,\;10^{-7}\right)}$
- Q: From the edge-share count $n_e$ (faces sharing an edge), derive the watertight 2-manifold condition over a shell. <br> <b>Formula:</b> $n_e=1$ free, $n_e=2$ contiguous, $n_e\ge3$ non-manifold; watertight 2-manifold $\iff$ every edge has $n_e=2$ <br> <i>Hint:</i> watertight bans $n_e=1$ (no boundary), 2-manifold bans $n_e\ge3$ (no triple edge); only $n_e=2$ survives. :: A: $n_e=1$ is a free (boundary) edge, $n_e=2$ a contiguous seam, $n_e\ge3$ non-manifold <br> watertight means no boundary survives, so no edge has $n_e=1$ <br> 2-manifold means no edge is shared by three or more faces, so no $n_e\ge3$ <br> $\boxed{\text{every edge has }n_e=2\iff\text{NbFreeEdges}=0\text{ and no multiple edges}}$
- Q: Specialize the Euler-Poincaré formula to a simply-connected single-shell solid (step 1 of 2). <br> <b>Formula:</b> $V-E+F=2(S-G)$, with $S=1,\ G=0$ <br> <i>Hint:</i> substitute one shell and zero handles into the right-hand side and evaluate $2(S-G)$. :: A: A single closed shell gives $S=1$; no through-handles gives genus $G=0$ <br> substitute into the right-hand side: $2(S-G)=2(1-0)$ <br> $\boxed{V-E+F=2}$
- Q: Specialize the Euler-Poincaré formula to a single-shell solid with one through-hole (step 2 of 2, torus regime). <br> <b>Formula:</b> $V-E+F=2(S-G)$, with $S=1,\ G=1$ <br> <i>Hint:</i> one through-handle means $G=1$; plug into $2(S-G)$. :: A: One through-handle means genus $G=1$; still a single shell $S=1$ <br> right-hand side: $2(S-G)=2(1-1)=0$ <br> $\boxed{V-E+F=0}$
- Q: Count $V,E,F$ for the unit cube prism and verify Euler-Poincaré. <br> <b>Formula:</b> $V-E+F=2(S-G)$; cube has $V=8,E=12,F=6,S=1,G=0$ <br> <i>Hint:</i> count corners, edges, faces, then check the left side equals $2(1-0)$. :: A: 8 corners so $V=8$, 12 edges so $E=12$, 6 faces so $F=6$ <br> $V-E+F=8-12+6=2$ <br> $S=1,G=0\Rightarrow 2(S-G)=2$ <br> $\boxed{8-12+6=2=2(1-0)\ \checkmark}$
- Q: Worked example (m site scale): set the sewing tolerance for a slab whose parts span $L=2$ m. <br> <b>Formula:</b> $\tau=\max(10^{-3}L,\,10^{-7})$ <br> <i>Hint:</i> compute the size-relative term $10^{-3}\cdot 2$ first, then compare it to the $10^{-7}$ floor. :: A: Use $\tau=\max(10^{-3}L,10^{-7})$ with $L=2$ <br> $10^{-3}\cdot 2=2\times10^{-3}$ m, well above the $10^{-7}$ floor <br> $\boxed{\tau=2\times10^{-3}\text{ m}=2\text{ mm}}$
- Q: Worked example (mm shop scale): a carved detail has characteristic size $L=5$ mm $=0.005$ m. Compute $\tau$. <br> <b>Formula:</b> $\tau=\max(10^{-3}L,\,10^{-7})$ <br> <i>Hint:</i> work in metres: evaluate $10^{-3}\cdot0.005$, then take the max with $10^{-7}$. :: A: $\tau=\max(10^{-3}\cdot0.005,\,10^{-7})=\max(5\times10^{-6},\,10^{-7})$ <br> $5\times10^{-6}>10^{-7}$, so the size-relative term wins <br> $\boxed{\tau=5\times10^{-6}\text{ m}=5\ \mu\text{m}}$
- Q: Worked example (site/large scale): a quarry block model has $L=40$ m. Compute $\tau$ and contrast with the mm-scale case. <br> <b>Formula:</b> $\tau=\max(10^{-3}L,\,10^{-7})$ <br> <i>Hint:</i> evaluate $10^{-3}\cdot40$; the floor $10^{-7}$ is irrelevant here. :: A: $\tau=\max(10^{-3}\cdot40,\,10^{-7})=\max(0.04,\,10^{-7})$ <br> $\boxed{\tau=0.04\text{ m}=40\text{ mm}}$ <br> 40 mm at site scale vs 5 µm at detail scale shows why $\tau$ must be size-relative, not a fixed constant
- Q: Worked example (boolean Cut): from $A=[0,1]^3$ cut $B=[0.5,1.5]^3$. Find the overlap region and verify Euler on the result. <br> <b>Formula:</b> overlap $=A\cap B$; check $V-E+F=2(S-G)$ with result counts $V=10,E=15,F=7$ <br> <i>Hint:</i> intersect the two boxes coordinate-by-coordinate first, then plug the post-cut counts into Euler. :: A: Overlap is the intersection $[0.5,1]^3$, a corner of $A$ <br> the notched solid has $V=10,E=15,F=7$ <br> $\boxed{10-15+7=2=2(1-0)\ \checkmark}$ still a valid genus-0 single shell
- Q: Worked example (prism thickness, m scale): sweep a fitted top face along $\vec d=(0,0,0.3)$ m. Give the slab thickness. <br> <b>Formula:</b> $h=\lVert\vec d\rVert$ <br> <i>Hint:</i> the vector is axis-aligned, so its length is just the magnitude of the nonzero component. :: A: Thickness is the sweep length $h=\lVert\vec d\rVert=\lVert(0,0,0.3)\rVert$ <br> $\boxed{h=0.3\text{ m}=300\text{ mm}}$
- Q: Worked example (oblique prism, mm scale): sweep a face along $\vec d=(3,4,12)$ mm. Compute the slab thickness $h$. <br> <b>Formula:</b> $h=\lVert\vec d\rVert=\sqrt{d_x^2+d_y^2+d_z^2}$ <br> <i>Hint:</i> square the three components, add, and take the square root (it is a Pythagorean triple). :: A: $h=\lVert\vec d\rVert=\sqrt{3^2+4^2+12^2}=\sqrt{9+16+144}=\sqrt{169}$ <br> $\boxed{h=13\text{ mm}}$
- Q: Closing test (close the loop): verify Euler-Poincaré for a triangular prism ($V=6,E=9,F=5$) as a single genus-0 shell. <br> <b>Formula:</b> $V-E+F=2(S-G)$, with $S=1,\ G=0$ <br> <i>Hint:</i> a pass means $V-E+F$ comes out to $2(1-0)=2$. :: A: $V-E+F=6-9+5=2$ <br> single shell $S=1$, no handle $G=0\Rightarrow2(S-G)=2$ <br> $\boxed{2=2\ \checkmark}$ the necessary condition holds (analyzer still needed for sufficiency)
- Q: Closing test (close the loop): cut a cylindrical through-hole in a slab so genus $G=1$. Predict $V-E+F$ and check the torus limit. <br> <b>Formula:</b> $V-E+F=2(S-G)$, with $S=1,\ G=1$ <br> <i>Hint:</i> a pass means the prediction equals the known one-handle (torus) Euler characteristic of $0$. :: A: From $V-E+F=2(S-G)$ with $S=1,G=1$: predicted value $=2(1-1)=0$ <br> verify against the known torus regime: a one-handle solid indeed gives Euler characteristic $0$ <br> $\boxed{V-E+F=0\ \checkmark}$ matches the genus-1 limiting case
- Q: Closing test (orientation sanity check): a face's $\sigma$ flag flips FORWARD$\to$REVERSED. Show the effect on $\hat n_{\text{out}}$ and why Cut then fails. <br> <b>Formula:</b> $\hat n_{\text{out}}=\sigma\,\hat n_{uv}$, with $\sigma:+1\to-1$ <br> <i>Hint:</i> a pass is realizing the normal negates, so the inside/outside classification inverts. :: A: $\hat n_{\text{out}}=\sigma\,\hat n_{uv}$, so $\sigma:+1\to-1$ sends $\hat n_{\text{out}}\to-\hat n_{\text{out}}$ <br> the boolean classifies face centers by the outward normal; a reversed normal labels "inside" as "outside" <br> $\boxed{\text{Cut removes the wrong region (inside-out solid)}}$ confirming $\sigma$ must be consistent before any boolean
- Q: When do you tessellate in the pipeline, and why not earlier? <br> <i>Hint:</i> think about which gate must pass first, and what a mesh can hide. :: A: Last, only after `BRepCheck_Analyzer.IsValid()` passes; meshing a broken B-Rep hides the defect behind plausible triangles that still reach IFC.
- Q: Is Euler-Poincaré ($V-E+F=2(S-G)$) sufficient to prove a B-Rep is valid? <br> <i>Hint:</i> recall it is a one-way test; ask what it cannot see. :: A: No. It is necessary only; a self-intersecting or misoriented solid can still satisfy it. The real gate is `BRepCheck_Analyzer`.
- Q: Name the OCCT topological hierarchy from highest to lowest cell. <br> <i>Hint:</i> eight names, starting at `Compound` and ending at `Vertex`. :: A: Compound, CompSolid, Solid, Shell, Face, Wire, Edge, Vertex.
- Q: What are the three set-operation booleans and their OCCT classes? <br> <i>Hint:</i> union / difference / intersection map to Fuse / Cut / Common. :: A: Union `BRepAlgoAPI_Fuse`, difference `BRepAlgoAPI_Cut`, intersection `BRepAlgoAPI_Common` (and `_Section` for intersection curves only).
- Q: Pitfall: why is `BRepAlgoAPI_Cut` order-sensitive? <br> <b>Formula:</b> $A\setminus B\ne B\setminus A$ <br> <i>Hint:</i> set difference is not commutative; identify which argument is kept. :: A: Difference is non-commutative, $A\setminus B\ne B\setminus A$; always pass (shape-to-keep, shape-to-remove).
- Q: What does `SetFuzzyValue(f)` do and what is its default? <br> <i>Hint:</i> it is an extra coincidence tolerance; recall the default is zero. :: A: Adds tolerance $f$ so near-coincident entities are treated as touching; default is $0$ (no fuzzing). Raise it near the real gap size for imported/scan geometry.
- CLOZE: An OCCT `Edge` carries a 3D curve plus one or two {{c1::p-curves}} (the curve in each adjacent face's (u,v) space).
- CLOZE: The kernel requires the tolerance nesting {{c1::$t_V \ge t_E \ge t_F$}} so the vertex sphere contains the edge tube contains the face slab.
- CLOZE: After sewing, a watertight 2-manifold solid has {{c1::zero}} free edges and every edge is contiguous (shared by two faces).
- CLOZE: OCCT's default angular deflection is {{c1::0.5}} rad, while the linear deflection has {{c2::no default}} (the caller must supply it). <br> <i>hint:</i> one is a fixed radian value, the other the caller must always pass.
- CLOZE: Default `Precision::Confusion()` in OCCT model units is {{c1::$10^{-7}$}}.
