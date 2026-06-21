# IFC serialization and BIM semantics (IfcOpenShell)

Placement: [[24_brep-boolean-occt]] => **IFC / BIM** => (handoff to BIM and fabrication)

## Why it matters for the stone scan->CAD->IFC / CADFit / robotic-fabrication workflow

IFC is the last stage of the pipeline: it packages the reconstructed geometry as a semantic building model that other tools can open, query, and fabricate from. A wall, a vault stone, or a slab becomes a named IFC element with geometry, a placement, materials, and relationships, not just a loose mesh. IfcOpenShell is the open-source bridge that reads and writes IFC and turns its geometric representations into triangles or OpenCASCADE B-Reps. Getting placement transforms, units, and the B-Rep-vs-tessellation choice right here is what makes the whole scan-to-BIM effort usable downstream.

## Intuition

A mesh says "here are some triangles." IFC says "here is a *wall*, made of *limestone*, located *here* relative to *this storey*, with *this* shape, *contained in* this building." It is a structured database of building objects (an entity graph) plus their geometry and relationships. Think of it as the difference between a photo of a car and the car's parts list with where each part bolts on.

## The math (first principles)

IFC is mostly data modeling, but two pieces are genuinely mathematical.

**Placement is composed rigid transforms.** Each element carries an `IfcLocalPlacement` that is relative to its parent's placement, forming a chain. The world transform of an element is the product of 4x4 rigid transforms up the chain (the same $SE(3)$ from [[05_transforms-se3]]):

$$ T_\text{world} = T_\text{site}\,T_\text{building}\,T_\text{storey}\,T_\text{element}. $$

An `IfcAxis2Placement3D` defines a local frame by an origin, a $Z$ axis (the "Axis"), and an $X$ axis (the "RefDirection"); $Y = Z \times X$ rebuilds an orthonormal rotation $R$, and with the origin $t$ you get $T = \begin{bmatrix} R & t \\ 0 & 1\end{bmatrix}$.

**Units scale everything.** IFC stores a length unit (often millimetres). Every coordinate must be multiplied by the model's length-unit scale to reach SI metres. A unit mismatch is the single most common scan-to-BIM bug.

**Geometry representation choice.** An element's shape can be an exact curved boundary (`IfcAdvancedBrep`, carrying the NURBS/B-Rep from [[24_brep-boolean-occt]]) or a faceted approximation (`IfcTriangulatedFaceSet` / tessellation). Exactness trades off against interoperability: advanced B-Reps preserve curvature but not every viewer supports them; tessellated geometry is universally readable but approximate. IfcOpenShell can emit either.

**Bulk processing.** For large models, IfcOpenShell exposes a geometry *iterator* that processes elements with multicore support, caching, and reuse of repeated geometry, returning per-element vertex/face arrays plus the 4x4 placement matrix.

## Worked example

A stone block sits at local origin $(0,0,0)$ in a storey whose placement translates it by $(10, 4, 3)$ m and the building adds a $90^\circ$ rotation about $Z$. The block's world position is $T_\text{building}\,T_\text{storey}\,p$. Take $p = (1,0,0)$ m in block coordinates: the storey moves it to $(11, 4, 3)$, then the building's $90^\circ$ $Z$-rotation about the origin sends $(11,4,3) \to (-4, 11, 3)$. If the model's length unit were millimetres and you forgot to scale, the same point would land at $(-4000, 11000, 3000)$, off by a factor of 1000. This is the units pitfall in one line.

## Code

IfcOpenShell C++ geometry iterator: pull triangulated geometry and the 4x4 placement for every element. (Python is common for authoring; the iterator pattern is the same.)

```cpp
#include <ifcgeom/IfcGeomIterator.h>
#include <ifcparse/IfcFile.h>

void exportTriangles(const std::string& ifcPath) {
    IfcParse::IfcFile file(ifcPath);

    ifcopenshell::geometry::Settings settings;          // request world coords, welding, etc.
    settings.set("use-world-coords", true);             // bake placements into vertices

    IfcGeom::Iterator it(settings, &file);              // multicore, cached iterator
    if (!it.initialize()) return;

    do {
        const IfcGeom::TriangulationElement* e =
            static_cast<const IfcGeom::TriangulationElement*>(it.get());

        const auto& mesh   = e->geometry();
        const auto& verts  = mesh.verts();              // flat [x,y,z, x,y,z, ...] in metres
        const auto& faces  = mesh.faces();              // index triples
        const auto& M      = e->transformation().matrix(); // 4x4 placement (if not world-coords)

        // ... write verts/faces to your target format, applying M if needed ...
    } while (it.next());
}
// Alternatively, create_shape() can yield an OpenCASCADE TopoDS_Shape (advanced B-Rep path).
```

## Connections

- [[05_transforms-se3]] is the placement algebra: IFC placements are chained 4x4 rigid transforms.
- [[24_brep-boolean-occt]] produces the B-Rep that becomes an `IfcAdvancedBrep`; the tessellated path uses OCCT meshing.
- [[25_cadfit-program-synthesis]] recovers editable operations whose result is exported here as semantic elements.
- [[17_mesh-topology]] is the triangulated representation IfcOpenShell returns for viewers.

## Pitfalls and failure modes

- **Units.** Forgetting the model length-unit scale (often mm) corrupts every coordinate by a constant factor. Read the unit assignment and scale to metres.
- **Placement chains.** Each placement is relative to its parent; you must compose the whole chain (or ask IfcOpenShell for world coordinates).
- **Schema version.** IFC2x3 vs IFC4 differ in entities (e.g. tessellation and advanced B-Reps are IFC4). Target the right schema for your consumers.
- **Advanced B-Rep support.** Some viewers cannot render `IfcAdvancedBrep`; ship a tessellated fallback when interoperability matters more than exactness.
- **Performance.** Building geometry element-by-element without the iterator (and without caching repeated geometry) is slow on large models.
- **Validity.** Export invalid B-Reps and consumers choke; validate with `BRepCheck_Analyzer` before serialization ([[24_brep-boolean-occt]]).

## Practice

1. Open a sample IFC with IfcOpenShell, read the length unit, and print the world bounding box in metres.
2. Iterate all elements and report their type counts (walls, slabs, etc.) and triangle totals.
3. Take one element and recompute its world placement by composing the `IfcLocalPlacement` chain by hand; compare to the iterator's 4x4 matrix.
4. Export the same stone solid twice, once as a triangulated face set and once as an advanced B-Rep; compare file size and which viewers open each.
5. Author a minimal IFC4 file with one `IfcBuildingElementProxy` carrying a tessellated stone block at a non-trivial placement.

## References

- IfcOpenShell documentation, geometry processing and iterator: https://docs.ifcopenshell.org/
- buildingSMART IFC4 specification: https://standards.buildingsmart.org/IFC/RELEASE/IFC4/
- Thomas Krijnen, IfcOpenShell academic and technical write-ups.
- OpenCASCADE reference manual (B-Rep meshing for the tessellated path).

## Flashcard seeds

- Q: Build the rotation matrix $R$ of an `IfcAxis2Placement3D` from the Axis $\hat z$ and RefDirection $x$ (derive, step 1 of 2). <br> <b>Formula:</b> $x'=x-(x\cdot\hat z)\hat z,\ \ \hat x=x'/\lVert x'\rVert,\ \ \hat y=\hat z\times\hat x$ <br> <i>Hint:</i> first project $x$ onto the plane perpendicular to $\hat z$ (subtract the part along $\hat z$), then normalize. :: A: $\hat z$ is the unit Axis <br> remove any component of $x$ along $\hat z$: $x' = x-(x\cdot\hat z)\hat z$, then $\hat x = x'/\lVert x'\rVert$ <br> the third axis closes the right-handed frame: $\hat y=\hat z\times\hat x$ <br> $\boxed{R=[\,\hat x\ \ \hat y\ \ \hat z\,]}$ (columns)
- Q: From $R=[\hat x\ \hat y\ \hat z]$ and origin $t$, assemble the homogeneous local placement $T$ (derive, step 2 of 2). <br> <b>Formula:</b> $p_\text{parent}=Rp+t$, target $T=\begin{bmatrix}R & t\\ 0 & 1\end{bmatrix}$ <br> <i>Hint:</i> write the affine map $Rp+t$ as a single $4\times4$ times $(p,1)$. :: A: a local point $p$ maps as $p_\text{parent}=Rp+t$ <br> stack into homogeneous form so one matrix carries rotation and translation <br> $\boxed{T=\begin{bmatrix}R & t\\ 0 & 1\end{bmatrix}}$
- Q: Each frame's placement is relative to its parent. Write the world transform of an element (derive, step 1 of 2). <br> <b>Formula:</b> $T_\text{world}=T_\text{site}\,T_\text{building}\,T_\text{storey}\,T_\text{element}$ <br> <i>Hint:</i> start from $T_\text{site}\,T_\text{building}$ and left-multiply each child placement by the accumulated parent transform. :: A: composition runs parent-to-child: $T_\text{world\leftarrow building}=T_\text{site}\,T_\text{building}$ <br> append storey then element, each left-multiplied by its parent's accumulated transform <br> $\boxed{T_\text{world}=T_\text{site}\,T_\text{building}\,T_\text{storey}\,T_\text{element}}$
- Q: Given $T_\text{world}=T_\text{site}T_\text{building}T_\text{storey}T_\text{element}$, derive the world position of a local point $p$ (derive, step 2 of 2). <br> <b>Formula:</b> $\tilde p=(p,1),\ \ \tilde p_w=T_\text{world}\,\tilde p$ <br> <i>Hint:</i> lift $p$ to homogeneous coordinates, apply $T_\text{world}$ once, then drop the trailing 1. :: A: lift $p$ to homogeneous $\tilde p=(p,1)$ <br> apply the composed transform once: $\tilde p_w=T_\text{world}\,\tilde p$ <br> drop the homogeneous coordinate <br> $\boxed{p_w=R_\text{world}\,p+t_\text{world}}$
- Q: Coordinates are stored in model units; you need SI metres. Derive how the length-unit scale acts on a coordinate (derive, step 1 of 2). <br> <b>Formula:</b> $s=$ metres per model unit (mm $\Rightarrow s=10^{-3}$), $p_\text{SI}=s\,p$ <br> <i>Hint:</i> each stored value is a count of units, so multiply componentwise by $s$. :: A: let $s$ = metres per model unit (mm $\Rightarrow s=10^{-3}$) <br> each stored axis value is a count of units, so multiply: $x_\text{SI}=s\,x$ <br> $\boxed{p_\text{SI}=s\,p}$
- Q: Starting from $p_\text{SI}=s\,p$, derive the error factor from forgetting the unit scale (derive, step 2 of 2). <br> <b>Formula:</b> $p_\text{wrong}=p$, ratio $p_\text{wrong}/p_\text{SI}=1/s$ <br> <i>Hint:</i> set the applied scale to $1$ and take the ratio to the true $p_\text{SI}$; plug $s=10^{-3}$. :: A: skipping the scale uses $p$ directly as if metres: $p_\text{wrong}=p$ <br> ratio to truth: $p_\text{wrong}/p_\text{SI}=1/s$ <br> for mm, $1/s=1000$ <br> $\boxed{p_\text{wrong}=1000\,p_\text{SI}\ \text{(mm model)}}$
- Q: Derive the matrix of a $90^\circ$ rotation about $Z$ to plug into a building placement (derive). <br> <b>Formula:</b> $R_z(\theta)=\begin{bmatrix}\cos\theta & -\sin\theta & 0\\ \sin\theta & \cos\theta & 0\\ 0&0&1\end{bmatrix}$ <br> <i>Hint:</i> substitute $\theta=90^\circ$ so $\cos\theta=0,\ \sin\theta=1$. :: A: $R_z(\theta)=\begin{bmatrix}\cos\theta & -\sin\theta & 0\\ \sin\theta & \cos\theta & 0\\ 0&0&1\end{bmatrix}$ <br> set $\theta=90^\circ$: $\cos=0,\ \sin=1$ <br> $\boxed{R_z(90^\circ)=\begin{bmatrix}0&-1&0\\1&0&0\\0&0&1\end{bmatrix}}$, i.e. $(x,y,z)\to(-y,x,z)$
- Q: Block point $p=(1,0,0)$ m; storey translates by $(10,4,3)$; building applies $90^\circ$ about $Z$. Find $p_\text{world}$ (worked, m scale). <br> <b>Formula:</b> $p_1=p+(10,4,3)$, then $R_z(90^\circ):(x,y,z)\to(-y,x,z)$ <br> <i>Hint:</i> apply the storey translation first, then the building rotation. :: A: storey: $(1,0,0)+(10,4,3)=(11,4,3)$ <br> building $90^\circ$ Z: $(x,y,z)\to(-y,x,z)$ <br> $(11,4,3)\to(-4,11,3)$ <br> $\boxed{p_\text{world}=(-4,11,3)\ \text{m}}$
- Q: Same setup as the m-scale example, but the model length unit is millimetres and you forget to scale. Where does $p$ land? (worked, mm regime). <br> <b>Formula:</b> true answer $(-4,11,3)$ m; unscaled error factor $1/s=1000$ <br> <i>Hint:</i> multiply the correct metres answer by $1000$. :: A: true metres answer is $(-4,11,3)$ <br> unscaled treats mm-counts as metres, factor $1/s=1000$ <br> $(-4,11,3)\times1000$ <br> $\boxed{p_\text{wrong}=(-4000,11000,3000)\ \text{(off by 1000)}}$
- Q: A 2 mm scan feature is stored in a metre-unit IFC ($s=1$) but you wrongly assume mm and multiply by $s_\text{assumed}=10^{-3}$. What size do you report? (worked, shop scale). <br> <b>Formula:</b> $p_\text{reported}=p_\text{stored}\cdot s_\text{assumed}$, with $p_\text{stored}=2\times10^{-3}$ m <br> <i>Hint:</i> the stored value is already in metres; just multiply it by the wrong $10^{-3}$. :: A: stored value is already metres: $2\times10^{-3}$ m <br> wrong assumption multiplies by $s_\text{assumed}=10^{-3}$: $2\times10^{-3}\cdot10^{-3}$ <br> $\boxed{2\times10^{-6}\ \text{m}=2\,\mu\text{m}}$ (1000x too small — scale bugs cut both ways)
- Q: Axis $\hat z=(0,0,1)$, RefDirection $x=(2,0,0)$, origin $t=(5,5,0)$. Build the $4\times4$ placement $T$ (worked, 3D build). <br> <b>Formula:</b> $\hat x=x/\lVert x\rVert,\ \ \hat y=\hat z\times\hat x,\ \ T=\begin{bmatrix}R & t\\ 0 & 1\end{bmatrix}$ <br> <i>Hint:</i> normalize $x$ first (it is already perpendicular to $\hat z$), then form $\hat y$ by the cross product. :: A: normalize $\hat x=(2,0,0)/2=(1,0,0)$ <br> $\hat y=\hat z\times\hat x=(0,1,0)$ <br> $R=I_3$ <br> $\boxed{T=\begin{bmatrix}1&0&0&5\\0&1&0&5\\0&0&1&0\\0&0&0&1\end{bmatrix}}$
- Q: Axis $\hat z=(0,0,1)$, RefDirection $x=(1,1,0)$. Find $\hat x,\hat y$ (worked, non-axis-aligned X). <br> <b>Formula:</b> $x'=x-(x\cdot\hat z)\hat z,\ \ \hat x=x'/\lVert x'\rVert,\ \ \hat y=\hat z\times\hat x$ <br> <i>Hint:</i> check $x\cdot\hat z$ first (it is $0$ here), then divide $x$ by its length $\sqrt2$. :: A: $x\cdot\hat z=0$, so $x'=x$ <br> $\hat x=(1,1,0)/\sqrt2=(\tfrac{1}{\sqrt2},\tfrac{1}{\sqrt2},0)$ <br> $\hat y=\hat z\times\hat x=(-\tfrac{1}{\sqrt2},\tfrac{1}{\sqrt2},0)$ <br> $\boxed{\hat x=(\tfrac{1}{\sqrt2},\tfrac{1}{\sqrt2},0),\ \hat y=(-\tfrac{1}{\sqrt2},\tfrac{1}{\sqrt2},0)}$
- Q: Check that $\hat x=(\tfrac{1}{\sqrt2},\tfrac{1}{\sqrt2},0)$, $\hat y=(-\tfrac{1}{\sqrt2},\tfrac{1}{\sqrt2},0)$, $\hat z=(0,0,1)$ form a valid rotation (closing test). <br> <b>Formula:</b> $R\in SO(3)\iff \hat x\cdot\hat y=0,\ \lVert\hat x\rVert=1,\ \det[\hat x\,\hat y\,\hat z]=+1$ <br> <i>Hint:</i> a pass means all three hold; compute the dot product, one norm, and the determinant. :: A: $\hat x\cdot\hat y=-\tfrac12+\tfrac12+0=0$ (orthogonal) <br> $\lVert\hat x\rVert^2=\tfrac12+\tfrac12=1$ (unit) <br> $\det[\hat x\,\hat y\,\hat z]=+1$ (right-handed) <br> $\boxed{R\in SO(3)\ \checkmark}$
- Q: World point $(-4,11,3)$ m. Invert the chain to recover $p_\text{local}$ and confirm it is $(1,0,0)$ (closing test). <br> <b>Formula:</b> inverse rotation $(x,y,z)\to(y,-x,z)$, then subtract storey shift $(10,4,3)$ <br> <i>Hint:</i> undo in reverse order — rotation first, then the translation; a pass returns exactly $(1,0,0)$. :: A: inverse rotation $(x,y,z)\to(y,-x,z)$: $(-4,11,3)\to(11,4,3)$ <br> subtract storey shift $(10,4,3)$: $(1,0,0)$ <br> matches the input <br> $\boxed{p_\text{local}=(1,0,0)\ \checkmark}$
- Q: A wall has width $200$ in an IFC whose length unit is millimetre. Convert to metres and sanity-check it (closing test). <br> <b>Formula:</b> $w_\text{SI}=s\,w$ with $s=10^{-3}$ m/unit <br> <i>Hint:</i> multiply by $10^{-3}$; a pass is a plausible wall thickness (tens of cm), not metres. :: A: $s=10^{-3}$ m/unit <br> $200\times10^{-3}=0.2$ m <br> $0.2$ m $=20$ cm, a sensible wall thickness <br> $\boxed{0.2\ \text{m}\ \checkmark}$ (had we read it as metres, $200$ m is absurd — the check catches it)
- CLOZE: An `IfcAxis2Placement3D` is defined by an origin, a {{c1::Z axis}} (Axis) and an {{c2::X axis}} (RefDirection); the remaining axis is rebuilt as $Y=Z\times X$.
- CLOZE: The world transform is the product of relative placements up the chain, $T_\text{world}=T_\text{site}\,T_\text{building}\,T_\text{storey}\,T_\text{element}$, composed {{c1::parent-to-child}} (each parent left-multiplies).
- CLOZE: To reach SI metres, multiply every coordinate by the length-unit scale $s$; for a millimetre model $s={{c1::10^{-3}}}$, and forgetting it scales coordinates by {{c2::1000}}.
- CLOZE: IfcOpenShell's `create_shape` can return triangulated geometry or an {{c1::OpenCASCADE}} B-Rep (`TopoDS_Shape`).
- Q: When do you export tessellation vs an advanced B-Rep? <br> <i>Hint:</i> weigh universal viewer support against exact curvature. :: A: Tessellated (`IfcTriangulatedFaceSet`) when universal viewer support matters; `IfcAdvancedBrep` when exact curvature must be preserved and consumers support IFC4. <br> Ship a tessellated fallback if interoperability outweighs exactness.
- Q: What does IFC add beyond a mesh? <br> <i>Hint:</i> think semantics, not just triangles. :: A: Semantics — named building elements, materials, placements, and parent/containment relationships, i.e. an entity graph, not just loose triangles.
- Q: Why must you compose the whole placement chain rather than read one element's frame? <br> <i>Hint:</i> recall what `IfcLocalPlacement` is relative to. :: A: Each `IfcLocalPlacement` is relative to its parent, so a single frame is incomplete; compose the full chain (or request `use-world-coords` to bake placements into vertices).
- Q: Why validate B-Reps before IFC export? <br> <i>Hint:</i> think about what a consumer does with invalid topology. :: A: Invalid topology or out-of-tolerance edges make downstream consumers choke; check with `BRepCheck_Analyzer` first.
- Q: When do you use IfcOpenShell's geometry iterator? <br> <i>Hint:</i> consider model size and repeated geometry. :: A: On large models — it gives multicore processing, caching, and reuse of repeated geometry, returning per-element vertex/face arrays plus the 4x4 placement.
