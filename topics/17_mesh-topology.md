# Mesh data structures and topology

Prereqs: [[01_vectors]] => this => [[18_discrete-curvature]], [[19_laplacian-poisson]]

## Why it matters for the stone scan->CAD->IFC / CADFit workflow

A raw stone scan arrives as a point cloud, but every downstream step (remeshing, curvature estimation, segmentation, B-Rep fitting, IFC export) needs a *connected* triangle mesh with known neighbor relations, not a bag of triangles. Topology is the part of the mesh that survives bending and stretching: it tells you which faces share an edge, where the boundary is, and whether the surface is closed, and that connectivity is exactly what feeds the cotangent Laplacian in [[19_laplacian-poisson]] and the angle defect in [[18_discrete-curvature]]. If the connectivity is wrong (duplicated vertices that should be welded, non-manifold edges shared by three faces, inconsistent winding), every later operator silently produces garbage: the Laplacian is asymmetric, normals flip, and OpenCASCADE refuses to sew the shells into a solid. The Euler characteristic $\chi = V - E + F$ is a one-integer health check you can run before committing compute: a clean scanned stone slab should be a topological disk ($\chi=1$) and a closed scanned boulder a topological sphere ($\chi=2$), so any other value means holes, handles, or non-manifold junk to repair. In a CADFit-style synthesis layer that emits CAD programs from geometry, the program must reference *faces and their adjacency*, so the half-edge structure is the substrate the synthesizer reads.

## Intuition

Think of a mesh as a city map, not a photograph. The photograph (a rendered picture) only shows you colors and silhouettes. The map encodes which streets connect to which intersections, and that connection data is what lets you compute routes. A mesh is the map: vertices are intersections, edges are streets, faces are city blocks, and the topology is the street-connection table.

The key mental shift for a Grasshopper user is this: a mesh is *data*, not a picture. Rhino draws it, but the value lives in the index arrays that say "face 7 uses vertices 3, 9, 12" and "edge 3-9 is also used by face 4". The half-edge structure is just a very disciplined version of that connection table. It splits every undirected edge into two oppositely pointing arrows (half-edges), one for each side. Each arrow remembers three things: the vertex it starts from, the next arrow going counterclockwise around its face, and its twin (the arrow on the other side of the same edge). With only "next" and "twin", you can walk the entire neighborhood of any vertex by spinning around it like a turnstile.

## The math (first principles)

### Mesh as a graph plus faces

A triangle mesh is a triple $M = (V, E, F)$.

- $V = \{1, \dots, n\}$ is the vertex set. Geometry attaches a position $\mathbf{p}_i \in \mathbb{R}^3$ to each vertex.
- $E \subseteq V \times V$ is the edge set (unordered pairs).
- $F$ is the face set. For a triangle mesh each face is an ordered triple $(i, j, k)$.

Separate two notions. **Connectivity (topology)** is $(V, E, F)$ as a combinatorial object: who is adjacent to whom. **Geometry** is the embedding $\mathbf{p}: V \to \mathbb{R}^3$. Topology answers "is this surface closed?"; geometry answers "how curved is it?". This document is about topology; the positions only enter in the worked example and code as labels.

Convention: faces are stored **counterclockwise (CCW)** when viewed from the outside (the side the outward normal points toward). For face $(i,j,k)$ the outward normal is
$$\mathbf{n} = \frac{(\mathbf{p}_j - \mathbf{p}_i) \times (\mathbf{p}_k - \mathbf{p}_i)}{\lVert (\mathbf{p}_j - \mathbf{p}_i) \times (\mathbf{p}_k - \mathbf{p}_i) \rVert}.$$
This is the right-hand-rule convention shared by libigl, OpenMesh, and Rhino. Consistent winding across all faces is what *orientation* means.

### The half-edge data structure

An undirected edge $\{i,j\}$ is split into two directed **half-edges**: $h$ going $i \to j$ and its twin $\bar{h}$ going $j \to i$. Each half-edge $h$ stores pointers:

- $\mathrm{vertex}(h)$: the vertex at the *base* (tail) of $h$.
- $\mathrm{twin}(h)$: the oppositely oriented half-edge on the same edge.
- $\mathrm{next}(h)$: the next half-edge CCW around the same face.
- $\mathrm{face}(h)$: the face to the left of $h$ (the face $h$ circulates).
- $\mathrm{edge}(h)$: the undirected edge $h$ lies on.

Each vertex, edge, and face stores one representative half-edge. The two combinatorial axioms are:
$$\mathrm{twin}(\mathrm{twin}(h)) = h, \qquad \mathrm{twin}(h) \neq h.$$
The first says twin is an involution. The second forbids a half-edge from being its own twin. Walking $\mathrm{next}$ repeatedly cycles around a face; for a triangle, $\mathrm{next}^3(h) = h$.

The **previous** half-edge and the rotation around a vertex are derived, not stored:
$$\mathrm{prev}(h) = \mathrm{next}(\mathrm{next}(h)) \ \text{(triangle)}, \qquad \mathrm{rotateCW}(h) = \mathrm{next}(\mathrm{twin}(h)).$$
Applying $\mathrm{rotateCW}$ from an outgoing half-edge of vertex $v$ visits all outgoing half-edges of $v$ in order. That circulation *is* the one-ring traversal.

### One-ring neighborhood

The **one-ring** of vertex $i$ is its set of immediate neighbors
$$N(i) = \{\, j \in V : \{i,j\} \in E \,\}.$$
The one-ring is the support of every local mesh operator. The cotangent Laplacian (next topics) sums only over $N(i)$:
$$(\Delta f)_i = \frac{1}{2A_i} \sum_{j \in N(i)} (\cot \alpha_{ij} + \cot \beta_{ij})\,(f_j - f_i).$$
You cannot evaluate this without an efficient way to enumerate $N(i)$, which is precisely what $\mathrm{rotateCW}$ around $i$ gives you in time proportional to the vertex degree, not to mesh size.

### Manifoldness, boundary, orientation

A mesh is **edge-manifold** if every edge is incident to exactly one face (a boundary edge) or exactly two faces (an interior edge). Three or more faces meeting at one edge is non-manifold.

A mesh is **vertex-manifold** if the faces around every vertex form a single fan (interior vertex) or a single open fan (boundary vertex). Two cones touching at a tip share a vertex but are not vertex-manifold (a "bowtie").

A **boundary edge** is a half-edge whose twin has no face (or, equivalently, an edge used by only one face). Boundary half-edges chain together into closed **boundary loops**. A closed surface has no boundary.

**Orientability** means winding can be made consistent so that twin half-edges always run in opposite directions. A Mobius band is the standard non-orientable counterexample. Scanned solids are orientable; a careless triangle-soup import may not be until you fix winding.

### Euler characteristic and genus

For any mesh, define the **Euler characteristic**
$$\chi = V - E + F.$$
This is a topological invariant: it does not change under any deformation, subdivision, or remesh that preserves connectivity type. For a **closed orientable** surface of **genus** $g$ (number of handles or "holes through" the solid),
$$\chi = 2 - 2g \quad\Longrightarrow\quad g = \frac{2 - \chi}{2} = \frac{2 - (V - E + F)}{2}.$$
Examples: sphere $g=0,\ \chi=2$; torus (one handle) $g=1,\ \chi=0$; double torus $g=2,\ \chi=-2$. For a surface with $b$ boundary loops the generalized relation is
$$\chi = 2 - 2g - b,$$
so a topological disk (a flat slab patch) has $g=0,\ b=1,\ \chi = 1$.

**Triangle-mesh edge count.** For a closed triangle mesh, every face has 3 edges and every edge is shared by 2 faces, so counting face-edge incidences two ways gives $3F = 2E$, hence
$$E = \tfrac{3}{2}F, \qquad F = 2V - 4 \ \ (\text{closed, genus } 0).$$
Derivation of the last: substitute $E = \tfrac{3}{2}F$ into $\chi = V - E + F = 2$ to get $V - \tfrac{1}{2}F = 2$, so $F = 2V - 4$ and $E = 3V - 6$. These ratios are a fast sanity check on mesh size and on triangulation quality.

**Discrete Gauss-Bonnet (the bridge to curvature).** The **angle defect** at an interior vertex $i$ is
$$d_i = 2\pi - \sum_{f \ni i} \theta_i^f,$$
where $\theta_i^f$ is the interior angle of face $f$ at vertex $i$. Summing over all vertices of a closed mesh recovers topology:
$$\sum_{i \in V} d_i = 2\pi\,\chi(M).$$
This is the discrete Gauss-Bonnet theorem. Total curvature is fixed by topology alone, which is why $\chi$ matters for the curvature topic that this one unlocks.

## Worked example

Take a regular **tetrahedron**: 4 triangular faces, 4 vertices. It is a closed mesh, topologically a sphere.

Count by hand.
- $V = 4$.
- $F = 4$.
- Edges: each face contributes 3 edge-incidences, total $3 \times 4 = 12$, each interior edge shared by 2 faces, so $E = 12 / 2 = 6$. List them to confirm: with vertices $\{1,2,3,4\}$ the edges are $12,13,14,23,24,34$, which is $6$.

Euler characteristic:
$$\chi = V - E + F = 4 - 6 + 4 = 2. \checkmark$$
Genus:
$$g = \frac{2 - \chi}{2} = \frac{2 - 2}{2} = 0. \checkmark$$
So the tetrahedron is a genus-0 closed surface (a sphere), as expected. Check the closed-mesh identities: $F = 2V - 4 = 2(4) - 4 = 4 \checkmark$ and $E = 3V - 6 = 3(4) - 6 = 6 \checkmark$.

**Angle defect.** Each face is equilateral, so every interior angle is $\theta = 60^\circ = \pi/3$. Each vertex of a tetrahedron touches 3 faces, so the angle sum at a vertex is $3 \times \pi/3 = \pi$. The defect at each vertex is
$$d_i = 2\pi - \pi = \pi.$$
Sum over all 4 vertices:
$$\sum_i d_i = 4\pi.$$
Gauss-Bonnet predicts $2\pi\chi = 2\pi(2) = 4\pi$. They match. The total curvature equals $4\pi$ for *any* closed genus-0 mesh, tetrahedron or fine-scanned sphere alike, because it is set by topology.

**Now open it.** Delete one face (say $234$). $V=4,\ E=6,\ F=3$, so $\chi = 4 - 6 + 3 = 1$. With $g=0$ and one boundary loop, $\chi = 2 - 0 - 1 = 1 \checkmark$. The deleted face exposed boundary edges $23,24,34$, which form one closed boundary loop. This is the disk case, exactly what a scanned slab patch looks like.

## Code

C++ with **libigl** (header-only, Eigen-based). This is the fastest way to compute topology invariants and run manifold checks in a scan pipeline. Mesh is `V` (n x 3 positions) and `F` (m x 3 vertex indices, CCW).

```cpp
#include <igl/euler_characteristic.h>
#include <igl/is_edge_manifold.h>
#include <igl/is_vertex_manifold.h>
#include <igl/boundary_loop.h>
#include <igl/edges.h>
#include <Eigen/Core>
#include <iostream>
#include <vector>

int main() {
    // A tetrahedron: 4 vertices (V), 4 CCW triangle faces (F).
    Eigen::MatrixXd V(4, 3);
    V << 0, 0, 0,
         1, 0, 0,
         0, 1, 0,
         0, 0, 1;
    Eigen::MatrixXi F(4, 3);
    F << 0, 2, 1,   // each row is a face (i,j,k), wound CCW (outward normal)
         0, 1, 3,
         0, 3, 2,
         1, 2, 3;

    // chi = V - E + F, computed directly from connectivity.
    int chi = igl::euler_characteristic(F);          // returns 2 here
    std::cout << "chi = " << chi << "\n";
    int genus = (2 - chi) / 2;                        // g = (2 - chi)/2 (closed)
    std::cout << "genus = " << genus << "\n";

    // E: count unique undirected edges to confirm 3F = 2E for a closed mesh.
    Eigen::MatrixXi E;
    igl::edges(F, E);                                 // one row per undirected edge
    std::cout << "V=" << V.rows() << " E=" << E.rows()
              << " F=" << F.rows() << "\n";           // 4, 6, 4

    // Manifold health checks BEFORE any Laplacian / fitting step.
    bool edge_ok = igl::is_edge_manifold(F);          // every edge has 1 or 2 faces
    Eigen::VectorXi vmark;                             // per-vertex manifold flag
    bool vert_ok = igl::is_vertex_manifold(F, vmark); // single fan per vertex
    std::cout << "edge_manifold=" << edge_ok
              << " vertex_manifold=" << vert_ok << "\n";

    // Boundary: empty loops => closed surface. (Tetrahedron is closed.)
    std::vector<std::vector<int>> loops;
    igl::boundary_loop(F, loops);                     // ordered boundary vertex loops
    std::cout << "boundary_loops = " << loops.size() << "\n"; // 0

    return 0;
}
```

Symbol map: `igl::euler_characteristic(F)` is $\chi = V - E + F$; `igl::edges` materializes $E$; `is_edge_manifold` enforces the "1 or 2 faces per edge" rule; `boundary_loop` returns the boundary half-edge chains as vertex loops. Run the same code on a scanned mesh: `chi != 2` plus nonzero boundary loops flags holes to fill before fitting.

**One-ring traversal with OpenMesh** (explicit half-edge library). Use this when you need to *walk* connectivity, not just count it.

```cpp
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <iostream>
typedef OpenMesh::TriMesh_ArrayKernelT<> Mesh;

int main() {
    Mesh mesh;
    // ... mesh.add_vertex(...) / mesh.add_face(...) or OpenMesh::IO::read_mesh ...

    for (auto vh : mesh.vertices()) {                 // every vertex i
        std::cout << "vertex " << vh.idx() << " one-ring:";
        // vv_iter = vertex-vertex circulator = the one-ring N(i).
        for (auto vv = mesh.vv_iter(vh); vv.is_valid(); ++vv)
            std::cout << " " << vv->idx();            // each j in N(i)
        std::cout << "  (degree = " << mesh.valence(vh) << ")\n";

        // Boundary test via half-edge: a boundary vertex has a halfedge
        // whose face handle is invalid (twin has no face).
        if (mesh.is_boundary(vh))
            std::cout << "  (on boundary)\n";
    }
    return 0;
}
```

`vv_iter(vh)` circulates the one-ring by chasing $\mathrm{rotateCW}(h) = \mathrm{next}(\mathrm{twin}(h))$ internally. `is_boundary` checks for a half-edge with no incident face, the operational definition of a boundary edge.

A short Python check (libigl bindings) clarifies the invariant in one line:

```python
import igl, numpy as np
F = np.array([[0,2,1],[0,1,3],[0,3,2],[1,2,3]])   # tetrahedron faces
print(igl.euler_characteristic(F))                 # -> 2  (chi)
print((2 - igl.euler_characteristic(F)) // 2)      # -> 0  (genus)
```

## Connections

- [[01_vectors]]: vertex positions are vectors in $\mathbb{R}^3$; face normals come from the cross product of edge vectors, and CCW winding is fixed by the right-hand rule. Topology indexes vectors; without the vector layer there is no geometry to attach.
- [[18_discrete-curvature]]: the angle defect $d_i = 2\pi - \sum \theta_i^f$ is computed by circulating the one-ring, and its global sum is locked to $\chi$ by discrete Gauss-Bonnet. This topic supplies the neighborhood traversal and the invariant that topic consumes.
- [[19_laplacian-poisson]]: the cotangent Laplacian is a sum over $N(i)$ assembled by walking half-edges; manifoldness and consistent orientation are what make the matrix symmetric and the Poisson solve well-posed.
- [[03_eig-svd-pca]]: per-vertex normals and local frames are often estimated by PCA over the one-ring, so the neighborhood structure defined here feeds the eigen-decomposition there.
- [[05_transforms-se3]]: aligning two scans before merging meshes is an $\mathrm{SE}(3)$ transform; topology must be re-stitched (welded, re-indexed) after the rigid alignment of [[05_transforms-se3]].

## Pitfalls and failure modes

- **Duplicated (unwelded) vertices.** OBJ/STL exporters often write each triangle independently, so coincident corners are distinct indices. The mesh then looks watertight but has $\chi$ far from the true value and a non-manifold structure. Weld by spatial hashing with a tolerance ($\sim 10^{-6}$ of the bounding-box diagonal) before computing topology. Too large a tolerance collapses thin features; too small leaves cracks.
- **Inconsistent winding.** A triangle-soup import may have neighboring faces wound in opposite directions. Twin half-edges then point the *same* way, the structure is unorientable in software even when the surface is orientable, and normals flip face to face. Run an orientation pass (BFS flood from a seed face, flipping neighbors to agree) before anything else.
- **Non-manifold edges and bowtie vertices.** Three faces on one edge, or two fans touching at one vertex, break the half-edge invariant $\mathrm{twin}(\mathrm{twin}(h)) = h$. Detect with `is_edge_manifold` and `is_vertex_manifold`; many half-edge builders simply refuse to load such input. Scanned stone with self-contact is a common source.
- **Confusing genus from $\chi$ on an open mesh.** $g = (2-\chi)/2$ is only valid for closed surfaces. On a slab patch with boundary you must use $\chi = 2 - 2g - b$ and count boundary loops, or you will report a fractional or negative "genus" that is meaningless.
- **Holes vs handles.** Both lower $\chi$, but a hole (missing boundary patch from incomplete scan coverage) is repaired by hole-filling, while a handle (a genuine tunnel through the stone) is real geometry. Check boundary loops: a hole has a boundary loop, a handle does not. Misclassifying one as the other either destroys a real tunnel or leaves a gap unfilled.
- **Float tolerance leaking into topology.** Topology is combinatorial and must be decided by *integer* index equality after welding, never by repeated floating-point distance tests downstream. Decide adjacency once, store indices, and stop comparing coordinates.
- **Degenerate faces.** Zero-area triangles (collinear vertices) have undefined normals and infinite cotangents later. Remove them during cleanup; they pass a naive $\chi$ count but poison every geometric operator.

## Practice

1. **Hand count.** Take a cube drawn as 6 quad faces, then as 12 triangles (each quad split). Compute $V, E, F$ and $\chi$ for both. Confirm $\chi = 2$ in both cases and that $3F = 2E$ holds only for the triangulated version. Explain why $\chi$ is unchanged by triangulation.
2. **Genus from a count.** A closed triangle mesh reports $V = 1000$, $F = 2000$. Compute $E$, then $\chi$, then $g$. (Use $3F = 2E$.) State what kind of solid this could be.
3. **One-ring code.** Load any closed mesh in libigl or OpenMesh and print, for 5 vertices, their one-ring neighbor lists and degrees. Verify that $\sum_i \deg(i) = 2E$ (handshake lemma).
4. **Gauss-Bonnet check.** For an icosahedron (20 faces, 12 vertices), compute $\chi$, then sum the angle defects numerically and confirm the total equals $2\pi\chi$. The angle defect per vertex should be $2\pi - 5(\pi/3)$.
5. **Boundary detection.** Delete a strip of faces from a closed scanned mesh, recompute $\chi$ and the boundary loops, and confirm $\chi$ rose by the number of new boundary loops minus handles, consistent with $\chi = 2 - 2g - b$.
6. **Repair decision.** Write a function that, given a mesh, reports: edge-manifold? vertex-manifold? number of boundary loops? estimated genus? Then decide programmatically whether the mesh is "ready to fit" (closed, manifold, genus 0) or needs repair, and print the reason.

## References

- Keenan Crane, "Discrete Differential Geometry: An Applied Introduction" (CMU course notes), section 2.5 on the half-edge mesh and the curvature/Gauss-Bonnet chapters. https://www.cs.cmu.edu/~kmcrane/Projects/DDG/paper.pdf
- CS184/284A (UC Berkeley), "An Introduction to the Half-Edge Data Structure." https://cs184.eecs.berkeley.edu/sp24/docs/half-edge-intro
- libigl tutorial and headers: `euler_characteristic`, `is_edge_manifold`, `is_vertex_manifold`, `boundary_loop`, `edges`. https://libigl.github.io/ and https://github.com/libigl/libigl
- OpenMesh documentation, "Mesh Iterators and Circulators" (vertex-vertex circulator `vv_iter`, `is_boundary`). https://www.openmesh.org/media/Documentations/OpenMesh-6.3-Documentation/a00015.html
- Geogram geometric-algorithm library (Delaunay, mesh, remeshing). https://github.com/BrunoLevy/geogram
- Wolfram MathWorld, "Euler Characteristic" ($\chi = 2 - 2g$ for closed orientable surfaces). https://mathworld.wolfram.com/EulerCharacteristic.html
- Sohini Upadhyay, "Gauss-Bonnet for Discrete Surfaces" (UChicago REU 2015), $\sum_v d_v = 2\pi\chi$. http://math.uchicago.edu/~may/REU2015/REUPapers/Upadhyay.pdf
- OpenCASCADE reference manual, `BRepMesh_IncrementalMesh` (triangulating B-Rep shapes for exchange). https://dev.opencascade.org/doc/overview/html/
- IfcOpenShell documentation (IFC geometry serialization). https://docs.ifcopenshell.org/
- CADFit, CAD-program synthesis from geometry. https://arxiv.org/abs/2605.01171 and https://github.com/ghadinehme/CADFit

## Flashcard seeds

- Q: Derive the edge count $E$ for a closed triangle mesh by counting face-edge incidences two ways. <br> <b>Formula:</b> $3F = 2E$ (each face has 3 edges; each interior edge is shared by 2 faces) <br> <i>Hint:</i> Write the total incidences as $3F$ from the face side and $2E$ from the edge side, then equate. :: A: Each triangular face has 3 edges, so total face-edge incidences $= 3F$ <br> Every interior edge is shared by exactly 2 faces, so the same incidences $= 2E$ <br> Equate the two counts: $3F = 2E$ <br> $\boxed{E = \tfrac{3}{2}F}$
- Q: For a closed genus-0 triangle mesh, derive $F$ and $E$ in terms of $V$. <br> <b>Formula:</b> $E=\tfrac32 F$ and $\chi=V-E+F=2$; target $F=2V-4,\ E=3V-6$ <br> <i>Hint:</i> Substitute $E=\tfrac32 F$ into the Euler equation first, then solve for $F$. :: A: Substitute $E=\tfrac32 F$ into $\chi=V-E+F=2$: $V-\tfrac32 F+F=2$ <br> $V-\tfrac12 F=2 \Rightarrow F=2V-4$ <br> Then $E=\tfrac32 F=\tfrac32(2V-4)=3V-6$ <br> $\boxed{F=2V-4,\quad E=3V-6}$
- Q: Derive the genus $g$ as a function of $\chi$ for a closed orientable surface. <br> <b>Formula:</b> $\chi = 2 - 2g$; target $g=\tfrac{2-\chi}{2}$ <br> <i>Hint:</i> Isolate the $2g$ term, then divide by 2. :: A: Start from the closed-surface relation $\chi = 2 - 2g$ <br> Solve for $g$: $2g = 2 - \chi$ <br> $\boxed{g = \dfrac{2 - \chi}{2}}$
- Q: A tetrahedron has $V=4$, $F=4$ and is closed. Compute $\chi$. <br> <b>Formula:</b> $E=\tfrac32 F$, then $\chi=V-E+F$ <br> <i>Hint:</i> Get $E$ from $\tfrac32 F$ first, then plug $V,E,F$ into the Euler sum. :: A: First get $E$ from $E=\tfrac32 F=\tfrac32(4)=6$ <br> Apply $\chi=V-E+F=4-6+4$ <br> $\boxed{\chi=2}$
- Q: The tetrahedron has $\chi=2$. Compute its genus. <br> <b>Formula:</b> $g=\tfrac{2-\chi}{2}$ <br> <i>Hint:</i> Substitute $\chi=2$ directly. :: A: Use $g=(2-\chi)/2$ with $\chi=2$ <br> $g=(2-2)/2=0/2$ <br> $\boxed{g=0}$ (a topological sphere)
- Q: A regular tetrahedron has equilateral faces with 3 faces meeting at each vertex. Compute the angle defect $d_i$ at one vertex. <br> <b>Formula:</b> $d_i = 2\pi - \sum_{f\ni i}\theta_i^f$, with each $\theta_i^f=\tfrac{\pi}{3}$ <br> <i>Hint:</i> Sum the three $\tfrac{\pi}{3}$ interior angles first, then subtract from $2\pi$. :: A: Each interior angle of an equilateral triangle is $\theta=\pi/3$ <br> The angle sum at the vertex is $\sum_{f\ni i}\theta_i^f = 3\cdot\tfrac{\pi}{3}=\pi$ <br> Defect $d_i=2\pi-\sum\theta_i^f=2\pi-\pi$ <br> $\boxed{d_i=\pi}$
- Q: Each of the 4 tetrahedron vertices has defect $d_i=\pi$ and $\chi=2$. Find the total curvature and check Gauss-Bonnet. <br> <b>Formula:</b> $\sum_i d_i = 2\pi\chi$ <br> <i>Hint:</i> Add the four defects, then compare with $2\pi\chi$. :: A: Sum over all 4 vertices: $\sum_i d_i = 4\cdot\pi=4\pi$ <br> Gauss-Bonnet predicts $2\pi\chi=2\pi(2)=4\pi$ <br> They agree, confirming total curvature is set by topology <br> $\boxed{\sum_i d_i = 2\pi\chi = 4\pi}$
- Q: A regular icosahedron has 5 equilateral faces meeting at each vertex. Derive the angle defect per vertex. <br> <b>Formula:</b> $d_i = 2\pi - \sum_{f\ni i}\theta_i^f$, with each $\theta_i^f=\tfrac{\pi}{3}$ <br> <i>Hint:</i> Compute the angle sum $5\cdot\tfrac{\pi}{3}$, then subtract from $2\pi$ over a common denominator. :: A: Interior angle of an equilateral face $\theta=\pi/3$ <br> Five faces meet, so angle sum $=5\cdot\tfrac{\pi}{3}=\tfrac{5\pi}{3}$ <br> $d_i=2\pi-\tfrac{5\pi}{3}=\tfrac{6\pi-5\pi}{3}$ <br> $\boxed{d_i=\tfrac{\pi}{3}}$
- Q: Delete one face of the tetrahedron (open/disk case), leaving $V=4,\ E=6,\ F=3$. Compute $\chi$. <br> <b>Formula:</b> $\chi=V-E+F$; cross-check with $\chi=2-2g-b$ <br> <i>Hint:</i> Plug $V,E,F$ into the Euler sum, then verify with $g=0,\ b=1$. :: A: Deleting a face removes a face but keeps all vertices and edges: $V=4,\ E=6,\ F=3$ <br> $\chi=V-E+F=4-6+3$ <br> Cross-check with $\chi=2-2g-b$ for $g=0$, one boundary loop $b=1$: $2-0-1=1$ <br> $\boxed{\chi=1}$ (a topological disk)
- Q: Generalize the closed-surface Euler relation to surfaces with boundary. <br> <b>Formula:</b> closed: $\chi=2-2g$; target with boundary: $\chi=2-2g-b$ <br> <i>Hint:</i> Each open boundary loop is one missing disk-cap; a disk contributes $\chi=1$, so removing $b$ caps subtracts $b$. :: A: A closed genus-$g$ surface has $\chi=2-2g$ <br> Each open boundary loop is a missing cap; removing one disk lowers $\chi$ by 1 (a disk has $\chi=1$) <br> Removing $b$ caps subtracts $b$ <br> $\boxed{\chi=2-2g-b}$
- Q: Worked example (large mesh). A closed triangle mesh reports $V=1000,\ F=2000$. Find $E$, $\chi$, $g$. <br> <b>Formula:</b> $E=\tfrac32 F$, $\chi=V-E+F$, $g=\tfrac{2-\chi}{2}$ <br> <i>Hint:</i> Compute $E$ from $\tfrac32 F$ first, then $\chi$, then $g$. :: A: $E=\tfrac32 F=\tfrac32(2000)=3000$ <br> $\chi=V-E+F=1000-3000+2000=0$ <br> $g=(2-\chi)/2=(2-0)/2=1$ <br> $\boxed{E=3000,\ \chi=0,\ g=1}$ (a one-handle solid / torus)
- Q: Worked example (small mesh, octahedron). $V=6,\ F=8$, closed. Find $E$, $\chi$, $g$. <br> <b>Formula:</b> $E=\tfrac32 F$, $\chi=V-E+F$, $g=\tfrac{2-\chi}{2}$ <br> <i>Hint:</i> Get $E=\tfrac32(8)$ first, then the Euler sum, then the genus. :: A: $E=\tfrac32 F=\tfrac32(8)=12$ <br> $\chi=V-E+F=6-12+8=2$ <br> $g=(2-2)/2=0$ <br> $\boxed{E=12,\ \chi=2,\ g=0}$ (a sphere)
- Q: Worked example (open slab patch). A scanned slab triangulates to $V=50,\ E=129,\ F=80$ with one boundary loop. Find $\chi$ and the implied genus. <br> <b>Formula:</b> $\chi=V-E+F$, then open-mesh $\chi=2-2g-b$ with $b=1$ <br> <i>Hint:</i> Compute $\chi$ from the counts, then solve $1=2-2g-1$ for $g$. :: A: $\chi=V-E+F=50-129+80=1$ <br> Open mesh, so use $\chi=2-2g-b$ with $b=1$: $1=2-2g-1$ <br> $2g=0 \Rightarrow g=0$ <br> $\boxed{\chi=1,\ g=0}$ (a flat disk, as expected for a slab)
- Q: Worked example (regime check, fine sphere scan). A scanned boulder mesh has $V=100000$ and is genus 0. Predict $F$, $E$, and the total curvature. <br> <b>Formula:</b> $F=2V-4$, $E=3V-6$, $\sum d_i=2\pi\chi$ with $\chi=2$ <br> <i>Hint:</i> Plug $V$ into the closed genus-0 size formulas; note the curvature does not depend on $V$. :: A: Closed genus-0: $F=2V-4=2(100000)-4=199996$ <br> $E=3V-6=3(100000)-6=299994$ <br> Total curvature $\sum d_i=2\pi\chi=2\pi(2)=4\pi$ regardless of resolution <br> $\boxed{F=199996,\ E=299994,\ \sum d_i=4\pi}$
- Q: Test yourself (close the loop). You computed $g=1$ for $V=1000,\ F=2000$. Verify it. <br> <b>Formula:</b> $\chi=2-2g$ should match $\chi=V-E+F$ (with $E=3000$) <br> <i>Hint:</i> A pass means both routes give the same $\chi$. Compute each side and compare. :: A: Plug $g=1$ into $\chi=2-2g=2-2(1)=0$ <br> Independently $\chi=V-E+F=1000-3000+2000=0$ <br> Both give $\chi=0$, so the genus-1 (torus) classification is consistent <br> $\boxed{\chi=0 \text{ both ways} \Rightarrow g=1 \text{ verified}}$
- Q: Test yourself (handshake lemma). For the tetrahedron with $E=6$, verify the degree sum. <br> <b>Formula:</b> $\sum_i \deg(i)=2E$ <br> <i>Hint:</i> A pass means the degree sum equals $2E$. Each vertex connects to the other 3, so add up the degrees. :: A: Each of the 4 vertices connects to the other 3, so $\deg(i)=3$ for all $i$ <br> $\sum_i \deg(i)=4\cdot 3=12$ <br> Check: $2E=2\cdot 6=12$ <br> $\boxed{\sum_i\deg(i)=2E=12 \checkmark}$
- Q: Test yourself (limit check on Gauss-Bonnet). For a closed genus-0 mesh, confirm $\sum d_i$ is dimensionless and equals $4\pi$ regardless of $V$. <br> <b>Formula:</b> $d_i=2\pi-\sum\theta_i^f$ and $\sum_i d_i=2\pi\chi$ with $\chi=2-2g$ <br> <i>Hint:</i> A pass means no length/count symbol survives. Check the units of $d_i$, then evaluate $2\pi\chi$ at $g=0$. :: A: $d_i=2\pi-\sum\theta_i^f$ is a sum of radians (dimensionless angles) <br> $\sum_i d_i=2\pi\chi$ with $\chi=2-2(0)=2$ gives $4\pi$ <br> No $V$, $E$, or length appears, so refinement leaves it fixed <br> $\boxed{\sum d_i=4\pi \text{ for all genus-0 closed meshes}}$
- Q: Construct the operator that visits all outgoing half-edges of vertex $v$, given only $\mathrm{next}$ and $\mathrm{twin}$. <br> <b>Formula:</b> target $\mathrm{rotateCW}(h)=\mathrm{next}(\mathrm{twin}(h))$ <br> <i>Hint:</i> Start from an outgoing $h$; apply $\mathrm{twin}$ to point back into $v$, then $\mathrm{next}$ to step CCW around the next face. :: A: From an outgoing half-edge $h$ (tail at $v$), $\mathrm{twin}(h)$ points back into $v$ along the same edge <br> $\mathrm{next}(\mathrm{twin}(h))$ is the next half-edge CCW around the neighboring face, again outgoing from $v$ <br> Iterating returns to $h$ after $\deg(v)$ steps <br> $\boxed{\mathrm{rotateCW}(h)=\mathrm{next}(\mathrm{twin}(h))}$
- Q: What three sets define a triangle mesh as a combinatorial object? <br> <i>Hint:</i> The triple is written $M=(V,E,F)$. :: A: Vertices $V$, edges $E$, and faces $F$ (the triple $M=(V,E,F)$).
- Q: When is the formula $g=(2-\chi)/2$ invalid, and what replaces it? <br> <b>Formula:</b> replacement $\chi=2-2g-b$ <br> <i>Hint:</i> Think about what changes when the mesh has an edge that touches only one face. :: A: It is invalid on an open mesh (with boundary); use $\chi=2-2g-b$ where $b$ is the number of boundary loops.
- Q: How do you distinguish a scan hole from a real handle, since both lower $\chi$? <br> <i>Hint:</i> Check the boundary loops returned by `boundary_loop`. :: A: A hole has a boundary loop (fill it); a handle has no boundary loop (real tunnel, keep it).
- Q: Define an edge-manifold mesh, and say why it matters before a Laplacian solve. <br> <i>Hint:</i> State the allowed face-counts per edge, then link to twin pairing. :: A: Every edge is incident to exactly 1 face (boundary) or exactly 2 faces (interior); non-manifold edges break twin pairing and make the Laplacian asymmetric.
- CLOZE: A mesh is data not a picture: its value lives in the {{c1::connectivity}} (who is adjacent to whom), separate from the geometric embedding $\mathbf{p}: V \to \mathbb{R}^3$.
- CLOZE: The twin pointer is an involution: $\mathrm{twin}(\mathrm{twin}(h))={{c1::h}}$ with $\mathrm{twin}(h)\neq h$.
- CLOZE: A boundary edge is a half-edge whose {{c1::twin has no incident face}}. <br> <i>hint:</i> equivalently, the edge is used by only one face.
- CLOZE: For a surface with boundary the generalized Euler relation is $\chi = {{c1::2 - 2g - b}}$, where $b$ counts boundary loops.
