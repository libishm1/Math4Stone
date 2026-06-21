# CADFit: recovering an editable CAD program (program synthesis)

Prereqs: [[08_nonlinear-least-squares]], [[09_robust-ransac]], [[24_brep-boolean-occt]] => this => unlocks [[26_ifc-bim]].

## Why it matters for the stone scan->CAD->IFC / CADFit workflow

A watertight stone mesh from the scan stage is geometry without intent. You cannot tell a fabricator "make this slot 3 mm deeper" because a mesh has no slot, only triangles. CADFit recovers an executable parametric program `P = [op_1 .. op_n]` whose operations (extrude, revolve, fillet, chamfer, union, cut) replay in a CAD kernel to reproduce the shape, so every machined feature becomes an editable parameter. In the stone pipeline you apply it selectively: the rough split face of a granite block stays a freeform SDF or mesh, but the sawn back face, the drilled dowel holes, and the template-cut profile become real CAD operations you can dimension, tolerance, and emit as IFC `IfcExtrudedAreaSolid` swept solids. CADFit needs a watertight mesh, so it runs after Signed Heat / reconstruction has closed the surface, and it feeds the IFC/BIM stage that needs parametric solids, not soup.

## Intuition

Imagine you are handed a finished wooden box and told to write the assembly instructions that built it. You do not start from atoms. You guess "this looks like a rectangle pushed up 40 mm," build that guess in a real CAD kernel, and overlay it on the target. If it covers 80 percent of the true volume you keep it; then you look at what is still missing or sticking out and add the next operation. CADFit is exactly this game played by a computer, scored by how much solid the program and the target share.

The one analogy to hold: it is sculpture by accretion and correction, not by direct measurement. You add a primitive, measure overlap, add another to fill the gap, subtract one to trim the excess, and stop when the overlap stops improving. The score is volumetric IoU: shared volume over combined volume. A perfect program scores 1.

## The math (first principles)

### Conventions

- Frame: right-handed, units in millimetres for shop geometry, metres only at site scale. Normalize the target to the box $[-1,1]^3$ for scale-free fitting, then rescale parameters back.
- A solid is a watertight B-Rep (boundary representation): faces, edges, vertices, plus tolerances. The target mesh $\mathcal{M}\subset\mathbb{R}^3$ must be watertight, or volume is undefined.
- Sketches live on explicit planes built from line, arc, and circle primitives. A closed sketch loop is a profile.

### The program and its operations

A CAD program is an ordered sequence

$$\Pi = (o_1, o_2, \dots, o_T).$$

Each operation is a tuple

$$o_t = (\tau_t,\ \theta_t,\ R_t),$$

where $\tau_t$ is the operation type, $\theta_t$ are the continuous parameters, and $R_t$ are the geometric references (which sketch, which edges, which axis). The six types are:

| Type $\tau$ | Parameters $\theta$ | References $R$ | Creates |
|---|---|---|---|
| Extrude | depth $d$ | sketch profile, direction | prism |
| Revolve | angle $\alpha$ | profile, axis | revolved solid |
| Fillet | radius $r$ | edge set | rounded edges |
| Chamfer | distance $c$ | edge set | beveled edges |
| Union | none | two solids | fused solid |
| Cut | none | two solids | subtracted solid |

Executing the program yields a solid $\mathrm{Solid}(\Pi) = \mathrm{Exec}(\Pi)$. The kernel validates each step; an operation that produces an invalid B-Rep is rejected before it ever enters the program.

### The objective

CADFit maximizes volumetric Intersection-over-Union between the executed solid and the target. For two solids $A,B$,

$$\mathrm{IoU}(A,B) = \frac{\mathrm{Vol}(A \cap B)}{\mathrm{Vol}(A \cup B)}.$$

The reconstruction objective searches over the space $\mathcal{S}$ of syntactically and geometrically valid programs:

$$\Pi^\* = \arg\max_{\Pi \in \mathcal{S}}\ \mathrm{IoU}\big(\mathrm{Solid}(\Pi),\ \mathcal{M}\big).$$

$\mathrm{IoU}\in[0,1]$, monotone in shared volume, scale invariant after the $[-1,1]^3$ normalization, and worth 1 only at exact coincidence. It is non-differentiable in the discrete program structure, which is why the search is combinatorial, not gradient descent over $\Pi$.

A second metric reports surface accuracy: the symmetric Chamfer distance between point sets $P$ (sampled from the target) and $\hat P$ (sampled from the result),

$$\mathrm{CD}(P,\hat P) = \tfrac12\!\left[\frac{1}{|P|}\sum_{x\in P}\min_{y\in\hat P}\lVert x-y\rVert_2 \;+\; \frac{1}{|\hat P|}\sum_{y\in\hat P}\min_{x\in P}\lVert y-x\rVert_2\right].$$

Lower is better; CD penalizes surface drift that volumetric IoU can hide.

### The stone-specific extended objective (curriculum extension, not in base CADFit)

Base CADFit maximizes IoU alone. For stone fabrication you want a program that is also simple, fabricable, assemblable, and heritage-respecting. Define a cost to minimize:

$$E(\Pi) = -\,\mathrm{IoU}\big(\mathrm{Solid}(\Pi),\mathcal{M}\big) + \lambda_1\,C_{\text{complexity}}(\Pi) + \lambda_2\,C_{\text{fab}}(\Pi) + \lambda_3\,C_{\text{assembly}}(\Pi) + \lambda_4\,C_{\text{heritage}}(\Pi).$$

Read each term. $C_{\text{complexity}}$ is operation count $T$ (Occam: fewer ops, fewer ways to be wrong). $C_{\text{fab}}$ penalizes radii below the tool diameter or undercuts a saw cannot reach. $C_{\text{assembly}}$ penalizes interfaces that collide when blocks are stacked. $C_{\text{heritage}}$ penalizes CAD-izing a face that should stay freeform (tooled or weathered stone). The $\lambda_i\ge 0$ trade accuracy against manufacturability. Setting all $\lambda_i=0$ recovers base CADFit. Mark this clearly in any writeup: the extended $E(\Pi)$ is a designed objective for the stone workflow, not a formula from the CADFit paper.

### The hybrid optimization (how the argmax is actually solved)

Direct search over $\mathcal{S}$ is intractable, so CADFit assembles a program greedily from validated candidates.

1. Sketch extraction. Cluster planar faces and slice with axis-aligned planes to get candidate profiles.
2. Candidate generation. For each profile fit Extrude and Revolve by sweeping the continuous parameter (depth, angle) and selecting a stable interval via one-sided Chamfer fit. Each candidate $c$ is kernel-validated.
3. Forward greedy selection. Starting from current solid $S=\mathrm{Solid}(\Pi)$, add the candidate with the largest positive IoU gain:

$$c^\* = \arg\max_{c \in \mathcal{C}\setminus\Pi}\ \big[\mathrm{IoU}(S\cup c,\ \mathcal{M}) - \mathrm{IoU}(S,\ \mathcal{M})\big].$$

4. Backward marginal pruning. Define the marginal effect of removing $c$:

$$\Delta^-(c\mid S) = \mathrm{IoU}\big(\mathrm{Solid}(\Pi\setminus\{c\}),\ \mathcal{M}\big) - \mathrm{IoU}(S,\ \mathcal{M}).$$

Remove the candidate with the largest $\Delta^-$ (least damage) and repeat until removing any candidate would lower IoU. This deletes redundant ops the greedy pass over-added.

5. Residual refinement. Compute the under- and over-reconstructed regions:

$$R^+ = \mathcal{M}\setminus S \quad(\text{missing}), \qquad R^- = S\setminus \mathcal{M}\quad(\text{excess}).$$

Reconstruct $R^+$ and add it by Union; reconstruct $R^-$ and remove it by Cut. Iterate until both residual volumes fall below threshold. Finishing fillet and chamfer passes round or bevel edges last, because they depend on edges that only exist after the solids are built.

This is residual refinement in the [[08_nonlinear-least-squares]] sense: fit the dominant structure, then fit what is left over, never the whole thing at once.

## Worked example

Target $\mathcal{M}$: an L-shaped stone bracket. In the $[-1,1]^3$ frame take a $2\times2\times1$ base slab minus a $1\times1\times1$ corner notch. True volume $= 2\cdot2\cdot1 - 1\cdot1\cdot1 = 4 - 1 = 3$ cubic units.

Step 1, first candidate. Extrude the $2\times2$ rectangle by depth 1. $S_1$ is a $2\times2\times1$ box, volume 4. Overlap with $\mathcal{M}$ is the whole L (the box contains it), so $\mathrm{Vol}(S_1\cap\mathcal{M})=3$. Union volume is $\mathrm{Vol}(S_1\cup\mathcal{M}) = 4$ (box already contains the L).

$$\mathrm{IoU}(S_1,\mathcal{M}) = \frac{3}{4} = 0.75.$$

Step 2, residual. $R^- = S_1\setminus\mathcal{M}$ is the $1\times1\times1$ corner the L does not have, volume 1. $R^+ = \mathcal{M}\setminus S_1 = \varnothing$. So add a Cut: extrude the $1\times1$ corner square by depth 1 and subtract it.

Step 3, after the cut. $S_2 = S_1 \setminus (\text{corner cube})$ has volume $4-1 = 3$ and matches $\mathcal{M}$ exactly. Intersection 3, union 3.

$$\mathrm{IoU}(S_2,\mathcal{M}) = \frac{3}{3} = 1.0.$$

The greedy gain of the Cut was $1.0 - 0.75 = 0.25 > 0$, so it is kept. Backward pruning tries removing either op: removing the Extrude drops IoU to 0, removing the Cut drops it to 0.75, so both stay. Final program, two operations:

```
Extrude(rect 2x2, depth=1)
Cut( Extrude(square 1x1 at corner, depth=1) )
```

Check the arithmetic by hand: $0.75 \to 1.0$, two ops, no residual. That is the entire loop in miniature.

## Code

Idiomatic C++ with OpenCASCADE (OCCT). This builds the worked L-bracket with the same op mapping CADFit uses under CadQuery, then estimates IoU by Monte Carlo point sampling. Symbols map to the math: `prism` = Extrude, `cutter` = the Cut operand, `BRepAlgoAPI_Cut` = the Cut op.

```cpp
// g++ lbracket.cpp -o lbracket \
//   $(pkg-config --cflags --libs opencascade) -std=c++17
#include <BRepPrimAPI_MakeBox.hxx>      // box helper for the example
#include <BRepPrimAPI_MakePrism.hxx>    // Extrude  -> MakePrism
#include <BRepPrimAPI_MakeRevol.hxx>    // Revolve  -> MakeRevol
#include <BRepFilletAPI_MakeFillet.hxx> // Fillet   -> BRepFilletAPI
#include <BRepFilletAPI_MakeChamfer.hxx>// Chamfer  -> BRepFilletAPI
#include <BRepAlgoAPI_Fuse.hxx>         // Union    -> Fuse
#include <BRepAlgoAPI_Cut.hxx>          // Cut      -> Cut
#include <BRepGProp.hxx>
#include <GProp_GProps.hxx>
#include <BRepClass3d_SolidClassifier.hxx>
#include <Bnd_Box.hxx>
#include <BRepBndLib.hxx>
#include <gp_Pnt.hxx>
#include <random>
#include <iostream>

// Volume of a solid via OCCT mass properties (Vol(.) in the math).
static double Volume(const TopoDS_Shape& s) {
    GProp_GProps props;
    BRepGProp::VolumeProperties(s, props);
    return props.Mass();
}

// Monte Carlo IoU(A,B) = Vol(A∩B) / Vol(A∪B), no boolean needed:
// sample a box covering both, classify each point in/out of A and B.
static double IoU(const TopoDS_Shape& A, const TopoDS_Shape& B, int N = 200000) {
    Bnd_Box bb; BRepBndLib::Add(A, bb); BRepBndLib::Add(B, bb);
    double xmin,ymin,zmin,xmax,ymax,zmax;
    bb.Get(xmin,ymin,zmin,xmax,ymax,zmax);
    BRepClass3d_SolidClassifier ca(A), cb(B);
    std::mt19937 rng(0);
    std::uniform_real_distribution<double> ux(xmin,xmax), uy(ymin,ymax), uz(zmin,zmax);
    long inter = 0, uni = 0;
    const double tol = 1e-7;
    for (int i = 0; i < N; ++i) {
        gp_Pnt p(ux(rng), uy(rng), uz(rng));
        ca.Perform(p, tol); bool inA = (ca.State() == TopAbs_IN);
        cb.Perform(p, tol); bool inB = (cb.State() == TopAbs_IN);
        if (inA && inB) ++inter;       // A ∩ B
        if (inA || inB) ++uni;         // A ∪ B
    }
    return uni ? double(inter) / double(uni) : 0.0;
}

int main() {
    // Target M: 2x2x1 slab minus a 1x1x1 corner cube (true volume 3).
    TopoDS_Shape slab   = BRepPrimAPI_MakeBox(gp_Pnt(-1,-1,0), 2, 2, 1).Shape();
    TopoDS_Shape notch  = BRepPrimAPI_MakeBox(gp_Pnt( 0, 0,0), 1, 1, 1).Shape();
    TopoDS_Shape M      = BRepAlgoAPI_Cut(slab, notch).Shape();

    // Step 1: first candidate S1 = the full slab (Extrude of the 2x2 rect).
    TopoDS_Shape S1 = slab;
    std::cout << "IoU(S1, M) = " << IoU(S1, M)
              << "  (expect ~0.75)\n";

    // Step 2: residual R- = S1 \ M is the corner cube -> apply Cut.
    TopoDS_Shape cutter = notch;                       // the Cut operand
    TopoDS_Shape S2 = BRepAlgoAPI_Cut(S1, cutter).Shape();
    std::cout << "IoU(S2, M) = " << IoU(S2, M)
              << "  (expect ~1.00)\n";

    std::cout << "Vol(M)=" << Volume(M)
              << "  Vol(S2)=" << Volume(S2) << "  (expect 3, 3)\n";
    return 0;
}
```

The OCCT op mapping, exactly as CADFit's CadQuery layer compiles down: Extrude -> `BRepPrimAPI_MakePrism`, Revolve -> `BRepPrimAPI_MakeRevol`, Fillet -> `BRepFilletAPI_MakeFillet`, Chamfer -> `BRepFilletAPI_MakeChamfer`, Union -> `BRepAlgoAPI_Fuse`, Cut -> `BRepAlgoAPI_Cut`.

The same program in CadQuery (the language CADFit emits), shown only because it is the output format you will read and edit:

```python
import cadquery as cq

# Extrude(rect 2x2, depth=1)  ->  Cut(Extrude(square 1x1, depth=1))
result = (
    cq.Workplane("XY")
      .rect(2, 2).extrude(1)                 # Extrude  -> MakePrism
      .faces(">Z").workplane()
      .moveTo(0.5, 0.5).rect(1, 1)
      .cutBlind(-1)                          # Cut      -> BRepAlgoAPI_Cut
)
# result.val()  is the OCCT TopoDS_Shape; IoU is measured against the mesh.
```

## Connections

- [[08_nonlinear-least-squares]]: parameter sweeps (depth, angle, radius) and residual refinement are nonlinear fitting; each candidate's continuous $\theta$ is a small least-squares fit, and $R^+/R^-$ is the residual you fit next.
- [[09_robust-ransac]]: the upstream primitive extraction labels planar/cylindrical consensus regions and rejects outlier points before any sketch is formed; clean primitives are CADFit's tokens.
- [[24_brep-boolean-occt]]: Union and Cut are kernel booleans; the whole assembly stage and IoU-by-boolean rely on robust B-Rep boolean evaluation in OCCT.
- [[26_ifc-bim]]: an editable program exports as parametric IFC solids (`IfcExtrudedAreaSolid`, `IfcRevolvedAreaSolid`), which a mesh cannot; CADFit is the bridge from geometry to BIM.
- Signed Heat / reconstruction (watertight-mesh stage): supplies the closed $\mathcal{M}$ that makes $\mathrm{Vol}$ and IoU well defined; CADFit must run after it.

## Pitfalls and failure modes

- Non-watertight input. IoU and Chamfer assume volume is defined. A single boundary gap makes $\mathrm{Vol}(\mathcal{M})$ meaningless and the classifier returns garbage. Close the mesh first.
- Non-commutative operations. Fillet-then-Cut is not Cut-then-Fillet. A fillet needs edges that the later cut may delete, so finishing ops must come last and reference edges that still exist.
- Frame and handedness. A revolve about the wrong-signed axis or a sketch on a flipped plane builds the mirror solid with high false IoU on symmetric parts. Pin the right-handed frame and the axis direction explicitly.
- Tolerance mismatch. OCCT booleans fail or leave sliver faces when two solids are coincident within kernel tolerance. Scale-relative epsilon, not a hardcoded 1e-3, or the boolean throws on small features.
- Ill-posed parameters. A near-degenerate profile (collinear points, zero-area loop) makes depth/angle unidentifiable; the sweep finds no stable interval. Filter invalid loops before fitting.
- Greedy local optima. Forward selection can pick an early op that blocks the true decomposition. Backward pruning recovers some of this, but a bad first primitive can cap final IoU. Diversify candidates.
- Over-CAD-izing stone. Forcing CAD ops onto a rough split face yields a high-op program that fits noise. Gate CADFit to sawn/slot/template faces; leave natural surfaces freeform (the $\lambda_4$ heritage term, or a hard mask).
- IoU plateaus. Monte Carlo IoU has sampling noise; tiny features below the sample density read as zero gain and get dropped. Raise $N$ near small features or switch to exact boolean volume for the final accept.
- Outlier-driven primitives. If upstream RANSAC let outliers into a plane fit, the extrude direction tilts and every dependent op inherits the bias.

## Practice

1. By hand, take a target that is a $3\times3\times2$ slab with a centered $1\times1\times2$ through-hole. Write the two-op program, compute IoU after the extrude only, then after the cut. Confirm the first IoU equals $18/18=1.0$ for the solid box vs the holed target's union, then recompute correctly: intersection $16$, union $18$, IoU $=16/18\approx0.889$, and $1.0$ after the cut.
2. Modify the C++ above to add a `Fillet` of radius 0.2 on the top edges via `BRepFilletAPI_MakeFillet`, then measure IoU against an un-filleted target. Explain why IoU drops slightly and why the fillet must run after the Cut.
3. Implement forward greedy selection over three candidate extrudes (only two are needed). Print the IoU gain of each at every step and confirm the redundant one is pruned by $\Delta^-$.
4. Replace Monte Carlo IoU with exact boolean volume: compute $\mathrm{Vol}(A\cap B)$ via `BRepAlgoAPI_Common` and $\mathrm{Vol}(A\cup B)$ via `BRepAlgoAPI_Fuse`. Compare runtime and accuracy against the sampling estimate.
5. Write the stone gate: given per-face roughness (e.g. normal-variation), mask faces above a threshold so CADFit only fits the smooth (sawn) faces. Report how many ops the gate removes.
6. Harder: take the extended objective $E(\Pi)$, set $\lambda_1>0$ only, and show that adding the complexity penalty makes the optimizer prefer the two-op L-bracket over an equivalent five-op program with the same IoU.

## References

- Ghadi Nehme, Eamon Whalen, Faez Ahmed. "CADFit: Precise Mesh-to-CAD Program Generation with Hybrid Optimization." arXiv:2605.01171, 2026. https://arxiv.org/abs/2605.01171
- CADFit source code (CC BY-NC 4.0; provisional patent filed). https://github.com/ghadinehme/CADFit
- CadQuery: a Python parametric CAD scripting framework based on OCCT. https://github.com/CadQuery/cadquery and https://cadquery.readthedocs.io
- Open CASCADE Technology reference manual: BRepPrimAPI_MakePrism, BRepPrimAPI_MakeRevol, BRepFilletAPI_MakeFillet, BRepFilletAPI_MakeChamfer, BRepAlgoAPI_Fuse, BRepAlgoAPI_Cut. https://dev.opencascade.org/doc/refman/html/
- Nicholas Sharp, Keenan Crane. "A Heat Method for Generalized Signed Distance" (Signed Heat method, watertight-surface stage). https://www.cs.cmu.edu/~kmcrane/Projects/SignedHeatMethod/
- IfcOpenShell documentation (IFC export of parametric solids). https://docs.ifcopenshell.org/
- DeepSDF: Park et al., "DeepSDF: Learning Continuous Signed Distance Functions for Shape Representation," CVPR 2019 (freeform alternative to CAD ops). https://arxiv.org/abs/1901.05103

## Flashcard seeds

- Q: Derive the volumetric IoU score (step 1 of 2): partition the union of two solids into shared, excess, and missing volumes. <br> <b>Formula:</b> $\mathrm{Vol}(A\cup B)=\mathrm{Vol}(A\cap B)+\mathrm{Vol}(A\setminus B)+\mathrm{Vol}(B\setminus A)$ <br> <i>Hint:</i> the three pieces are disjoint; rename $A\setminus B$ as excess $R^-$ and $B\setminus A$ as missing $R^+$. :: A: The union splits into three disjoint pieces: shared $A\cap B$, excess $A\setminus B$, missing $B\setminus A$ <br> rename $A\setminus B=R^-$ and $B\setminus A=R^+$ <br> $\boxed{\mathrm{Vol}(A\cup B)=\mathrm{Vol}(A\cap B)+\mathrm{Vol}(R^-)+\mathrm{Vol}(R^+)}$.
- Q: Derive the IoU score (step 2 of 2) from the partition. <br> <b>Formula:</b> $\mathrm{Vol}(A\cup B)=\mathrm{Vol}(A\cap B)+R^-+R^+$ and $\mathrm{IoU}=\dfrac{\text{shared}}{\text{combined}}$ <br> <i>Hint:</i> put shared volume over the partitioned union, then collapse the denominator back to $\mathrm{Vol}(A\cup B)$. :: A: IoU is the fraction of the combined region the solids share <br> $\mathrm{IoU}=\dfrac{\mathrm{Vol}(A\cap B)}{\mathrm{Vol}(A\cap B)+R^-+R^+}$ <br> collapse the denominator <br> $\boxed{\mathrm{IoU}(A,B)=\dfrac{\mathrm{Vol}(A\cap B)}{\mathrm{Vol}(A\cup B)}}$.
- Q: Derive the CADFit reconstruction objective from "find the program whose solid best matches the target." <br> <b>Formula:</b> score $=\mathrm{IoU}(\mathrm{Solid}(\Pi),\mathcal{M})$; target $\Pi^\*=\arg\max_{\Pi\in\mathcal{S}}\,\mathrm{IoU}(\mathrm{Solid}(\Pi),\mathcal{M})$ <br> <i>Hint:</i> wrap the IoU score in an $\arg\max$ over the space $\mathcal{S}$ of valid programs. :: A: A program's quality is $\mathrm{IoU}(\mathrm{Solid}(\Pi),\mathcal{M})$ <br> maximize it over $\mathcal{S}$, the syntactically and geometrically valid programs <br> $\boxed{\Pi^\*=\arg\max_{\Pi\in\mathcal{S}}\,\mathrm{IoU}(\mathrm{Solid}(\Pi),\mathcal{M})}$.
- Q: Derive the forward greedy gain (step 1 of 2): the value of adding candidate $c$ to current solid $S$. <br> <b>Formula:</b> $\mathrm{gain}(c\mid S)=\mathrm{IoU}(S\cup c,\mathcal{M})-\mathrm{IoU}(S,\mathcal{M})$ <br> <i>Hint:</i> adding $c$ updates the solid to $S\cup c$; gain is value-after minus value-before. :: A: Adding $c$ updates the solid to $S\cup c$, with value $\mathrm{IoU}(S\cup c,\mathcal{M})$ <br> marginal gain is value-after minus value-before <br> $\boxed{\mathrm{gain}(c\mid S)=\mathrm{IoU}(S\cup c,\mathcal{M})-\mathrm{IoU}(S,\mathcal{M})}$.
- Q: Derive which candidate the greedy step picks (step 2 of 2). <br> <b>Formula:</b> $c^\*=\arg\max_{c\in\mathcal{C}\setminus\Pi}\big[\mathrm{IoU}(S\cup c,\mathcal{M})-\mathrm{IoU}(S,\mathcal{M})\big]$ <br> <i>Hint:</i> take the $\arg\max$ of the gain over unused candidates $\mathcal{C}\setminus\Pi$, and accept only if it is positive. :: A: Greedy keeps the largest gain over unused candidates <br> pick $c$ maximizing $\mathrm{gain}(c\mid S)$, accept only if the max is $>0$ <br> $\boxed{c^\*=\arg\max_{c\in\mathcal{C}\setminus\Pi}\big[\mathrm{IoU}(S\cup c,\mathcal{M})-\mathrm{IoU}(S,\mathcal{M})\big]}$.
- Q: Derive the backward-pruning marginal $\Delta^-(c\mid S)$ for removing op $c$. <br> <b>Formula:</b> $\Delta^-(c\mid S)=\mathrm{IoU}(\mathrm{Solid}(\Pi\setminus\{c\}),\mathcal{M})-\mathrm{IoU}(S,\mathcal{M})$ <br> <i>Hint:</i> removing $c$ gives program $\Pi\setminus\{c\}$; the effect is value-after-removal minus current value. :: A: Removing $c$ gives program $\Pi\setminus\{c\}$ with solid $\mathrm{Solid}(\Pi\setminus\{c\})$ <br> effect is value-after-removal minus current value <br> $\boxed{\Delta^-(c\mid S)=\mathrm{IoU}(\mathrm{Solid}(\Pi\setminus\{c\}),\mathcal{M})-\mathrm{IoU}(S,\mathcal{M})}$ <br> remove the candidate with the largest (least negative) $\Delta^-$.
- Q: Derive the two residual sets used in refinement from "what is missing" and "what is excess." <br> <b>Formula:</b> $R^+=\mathcal{M}\setminus S$, $R^-=S\setminus\mathcal{M}$ <br> <i>Hint:</i> missing = target minus current solid (Union it in); excess = solid minus target (Cut it out). :: A: Missing = target not yet covered by $S$: $\mathcal{M}\setminus S$ <br> excess = solid not in the target: $S\setminus\mathcal{M}$ <br> $\boxed{R^+=\mathcal{M}\setminus S\ (\text{Union in}),\quad R^-=S\setminus\mathcal{M}\ (\text{Cut out})}$.
- Q: Compute $\mathrm{IoU}(S_1,\mathcal{M})$ for the first L-bracket candidate: a $2\times2\times1$ box $S_1$ over the L-shaped target (true volume 3). <br> <b>Formula:</b> $\mathrm{IoU}=\dfrac{\mathrm{Vol}(S_1\cap\mathcal{M})}{\mathrm{Vol}(S_1\cup\mathcal{M})}$ <br> <i>Hint:</i> the box contains the whole L, so $S_1\cap\mathcal{M}=\mathcal{M}$ and $S_1\cup\mathcal{M}=S_1$; compute $\mathrm{Vol}(S_1)$ first. :: A: $\mathrm{Vol}(S_1)=2\cdot2\cdot1=4$; box contains the L so $\mathrm{Vol}(S_1\cap\mathcal{M})=3$ <br> union: box already contains L, so $\mathrm{Vol}(S_1\cup\mathcal{M})=4$ <br> $\boxed{\mathrm{IoU}(S_1,\mathcal{M})=\tfrac34=0.75}$.
- Q: From $\mathrm{IoU}(S_1,\mathcal{M})=0.75$, compute the IoU after Cutting the corner cube, and the greedy gain. <br> <b>Formula:</b> $R^-=S_1\setminus\mathcal{M}$, $\mathrm{IoU}=\dfrac{\mathrm{Vol}(S_2\cap\mathcal{M})}{\mathrm{Vol}(S_2\cup\mathcal{M})}$, $\mathrm{gain}=\mathrm{IoU}(S_2)-\mathrm{IoU}(S_1)$ <br> <i>Hint:</i> the Cut removes the $1\times1\times1$ corner of volume 1, making $S_2=\mathcal{M}$ exactly. :: A: Cut removes $R^-=S_1\setminus\mathcal{M}$ (corner volume 1), so $\mathrm{Vol}(S_2)=4-1=3$ and $S_2=\mathcal{M}$ <br> intersection $=3$, union $=3$ <br> $\mathrm{IoU}(S_2,\mathcal{M})=\tfrac33=1.0$; gain $=1.0-0.75=0.25>0$ so kept <br> $\boxed{\mathrm{IoU}(S_2,\mathcal{M})=1.0}$.
- Q: Worked example (3D, different numbers): a $3\times3\times2$ slab with a centered $1\times1\times2$ through-hole. Compute $\mathrm{IoU}(S_1,\mathcal{M})$ after the extrude-only box $S_1$. <br> <b>Formula:</b> $\mathrm{IoU}=\dfrac{\mathrm{Vol}(S_1\cap\mathcal{M})}{\mathrm{Vol}(S_1\cup\mathcal{M})}$ <br> <i>Hint:</i> the solid box contains the holed target, so intersection $=\mathrm{Vol}(\mathcal{M})$ and union $=\mathrm{Vol}(S_1)$; subtract the hole to get $\mathrm{Vol}(\mathcal{M})$. :: A: $\mathrm{Vol}(S_1)=3\cdot3\cdot2=18$; hole $=1\cdot1\cdot2=2$ so $\mathrm{Vol}(\mathcal{M})=18-2=16$ <br> box contains target: intersection $=16$, union $=18$ <br> $\boxed{\mathrm{IoU}(S_1,\mathcal{M})=\tfrac{16}{18}\approx0.889}$, then $1.0$ after the Cut.
- Q: Worked example (mm shop scale): the true sawn slab is $300\times200\times40$ mm; candidate $S$ is over-swept to $300\times200\times50$ mm. Compute $\mathrm{IoU}(S,\mathcal{M})$. <br> <b>Formula:</b> $\mathrm{IoU}=\dfrac{\mathrm{Vol}(S\cap\mathcal{M})}{\mathrm{Vol}(S\cup\mathcal{M})}$ <br> <i>Hint:</i> $S$ contains $\mathcal{M}$ (same footprint, deeper), so intersection $=\mathrm{Vol}(\mathcal{M})$ and union $=\mathrm{Vol}(S)$; compute both volumes first. :: A: $\mathrm{Vol}(\mathcal{M})=300\cdot200\cdot40=2.4\times10^6$; $\mathrm{Vol}(S)=300\cdot200\cdot50=3.0\times10^6\ \mathrm{mm^3}$ <br> intersection $=2.4\times10^6$, union $=3.0\times10^6$ <br> $\boxed{\mathrm{IoU}=2.4/3.0=0.8}$, signalling depth is 10 mm too deep.
- Q: Worked example (greedy gain): current $\mathrm{IoU}(S,\mathcal{M})=0.62$; adding a dowel-hole Cut raises it to $0.74$. What is the gain, and is it accepted? <br> <b>Formula:</b> $\mathrm{gain}(c\mid S)=\mathrm{IoU}(S\cup c,\mathcal{M})-\mathrm{IoU}(S,\mathcal{M})$, accept iff gain $>0$ <br> <i>Hint:</i> subtract the before-IoU from the after-IoU, then check the sign. :: A: $\mathrm{gain}=0.74-0.62=0.12$ <br> $0.12>0$ <br> $\boxed{\text{accept the Cut (positive gain)}}$.
- Q: Closing test (close the loop): plug the final L-bracket program back and confirm Vol and IoU. <br> <b>Formula:</b> executed program Extrude($2\times2$, d=1) then Cut($1\times1$ corner, d=1); check $\mathrm{Vol}(\mathrm{Solid}(\Pi))=\mathrm{Vol}(\mathcal{M})$ and $\mathrm{IoU}=\dfrac{\mathrm{Vol}(\cap)}{\mathrm{Vol}(\cup)}$ <br> <i>Hint:</i> a pass means executed volume equals $\mathrm{Vol}(\mathcal{M})=3$ and IoU $=1$. :: A: Executed volume $=4-1=3=\mathrm{Vol}(\mathcal{M})$ (volumes match) <br> intersection $=3$, union $=3$ <br> $\boxed{\mathrm{IoU}=3/3=1.0}$, confirming exact reconstruction.
- Q: Closing test (limiting cases): show $\mathrm{IoU}(A,A)=1$, $\mathrm{IoU}(A,B)=0$ for disjoint solids, and that the range is $[0,1]$. <br> <b>Formula:</b> $\mathrm{IoU}(A,B)=\dfrac{\mathrm{Vol}(A\cap B)}{\mathrm{Vol}(A\cup B)}$ <br> <i>Hint:</i> substitute $B=A$, then $A\cap B=\varnothing$; a pass uses $\mathrm{Vol}(A\cap B)\le\mathrm{Vol}(A\cup B)$ for the bound. :: A: Identical: $A\cap A=A$, $A\cup A=A$, so $\mathrm{IoU}=\mathrm{Vol}(A)/\mathrm{Vol}(A)=1$ <br> disjoint: $A\cap B=\varnothing$ so numerator $=0$, $\mathrm{IoU}=0$ <br> $\mathrm{Vol}(A\cap B)\le\mathrm{Vol}(A\cup B)$ always <br> $\boxed{\mathrm{IoU}\in[0,1]}$, attaining 1 only at coincidence.
- Q: Write the symmetric Chamfer distance $\mathrm{CD}(P,\hat P)$ between sampled point sets $P$ (target) and $\hat P$ (result). <br> <b>Formula:</b> $\mathrm{CD}(P,\hat P)=\tfrac12\Big[\tfrac1{|P|}\sum_{x\in P}\min_{y\in\hat P}\lVert x-y\rVert_2+\tfrac1{|\hat P|}\sum_{y\in\hat P}\min_{x\in P}\lVert y-x\rVert_2\Big]$ <br> <i>Hint:</i> average each point's nearest-neighbour distance in both directions, then halve the sum. :: A: For each $x\in P$ take the nearest $y\in\hat P$; average over $P$ <br> do the same from $\hat P$ to $P$; halve the sum of the two means <br> $\boxed{\mathrm{CD}(P,\hat P)=\tfrac12\Big[\tfrac1{|P|}\!\sum_{x\in P}\min_{y\in\hat P}\lVert x-y\rVert_2+\tfrac1{|\hat P|}\!\sum_{y\in\hat P}\min_{x\in P}\lVert y-x\rVert_2\Big]}$ (lower is better; catches drift IoU hides).
- Q: Construct the stone extended objective $E(\Pi)$ from base CADFit plus four manufacturability penalties. <br> <b>Formula:</b> $E(\Pi)=-\mathrm{IoU}(\mathrm{Solid}(\Pi),\mathcal{M})+\lambda_1 C_{\text{complexity}}+\lambda_2 C_{\text{fab}}+\lambda_3 C_{\text{assembly}}+\lambda_4 C_{\text{heritage}}$ <br> <i>Hint:</i> turn "maximize IoU" into "minimize $-\mathrm{IoU}$," then add four $\lambda_i$-weighted cost terms. :: A: Maximizing IoU = minimizing $-\mathrm{IoU}$ <br> add weighted costs for complexity, fabricability, assembly, heritage <br> $\boxed{E(\Pi)=-\mathrm{IoU}(\mathrm{Solid}(\Pi),\mathcal{M})+\lambda_1 C_{\text{complexity}}+\lambda_2 C_{\text{fab}}+\lambda_3 C_{\text{assembly}}+\lambda_4 C_{\text{heritage}}}$ <br> all $\lambda_i=0$ recovers base CADFit (a designed objective, not in the paper).
- Q: How is each operation represented formally? <br> <b>Formula:</b> $o_t=(\tau_t,\theta_t,R_t)$ <br> <i>Hint:</i> name what each of the three slots holds. :: A: A tuple $o_t=(\tau_t,\theta_t,R_t)$: $\tau_t$ type, $\theta_t$ continuous params, $R_t$ geometric references <br> $\boxed{o_t=(\tau_t,\theta_t,R_t)}$.
- Q: What are the six CADFit operation types? <br> <i>Hint:</i> two add solids (sweep), two finish edges, two combine solids. :: A: Sweeps: Extrude, Revolve <br> edge finishers: Fillet, Chamfer <br> booleans: Union, Cut <br> $\boxed{\text{Extrude, Revolve, Fillet, Chamfer, Union, Cut}}$.
- Q: Pitfall: why must CADFit input be watertight? <br> <i>Hint:</i> think about what a single boundary gap does to $\mathrm{Vol}(\mathcal{M})$. :: A: Volume, hence IoU and Chamfer distance, is undefined on an open mesh <br> one boundary gap makes $\mathrm{Vol}(\mathcal{M})$ meaningless and the solid classifier returns garbage <br> $\boxed{\text{close the mesh before fitting}}$.
- Q: Pitfall: why must fillet and chamfer run last? <br> <i>Hint:</i> recall that operations are non-commutative and which entities a fillet references. :: A: Fillet/chamfer reference edges that earlier booleans create or destroy <br> operations are non-commutative: fillet-then-cut $\ne$ cut-then-fillet <br> $\boxed{\text{finish on edges that still exist, so do them last}}$.
- Q: When-to-use: apply CADFit vs leave geometry freeform on stone? <br> <i>Hint:</i> sort faces by whether you must dimension and tolerance them; the $\lambda_4$ term enforces this. :: A: CAD-ize sawn / slot / template / drilled faces you must dimension and tolerance <br> keep rough natural (tooled or weathered) surfaces freeform <br> $\boxed{\text{gate via the } \lambda_4 \text{ heritage term or a hard mask}}$.
- Q: Map Extrude and Cut to OCCT classes. <br> <i>Hint:</i> Extrude is a linear sweep (prism); Cut is a boolean subtraction. :: A: Extrude $\to$ `BRepPrimAPI_MakePrism` <br> Cut $\to$ `BRepAlgoAPI_Cut`.
- Q: Map Revolve, Fillet, Union to OCCT classes. <br> <i>Hint:</i> Revolve sweeps about an axis; Fillet rounds edges; Union fuses. :: A: Revolve $\to$ `BRepPrimAPI_MakeRevol` <br> Fillet $\to$ `BRepFilletAPI_MakeFillet` <br> Union $\to$ `BRepAlgoAPI_Fuse`.
- CLOZE: The union partitions as $\mathrm{Vol}(A\cup B)=\mathrm{Vol}(A\cap B)+{{c1::R^-}}+{{c2::R^+}}$, which is the denominator of IoU.
- CLOZE: CADFit emits executable {{c1::CadQuery}} programs over an {{c2::OpenCASCADE}} kernel.
- CLOZE: Targets are normalized to the box {{c1::$[-1,1]^3$}} so IoU fitting is {{c2::scale invariant}}.
- CLOZE: CADFit is released under {{c1::CC BY-NC 4.0}} with a {{c2::provisional patent}} filed, so commercial use needs licensing.
