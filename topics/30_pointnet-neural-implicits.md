# Point clouds and neural implicits: PointNet, DeepSDF, Occupancy Networks

Placement: [[27_supervised-learning-optimization]], [[20_sdf-signed-heat]] => **PointNet / neural implicits** => (learned shape from scans)

## Why it matters for the stone scan->CAD->IFC / CADFit / robotic-fabrication workflow

Scans arrive as raw, unordered point clouds, and the surfaces you want are implicit (defined by where a function crosses zero). PointNet is the canonical way to learn directly on unordered points without forcing them into a voxel grid, which is exactly the structure of stone scan data. Neural implicits (DeepSDF, Occupancy Networks) represent a shape as the zero level set of a learned function, sitting right next to the classical signed distance field from [[20_sdf-signed-heat]]. Together they let AI propose shape completions, classify points, or fill holes, while the deterministic geometry pipeline still does the exact reconstruction and CAD work.

## Intuition

A point cloud has no natural order: shuffling the points does not change the shape, so any function of it must give the same answer regardless of order. PointNet enforces this by processing each point identically and then pooling with a symmetric operation (max), which forgets order on purpose. Neural implicits flip the usual representation: instead of storing a surface, they store a function that, asked "are you inside, outside, or on the surface at point $x$?", answers for any $x$. The surface is wherever the answer flips. Analogy: a topographic function over a landscape; the shoreline (zero level set) is the shape, defined implicitly by water level.

## The math (first principles)

**Permutation invariance.** A function $f$ on a set $\{x_1,\dots,x_n\}$ is permutation-invariant if reordering the points leaves $f$ unchanged. The general symmetric form is

$$ f(\{x_i\}) = \gamma\!\left(\operatorname*{POOL}_{i=1}^{n}\, h(x_i)\right), $$

where $h$ is a shared per-point map (an MLP applied identically to every point), $\operatorname{POOL}$ is a symmetric aggregator, and $\gamma$ is a final map. **PointNet** uses a shared MLP for $h$ and **max-pooling** for $\operatorname{POOL}$:

$$ g = \max_{i} \, h(x_i), \qquad \text{output} = \gamma(g). $$

Max-pooling makes the global feature depend only on the set, not its order. Optional learned input/feature alignment (a small T-Net predicting a transform) adds rigid-pose robustness.

**Neural signed distance (DeepSDF).** Represent a shape by a network $f_\theta(x)$ trained to output the signed distance to the surface. The surface is the zero level set:

$$ S = \{\,x \in \mathbb{R}^3 : f_\theta(x) = 0\,\}. $$

Training minimizes $\sum_i |f_\theta(x_i) - s_i|$ (often with a clamp) over sampled points $x_i$ with known signed distances $s_i$. A shape code $z$ can condition $f_\theta(x; z)$ so one network represents a family of shapes.

**Occupancy networks.** Instead of distance, predict occupancy probability $o_\theta(x) \in [0,1]$ (inside vs outside); the surface is the decision boundary $o_\theta(x) = 0.5$. Training is binary cross-entropy on sampled inside/outside points.

**Extracting the mesh.** A learned implicit is evaluated on a grid and meshed with Marching Cubes to recover an explicit surface for downstream CAD ([[24_brep-boolean-occt]]).

**Relation to classical SDF.** [[20_sdf-signed-heat]] *solves* for a signed distance field from a PDE; neural implicits *learn* one from data. Same object (a scalar field whose zero set is the surface), different machinery, and they compose: learn a prior, then refine deterministically.

## Worked example

Four 2-D points $x_1=(1,0)$, $x_2=(0,1)$, $x_3=(-1,0)$, $x_4=(0,-1)$. Let the shared map be the identity $h(x)=x$ and pool with element-wise max. The global feature is $\big(\max(1,0,-1,0),\ \max(0,1,0,-1)\big) = (1,1)$. Shuffle the points to any order and the max is still $(1,1)$: the output is permutation-invariant by construction. For a neural implicit, query the unit circle's SDF $f(x)=\lVert x\rVert - 1$: at $x_1=(1,0)$, $f=0$ (on the surface); at $(2,0)$, $f=1$ (outside); at $(0,0)$, $f=-1$ (inside). The zero set is exactly the unit circle.

## Code

Idiomatic C++ with Eigen: the PointNet symmetric-function core (shared per-point map then max-pool), and evaluating a neural-implicit-style SDF on a grid for marching cubes.

```cpp
#include <Eigen/Dense>
#include <vector>
#include <functional>
using Vec = Eigen::VectorXd;

// PointNet global feature: apply a shared map h to each point, then max-pool.
Vec pointnetGlobal(const std::vector<Vec>& pts,
                   const std::function<Vec(const Vec&)>& h) {
    Vec g = h(pts[0]);                              // initialize with first feature
    for (size_t i = 1; i < pts.size(); ++i)
        g = g.cwiseMax(h(pts[i]));                  // element-wise max == symmetric pool
    return g;                                       // invariant to point order
}

// Sample a (learned or analytic) SDF f on a grid; the zero crossings are the surface.
// Feed the grid to marching cubes to get a mesh (see 24_brep-boolean-occt downstream).
std::vector<double> sampleSdfGrid(const std::function<double(const Vec&)>& f,
                                  const Vec& lo, const Vec& hi, int n) {
    std::vector<double> vals; vals.reserve(n*n*n);
    for (int k = 0; k < n; ++k)
      for (int j = 0; j < n; ++j)
        for (int i = 0; i < n; ++i) {
            Vec x(3);
            x << lo(0)+(hi(0)-lo(0))*i/(n-1),
                 lo(1)+(hi(1)-lo(1))*j/(n-1),
                 lo(2)+(hi(2)-lo(2))*k/(n-1);
            vals.push_back(f(x));                   // f(x)=0 is the surface
        }
    return vals;
}
```

## Connections

- [[27_supervised-learning-optimization]] trains the shared MLP and implicit networks by gradient descent.
- [[20_sdf-signed-heat]] is the classical, PDE-solved signed distance field; neural implicits learn the same scalar field from data.
- [[24_brep-boolean-occt]] consumes the meshed zero level set for CAD.
- [[28_transformers-attention]] and [[29_gnn-message-passing]] are alternative set/graph learners; point transformers and point GNNs extend PointNet.
- [[03_eig-svd-pca]] underlies the T-Net alignment idea (pose normalization).

## Pitfalls and failure modes

- **Forcing points into voxels.** Voxelizing a sparse cloud wastes memory and resolution; PointNet avoids it by operating on points directly.
- **Breaking permutation invariance.** Any order-dependent operation (e.g. feeding points as an ordered sequence to an order-sensitive net) is a bug for sets.
- **Max-pool loses local structure.** Vanilla PointNet captures global features but weak local geometry; use PointNet++ or local groupings for fine detail.
- **Implicit field needs good sampling.** DeepSDF/Occupancy training needs well-chosen inside/outside or near-surface samples; bad sampling blurs the surface.
- **Featureless convex shapes overfit.** A plain stone block has little distinctive geometry; learned per-object pose/shape can overfit when data is limited (prefer classical methods there).
- **Grid resolution vs cost.** Marching cubes resolution trades surface fidelity against an $O(n^3)$ sampling cost.

## Practice

1. Implement `pointnetGlobal` with identity $h$ and confirm the max-pool output is unchanged under random point shuffles.
2. Define the SDF of a sphere $f(x)=\lVert x-c\rVert - r$, sample it on a grid, and run any marching-cubes implementation to recover the mesh.
3. Train a tiny MLP to regress that sphere's SDF from sampled $(x, f(x))$ pairs; visualize its zero level set.
4. Compare an occupancy formulation (binary cross-entropy on inside/outside) to the SDF regression on the same shape.
5. Show that max-pooling discards local detail by reconstructing a global feature from two very different clouds that share extremes.

## References

- Qi et al., "PointNet: Deep Learning on Point Sets" (2017): https://arxiv.org/abs/1612.00593
- Qi et al., "PointNet++" (local features): https://arxiv.org/abs/1706.02413
- Park et al., "DeepSDF" (2019): https://arxiv.org/abs/1901.05103
- Mescheder et al., "Occupancy Networks" (2019): https://arxiv.org/abs/1812.03828

## Flashcard seeds

- Q: Derive the general permutation-invariant form for a function of an unordered set $\{x_1,\dots,x_n\}$ (step 1 of 2). <br> <b>Formula:</b> target is $f(\{x_i\}) = \operatorname{POOL}_{i=1}^{n} h(x_i)$ with $h$ a per-point map and $\operatorname{POOL}$ a symmetric (order-free) aggregator <br> <i>Hint:</i> first make every point go through the SAME map $h$ so no index is special; then ask what kind of aggregator ignores order. :: A: Process each point with the same map $h$, so per-point features carry no index identity <br> aggregate the multiset $\{h(x_i)\}$ with an operator POOL whose value is independent of order (symmetric) <br> any permutation $P$ reorders POOL's inputs but leaves its output fixed <br> $\boxed{f(\{x_i\}) = \operatorname{POOL}_{i=1}^{n} h(x_i)}$.
- Q: From the symmetric core $\operatorname{POOL}_i h(x_i)$, derive PointNet's concrete instance (step 2 of 2). <br> <b>Formula:</b> target is $f(\{x_i\}) = \gamma\!\big(\max_i h(x_i)\big)$ <br> <i>Hint:</i> pick the specific symmetric aggregator PointNet uses ($\max$), then wrap the pooled feature in a final learned readout $\gamma$. :: A: Choose POOL $=\max$, which is order-free, giving global feature $g=\max_i h(x_i)$ <br> append a final learned map $\gamma$ to read out the answer from $g$ <br> $\boxed{f(\{x_i\}) = \gamma\!\big(\max_i h(x_i)\big)}$.
- Q: Prove max-pooling is permutation-invariant for any permutation $\pi$ of the indices. <br> <b>Formula:</b> show $\max_i h(x_{\pi(i)}) = \max_i h(x_i)$ <br> <i>Hint:</i> note that $\{h(x_{\pi(1)}),\dots\}$ is the SAME multiset as $\{h(x_1),\dots\}$, then use that $\max$ of a multiset ignores listing order. :: A: $g_\pi = \max_i h(x_{\pi(i)})$ ranges over the same multiset $\{h(x_1),\dots,h(x_n)\}$ <br> the maximum of a multiset does not depend on the listing order <br> so $g_\pi = g$ for every $\pi$ <br> $\boxed{\max_i h(x_{\pi(i)}) = \max_i h(x_i)}$.
- Q: Derive the DeepSDF surface set from a network $f_\theta(x)$ that outputs signed distance ($f_\theta<0$ inside, $f_\theta>0$ outside). <br> <b>Formula:</b> target is $S=\{x\in\mathbb{R}^3 : f_\theta(x)=0\}$ <br> <i>Hint:</i> the sign flips at the boundary, so set the value to the level that is neither $+$ nor $-$. :: A: Signed distance is negative inside and positive outside, so it changes sign exactly at the boundary <br> by continuity the boundary is where the value is neither $+$ nor $-$, i.e. $f_\theta(x)=0$ <br> collect all such points <br> $\boxed{S=\{x\in\mathbb{R}^3 : f_\theta(x)=0\}}$.
- Q: Derive the DeepSDF training objective from samples $x_i$ with known signed distances $s_i$. <br> <b>Formula:</b> target is $\min_\theta \sum_i |f_\theta(x_i)-s_i|$ <br> <i>Hint:</i> write the per-sample mismatch you want near zero ($f_\theta(x_i)-s_i$), wrap it in an L1 penalty, then sum and minimize over $\theta$. :: A: We want the prediction to match the target at every sample: $f_\theta(x_i)\approx s_i$ <br> penalize the mismatch with an L1 residual $|f_\theta(x_i)-s_i|$ (robust, often clamped) <br> sum over samples and minimize over $\theta$ <br> $\boxed{\min_\theta \sum_i |f_\theta(x_i)-s_i|}$.
- Q: Derive where an Occupancy Network places the surface, given $o_\theta(x)\in[0,1]$ = probability the point is inside. <br> <b>Formula:</b> target is $S=\{x : o_\theta(x)=0.5\}$ <br> <i>Hint:</i> inside $\to 1$, outside $\to 0$; the boundary is the level halfway between the two classes. :: A: Inside has $o_\theta\to 1$, outside has $o_\theta\to 0$, so the boundary is the maximally-uncertain level <br> that is the decision threshold halfway between the two classes <br> set $o_\theta(x)=\tfrac12$ <br> $\boxed{S=\{x : o_\theta(x)=0.5\}}$.
- Q: Construct the analytic SDF of a sphere of center $c$ and radius $r$, and state its sign regime. <br> <b>Formula:</b> target is $f(x)=\lVert x-c\rVert - r$ <br> <i>Hint:</i> compute the distance from the query to the center first, then subtract the radius. :: A: Signed distance is the gap between the point's distance-to-center and the radius <br> distance to center is $\lVert x-c\rVert$; subtract $r$ <br> sign is $-$ for $\lVert x-c\rVert<r$ (inside), $+$ for $>r$ (outside), $0$ on the shell <br> $\boxed{f(x)=\lVert x-c\rVert - r}$.
- Q: Worked example (2D, unit circle). Evaluate $f$ at $(1,0)$, $(2,0)$, $(0,0)$ and classify each. <br> <b>Formula:</b> $f(x)=\lVert x\rVert-1$ <br> <i>Hint:</i> compute the norm $\lVert x\rVert$ first, then subtract 1 and read the sign. :: A: $(1,0)$: $\lVert x\rVert=1$, $f=0$ -> on surface <br> $(2,0)$: $\lVert x\rVert=2$, $f=1>0$ -> outside <br> $(0,0)$: $\lVert x\rVert=0$, $f=-1<0$ -> inside <br> $\boxed{f=(0,\,+1,\,-1)}$.
- Q: Worked example (3D, mm scale). Sphere center $c=(0,0,0)$, radius $r=50$ mm; query $x=(30,40,0)$ mm. Compute $f$ and classify. <br> <b>Formula:</b> $f(x)=\lVert x-c\rVert - r$ <br> <i>Hint:</i> compute $\lVert x-c\rVert=\sqrt{30^2+40^2+0^2}$ first, then subtract 50. :: A: $\lVert x-c\rVert=\sqrt{30^2+40^2+0^2}=\sqrt{2500}=50$ <br> $f=50-50=0$ <br> $\boxed{f=0\text{ mm} \Rightarrow \text{on the surface}}$.
- Q: Worked example (3D, m scale, far query). Block-pile bound by sphere $c=(2,2,2)$ m, $r=1$ m; query $x=(5,2,2)$ m. Compute $f$. <br> <b>Formula:</b> $f(x)=\lVert x-c\rVert - r$ <br> <i>Hint:</i> form the difference $x-c=(3,0,0)$ first, take its norm, then subtract 1. :: A: $\lVert x-c\rVert=\sqrt{3^2+0+0}=3$ m <br> $f=3-1=2$ m $>0$ <br> $\boxed{f=+2\text{ m} \Rightarrow \text{outside by 2 m}}$.
- Q: Worked example (PointNet max-pool, 2D). Points $x_1=(1,0)$, $x_2=(0,1)$, $x_3=(-1,0)$, $x_4=(0,-1)$, identity $h$. Compute the global feature $g$. <br> <b>Formula:</b> $g=\max_i h(x_i)$, element-wise across coordinates <br> <i>Hint:</i> take the max of all first coordinates, then the max of all second coordinates, separately. :: A: first coord: $\max(1,0,-1,0)=1$ <br> second coord: $\max(0,1,0,-1)=1$ <br> $\boxed{g=(1,1)}$.
- Q: Closing test (invariance). Reshuffle the four points above to order $x_3,x_1,x_4,x_2$ and recompute $g$; check it matches $(1,1)$. <br> <b>Formula:</b> $g_\pi=\max_i h(x_{\pi(i)})$ element-wise <br> <i>Hint:</i> recompute both per-coordinate maxes on the new order; a pass means $g_\pi=g=(1,1)$. :: A: first coord: $\max(-1,1,0,0)=1$ <br> second coord: $\max(0,0,-1,1)=1$ <br> $g=(1,1)$, identical to the unshuffled result <br> $\boxed{g_\pi=(1,1)=g \;\checkmark\text{ invariant}}$.
- Q: Closing test (SDF zero set). Take the claimed surface point $(\tfrac{\sqrt2}{2},\tfrac{\sqrt2}{2})$ and verify it lies on the unit circle. <br> <b>Formula:</b> $f(x)=\lVert x\rVert-1$; on-surface means $f=0$ <br> <i>Hint:</i> compute $\lVert x\rVert=\sqrt{\tfrac12+\tfrac12}$; a pass is $f=0$. :: A: $\lVert x\rVert=\sqrt{\tfrac12+\tfrac12}=\sqrt1=1$ <br> $f=1-1=0$ <br> matches the claimed zero level set <br> $\boxed{f=0 \;\checkmark\text{ point lies on }S}$.
- Q: Closing test (occupancy boundary). A network gives $o_\theta=0.5$ at $x^\star$. Is $x^\star$ on $S$, and which way does $o_\theta$ move just inside? <br> <b>Formula:</b> $S=\{x:o_\theta(x)=0.5\}$, with $o_\theta\to1$ inside, $o_\theta\to0$ outside <br> <i>Hint:</i> check whether $o_\theta(x^\star)$ equals the threshold $0.5$; a pass means $x^\star\in S$ and moving inward raises $o_\theta$ above $0.5$. :: A: inside means $o_\theta\to1$, so moving inward must push $o_\theta>0.5$ <br> at $x^\star$ exactly $o_\theta=0.5$, the threshold <br> $\boxed{x^\star\in S,\; o_\theta\text{ increases inward}\;\checkmark}$.
- Q: Derive the cost of sampling an implicit on an $n\times n\times n$ grid with one $f$ call per cell. <br> <b>Formula:</b> target is $\text{cost}=O(n^3)$ <br> <i>Hint:</i> multiply the cell counts along the three axes. :: A: grid has $n$ cells per axis in 3 axes <br> total cells $=n\cdot n\cdot n=n^3$, each one $f$ call <br> $\boxed{\text{cost}=O(n^3)}$ (resolution trades against compute).
- CLOZE: A neural implicit represents a surface as the {{c1::zero level set}} of a learned function $f_\theta(x)$.
- CLOZE: PointNet attains order-invariance by applying a {{c1::shared}} map $h$ to each point and then a {{c2::symmetric (max)}} pooling. <br> <i>hint:</i> first the per-point map, then the order-free aggregator.
- CLOZE: PointNet avoids forcing points into a {{c1::voxel grid}}, operating on the points directly.
- CLOZE: DeepSDF {{c1::learns}} the signed-distance field from data, whereas the Signed Heat method {{c2::solves}} a PDE for the same scalar field.
- Q: When should you use PointNet instead of voxelizing a scan? <br> <i>Hint:</i> think about sparse, unordered clouds and memory cost. :: A: Use PointNet when the cloud is sparse/unordered; voxelizing wastes memory and resolution.
- Q: What breaks if you feed set points to an order-sensitive net? <br> <i>Hint:</i> the defining property of a set is that order does not matter. :: A: You break permutation invariance; the output depends on listing order, which is a bug for sets.
- Q: What is the limitation of vanilla PointNet's max-pool, and the fix? <br> <i>Hint:</i> global vs local geometry. :: A: It captures global features but weak local geometry; use PointNet++ or local groupings.
- Q: Why can a plain stone block overfit a per-object learned implicit? <br> <i>Hint:</i> consider how much distinctive geometry a featureless convex shape carries. :: A: It is featureless/convex with little distinctive geometry, so limited data overfits (prefer classical methods).
- Q: How do you turn a learned implicit into an explicit CAD mesh? <br> <b>Formula:</b> surface is $\{x:f_\theta(x)=0\}$ <br> <i>Hint:</i> sample $f_\theta$ on a grid, then extract the zero crossings. :: A: Evaluate $f_\theta$ on a grid and run Marching Cubes on the zero crossings.
- Q: Code idiom: which Eigen op implements the symmetric max-pool? <br> <i>Hint:</i> an element-wise max accumulated over per-point features. :: A: Element-wise `cwiseMax` accumulated over the per-point features.
