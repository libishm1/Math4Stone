# Graph neural networks and message passing

Placement: [[27_supervised-learning-optimization]], [[19_laplacian-poisson]] => **GNNs / message passing** => (mesh and point-graph learning)

## Why it matters for the stone scan->CAD->IFC / CADFit / robotic-fabrication workflow

Meshes, point neighbourhoods, kinematic chains, and CAD topology are all graphs, so graph neural networks are the natural learner for them. A GNN can classify each mesh face as rough or sawn, predict per-vertex curvature or segmentation, or learn priors over CAD-element adjacency. Because GNNs are built on the same graph Laplacian you already met in geometry processing ([[19_laplacian-poisson]]), they are not a new world: they are learnable message passing on top of operators you recognize, trained by the optimization of [[27_supervised-learning-optimization]].

## Intuition

Each node holds a feature vector. In each round, every node sends a message to its neighbours, gathers the incoming messages, and updates its own state from them. After $k$ rounds a node's state summarizes its $k$-hop neighbourhood. It is rumor-spreading on a graph: information diffuses one ring per layer. The cotangent Laplacian smooths a signal on a mesh by the same neighbour-averaging logic; a GNN learns *what* to average and *how* to transform it.

## The math (first principles)

**Graph.** A graph $G=(V,E)$ has nodes $V$, edges $E$, an adjacency matrix $A$ ($A_{ij}=1$ if $i\sim j$), and a degree matrix $D$ (diagonal, $D_{ii}=\sum_j A_{ij}$). Node features are rows of $H^{(0)}=X$.

**Graph convolutional network (GCN).** Kipf and Welling's layer propagates features through a normalized adjacency:

$$ H^{(l+1)} = \sigma\!\left(\tilde D^{-1/2}\,\tilde A\,\tilde D^{-1/2}\,H^{(l)}\,W^{(l)}\right), $$

where $\tilde A = A + I$ adds self-loops, $\tilde D$ is its degree matrix, $W^{(l)}$ is a learned weight matrix, and $\sigma$ is a nonlinearity. The factor $\tilde D^{-1/2}\tilde A\tilde D^{-1/2}$ is a symmetric normalized aggregation closely related to the graph Laplacian $L = D - A$ (normalized $L_\text{sym}=I - D^{-1/2}AD^{-1/2}$). This is the direct bridge to [[19_laplacian-poisson]].

**General message passing (MPNN).** GCN is one instance of the message-passing framework. For node $u$ with neighbours $N(u)$:

$$ h_u^{(l+1)} = \phi\!\left(h_u^{(l)},\ \bigoplus_{v\in N(u)} \psi\big(h_u^{(l)}, h_v^{(l)}, e_{uv}\big)\right), $$

where $\psi$ builds a message from a neighbour (and edge feature $e_{uv}$), $\bigoplus$ is a permutation-invariant aggregator (sum, mean, or max), and $\phi$ updates the node. Choosing $\psi$, $\bigoplus$, and $\phi$ gives GCN, GraphSAGE, GAT (attention-weighted messages, link to [[28_transformers-attention]]), and others.

**Permutation invariance.** Neighbour order is meaningless, so the aggregator must be permutation-invariant (sum/mean/max). This is the same symmetry requirement as point clouds in [[30_pointnet-neural-implicits]].

**Receptive field.** $k$ stacked layers give each node a $k$-hop receptive field. Too many layers cause *over-smoothing*: all node features collapse toward the same value.

## Worked example

A path graph $1 - 2 - 3$ with scalar features $h=(1,0,0)$. Add self-loops: node 1 neighbours $\{1,2\}$ (degree 2), node 2 neighbours $\{1,2,3\}$ (degree 3), node 3 neighbours $\{2,3\}$ (degree 2). Use a mean aggregator and identity weight. New node 2 value $=\frac{1}{3}(h_1+h_2+h_3)=\frac{1}{3}(1+0+0)=0.333$. New node 1 value $=\frac{1}{2}(h_1+h_2)=\frac{1}{2}(1+0)=0.5$. New node 3 value $=\frac{1}{2}(h_2+h_3)=0$. After one round the "1" on node 1 has diffused: $(0.5, 0.333, 0)$. Another round spreads it further. That is message passing smoothing a signal, exactly like Laplacian diffusion.

## Code

Idiomatic C++ with Eigen: one GCN layer using a sparse normalized adjacency. Build $\hat A = \tilde D^{-1/2}\tilde A\tilde D^{-1/2}$ once, then each layer is $\sigma(\hat A H W)$.

```cpp
#include <Eigen/Sparse>
#include <Eigen/Dense>
using SpMat = Eigen::SparseMatrix<double>;
using Mat   = Eigen::MatrixXd;

// Build normalized adjacency  Ahat = Dt^{-1/2} (A + I) Dt^{-1/2}.
SpMat normalizedAdjacency(const SpMat& A) {
    int n = A.rows();
    SpMat At = A;
    for (int i = 0; i < n; ++i) At.coeffRef(i, i) += 1.0;     // add self-loops (A + I)
    Eigen::VectorXd d = Eigen::VectorXd::Zero(n);
    for (int k = 0; k < At.outerSize(); ++k)
        for (SpMat::InnerIterator it(At, k); it; ++it) d(it.row()) += it.value();
    Eigen::VectorXd dinv = d.array().rsqrt();                 // Dt^{-1/2}
    SpMat Dm(n, n);
    for (int i = 0; i < n; ++i) Dm.insert(i, i) = dinv(i);
    return Dm * At * Dm;
}

// One GCN layer: H' = relu(Ahat * H * W).
Mat gcnLayer(const SpMat& Ahat, const Mat& H, const Mat& W) {
    Mat Z = Ahat * (H * W);                                   // aggregate then transform
    return Z.cwiseMax(0.0);                                   // ReLU nonlinearity
}
```

## Connections

- [[19_laplacian-poisson]] is the graph/cotangent Laplacian that GCN's normalized aggregation is built from.
- [[27_supervised-learning-optimization]] trains the weights $W^{(l)}$ by gradient descent.
- [[28_transformers-attention]] is message passing with attention-weighted messages (Graph Attention Networks).
- [[30_pointnet-neural-implicits]] shares the permutation-invariant aggregation requirement.
- [[17_mesh-topology]] supplies the graph (vertices, edges, one-ring neighbourhoods) the GNN runs on.

## Pitfalls and failure modes

- **Over-smoothing.** Too many layers make all node features converge; keep depth modest or use residual/jumping connections.
- **Wrong aggregator.** A non-permutation-invariant aggregator breaks on reordered neighbours; use sum/mean/max.
- **Normalization.** Forgetting self-loops or degree normalization causes exploding/vanishing activations on high-degree nodes.
- **Graph construction.** For point clouds, the kNN graph choice (k, radius) changes results; too sparse loses context, too dense over-smooths.
- **Scalability.** Full-graph training is memory-heavy on large meshes; use neighbour sampling or mini-batched subgraphs.

## Practice

1. Build the path-graph example and reproduce $(0.5, 0.333, 0)$ after one mean-aggregation round.
2. Construct a kNN graph over a small point cloud and run two GCN layers; visualize how a one-hot signal diffuses.
3. Show that the normalized aggregation $\tilde D^{-1/2}\tilde A\tilde D^{-1/2}$ equals $I - L_\text{sym}$ with self-loops folded in.
4. Stack 2, 4, 8, 16 layers on a mesh and measure feature variance to demonstrate over-smoothing.
5. Replace the mean aggregator with a max aggregator and compare segmentation behavior on a toy graph.

## References

- Kipf and Welling, "Semi-Supervised Classification with Graph Convolutional Networks" (2017): https://arxiv.org/abs/1609.02907
- Gilmer et al., "Neural Message Passing for Quantum Chemistry" (the MPNN framework): https://arxiv.org/abs/1704.01212
- Velickovic et al., "Graph Attention Networks": https://arxiv.org/abs/1710.10903
- Hamilton, *Graph Representation Learning* (free): https://www.cs.mcgill.ca/~wlh/grl_book/

## Flashcard seeds

- Q: Derive the GCN aggregation matrix $\hat A$ (step 1 of 2): start from raw adjacency $A$, add self-loops, then symmetric-normalize. <br> <b>Formula:</b> self-loops $\tilde A=A+I$; degree $\tilde D_{ii}=\sum_j\tilde A_{ij}$; target $\hat A=\tilde D^{-1/2}\tilde A\tilde D^{-1/2}$ <br> <i>Hint:</i> first add $I$ to $A$ so a node keeps its own feature, then rescale each edge end by $1/\sqrt{\deg}$. :: A: add self-loops: $\tilde A=A+I$ <br> degree matrix $\tilde D_{ii}=\sum_j\tilde A_{ij}$ <br> sandwich with $\tilde D^{-1/2}$ on both sides <br> $\boxed{\hat A=\tilde D^{-1/2}(A+I)\tilde D^{-1/2}}$.
- Q: Derive the full GCN layer (step 2 of 2) from the aggregator $\hat A$, a learned weight $W^{(l)}$, and a nonlinearity $\sigma$. <br> <b>Formula:</b> start $\hat A=\tilde D^{-1/2}\tilde A\tilde D^{-1/2}$; target $H^{(l+1)}=\sigma(\hat A\,H^{(l)}W^{(l)})$ <br> <i>Hint:</i> apply the three operations in order — aggregate ($\hat A H^{(l)}$), transform ($\cdot\,W^{(l)}$), then squash ($\sigma$). :: A: aggregate neighbours: $\hat A H^{(l)}$ <br> apply the learned linear map: $\hat A H^{(l)}W^{(l)}$ <br> pass through the nonlinearity <br> $\boxed{H^{(l+1)}=\sigma(\tilde D^{-1/2}\tilde A\tilde D^{-1/2}H^{(l)}W^{(l)})}$.
- Q: Derive the GCN–Laplacian link (step 1 of 2): from the symmetric normalized Laplacian, write the normalized adjacency it implies. <br> <b>Formula:</b> $L_\text{sym}=I-D^{-1/2}AD^{-1/2}$ <br> <i>Hint:</i> just solve algebraically for the $D^{-1/2}AD^{-1/2}$ term. :: A: move the adjacency term to one side <br> $D^{-1/2}AD^{-1/2}=I-L_\text{sym}$ <br> $\boxed{D^{-1/2}AD^{-1/2}=I-L_\text{sym}}$.
- Q: Show the GCN aggregator is a Laplacian smoother (step 2 of 2) using the self-looped Laplacian $\tilde L_\text{sym}$. <br> <b>Formula:</b> $\hat A=\tilde D^{-1/2}\tilde A\tilde D^{-1/2}$ and $\tilde L_\text{sym}=I-\hat A$ <br> <i>Hint:</i> substitute the self-looped version of $D^{-1/2}AD^{-1/2}=I-L_\text{sym}$, then read off what one layer does to $H^{(l)}$. :: A: fold self-loops into the previous identity: $\hat A=I-\tilde L_\text{sym}$ <br> one layer applies $H^{(l)}\mapsto(I-\tilde L_\text{sym})H^{(l)}$ <br> that is one explicit diffusion (smoothing) step <br> $\boxed{\hat A=I-\tilde L_\text{sym}}$.
- Q: Derive the GCN message+update as a special case of the MPNN template. <br> <b>Formula:</b> $h_u^{(l+1)}=\phi\big(h_u^{(l)},\bigoplus_{v\in N(u)}\psi(h_u^{(l)},h_v^{(l)},e_{uv})\big)$; target $h_u^{(l+1)}=\sigma\big(W\sum_{v\in N(u)\cup\{u\}}\tfrac{1}{\sqrt{\tilde d_u\tilde d_v}}h_v^{(l)}\big)$ <br> <i>Hint:</i> set $\psi=\tfrac{1}{\sqrt{\tilde d_u\tilde d_v}}h_v^{(l)}$, $\bigoplus=\sum$ (over neighbours plus self), and $\phi(\cdot)=\sigma(\cdot\,W)$. :: A: choose message $\psi=\tfrac{1}{\sqrt{\tilde d_u\tilde d_v}}h_v^{(l)}$ <br> choose aggregator $\bigoplus=\sum$ over $N(u)\cup\{u\}$ <br> choose update $\phi(\cdot)=\sigma(\cdot\,W)$ <br> $\boxed{h_u^{(l+1)}=\sigma\!\Big(W\sum_{v\in N(u)\cup\{u\}}\tfrac{1}{\sqrt{\tilde d_u\tilde d_v}}h_v^{(l)}\Big)}$.
- Q: Derive the symmetric-normalized weight on a single edge $u\!-\!v$ from the aggregator matrix. <br> <b>Formula:</b> $\hat A=\tilde D^{-1/2}\tilde A\tilde D^{-1/2}$, with $\tilde A_{uv}=1$ and $(\tilde D^{-1/2})_{ii}=\tilde d_i^{-1/2}$ <br> <i>Hint:</i> the entry $\hat A_{uv}$ is the product of the row factor, the unit edge, and the column factor. :: A: row factor $(\tilde D^{-1/2})_{uu}=\tilde d_u^{-1/2}$, column factor $(\tilde D^{-1/2})_{vv}=\tilde d_v^{-1/2}$ <br> multiply through the unit entry $\tilde A_{uv}=1$ <br> $\boxed{\hat A_{uv}=\dfrac{1}{\sqrt{\tilde d_u\,\tilde d_v}}}$.
- Q: Worked example (small, mean aggregator). Path $1\!-\!2\!-\!3$, features $h=(1,0,0)$, self-loops added, identity weight. Compute the new node-2 value. <br> <b>Formula:</b> $h_u'=\tfrac{1}{|\tilde N(u)|}\sum_{v\in\tilde N(u)}h_v$ where $\tilde N(u)$ includes $u$ itself <br> <i>Hint:</i> node 2's self-looped neighbourhood is $\{1,2,3\}$, so divide the feature sum by 3. :: A: node 2 self-looped neighbours $\{1,2,3\}$, degree 3 <br> mean $=\tfrac13(h_1+h_2+h_3)=\tfrac13(1+0+0)$ <br> $\boxed{h_2'=0.333}$.
- Q: Worked example (same path graph). Compute the new node-1 and node-3 values after one mean round. <br> <b>Formula:</b> $h_u'=\tfrac{1}{|\tilde N(u)|}\sum_{v\in\tilde N(u)}h_v$, self-loops included <br> <i>Hint:</i> nodes 1 and 3 are endpoints with self-looped neighbourhoods $\{1,2\}$ and $\{2,3\}$, so divide by 2 each. :: A: node 1 neighbours $\{1,2\}$, degree 2: $\tfrac12(1+0)=0.5$ <br> node 3 neighbours $\{2,3\}$, degree 2: $\tfrac12(0+0)=0$ <br> $\boxed{(h_1',h_3')=(0.5,\,0)}$, full state $(0.5,0.333,0)$.
- Q: Worked example (symmetric GCN norm, contrast with mean). Same path $1\!-\!2\!-\!3$, $h=(1,0,0)$, degrees $\tilde d=(2,3,2)$. Compute node 2's value. <br> <b>Formula:</b> $h_u'=\sum_{v\in\tilde N(u)}\tfrac{1}{\sqrt{\tilde d_u\tilde d_v}}\,h_v$ <br> <i>Hint:</i> only $h_1=1$ is nonzero, so compute the single weight $1/\sqrt{\tilde d_2\tilde d_1}=1/\sqrt{3\cdot2}$. :: A: node 2 sums over $\{1,2,3\}$ <br> $\hat A_{21}h_1=\tfrac{1}{\sqrt{3\cdot2}}\cdot1=0.408$; $\hat A_{22}h_2=0$; $\hat A_{23}h_3=0$ <br> $\boxed{h_2'=0.408}$ (vs $0.333$ for the plain mean).
- Q: Worked example (large 3D mesh scale). A scan mesh has $V=2{,}000{,}000$ vertices, average one-ring degree $\bar d=6$. Count directed messages for one GCN layer and for $k=4$ layers. <br> <b>Formula:</b> messages/layer $\approx V\bar d$; for $k$ layers $\approx kV\bar d$ <br> <i>Hint:</i> compute $V\bar d$ first, then multiply by $k=4$. :: A: per layer $\approx V\bar d=2\times10^6\times6=1.2\times10^7$ <br> over $k=4$ layers $\approx 4\times1.2\times10^7=4.8\times10^7$ <br> $\boxed{\sim1.2\times10^7\text{/layer},\ \sim4.8\times10^7\text{ for }k=4}$ (why full-graph training is memory-heavy).
- Q: Worked example (receptive field, mm scale). A kNN point-graph has mean spacing $\delta=5$ mm. What physical radius does a node see after $k=3$ layers (one hop $\approx$ one spacing)? <br> <b>Formula:</b> receptive radius $\approx k\,\delta$ <br> <i>Hint:</i> the receptive field is $k$ hops; multiply hops by the per-hop spacing $\delta$. :: A: receptive field $=k$ hops $=3$ <br> physical reach $\approx k\,\delta=3\times5\text{ mm}$ <br> $\boxed{\approx 15\text{ mm}}$; for an $m$-scale block you would need more hops or a coarser graph.
- Q: Closing test (verify, mean aggregator). Check the claim that one mean round on path $1\!-\!2\!-\!3$, $h=(1,0,0)$ gives $(0.5,0.333,0)$ by testing total-mass conservation. <br> <b>Formula:</b> compare input sum $\sum_i h_i$ with output sum $\sum_i h_i'$ <br> <i>Hint:</i> a pass here is recognizing the sums need NOT match — mean rows are not column-stochastic, so a drop is expected, not a bug. :: A: input sum $=1+0+0=1$ <br> output sum $=0.5+0.333+0=0.833\ne1$ <br> mean aggregation is not mass-conserving, so the drop is expected <br> $\boxed{(0.5,0.333,0)\text{ confirmed; sum }0.833}$.
- Q: Closing test (verify limiting case). Confirm $\hat A=I-\tilde L_\text{sym}$ is a low-pass (Laplacian) smoother by checking its action on the smoothest eigenmode. <br> <b>Formula:</b> $\hat A=I-\tilde L_\text{sym}$; on eigenmode with eigenvalue $\lambda$, $\hat A$ scales by $1-\lambda$; the flat mode has $\lambda=0$ <br> <i>Hint:</i> a pass = the $\lambda=0$ (constant) mode is preserved while large-$\lambda$ modes shrink. :: A: flat mode has $\lambda=0$, so $\hat A$ scales it by $1-0=1$ (preserved) <br> high-frequency modes have large $\lambda$, scaled by $1-\lambda<1$ (shrunk) <br> $\boxed{\text{constant mode preserved, oscillations damped}=\text{low-pass}}$.
- Q: Closing test (over-smoothing limit). Predict node-feature variance as layers $k\to\infty$ on a connected graph, then state the fix. <br> <b>Formula:</b> each layer scales mode $\lambda$ by $(1-\lambda)$ with $0<1-\lambda<1$ for $\lambda>0$ <br> <i>Hint:</i> a pass = realizing only the flat mode (constant across nodes) survives repeated shrinking, so per-node variance collapses. :: A: each layer shrinks every non-flat mode by $1-\lambda<1$ <br> as $k\to\infty$ only the flat mode survives, all nodes converge <br> variance $\to0$ <br> $\boxed{\operatorname{Var}\to0\ \Rightarrow\ \text{over-smoothing; use residual/jumping connections or fewer layers}}$.
- Q: Write the general message-passing (MPNN) update for node $u$ (recall). <br> <b>Formula:</b> $h_u^{(l+1)}=\phi\big(h_u^{(l)},\bigoplus_{v\in N(u)}\psi(h_u^{(l)},h_v^{(l)},e_{uv})\big)$ <br> <i>Hint:</i> read it inside-out — message $\psi$ per neighbour, aggregate with $\bigoplus$, update with $\phi$. :: A: $\boxed{h_u^{(l+1)}=\phi\big(h_u^{(l)},\ \bigoplus_{v\in N(u)}\psi(h_u^{(l)},h_v^{(l)},e_{uv})\big)}$.
- Q: Definition: what are $\psi$, $\bigoplus$, and $\phi$ in an MPNN? <br> <i>Hint:</i> map each symbol to one of message / aggregate / update. :: A: $\psi$ = message function (builds a message from a neighbour + edge feature) <br> $\bigoplus$ = permutation-invariant aggregator (sum/mean/max) <br> $\phi$ = update function (refreshes the node state).
- CLOZE: The message aggregator $\bigoplus$ must be {{c1::permutation-invariant}} (sum, mean, or max) because neighbour order is meaningless.
- CLOZE: A Graph Attention Network is message passing with {{c1::attention-weighted}} messages, linking GNNs to transformers.
- CLOZE: In a GCN, $\tilde A=A+I$ adds {{c1::self-loops}} so a node keeps its own feature during aggregation.
- CLOZE: Stacking $k$ GNN layers gives each node a {{c1::$k$-hop}} receptive field. <br> <i>hint:</i> one layer = one ring of diffusion.
- Q: When to use a GNN here? <br> <i>Hint:</i> ask whether the data is naturally a graph and whether you want per-node/per-edge outputs. :: A: When the data is a graph (mesh faces/vertices, kNN point neighbourhoods, CAD-element adjacency, kinematic chains) and you want per-node or per-edge predictions.
- Q: Pitfall: what breaks if you forget degree normalization or self-loops? <br> <i>Hint:</i> think about high-degree nodes and a node's own feature. :: A: Exploding or vanishing activations, especially on high-degree nodes, and loss of the node's own feature.
- Q: Pitfall: for a point-cloud kNN graph, how do $k$ (or radius) choices fail? <br> <i>Hint:</i> consider the two extremes, too sparse and too dense. :: A: Too sparse loses context; too dense over-smooths and blurs distinct regions.
- Q: Code idiom: what do you precompute once for all GCN layers, and why? <br> <b>Formula:</b> $\hat A=\tilde D^{-1/2}\tilde A\tilde D^{-1/2}$, then each layer is $\sigma(\hat A H W)$ <br> <i>Hint:</i> ask which factor does not change between layers. :: A: The normalized adjacency $\hat A=\tilde D^{-1/2}\tilde A\tilde D^{-1/2}$, because it is constant across layers; each layer is then just $\sigma(\hat A H W)$.
