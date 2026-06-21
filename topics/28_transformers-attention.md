# Transformers and attention

Placement: [[27_supervised-learning-optimization]], [[02_matrices]] => **Transformers / attention** => (correspondence, sequence and set models)

## Why it matters for the stone scan->CAD->IFC / CADFit / robotic-fabrication workflow

Attention is the mechanism that lets a model weigh relationships between elements: which scan patch matches which template region, which CAD operation should come next, which points belong to the same face. Modern correspondence, point-cloud transformers, and CAD-program generators (the kind of model CADFit-style synthesis can use to propose operation sequences) are built on attention. You do not need to train a large language model; you need to read attention as "mostly matrix multiplications plus a softmax," so it stops being a black box and becomes another differentiable function trained by the optimization of [[27_supervised-learning-optimization]].

## Intuition

Attention answers: for each item, which other items should I look at, and how much? Each item emits a *query* (what it is looking for), a *key* (what it offers), and a *value* (what it passes on). An item compares its query to every key, turns the similarities into weights that sum to one, and reads a weighted blend of the values. It is a soft, learned lookup table: a content-addressable average. Analogy: in a room, each person asks a question (query); everyone wears a name tag (key); you listen more to the people whose tags match your question and average what they say (values).

## The math (first principles)

**Scaled dot-product attention.** Pack queries, keys, and values as rows of matrices $Q\in\mathbb{R}^{n\times d_k}$, $K\in\mathbb{R}^{m\times d_k}$, $V\in\mathbb{R}^{m\times d_v}$:

$$ \mathrm{Attention}(Q,K,V) = \mathrm{softmax}\!\left(\frac{QK^\top}{\sqrt{d_k}}\right) V. $$

$QK^\top$ is an $n\times m$ matrix of similarities; dividing by $\sqrt{d_k}$ keeps the scores from growing with dimension; the row-wise softmax turns each row into weights that sum to 1; multiplying by $V$ produces the weighted blend.

**Softmax.** For a vector $z$, $\mathrm{softmax}(z)_i = \dfrac{e^{z_i}}{\sum_j e^{z_j}}$. It maps real scores to a probability distribution.

**Self-attention.** When $Q$, $K$, $V$ are all linear projections of the *same* input $X$ ($Q=XW_Q$, $K=XW_K$, $V=XW_V$), each element attends to all elements of the same set, mixing information across the set. Self-attention is **permutation-equivariant**: reorder the inputs and the outputs reorder the same way, which is why it suits unordered sets like point clouds (link to [[30_pointnet-neural-implicits]]).

**Multi-head attention.** Run $h$ attention functions in parallel on different learned projections and concatenate:

$$ \mathrm{MultiHead}(X) = [\,\mathrm{head}_1,\dots,\mathrm{head}_h\,]\,W_O, \quad \mathrm{head}_i = \mathrm{Attention}(XW_Q^i, XW_K^i, XW_V^i). $$

Different heads learn different relationships.

**Positional encoding.** Attention itself ignores order, so when order matters (sequences) you add a positional signal to the inputs. For pure sets you omit it.

**Feed-forward block.** Each transformer layer follows attention with a position-wise MLP, two affine maps with a nonlinearity between: $\mathrm{FFN}(x) = W_2\,\sigma(W_1 x + b_1) + b_2$, plus residual connections and normalization. The whole layer is matrix multiplies, a softmax, and a nonlinearity, repeated.

## Worked example

Three tokens, $d_k=2$. Let $Q=K=V$ rows be $x_1=(1,0)$, $x_2=(0,1)$, $x_3=(1,1)$. Take the query $q_1=(1,0)$. Scores before scaling: $q_1\cdot x_1=1$, $q_1\cdot x_2=0$, $q_1\cdot x_3=1$. Scale by $\sqrt{2}\approx1.414$: $(0.707, 0, 0.707)$. Softmax: $e^{0.707}\approx2.028$, $e^{0}=1$, $e^{0.707}\approx2.028$; sum $5.056$; weights $\approx(0.401, 0.198, 0.401)$. Output is $0.401\,x_1 + 0.198\,x_2 + 0.401\,x_3 = (0.401+0.401,\ 0.198+0.401) = (0.802, 0.599)$. Token 1 blended mostly with tokens 1 and 3, which its query matched.

## Code

Idiomatic C++ with Eigen: single-head scaled dot-product attention with a numerically stable row softmax.

```cpp
#include <Eigen/Dense>
using Mat = Eigen::MatrixXd;

// Row-wise softmax, stabilized by subtracting each row's max.
Mat softmaxRows(const Mat& S) {
    Mat out(S.rows(), S.cols());
    for (int i = 0; i < S.rows(); ++i) {
        Eigen::RowVectorXd r = S.row(i);
        double m = r.maxCoeff();                 // stability: avoid overflow in exp
        Eigen::RowVectorXd e = (r.array() - m).exp();
        out.row(i) = e / e.sum();                // weights sum to 1
    }
    return out;
}

// Attention(Q,K,V) = softmax(Q K^T / sqrt(dk)) V
Mat attention(const Mat& Q, const Mat& K, const Mat& V) {
    double dk = static_cast<double>(Q.cols());
    Mat scores = (Q * K.transpose()) / std::sqrt(dk);   // n x m similarities
    Mat W = softmaxRows(scores);                        // attention weights
    return W * V;                                       // weighted blend of values
}
```

## Connections

- [[27_supervised-learning-optimization]] trains the projection matrices $W_Q,W_K,W_V,W_O$ by gradient descent on a loss.
- [[02_matrices]] is the whole computation: attention is matrix products plus a softmax.
- [[30_pointnet-neural-implicits]] shares the permutation-symmetry concern; set transformers attend over unordered points.
- [[29_gnn-message-passing]] is closely related: attention is message passing where the messages are softmax-weighted.

## Pitfalls and failure modes

- **Missing the $1/\sqrt{d_k}$ scale.** Without it, large $d_k$ makes scores huge, the softmax saturates, and gradients vanish.
- **Unstable softmax.** Exponentiating raw scores overflows; subtract the row max first.
- **Forgetting permutation behavior.** Self-attention is permutation-equivariant; if order matters you must add positional encodings.
- **Quadratic cost.** Attention is $O(n^2)$ in the number of elements; large point clouds need sampling or efficient-attention variants.
- **Masking.** Causal or padding masks must set disallowed scores to $-\infty$ before the softmax, not after.

## Practice

1. Implement `attention` and reproduce the worked example's output $(0.802, 0.599)$.
2. Show numerically that scaling by $1/\sqrt{d_k}$ keeps softmax from saturating as $d_k$ grows.
3. Permute the input rows and confirm self-attention outputs permute identically (equivariance).
4. Add a causal mask (token $i$ may attend only to $j\le i$) and verify the weight matrix becomes lower-triangular.
5. Extend to two heads with different random projections and concatenate; compare outputs to single-head.

## References

- Vaswani et al., "Attention Is All You Need" (2017): https://arxiv.org/abs/1706.03762
- Jay Alammar, "The Illustrated Transformer": https://jalammar.github.io/illustrated-transformer/
- Lee et al., "Set Transformer" (attention over unordered sets): https://arxiv.org/abs/1810.00825
- Zhao et al., "Point Transformer" (attention for point clouds): https://arxiv.org/abs/2012.09164

## Flashcard seeds

- Q: Derive the raw similarity score between one query $q$ and the $m$ keys (rows of $K$), starting from "compare query to key by content match." <br> <b>Formula:</b> dot-product similarity $s_j = q\cdot k_j = \sum_t q_t k_{jt}$; target $s = qK^\top$ <br> <i>Hint:</i> write the single dot product first, then stack the $m$ of them as one row times $K^\top$. :: A: An item compares its query $q$ to each key $k_j$ by inner product: $s_j = q\cdot k_j = \sum_t q_t k_{jt}$ <br> stacking the $m$ scores, $k_j$ are the rows of $K$ so $q\cdot k_j$ is the $j$-th entry of $qK^\top$ <br> $\boxed{s_j = q\cdot k_j,\quad s = qK^\top}$
- Q: From the score vector $s=qK^\top$, derive the scaled scores and justify the divisor. <br> <b>Formula:</b> for unit-variance entries in $\mathbb{R}^{d_k}$, $\mathrm{Var}(q\cdot k)=d_k$; target $\tilde s = \dfrac{qK^\top}{\sqrt{d_k}}$ <br> <i>Hint:</i> divide by the standard deviation (the square root of the variance) to renormalize to variance $\approx 1$. :: A: Each score $q\cdot k=\sum_{t=1}^{d_k} q_t k_t$ is a sum of $d_k$ unit-variance terms, so $\mathrm{Var}(q\cdot k)=d_k$ and the typical magnitude is $\sqrt{d_k}$ <br> divide by the standard deviation $\sqrt{d_k}$ so the scaled variance is $\approx 1$, stopping the softmax from saturating <br> $\boxed{\tilde s = \dfrac{qK^\top}{\sqrt{d_k}}}$
- Q: From scaled scores $\tilde s$, derive the weights and the final per-query output $o$. <br> <b>Formula:</b> $w_j = e^{\tilde s_j}/\sum_l e^{\tilde s_l}$ and $o = \sum_j w_j v_j$; target $o = \mathrm{softmax}\!\big(\tfrac{qK^\top}{\sqrt{d_k}}\big)V$ <br> <i>Hint:</i> turn scores into weights with softmax first, then take the weighted sum of the value rows $v_j$. :: A: Softmax the scaled scores into weights that sum to 1: $w_j = e^{\tilde s_j}/\sum_l e^{\tilde s_l}$ <br> read a weighted blend of the value rows $v_j$ of $V$: $o = \sum_j w_j v_j = wV$ <br> $\boxed{o = \mathrm{softmax}\!\big(\tfrac{qK^\top}{\sqrt{d_k}}\big)V}$
- Q: Promote the single-query attention result to the full matrix form for $n$ queries. <br> <b>Formula:</b> single query gives $o = \mathrm{softmax}\!\big(\tfrac{qK^\top}{\sqrt{d_k}}\big)V$; target $\mathrm{Attention}(Q,K,V)=\mathrm{softmax}\!\big(\tfrac{QK^\top}{\sqrt{d_k}}\big)V$ <br> <i>Hint:</i> stack the $n$ query rows into $Q$ and apply the softmax row-wise. :: A: Stack $n$ queries as rows of $Q$; the scaled scores become the $n\times m$ matrix $\tfrac{QK^\top}{\sqrt{d_k}}$ <br> apply softmax to each row (each row's weights sum to 1), then right-multiply by $V$ <br> $\boxed{\mathrm{Attention}(Q,K,V)=\mathrm{softmax}\!\big(\tfrac{QK^\top}{\sqrt{d_k}}\big)V}$
- Q: Derive the softmax of a vector $z$ as a normalized exponential. <br> <b>Formula:</b> requirements: $w_i>0$, $\sum_i w_i=1$, monotone in $z_i$; target $\mathrm{softmax}(z)_i = \dfrac{e^{z_i}}{\sum_j e^{z_j}}$ <br> <i>Hint:</i> exponentiate to force positivity, then divide by the total to normalize. :: A: To make weights positive and monotone in $z_i$, exponentiate: $e^{z_i}>0$ and increasing in $z_i$ <br> normalize by the sum so the weights add to 1 <br> $\boxed{\mathrm{softmax}(z)_i = \dfrac{e^{z_i}}{\sum_j e^{z_j}}}$
- Q: Derive the numerically stable softmax using shift invariance. <br> <b>Formula:</b> $\mathrm{softmax}(z)_i = \dfrac{e^{z_i}}{\sum_j e^{z_j}}$; target $\dfrac{e^{z_i-\max_j z_j}}{\sum_l e^{z_l-\max_j z_j}}$ <br> <i>Hint:</i> multiply top and bottom by $e^{-c}$ and pick $c=\max_j z_j$ so the largest exponent becomes $0$. :: A: For any constant $c$, multiply numerator and denominator by $e^{-c}$: $\dfrac{e^{z_i}}{\sum_j e^{z_j}}=\dfrac{e^{z_i-c}}{\sum_j e^{z_j-c}}$ (the factor cancels) <br> choose $c=\max_j z_j$ so the largest exponent is $0$ and every $e^{z_i-c}\le 1$ (no overflow) <br> $\boxed{\mathrm{softmax}(z)_i=\dfrac{e^{z_i-\max_j z_j}}{\sum_l e^{z_l-\max_j z_j}}}$
- Q: From $Q=XW_Q,\ K=XW_K,\ V=XW_V$, derive the self-attention output of input $X$. <br> <b>Formula:</b> $\mathrm{Attention}(Q,K,V)=\mathrm{softmax}\!\big(\tfrac{QK^\top}{\sqrt{d_k}}\big)V$ with the three projections above <br> <i>Hint:</i> just substitute the three projection expressions into the attention formula. :: A: All three of $Q,K,V$ are linear projections of the same input $X$ <br> substitute into scaled dot-product attention <br> $\boxed{\mathrm{SelfAttn}(X)=\mathrm{softmax}\!\big(\tfrac{(XW_Q)(XW_K)^\top}{\sqrt{d_k}}\big)XW_V}$
- Q: Show that self-attention is permutation-equivariant: reorder the inputs and the outputs reorder the same way. <br> <b>Formula:</b> permutation matrix $P$ with $P^\top P=I$ sends $X\to PX$; target $\mathrm{Attn}(PX)=P\,\mathrm{Attn}(X)$ <br> <i>Hint:</i> push $P$ through the projections, then simplify the scores $(PQ)(PK)^\top$ and use that softmax acts row-wise. :: A: Permuting inputs by $P$ gives $Q\to PQ$, $K\to PK$, $V\to PV$ <br> scores become $(PQ)(PK)^\top=P(QK^\top)P^\top$; row-wise softmax with the trailing $P^\top$ relabels columns, and the result times $PV$ relabels rows by $P$ <br> $\boxed{\mathrm{Attn}(PX)=P\,\mathrm{Attn}(X)}$
- Q: Compute the attention weights (2D, $d_k=2$): $q=(1,0)$, keys $k_1=(1,0),\,k_2=(0,1),\,k_3=(1,1)$. <br> <b>Formula:</b> $w=\mathrm{softmax}\!\big(\tfrac{q\cdot k_j}{\sqrt{d_k}}\big)$ with $e^{0.707}\approx2.028,\ e^{0}=1$ <br> <i>Hint:</i> compute the three dot products first, then divide each by $\sqrt2\approx1.414$. :: A: Raw scores $q\cdot k_j=(1,0,1)$ <br> scale by $\sqrt2\approx1.414$: $(0.707,0,0.707)$ <br> $\exp$: $(2.028,1,2.028)$, sum $5.056$; divide each by the sum <br> $\boxed{w\approx(0.401,0.198,0.401)}$
- Q: Close the loop: with $w\approx(0.401,0.198,0.401)$ and values equal to the keys $(1,0),(0,1),(1,1)$, compute the output and check the weights are a valid distribution. <br> <b>Formula:</b> $o=\sum_j w_j v_j$; valid iff $\sum_j w_j = 1$ <br> <i>Hint:</i> form the weighted sum component-wise, then add the weights — a pass means they total $1.000$. :: A: $o=0.401(1,0)+0.198(0,1)+0.401(1,1)=(0.802,\,0.599)$ <br> check: $0.401+0.198+0.401=1.000$ (valid distribution) and the aligned keys $1,3$ dominate, as expected <br> $\boxed{o=(0.802,0.599)\ \checkmark}$
- Q: Worked example (small scale, $d_k=1$): scaled scores $s=(2,0)$, values $v_1=10,\ v_2=0$. Compute the attention output. <br> <b>Formula:</b> $w=\mathrm{softmax}(s/\sqrt{d_k})$, $o=\sum_j w_j v_j$, with $e^{2}\approx7.389,\ e^{0}=1$ <br> <i>Hint:</i> $\sqrt{1}=1$ so the scores are unchanged; exponentiate and normalize first. :: A: Scale by $\sqrt1=1$: scores stay $(2,0)$ <br> $\exp$: $(7.389,1)$, sum $8.389$; weights $(0.881,0.119)$ <br> $o=0.881\cdot10+0.119\cdot0=8.81$ <br> $\boxed{o\approx 8.81}$
- Q: Worked example (large-$d_k$, why scaling matters): two raw dot products are $20$ and $0$ at $d_k=64$. Compare the top weight with vs without the $1/\sqrt{d_k}$ scale. <br> <b>Formula:</b> $w_1=\dfrac{e^{s_1}}{e^{s_1}+e^{s_2}}$; scaled uses $s/\sqrt{64}$, with $\sqrt{64}=8$ and $e^{2.5}\approx12.18$ <br> <i>Hint:</i> compute $w_1$ for the raw scores $(20,0)$, then for the scaled scores $(2.5,0)$. :: A: Unscaled: $w_1=e^{20}/(e^{20}+1)\approx 1-2\times10^{-9}$ — saturated, near-zero gradient <br> scaled by $\sqrt{64}=8$: scores $(2.5,0)$, $w_1=e^{2.5}/(e^{2.5}+1)\approx0.924$ — soft, useful gradient <br> $\boxed{\text{scaling keeps the softmax in its informative range}}$
- Q: Worked example (uniform case, 3 keys): scaled scores are all equal, $\tilde s=(c,c,c)$. Find the weights and the output. <br> <b>Formula:</b> $w_j=e^{\tilde s_j}/\sum_l e^{\tilde s_l}$, $o=\sum_j w_j v_j$ <br> <i>Hint:</i> equal numerators over their sum — what does each $w_j$ reduce to? :: A: $w_j=e^{c}/(3e^{c})=1/3$ for every $j$ <br> the output is the plain mean of the values <br> $\boxed{o=\tfrac13(v_1+v_2+v_3)}$
- Q: Limiting case: as inverse temperature $\beta\to\infty$ (scores scaled by $\beta$), what does attention compute? Verify on two scores $\tilde s=(3,1)$. <br> <b>Formula:</b> $w_1=\dfrac{1}{1+e^{-\beta(\tilde s_1-\tilde s_2)}}$ <br> <i>Hint:</i> plug the score gap $3-1=2$ in; let $\beta\to\infty$ to see the limit, then set $\beta=1$ to verify it is already leaning. :: A: With gap $\Delta=3-1=2$: $w_1=1/(1+e^{-2\beta})\to 1$ as $\beta\to\infty$, so all weight goes to the max-score key and $o\to v_{\arg\max}$ (hard lookup) <br> check at $\beta=1$: $w_1=1/(1+e^{-2})=0.881$, already leaning to $v_1$; larger $\beta$ pushes it toward $1$ <br> $\boxed{\lim_{\beta\to\infty}\mathrm{Attn}=v_{\arg\max s}\ \checkmark}$
- Q: Derive multi-head attention from the single-head function. <br> <b>Formula:</b> $\mathrm{head}_i=\mathrm{Attention}(XW_Q^i,XW_K^i,XW_V^i)$; target $\mathrm{MultiHead}(X)=[\mathrm{head}_1,\dots,\mathrm{head}_h]\,W_O$ <br> <i>Hint:</i> run $h$ heads on their own projections, concatenate the outputs, then mix with one output matrix $W_O$. :: A: Run $h$ attention heads, each with its own learned projections: $\mathrm{head}_i=\mathrm{Attention}(XW_Q^i,XW_K^i,XW_V^i)$ <br> concatenate the $h$ head outputs side by side and mix with $W_O$ <br> $\boxed{\mathrm{MultiHead}(X)=[\mathrm{head}_1,\dots,\mathrm{head}_h]\,W_O}$
- Q: Derive how a causal mask makes token $i$ attend only to tokens $j\le i$. <br> <b>Formula:</b> mask $M$ with $M_{ij}=-\infty$ for $j>i$ else $0$; target $\mathrm{softmax}\!\big(\tfrac{QK^\top}{\sqrt{d_k}}+M\big)V$ <br> <i>Hint:</i> add $M$ to the scores before the softmax and use $e^{-\infty}=0$. :: A: Add the mask $M$ ($M_{ij}=-\infty$ for $j>i$, else $0$) to the scaled scores before softmax <br> $e^{-\infty}=0$ zeroes the future weights, leaving a lower-triangular weight matrix <br> $\boxed{\mathrm{softmax}\!\big(\tfrac{QK^\top}{\sqrt{d_k}}+M\big)V}$
- Q: Closing test (equivariance check): permute two input rows of $X$ and confirm self-attention outputs permute identically. <br> <b>Formula:</b> $\mathrm{Attn}(PX)=P\,\mathrm{Attn}(X)$ for a permutation matrix $P$ <br> <i>Hint:</i> a pass means the output rows swap in exactly the same order as the inputs, with values unchanged. :: A: Apply $P$ that swaps rows $1$ and $2$; recompute attention to get $\mathrm{Attn}(PX)$ <br> compare to swapping the same two rows of the original output $P\,\mathrm{Attn}(X)$: they match entry-for-entry <br> $\boxed{\mathrm{Attn}(PX)=P\,\mathrm{Attn}(X)\ \checkmark}$
- Q: Closing test (causal mask check): after adding the mask $M$, verify the $4\times4$ weight matrix is lower-triangular and each row still sums to $1$. <br> <b>Formula:</b> $W=\mathrm{softmax}\!\big(\tfrac{QK^\top}{\sqrt{d_k}}+M\big)$, with $W_{ij}=0$ for $j>i$ <br> <i>Hint:</i> a pass means every above-diagonal entry is $0$ and each row of the surviving entries sums to $1$. :: A: For each $j>i$ the masked score is $-\infty$ so $W_{ij}=e^{-\infty}/(\cdot)=0$ — the matrix is lower-triangular <br> within each row $i$ the softmax over the $j\le i$ entries still normalizes to $1$ <br> $\boxed{W\ \text{lower-triangular},\ \textstyle\sum_j W_{ij}=1\ \checkmark}$
- Q: What roles do the query, key, and value play in attention? <br> <i>Hint:</i> think of asking a question, wearing a name tag, and the answer you pass on. :: A: Query = what an item seeks; key = what an item offers; value = what it passes on <br> $\boxed{\text{query asks, key advertises, value is read out}}$
- Q: When do you add positional encodings, and when do you omit them? <br> <i>Hint:</i> attention alone is permutation-equivariant, so think about whether input order carries meaning. :: A: Add a positional signal when input order matters (sequences) so the model can tell positions apart <br> omit it for pure sets such as point clouds, where self-attention's permutation-equivariance is exactly what you want <br> $\boxed{\text{add for ordered sequences, omit for unordered sets}}$
- Q: Pitfall: what breaks if you omit the $1/\sqrt{d_k}$ scale at large $d_k$? <br> <b>Formula:</b> $\mathrm{Var}(q\cdot k)=d_k$, so scores scale like $\sqrt{d_k}$ <br> <i>Hint:</i> trace large scores through the softmax to its gradient. :: A: Scores grow like $\sqrt{d_k}$, so at large $d_k$ they become huge <br> the softmax saturates (one weight $\to1$, the rest $\to0$) and its gradient vanishes, so the projection matrices stop learning <br> $\boxed{\text{saturated softmax} \Rightarrow \text{vanishing gradients}}$
- Q: State the computational cost of attention in the number of elements $n$ and explain why. <br> <b>Formula:</b> $QK^\top$ is $n\times n$ <br> <i>Hint:</i> the cost is dominated by forming and storing the score matrix. :: A: Forming $QK^\top$ produces an $n\times n$ score matrix, so both compute and memory scale as $O(n^2)$ <br> large point clouds therefore need sampling or efficient-attention variants <br> $\boxed{O(n^2)\ \text{in compute and memory}}$
- Q: How does attention relate to message passing on a graph? <br> <i>Hint:</i> picture a fully connected graph and ask what sets each edge's weight. :: A: Treat the elements as nodes of a fully connected graph; each node sends a message (its value) along every edge <br> attention is message passing where each message is softmax-weighted by query-key similarity <br> $\boxed{\text{attention} = \text{softmax-weighted message passing}}$
- CLOZE: In self-attention, $Q$, $K$, $V$ are all {{c1::linear projections of the same input}} ($Q=XW_Q$, etc.).
- CLOZE: Dividing scores by {{c1::$\sqrt{d_k}$}} keeps their variance near 1 so the softmax does not saturate. <br> <i>hint:</i> it is the square root of the key dimension.
- CLOZE: A transformer layer follows attention with a position-wise {{c1::feed-forward (MLP)}} block, plus residuals and normalization.
- CLOZE: A causal/padding mask sets disallowed scores to {{c1::$-\infty$}} before the softmax, not after. <br> <i>hint:</i> the value whose exponential is $0$.
