# Corpus status and QA notes

Status as of this build. Honest accounting of what was produced and what still needs a verification pass.

## Coverage

- **35 / 35 topic files present** (`01`–`35`), every node of the master dependency graph covered.
- **621 flashcard seeds** total, all in the required `Q: … :: A: …` or `CLOZE: {{c1::…}}` formats (13–20 per file).
- Master index at [00_INDEX.md](00_INDEX.md) with module map, mermaid dependency graph, learning order, and the flashcard build plan.

## How each file was produced (provenance)

- **Agent-written with web research (26 files):** `01`–`25` and `27`. Each was researched against the two `D:\math` reports plus live web sources and returned with a structured summary.
- **Authored directly in the main thread (9 files):** `26_ifc-bim`, `28_transformers-attention`, `29_gnn-message-passing`, `30_pointnet-neural-implicits`, and the five extensions `31`–`35`. These were written without a fresh web round-trip because the subagent batch hit an account **session limit** (`resets ~2:10am Asia/Riyadh`). The math is standard and drawn from canonical sources (cited in each file), but these nine have **not been independently web-verified this session**.

## Recommended verification pass (after the session limit resets)

1. **Equation spot-check** on the 9 directly-authored files against their cited sources (esp. RK4 coefficients in `34`, the manipulator equation and computed-torque law in `35`, the GCN normalization in `29`, attention scaling in `28`).
2. **Compile-check the C++** snippets (Eigen) across all files; confirm they build and the worked-example outputs match (e.g. attention worked example `(0.802, 0.599)`, gradient check in `31`).
3. **Link integrity:** confirm every `[[NN_slug]]` wiki-link resolves to an existing file (all 35 slugs now exist).
4. **Automated critic:** re-run the completeness/correctness critic over all 35 once subagents are available, writing findings here.

## Known minor gaps to consider

- **Geodesics depth.** The master graph's "curvature, smoothing, geodesics" node is covered via the heat method inside `19`, `20`, `34`, but there is no standalone geodesic-distance file. Add one if you want fuller treatment (MMP exact geodesics, heat method, fast marching).
- **General "Optimization" overlap.** `33` (optimization foundations) and `08` (nonlinear least squares) overlap by design; `33` is the general node, `08` the least-squares specialization. Cross-links are in place.
- **Calculus split.** `31` (foundations) vs `07` (Jacobians/linearization for solvers) is an intentional split; verify they do not duplicate the chain-rule treatment.

## Not blocking

None of the above blocks flashcard generation. The seeds are atomic, self-contained, and correctly formatted, and can be imported now; the verification pass improves the long-form teaching text and code, not the card fronts/backs.
