# Anki 3Blue1Brown rollout — design, recipe, progress

Goal: every topic's flashcards rendered in a 3Blue1Brown style — dark Manim canvas, a hero SVG per topic, MathJax, and a link to the relevant 3b1b (or closest) video. Deck tree `Math4Stone::<Module>`.

## Connection
- Local Anki MCP (no auth): `http://127.0.0.1:3141/` (Streamable-HTTP; send `Accept: application/json, text/event-stream`; result arrives in the SSE `data:` line).
- Tunnel (OAuth, registered as `anki` in `.mcp.json`): `https://tunnel.ankimcp.ai/a940fefe-...` — needs session reconnect + sign-in for native tools.

## Card design
- **Basic** → custom note type **`3B1B Basic`**, fields `Front, Back, Visual, Video, Topic`.
  - Front template: topic chip + question. Back: + answer + `{{Visual}}` (hero SVG `<img>`) + `{{Video}}` (link).
- **Cloze** → **stock `Cloze`** model (this Anki build can't create a custom cloze type — `create_model is_cloze` errors on `new_cloze`). Restyled its CSS to 3b1b. Carry the hero SVG + video + topic chip inline in **`Back Extra`**.
- **CSS** (applied to both `3B1B Basic` and `Cloze`): dark `#0e1116`, blue topic chip `#58C4DD`, white question, teal answer `#9CE8D6`, yellow video pill `#FFD54A`, `.cloze` blue. (Full CSS stored in the models.)
- **Hero SVGs**: `D:\math\topics\media\m4s_<slug>.svg`, dark `#0e1116` + blue grid, yellow vectors, green î / red ĵ. Stored in Anki media via `store_media_file`, referenced `<img src="m4s_<slug>.svg" class="viz">`.
- **Tags**: `module::<name>`, `topic::<NN_slug>`.

## Build recipe (per module)
1. Author one hero SVG per topic → `media/m4s_<slug>.svg`.
2. `store_media_file {filename, path}` for each.
3. Parse each topic file's `## Flashcard seeds`: `- Q: .. :: A: ..` → Basic; `- CLOZE: ..` → Cloze.
4. Convert math delimiters: `$$..$$`→`\[..\]`, `$..$`→`\(..\)`.
5. `create_deck Math4Stone::<Module>`.
6. `add_notes` (model `3B1B Basic` for Q/A; model `Cloze` with `Back Extra` for cloze), `allow_duplicate=true`, chunk ≤100.

## Gotchas (learned)
- `find_notes` caps at **limit 50** (passing 2000 errors). To clear a deck, paginate-delete in a loop at limit 50 until count 0.
- Cloze notes with two distinct indices ({{c1}}, {{c2}}) generate 2 cards — so card count > note count is normal.
- `delete_notes` needs `confirmDeletion=true`.
- Do not sync to AnkiWeb without the user's OK.

## Video map (3b1b proper, or closest; verify IDs before embed)
- 01 vectors → youtube fNk_zzaMoSs (EoLA Ch.1) ✓
- 02 matrices → 3blue1brown.com/lessons/linear-transformations (EoLA Ch.3) ✓
- 31 calculus → youtube WUvTyaaNkzM (EoC Ch.1) ✓
- 03 eigen/svd/pca → 3blue1brown.com/lessons/eigenvalues (EoLA Ch.14) ✓
- 04 trig/rotations → 3blue1brown.com/lessons/quaternions-and-3d-rotation ✓
- 05 SE(3) → 3blue1brown.com/lessons/3d-transformations (EoLA Ch.5) ✓
- 32 probability → 3blue1brown.com/lessons/clt ✓
- 06 linear-least-squares → EoLA dot-products (lessons/dot-products) [closest]
- 07 calculus-jacobians → EoC derivative/chain-rule, or Khan "Jacobian matrix" (Grant) [closest]
- 33 optimization → youtube IHZwWFHWa-w "Gradient descent" (NN Ch.2)
- 08 nonlinear-least-squares → IHZwWFHWa-w [closest]
- 09 robust-ransac → none — omit button
- 34 differential-equations → youtube p_di4Zn4wz4 (DE Ch.1)
- 10 camera-pinhole → none / EoLA transformations [weak] — omit
- 11 calibration-pnp → none — omit
- 12 twists-screws-exp → youtube O85OWBJ2ayo "e to the matrix"
- 13 forward-kinematics → none — omit
- 14 manipulator-jacobian → Khan "Jacobian" (Grant) [closest]
- 15 inverse-kinematics → IHZwWFHWa-w gradient descent [closest]
- 35 dynamics-control → p_di4Zn4wz4 DE Ch.1 (pendulum)
- 16 registration-icp → none — omit
- 17 mesh-topology → none — omit
- 18 discrete-curvature → none — omit
- 19 laplacian-poisson → youtube rB83DpBJQsE "Divergence and curl"
- 20 sdf-signed-heat → youtube ly4S0oi3Yz8 "What is a PDE? (heat equation)"
- 21 remeshing-segmentation → none — omit
- 22 bsplines-nurbs → youtube aVwxzDHniEw Freya Holmér "Bézier curves" (3b1b-style, not 3b1b) 
- 23 nurbs-surfaces-fitting → aVwxzDHniEw [reuse] / omit
- 24 brep-boolean-occt → none — omit
- 25 cadfit → none — omit
- 26 ifc-bim → none — omit
- 27 supervised-learning → youtube aircAruvnKk (NN Ch.1)
- 28 transformers-attention → youtube eMlx5fFNoYc "Attention in transformers"
- 29 gnn-message-passing → aircAruvnKk NN Ch.1 [closest]
- 30 pointnet-neural-implicits → aircAruvnKk NN Ch.1 [closest]

## Progress
- [x] Models `3B1B Basic` (Visual now on the FRONT — diagram shows with the question) + restyled `Cloze`; 3b1b CSS.
- [x] **Foundations** — 7 hero SVGs, 127 notes / 132 cards, videos. DONE.
- [x] **Fab Loop** (36 state-estimation, 37 trajectory-planning, 38 toolpath/CAM) — 3 derivation SVGs, 43 notes / 46 cards, "construct the equation from the diagram" style. DONE. (New corpus topics 36-38 written; add to 00_INDEX learning order.)
- [x] **ALL 8 modules + Fab Loop rebuilt DERIVATION-FIRST** — front: "derive from the diagram"; back: step chain ending in `\boxed{final}` + the diagram. Plus diversified examples + closing test-loop cards. 28 new derivation SVGs authored (workflow w1jh1knq0).
- [x] TOTAL **1043 cards**: Foundations 201, Optimization 182, Robotics 142, CAD 150, Geometry 138, AI 102, Vision 54, Registration 28, Fab Loop 46.
- [x] Video map applied (3b1b proper for Foundations/diff-eq/AI; closest 3b1b for others; omitted where Grant has no match).
- Note: asserted YouTube IDs were link-checked and confirmed; lesson-page URLs are stable.
- [x] **FINAL: all 38 topics re-scaffolded FORMULA + HINT on every card front** (task -> Formula -> Hint -> solve; back = worked steps to boxed answer). Gentle, always-solvable difficulty.
- [x] **Rebuilt into NUMBERED subdecks** `Math4Stone::1 Foundations` … `9 Fab Loop` so the deck list + parent-deck queue read foundation-first (verified: first card across `Math4Stone` is from `1 Foundations`). Cards added per-topic (kept together, pedagogical order within). **TOTAL 1074 cards.**
- Cleanup: the 9 old unprefixed subdecks (Math4Stone::Foundations etc.) are now empty — delete them in Anki (right-click → Delete) or ignore.
