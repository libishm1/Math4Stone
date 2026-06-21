# Math4Stone

**3Blue1Brown-style derivation flashcards for the mathematics behind scan → CAD → IFC and robotic stone fabrication.**

A complete, dependency-ordered Anki deck (≈1,000 cards across 38 topics) plus the full teaching corpus it is built from. Every card hands you the **formula** and a **hint**, asks you to **derive or compute** the result, and shows the worked steps to a boxed answer next to a diagram — in a dark, Manim-inspired style.

The math is the dependency graph that underlies an irregular-stone pipeline: point-cloud → registration → geometry processing → CAD/B-Rep → IFC, plus the robotics (kinematics, dynamics, control), the optimization/numerics spine, and the AI layer on top.

## What's here

```
topics/                     the teaching corpus (source of truth)
  00_INDEX.md               module map, dependency graph, learning order
  00_ANKI_ROLLOUT.md        how the Anki deck is built (design + recipe)
  00_GAPS.md                QA notes
  NN_slug.md                38 topic files: intuition · math · C++ code ·
                            pitfalls · references · derivation flashcard seeds
  media/                    38 hero diagrams (SVG, Manim palette)
deep-research-report*.md    the two research reports that seeded the curriculum
crowdanki/                  the Anki deck in CrowdAnki JSON (git-friendly, importable)
```

## The 8 + 1 modules (dependency order)

1. **Foundations** — vectors · matrices · eigen/SVD/PCA · trig & rotations · SE(3) · calculus · probability
2. **Optimization & numerics** — least squares · Jacobians · optimization · nonlinear least squares · RANSAC · differential equations
3. **Vision** — pinhole camera · calibration & PnP
4. **Robotics** — twists/screws/exp · forward kinematics · manipulator Jacobian · inverse kinematics · dynamics & control
5. **Registration** — ICP / Kabsch-Umeyama
6. **Geometry processing** — mesh topology · discrete curvature · Laplacian/Poisson · SDF / Signed Heat · remeshing & segmentation
7. **CAD** — B-splines/NURBS · NURBS surfaces & fitting · B-Rep & booleans · CADFit program synthesis · IFC/BIM
8. **AI layer** — supervised learning · transformers · GNNs · PointNet / neural implicits
9. **Fab Loop** — state estimation & filtering · trajectory planning & collision · toolpath / CAM

## Card design

- **Front:** a topic diagram, a task, the **Formula** you need, and a **Hint**.
- **Back:** the step-by-step derivation ending in a `□ boxed` result, the diagram again, and a link to the matching 3Blue1Brown video where one exists.
- Difficulty is tuned to be **solvable when given the formula** — derivation chains, diversified examples (different scales/regimes), and closing "test-yourself" cards that verify the result.

## Use the deck (CrowdAnki import)

1. In Anki, install **CrowdAnki** (add-on code `1788670778`) and restart.
2. `File → Import…` and select the `crowdanki/Math4Stone` folder.
3. Study the **`Math4Stone`** deck top-down — the subdecks are numbered `1 Foundations … 9 Fab Loop` so the order follows the dependency graph. Math renders via MathJax in the Anki app.

## Rebuild the deck from source

The deck is generated from the `## Flashcard seeds` sections of the `topics/*.md` files (each line is one card) plus the SVGs in `topics/media/`. See [topics/00_ANKI_ROLLOUT.md](topics/00_ANKI_ROLLOUT.md) for the exact build recipe (note types, CSS, `$…$ → \(…\)` MathJax conversion, deck routing).

## Notes

- Educational / research material. Equations are sourced (see each topic's References section); the directly-authored topics carry a verification note in `00_GAPS.md`.
- 3Blue1Brown videos are linked for reference and belong to Grant Sanderson; this project is not affiliated with 3Blue1Brown.
