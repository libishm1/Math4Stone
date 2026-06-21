# Trajectory planning and collision checking

Placement: [[15_inverse-kinematics]], [[33_optimization-foundations]] => **Trajectory planning & collision** => (collision-free tool motion)

## Why it matters for the stone scan->CAD->IFC / CADFit / robotic-fabrication workflow

Once you know the tool poses a job needs, the robot still has to *get there* without driving the tool, arm, or fixture into the stone, the table, or itself. Planning finds a continuous, collision-free, joint-limit-respecting path from the current configuration to the target, and then times it smoothly. This is the layer between [[15_inverse-kinematics]] (which gives one pose) and [[35_dynamics-control]] (which executes the motion).

## Intuition

Stop thinking in the workspace and start thinking in **configuration space**: one point is one whole pose of the robot (all joint angles). Obstacles in the world become forbidden regions in that space. Planning is then drawing a curve from start to goal that never enters a forbidden region. The hero diagram shows it: a straight line collides, so a planner weaves a path through the free space around the obstacles.

## The math (first principles)

Let $\mathcal{C} \subset \mathbb{R}^n$ be the **configuration space** ($n$ joints). It splits into the free space $\mathcal{C}_{\text{free}}$ and the obstacle region $\mathcal{C}_{\text{obs}}$. A **path** is a continuous map

$$ \tau:[0,1] \to \mathcal{C}, \qquad \tau(0)=q_{\text{start}},\ \ \tau(1)=q_{\text{goal}}. $$

It is **valid** (collision-free) when

$$ \tau(t) \in \mathcal{C}_{\text{free}} \ \ \text{for all } t \in [0,1], \quad\text{i.e.}\quad \tau \cap \mathcal{C}_{\text{obs}} = \varnothing. $$

**Collision checking** rarely has a closed form, so a path is checked discretely: sample it at a fine resolution $\Delta$ and test each $\tau(t_i)$ for overlap between robot geometry and obstacles (broad-phase bounding volumes, then narrow-phase mesh tests, link to [[17_mesh-topology]]).

**Sampling-based planners** avoid building $\mathcal{C}_{\text{obs}}$ explicitly:
- **RRT** grows a tree from $q_{\text{start}}$. Each iteration: sample $q_{\text{rand}} \in \mathcal{C}$; find the nearest tree node $q_{\text{near}}$; **steer** a step $\eta$ toward it,
$$ q_{\text{new}} = q_{\text{near}} + \eta\,\frac{q_{\text{rand}} - q_{\text{near}}}{\lVert q_{\text{rand}} - q_{\text{near}}\rVert}; $$
add the edge only if $q_{\text{new}}$ and the segment to it are collision-free. RRT is *probabilistically complete*: the probability of finding a path goes to 1 as samples increase.
- **PRM** samples many configurations, connects nearby valid pairs into a roadmap graph, then searches it (good for many queries in a static scene).

**Trajectory optimization** turns a path into a smooth, time-parameterized trajectory by minimizing effort subject to staying free and within limits:

$$ \min_{q(t)} \int_0^T \lVert \ddot q(t) \rVert^2\, dt \quad \text{s.t.}\quad q(t)\in\mathcal{C}_{\text{free}},\ \ \lvert\dot q\rvert\le \dot q_{\max},\ \ \lvert\ddot q\rvert\le \ddot q_{\max}, $$

solved with the tools of [[33_optimization-foundations]] (CHOMP, TrajOpt, STOMP are instances). Finally, **time scaling** assigns timestamps so velocity and acceleration limits hold.

## Worked example

A 2-joint arm, so $\mathcal{C} = \mathbb{R}^2$ (radians). Tree currently has $q_{\text{near}} = (0.2, 0.5)$. Sample $q_{\text{rand}} = (1.2, 1.5)$. Direction $q_{\text{rand}} - q_{\text{near}} = (1.0, 1.0)$, norm $\sqrt 2 \approx 1.414$, unit $(0.707, 0.707)$. With step $\eta = 0.3$: $q_{\text{new}} = (0.2, 0.5) + 0.3(0.707, 0.707) = (0.412, 0.712)$. If the segment from $q_{\text{near}}$ to $q_{\text{new}}$ is collision-free, add it; otherwise discard and resample.

## Code

Idiomatic C++ with Eigen: the core RRT extend step (sample / nearest / steer / collision-free check).

```cpp
#include <Eigen/Dense>
#include <vector>
#include <functional>
using Vec = Eigen::VectorXd;

struct Node { Vec q; int parent; };

// One RRT iteration. collisionFree(a,b) tests the segment a->b.
bool rrtExtend(std::vector<Node>& tree, const Vec& qRand, double eta,
               const std::function<bool(const Vec&, const Vec&)>& collisionFree) {
    // nearest node
    int near = 0; double best = 1e18;
    for (int i = 0; i < (int)tree.size(); ++i) {
        double d = (tree[i].q - qRand).squaredNorm();
        if (d < best) { best = d; near = i; }
    }
    // steer a step eta toward qRand
    Vec dir = qRand - tree[near].q;
    double n = dir.norm();
    if (n < 1e-9) return false;
    Vec qNew = tree[near].q + eta * dir / n;            // q_new = q_near + eta * unit(q_rand - q_near)

    if (!collisionFree(tree[near].q, qNew)) return false;
    tree.push_back({qNew, near});                       // add edge
    return true;
}
```

## Connections

- [[15_inverse-kinematics]] supplies the goal configuration(s); planning connects start to goal.
- [[33_optimization-foundations]] is the engine of trajectory optimization (smooth, constrained).
- [[35_dynamics-control]] executes the planned trajectory and respects torque limits.
- [[17_mesh-topology]] supplies the collision geometry (the stone, fixtures, the arm).
- [[09_robust-ransac]] shares the random-sampling idea that powers RRT/PRM.

## Pitfalls and failure modes

- **Narrow passages.** Sampling planners struggle to thread tight gaps (low sample probability); bias sampling or use connectors.
- **Collision-check resolution.** Too coarse $\Delta$ can step *through* a thin obstacle; match $\Delta$ to the smallest obstacle feature.
- **Probabilistic completeness, not optimality.** Plain RRT finds *a* path, not a short one; use RRT* or shortcut/smooth the result.
- **Joint limits and self-collision.** The robot's own body is part of $\mathcal{C}_{\text{obs}}$; do not forget self-collision.
- **Jerky paths.** Raw planner output has corners; smooth and time-scale before execution.

## Practice

1. Implement `rrtExtend` and grow a 2-D tree around a rectangular obstacle from start to goal; plot the tree and path.
2. Add segment collision checking by sampling the segment at resolution $\Delta$; show too-large $\Delta$ tunnels through a thin wall.
3. Shortcut the found path (try to connect non-adjacent nodes directly) and measure the length reduction.
4. Swap in goal-biased sampling (with probability $p$, sample the goal) and measure the speed-up.
5. Add a velocity/acceleration limit and time-parameterize a straight segment with a trapezoidal profile (link to [[35_dynamics-control]]).

## References

- LaValle, *Planning Algorithms* (free): http://lavalle.pl/planning/
- Karaman and Frazzoli, "Sampling-based Algorithms for Optimal Motion Planning" (RRT*): https://arxiv.org/abs/1105.1186
- Schulman et al., "Motion planning with sequential convex optimization" (TrajOpt).
- Choset et al., *Principles of Robot Motion*.

## Flashcard seeds

- Q: A 2-joint arm has tree node $q_{\text{near}}=(0.2,0.5)$ and samples $q_{\text{rand}}=(1.2,1.5)$ rad. Take one RRT step with $\eta=0.3$. <br> <b>Formula:</b> $q_{\text{new}} = q_{\text{near}} + \eta\,\frac{q_{\text{rand}}-q_{\text{near}}}{\lVert q_{\text{rand}}-q_{\text{near}}\rVert}$ <br> <i>Hint:</i> First compute the difference vector $q_{\text{rand}}-q_{\text{near}}$ and its norm, then scale the unit direction by $\eta$. :: A: Difference $= (1.0,\,1.0)$ <br> Norm $= \sqrt{1^2+1^2}=\sqrt2\approx1.414$ <br> Unit $= (0.707,\,0.707)$ <br> Step $= 0.3\,(0.707,0.707)=(0.212,0.212)$ <br> $q_{\text{new}} = (0.2,0.5)+(0.212,0.212)$ <br> $\boxed{q_{\text{new}}=(0.412,\,0.712)}$
- Q: A 2-joint arm at $q_{\text{near}}=(0,0)$ samples $q_{\text{rand}}=(0.6,0.8)$ rad. Take one RRT step with $\eta=0.5$. <br> <b>Formula:</b> $q_{\text{new}} = q_{\text{near}} + \eta\,\frac{q_{\text{rand}}-q_{\text{near}}}{\lVert q_{\text{rand}}-q_{\text{near}}\rVert}$ <br> <i>Hint:</i> The difference is $(0.6,0.8)$ — a 3-4-5 vector, so the norm is exactly 1.0; the unit direction is the difference itself. :: A: Difference $=(0.6,0.8)$ <br> Norm $=\sqrt{0.36+0.64}=\sqrt{1.0}=1.0$ <br> Unit $=(0.6,0.8)$ <br> Step $=0.5\,(0.6,0.8)=(0.3,0.4)$ <br> $q_{\text{new}}=(0,0)+(0.3,0.4)$ <br> $\boxed{q_{\text{new}}=(0.3,\,0.4)}$
- Q: A planner samples a straight segment of length $L=120$ mm. The thinnest obstacle feature is $5$ mm. How many interior collision-check samples are needed if you set $\Delta = 2$ mm? <br> <b>Formula:</b> $\#\text{samples} = \left\lceil \dfrac{L}{\Delta}\right\rceil$, with the rule $\Delta <$ smallest obstacle feature. <br> <i>Hint:</i> First confirm $\Delta=2$ mm is finer than the $5$ mm feature (so no tunneling), then divide $L$ by $\Delta$. :: A: Check: $\Delta=2 < 5$ mm, so resolution is safe (no tunneling) <br> $\#\text{samples} = \lceil 120/2\rceil = 60$ <br> $\boxed{60 \text{ samples}}$
- Q: A 3-joint arm moves from $q_{\text{start}}=(0,0,0)$ to $q_{\text{goal}}=(0.9,1.2,0.8)$ rad along a straight C-space line. What is the Euclidean path length in C-space? <br> <b>Formula:</b> $\text{length} = \lVert q_{\text{goal}}-q_{\text{start}}\rVert = \sqrt{\sum_i (q_{\text{goal},i}-q_{\text{start},i})^2}$ <br> <i>Hint:</i> Subtract componentwise to get $(0.9,1.2,0.8)$, then take the root of the sum of squares. :: A: Difference $=(0.9,1.2,0.8)$ <br> Sum of squares $=0.81+1.44+0.64=2.89$ <br> $\text{length}=\sqrt{2.89}=1.7$ <br> $\boxed{1.7 \text{ rad}}$
- Q: Goal-biased sampling: with probability $p=0.1$ the planner samples the goal directly. Over $N=50$ iterations, what is the expected number of goal-biased samples? <br> <b>Formula:</b> $\mathbb{E}[\text{goal samples}] = pN$ <br> <i>Hint:</i> This is just the mean of a binomial — multiply the bias probability by the iteration count. :: A: $\mathbb{E} = pN = 0.1 \times 50$ <br> $\boxed{5 \text{ goal-biased samples}}$
- Q: Derive the RRT steer step: show that moving a fixed distance $\eta$ from $q_{\text{near}}$ toward $q_{\text{rand}}$ gives the boxed formula. <br> <b>Formula:</b> start from "displacement = (step length) $\times$ (unit direction)"; target $q_{\text{new}} = q_{\text{near}} + \eta\,\frac{q_{\text{rand}}-q_{\text{near}}}{\lVert q_{\text{rand}}-q_{\text{near}}\rVert}$ <br> <i>Hint:</i> First write the unit direction from $q_{\text{near}}$ to $q_{\text{rand}}$, then multiply by step length $\eta$ and add to $q_{\text{near}}$. :: A: Direction vector $= q_{\text{rand}}-q_{\text{near}}$ <br> Unit direction $= \dfrac{q_{\text{rand}}-q_{\text{near}}}{\lVert q_{\text{rand}}-q_{\text{near}}\rVert}$ <br> Displacement $= \eta \times$ unit direction <br> Add to start: $q_{\text{new}} = q_{\text{near}} + \text{displacement}$ <br> $\boxed{q_{\text{new}} = q_{\text{near}} + \eta\,\frac{q_{\text{rand}}-q_{\text{near}}}{\lVert q_{\text{rand}}-q_{\text{near}}\rVert}}$
- Q: From the C-space picture, write the validity (collision-free) condition for a path $\tau:[0,1]\to\mathcal{C}$. <br> <b>Formula:</b> $\tau(t)\in\mathcal{C}_{\text{free}}\ \forall t\in[0,1]$, equivalently $\tau\cap\mathcal{C}_{\text{obs}}=\varnothing$ <br> <i>Hint:</i> State the endpoint conditions too — $\tau(0)$ and $\tau(1)$ are pinned. :: A: Every point on the path must lie in free space <br> Endpoints: $\tau(0)=q_{\text{start}},\ \tau(1)=q_{\text{goal}}$ <br> $\boxed{\tau(t)\in\mathcal{C}_{\text{free}}\ \forall t\in[0,1],\quad \tau\cap\mathcal{C}_{\text{obs}}=\varnothing}$
- Q: Write the trajectory-optimization objective that turns a path into a smooth timed trajectory. <br> <b>Formula:</b> $\min_{q(t)}\int_0^T \lVert\ddot q(t)\rVert^2\,dt$ subject to $q(t)\in\mathcal{C}_{\text{free}}$, $\lvert\dot q\rvert\le\dot q_{\max}$, $\lvert\ddot q\rvert\le\ddot q_{\max}$ <br> <i>Hint:</i> The integrand penalizes acceleration (effort); the constraints keep it free and within velocity/acceleration limits. :: A: Minimize integrated squared acceleration (effort) <br> Subject to staying collision-free and within velocity/acceleration limits <br> $\boxed{\min_{q(t)}\int_0^T\lVert\ddot q\rVert^2\,dt \ \text{ s.t. } q(t)\in\mathcal{C}_{\text{free}},\ \lvert\dot q\rvert\le\dot q_{\max},\ \lvert\ddot q\rvert\le\ddot q_{\max}}$
- Q: Closing test — verify a candidate RRT step lands the right distance from $q_{\text{near}}=(1,1)$: is $q_{\text{new}}=(1.3,1.4)$ a valid $\eta=0.5$ step? <br> <b>Formula:</b> a valid step satisfies $\lVert q_{\text{new}}-q_{\text{near}}\rVert = \eta$ <br> <i>Hint:</i> Compute the displacement norm; a pass means it equals $\eta=0.5$. :: A: Displacement $=(0.3,0.4)$ <br> Norm $=\sqrt{0.09+0.16}=\sqrt{0.25}=0.5$ <br> $0.5 = \eta$, so the step length checks out <br> $\boxed{\text{valid: } \lVert q_{\text{new}}-q_{\text{near}}\rVert = 0.5 = \eta}$
- Q: Closing test — collision-resolution check. A wall is $4$ mm thick and you set $\Delta=6$ mm. Will the discrete check reliably catch the wall? <br> <b>Formula:</b> safe iff $\Delta < $ (smallest obstacle feature) <br> <i>Hint:</i> Compare $\Delta$ against the $4$ mm feature; a pass requires $\Delta$ strictly smaller. :: A: Compare: $\Delta=6$ mm vs feature $=4$ mm <br> $6 \not< 4$, so a sample can step across the wall <br> $\boxed{\text{fails — tunneling possible; need } \Delta < 4 \text{ mm}}$
- Q: What is a configuration space $\mathcal{C}$? <br> <i>Hint:</i> Think one point = one whole pose. :: A: The space of all robot configurations; one point is one full pose (all joint values). <br> $\boxed{\mathcal{C}\subset\mathbb{R}^n,\ n=\#\text{joints}}$
- Q: Why plan in configuration space instead of the workspace? <br> <i>Hint:</i> Consider how a whole robot pose maps to a single object. :: A: One C-space point is a whole robot pose <br> So planning is just drawing a curve that avoids the obstacle region <br> $\boxed{\text{path = curve in } \mathcal{C}_{\text{free}} \text{ from } q_{\text{start}} \text{ to } q_{\text{goal}}}$
- Q: How is collision along a path segment usually checked (no closed form)? <br> <b>Formula:</b> sample at resolution $\Delta$ and test each $\tau(t_i)$ for robot–obstacle overlap <br> <i>Hint:</i> Broad-phase bounding volumes first, then narrow-phase mesh tests. :: A: Discretize the segment at resolution $\Delta$ <br> Test each sample for overlap (broad-phase, then narrow-phase) <br> $\boxed{\text{discrete check at spacing } \Delta}$
- Q: Difference between RRT and PRM? <br> <i>Hint:</i> One is single-query (tree), one is multi-query (roadmap). :: A: RRT grows a tree for a single query <br> PRM builds a reusable roadmap for many queries in a static scene <br> $\boxed{\text{RRT = single-query tree; PRM = multi-query roadmap}}$
- Q: What turns RRT into an asymptotically optimal planner? <br> <i>Hint:</i> Think about rewiring, or post-processing the raw path. :: A: RRT* (rewiring) gives asymptotic optimality <br> Or shortcut/smooth the plain-RRT path <br> $\boxed{\text{RRT* or shortcutting/smoothing}}$
- Q: Why is the robot's own body part of $\mathcal{C}_{\text{obs}}$? <br> <i>Hint:</i> Consider configurations where two links overlap. :: A: Self-collision must be avoided <br> So the robot's own geometry forbids some configurations <br> $\boxed{\text{self-collision} \Rightarrow \text{those } q \in \mathcal{C}_{\text{obs}}}$
- Q: Final step before execution after planning a geometric path? <br> <i>Hint:</i> Raw paths have corners and no timestamps. :: A: Smooth the corners <br> Time-parameterize (assign timestamps within velocity/acceleration limits) <br> $\boxed{\text{smooth + time-scale}}$
- CLOZE: RRT is {{c1::probabilistically complete}}: the chance of finding a path tends to 1 as samples increase, but it is not {{c2::optimal}}.
- CLOZE: Sampling planners struggle with {{c1::narrow passages}} because the chance of sampling inside them is low. <br> <i>hint:</i> bias sampling or use connectors to help.
- CLOZE: The RRT inner loop first finds the {{c1::nearest}} tree node to $q_{\text{rand}}$, then steers a step $\eta$ toward it.
- CLOZE: Set the collision-check resolution $\Delta$ smaller than the {{c1::smallest obstacle feature}} to avoid tunneling through thin obstacles.
