# Unification framework

## Notation and stationarity conditions

### Notation

We start by defining the scaled Lagrangian or Fritz John function of (NLP) at $(x, y, z, \pi)$:

$$\mathcal{L}_\pi(x, y, z) \stackrel{\text{def}}{=} \pi f(x) - y^T c(x) - z^T x = \pi f(x) - \sum_{j=1}^m y_j c_j(x) - \sum_{i=1}^n z_i x_i,$$

where $y$ and $z \ge 0$ are the Lagrange multipliers of the general constraints $c(x) = 0$ and the bound constraints $x \ge 0$, respectively, and $\pi \ge 0$ is an objective multiplier that is introduced to handle infeasibility or lack of constraint qualification (CQ) in a consistent way.

$\nabla_x \mathcal{L}_\pi(x, y, z)$ is the gradient of the scaled Lagrangian with respect to $x$:

$$\nabla_x \mathcal{L}_\pi(x, y, z) \stackrel{\text{def}}{=} \pi \nabla f(x) - \sum_{j=1}^m y_j \nabla c_j(x) - z.$$

$\nabla^2_{xx} \mathcal{L}_\pi(x, y)$ is the Hessian of the scaled Lagrangian with respect to $x$:

$$\nabla^2_{xx} \mathcal{L}_\pi(x, y) = \pi \nabla^2 f(x) - \sum_{j=1}^m y_j \nabla^2 c_j(x).$$

### First-order stationarity conditions

We are primarily concerned with first-order stationary points. The first-order optimality conditions (aka Fritz John conditions) of problem (NLP) at a stationary point $x^*$ state that there exist $(\pi^*, y^*, z^*)$ such that:

$$\text{(stationarity)} \quad \nabla_x \mathcal{L}_{\pi^*}(x^*, y^*, z^*) = 0$$

$$\text{(primal feasibility)} \quad c(x^*) = 0, \quad x^* \ge 0$$

$$\text{(dual feasibility)} \quad \pi^* \ge 0, \quad z_i^* \geq 0, \quad (\pi^*, y^*, z^*) \not= (0, 0, 0)$$

$$\text{(complementarity)} \quad x_i^* z_i^* = 0, \quad \forall i \in \{1, \ldots, n\}.$$

If $\pi^* > 0$, the optimality conditions are equivalent to the KKT conditions, which can be recovered by scaling the stationarity equation by $1/\pi^*$. If $\pi^* = 0$, they characterize Fritz John points, that is, feasible points at which a constraint qualification is violated.

## A unifying framework for nonlinearly constrained optimization

A local quadratic approximation of (NLP) at iteration $k$ about the current primal-dual point is given by:

$$
\begin{array}{ll} \displaystyle
\min_x      & \frac{1}{2} (x - x^{(k)})^T H^{(k)} (x - x^{(k)}) + (\nabla f^{(k)})^T (x - x^{(k)}) \\
\mbox{s.t.} & c^{(k)} + (\nabla c^{(k)})^T (x - x^{(k)}) = 0 \\
				& x \ge 0,
$$

where $H(x, y)$ is any approximation of $\nabla^2_{xx} \mathcal{L}_\pi(x, y)$, and the superscript $(k)$ denotes evaluation at the current iterate. For convex equality-constrained problems, an iteration can be interpreted as taking a Newton step on the first-order optimality conditions of the problem, hence the name *Lagrange-Newton methods*.

In this section we introduce a unified view for describing Lagrange-Newton methods and argue that they can be assembled by combining **eight generic building blocks or *ingredients***. These ingredients are organized in four layers:

- the **reformulation** layer:
  - a **constraint relaxation strategy** constructs feasible subproblems by relaxing the general constraints;
  - an **inequality handling method** handles the inequality constraints;
- the **subproblem** layer:
  - a **Hessian model** determines the approximation of the Lagrangian Hessian;
  - an **inertia correction strategy** corrects the inertia of the Lagrangian Hessian or of the KKT matrix;
- the **subproblem solver** layer:
  - a **subproblem solver** approximately solves the subproblem by exploiting its characteristics, such as convexity, presence of constraints, presence of bounds, and presence of inequality constraints;
- the **globalization** layer:
  - a **globalization strategy** determines whether a trial iterate makes sufficient progress toward a solution and accepts or rejects it; and
  - a **globalization mechanism** controls the length of the direction and defines the recourse action taken when a trial iterate is rejected.

The role of each ingredient is shown in the following abstract algorithm.

**Algorithm 1 (Abstract algorithm).**

> **Input:** initial point $(x^{(0)}, y^{(0)}, z^{(0)})$  
> Set $k \gets 0$  
> **while** termination criteria at $(x^{(k)}, y^{(k)}, z^{(k)})$ not met **do**  
> &emsp;**repeat**  
> &emsp;&emsp;Solve a (sequence of) feasible subproblem(s) that approximate(s) (NLP) at $(x^{(k)}, y^{(k)}, z^{(k)})$  
> &emsp;&emsp;Assemble trial iterate $(\hat{x}^{(k+1)}, \hat{y}^{(k+1)}, \hat{z}^{(k+1)})$  
> &emsp;**until** $(\hat{x}^{(k+1)}, \hat{y}^{(k+1)}, \hat{z}^{(k+1)})$ is acceptable  
> &emsp;Update $(x^{(k+1)}, y^{(k+1)}, z^{(k+1)}) \gets (\hat{x}^{(k)}, \hat{y}^{(k)}, \hat{z}^{(k)})$  
> &emsp;$k \gets k + 1$  
> **end while**  
> **Output:** $(x^{(k)}, y^{(k)}, z^{(k)})$

The inner loop (**repeat**) generates and solves a feasible subproblem (possibly a sequence of feasible subproblems) until a trial iterate is accepted by the globalization strategy, and the outer loop (**while**) generates a sequence of acceptable iterates until termination.

The "wheel of strategies" (Figure 1) organizes a wide range of strategies into a coherent hierarchy. The outer layer represents the optimization layers, the middle layer represents the ingredients, and the inner layer represents the strategies. Strategies that perform similar tasks within an optimization method are listed under the same ingredient (e.g., a line search and a trust-region method are both globalization mechanisms).

![Unifying framework: the "wheel of strategies".](figures/wheel.png)

**Figure 1.** Unifying framework: the "wheel of strategies".