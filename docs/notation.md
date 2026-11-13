# Notation

## Nonlinear optimization problems

We consider nonlinearly constrained optimization problems of the form

$$
\begin{array}{ll} \displaystyle
\min_x      & f(x) \\
\mbox{s.t.} & l \le
\begin{Bmatrix}
c(x) \\
Ax \\
x
\end{Bmatrix}
\le u,
\end{array}
$$

where $x \in \mathbb{R}^n$, $f : \mathbb{R}^n \to \mathbb{R}$, $c: \mathbb{R}^n \to \mathbb{R}^{m_c}$, $A \in \mathbb{R}^{m_A \times n}$, and $l \in (\mathbb{R} \cup \{-\infty\})^{m_c + m_A + n}$ and $u \in (\mathbb{R} \cup \{+\infty\})^{m_c + m_A + n}$. $f$ and $c$ may be nonconvex, which results in a nonconvex optimization problem. This formulation allows for unbounded variables and equality constraints and explicitly separates general nonlinear, linear, and bound constraints, enabling solvers to readily exploit this structure. However, for the sake of simplicity of this presentation and without loss of generality, we consider the problem in the following form:

$$
\tag{NLP}
\begin{array}{ll} \displaystyle
\min_x      & f(x) \\
\mbox{s.t.} & c(x) = 0 \\
				& x \ge 0,
\end{array}
$$

where $x \in \mathbb{R}^n$, $f : \mathbb{R}^n \to \mathbb{R}$, and $c: \mathbb{R}^n \to \mathbb{R}^m$.

## Stationarity conditions

We define the scaled Lagrangian or Fritz John function of (NLP) at $(x, y, z, \pi)$:

$$\mathcal{L}_\pi(x, y, z) \stackrel{\text{def}}{=} \pi f(x) - y^T c(x) - z^T x = \pi f(x) - \sum_{j=1}^m y_j c_j(x) - \sum_{i=1}^n z_i x_i,$$

where $y$ and $z \ge 0$ are the Lagrange multipliers of the general constraints $c(x) = 0$ and the bound constraints $x \ge 0$, respectively, and $\pi \ge 0$ is an objective multiplier that is introduced to handle infeasibility or lack of constraint qualification (CQ) in a consistent way.

$\nabla_x \mathcal{L}_\pi(x, y, z)$ is the gradient of the scaled Lagrangian with respect to $x$:

$$\nabla_x \mathcal{L}_\pi(x, y, z) \stackrel{\text{def}}{=} \pi \nabla f(x) - \sum_{j=1}^m y_j \nabla c_j(x) - z.$$

$\nabla^2_{xx} \mathcal{L}_\pi(x, y)$ is the Hessian of the scaled Lagrangian with respect to $x$:

$$\nabla^2_{xx} \mathcal{L}_\pi(x, y) = \pi \nabla^2 f(x) - \sum_{j=1}^m y_j \nabla^2 c_j(x).$$

We are primarily concerned with first-order stationary points. The first-order optimality conditions (aka Fritz John conditions) of problem (NLP) at a stationary point $x^*$ state that there exist $(\pi^*, y^*, z^*)$ such that:

$$\text{(stationarity)} \quad \nabla_x \mathcal{L}_{\pi^*}(x^*, y^*, z^*) = 0$$

$$\text{(primal feasibility)} \quad c(x^*) = 0, \quad x^* \ge 0$$

$$\text{(dual feasibility)} \quad \pi^* \ge 0, \quad z_i^* \geq 0, \quad (\pi^*, y^*, z^*) \not= (0, 0, 0)$$

$$\text{(complementarity)} \quad x_i^* z_i^* = 0, \quad \forall i \in \{1, \ldots, n\}.$$

If $\pi^* > 0$, the optimality conditions are equivalent to the KKT conditions, which can be recovered by scaling the stationarity equation by $1/\pi^*$. If $\pi^* = 0$, they characterize Fritz John points, that is, feasible points at which a constraint qualification is violated.