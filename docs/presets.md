## Uno presets

Uno implements presets, that is combinations of ingredients that correspond to existing solvers. At the moment, the available presets are `filtersqp` (after the trust-region restoration filter SQP solver filterSQP) and `ipopt` (after the line-search filter restoration infeasible interior-point solver IPOPT). We show below how the eight ingredients naturally arise in these two methods.

### Trust-region restoration filter SQP

Convergence of SQP filter methods has been proven under mild conditions in the context of trust-region methods and of line-search methods. The trust-region optimality QP subproblem about $(x^{(k)}, y^{(k)})$ is defined as:

$$
\tag{$QP^{(k)}(\Delta^{(l)})$}
\begin{array}{ll} \displaystyle
\min_{d_x}  & \frac{1}{2} d_x^T H_1^{(k)} d_x + (\nabla f^{(k)})^T d_x \\
\mbox{s.t.} & c^{(k)} + (\nabla c^{(k)})^T d_x = 0 \\
				& x^{(k)} + d_x \ge 0 \\
				& \|d_x\|_\infty \le \Delta^{(l)},
\end{array}
$$

where $\Delta^{(l)} > 0$ is the current trust-region radius. If the QP subproblem is infeasible ($\Delta^{(l)}$ is too small or the linearized constraints are inconsistent), we switch to feasibility restoration and solve a smooth reformulation of the $\ell_1$ feasibility problem with elastic variables $u^+ \in \mathbb{R}^m$ and $u^- \in \mathbb{R}^m$:

$$
\tag{$FQP^{(k)}(\Delta^{(l)})$}
\begin{array}{ll} \displaystyle
\min_{d_x, u^+, u^-} & \frac{1}{2} d_x^T H_0^{(k)} d_x + e^T u^+ + e^T u^- \\
\mbox{s.t.} 			& c^{(k)} + (\nabla c^{(k)})^T d_x - u^+ + u^- = 0 \\
							& x^{(k)} + d_x \ge 0 \\
							& \|d_x\|_\infty \le \Delta^{(l)} \\
							& u^+ \ge 0, \; u^- \ge 0.
\end{array}
$$

If the trial iterate $x^{(k)} + d_x$ makes sufficient progress with respect to the filter method, it is accepted. If the trust region was active at the solution of the QP ($\|d_x^*\|_\infty = \Delta^{(l)}$), we enlarge the radius. If the trial iterate is rejected, we resolve the trust-region subproblem with a smaller trust-region radius.

### Line-search filter restoration infeasible interior-point method

An infeasible interior-point method does not require feasibility with respect to the general constraints. A prerequisite is to turn inequality constraints into equality constraints using slack variables:

$$\min_{x,s} \quad f(x) \quad \text{s.t.} \quad c(x) - s = 0, \quad l_x \le x \le u_x, \quad l_c \le s \le u_c.$$

Provided that the subproblem is convex, the primal-dual direction is the solution of the primal-dual system:

$$
\tag{$IPSP_\mu$}
\begin{pmatrix}
H_1^{(k)} + (X^{(k)})^{-1} Z^{(k)} + \delta_w I & \nabla c^{(k)} \\
(\nabla c^{(k)})^T & -\delta_c I
\end{pmatrix}
\begin{pmatrix} d_x \\ -d_y \end{pmatrix}
= -
\begin{pmatrix}
\nabla f^{(k)} - \nabla c^{(k)} y^{(k)} - \mu (X^{(k)})^{-1} e \\
c^{(k)}
\end{pmatrix},
$$

where $X^{(k)} = \text{diag}(x^{(k)})$, $Z^{(k)} = \text{diag}(z^{(k)})$, $e$ is a vector of ones of appropriate size, and $\delta_w$ and $\delta_c$ are primal and dual inertia correction coefficients. The dual direction for the bound constraints is given by $d_z = (X^{(k)})^{-1} (\mu e - Z^{(k)} d_x) - z^{(k)}$. The fraction-to-boundary rule determines primal and dual step lengths that maintain positivity of $x$ and $z$:

$$\alpha_x^{(k)} \stackrel{\text{def}}{=} \max\{\alpha \in (0, 1] \mid x^{(k)} + \alpha d_x \ge (1 - \tau) x^{(k)}\}$$
$$\alpha_z^{(k)} \stackrel{\text{def}}{=} \max\{\alpha \in (0, 1] \mid z^{(k)} + \alpha d_z \ge (1 - \tau) z^{(k)}\},$$

where $\tau$ is a parameter close to 1. A filter line search assesses whether the trial iterate makes sufficient progress with respect to the filter method. If the step length ultimately falls below a given threshold (e.g., $10^{-7}$), we switch to feasibility restoration. Note that by construction the filter entries depend on the barrier parameter $\mu$ through the auxiliary measure $\xi$. Consequently, the filter must be flushed whenever $\mu$ is updated.