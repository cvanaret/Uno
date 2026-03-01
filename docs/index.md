Uno, a unified solver for nonlinearly constrained optimization
==============================================================

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
$$

where $x \in \mathbb{R}^n$, $f : \mathbb{R}^n \to \mathbb{R}$, and $c: \mathbb{R}^n \to \mathbb{R}^m$.

Most derivative-based iterative methods for nonlinearly constrained nonconvex optimization share common algorithmic components. Building on this insight, we introduce an abstract framework structured around generic building blocks that describes these methods in a unified fashion. We then present **Uno** (Unifying Nonlinear Optimization)[^1], a modular open-source solver for nonlinearly constrained optimization that unifies numerous state-of-the-art methods and organizes existing strategies into a coherent hierarchy. Optimization strategies are implemented as unique software components (e.g., all line-search methods may be generated based on a single `LineSearch` class); these efficient and robust implementations can be combined at will and compared on a given instance.

Uno allows the automatic generation of various strategy combinations on the fly with no programming effort from the user. While all combinations do not lead to convergent methods, some of them result in efficient solvers that may not exist as software implementations. We demonstrate that Uno is competitive against state-of-the-art solvers on a subset of 429 CUTE test problems translated into AMPL, while being extensible and lightweight. We believe that Uno has the potential to serve as an experimental laboratory for practitioners and optimizers and to accelerate research in nonlinearly constrained optimization. Our ultimate goal is to promote the extension of state-of-the-art nonlinear optimization techniques to new classes of problems such as problems with equilibrium constraints and nonlinear robust optimization.

[^1]: Uno was first introduced at the ISMP 2018 conference under the name Argonot.