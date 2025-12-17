## Unification theory

We argue that most Lagrange-Newton methods can be broken down into the following generic ingredients:
* a **constraint relaxation strategy**: a systematic way to relax the general constraints;
* an **inequality handling method**: a systematic way to handle the inequality constraints;
* a **Lagrange-Newton subproblem**: a local Lagrange-Newton approximation of the reformulated problem, composed of:
	* a **Hessian model**: a model of the Lagrangian Hessian of the original problem;
	* an **inertia control strategy**: a strategy to correct the inertia of the Lagrangian Hessian or the augmented system of the reformulated problem;
* a **globalization strategy**: an acceptance test of the trial iterate;
* a **globalization mechanism**: a recourse action upon rejection of the trial iterate.


The following graph gives an overview of state-of-the-art strategies:
<p align="center">
   <img src="docs/figures/wheel.png" alt="Uno hypergraph" width="55%" />
</p>
