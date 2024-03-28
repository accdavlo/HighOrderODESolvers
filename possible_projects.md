# Possible projects
  * Stability (Minion) and/or code DeC IMEX for FD advection/diffusion (with a fixed coefficient, letting only dt/dx change?)
  * Relaxation ADER/DeC applications on Burgers smooth
  * Implicit ADER as RK and IMEX ADER
  * Modified Patankar DeC for a simple Shallow Water problem
  * Efficient DeC: code, write RK matrix
  * Modified Patankar - ADER
  * ADER with L1 that is diagonal: stability and accuracy wrt classical one, explicit, implicit, IMEX
  * IMEX application for complicated ODE: inverse pendulum with forcing term to stabilize the problem
  * Prove positivity of implicit Euler for nonlinear PDS + ...

## ODE solvers and time performances for generative models -> Guglielmo Padula 

## Relaxation ADER/DeC applications on Burgers smooth -> Dario Coscia/Rashid Ashraf

Study the discretization of the (in)viscid Burgers' equation

$\partial_t u + \partial_x (u^2/2) - \varepsilon \partial_{xx} u=0,$

given by a FD formula on a 1D periodic domain

$\partial_t u_j = \frac{F_{j+1/2}-F_{j-1/2}}{\Delta x}$

with

$F_{j+1/2} = \frac{u_{j+1}^2 + u_{j}u_{j+1} + u_j^2}{6} -\frac{\varepsilon}{\Delta x} (u_{j+1}-u_j).$

1. Check that the discretization is consistent.
2. Chech that for the inviscid case ($\varepsilon =0$) the spatial discretization is energy preserving and in the viscous case ($\varepsilon >0$) the spatial discretization is energy dissipative.
3. Implement the method for some (not too) high order time discretization with/without relaxation (ADER/DeC) and check the energy production.
4. Test on some problems: in the inviscid case try to stop before the shock formation and see what happens when you go further and try to understand how much dissipation you need to make the scheme stable.

## Stability for IMEX schemes with spectral methods -> Nicola Clinco

Take the advection diffusion equation
$\partial_t u + a \partial_x u + d \partial_{xx} u = 0$
consider the Fourier spectral method on it, considering $u = \sum_k \hat{u}_k e^{ixk}$. 
The weak formulation generated has a very simple matrix structure with diagonal matrices.

1. Study the stability of various IMEX methods with the Minion approach.
2. Apply the analysis to problem using as parameters $\Delta t a$ and $\Delta t d$ varying $N/2$ the maximum $k$ (theoretically and/or computationally).
3. Compare the numerical solutions with the found results to check whether there is a correspondence.
4. Nonlinear? I have an option, but we can discuss if you have something in mind.
 

## ODE solvers comparison for Entropic Wasserstein Gradient Flows/JKO algorithm -> Isabella Gonnella

