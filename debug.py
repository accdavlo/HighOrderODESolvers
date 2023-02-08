# Loading/installing packages
import numpy as np
# This package allows to  plot
import matplotlib.pyplot as plt 
#This package already implemented some functions for Runge Kutta and multistep methods
from nodepy import rk
#This package contains functions to solve nonlinear problems
from scipy import optimize
from ODEproblems import ODEproblem

## explicit RK method
def explicitRK(flux, tspan, y_0, A, b, c):
    # Solving u'=F(u,t)
    # input: flux=F, tspan is a vector of times determining the RK steps
    # input: y_0 the initial condition
    # input: A,b,c are matrix and vectors of RK methods
    N_time=len(tspan)  # N+1
    dim=len(y_0)          # S
    y=np.zeros((dim,N_time))    # initializing the variable of solutions    
    y[:,0]=y_0                 # first timestep 
    S=np.shape(A)[0]
    u=np.zeros((dim,S))       # Internal stages
    Fu=np.zeros((dim,S))       # Flux at internal stages
    for n in range(N_time-1):    # n=0,..., N-1
        delta_t=tspan[n+1]-tspan[n]
        for k in range(S):
            u[:,k]=y[:,n] 
            for j in range(k):
                u[:,k] =u[:,k]+ delta_t*A[k,j]*Fu[:,j]
            Fu[:,k] = flux(u[:,k],tspan[n]+delta_t*c[k])
        y[:,n+1]=y[:,n]
        for j in range(S):
            y[:,n+1]=y[:,n+1]+delta_t*b[j]*Fu[:,j]
    return tspan, y 

def CrankNicolson(func, jac_func, tspan, y_0):
    '''
    Crank-Nicolson method with a nonlinear solver
    Input:
    func (nonlinear) function of the ODE, takes input u, t
    jac_func jacobian wrt to u of func, takes input u, t
    tspan vector of timesteps (t^0,...,t^N)
    y_0 initial value    
    '''
    N_time=len(tspan)  # N+1
    dim=len(y_0)          # S
    y=np.zeros((dim,N_time))    # initializing the variable of solutions    
    y[:,0]=y_0                 # first timestep 
    for n in range(N_time-1):    # loop through timesteps n=0,..., N-1
        # define the nonlinear function to be solved
        res = lambda yn1: yn1 -y[:,n] -(tspan[n+1]-tspan[n])*0.5*(func(yn1,tspan[n+1])+func(y[:,n],tspan[n]))
        jacRes = lambda yn1: np.eye(dim) - 0.5*(tspan[n+1]-tspan[n])*jac_func(yn1,tspan[n+1])
        z = optimize.root(res, y[:,n], jac=jacRes) # using Newton's method from scipy.optimize
        y[:,n+1] = z.x
    return tspan, y 
import matplotlib.pyplot as plt
pr=ODEproblem("stiff_scalar")
t_span=np.linspace(0,pr.T_fin,500)
rk44=rk.loadRKM('Heun33')
rk.loadRKM()



pr=ODEproblem("nonLinearOscillator")
t_span=np.linspace(0,pr.T_fin,200)
rk44=rk.loadRKM('Heun33')
rk.loadRKM()
tt,uu=explicitRK(pr.flux,t_span,pr.u0,np.float64(rk44.A),np.float64(rk44.b),np.float64(rk44.c))

uu_ex = pr.exact_solution_times(pr.u0,tt)
for j in range(uu.shape[0]):
    plt.plot(tt,uu[j,:],label="Heun")
    plt.plot(tt,uu_ex[j,:],"x",label="exact")

plt.legend()
pr.jacobian(uu[:,7])
plt.show()
