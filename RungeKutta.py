## explicit RK method
def explicitRK(flux, tspan, y_0, A, b, c):
    """
    # Solving u'=F(u,t)
    # input: flux=F, tspan is a vector of times determining the RK steps
    # input: y_0 the initial condition
    # input: A,b,c are matrix and vectors of RK methods
    """
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

def RRK(rkm, dt, f, w0=[1.,0], t_final=1., relaxation=True, 
        rescale_step=True, debug=False, gammatol=0.1, print_gamma=False,
        one_step=False, dissip_factor=1.0):
    """
    Relaxation Runge-Kutta method implementation.
    
    Options:
    
        rkm: Base Runge-Kutta method, in Nodepy format
        dt: time step size
        f: RHS of ODE system
        w0: Initial data
        t_final: final solution time
        relaxation: if True, use relaxation method.  Otherwise, use vanilla RK method.
        rescale_step: if True, new time step is t_n + \gamma dt
        debug: output some additional diagnostics
        gammatol: Fail if abs(1-gamma) exceeds this value
        
    """
    w = np.array(w0)
    t = 0
    # We pre-allocate extra space because if rescale_step==True then
    # we don't know exactly how many steps we will take.
    ww = np.zeros([len(w0),int((t_final-t)/dt*2.5)+10000])
    ww[:,0] = w.copy()
    tt = [t]
    ii = 0
    s = len(rkm)
    b = rkm.b
    y = np.zeros((s,len(w0)))
    max_gammam1 = 0.
    gams = []
    
    while t < t_final:
        if t + dt >= t_final:
            dt = t_final - t # Hit final time exactly
        
        for i in range(s):
            y[i,:] = w.copy()
            for j in range(i):
                y[i,:] += rkm.A[i,j]*dt*f(y[j,:])
                
        F = np.array([f(y[i,:]) for i in range(s)])
        
        if relaxation:
            numer = 2*sum(b[i]*rkm.A[i,j]*np.dot(F[i],F[j]) \
                                for i in range(s) for j in range(s))
            denom = sum(b[i]*b[j]*np.dot(F[i],F[j]) for i in range(s) for j in range(s))
            if denom != 0:
                gam = numer/denom
            else:
                gam = 1.
        else:  # Use standard RK method
            gam = 1.
           
        if print_gamma:
            print(gam)
        
        if np.abs(gam-1.) > gammatol:
            print(gam)
            raise Exception("The time step is probably too large.")
        
        w = w + dissip_factor*gam*dt*sum([b[j]*F[j] for j in range(s)])
        if (t+dt < t_final) and rescale_step:
            t += dissip_factor*gam*dt
        else:
            t += dt
        ii += 1
        tt.append(t)
        ww[:,ii] = w.copy()
        if debug:
            gm1 = np.abs(1.-gam)
            max_gammam1 = max(max_gammam1,gm1)
            gams.append(gam)
            
        if one_step:
            return w, gam
            
    if debug:
        print(max_gammam1)
        return tt, ww[:, :ii+1], np.array(gams)
    else:
        return tt, ww[:,:ii+1]