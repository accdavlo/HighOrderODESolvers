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

## explicit RK method
def explicitRelaxRK(flux, y_0, dt0, T_fin, KtMax,  A, b, c):
    # Solving u'=F(u,t)
    # input: flux=F, tspan is a vector of times determining the RK steps
    # input: y_0 the initial condition
    # dt0 is the basic time interval, that will be modified along the steps
    # T_fin is the final time
    # KtMax is maximum number of timesteps
    # input: A,b,c are matrix and vectors of RK methods
    dim=len(y_0)          # S
    y=np.zeros((dim,KtMax))    # initializing the variable of solutions
    tspan=np.zeros(KtMax)      # times will be stored here
    gammas = np.zeros(KtMax)   # Gamma relaxation coefficients
    time= 0.
    gammas[0] = 1
    n=0                        # Time step index
    tspan[0] = time
    y[:,0]=y_0                 # first timestep 
    S=np.shape(A)[0]
    u=np.zeros((dim,S))       # Internal stages
    Fu=np.zeros((dim,S))       # Flux at internal stages
    while(time<T_fin and n<KtMax):    # n=0,..., N-1
        delta_t=min(dt0,T_fin-time)
        #Classic RK step
        for k in range(S):
            u[:,k]=y[:,n] 
            for j in range(k):
                u[:,k] = u[:,k]+ delta_t*A[k,j]*Fu[:,j]
            Fu[:,k] = flux(u[:,k],tspan[n]+delta_t*c[k])
        yn1=y[:,n]
        for j in range(S):
            yn1=yn1+delta_t*b[j]*Fu[:,j]
        # Compute the relaxation gamma
        deltay = yn1-y[:,n]
        sumBScal=0.
        for j in range(S):
            sumBScal=sumBScal + b[j]* np.dot(u[:,j]-y[:,n],Fu[:,j])
        gamma = 2* delta_t* sumBScal/np.dot(deltay,deltay)
        # Update the n+1 values
        y[:,n+1]=y[:,n] +gamma*deltay
        if (time+delta_t<T_fin -10**-16):
            time = time + gamma*delta_t
        else:
            time=T_fin
        tspan[n+1]=time
        gammas[n+1]=gamma
        n=n+1
        
    return tspan[:n+1], y[:,:n+1] , gammas[:n+1]