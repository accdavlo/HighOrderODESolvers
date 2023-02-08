## Modified Patankar 1st order scheme

def patankar(prod_dest, tspan, u0):
    '''
    Input: prod_dest is the function that returns the matrices p_{i,j}(c) and d_{i,j}(c)
    tspan is the time vector
    u0 is the initial condition
    ''' 
    dim=len(u0)            # Dimension of the problem
    Nt=len(tspan)          # Length of time span
    U=np.zeros((dim,Nt))   # Solution vector
    p=np.zeros((dim,dim))  # Temporary production matrix
    d=np.zeros((dim,dim))  # Temporary destruction matrix
    U[:,0]=u0
    for it in range(1,Nt):   # Loop over timesteps
        dt=tspan[it]-tspan[it-1]    
        p,d =prod_dest(U[:,it-1])   # Computing the production and destruction at the previous timestep
        for i in range(dim):        # Adding all the terms
            lhs = 1.           # Initializing the lhs coefficients
            rhs = U[i,it-1]           # Initializing the rhs
            for j in range(dim):
                lhs = lhs + dt*d[i,j]/U[i,it-1]
                rhs = rhs + dt*p[i,j]
            U[i,it] = rhs/lhs    # Solve the final system
    return tspan, U


## Modified Patankar 1st order scheme

def mPEuler(prod_dest, tspan, u0):
    '''
    Input: prod_dest is the function that returns the matrices p_{i,j}(c) and d_{i,j}(c)
    tspan is the time vector
    u0 is the initial condition
    ''' 
    dim=len(u0)            # Dimension of the problem
    Nt=len(tspan)          # Length of time span
    U=np.zeros((dim,Nt))   # Solution vector
    p=np.zeros((dim,dim))  # Temporary production matrix
    d=np.zeros((dim,dim))  # Temporary destruction matrix
    U[:,0]=u0
    for it in range(1,Nt):   # Loop over timesteps
        dt=tspan[it]-tspan[it-1]    
        p,d =prod_dest(U[:,it-1])   # Computing the production and destruction at the previous timestep
        MM = np.eye(dim)            # Initializing the mass matrix
        for i in range(dim):        # Adding all the terms
            for j in range(dim):
                MM[i,j] = MM[i,j] - dt*p[i,j]/U[j,it-1]
                MM[i,i] = MM[i,i] + dt*d[i,j]/U[i,it-1]
        U[:,it] = np.linalg.solve(MM,U[:,it-1])   # Solve the final system
    return tspan, U