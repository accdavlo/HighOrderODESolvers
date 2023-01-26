import numpy as np
import matplotlib.pyplot as plt
from numpy.polynomial.legendre import leggauss
from quadr import lglnodes,equispaced

def lagrange_basis(nodes,x,k):
    """
    Compute the Lagrange polynomial basis for the given set of nodes and x values
    Parameters:
        nodes (np.array): set of nodes
        x (np.array): x values
        k (int): index of node
    Returns:
        y (np.array): y values of the Lagrange polynomial basis
    """
    y=np.zeros(x.size)
    for ix, xi in enumerate(x):
        tmp=[(xi-nodes[j])/(nodes[k]-nodes[j])  for j in range(len(nodes)) if j!=k]
        y[ix]=np.prod(tmp)
    return y

def get_nodes(n_nodes,nodes_type):
    """
    Compute the quadrature nodes and weights for the specified type of nodes
    Parameters:
        n_nodes (int): number of the nodes
        nodes_type (str): type of nodes. Can be "equispaced", "gaussLegendre", "gaussLobatto"
    Returns:
        nodes (np.array): the computed nodes on [0,1]
        w (np.array): the computed weights
    """
    if nodes_type=="equispaced":
        nodes,w = equispaced(n_nodes)
    elif nodes_type == "gaussLegendre":
        nodes,w = leggauss(n_nodes)
    elif nodes_type == "gaussLobatto":
        nodes, w = lglnodes(n_nodes-1,10**-15)
    nodes=nodes*0.5+0.5
    w = w*0.5
    return nodes, w
        
def compute_theta_DeC(n_nodes, nodes_type):
    """
    Compute the DeC theta coefficients and beta values for the given number of nodes n_nodes and type of nodes
    Parameters:
        n_nodes (int): number of the nodes
        nodes_type (str): type of nodes. Can be "equispaced", "gaussLegendre", "gaussLobatto"
    Returns:
        theta (np.array): the computed theta coefficients of shape (n_nodes, n_nodes)
        beta (np.array): the computed beta values of shape (n_nodes)
    """
    nodes, w = get_nodes(n_nodes,nodes_type)
    int_nodes, int_w = get_nodes(n_nodes,"gaussLobatto")
    # generate theta coefficients 
    theta = np.zeros((n_nodes,n_nodes))
    beta = np.zeros(n_nodes)
    for m in range(n_nodes):
        beta[m] = nodes[m]
        nodes_m = int_nodes*(nodes[m])
        w_m = int_w*(nodes[m])
        for r in range(n_nodes):
            theta[r,m] = sum(lagrange_basis(nodes,nodes_m,r)*w_m)
    return theta, beta


def compute_RK_from_DeC(M_sub,K_iter,nodes_type):
    """
    Compute the Runge--Kutta matrices for a certain DeC method given by the number of subtimesteps, the number of iterations K_iter, and type of nodes
    Parameters:
        M_sub (int): number subtimesteps = n_nodes-1
        K_iter (int): number of DeC iterations
        nodes_type (str): type of nodes. Can be "equispaced", "gaussLegendre", "gaussLobatto"
    Returns:
        A (np.array): the A RK matrix of shape (NRK, NRK)
        beta (np.array): the b RK coefficients of shape (NRK)
        c (np.array): the c RK coefficients of shape (NRK)
    """
    n_nodes=M_sub+1
    [theta,beta]=compute_theta_DeC(n_nodes,nodes_type)
    bar_beta=beta[1:]  # M_sub
    bar_theta=theta[:,1:].transpose() # M_sub x (M_sub +1)
    theta0= bar_theta[:,0]  # M_sub x 1
    bar_theta= bar_theta[:,1:] #M_sub x M_sub
    A=np.zeros((M_sub*(K_iter-1)+1,M_sub*(K_iter-1)+1))  # (M_sub x K_iter +1)^2
    b=np.zeros(M_sub*(K_iter-1)+1)
    c=np.zeros(M_sub*(K_iter-1)+1)

    c[1:M_sub+1]=bar_beta
    A[1:M_sub+1,0]=bar_beta
    for k in range(1,K_iter-1):
        r0=1+M_sub*k
        r1=1+M_sub*(k+1)
        c0=1+M_sub*(k-1)
        c1=1+M_sub*(k)
        c[r0:r1]=bar_beta
        A[r0:r1,0]=theta0
        A[r0:r1,c0:c1]=bar_theta
    b[0]=theta0[-1]
    b[-M_sub:]=bar_theta[M_sub-1,:]
    return A,b,c


def dec(func, tspan, y_0, M_sub, K_iter, nodes_type):
    """
    Deferred Correction time integrator on the ode u'=func(u)
    The DeC is defined with L^{1,m}(U)=U^m-U^0-dt*func(U^0) and L^{2,m} = U^m-U^0-dt*\sum_r theta^m_r func(U^r)
    Parameters:
        func: right hand side of the ode, it is a function of u
        tspan (np.array): the array of the times at which we want to compute the solution of size N_time
        y0 (np.array): initial conditions of size dim
        M_sub (int): number subtimesteps = n_nodes-1
        K_iter (int): number of DeC iterations
        nodes_type (str): type of nodes. Can be "equispaced", "gaussLegendre", "gaussLobatto"
    Returns:        
        tspan (np.array): the array of the times at which we want to compute the solution of size N_time
        U (np.array): the solution at times tspan of shape (dim, N_time)
    """
    N_time=len(tspan)
    dim=len(y_0)
    U=np.zeros((dim, N_time))
    u_p=np.zeros((dim, M_sub+1))
    u_a=np.zeros((dim, M_sub+1))
    rhs= np.zeros((dim,M_sub+1))
    Theta, beta = compute_theta_DeC(M_sub+1,nodes_type)
    U[:,0]=y_0
    for it in range(1, N_time): # Iterating on time
        delta_t=(tspan[it]-tspan[it-1])
        for m in range(M_sub+1): # Initializing u at all subtimenodes
            u_a[:,m]=U[:,it-1]
            u_p[:,m]=U[:,it-1]
        for k in range(1,K_iter+1): # Iterating on the DeC corrections
            u_p=np.copy(u_a)
            for r in range(M_sub+1): # Computing the rhs function at the previous iteration
                rhs[:,r]=func(u_p[:,r])
            for m in range(1,M_sub+1): # updating the value at all subtimenodes
                u_a[:,m]= U[:,it-1]+delta_t*sum([Theta[r,m]*rhs[:,r] for r in range(M_sub+1)])
        U[:,it]=u_a[:,M_sub]
    return tspan, U
            
def decImplicit(func,jac_stiff, tspan, y_0, M_sub, K_iter, nodes_type):
    """
    Deferred Correction time integrator on the stiff ode u'=func(u) with jac_stiff is the jacobian of the stiff part of func
    The DeC is defined with L^{1,m}(U)=U^m-U^0-dt*jac_stiff(U^0)*(U^m-U^0) and L^{2,m} = U^m-U^0-dt*\sum_r theta^m_r func(U^r)
    Parameters:
        func: right hand side of the ode, it is a function of u
        jac_stiff: jacobian of the stiff part of func, it is a function of u
        tspan (np.array): the array of the times at which we want to compute the solution of size N_time
        y0 (np.array): initial conditions of size dim
        M_sub (int): number subtimesteps = n_nodes-1
        K_iter (int): number of DeC iterations
        nodes_type (str): type of nodes. Can be "equispaced", "gaussLegendre", "gaussLobatto"
    Returns:        
        tspan (np.array): the array of the times at which we want to compute the solution of size N_time
        U (np.array): the solution at times tspan of shape (dim, N_time)
    """
    N_time=len(tspan)
    dim=len(y_0)
    U=np.zeros((dim, N_time))
    u_p=np.zeros((dim, M_sub+1))
    u_a=np.zeros((dim, M_sub+1))
    u_help= np.zeros(dim)
    rhs= np.zeros((dim,M_sub+1))
    Theta, beta = compute_theta_DeC(M_sub+1,nodes_type)
    invJac=np.zeros((M_sub+1,dim,dim))
    U[:,0]=y_0
    for it in range(1, N_time):
        delta_t=(tspan[it]-tspan[it-1])
        for m in range(M_sub+1):
            u_a[:,m]=U[:,it-1]
            u_p[:,m]=U[:,it-1]
        SS=jac_stiff(u_p[:,0])
        for m in range(1,M_sub+1):
            invJac[m,:,:]=np.linalg.inv(np.eye(dim) - delta_t*beta[m]*SS)
        for k in range(1,K_iter+1):
            u_p=np.copy(u_a)
            for r in range(M_sub+1):
                rhs[:,r]=func(u_p[:,r])
            for m in range(1,M_sub+1):
                u_a[:,m]= u_p[:,m]+delta_t*np.matmul(invJac[m,:,:],\
                (-(u_p[:,m]-u_p[:,0])/delta_t\
                 +sum([Theta[r,m]*rhs[:,r] for r in range(M_sub+1)])))
        U[:,it]=u_a[:,M_sub]
    return tspan, U



def decMPatankar(prod_dest, rhs, tspan, y_0, M_sub, K_iter, nodes_type):
    """
    Modified Patankar Deferred Correction time integrator on the Production-Destruction ODE (u_i)'=sum_j p_ij(u)-d_ij(u) + rhs_i(u)
    with nonnegative production and destruction terms p_ij,d_ij >= 0 and nonnegative rest term rhs_i(u)
    Parameters:
        prod_dest: production and destruction functions depending on u, it returns two matrices of shape (dim,dim)
        rhs: nonnegative rest term function of u, which returns a term of shape (dim)
        tspan (np.array): the array of the times at which we want to compute the solution of size N_time
        y0 (np.array): initial conditions of size dim
        M_sub (int): number subtimesteps = n_nodes-1
        K_iter (int): number of DeC iterations
        nodes_type (str): type of nodes. Can be "equispaced", "gaussLegendre", "gaussLobatto"
    Returns:        
        tspan (np.array): the array of the times at which we want to compute the solution of size N_time
        U (np.array): the solution at times tspan of shape (dim, N_time)
    """
    N_time=len(tspan)
    dim=len(y_0)
    U=np.zeros((dim, N_time))
    u_p=np.zeros((dim, M_sub+1))
    u_a=np.zeros((dim, M_sub+1))
    prod_p = np.zeros((dim,dim,M_sub+1))
    dest_p = np.zeros((dim,dim,M_sub+1))
    rhs_p= np.zeros((dim,M_sub+1))
    Theta, beta = compute_theta_DeC(M_sub+1,nodes_type)
    U[:,0]=y_0
    for it in range(1, N_time): # loop on time steps
        delta_t=(tspan[it]-tspan[it-1])
        for m in range(M_sub+1):
            u_a[:,m]=U[:,it-1]
            u_p[:,m]=U[:,it-1]
        for k in range(1,K_iter+1): # loop on DeC iterations
            u_p=np.copy(u_a)
            for r in range(M_sub+1): # compute production destruction and rest terms on the previous iteration
                prod_p[:,:,r], dest_p[:,:,r]=prod_dest(u_p[:,r])
                rhs_p[:,r]=rhs(u_p[:,r])
            for m in range(1,M_sub+1): # updates of u at all subtimenodes using modified Patankar technique
                u_a[:,m]= patankar_type_dec(prod_p,dest_p,rhs_p,delta_t,m,M_sub,Theta,u_p,dim)
        U[:,it]=u_a[:,M_sub]
    return tspan, U


def patankar_type_dec(prod_p,dest_p,rhs_p,delta_t,m,M_sub,Theta,u_p,dim):
    """
    Modified Patankar Deferred Correction step: assembles the matrix and solves the system
    Parameters:
        prod_p (np.array): production matrices of previous iteration at all subtimenodes of shape (dim,dim,M_sub+1)
        dest_p (np.array): destruction matrices of previous iteration at all subtimenodes of shape (dim,dim,M_sub+1)
        rhs_p: nonnegative rest term function of u at all subtimenodes of shape (dim,M_sub+1)
        delta_t: timestep
        m (int): current subtimenode
        M_sub (int): number subtimesteps = n_nodes-1
        Theta (np.array): DeC theta coefficients
        u_p (np.array): u at previous itertion at all subtimenodes of shape (dim,M_sub+1)
        dim (int): dimension of u
    Returns:        
        u_a^m (np.array): the solution of the next iteration at subtimenode m of shape (dim)
    """
    mass = np.eye(dim) # setting up the matrix of the linear system
    RHS  = u_p[:,0]    # setting up the RHS of the linear system
    for i in range(dim):
        for r in range(M_sub+1): #RHS = u^0 + sum_r dt*theta^m_r rest(u^r)
            RHS[i]=RHS[i]+delta_t*Theta[r,m]*rhs_p[i,r]
            if Theta[r,m]>0:  # matrix assembled according to the sign of theta_r^m
                for j in range(dim):
                    mass[i,j]=mass[i,j]-delta_t*Theta[r,m]*(prod_p[i,j,r]/u_p[j,m])
                    mass[i,i]=mass[i,i]+ delta_t*Theta[r,m]*(dest_p[i,j,r]/u_p[i,m])
            elif Theta[r,m]<0:
                for j in range(dim):
                    mass[i,i]=mass[i,i]- delta_t*Theta[r,m]*(prod_p[i,j,r]/u_p[i,m])
                    mass[i,j]=mass[i,j]+ delta_t*Theta[r,m]*(dest_p[i,j,r]/u_p[j,m])
    return np.linalg.solve(mass,RHS) # solution of linear system 
