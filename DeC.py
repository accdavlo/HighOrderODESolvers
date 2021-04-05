import numpy as np
import matplotlib.pyplot as plt
from numpy.polynomial.legendre import leggauss
from quadr import lglnodes,equispaced
def lagrange_basis(nodes,x,k):
    y=np.zeros(x.size)
    for ix, xi in enumerate(x):
        tmp=[(xi-nodes[j])/(nodes[k]-nodes[j])  for j in range(len(nodes)) if j!=k]
        y[ix]=np.prod(tmp)
    return y

def get_nodes(order,nodes_type):
    if nodes_type=="equispaced":
        nodes,w = equispaced(order)
    elif nodes_type == "gaussLegendre":
        nodes,w = leggauss(order)
    elif nodes_type == "gaussLobatto":
        nodes, w = lglnodes(order-1,10**-15)
    nodes=nodes*0.5+0.5
    w = w*0.5
    return nodes, w
        
def compute_theta_DeC(order, nodes_type):
    nodes, w = get_nodes(order,nodes_type)
    int_nodes, int_w = get_nodes(order,"gaussLobatto")
    # generate theta coefficients 
    theta = np.zeros((order,order))
    beta = np.zeros(order)
    for m in range(order):
        beta[m] = nodes[m]
        nodes_m = int_nodes*(nodes[m])
        w_m = int_w*(nodes[m])
        for r in range(order):
            theta[r,m] = sum(lagrange_basis(nodes,nodes_m,r)*w_m)
    return theta, beta


def compute_RK_from_DeC(M_sub,K_corr,nodes_type):
    order=M_sub+1;
    [theta,beta]=compute_theta_DeC(order,nodes_type)
    bar_beta=beta[1:]  # M_sub
    bar_theta=theta[:,1:].transpose() # M_sub x (M_sub +1)
    theta0= bar_theta[:,0]  # M_sub x 1
    bar_theta= bar_theta[:,1:] #M_sub x M_sub
    A=np.zeros((M_sub*(K_corr-1)+1,M_sub*(K_corr-1)+1))  # (M_sub x K_corr +1)^2
    b=np.zeros(M_sub*(K_corr-1)+1)
    c=np.zeros(M_sub*(K_corr-1)+1)

    c[1:M_sub+1]=bar_beta
    A[1:M_sub+1,0]=bar_beta
    for k in range(1,K_corr-1):
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


def dec(func, tspan, y_0, M_sub, K_corr, distribution):
    N_time=len(tspan)
    dim=len(y_0)
    U=np.zeros((dim, N_time))
    u_p=np.zeros((dim, M_sub+1))
    u_a=np.zeros((dim, M_sub+1))
    rhs= np.zeros((dim,M_sub+1))
    Theta, beta = compute_theta_DeC(M_sub+1,distribution)
    U[:,0]=y_0
    for it in range(1, N_time):
        delta_t=(tspan[it]-tspan[it-1])
        for m in range(M_sub+1):
            u_a[:,m]=U[:,it-1]
            u_p[:,m]=U[:,it-1]
        for k in range(1,K_corr+1):
            u_p=np.copy(u_a)
            for r in range(M_sub+1):
                rhs[:,r]=func(u_p[:,r])
            for m in range(1,M_sub+1):
                u_a[:,m]= U[:,it-1]+delta_t*sum([Theta[r,m]*rhs[:,r] for r in range(M_sub+1)])
        U[:,it]=u_a[:,M_sub]
    return tspan, U
            
def decImplicit(func,jac_stiff, tspan, y_0, M_sub, K_corr, distribution):
    N_time=len(tspan)
    dim=len(y_0)
    U=np.zeros((dim, N_time))
    u_p=np.zeros((dim, M_sub+1))
    u_a=np.zeros((dim, M_sub+1))
    u_help= np.zeros(dim)
    rhs= np.zeros((dim,M_sub+1))
    Theta, beta = compute_theta_DeC(M_sub+1,distribution)
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
        for k in range(1,K_corr+1):
            u_p=np.copy(u_a)
            for r in range(M_sub+1):
                rhs[:,r]=func(u_p[:,r])
            for m in range(1,M_sub+1):
                u_a[:,m]= u_p[:,m]+delta_t*np.matmul(invJac[m,:,:],\
                (-(u_p[:,m]-u_p[:,0])/delta_t\
                 +sum([Theta[r,m]*rhs[:,r] for r in range(M_sub+1)])))
        U[:,it]=u_a[:,M_sub]
    return tspan, U



def decMPatankar(prod_dest, rhs, tspan, y_0, M_sub, K_corr, distribution):
    N_time=len(tspan)
    dim=len(y_0)
    U=np.zeros((dim, N_time))
    u_p=np.zeros((dim, M_sub+1))
    u_a=np.zeros((dim, M_sub+1))
    prod_p = np.zeros((dim,dim,M_sub+1))
    dest_p = np.zeros((dim,dim,M_sub+1))
    rhs_p= np.zeros((dim,M_sub+1))
    Theta, beta = compute_theta_DeC(M_sub+1,distribution)
    U[:,0]=y_0
    for it in range(1, N_time):
        delta_t=(tspan[it]-tspan[it-1])
        for m in range(M_sub+1):
            u_a[:,m]=U[:,it-1]
            u_p[:,m]=U[:,it-1]
        for k in range(1,K_corr+1):
            u_p=np.copy(u_a)
            for r in range(M_sub+1):
                prod_p[:,:,r], dest_p[:,:,r]=prod_dest(u_p[:,r])
                rhs_p[:,r]=rhs(u_p[:,r])
            for m in range(1,M_sub+1):
                u_a[:,m]= patankar_type_dec(prod_p,dest_p,rhs_p,delta_t,m,M_sub,Theta,u_p,dim)
        U[:,it]=u_a[:,M_sub]
    return tspan, U


def patankar_type_dec(prod_p,dest_p,rhs_p,delta_t,m,M_sub,Theta,u_p,dim):
    mass= np.eye(dim)
    RHS= u_p[:,0]
    for i in range(dim):
        for r in range(M_sub+1):
            RHS[i]=RHS[i]+delta_t*Theta[r,m]*rhs_p[i,r]
            if Theta[r,m]>0:
                for j in range(dim):
                    mass[i,j]=mass[i,j]-delta_t*Theta[r,m]*(prod_p[i,j,r]/u_p[j,m])
                    mass[i,i]=mass[i,i]+ delta_t*Theta[r,m]*(dest_p[i,j,r]/u_p[i,m])
            elif Theta[r,m]<0:
                for j in range(dim):
                    mass[i,i]=mass[i,i]- delta_t*Theta[r,m]*(prod_p[i,j,r]/u_p[i,m])
                    mass[i,j]=mass[i,j]+ delta_t*Theta[r,m]*(dest_p[i,j,r]/u_p[j,m])
    return np.linalg.solve(mass,RHS)
