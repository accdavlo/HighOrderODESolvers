import numpy as np
import matplotlib.pyplot as plt
from numpy.polynomial.legendre import leggauss
from quadr import lglnodes,equispaced
from scipy.interpolate import lagrange

def lagrange_poly(nodes,k):
    interpVal=np.zeros(np.size(nodes))
    interpVal[k] = 1.
    pp=lagrange(nodes,interpVal)
    return pp

def lagrange_basis(nodes,x,k):
    pp=lagrange_poly(nodes,k)
    return pp(x)

def lagrange_deriv(nodes,x,k):
    pp=lagrange_poly(nodes,k)
    dd=pp.deriv()
    return dd(x)

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
        
def getADER_matrix(order, nodes_type):
    nodes_poly, w_poly = get_nodes(order,nodes_type)
    if nodes_type=="equispaced":
        quad_order=order
        nodes_quad, w = get_nodes(quad_order,"gaussLegendre")
    else:
        quad_order=order
        nodes_quad, w = get_nodes(quad_order,nodes_type)
                    
    # generate mass matrix
    M = np.zeros((order,order))
    for i in range(order):
        for j in range(order):
            M[i,j] = lagrange_basis(nodes_poly,1.0,i) *lagrange_basis(nodes_poly,1.0,j)\
                  -sum([lagrange_deriv(nodes_poly,nodes_quad[q],i)\
                  *lagrange_basis(nodes_poly,nodes_quad[q],j)\
                  *w[q] for q in range(quad_order)])
    # generate mass matrix
    RHSmat = np.zeros((order,order))
    for i in range(order):
        for j in range(order):
            RHSmat[i,j] = sum([lagrange_basis(nodes_poly,nodes_quad[q],i)*\
                               lagrange_basis(nodes_poly,nodes_quad[q],j)*\
                               w[q] for q in range(quad_order)])
    return nodes_poly, w_poly, M, RHSmat

def ader(func, tspan, y_0, M_sub, K_corr, distribution):
    N_time=len(tspan)
    dim=len(y_0)
    U=np.zeros((dim, N_time))
    u_p=np.zeros((dim, M_sub+1))
    u_a=np.zeros((dim, M_sub+1))
    u_tn=np.zeros((dim, M_sub+1))
    rhs= np.zeros((dim,M_sub+1))
    
    x_poly, w_poly, ADER, RHS_mat = getADER_matrix(M_sub+1, distribution)
    invader = np.linalg.inv(ADER)
    evolMatrix=np.matmul(invader,RHS_mat)
    
    U[:,0]=y_0
    
    for it in range(1, N_time):
        delta_t=(tspan[it]-tspan[it-1])
        for m in range(M_sub+1):
            u_a[:,m]=U[:,it-1]
            u_p[:,m]=U[:,it-1]
            u_tn[:,m]=U[:,it-1]
        for k in range(1,K_corr+1):
            u_p=np.copy(u_a)
            for r in range(M_sub+1):
                rhs[:,r]=func(u_p[:,r])
            for d in range(dim):
                u_a[d,:] = u_tn[d,:] + delta_t*np.matmul(evolMatrix,rhs[d,:])
        for d in range(dim):
            U[d,it]=sum(u_a[d,:]*[lagrange_basis(x_poly,1.0,i) for i in range(M_sub+1)])
    return tspan, U