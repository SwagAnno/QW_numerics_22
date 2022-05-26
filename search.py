import numpy as np
import matplotlib.pyplot as plt
import qutip as qt

from scipy import optimize as opt

def get_projector(l, i):
    out_mat = np.zeros((l,l))
    out_mat[i,i] = 1

    out = qt.Qobj(out_mat)

    return out

def get_even_state(N):

    out = np.empty(N)

    for i in range(N):
        out[i] = 1/np.sqrt(N)

    return qt.Qobj(out)


def complete_search( N , gamma= None, target = 1, adjacency = True, end = None, show = True) :

    if gamma == None:
        gamma = -1/N

    if end == None:
        end = N

    adj = np.ones((N,N))

    for i in range(N):
        if not adjacency:
            adj[i,i] = 1-N
        else :
            adj[i,i] = 0

    offset = get_projector(N, target)

    mat = gamma*adj - offset

    print(mat)

    psi = get_even_state(N)
    t = np.arange(0,end,.1)

    H = qt.Qobj(mat)

    E_list = []
    E_list.append(get_projector(N,target))

    res = qt.sesolve(H, psi, t, E_list)

    data = res.expect[0]

    if show :  
        fig, ax = plt.subplots()
        ax.plot( t, data)

        ax.vlines( np.pi/2*np.sqrt(N),0,1, colors = "red")
        ax.set_xlabel('t')
        ax.set_ylabel('P_target')
        
        ax.legend("Quantum")

        #display finished plot or pass parameter for firther additions
        plt.show()
        
    else :
        return data

def star_search( N , gamma = None, target = 1, adjacency = False, end = None, show = True) :
    
    if gamma == None:
        gamma = 1/N

    if end == None:
        end = N

    adj = np.zeros((N,N))

    for i in range(N):

        adj[0,i] = -1
        adj[i,0] = -1

        if not adjacency:
            adj[i,i] = 1
        else :
            adj[i,i] = 0

    if not adjacency:
        adj[0,0] = N-1

    offset = get_projector(N, target)

    mat = gamma*adj - offset

    #print(mat)

    psi = get_even_state(N)
    t = np.arange(0,end,.1)
    
    H = qt.Qobj(mat)

    E_list = []
    E_list.append(get_projector(N,target))

    res = qt.sesolve(H, psi, t, E_list)

    data = res.expect[0]
    
    if show :  
        fig, ax = plt.subplots()
        ax.plot( t, data)

        ax.vlines( np.pi/2*np.sqrt(N),0,1, colors = "red")
        ax.set_xlabel('t')
        ax.set_ylabel('P_target')
        
        ax.legend("Quantum")

        #display finished plot or pass parameter for firther additions
        plt.show()
        
    else :
        return data

def star_optimize(N,target, end = None, adjacency = False, show = True) :
    if end == None:
        end = 2*N

    gamma_grid = np.arange(0,.6, .01)
    t_grid = np.arange(0,end, .1)

    data = np.empty( (len(gamma_grid), len(t_grid)))

    for i in range(len(gamma_grid)):
        cur = star_search(N, target = target, gamma = gamma_grid[i],end = end, adjacency = adjacency,show = False)
        data[i,:] = cur


    fig, ax = plt.subplots()

    c = ax.pcolormesh(t_grid, gamma_grid, data)


    ax.set_xlabel('t')
    ax.set_ylabel('gamma')

    fig.colorbar(c, ax=ax)

    plt.show()

def star_lap_vs_adjacency(N,target, gamma = None, end = None):

    if end == None:
        end = N
    
    data_adj = star_search(N, target = target, gamma = .5  ,end = end, adjacency = True ,show = False)
    data_lap = star_search(N, target = target, gamma = 1/N ,end = end, adjacency = False,show = False)

    t = np.arange(0,end, .1)

    fig, ax = plt.subplots()
    ax.plot( t, data_adj, color = "green")
    ax.plot( t, data_lap, color = "orange")

    ax.vlines( np.pi/2*np.sqrt(N),0,1, colors = "red")
    ax.set_xlabel('t')
    ax.set_ylabel('P_target')
    
    ax.legend(["adj","lap"])

    #display finished plot or pass parameter for firther additions
    plt.show()


def t_grover_star(N,target, gamma = None, adjacency = False):

    if gamma == None:
        if adjacency:
            gamma = .5
        else : 
            gamma = 1/N

    data = star_search(N, target = target, gamma = gamma, adjacency = adjacency, end = N, show = False)

    t = np.arange(0,N, .1)
    i = 0

    while  data[i]< data[i+1] :
        i = i+1

    return (t[i]+ t[i+1])/2

def star_search_time_trend(target,bounds = (5,20)):

    N_grid = np.linspace(bounds[0],bounds[1]+1)

    data_adj = np.empty( len(N_grid))
    data_lap = np.empty( len(N_grid))
    data_grover = np.empty( len(N_grid))

    for i in range(len(N_grid)):
        data_adj[i] = t_grover_star(int(N_grid[i]), target = target, adjacency = True)
        data_lap[i] = t_grover_star(int(N_grid[i]), target = target, adjacency = False)
        data_grover[i] = np.pi/2*np.sqrt(N_grid[i])


    fig, ax = plt.subplots()
    ax.plot( N_grid, data_adj, color = "green")
    ax.plot( N_grid, data_lap, color = "orange")
    ax.plot( N_grid, data_grover, color = "red")

    ax.set_xlabel('N')
    ax.set_ylabel('T Grover')
    
    ax.legend(["adj","lap"])

    #display finished plot or pass parameter for firther additions
    plt.show()

star_optimize( 5, target = 1, adjacency = False)
#star_search(50, target = 0, adjacency = False)


#star_lap_vs_adjacency(30, target = 0)
#star_search_time_trend(0)





