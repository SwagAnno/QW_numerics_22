import numpy as np
import matplotlib.pyplot as plt
import qutip as qt

def get_projector(l, i):
    out_mat = np.zeros((l,l))
    out_mat[i,i] = 1

    out = qt.Qobj(out_mat)

    return out

def classic_line(step, show = False, mode = None):
    prob = np.zeros(step*2 +1)
    prob[step] = 1

    coin = [.5,.5]

    for j in range(step):

        next_p = np.zeros(step*2 +1)

        for m in range(2*step):
            next_p[m] += coin[0]* prob[m + 1]
            next_p[m+1] += coin[1]* prob[m]

        prob = next_p

    #print( prob[0], prob[-1], prob[1])
    print(np.sum(prob))

    if show :
        fig, ax = plt.subplots()
        ax.plot( np.arange(-step, step+1,1), prob)
        ax.set_xlabel('N')
        ax.set_ylabel('p_n')
        
        ax.legend("Classic")

        #display finished plot or pass parameter for firther additions
        plt.show()
    if mode == "data":
        return data
    else :
        return line_momenta(step,prob)       

def quantum_line(step, show = False, mode = None):
    A = np.zeros((step*2 + 1,2), dtype = complex)
    A[step,0] = 1/np.sqrt(2)
    A[step,1] = 1j*1/np.sqrt(2)

    coin = np.ones((2,2), dtype = complex)
    coin[1,1] = -1
    coin = coin * 1/np.sqrt(2)

    #print(coin)
    #print( np.dot(coin,A[step,]))

    for j in range(step):

        next_A = np.zeros((step*2 +1,2), dtype = complex)

        for m in range(2*step):
            next_A[m,0] += np.dot(coin,A[m+1,])[0]
            next_A[m+1,1] += np.dot(coin,A[m,])[1]
            
        A = next_A

    #print( prob[0], prob[-1], prob[1])
    prob = A.conjugate()
    np.multiply(A,prob,out=prob)
    prob = np.abs(prob)
    #print(prob)
    print(np.sum(prob))

    data = prob[:,0]+prob[:,1]

    if show :  
        fig, ax = plt.subplots()
        ax.plot( np.arange(-step, step+1,1), data)
        ax.set_xlabel('N')
        ax.set_ylabel('p_n')
        
        ax.legend("Quantum")

        #display finished plot or pass parameter for firther additions
        plt.show()
        
    if mode == "data":
        return data
    else :
        return line_momenta(step,data)

def CTQW_line(t, length, show = False):
    psi = qt.basis(2*length +1,length)
    t = [0.,t]

    H = np.zeros( (2*length +1,2*length +1) )
    for i in range(2*length):
        H[i,i+1] = -1
        H[i+1,i] = -1
    gamma = .35
    H = H*gamma
    H = qt.Qobj(H)

    E_list = []
    for i in range(2*length+1):
        E_list.append(get_projector(2*length+1,i))

    res = qt.sesolve(H, psi, t, E_list)

    #print(res.expect[0])

    data = np.ndarray(2*length +1)

    for i in range(2*length +1):
        data[i] = res.expect[i][1]
    
    if show :  
        fig, ax = plt.subplots()
        ax.plot( np.arange(-length, length+1,1), data)
        ax.set_xlabel('N')
        ax.set_ylabel('p_n')
        
        ax.legend("Quantum")

        #display finished plot or pass parameter for firther additions
        plt.show()
        
    else :
        return data

def line_momenta(step, p):

    mean = 0
    var = 0
    for i in range(2*step +1):
            mean += (i - step)*p[i]
            var += (i - step)*(i - step)*p[i]

    return mean,var

def plot_var(step):

    quantum_m = np.ndarray( (step,2))
    class_m = np.ndarray( (step,2))

    for i in range(step) :
        quantum_m[i,:] = quantum_line(i)
        class_m[i,:] = classic_line(i)

    fig, ax = plt.subplots()
    ax.plot( np.arange(0, step,1), quantum_m[:,1])
    ax.plot( np.arange(0, step,1), class_m[:,1])
    ax.set_xlabel('step')
    ax.set_ylabel('var')
    
    ax.legend(["Quantum","Class"])

    plt.show()

def compare_dc(step):

    disc = quantum_line(step, mode = "data")
    cont = CTQW_line(step,step)

    fig, ax = plt.subplots()
    ax.plot( np.arange(-step, step+1,1), disc)
    ax.plot( np.arange(-step, step+1,1), cont)
    ax.set_xlabel('n')
    ax.set_ylabel('p')
    
    ax.legend(["Discrete","CTQW"])

    plt.show()
    

#print(line_momenta(1,[1,1,1]))
compare_dc(80)

    
