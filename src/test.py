import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from scipy.linalg import eigh, det

def jacobi_transform(m_list):
    dim=m_list.shape[0]
    J=np.zeros((dim,dim))
    for i in range(dim):
        sum_m=np.sum(m_list[:i+1])
        for j in range(dim):
            if j==i+1:
                J[i,j]=-1
            elif i+1<j:
                J[i,j]=0
            else:
                J[i,j]=m_list[j]/sum_m
            if np.isnan(J[i,j]):
                J[i,j]=1
    U=np.linalg.inv(J)
    if 1<dim:
        U=np.delete(U,dim-1,1)
        J=np.delete(J,dim-1,0)
    return J,U

def w_gen_3():
    w1=np.zeros((3))
    w2=np.zeros((3))
    w3=np.zeros((3))
    w1[0]=1
    w1[1]=-1
    w2[0]=1
    w2[2]=-1
    w3[1]=1
    w3[2]=-1
    w_list=[w1,w2,w3]
    return w_list

def corput(n, b=3):
    q=0
    bk=1/b
    while n>0:
        n, rem= np.divmod(n,b)
        q += rem*bk
        bk /= b
    return q

def halton(n,d):
    x=np.zeros(d)
    base=np.array([2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,101,103,107,109,113,127,131,137,139,149,151,157,163,167,173,179,181,191,193,197,199,211,223,227,233,239,241,251,257,263,269,271,277,281])
    assert base.shape[0]>d, "Error: d exceeds the number of basis elements."
    for i in range(d):
        x[i]=corput(n,base[i])
    return x

def A_generate(bij,w_list):
    if type(bij)==list or type(bij)==np.ndarray:
        dim=len(w_list)
        mat_list=[np.outer(w_list[i],w_list[i]) for i in range(dim)]
        for i in range(dim):
            mat_list[i]=mat_list[i]/(bij[i]**2)
        A=sum(mat_list)
        return A
    else:
        dim=len(w_list)
        mat_list=[np.outer(w_list[i],w_list[i]) for i in range(dim)]
        for i in range(dim):
            mat_list[i]=mat_list[i]/(bij**2)
        A=sum(mat_list)
        return A

def shift_dot(a,b,mat=None):
    n=a.shape[1]
    sum=0
    if not hasattr(mat, "__len__"):
        mat=np.identity(n)
    assert n==np.shape[0], "ERROR! Matrix shape does not match number of shift vectors."
    for i in range(n):
        for j in range(n):
            dot=a[:,i]@b[:,j]
            sum +=mat[i,j]*dot
    return sum

def transform_list(alphas):
    g_new=[np.ones((1,1))*alphas[i] for i in range(len(alphas))]
    print("size=",np.shape(g_new))
    return g_new

def w_gen(dim, i,j):  
    if dim==1:
        w=np.ones((1,1))
        return w
    else:
        w=np.zeros((1,dim))
        w[i]=1
        w[j]=-1
        return w

def S_elem(A,B,K,w=None):
    dim=A.shape[0]
    coul=0
    D=A+B
    R=np.linalg.inv(D)
    M0=(np.pi**(dim)/np.linalg.det(D))**(3.0/2)
    trace=np.trace(B@K@A@R)
    if w!=None:
        for k in range(len(w)): ##Included to make the anion test run properly. Basically, script now handles taking a list of w's instead of just a single one. 
#The implementation is a bit hacky, since it only works if the list of w's has a length of two. 
            beta=1/(w[k].T@R@w[k])
            if k==2:
                coul+=2*np.sqrt(beta/np.pi)*M0 #Plus since we have two negative terms and one positive, meaning net one negative, that cancels minus.
            else:
                coul-=2*np.sqrt(beta/np.pi)*M0
        return M0,trace,coul
    else:
        return M0, trace

def transform_list(alphas):
    g_new = [np.ones((1, 1)) * alphas[i] for i in range(len(alphas))]
    return g_new

def S_wave(alphas, K, w=None):
    length = len(alphas)
    alphas = transform_list(alphas)
    kinetic = np.zeros((length, length))
    overlap = np.zeros((length, length))
    coulomb = np.zeros((length, length))
    for i in range(length):
        for j in range(length):
            if j <= i:
                A = alphas[i]
                B = alphas[j]
                M0, trace, coul = S_elem(A, B, K, w)
                R = np.linalg.inv(A + B)
                overlap[i, j] = M0
                overlap[j, i] = overlap[i, j]
                kinetic[i, j] = 6 * trace * M0
                kinetic[j, i] = kinetic[i, j]
                coulomb[i, j] = coul
                coulomb[j, i] = coulomb[i, j]
    return overlap, kinetic, coulomb

w_list=w_gen_3()
m_list=np.array([1, 1, 1])
K=np.array([[0,0,0],[0,1/2,0],[0,0,1/2]])
J,U=jacobi_transform(m_list)
K_trans=J@K@J.T
w_trans=[U.T @ w_list[i] for i in range(len(w_list))]

def energyS(bij,K,w):
    alphas=[]
    dim=len(w)
    for i in range(0,len(bij),dim):
        A=A_generate(bij[i:i+dim],w)
        alphas.append(A)
    N, kinetic, coulomb=S_wave(alphas,K,w)
    H=kinetic+coulomb
    E=eigh(H,N, eigvals_only='true')
    E0=np.amin(E)
    return E0

b1=7
E_list=[]
gaussians=[]
E_theoS=[]
bij=np.array([])
E_S=-0.527

#print("---------QUASI-RANDOM METHOD---------")

E_low=np.inf
bases=np.array([])
base_test=np.array([])
for i in range(50):
   hal=halton(i+1,15*len(w_trans))
   bij=-np.log(hal)*b1
   for j in range(0,len(hal),len(w_trans)):
       base_test=np.append(base_test,bij[j:j+len(w_trans)])
       E0=energyS(base_test,K_trans,w_trans)
       if E0<=E_low:
           E_low=E0
           base_curr=np.copy(bij[j:j+len(w_trans)])
       base_test=base_test[:-len(w_trans)]
   bases=np.append(bases,base_curr)
   base_test=np.append(base_test,base_curr)
   E_list.append(E_low)
   print(E_low)
   gaussians.append(i)
   E_theoS.append(E_S)

print("Best convergent numerical value:", E_list[-1])
print("Theoretical value:", E_S)
print("Difference:", np.abs(E_list[-1]-E_S))

plt.figure(2)
plt.plot(gaussians,E_list, marker='.')
plt.plot(gaussians,E_theoS, '--')
plt.title('S-wave convergence of Positron and two Electron System')
plt.xlabel('Number of Gaussians')
plt.ylabel('Energy [Hartree]')
plt.legend(['Numerical result', 'Theoretical value'])
plt.show()