
import numpy as np
import matplotlib.pyplot as plt
x=np.linspace(-15,10,100)
s=1
m=1
u=2
k=1
h=1
def psi(x):
    return (1/(s*((2*np.pi)**0.5))**0.5)*np.exp(-1*((x-u)/2/s)**2)

#defining x operator in position basis

x0=np.zeros((len(x),len(x)))
for i in range(len(x)):
    x0[i][i]=x[i]
    

potentialenergy=0.5*k*np.dot(x0,x0)

#defining p operator in x basis

d2=np.zeros((len(x),len(x)))

d2[0][0]=-2
d2[1][0]=1

d2[len(x)-1][len(x)-1]=-2
d2[len(x)-2][len(x)-1]=1

for i in range(1,len(x)-1):
    d2[i-1][i]=1
    d2[i][i]=-2
    d2[i+1][i]=1

print(d2)
dx=x[1]-x[0]
    
kineticenergy=(-1*(h**2)*d2)/(2*m*dx**2)

hamiltonian=kineticenergy + potentialenergy

Time=np.linspace(0,4*np.pi,50)


from scipy.linalg import expm

def time(t):
    timeevolution=expm((hamiltonian*-1j*t)/h)
    return timeevolution

def phi(x,t):
    z=np.array(np.matmul(time(t),psi(x)),dtype=complex)
    return z

plt.plot(x,np.abs(phi(x,1))**2)


for j in Time:
    plt.plot(x,abs(psi(x))**2,label="$|\Psi(X,0)|^2$")
    plt.plot(x,abs(phi(x,j))**2,label="$|\Psi(X,T)|^2$")
    plt.ylim(-0.5,1)
    plt.title('TIME EVOLUTION OF PROBABLITY DISTRIBUTION')
    plt.ylabel('PROBABLITY DISTRIBUTION')
    plt.xlabel('POSITION')
    plt.legend()
    plt.pause(0.1)
    plt.clf()
plt.close()
plt.show()

