import numpy as np
import matplotlib.pyplot as plt

x=np.linspace(-0.5,0.5,300)
dx=x[1]-x[0]
l=1
h=1
m=1


def eigen(x,n):
    return ((2/l)**0.5)*np.sin(((n*np.pi*x)/l)+(n*np.pi/2))

def psi(x):
    return (1/(6**0.5))*eigen(x,3)   + (2/(6**0.5))*eigen(x,4)     +(1/(6**0.5))*eigen(x,6)                        
    
def c(n):
   return (sum(psi(x)*eigen(x,n)))*dx

def energy(n):
    return ((n**2)*(h**2))/(8*m*(l**2))

def time(t,n):
    return np.exp((-1j*energy(n)*t*2*np.pi)/h)

def phi(x,t):
    z=np.zeros(len(x),dtype='complex')
    for i in range(1,100):
        z+=(c(i)*eigen(x,i)*time(t,i))
    return z




Time=np.linspace(0,10,len(x))


for j in Time:
    plt.plot(x,abs(psi(x))**2,label="$|\Psi(X,0)|^2$")
    plt.plot(x,abs(phi(x,j))**2,label="$|\Psi(X,T)|^2$")
    plt.ylim(-0.5,10)
    plt.title('TIME EVOLUTION OF PROBABLITY DISTRIBUTION')
    plt.ylabel('PROBABLITY DISTRIBUTION')
    plt.xlabel('POSITION')
    plt.legend()
    plt.pause(0.01)
    plt.clf()
plt.close()


