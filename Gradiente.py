from sympy.core import symbol
import numpy as np
from sympy.abc import x, y
from sympy import *
import time
import matplotlib.pyplot as plt


class NumericalMethods:
    
    def newtonMethod(self, tol, maxIter, x0, funcOr):
        raiz = np.inf
        func= diff(funcOr,h)
        # Calcular derivada de la función original
        firstDer = diff(func, h)

        # Inicializar error e iteraciones
        error = np.inf
        iter = 0

        while iter < maxIter:
            # Validación de división entre 0
            if  firstDer.subs(h,x0)!=0:
                X1 = x0 - (func.subs(h, x0)/ firstDer.subs(h,x0))
                x0 = X1 
            # calcular error
            if abs(func.subs(h,X1))<=tol:
                return X1
            iter += 1
            func.simplify()
            firstDer.simplify
        
        return raiz

def grad(f,n,X): #Cálculo de gradiente
    g=list()
    for i in range(n):
        g.append(-1*diff(f,X[i]))
    return g

def maxGrad(f,g,n,X,v): #Determinación de la distancia óptima en la dirección del gradiente
    tol=0.000001
    maxIter=1000000
    h=symbols('h')
    k=list()
    for i in range(n):
        k.append(v[i]+h*g[i]) #Se crean las nuevas expresiones para las variables a evaluar
    objNM = NumericalMethods()
    return objNM.newtonMethod(tol, maxIter, 0.1, f.subs(dict(zip(X,k))))


def hess(f,n,X):
    M=Matrix.zeros(n, n)
    #print(M)
    for i in range(n):
        #M.append(list())
        for j in range(n):
            w=diff(f,X[i])
            M[i,j]=diff(w,X[j])
    return M


n=1 #Número de variables
X=symbols('x0:4') #Declaración de variables simbólicas
h=symbols('h')
if n==1: #Definición de las funciones a maximizar
    f=-(X[0]-10)**4
elif n==2:
    f=-(X[0]-10)**4-5*(X[1]+6)**6-3
elif n==3:
    f=-(X[0]-10)**4-5*(X[1]+6)**6-3-5*(3*X[2]+2)**8
else:
    f=-(X[0]-10)**4-5*(X[1]+6)**6-3-5*(3*X[2]+2)**8-2*(X[3]-5)**2

v=[0.1,0.1,0.1,0.1] #Valores iniciales

m=hess(f,n,X)

p=1
gr=grad(f,n,X)
dg=[10,10,10,10]
mag=1

start=time.time()
while sqrt(mag)>0.0001: #Se usa la longitud del último "paso" dado como criterio de cercanía a la solución
    mag=0
    for i in range(n):
        dg[i]=gr[i].subs(dict(zip(X,v)))
        
    p=maxGrad(f,dg,n,X,v)
    for i in range(n):
        v[i]+=dg[i]*p
        mag+=(dg[i]*p)**2

end=time.time()

print("El máximo se enuentra en") #Impresión de resultados
for i in range(n):
    print(X[i],"=",v[i])
print("El tiempo de ejecución fue de",end-start,"s")

if n==1: #Generación de gráfica para el caso de una variable
    Xg = np.arange(v[0]-0.1, v[0]+0.1,step=0.01)
    Yg = np.zeros(len(Xg))
    for i in range(len(Xg)):
        Yg[i] = f.subs(X[0],Xg[i])
    plt.plot(Xg,Yg)
    plt.grid()
    plt.savefig("graph",dpi=250)
