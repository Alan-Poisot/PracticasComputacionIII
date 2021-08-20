import numpy as np
from sympy.abc import x, y
from sympy import *
import time
import math
import matplotlib.pyplot as plt


class NumericalMethods:
    
    def newtonMethod(self, tol, maxIter, x0, funcOr):
        raiz = np.inf
        func= diff(funcOr,x)
        # Calcular derivada de la función original
        firstDer = diff(func, x)

        # Inicializar error e iteraciones
        error = np.inf
        iter = 0

        while iter < maxIter:
            # Validación de división entre 0
            if  firstDer.subs(x,x0)!=0:
                X1 = x0 - (func.subs(x, x0)/ firstDer.subs(x,x0))
                x0 = X1 
            # calcular error
            if abs(func.subs(x,X1))<=tol:
                return X1
            iter += 1
        
        return raiz


        
class NaiveMethods:

    """
        Este método efectua una búsqueda incremental sobre el intervalo [a,b] con pasos de 0.1
        Entradas: Intervalo [a,b], función
        Salida: Máximo y mínimo
    """
    def incrementalSearch(self, a, b, func):
        #X = np.linspace(a, b, 10000)
        #Y = np.zeros_like(X)
        i=a+0.1
        maxC=-1*np.inf
        minC=np.inf
        num1=func.subs(x, i-0.1)
        num2=func.subs(x, i)
        num3=func.subs(x, i+0.1)
        while i < b: #Se comparan los valores adyacentes para determinar la existenia de un mínimo o máximo
            if  num1>=num2 and num3>=num2 and num2 < minC:
                minC=num2
            if num1<=num2 and num3<=num2 and num2 > maxC:
                maxC=num2
            i+=0.1
            num1=num2
            num2=num3
            num3=func.subs(x, i+0.1)
        return minC,maxC
            
class Graph:

    def graph2DPlot(self, start, end, equation):
        X = np.arange(start, end,step=10)
        Y = np.zeros(len(X))
        for i in range(len(X)):
            Y[i] = equation.subs(x, X[i])
        plt.plot(X,Y)
        plt.grid()
        plt.savefig("graph",dpi=250)
     

def main():

    maxIter = 1000000
    tol = 0.0001
    func = x**3-670*x**2 + 790*x + 3000
    objG = Graph()
    #Se tuvieron que usar límites de -1000 y 1000 debido a que con 10000 no se podían apreciar correctamente los mínimos y máximos
    linf=-1000
    lsup=1000
    x0=np.linspace(linf,lsup,10)

    objG.graph2DPlot(linf, lsup, func)

    # Crear objetos para el método ingenuo
    objN = NaiveMethods()

    start = time.time()
    mni,mxi=objN.incrementalSearch(linf, lsup, func)
    end = time.time()
    print("Tiempo Búsqueda Exhaustiva:", (end-start))

    # Objeto para métodos numéticos
    objNM = NumericalMethods()
    
    
    # Método de Newton-Raphson
    criticalP=list()
    start = time.time()
    for i in x0:
        criticalP.append(func.subs(x,objNM.newtonMethod(tol, maxIter, i, func)))
    end = time.time()
    print("Tiempo Newton-Raphson:", end - start)
    mn=min(criticalP)
    mx=max(criticalP)

    print("Mínimo y máximo encontrados por el método ingenuo: ",float(mni),float(mxi))
    print("Mínimo y máximo encontrados por el método de Newton-Raphson: ",float(mn),float( mx))
    


if __name__ == "__main__":
    main()
