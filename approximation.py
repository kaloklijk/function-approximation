# This is a program for guessing function expressions, the guess is done by
# inputting a sequence of sample values(with distance t) of the function, and
#then the program guess the expression of the functions using taylor's expansion
# It is forced that the value have to be >= 0
import sys
import math
import pyomo
import scipy
import matplotlib.pyplot as plt
from fractions import Fraction
import numpy as np


class differentiation:
    # A class for getting derivatives of functions
    def __init__(self, distance):
        self.distance = distance

    def first_derivatives(self, array):
        # return a first derivatives group 
        a = []

        for ct in range(len(array)-1):
            count = ct + 1
            for i in range(len(array)-count):
                k = i + count        
                #res = Fraction((array[k] - array[i])/t)
                res = (array[k]-array[i])/t
                try:
                    a[ct].append(res)
                except:
                    a.append([res])
        return a          

    def derivatives(self, array, g = 0):
        # return a full derivatives group

        if g == 0:
            g = len(array)
        a = []
        
        for i in range(g - 1):
            j = i+1
            if j == 1:
                a.append(differentiation(self.distance).first_derivatives(array)[0])
            else:
                a.append(differentiation(self.distance).first_derivatives(a[i-1])[0])
        return a


        
def funcguess(s, t, g):
    # guessing the expression of the function & it's exponential approximation
    a = differentiation(t).derivatives(s, g)
    x0 = np.linspace(0, (len(s)-1)*t, len(s))
    plt.plot(x0, s, "y")
    # guessing the expression of the function
    taylcoeff = [s[0]]    #taylor coefficients
    for i in range(len(s)-1):
        taylcoeff.append(a[i][0]/(math.factorial(i+1) ** 2))

    for j in range(len(s) ** 4):
        for i in range(len(s)-1):   # Compare the value to adjust the result
            taylcoeff.reverse()
            fe = np.poly1d(taylcoeff)
            taylcoeff.reverse()
            if i != len(s)-2 or len(taylcoeff) > len(s):
                try:
                    if s[i] != 0:
                        dis = (s[i] - fe(i*t))/s[i]
                        if dis > 0.2:
                            taylcoeff[i+1] += (abs(a[i][0]/(math.factorial(i+2) ** 3) * ((3*dis) ** 3)))
                        if dis < 1:
                            taylcoeff[i+1] -= (a[i][0]/(math.factorial(i+2) ** 3) * (abs(dis) ** 3))
                    else:
                        dis = fe(i*t)/3
                        if dis > 2:
                            taylcoeff[i+1] += (a[i][0]/(math.factorial(i+2) ** 3) * (dis ** 3))
                        if dis < 2:
                            taylcoeff[i+1] -= (a[i][0]/(math.factorial(i+2) ** 3) * (abs(dis) ** 3))                   
                except:
                    pass
            else:
                if s[i] != 0:
                    dis = (s[i] - fe(i*t))/s[i]
                    if dis > 0.2:
                        taylcoeff.append(abs(a[i][0]/(math.factorial(i+2) ** 3) * ((3*dis) ** 3)))
                    if dis < 1:
                        taylcoeff.append(-a[i][0]/(math.factorial(i+2) ** 3) * (dis ** 3))
                else:
                    dis = fe(i*t)/3
                    if dis > 2:
                        taylcoeff.append(a[i][0]/(math.factorial(i+2) ** 3) * (dis ** 3))
                    if dis < 2:
                        taylcoeff.append(-a[i][0]/(math.factorial(i+2) ** 3) * (abs(dis) ** 3))             
              
    taylcoeff.reverse()
    f = np.poly1d(taylcoeff)
    print(f)
    x1 = np.linspace(0, (len(s)-1)*t+2, len(s)*10)
    value = []
    for j in x1:
        value.append(f(j))
    plt.plot(x1, value, "m")
    # math.log() the samples
    ss = []

    for i in s:
        if i == 0:
            i = 0.1
        ss.append(math.log(i))
        
    a = differentiation(t).derivatives(ss, g)
    # guessing function's exponential approximation
    taylcoeff = [s[0]]    #taylor coefficients
    for i in range(len(s)-1):
        taylcoeff.append(a[i][0]/math.factorial(i+1))
    taylcoeff.reverse()
    f = np.poly1d(taylcoeff)
    print(f"\nThe exponential approximation of the function:\nexp(\n{f}\n)")
    value = []
    for j in x1:
        try:
            value.append(math.exp(f(j)))
        except:
            value.append(-1)
    plt.plot(x1, value, "r")
    plt.ylim(min(s)-1, max(s)*2)
    plt.show()

    



#driver code 
while True:
    try:
        plt.clf()
    except:
        pass
    try:
        t = int(input("\n\ninput the distances between values\n"))
    except:
        t = 1
    print("\ninput the value of the function:")

    try:
        s = []
        while True:
            s.append(int(input()))
    except:
        try:
            g = int(input("input number of derivatives for approximation(max as default)\n "))
        except:
            g = 0
        print(differentiation(t).derivatives(s, g))
        funcguess(s, t, g)

    

