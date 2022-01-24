import numpy as np
import matplotlib.pyplot as plt 

def f(x):
    return (x-2)**9

def f2(x):
    return x**9 - 18*x**8 + 144*x**7 - 672*x**6 + 2016*x**5 - 4032*x**4 + 5376*x**3 - 4608*x**2 + 2304*x - 512

f = np.vectorize(f)
f2 = np.vectorize(f2)

x = np.arange(1.92, 2.08, 0.001)
plt.plot(x, f(x), c='Blue', label=r'$(x-2)^9$')
plt.plot(x, f2(x), c='Red', label=r'Expanded form')
plt.legend()
plt.show()
