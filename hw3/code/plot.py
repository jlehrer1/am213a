import matplotlib.pyplot as plt 
import numpy as np

f = np.loadtxt('least_squares_data.dat')
x = [x[0] for x in f]
y = [x[1] for x in f]
xs = np.linspace(0, 1, 100)
three_a = np.loadtxt('vand1.txt')
five_a = np.loadtxt('vand2.txt')

def inter_1(x):
    y = three_a[0]
    for exp, a in enumerate(three_a[1:]):
        y += a*(x**(exp + 1))
    return y

def inter_2(x):
    y = five_a[0]
    for exp, a in enumerate(five_a[1:]):
        y += a*(x**(exp + 1))
    return y

inter_1 = np.vectorize(inter_1)
inter_2 = np.vectorize(inter_2)

plt.scatter(x, y)
plt.plot(x, inter_1(x))
plt.title('3rd Degree Polynomial Fitting via QR Decomposition')
plt.savefig('3d_degree.png')
plt.clf()

plt.scatter(x, y)
plt.plot(xs, inter_2(xs))

plt.title('5th Degree Polynomial Fitting via QR Decomposition')
plt.savefig('5th-degree.png')