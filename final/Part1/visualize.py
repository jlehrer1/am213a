import matplotlib.pyplot as plt 
import numpy as np 

data = np.loadtxt('dog_bw_data.dat')
plt.gray()
plt.imshow(data)
plt.show()