import matplotlib.pyplot as plt 
import numpy as np 
import os 
import plotly.express as px 

data = np.loadtxt('dog_bw_data.dat')
plt.gray()
plt.imshow(data)
plt.savefig('dog_original.png')

singulars = [20, 40, 80, 160, 320, 640, 1280, 2560, 3355]

print('Generating compressed images')
for i, num in enumerate(singulars):
    if not os.path.isfile(f'compressed_{num}.png'):
        print(f'Reading in {i}')
        data = np.loadtxt(f'compr {i+1}')
        plt.gray()
        plt.imshow(data)
        plt.savefig(f'compressed_{num}.png')
    

print('Generating error plot')
# plt.clf()
errors = np.loadtxt('frobenius_errors.dat')
fig = px.line(x=singulars, y=errors)
fig.show()
# plt.plot(singulars, errors)
# plt.xlabel('Number of Singular Values')
# plt.xticks(singulars)
# plt.ylabel('Averaged Frobenius Error')
# plt.title('Frobenius Error vs. Number of Singular Values')
# plt.savefig('frobenius_errors.png', dpi=300)