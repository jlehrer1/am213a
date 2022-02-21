import numpy as np 
from numpy import array, zeros, fabs, linalg

A = np.array(
    [[2,1,1,0],
    [4,3,3,1],
    [8,7,9,5],
    [6,7,9,8]]
).astype(float)

B = np.array(
    [[3, 0, 1, 0.9, 2.1, 3.141592653589793],
    [6.0, -2.0, 1.0, 10.4, -491.2, -4.712388980384690],
    [10.0, 2.0,	0.0, -20.2, 0.12, 2.243994752564138],
    [1.0, 10.0, -5.0, -5.12, -51.3, 2.356194490192345]]
)

print(linalg.solve(A, B))

# a = A.copy()
# b = B.copy()[:, 0]

# n = len(b)
# x = zeros(n, float)

# #first loop specifys the fixed row
# for k in range(n-1):
#     if fabs(a[k,k]) < 1.0e-12:
        
#         for i in range(k+1, n):
#             if fabs(a[i,k]) > fabs(a[k,k]):
#                 a[[k,i]] = a[[i,k]]
#                 b[[k,i]] = b[[i,k]]
#                 break

#  #applies the elimination below the fixed row

#     for i in range(k+1,n):
#         if a[i,k] == 0:continue

#         factor = a[k,k]/a[i,k]
#         for j in range(k,n):
#             a[i,j] = a[k,j] - a[i,j]*factor
#             #we also calculate the b vector of each row
#         b[i] = b[k] - b[i]*factor
# print(a)
# print(b)

