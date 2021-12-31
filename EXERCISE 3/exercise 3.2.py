import math
import numpy as np

N = 3


def CholeskyDecomposition(A):
    L = np.zeros((N, N), dtype=float)

    for i in range(N):
        for j in range(i + 1):
            if j == i:
                sum1 = 0
                for k in range(j):
                    sum1 += (L[j][k] ** 2)
                L[j][j] = math.sqrt(A[j][j] - sum1)
            elif i > j:
                sum2 = 0
                for k in range(j):
                    sum2 += (L[i][k] * L[j][k])
                L[i][j] = (A[i][j] - sum2) / L[j][j]
    return L


A = np.array([[4, 12, -16], [12, 37, -43], [-16, -43, 98]])

L = CholeskyDecomposition(A)

print(L)