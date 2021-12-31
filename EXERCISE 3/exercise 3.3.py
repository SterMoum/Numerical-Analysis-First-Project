import numpy as np
import math


def GaussSeidel(A, b, x, N):
    for k in range(N):
        sum1 = 0
        sum2 = 0
        for j in range(k):
            sum1 += A[k][j] * x[j]
        for j in range(k + 1, N):
            sum2 += A[k][j] * x[j]
        x[k] = (b[k] - sum1 - sum2) / A[k][k]
    return x


def start(size):
    array = np.zeros((size, size), dtype=int)
    array[0][1] = -2
    array[size - 1][size - 2] = -2
    for i in range(size):
        array[i][i] = 5
        if i == size - 1:
            break
        array[i + 1, i] = -2
        array[i, i + 1] = -2

    x = np.zeros(size, dtype=float)
    vector = np.empty(size)
    for i in range(size):
        if i == 0 or i == size - 1:
            vector[i] = 3
        else:
            vector[i] = 1
    vector = np.transpose(vector)
    k = 0
    while True:
        x_previous = np.copy(x)
        x = GaussSeidel(array, vector, x, 10)
        if abs(max(abs(x)) - max(abs(x_previous))) < 10 ** (-4):
            break
        k += 1

    return x, k


result, repeats = start(10)
print("For N = 10")
print("REPEATS:", repeats)
print("Solution : ", result)

result, repeats = start(10000)
print("For N = 10000")
print("REPEATS:", repeats)
print("Solution : ", result)
