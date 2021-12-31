import numpy as np
import sys
import math

N = 3
P = np.zeros((N, N), dtype=int)
for m in range(N):
    P[m][m] = 1
L = np.zeros((N, N), dtype=float)
for m in range(N):
    L[m][m] = 1


def arrayInitialize():
    main_array = np.array([[2, 1, 5], [4, 4, -4], [1, 3, 1]])
    return main_array


def vectorInitialize():
    b = np.array([[5, 3, 2]])
    b = np.transpose(b)
    return b


def Pfactorize(A, rank):
    # rank must be smaller than N to manage the swap between the rows if necessary
    if rank >= N:
        return A
    # we need to find the maximum value of the element dependent on rank
    max = abs(A[rank][rank])
    pos = 0
    l_counter = 0
    for i in range(rank, N):
        if abs(A[i][rank]) > max:
            max = abs(A[i][rank])
            pos = i
            l_counter += 1
    # if changes happened above, make the necessary swaps
    if pos != 0:
        A[[rank, pos]] = A[[pos, rank]]
        P[[rank, pos]] = P[[pos, rank]]
        if 1 <= rank <= N:
            temp = L[pos, rank - 1]
            L[pos, rank - 1] = L[l_counter, rank - 1]
            L[l_counter, rank - 1] = temp

    return A


# Just a function that implements the multiplication between two np.arrays
def Multiply(A, B):
    result = np.zeros((N, N))

    for i in range(N):
        for j in range(len(B[0])):
            for k in range(len(B)):
                result[i][j] += A[i][k] * B[k][j]

    return result


def forwardSubstitution(b):
    x = np.empty(len(b))

    for v in range(len(b)):
        if L[v][v] == 0:
            # the value on the diagonal is zero
            x[v] = 0
            # the current variable's value is irrelevant
            continue
        # calculate v-th variable value
        value = b[v]
        # subtract linear combination of the already known variables
        # in the bottom left corner of the matrix
        for i in range(v):
            value = value - (L[v][i] * x[i])
        # divide by the coefficient by the v-th variable to get it's value
        value = value / L[v][v]
        # store the value in the resulting vector
        x[v] = value

    return x


def backwardSubstitution(U, b):
    x = np.empty(len(b))

    for v in range(len(b) - 1, -1, -1):
        if U[v][v] == 0:
            # the value on the diagonal is zero
            x[v] = 0
            continue
        # calculate v-th variable value
        value = b[v]
        # subtract linear combination of the already known variables
        # in the top right corner of the matrix
        for i in range(v + 1, m, 1):
            value -= U[v][i] * x[i]
        # divide by the coefficient before the i-th variable to get it's value
        value /= U[v][v]
        # store the value in the resulting vector
        x[v] = value

    return x


def SolveSystem(U):
    b = np.matmul(P, vectorInitialize())
    y = forwardSubstitution(b)
    x = backwardSubstitution(U, y)
    return x


# implements Gauss elimination
def GaussElimination():
    rank = 0
    ratio_list = []
    U = arrayInitialize()
    A = arrayInitialize()
    for i in range(N - 1):
        U = Pfactorize(U, rank)
        rank += 1
        if U[i][i] == 0:
            sys.exit('Divide by zero')
        for j in range(i + 1, N):
            ratio = U[j][i] / U[i][i]
            for k in range(N):
                U[j][k] = U[j][k] - ratio * U[i][k]
            L[j][i] = ratio
    A = Multiply(P, arrayInitialize())
    B = Multiply(L, U)
    if np.array_equal(A, B):
        print("PA = LU analysis was succesfull")
    X = SolveSystem(U)
    for i in range(len(X)):
        print('x' + str(i), ' = ', X[i])
    print("P = ", P)
    print("A = ", arrayInitialize())
    print("L = ", L)
    print("U = ", U)


def start():
    array = arrayInitialize()
    if np.linalg.det(array) != 0:
        GaussElimination()


# main:
start()
