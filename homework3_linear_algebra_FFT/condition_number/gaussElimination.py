"""

    This module calculates a linear system by Gaussian elimination with pivoting.
    Almost a copy of on Mike Zingale's code, spring 2013.

"""

import numpy as npy

def gaussElim(A, b):
    """ perform gaussian elimination with pivoting, solving A x = b A
        is an NxN matrix, x and b are an N-element vectors.  Note: A
        and b are changed upon exit to be in upper triangular (row
        echelon) form """

    # b is a vector
    if not b.ndim == 1:
        print "ERROR: b should be a vector"
        return None

    N = len(b)

    # A is square, with each dimension of length N
    if not (A.shape[0] == N and A.shape[1] == N):
        print "ERROR: A should be square with each dim of same length as b"
        return None

    # allocation the solution array
    x = npy.zeros((N), dtype=A.dtype)

    # find the scale factors for each row -- this is used when pivoting
    scales = npy.max(npy.abs(A), 1)

    # keep track of the number of times we swapped rows
    numRowSwap = 0

    # main loop over rows
    for k in range(N):
        
        # find the pivot row based on the size of column k -- only consider
        # the rows beyond the current row
        rowMax = npy.argmax(A[k:, k]/scales[k:]) 
        if (k > 0): rowMax += k  # we sliced A from k:, correct for total rows

        # swap the row with the largest scaled element in the current column
        # with the current row (pivot) -- do this with b too!
        if not rowMax == k:
            A[[k, rowMax],:] = A[[rowMax, k],:]
            b[[k, rowMax]] = b[[rowMax, k]]
            numRowSwap += 1

        # do the forward-elimination for all rows below the current
        for i in range(k+1, N):
            coeff = A[i,k]/A[k,k]

            for j in range(k+1, N):
                A[i,j] += -A[k,j]*coeff

            A[i,k] = 0.0
            b[i] += -coeff*b[k]
       
    
    # last solution is easy
    x[N-1] = b[N-1]/A[N-1,N-1]

    for i in reversed(range(N-1)):
        isum = b[i]
        for j in range(i+1,N):
            isum += -A[i,j]*x[j]
        x[i] = isum/A[i,i]


    return x
    


