"""

    This module finds the inverse of a matrix using Gauss method.
    Almost a copy of on Mike Zingale's code, spring 2013.

"""


import numpy


def inverseGauss(AInput):
    """ return the inverse of AInput """

    A = AInput.copy()
    N = len(A[:,0])

    # A is square, with each dimension of length N
    if not (A.shape[0] == A.shape[1]):
        print "ERROR: A should be square"
        return None

    # create an identity matrix
    I = numpy.identity(N)

    # allocation for the inverse
    Ainv = numpy.zeros((N,N), dtype=A.dtype)

    # find the scale factors for each row -- this is used when pivoting
    scales = numpy.max(numpy.abs(A), 1)

    # keep track of the number of times we swapped rows
    numRowSwap = 0

    # we are essentially doing Gaussian elimination, but with A Ainv = I
    # each column of I represents a separate righthand side to an Ax = b
    # linear system

    # main loop over rows
    for k in range(N):
        
        # find the pivot row based on the size of column k -- only consider
        # the rows beyond the current row
        rowMax = numpy.argmax(A[k:, k]/scales[k:]) 
        if (k > 0): rowMax += k  # we sliced A from k:, correct for total rows

        # swap the row with the largest scaled element in the current column
        # with the current row (pivot) -- do this with b too!
        if not rowMax == k:
            A[[k, rowMax],:] = A[[rowMax, k],:]
            I[[k, rowMax],:] = I[[rowMax, k],:]
            numRowSwap += 1

        # do the forward-elimination for all rows below the current
        for i in range(k+1, N):
            coeff = A[i,k]/A[k,k]

            for j in range(k+1, N):
                A[i,j] += -A[k,j]*coeff

            A[i,k] = 0.0
            I[i,:] += -coeff*I[k,:]
    
    

    # back-substitution -- once for each column in the I matrix
    
    for c in range(N):

        # last solution is easy
        Ainv[N-1,c] = I[N-1,c]/A[N-1,N-1]

        for i in reversed(range(N-1)):
            isum = I[i,c]
            for j in range(i+1,N):
                isum += -A[i,j]*Ainv[j,c]
            Ainv[i,c] = isum/A[i,i]


    # determinant
    #idet = numpy.prod(numpy.diagonal(A))*(-1.0)**numRowSwap
    
    return Ainv
