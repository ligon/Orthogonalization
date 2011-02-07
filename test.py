#!/usr/bin/env python

import math
import numpy
from numpy.linalg import norm
import random

import qrfact


#define a few functions we might want
def test_ortho( Q, ret=0 ):
    if ret == 0 :
        print "\nangle between columns of Q"
    max_dot = 0;
    m,n = Q.shape
    for k in range(0,n) :
        for j in range(k+1,n) :
            dot = numpy.dot( Q[:,k], Q[:,j] )
            angle = math.acos(dot/(norm(Q[:,k])*norm(Q[:,j])))*180/math.pi
            if ret == 0 :
                print "[", k, ",", j, "]:", angle, "deg; dot =", dot
            if abs(dot) > max_dot :
                max_dot = abs(dot) ;
    if ret == 1 :
        return max_dot

def generate_A((m,n),method='range'):
    """
    Generate a matrix (2-dimensional numpy array) A.
    
    The variable 'method' describes the algorithm.
    If equal to 'range', then an mxn array will be generated as a range of
    integers.  If equal to 'random', as an mxn array of random numbers.
    """

    if method=='range':
        # create matrix A with range of numbers, to test MGS
        A = numpy.arange(m*n).reshape(m,n)
        A = A + 1000
    elif method=='random':
        # create matrix A with random numbers
        A = numpy.zeros( (m,n) )
        for k in range(0,m) :
            for j in range(0,n) :
                A[k,j] = round(random.uniform(1,100))

    return A


def main(qr=qrfact.qr_mgs,(m,n)=(5,3)):
    """
    Test decomposition of a matrix (2-dimensional numpy array) A.
    """

    A=generate_A((5,3))

    if max(m,n)<10:
        print "A = \n", A

    Q,R = qr(A)[:2]  # Some QR routines return a third permutation P solving AP=QR.

    # repeat QR factorization to ensure Q is orthogonal, just testing
    # for count in range(0,0) :
    Q1,R1 = qr(Q)

    try:
        e=norm(Q1-Q)
        assert(e<0.001)
    except AssertionError:
        print "Orthogonality error in Q: %g." % e
        
    try:
        e=norm(R-numpy.dot(R1,R))
        assert(e<0.001)
    except AssertionError:
        print "Error in R: %g." % e


    # display Q, R, and confirm A is correct
    print "\nQ = \n", Q
    print "\nR = \n", R
    print "\nconfirm A = \n", numpy.dot(Q,R)

    # we also check how close Q is to orthogonal 
    test_ortho(Q)
    print "max dot =", test_ortho(Q,1)


    '''
    ############################################
    # this uses the built-in QR factorization method
    Q, R = numpy.linalg.qr(A)
    print "\nq = \n", Q
    print "\nr = \n", R
    print "\nconfirm a = \n", numpy.dot(Q,R)
    
    # we also check how close Q is to orthogonal 
    test_ortho(Q)
    '''

if __name__=='__main__':
    main()

