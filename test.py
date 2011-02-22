#!/usr/bin/env python

import sys
import math
import numpy
from numpy.linalg import norm
import random

import qrfact
import leastsquares
import rank


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
        A = A + 10
    elif method=='random':
        # create matrix A with random numbers
        A = numpy.zeros( (m,n) ) 
	for k in range(0,m) :
            for j in range(0,n) :
                A[k,j] = round(random.uniform(1,100))

    return A


def main(qr=qrfact.qr_mgs,(m,n)=(5,3),alpha=0.5):
    """
    Test decomposition of a matrix (2-dimensional numpy array) A.
    """

    A=generate_A((5,3))

    if max(m,n)<10:
        print "A = \n", A


    try:
        Q,R = qr(A)[:2]  # Some QR routines return a third permutation P solving AP=QR.
    except TypeError:
        Q,R = qr(A,alpha)[:2]  # Some QR routines return a third permutation P solving AP=QR.


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

def test_ls( qr=qrfact.qri_mgs_piv, (m,n)=(5,5), alpha=0.5 ) :
    """
    Test least squares routine
    """

    A = generate_A( (m,n), 'random' )
    b = generate_A( (m,1), 'random' )

    if max(m,n) < 10:
        print "A = \n", A
        print "b = \n", b

    x,AP = leastsquares.leastsquares( A, b, qr=qr, alpha=alpha )
    
    print "AP = \n", AP
    print "x = ", x
    print "APx = \n", numpy.dot( AP, x )


def test_rank((m,n)=(4,4), alpha=0.5 ) :
    """
    Test rank finding routine
    """

    A = generate_A( (m,n) )
    #err = generate_A( (1,m), 'random' )
    
    #A[:,n-1] = A[:,0] + 0.1 * A[:,1]
    #A[:,n-1] = A[:,0] + 0.000001 * err
    print "A = \n", A

    therank,R = rank.rank(A)
    
    print "R = \n", R
    print "rank = ", therank
    
   

if __name__=='__main__':
    arg = sys.argv[1:]

    if (len(arg) > 0) and (arg[0] == 'ls') :
        print "doing leastsquares test"
        test_ls()
    if (len(arg) > 0) and (arg[0] == 'rank') :
        print "doing rank test"
        test_rank()
    else :
        print "doing main loop"
        main(qr=qrfact.qri_mgs_piv)


