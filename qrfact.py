import numpy

def qr_mgs( A ):
    """QR decomposition of A.
    
    Modified Gram-Schmidt algorithm, row version (Bjorck alg. 2.1)"""

    A = numpy.array(A, dtype=float)
    m,n = A.shape
    Q = numpy.zeros( (m,n) )
    R = numpy.zeros( (n,n) )

    for k in range( 0, n ) :
        R[k,k] = numpy.linalg.norm( A[:,k] )
        Q[:,k] = A[:,k] / R[k,k]

        for j in range( k+1, n ) :
            R[k,j] = numpy.dot( Q[:,k], A[:,j] )
            A[:,j] = A[:,j] - Q[:,k] * R[k,j]
    
    return Q,R


def qri_cgs( A, alpha ):
    """QR decomposition of A.

    Takes as arguments the matrix A and also a tolerance parameter alpha.

    Iterated CGS algorithm (Bjorck alg. 6.1(1); Hoffman section 3)"""
    
    A = numpy.array(A, dtype=float)
    m,n = A.shape
    Q1 = numpy.zeros( (m,n) )
    R = numpy.zeros( (n,n) )

    for k in range( 0, n ) :
        Qhat = A[:,k]
        #ii = 0
        while True : # iterate
            s = numpy.dot( numpy.transpose( Q1 ), Qhat )
            Qhat2 = Qhat - numpy.dot( Q1, s )
            R[:,k] = R[:,k] + s

            #ii = ii + 1
            Qhat2len = numpy.linalg.norm( Qhat2 )
            Qhatlen = numpy.linalg.norm( Qhat )
            if ( Qhat2len > alpha * Qhatlen ) :
                Qhat = Qhat2
                #print ii
                break
            Qhat = Qhat2

        R[k,k] = numpy.linalg.norm( Qhat )
        Q1[:,k] = Qhat / R[k,k]
    
    return Q1,R


def qri_mgs( A, alpha ): 
    """QR decomposition of A.

    Takes as arguments the matrix A and also a tolerance parameter alpha.

    Iterated MGS algorithm, column version (Bjorck alg. 6.1(2); Hoffman section 3)"""
    A = numpy.array(A, dtype=float)
    m,n = A.shape
    Q1 = numpy.zeros( (m,n) )
    R = numpy.zeros( (n,n) )

    for k in range( 0, n ) :
        Qhat = A[:,k]
        Qhat2 = Qhat
        ii = 0
        while True : # iterate
            for i in range( 0, k ) :
                s = numpy.dot( Q1[:,i], Qhat )
                Qhat2 = Qhat2 - s * Q1[:,i]
                R[i,k] = R[i,k] + s

            ii = ii + 1
            Qhat2len = numpy.linalg.norm( Qhat2 )
            Qhatlen = numpy.linalg.norm( Qhat )
            if (Qhat2len > alpha * Qhatlen) :
                Qhat = Qhat2
                print ii
                break
            Qhat = Qhat2
            
        R[k,k] = numpy.linalg.norm( Qhat )
        Q1[:,k] = Qhat / R[k,k]
    
    return Q1,R


def qrtest( A ):
    """QR decomposition of A.

    Bootleg reorthogonalization.
    Simply re-do the MGS QR factorization on Q to reorthogonalize."""
    
    A = numpy.array(A, dtype=float)
    m,n = A.shape
    Q = numpy.zeros( (m,n) )
    R1 = numpy.zeros( (m,n) )
    R = numpy.zeros( (n,n) )

    Q,R = qr_mgs(A)

    for ii in range(0,1) : # repeat how many times
        Q,R1 = qr_mgs(Q)
        R = numpy.dot(R1,R)

    return Q,R
    


def qri_mgs_piv( Q, alpha ):
    """QR decomposition of A, with column pivoting; returns Q,R,P.

    Takes an optional tolerance parameter alpha.
    
    Iterated MGS with column pivoting (Dax 1999)."""
    
    Q = numpy.array(Q, dtype=float)
    m,n = Q.shape
    R = numpy.zeros( (n,n) )
    Qnorms = numpy.zeros( n )
    piv = numpy.zeros( n )
    P = numpy.eye( n )

    for k in range( 0, n ) :
        print "column", k

        # step 0
        for j in range ( k, n ) :
            Qnorms[j] = numpy.linalg.norm( Q[:,j] )
        #print Qnorms
        j = numpy.where(Qnorms == max(Qnorms[k:n]))[0][0]
        Qnorms[k] = 0
        #print Q
        #print R
        #piv[k] = j
        if (j != k) :
            #print "switching columns", k, "and", j
            P[:, [j, k]] = P[:, [k, j]]
            Q[:, [j, k]] = Q[:, [k, j]]
            #if (k > 0) :
            #    R[0:k, [j, k]] = R[0:k, [k, j]]
            R[:, [j, k]] = R[:, [k, j]]
        #print Q
        #print R

        # step 1
        vl2norm = numpy.linalg.norm( Q[:,k] )
        ii = 0
        while True : # iterate
            for i in range( 0, k ) :
                s = numpy.dot( Q[:,i], Q[:,k] )
                Q[:,k] = Q[:,k] - s * Q[:,i]
                R[i,k] = R[i,k] + s

            ii = ii + 1
            vlnorm = vl2norm
            vl2norm = numpy.linalg.norm( Q[:,k] )
            if (vl2norm > alpha * vlnorm) :
                #print "on column", k, "used", ii, "orthogonalizations"
                break
            
        # step 2
        R[k,k] = numpy.linalg.norm( Q[:,k] )
        Q[:,k] = Q[:,k] / R[k,k]

        # step 3
        if (k == n) :
            break
        else :
            for j in range( k+1, n ) :
                R[k,j] = numpy.dot( Q[:,k], Q[:,j] )
                Q[:,j] = Q[:,j] - R[k,j] * Q[:,k]

        # step 4
        #Qhat = Q[:,k]
        #Qhat2 = Qhat
        for j in range( k+1, n ) :
            ii = 0
            vl2norm = numpy.linalg.norm( Q[:,j] )
            while True : # iterate
                s = numpy.dot( Q[:,j], Q[:,k] )
                R[k,j] = R[k,j] + s
                Q[:,j] = Q[:,j] - s * Q[:,k]
                
                ii = ii + 1
                vlnorm = vl2norm
                vl2norm = numpy.linalg.norm( Q[:,j] )
                if (vl2norm > alpha * vlnorm) :
                    #print "on column", j, "used", ii, "orthogonalizations"
                    break
            
    return Q,R,P






