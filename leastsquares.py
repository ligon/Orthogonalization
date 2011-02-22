import numpy
import qrfact


def leastsquares(A,b,qr=qrfact.qri_mgs_piv,alpha=0.5):
    """
    Modified least-squares algorithm, Dax section 7
    """
    

    A = numpy.array(A, dtype=float)
    m,n = A.shape
    z = numpy.zeros( n )
    a = numpy.zeros( n )
    x = numpy.zeros( n )
    b = numpy.transpose(b)[0]

    # do the QR factorization
    try:
        Q,R = qr(A)[:2]  # Some QR routines return a third permutation P solving AP=QR.
        PA = A
    except TypeError:
        Q,R,P = qr(A,alpha)[:3]  # Some QR routines return a third permutation P solving AP=QR.
        AP = numpy.dot( A, P )

    # Step 1'': orthogonalization of b against Q
    u = b
    for j in range( 0, n ) :
        # print "Qj = ", Q[:,j]
        # print "u = ", u
        # print "dot = ", numpy.dot( Q[:,j], u )
        z[j] = numpy.dot( Q[:,j], u )
        u = u - z[j] * Q[:,j]

    # Step 2'': iterative orthogonalization of u
    ul2norm = numpy.linalg.norm( u )
    ii = 0
    while True : # iterate
        for j in range( 0, n ) :
            a[j] = numpy.dot( Q[:,j], u )
            z[j] = z[j] + a[j]
            u = u - a[j] * Q[:,j]

        ii = ii + 1
        ulnorm = ul2norm
        ul2norm = numpy.linalg.norm( u )

        #print ul2norm, ulnorm
        
        if (ul2norm > alpha * ulnorm) or ul2norm == 0 :
            # print "used", ii, "orthogonalizations"
            break

    #print z
    #print R

    # Step 3'': use back substitution to solve Rx = z
    for i in range( n-1, -1, -1 ) :
        x[i] = z[i]
        for j in range( i+1, n ) :
            x[i] = x[i] - R[i,j] * x[j]
        x[i] = x[i] / R[i,i]
        #print x
    
    return x, AP
    



