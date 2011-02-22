import numpy
import qrfact


def rank(A,eps=(10^-15),alpha=0.5):

    m,n = A.shape
    Q,R,P = qrfact.qri_mgs_piv(A, alpha)

    normA = numpy.linalg.norm( A )
    # print normA
    for i in range( 0, n ) :
        normR22 = numpy.linalg.norm( R[i:n,i:n] )
        normratio = normR22 / normA
        # print normratio
        
        # lots of questions about this limit
        if ( normratio < eps ) :
            break
    
    rank = i
    return rank, R
    



