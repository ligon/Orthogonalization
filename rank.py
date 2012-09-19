import numpy
import qrfact


def rank(A,eps=1e-15,alpha=0.5):

    m,n = A.shape
    R = qrfact.qri_mgs_piv(A, alpha)[1]

    normA = numpy.linalg.norm( A )
    # print normA

    rank = -1

    for i in range( 0, n ) :
        normR22 = numpy.linalg.norm( R[i:n,i:n] )
        normratio = normR22 / normA
        # print normratio
        
        # lots of questions about this limit
        if ( normratio < eps ) :
            rank = i
            break
    
    if rank == -1 :
        rank = i + 1
    
    return rank, R
    



