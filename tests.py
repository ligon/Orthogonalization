#!/usr/bin/env python

import numpy
import qrfact
import unittest
import math
from numpy.linalg import norm
from numpy import matrix, eye

class TestQRDecompositions(unittest.TestCase):

    tol = 0.001
    
    def setUp(self,(m,n)=(3,5)):
        self.A = numpy.arange(m*n).reshape(m,n)
        self.A += 1000
        self.result=qrfact.qr_mgs(self.A)

    def test_A_construction(self):
        Q=self.result[0]
        R=self.result[1]
        Ahat=matrix(Q)*matrix(R)
        self.assertTrue(norm(self.A-Ahat)<0.00001)
        
    def test_orthogonality(self):
            Q,R=self.result[0:2]
            Q=matrix(Q)
            self.assertAlmostEqual(norm(Q.T*Q-eye(Q.shape[1])),0.)

if __name__ == '__main__':
    unittest.main()
