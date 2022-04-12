import numpy as np
if __name__=="__main__":

#    def EtaB(self):
#        # intial conditions in the order RN11, RN12, RN21, RN22, RNb11, RNb12, RNb21, RNb22, mudelta1, mudelta2,  mudelta3
#        y0       = np.array([1+0j,0+0j, 0+0j, 1+0j, 1+0j, 0+0j, 0+0j, 1+0j, 0+0j, 0+0j, 0+0j], dtype=np.complex128)

#        params   = np.array([epstt,epsmm,epsee,k], dtype=np.complex128)
        
        # CI parameters for test
        th12val  = np.arcsin(0.557)
        th13val  = np.arcsin(0.1497)
        th23val  = np.arcsin(0.75)
        omegaval = np.pi/4 - 0.7 * 1j
        m1v1l    = 0.
        m2val    = 8.6e-12
        m3val    = 58e-12
        M1val    = 40.0
        dMval    = 0.2e-10
        M2val    = M1val + dMval
        deltaval = 221. * np.pi/180.
        alpha1val= np.pi/3 - deltaval
        vev      = 246.
        

        
        # construct the CI parametrsition for internal tests
        
        mnu      = matrix_diag3(m1val, m2val, m3val)
        mM       = matrix_diag2(M1val, M2val)
        R_mat    = np.matrix([ [0,0] , [np.cos(omegaval), np.sin(omegaval)] , [-np.sin(omegaval), np.cos(omegaval)] ])
        
        
        Fmat     = (np.sqrt(2)/ vev) * matrix_pmns(th12val, th13val, th23val, deltaval, alpha1val).conjugate() @ np.sqrt( mnu ) @ R_mat.conjugate() @ np.sqrt(mM)
        
        print(Fmat)
       

       

