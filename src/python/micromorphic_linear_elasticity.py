import numpy as np
import hex8
import finite_difference as fd
import unittest
import os

from hex8 import T_to_V_mapping as T2V
from hex8 import V_to_T_mapping as V2T

"""
==============================================
|                                            |
|       Micromorphic Linear Elasticity       |
|                                            |
==============================================
|                                            |
| A constitutive model for micromorphic      |
| linear elasticity. Written in a style      |
| which will be simple to port to Abaqus     |
| for use in finite element simulations.     |
|                                            |
| Model derived from a quadratic form of     |
| the Helmholtz free energy detailed in      |
| "Nonlinear Micromorphic Continuum          |
| Mechanics and Finite Strain                |
| Elastoplasticity" by Dr. Richard Regueiro  |
| with expansions by Nathan Miller           |
==============================================
"""

#Compute strain measures

def compute_strain_measures(C,Psi): #Test function written
    """Compute the strain measures"""
    Iten = np.eye(3)
    E_macro = np.zeros([9,])
    E_micro = np.zeros([9,])
    for index in range(len(E_macro)):
        I,J = V2T(index,[3,3])
        E_macro[index] = 0.5*(  C[index]-Iten[I,J])
        E_micro[index] =      Psi[index]-Iten[I,J]
    return E_macro,E_micro
    
def compute_dCinvdC(Cinv): #Test function written
    """Compute the derivative of Cinv w.r.t. C"""
    dCinvdC = np.zeros([3*3*3*3])
    for I in range(3):
        for J in range(3):
            for K in range(3):
                for L in range(3):
                    IJKL = T2V([I,J,K,L],[3,3,3,3])
                    IK   = T2V([I,K],[3,3])
                    LJ   = T2V([L,J],[3,3])
                    dCinvdC[IJKL] = -Cinv[IK]*Cinv[LJ]
    return dCinvdC

###### Form Stiffness Tensors ######

def form_stiffness_tensors(PARAMS): #Test function written
    """Form the stiffness tensors"""
    LAMBDA,MU,ETA,TAU,KAPPA,NU,SIGMA,TAUS = PARAMS
    AFOT = form_A(LAMBDA,MU)
    BFOT = form_B(ETA,TAU,KAPPA,NU,SIGMA)
    CFOT = form_C(TAUS)
    DFOT = form_D(TAU,SIGMA)
    return AFOT,BFOT,CFOT,DFOT
    
def form_A(LAMBDA,MU): #Test function written
    """Form the A stiffness tensor"""
    A = np.zeros([3*3*3*3])
    I = np.eye(3)
    for k in range(3):
        for l in range(3):
            for m in range(3):
                for n in range(3):
                    KLMN = T2V([k,l,m,n],[3,3,3,3])
                    A[KLMN] = LAMBDA*I[k,l]*I[m,n] + MU*(I[k,m]*I[l,n]+I[k,n]*I[l,m])
    return A
    
def form_B(ETA,TAU,KAPPA,NU,SIGMA): #Test function written
    """Form the B stiffness tensor"""
    B = np.zeros([3*3*3*3])
    I = np.eye(3)
    for k in range(3):
        for l in range(3):
            for m in range(3):
                for n in range(3):
                    KLMN = T2V([k,l,m,n],[3,3,3,3])
                    B[KLMN] = (ETA-TAU)*I[k,l]*I[m,n]+KAPPA*I[k,m]*I[l,n]+NU*I[k,n]*I[l,m]\
                              -SIGMA*(I[k,m]*I[l,n]+I[k,n]*I[l,m])
    return B
                                 
def form_C(TAUS): #Test function written
    """Form the C stiffness tensor"""
    TAU1,TAU2,TAU3,TAU4,TAU5,TAU6,TAU7,TAU8,TAU9,TAU10,TAU11 = TAUS
    C = np.zeros([3*3*3*3*3*3])
    I = np.eye(3)
    for k in range(3):
        for l in range(3):
            for m in range(3):
                for n in range(3):
                    for p in range(3):
                        for q in range(3):
                            KLMNPQ = T2V([k,l,m,n,p,q],[3,3,3,3,3,3])
                            C[KLMNPQ] = TAU1*(I[k,l]*I[m,n]*I[p,q]+I[k,q]*I[l,m]*I[n,p])\
                                       +TAU2*(I[k,l]*I[m,p]*I[n,q]+I[k,m]*I[l,q]*I[n,p])\
                                       +TAU3*I[k,l]*I[m,q]*I[n,p]+TAU4*I[k,n]*I[l,m]*I[p,q]\
                                       +TAU5*(I[k,m]*I[l,n]*I[p,q]+I[k,p]*I[l,m]*I[n,q])\
                                       +TAU6*I[k,m]*I[l,p]*I[n,q]+TAU7*I[k,n]*I[l,p]*I[m,q]\
                                       +TAU8*(I[k,p]*I[l,q]*I[m,n]+I[k,q]*I[l,n]*I[m,p])+TAU9*I[k,n]*I[l,q]*I[m,p]\
                                       +TAU10*I[k,p]*I[l,n]*I[m,q]+TAU11*I[k,q]*I[l,p]*I[m,n]
    return C

def form_D(TAU,SIGMA): #Test function written
    """Form the D stiffness Matrix"""
    D = np.zeros([3*3*3*3])
    I = np.eye(3)
    for k in range(3):
        for l in range(3):
            for m in range(3):
                for n in range(3):
                    KLMN = T2V([k,l,m,n],[3,3,3,3])
                    D[KLMN] = TAU*I[k,l]*I[m,n]+SIGMA*(I[k,m]*I[l,n]+I[k,n]*I[l,m])
    return D
    
###### Compute Stresses ######
    
def compute_pk2_stress(AFOT,BFOT,CFOT,DFOT,E_macro,E_micro,C,Gamma): #Test function written
    """Compute the Second Piola Kirchhoff Stress returns PK2 and other useful terms"""
    PK2    = np.zeros([9,])
    I      = np.array([1,0,0,0,1,0,0,0,1]).astype(float)
    Cinv   = hex8.invert_3x3_matrix_V(C)[1]
    STRAIN_TERM = E_micro+I#hex8.matrix_dot_V(Cinv,E_micro)+I
    
    TERM1  = np.zeros([9,])
    TERM2  = np.zeros([9,])
    
    for i in range(3):
        for j in range(3):
            IJ = T2V([i,j],[3,3])
            for k in range(3):
                for l in range(3):
                    IJKL = T2V([i,j,k,l],[3,3,3,3])
                    KL   = T2V([k,l],[3,3])
                    TERM1[IJ] += AFOT[IJKL]*E_macro[KL] + DFOT[IJKL]*E_micro[KL]
    
    for i in range(3):
        for j in range(3):
            IJ = T2V([i,j],[3,3])
            for k in range(3):
                for l in range(3):
                    KL   = T2V([k,l],[3,3])
                    for q in range(3):
                        for r in range(3):
                            IQKL = T2V([i,q,k,l],[3,3,3,3])
                            KL   = T2V([k,l],[3,3])
                            RQ   = T2V([r,q],[3,3])
                            JR   = T2V([j,r],[3,3])
                            TERM2[IJ] += (BFOT[IQKL]*E_micro[KL]+DFOT[IQKL]*E_macro[KL])*STRAIN_TERM[RQ]*Cinv[JR]
    
    for i in range(3):
        for j in range(3):
            IJ = T2V([i,j],[3,3])
            for l in range(3):
                for m in range(3):
                    for n in range(3):
                        LMN = T2V([l,m,n],[3,3,3])
                        for q in range(3):
                            for r in range(3):
                                for s in range(3):
                                    IQRLMN = T2V([i,q,r,l,m,n],[3,3,3,3,3,3])
                                    JS     = T2V([j,s],[3,3])
                                    SQR    = T2V([s,q,r],[3,3,3])
                                    TERM2[IJ] += CFOT[IQRLMN]*Gamma[LMN]*Cinv[JS]*Gamma[SQR]
    
    for i in range(3):
        for j in range(3):
            IJ = T2V([i,j],[3,3])
            PK2[IJ] = TERM1[IJ]+TERM2[IJ]
    return PK2,TERM1,TERM2
    
def compute_symmetric_stress(TERM1,TERM2): #Test function written
    """Compute the symmetric stress from the given terms"""
    SIGMA = np.zeros([3*3,])
    TERM2_SYMM  = hex8.get_symm_matrix_V(TERM2)
    for i in range(3):
        for j in range(3):
            IJ = T2V([i,j],[3,3])
            SIGMA[IJ] = TERM1[IJ]+2*TERM2_SYMM[IJ]
    return SIGMA
    
def compute_ho_stress(CFOT,Gamma):
    """Get the higher order stress"""
    M = np.zeros([27,])
    for k in range(3):
        for l in range(3):
            for m in range(3):
                KLM = T2V([k,l,m],[3,3,3])
                for n in range(3):
                    for p in range(3):
                        for q in range(3):
                            KLMNPQ = T2V([k,l,m,n,p,q],[3,3,3,3,3,3])
                            NPQ    = T2V([n,p,q],[3,3,3])
                            M[KLM] += CFOT[KLMNPQ]*Gamma[NPQ]
    return M
    
###### Compute the tangents ######

def compute_tangents(AFOT,BFOT,CFOT,DFOT,E_macro,E_micro,C,Gamma):
    """Compute the tangents for a micromorphic linear elastic
    constitutive model"""
    
    #Common Terms
    I    = np.eye(3)
    Cinv = hex8.invert_3x3_matrix(C)
    
    #Compute dCinvdC
    dCinvdC = compute_dCinvdC(Cinv)
    
    #Compute tangents
    dSdC,    dSigmadC,    dMdC     = compute_stress_derivatives_wrt_C(E_macro,E_micro,AFOT,BFOT,CFOT,DFOT,Cinv,I,dCinvdC)
    dSdPhi,  dSigmadPhi,  dMdPhi   = compute_stress_derivatives_wrt_Phi(E_macro,E_micro,AFOT,BFOT,CFOT,DFOT,Cinv,I)
    dSdGamma,dSigmadGamma,dMdGamma = compute_stress_derivatives_wrt_Gamma(E_macro,E_micro,AFOT,BFOT,CFOT,DFOT,Cinv,I)
    
    return dSdC,dSdPhi,dSdGamma,dSigmadC,dSigmadPhi,dSigmadGamma,dMdC,dMdPhi,dMdGamma

def compute_stress_derivatives_wrt_C(E_macro,E_micro,AFOT,BFOT,CFOT,DFOT,Cinv,I,dCinvdC):
    """Compute the stress derivatives w.r.t. C"""
    
    #Initialize tangents with respect to C
    dSdC     = np.zeros([3,3,3,3])
    dSigmadC = np.zeros([3,3,3,3])
    dMdC     = np.zeros([3,3,3,3,3])
    
    #Useful terms for computation to prevent repetition
    TERM1 = np.zeros([3,3,3,3])
    TERM2 = np.zeros([3,3,3,3])
    
    #Compute dCinvdC
    dCinvdC = compute_dCinvdC(Cinv)
    
    #Compute dSdC, dSigmadC, dMdC
    for I in range(3):
        for J in range(3):
            for O in range(3):
                for P in range(3):
                    TERM1[I,J,O,P] = AFOT[I,J,O,P]
                    
                    for Q in range(3):
                        for R in range(3):
                            TERM1[I,J,O,P] += DFOT[I,Q,O,P]*(E_micro[R,Q]+I[R,Q])*Cinv[J,R]
                    
                            for K in range(3):
                                for L in range(3):
                                    TERM2[I,J,O,P] += (BFOT[I,Q,K,L]*E_micro[K,L]+DFOT[I,Q,K,L]*E_macro[K,L])*(E_micro[R,Q]+I[R,Q])*dCinvdC[J,R,O,P]
                                    
                            for L in range(3):
                                for M in range(3):
                                    for N in range(3):
                                        TERM2[I,J,O,P] += CFOT[I,Q,R,L,M,N]*Gamma[L,M,N]*Gamma[S,Q,R]*dCinvdC[J,S,O,P]
                                        
                    dSdC[I,J,O,P]     = TERM1[I,J,O,P] + TERM2[I,J,O,P]
                    dSigmadC[I,J,O,P] = TERM1[I,J,O,P] + (TERM2[I,J,O,P]+TERM2[J,I,O,P])
    
    return dSdC,dSigmadC,dMdC
    
def compute_stress_derivatives_wrt_Phi(E_macro,E_micro,AFOT,BFOT,CFOT,DFOT,Cinv,I):
    """Compute the stress derivatives w.r.t. Phi"""
    
    #Initialize tangents with respect to Phi
    dSdPhi     = np.zeros([3,3,3,3])
    dSigmadPhi = np.zeros([3,3,3,3])
    dMdPhi     = np.zeros([3,3,3,3,3])
    
    TERM1 = np.zeros([4,4,4,4])
    TERM2 = np.zeros([4,4,4,4])
    
    #Compute dSdPhi, dSigmadPhi, and dMdPhi
    for I in range(3):
        for J in range(3):
            for O in range(3):
                for P in range(3):
                    dSdPhi[I,J,O,P]     = DFOT[I,J,O,P]
                    dSigmadPhi[I,J,O,P] = DFOT[I,J,O,P]
                    
                    for Q in range(3):
                        for R in range(3):
                            TERM2[I,J,O,P] += BFOT[I,Q,O,P]*(E_micro[R,Q]+I[R,Q])*Cinv[J,R]
                            
                            for K in range(3):
                                for L in range(3):
                                    TERM2[I,J,O,P] += (BFOT[I,P,K,L]*E_micro[K,L]+DFOT[I,P,K,L]*E_macro[K,L])*Cinv[J,O]
                                    
                    dSdPhi[I,J,O,P]     += TERM2[I,J,O,P]
                    dSigmadPhi[I,J,O,P] += TERM2[I,J,O,P] + TERM2[J,I,O,P]
                    
    return dSdPhi,dSigmadPhi,dMdPhi
    
def compute_stress_derivatives_wrt_Gamma(E_macro,E_micro,AFOT,BFOT,CFOT,DFOT,Cinv,I):
    """Compute the stress derivatives w.r.t. Gamma"""
    #Initialize tangents with respect to Gamma
    dSdGamma     = np.zeros([3,3,3,3,3])
    dSigmadGamma = np.zeros([3,3,3,3,3])
    
    #Compute dSdGamma, dSigmadGamma, and dMdGamma
    TERM = np.zeros([3,3,3,3,3])
    for I in range(3):
        for J in range(3):
            for T in range(3):
                for U in range(3):
                    for V in range(3):
                        for S in range(3):
                            for Q in range(3):
                                for R in range(3):
                                    dSdGamma[I,J,T,U,V] += CFOT[I,Q,R,T,U,]*Cinv[J,S]*Gamma[S,Q,R] + CFOT[I,U,V,S,Q,R]*Gamma[S,Q,R]*Cinv[J,T]
                        dSigmadGamma[I,J,T,U,V] = dSdGamma[I,J,T,U,V] + dSdGamma[J,I,T,U,V]
    
    return dSdGamma,dSigmadGamma,CFOT
    
    
###### Call the model ######

def micromorphic_linear_elasticity(F,chi,grad_chi,params):
    """A constitutive model for micromorphic linear elasticity"""
    C,Phi,Gamma         = get_deformation_measures(F,chi,grad_chi)
    AFOT,BFOT,CFOT,DFOT = form_stiffness_tensors(params)
    CAUCHY,TERM1,TERM2  = compute_pk2_stress(AFOT,BFOT,CFOT,DFOT,E_macro,E_micro,C,Gamma)
    SIGMA               = compute_symmetric_stress(TERM1,TERM2)
    M                   = compute_higher_order_stress(CFOT,Gamma)
    dSdC,dSdPhi,dSdGamma,dSigmadC,dSigmadPhi,dSigmadGamma,dMdC,dMdPhi,dMdGamma = compute_tangents(AFOT,BFOT,CFOT,DFOT,E_macro,E_micro,C,Gamma)
    return CAUCHY,SIGMA,M,dSdC,dSdPhi,dSdGamma,dSigmadC,dSigmadPhi,dSigmadGamma,dMdC,dMdPhi,dMdGamma
    
class TestMicro_LE(unittest.TestCase):

    f                    = None
    original_directory   = ""
    output_file_name     = r"micromorphic_linear_elasticity.txt"
    output_file_location = r".\tests\unittests"
    currentResult        = None
    @classmethod
    def setUpClass(self):
        """Setup method"""
        output_file = os.path.join(self.output_file_location,self.output_file_name)
        if(os.path.isfile(output_file)):
            os.remove(output_file)
        self.f = open(output_file,"w+")
    @classmethod
    def tearDownClass(self):
        """Teardown method"""
        self.f.close()
        
    def setUp(self):
        pass
        
    def tearDown(self):
        ok = self.currentResult.wasSuccessful()
        tname = self.id().split(".")[-1]
        self.f.write(tname+"\t&\t"+str(ok)+"\n")
        
    def run(self, result=None):
        """Redefine run to keep track of results"""
        self.currentResult = result
        unittest.TestCase.run(self,result)
        
    def test_form_A(self):
        """Test forming the A stiffness tensor"""
        LAMBDA = 2.4
        MU     = 6.7
        
        A = form_A(LAMBDA,MU)
        
        I = np.eye(3)
        At = np.zeros([3,3,3,3])
        
        for K in range(3):
            for L in range(3):
                for M in range(3):
                    for N in range(3):
                        At[K,L,M,N] = LAMBDA*I[K,L]*I[M,N] + MU*(I[K,M]*I[L,N] + I[K,N]*I[L,M])
        At = hex8.reduce_tensor_to_vector_form(At)
        
        self.assertEqual(np.allclose(A,At),True)
        
    def test_form_B(self):
        """Test forming the B stiffness tensor"""
        ETA,TAU,KAPPA,NU,SIGMA = 2.4,5.1,5.6,8.2,2.
        
        B  = form_B(ETA,TAU,KAPPA,NU,SIGMA)
        
        I  = np.eye(3)
        Bt = np.zeros([3,3,3,3])
        
        for K in range(3):
            for L in range(3):
                for M in range(3):
                    for N in range(3):
                        Bt[K,L,M,N] = (ETA-TAU)*I[K,L]*I[M,N] + KAPPA*I[K,M]*I[L,N] + NU*I[K,N]*I[L,M]\
                                      -SIGMA*(I[K,M]*I[L,N]+I[K,N]*I[L,M])
        Bt = hex8.reduce_tensor_to_vector_form(Bt)                        
        self.assertEqual(np.allclose(B,Bt),True)
        
    def test_form_C(self):
        """Test forming the C stiffness tensor"""
        TAUS = [4.5,1.3,9.2,1.1,6.4,2.4,7.11,5.5,1.5,3.8,2.7]
         
        C = form_C(TAUS)
        I = np.eye(3)
        
        Ct = np.zeros([3,3,3,3,3,3])
        
        TAU1,TAU2,TAU3,TAU4,TAU5,TAU6,TAU7,TAU8,TAU9,TAU10,TAU11 = TAUS
        
        for K in range(3):
            for L in range(3):
                for M in range(3):
                    for N in range(3):
                        for P in range(3):
                            for Q in range(3):
                                Ct[K,L,M,N,P,Q] = TAU1*(I[K,L]*I[M,N]*I[P,Q] + I[K,Q]*I[L,M]*I[N,P])\
                                                 +TAU2*(I[K,L]*I[M,P]*I[N,Q] + I[K,M]*I[L,Q]*I[N,P])\
                                                 +TAU3*I[K,L]*I[M,Q]*I[N,P]  + TAU4*I[K,N]*I[L,M]*I[P,Q]\
                                                 +TAU5*(I[K,M]*I[L,N]*I[P,Q] + I[K,P]*I[L,M]*I[N,Q])\
                                                 +TAU6*I[K,M]*I[L,P]*I[N,Q]  + TAU7*I[K,N]*I[L,P]*I[M,Q]\
                                                 +TAU8*(I[K,P]*I[L,Q]*I[M,N] + I[K,Q]*I[L,N]*I[M,P])\
                                                 +TAU9*I[K,N]*I[L,Q]*I[M,P]  + TAU10*I[K,P]*I[L,N]*I[M,Q]\
                                                 +TAU11*I[K,Q]*I[L,P]*I[M,N]
        Ct = hex8.reduce_tensor_to_vector_form(Ct)
        
        self.assertEqual(np.allclose(C,Ct),True)

    def test_form_D(self):
        """Test forming the D stiffness tensor"""
        TAU   = 5.1
        SIGMA = 2.
        
        I = np.eye(3)
        
        D = form_D(TAU,SIGMA)
        
        Dt = np.zeros([3,3,3,3])
        
        for K in range(3):
            for L in range(3):
                for M in range(3):
                    for N in range(3):
                        Dt[K,L,M,N] = TAU*I[K,L]*I[M,N]+SIGMA*(I[K,M]*I[L,N]+I[K,N]*I[L,M])
        Dt = hex8.reduce_tensor_to_vector_form(Dt)
        
        self.assertEqual(np.allclose(D,Dt),True)
        
    def test_form_stiffness_tensors(self):
        """Test forming the stiffness tensors"""
        LAMBDA = 2.4
        MU     = 6.7
        ETA    = 2.4
        TAU    = 5.1
        KAPPA  = 5.6
        NU     = 8.2
        SIGMA  = 2.
        TAUS   = [4.5,1.3,9.2,1.1,6.4,2.4,7.11,5.5,1.5,3.8,2.7]
        PARAMS = LAMBDA,MU,ETA,TAU,KAPPA,NU,SIGMA,TAUS
        
        A,B,C,D = form_stiffness_tensors(PARAMS)
        
        At = form_A(LAMBDA,MU)
        Bt = form_B(ETA,TAU,KAPPA,NU,SIGMA)
        Ct = form_C(TAUS)
        Dt = form_D(TAU,SIGMA)
        
        self.assertTrue(np.allclose(A,At))
        self.assertTrue(np.allclose(B,Bt))
        self.assertTrue(np.allclose(C,Ct))
        self.assertTrue(np.allclose(D,Dt))
        
    def test_compute_strain_measures(self):
        """Test the computation of the strain measures"""
        F   = np.array(range(9))
        chi = np.array(range(9,18))
        
        C   = hex8.matrix_Tdot_V(F,F)
        Phi = hex8.matrix_Tdot_V(F,chi)
        I   = hex8.reduce_tensor_to_vector_form(np.eye(3))
        
        E_macro,E_micro = compute_strain_measures(C,Phi)
        
        E_macrot = 0.5*(C-I)
        E_microt = Phi-I
        
        self.assertTrue(np.allclose(E_macro,E_macrot))
        self.assertTrue(np.allclose(E_micro,E_microt))
        
    def test_compute_dCinvdC(self):
        """Test computinig the derivative of Cinv w.r.t. C"""
        Cinv  = np.array(range(9))
        Cinvt = hex8.convert_V_to_T(Cinv,[3,3])
        
        dCinvdC = compute_dCinvdC(Cinv)
        
        dCinvdCt = np.zeros([3,3,3,3])
        for I in range(3):
            for J in range(3):
                for K in range(3):
                    for L in range(3):
                        dCinvdCt[I,J,K,L] = -Cinvt[I,K]*Cinvt[L,J]
        dCinvdCt = hex8.reduce_tensor_to_vector_form(dCinvdCt)
        
        self.assertTrue(np.allclose(dCinvdC,dCinvdCt))
        
    def test_compute_pk2_stress(self):
        """Test the computation of the pk2 stress"""
        LAMBDA = 2.4
        MU     = 6.7
        ETA    = 2.4
        TAU    = 5.1
        KAPPA  = 5.6
        NU     = 8.2
        SIGMA  = 2.
        TAUS   = [4.5,1.3,9.2,1.1,6.4,2.4,7.11,5.5,1.5,3.8,2.7]
        PARAMS = LAMBDA,MU,ETA,TAU,KAPPA,NU,SIGMA,TAUS
        
        AS,BS,CS,DS = form_stiffness_tensors(PARAMS)
        
        F        = np.array([1,4,2,2,2,3,3,2,1]).astype(float)
        chi      = np.array(range(9,18)).astype(float)
        grad_chi = np.array(range(18,45)).astype(float)
        
        C   = hex8.matrix_Tdot_V(F,F)
        Cinv = hex8.invert_3x3_matrix_V(C)[1]
        Phi = hex8.matrix_Tdot_V(F,chi)
        
        Gamma = np.zeros([27,])
        
        for I in range(3):
            for J in range(3):
                for K in range(3):
                    IJK = T2V([I,J,K],[3,3,3])
                    for i in range(3):
                        iI  = T2V([i,I],[3,3])
                        iJK = T2V([i,J,K],[3,3,3])
                        Gamma[IJK] += F[iI]*grad_chi[iJK]
        
        Iten = hex8.reduce_tensor_to_vector_form(np.eye(3))
        
        E_micro,E_macro = compute_strain_measures(C,Phi)
        
        pk2   = compute_pk2_stress(AS,BS,CS,DS,E_macro,E_micro,C,Gamma)[0]
        
        St  = np.zeros([9,])
        
        for I in range(3):
            for J in range(3):
                IJ = T2V([I,J],[3,3])
                for K in range(3):
                    for L in range(3):
                        IJKL = T2V([I,J,K,L],[3,3,3,3])
                        KL   = T2V([K,L],        [3,3])
                        St[IJ] +=  AS[IJKL]*E_macro[KL] + DS[IJKL]*E_micro[KL]
                        
        for I in range(3):
            for J in range(3):
                IJ = T2V([I,J],[3,3])
                for K in range(3):
                    for L in range(3):
                        KL   = T2V([K,L],[3,3])
                        for Q in range(3):
                            for R in range(3):
                                IQKL = T2V([I,Q,K,L],[3,3,3,3])
                                RQ   = T2V([R,Q],[3,3])
                                JR   = T2V([J,R],[3,3])
                                St[IJ] += (BS[IQKL]*E_micro[KL] + DS[IQKL]*E_macro[KL])*(E_micro[RQ] + Iten[RQ])*Cinv[JR]
        
        for I in range(3):
            for J in range(3):
                IJ = T2V([I,J],[3,3])
                for L in range(3):
                    for M in range(3):
                        for N in range(3):
                            LMN = T2V([L,M,N],[3,3,3])
                            for Q in range(3):
                                for R in range(3):
                                    for S in range(3):
                                        IQRLMN = T2V([I,Q,R,L,M,N],[3,3,3,3,3,3])
                                        JS     = T2V([J,S],[3,3])
                                        SQR    = T2V([S,Q,R],[3,3,3])
                                        St[IJ] += CS[IQRLMN]*Gamma[LMN]*Cinv[JS]*Gamma[SQR]
        
        self.assertTrue(np.allclose(pk2,St))

    def test_compute_symmetric_stress(self):
        """Test the computation of the symmetric stress"""
        LAMBDA = 2.4
        MU     = 6.7
        ETA    = 2.4
        TAU    = 5.1
        KAPPA  = 5.6
        NU     = 8.2
        SIGMA  = 2.
        TAUS   = [4.5,1.3,9.2,1.1,6.4,2.4,7.11,5.5,1.5,3.8,2.7]
        PARAMS = LAMBDA,MU,ETA,TAU,KAPPA,NU,SIGMA,TAUS
        
        AS,BS,CS,DS = form_stiffness_tensors(PARAMS)
        
        F        = np.array([1,4,2,2,2,3,3,2,1]).astype(float)
        chi      = np.array(range(9,18)).astype(float)
        grad_chi = np.array(range(18,45)).astype(float)
        
        C   = hex8.matrix_Tdot_V(F,F)
        Cinv = hex8.invert_3x3_matrix_V(C)[1]
        Phi = hex8.matrix_Tdot_V(F,chi)
        
        Gamma = np.zeros([27,])
        
        for I in range(3):
            for J in range(3):
                for K in range(3):
                    IJK = T2V([I,J,K],[3,3,3])
                    for i in range(3):
                        iI  = T2V([i,I],[3,3])
                        iJK = T2V([i,J,K],[3,3,3])
                        Gamma[IJK] += F[iI]*grad_chi[iJK]
        
        Iten = hex8.reduce_tensor_to_vector_form(np.eye(3))
        
        E_micro,E_macro = compute_strain_measures(C,Phi)
        
        _,TEMP1,TEMP2   = compute_pk2_stress(AS,BS,CS,DS,E_macro,E_micro,C,Gamma)
        symm_stress     = compute_symmetric_stress(TEMP1,TEMP2)
        
        symm_stresst = np.zeros([9,])
        
        for I in range(3):
            for J in range(3):
                IJ = T2V([I,J],[3,3])
                for K in range(3):
                    for L in range(3):
                        IJKL = T2V([I,J,K,L],[3,3,3,3])
                        KL   = T2V([K,L],        [3,3])
                        symm_stresst[IJ] += AS[IJKL]*E_macro[KL]+DS[IJKL]*E_micro[KL]
                        for Q in range(3):
                            for R in range(3):
                                IQKL = T2V([I,Q,K,L],[3,3,3,3])
                                JQKL = T2V([J,Q,K,L],[3,3,3,3])
                                RQ = T2V([R,Q],[3,3])
                                JR = T2V([J,R],[3,3])
                                IR = T2V([I,R],[3,3])
                                tmp1 = (BS[IQKL]*E_micro[KL] + DS[IQKL]*E_macro[KL])*(E_micro[RQ]+Iten[RQ])*Cinv[JR]
                                tmp2 = (BS[JQKL]*E_micro[KL] + DS[JQKL]*E_macro[KL])*(E_micro[RQ]+Iten[RQ])*Cinv[IR]
                                symm_stresst[IJ] += tmp1+tmp2
        for I in range(3):
            for J in range(3):
                IJ = T2V([I,J],[3,3])
                for L in range(3):
                    for M in range(3):
                        for N in range(3):
                            LMN = T2V([L,M,N],[3,3,3])
                            for Q in range(3):
                                for R in range(3):
                                    for S in range(3):
                                        IQRLMN = T2V([I,Q,R,L,M,N],[3,3,3,3,3,3])
                                        JQRLMN = T2V([J,Q,R,L,M,N],[3,3,3,3,3,3])
                                        SQR    = T2V([S,Q,R],[3,3,3])
                                        JS     = T2V([J,S],[3,3])
                                        IS     = T2V([I,S],[3,3])
                                        tmp1   = CS[IQRLMN]*Gamma[LMN]*Gamma[SQR]*Cinv[JS]
                                        tmp2   = CS[JQRLMN]*Gamma[LMN]*Gamma[SQR]*Cinv[IS]
                                        symm_stresst[IJ] += tmp1 + tmp2
        self.assertTrue(np.allclose(symm_stress,symm_stresst))
    
    def test_compute_ho_stress(self):
        """Test the computation of the higher order stress"""
        LAMBDA = 2.4
        MU     = 6.7
        ETA    = 2.4
        TAU    = 5.1
        KAPPA  = 5.6
        NU     = 8.2
        SIGMA  = 2.
        TAUS   = [4.5,1.3,9.2,1.1,6.4,2.4,7.11,5.5,1.5,3.8,2.7]
        PARAMS = LAMBDA,MU,ETA,TAU,KAPPA,NU,SIGMA,TAUS
        
        AS,BS,CS,DS = form_stiffness_tensors(PARAMS)
        
        F        = np.array([1,4,2,2,2,3,3,2,1]).astype(float)
        chi      = np.array(range(9,18)).astype(float)
        grad_chi = np.array(range(18,45)).astype(float)
        
        C   = hex8.matrix_Tdot_V(F,F)
        Cinv = hex8.invert_3x3_matrix_V(C)[1]
        Phi = hex8.matrix_Tdot_V(F,chi)
        
        Gamma = np.zeros([27,])
        
        for I in range(3):
            for J in range(3):
                for K in range(3):
                    IJK = T2V([I,J,K],[3,3,3])
                    for i in range(3):
                        iI  = T2V([i,I],[3,3])
                        iJK = T2V([i,J,K],[3,3,3])
                        Gamma[IJK] += F[iI]*grad_chi[iJK]
        Mten  = compute_ho_stress(CS,Gamma)
        Mt = np.zeros([27,])
        for I in range(3):
            for J in range(3):
                for K in range(3):
                    IJK = T2V([I,J,K],[3,3,3])
                    for L in range(3):
                        for M in range(3):
                            for N in range(3):
                                IJKLMN = T2V([I,J,K,L,M,N],[3,3,3,3,3,3])
                                LMN    = T2V([L,M,N],[3,3,3])
                                Mt[IJK] += CS[IJKLMN]*Gamma[LMN]
        self.assertTrue(np.allclose(Mten,Mt))

if __name__ == '__main__':
    unittest.main()