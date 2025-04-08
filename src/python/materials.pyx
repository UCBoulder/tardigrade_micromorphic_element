from libcpp.vector cimport vector
from libcpp.string cimport string

import numpy as np
cimport numpy as np

cimport materials

# Function definitions

cdef map_array_to_vector(np.ndarray A, shape):
    """
    Map a numpy nd array to a vector

    :param np.ndarray A: The array to be mapped
    """

    cdef vector[double] V

    if (shape > 1):
        return [map_array_to_vector(a, shape-1) for a in A]
    
    return [a for a in A]

cdef map_vector_to_array(vector[double] V):
    """
    Map a vector to a numpy array

    :param vector[double] V: The vector to be mapped
    """

    cdef np.ndarray array = np.zeros(V.size())
    
    for i in range(array.size):
        array[i] = V[i]

    return array

cdef map_matrix_to_array(vector[vector[double]] M):

    return np.array([map_vector_to_array(V) for V in M])

cdef map_3matrix_to_array(vector[vector[vector[double]]] M3):

    return [map_matrix_to_array(M) for M in M3]


def evaluate_model(object model_name,
                   np.ndarray time, np.ndarray fparams,
                   np.ndarray current_grad_u,  np.ndarray current_phi,  np.ndarray current_grad_phi,
                   np.ndarray previous_grad_u, np.ndarray previous_phi, np.ndarray previous_grad_phi,
                   np.ndarray SDVS,
                   np.ndarray current_ADD_DOF,
                   np.ndarray current_ADD_grad_DOF,
                   np.ndarray previous_ADD_DOF,
                   np.ndarray previous_ADD_grad_DOF):
    """
    The python wrapper for the evaluate model function

    A wrapper for the evaluate model method of micromorphic_material_library::IMaterial which also includes the call to the material name
    
    Evaluate the jacobian of the material model using a numeric gradient

    :param str model_name: The name of the model to be evaluated
    :param np.ndarray time: The current time and the timestep
    :param np.ndarray fparams: The parameters for the constitutive model
    :param np.ndarray current_grad_u: The current displacement gradient
        Assumed to be of the form [ [ u_{1,1}, u_{1,2}, u_{1,3} ],
                                    [ u_{2,1}, u_{2,2}, u_{2,3} ],
                                    [ u_{3,1}, u_{3,2}, u_{3,3} ] ]
    :param np.ndarray current_phi: The current micro displacment values.
        Assumed to be of the form [ \phi_{11}, \phi_{12}, \phi_{13}, \phi_{21}, \phi_{22}, \phi_{23}, \phi_{31}, \phi_{32}, \phi_{33} ]
    :param np.ndarray current_grad_phi: The current micro displacement gradient
        Assumed to be of the form [ [ \phi_{11,1}, \phi_{11,2}, \phi_{11,3} ],
                                    [ \phi_{12,1}, \phi_{12,2}, \phi_{12,3} ],
                                    [ \phi_{13,1}, \phi_{13,2}, \phi_{13,3} ],
                                    [ \phi_{21,1}, \phi_{21,2}, \phi_{21,3} ],
                                    [ \phi_{22,1}, \phi_{22,2}, \phi_{22,3} ],
                                    [ \phi_{23,1}, \phi_{23,2}, \phi_{23,3} ],
                                    [ \phi_{31,1}, \phi_{31,2}, \phi_{31,3} ],
                                    [ \phi_{32,1}, \phi_{32,2}, \phi_{32,3} ],
                                    [ \phi_{33,1}, \phi_{33,2}, \phi_{33,3} ] ]
    :param np.ndarray previous_grad_u: The previous displacement gradient.
    :param np.ndarray previous_phi: The previous micro displacement.
    :param np.ndarray previous_grad_phi: The previous micro displacement gradient.
    :param np.ndarray SDVS: The previously converged values of the state variables
    :param np.ndarray current_ADD_DOF: The current values of the additional degrees of freedom
    :param np.ndarray current_ADD_grad_DOF: The current values of the gradients of the
        additional degrees of freedom
    
    :returns:
        0: No errors. Solution converged.
        1: Convergence Error. Request timestep cutback.
        2: Fatal Errors encountered. Terminate the simulation.
    """

    cdef string c_model_name = model_name.encode('UTF-8')

    cdef vector[double] c_time                          = map_array_to_vector(time, 1)
    cdef vector[double] c_fparams                       = map_array_to_vector(fparams, 1)
    cdef double[3][3]   c_current_grad_u
    cdef double[9]      c_current_phi                   = current_phi
    cdef double[9][3]   c_current_grad_phi
    cdef double[3][3]   c_previous_grad_u
    cdef double[9]      c_previous_phi                  = previous_phi
    cdef double[9][3]   c_previous_grad_phi
    cdef vector[double] c_SDVS                          = map_array_to_vector(SDVS, 1)
    cdef vector[double] c_current_ADD_DOF               = map_array_to_vector(current_ADD_DOF, 1)
    cdef vector[vector[double]] c_current_ADD_grad_DOF  = map_array_to_vector(current_ADD_grad_DOF, 2)
    cdef vector[double] c_previous_ADD_DOF              = map_array_to_vector(current_ADD_DOF, 1)
    cdef vector[vector[double]] c_previous_ADD_grad_DOF = map_array_to_vector(previous_ADD_grad_DOF, 2)
    cdef vector[double] c_PK2
    cdef vector[double] c_SIGMA
    cdef vector[double] c_M
    cdef vector[vector[double]] c_DPK2Dgrad_u
    cdef vector[vector[double]] c_DPK2Dphi
    cdef vector[vector[double]] c_DPK2Dgrad_phi
    cdef vector[vector[double]] c_DSIGMADgrad_u
    cdef vector[vector[double]] c_DSIGMADphi
    cdef vector[vector[double]] c_DSIGMADgrad_phi
    cdef vector[vector[double]] c_DMDgrad_u
    cdef vector[vector[double]] c_DMDphi
    cdef vector[vector[double]] c_DMDgrad_phi
    cdef vector[vector[double]] c_ADD_TERMS
    cdef vector[vector[vector[double]]] c_ADD_JACOBIANS
    cdef string c_output_message
    cdef int errorCode
    cdef np.ndarray PK2
    cdef np.ndarray SIGMA
    cdef np.ndarray M
    cdef np.ndarray DPK2Dgrad_u
    cdef np.ndarray DPK2Dphi
    cdef np.ndarray DPK2Dgrad_phi
    cdef np.ndarray DSIGMADgrad_u
    cdef np.ndarray DSIGMADphi
    cdef np.ndarray DSIGMADgrad_phi
    cdef np.ndarray DMDgrad_u
    cdef np.ndarray DMDphi
    cdef np.ndarray DMDgrad_phi
    cdef np.ndarray ADD_TERMS
    cdef list ADD_JACOBIANS
    cdef object output_message

    # Map the grad u variables
    for i in range(3):
        for j in range(3):
            c_current_grad_u[i][j]  = current_grad_u[i,j]
            c_previous_grad_u[i][j] = previous_grad_u[i,j]

    # Map the grad phi variables
    for i in range(9):
        for j in range(3):
            c_current_grad_phi[i][j]  = current_grad_phi[i,j]
            c_previous_grad_phi[i][j] = previous_grad_phi[i,j]

    errorCode = materials.evaluate_model(c_model_name,\
                                         c_time, c_fparams,\
                                         c_current_grad_u,  c_current_phi,  c_current_grad_phi,\
                                         c_previous_grad_u, c_previous_phi, c_previous_grad_phi,\
                                         c_SDVS,\
                                         c_current_ADD_DOF,\
                                         c_current_ADD_grad_DOF,\
                                         c_previous_ADD_DOF,\
                                         c_previous_ADD_grad_DOF,\
                                         c_PK2, c_SIGMA, c_M,\
                                         c_DPK2Dgrad_u,   c_DPK2Dphi,   c_DPK2Dgrad_phi,\
                                         c_DSIGMADgrad_u, c_DSIGMADphi, c_DSIGMADgrad_phi,\
                                         c_DMDgrad_u, c_DMDphi, c_DMDgrad_phi,\
                                         c_ADD_TERMS, c_ADD_JACOBIANS, c_output_message);

    PK2             = map_vector_to_array(c_PK2)
    SIGMA           = map_vector_to_array(c_SIGMA)
    M               = map_vector_to_array(c_M)

    SDVS            = map_vector_to_array(c_SDVS)

    DPK2Dgrad_u     = map_matrix_to_array(c_DPK2Dgrad_u)
    DPK2Dphi        = map_matrix_to_array(c_DPK2Dphi)
    DPK2Dgrad_phi   = map_matrix_to_array(c_DPK2Dgrad_phi)

    DSIGMADgrad_u   = map_matrix_to_array(c_DSIGMADgrad_u)
    DSIGMADphi      = map_matrix_to_array(c_DSIGMADphi)
    DSIGMADgrad_phi = map_matrix_to_array(c_DSIGMADgrad_phi)

    DMDgrad_u       = map_matrix_to_array(c_DMDgrad_u)
    DMDphi          = map_matrix_to_array(c_DMDphi)
    DMDgrad_phi     = map_matrix_to_array(c_DMDgrad_phi)

    ADD_TERMS       = map_matrix_to_array(c_ADD_TERMS)

    ADD_JACOBIANS   = map_3matrix_to_array(c_ADD_JACOBIANS)

    output_message  = c_output_message

    return errorCode, PK2, SIGMA, M, SDVS,\
           DPK2Dgrad_u, DPK2Dphi, DPK2Dgrad_phi,\
           DSIGMADgrad_u, DSIGMADphi, DSIGMADgrad_phi,\
           DMDgrad_u, DMDphi, DMDgrad_phi,\
           ADD_TERMS, ADD_JACOBIANS, output_message
