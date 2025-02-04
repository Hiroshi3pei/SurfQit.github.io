import numpy as np
cimport numpy as cnp
cimport cython
from cython import boundscheck, wraparound

ctypedef cnp.float64_t DTYPE_t
ctypedef cnp.int32_t RETURN_t

cpdef inline cnp.ndarray[RETURN_t, ndim=2] C_position_to_index(cnp.ndarray[DTYPE_t, ndim=2] host, cnp.ndarray[DTYPE_t, ndim=3] positions):
    cdef cnp.ndarray[RETURN_t, ndim=2] indice
    with boundscheck(False), wraparound(False):
        indice = np.empty((positions.shape[0],host.shape[0]),np.int32)
        for operation in range(positions.shape[0]):
            indice[operation] = np.argmin(np.sum((positions[operation][np.newaxis,:,:]%1 - host[:,np.newaxis,:]) ** 2, axis=2),axis=1)
    return indice.T

