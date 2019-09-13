"""
This module contains functions to calculate an (euclidean) and influence weighted distance transform 
"""

# cython: boundscheck= False, wraparound= False, cdivision= True, language_level=3, binding= True, embedsignature= True

import numpy as np
import netsim.chamfer as ch
cimport numpy as np
cimport cython
from cpython cimport bool

cdef double __poly(double[:] coef, double x):
    '''
    Internal function used to calculate a polynomial

    Paremeters
    ----------
        coef: array of doubles
            polynomial coefficients
        
        x: double
            value used at which the polynomial is calculated
    
    Returns
    -------
        y: double
            p(x) for a given set of coefficients
    '''
    cdef:
        double y = 0.0
        Py_ssize_t ncoef = coef.shape[0]
        Py_ssize_t i = 0
    while i < ncoef:
        y += x**(ncoef - i -1) * coef[i]
        i += 1
    return y

cdef cy_calculate_iwdt(double[:,:] iwdt_tmp, dict iwdt_dict, int option = 1):
    
    # convert iwdt to continguous
    cdef double[:,::1] iwdt = np.ascontiguousarray(iwdt_tmp)
    
    # unpack iwdt dict
    cdef:
        double [:,::1] dem = np.ascontiguousarray(iwdt_dict['dem'])
        double [:,::1] netcost = np.ascontiguousarray(iwdt_dict['netcost'])
        double [:] coef = iwdt_dict['coef']
        float cellsize = iwdt_dict['cellsize']
        float weight = iwdt_dict['weight']
       
    # initialize old_iwdt and backlinks
    cdef:
        double[:,::1] old_iwdt = np.ascontiguousarray(np.zeros((iwdt.shape[0],iwdt.shape[1]), np.float64))
        int[:,::1] blx = np.ascontiguousarray(np.zeros((iwdt.shape[0],iwdt.shape[1]), np.intc))
        int[:,::1] bly = np.ascontiguousarray(np.zeros((iwdt.shape[0],iwdt.shape[1]), np.intc))

    # calculate max_cost
    cdef:
        double[:,::1] deltax, deltay, all_gradients
        float max_gradient
        double max_cost

    deltay, deltax = np.gradient(dem, cellsize)
    all_gradients = np.arctan(np.sqrt(np.square(deltax) + np.square(deltay))) 
    max_gradient = np.float(np.tan(np.percentile(all_gradients, 75)))
    max_cost = __poly(coef, max_gradient)
          
    # calculate loop limits
    cdef:
        Py_ssize_t nsize = ch.nsize
        Py_ssize_t nrows = iwdt.shape[0] - nsize
        Py_ssize_t ncols = iwdt.shape[1] - nsize    
    
    # initialize convergence flag
    cdef bool converged = False

    # initialize variables inside all loops
    cdef:
        Py_ssize_t scan_dir
        Py_ssize_t r, c
        Py_ssize_t nr, nc
        Py_ssize_t n
        Py_ssize_t ncoef = ch.ncoef
        double[:,::1] ldm = ch.ldm.astype(np.float64)
        double threshold = ch.threshold
        int[:,::1] dx = ch.dx
        int[:,::1] dy = ch.dy       
        double d0, d1
        double current_gradient, gradient_cost, final_cost 

    while not converged:
        
        # forward
        scan_dir= 0
        
        for r in range(nsize, nrows):
            for c in range(nsize, ncols):
                
                    # initialize to current cost @ source
                    d0 = iwdt[r, c]
                    
                    for n in range(ncoef):
                        
                        # generate neighbor
                        nr = r + dx[scan_dir, n]
                        nc = c + dy[scan_dir, n]

                        # calculate cost to neighbor
                        current_gradient = (dem[nr,nc] - dem[r,c]) / (ldm[scan_dir, n] * cellsize)
                        gradient_cost = __poly(coef, current_gradient) / max_cost
                        final_cost = netcost[nr,nc]*weight + gradient_cost*(1.0 - weight)

                        # new accumulated cost (cost existing at neighbor plus cost to go to neighbor)
                        d1 = iwdt[nr,nc] + final_cost * ldm[scan_dir, n]

                        # update w new cost?
                        if d1 < d0:
                            d0 =  d1
                            blx[r, c] = dx[scan_dir, n]
                            bly[r, c] = dy[scan_dir, n]
               
                    # update with new cost
                    iwdt[r,c] = d0                  
        
        # backward
        scan_dir= 1
        
        for r in range(nrows-1, nsize-1, -1):
            for c in range(ncols-1, nsize-1, -1):
            
                    # initialize to current cost @ source
                    d0 = iwdt[r, c]
            
                    for n in range(ncoef):
                        
                        # generate neighbor
                        nr = r + dx[scan_dir, n]
                        nc = c + dy[scan_dir, n]

                         # calculate cost to neighbor
                        current_gradient = (dem[nr,nc] - dem[r,c]) / (ldm[scan_dir, n] * cellsize)
                        gradient_cost = __poly(coef, current_gradient) / max_cost
                        final_cost = netcost[nr,nc]*weight + gradient_cost*(1.0 - weight)

                        # new accumulated cost (cost existing at neighbor plus cost to neighbor)
                        d1 = iwdt[nr,nc] + final_cost  * ldm[scan_dir, n]

                        # update w new cost?
                        if d1 < d0:
                            d0 =  d1
                            blx[r, c] = dx[scan_dir, n]
                            bly[r, c] = dy[scan_dir, n]

                    # update with new cost
                    iwdt[r,c] = d0
        
        # convergence?
        converged = <bool> np.all(np.abs(np.subtract(old_iwdt, iwdt)) < threshold)
        
        # copy old iwdt values
        old_iwdt[:,:] = iwdt
        
    # end message
    #print('\n end of iwdt!')

    # add padding to edges
    iwdt = np.pad(iwdt[nsize:nrows, nsize:ncols], pad_width=((nsize, nsize), (nsize,nsize)), mode='mean')
    blx = np.pad(blx[nsize:nrows, nsize:ncols], pad_width=((nsize, nsize), (nsize,nsize)), mode='edge')
    bly = np.pad(bly[nsize:nrows, nsize:ncols], pad_width=((nsize, nsize), (nsize,nsize)), mode='edge')

   
    if option == 1:
        return np.asarray(iwdt), np.asarray(blx), np.asarray(bly)
    elif option == 2:
        return np.asarray(iwdt)
    else:
        return np.asarray(blx), np.asarray(bly)

    
def calculate_iwdt(double[:,:] iwdt, dict iwdt_dict, int option = 1):
    '''
    calculates influence weighted distance transform

    Parameters
    ----------
    iwdt: 2D numpy array
        array with source cells initiated to zero and remaining cells to 
        a large none-zero value (typically 999999.0)
    
    iwdt_dict: dictionary
        dictionary containing information needed to calculate iwdt (see *Notes* below)
    
    option: int
        parameter used to determine the output:

        1. Returns *iwdt*, *blx*, *bly* (*backlinks*)
        2. Return *iwdt* only
        3. Returns *blx*, *bly* (*backlinks*) 

    Returns
    -------
    iwdt: 2D numpy array, *optional*
        influence weighted distance transform
    
    blx, bly: 2D numpy arrays, *optional*
        influence weighted distance transform backlinks

    Notes
    -----

    Wrapper function for the cython function ``cy_calculate_iwdt()``. To calculate
    *influence weighted distance transform*, it is necessary to generate a dictionary (``iwdt_dict{}``) with 
    the following entries:
    
    - **'dem:'** - 2D numpy array
      array with elevation data
    - **'netcost:'** -  2D numpy array
      array with values between 0 and 1. **Must have the same dimensions as *dem**
    - **'cellsize:'** - float
      cellsize
    - **'weight:'**- float
      weight associated with netcost. Must be a value between 0 and 1
    - **'coef:'**- numpy array 
      coefficients for polynomial mapping gradient to cost

    '''

    return cy_calculate_iwdt(iwdt, iwdt_dict, option)



cdef cy_calculate_dt(double[:,:] dt_tmp, float cellsize, int option = 1):

     # convert iwdt to continguous
    cdef double[:,::1] dt = np.ascontiguousarray(dt_tmp)  
    
    # initialize backlinks
    cdef:
        int[:,::1] blx = np.ascontiguousarray(np.zeros((dt.shape[0],dt.shape[1]), np.intc))
        int[:,::1] bly = np.ascontiguousarray(np.zeros((dt.shape[0],dt.shape[1]), np.intc))        
    
    # calculate loop limits
    cdef:
        Py_ssize_t nsize = ch.nsize
        Py_ssize_t nrows = dt.shape[0] - nsize
        Py_ssize_t ncols = dt.shape[1] - nsize    
    

    # initialize variables inside all loops
    cdef:
        Py_ssize_t scan_dir
        Py_ssize_t r, c
        Py_ssize_t nr, nc
        Py_ssize_t n
        Py_ssize_t ncoef = ch.ncoef
        double[:,::1] ldm = ch.ldm.astype(np.float64)
        int[:,::1] dx = ch.dx
        int[:,::1] dy = ch.dy       
        double d0, d1

    # forward
    scan_dir= 0

    for r in range(nsize, nrows):
        for c in range(nsize, ncols):

                # initialize to current distance @ source
                d0 = dt[r, c]

                for n in range(ncoef):

                    # generate neighbor
                    nr = r + dx[scan_dir, n]
                    nc = c + dy[scan_dir, n]

                    # new accumulated distance to neighbor
                    d1 = dt[nr,nc] + ldm[scan_dir, n]*cellsize

                    # update w new distance?
                    if d1 < d0:
                        d0 =  d1
                        blx[r, c] = dx[scan_dir, n]
                        bly[r, c] = dy[scan_dir, n]

                # update with new distance
                dt[r,c] = d0                  

    # backward
    scan_dir= 1

    for r in range(nrows-1, nsize-1, -1):
        for c in range(ncols-1, nsize-1, -1):

                # initialize to current distance @ source
                d0 = dt[r, c]

                for n in range(ncoef):

                    # generate neighbor
                    nr = r + dx[scan_dir, n]
                    nc = c + dy[scan_dir, n]

                    # new accumulated distance to neighbor
                    d1 = dt[nr,nc] + ldm[scan_dir, n]*cellsize

                    # update w new distance?
                    if d1 < d0:
                        d0 =  d1
                        blx[r, c] = dx[scan_dir, n]
                        bly[r, c] = dy[scan_dir, n]

                # update with new distance
                dt[r,c] = d0                  

        
    # end message
    #print('\nend of dt!')

    # add padding to edges
    dt =  np.pad( dt[nsize:nrows, nsize:ncols], pad_width=((nsize, nsize), (nsize,nsize)), mode='mean')
    blx = np.pad(blx[nsize:nrows, nsize:ncols], pad_width=((nsize, nsize), (nsize,nsize)), mode='edge')
    bly = np.pad(bly[nsize:nrows, nsize:ncols], pad_width=((nsize, nsize), (nsize,nsize)), mode='edge')
    
    if option == 1:
        return np.asarray(dt), np.asarray(blx), np.asarray(bly)
    elif option == 2:
        return np.asarray(dt)
    else:
        return np.asarray(blx), np.asarray(bly)

def calculate_dt(double[:,:] dt, float cellsize, int option = 1):
    '''
    calculates euclidean distance transform. 
    
    Parameters
    ----------

    dt: 2D numpy array
       array with source cells initiated to zero and remaining cells to 
       a large none-zero value (typically 999999.0)

    cellsize: float
        cellsize
    
    option: int
        parameter used to determine the output:

        1. Returns *dt*, *blx*, *bly* (*backlinks*)
        2. Return *dt* only
        3. Returns *blx*, *bly* (*backlinks*) 
        
    Returns
    -------
    dt: 2D numpy array, *optional*
        distance transform
    
    blx, bly: 2D numpy arrays, *optional*
        distance transform backlinks
    
    Notes
    -----
    This is a wrapper for the cython function ``cy_calculate_dt()``

    '''
    return cy_calculate_dt(dt, cellsize, option)