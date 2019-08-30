"""
This module contains the main function to run a Network Simulation ``simulation()``

"""

import numpy as np
import netsim.path_tools as pt
from collections import OrderedDict
from netsim.iwdt import calculate_iwdt, calculate_dt


def simulation(pts, net_layout, cost_dict, netsim_dict):
    '''
    Generates a network of paths.
    
    Parameters
    ----------
    
    pts: dataframe
        contains the identifier, row and column of each location
    
    net_layout: dataframe
        dataframe specifying origin and destination of each path in the network
    
    cost_dict: dictionary
        contains parameters used for ``calculate_iwdt()``
    
    netsim_dict: dictionary
        contains parameters needed to execute network simulation


    Returns
    -------
    
    Gt: 2D numpy array
        final ground potential
    
    paths: 2D numpy array
        sum of all network paths
        
    paths_dict: dictionary
        dictionary with the track of every path in network


    Notes
    -----
    
    The ``netsim_dict{}`` must contain the following entries:

    - 'i: ' - float
      effect of new path. *Default:* 1.0
    - 'Gmax: ' - float
      maximum ground potential.
    - 'T: ' - float
      1/T represents the residuality of a path.
    - 'alpha: ' - float
      coefficient calculated from :math:`\\alpha_1 = \\frac {d_0} {ln(1 - NC_0)}`
    

    While the ``cost_dict{}`` must contain the same entries as ``iwdt_dict{}`` see **IWDT** module

        
    '''
    
    # unpack netsim dictionary
    i = netsim_dict['i']
    Gmax = netsim_dict['Gmax']
    T = netsim_dict['T']
    alpha = netsim_dict['alpha']
    
    # initialize netsim outputs
    Gt = np.zeros_like(cost_dict['dem'])
    Gt_1 = np.zeros_like(Gt)
    paths = np.zeros_like(Gt)
    
    # initialize paths dictionary
    path_lst = [] #paths_dict= {}#OrderedDict()
    
    for pth_id, pth_def in net_layout.iterrows():
        
        # origin & destination?
        o = pth_def['origin']
        d = pth_def['destination']

        # retrieve location @ origin & destination
        origin      = [[ pts.loc[pts['id'] == o, 'r'].values[0], pts.loc[pts['id'] == o,'c'].values[0] ]]
        destination = [[ pts.loc[pts['id'] == d, 'r'].values[0], pts.loc[pts['id'] == d,'c'].values[0] ]]
                
        # initialize ACS
        acs = np.full_like(Gt, 999999.0)
        for r,c in origin:
            acs[r,c] = 0.0
        
        # calculated influence weighted distance transform
        acs, blx, bly = calculate_iwdt(acs, cost_dict)
        
        # create new path to destination
        path_t, path_info = pt.create_paths(blx, bly, origin, destination, start_path=pth_id)

        # update paths
        paths += path_t
        path_lst.append(path_info) #paths_dict.update(path_info)
        
        # update ground potential
        Gt = Gt_1 - (Gt_1/T) + path_t * i * (1 - (Gt_1 / Gmax))
              
        # update network cost
        d = np.full_like(Gt, 99999.0)
        d[Gt >= 1.0] = 0.0
        d = calculate_dt(d, cost_dict['cellsize'], option=2)
        cost_dict['netcost'] = 1.0 - np.exp(d / alpha)
        
        # update Gt_1
        Gt_1 = np.copy(Gt)
        
    return Gt, paths, path_lst #paths_dict