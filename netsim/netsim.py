import numpy as np
import netsim.path_tools as pt
from netsim.iwdt import calculate_iwdt, calculate_dt


def netsim(pts, net_layout, cost_dict, netsim_dict):
    '''
    Generates a path network.
    
    Parameters
    ----------
    
    pts: dataframe
        contains the identifier, and row and column of each location
    
    net_layout: dataframe
        dataframe specifying origin and destination of each path network
    
    cost_dict: dictionary
        contains parameters used for *calculate_iwdt()*
    
    netsim_dict: dictionary
        contains parameters needed to execute netsim.
        
    Returns
    -------
    
    G_t: 2D numpy array
        final ground potential
    
    paths: 2D numpy array
        sum of all network paths
        
    paths_dict: dictionary
        dictionary with the track of every path
        
    Notes
    -----
    
    Have a function that runs checks netsim
    
    '''
    
    # unpack netsim dictionary
    i = netsim_dict['i']
    G_max = netsim_dict['G_max']
    ri = netsim_dict['ri']
    alpha = netsim_dict['alpha']
    
    # initialize netsim outputs
    G_t = np.zeros_like(cost_dict['DEM'])
    G_t_1 = np.zeros_like(G_t)
    paths = np.zeros_like(G_t)
    
    # initialize paths dictionary
    paths_dict= OrderedDict()
    
    for pth_id, pth_def in net_layout.iterrows():
        
        # origin & destination?
        o = pth_def['origin']
        d = pth_def['destination']

        # retrieve location @ origin & destination
        origin      = [[ pts.loc[pts['id'] == o, 'r'].values[0], pts.loc[pts['id'] == o,'c'].values[0] ]]
        destination = [[ pts.loc[pts['id'] == d, 'r'].values[0], pts.loc[pts['id'] == d,'c'].values[0] ]]
                
        # initialize ACS
        acs = np.full_like(G_t, 999999.0)
        for r,c in origin:
            acs[r,c] = 0.0
        
        # calculated influence weighted distance transform
        acs, blx, bly = calculate_iwdt(acs, cost_dict)
        
        # create new path to destination
        path_t, track = pt.create_paths(blx, bly, origin, destination, start_path=pth_id)

        # update paths
        paths += path_t
        
        # update paths_dict
        paths_dict.update(track)
        
        # update ground potential
        G_t = ri * G_t_1 + path_t * i * (1 - (G_t_1 / G_max))
        
        # create current network cost
        d = np.full_like(G_t, 99999.0)
        d[G_t > 0] = 0.0
        d = calculate_dt(d, cellsize, option=2)
        netcost = (G_t + 1.0) * (1.0 - np.exp(d / alpha)) # G_t *
        
        #update netcost
        cost_dict['NETCOST'] = netcost
        
        # update G_t_1
        G_t_1 = np.copy(G_t)
    
    return G_t, paths, paths_dict