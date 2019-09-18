'''
This module contains functions used to generate paths and explore them.
'''

import numpy as np
import pandas as pd
import geopandas as gpd

def __segment(r0, c0, r1, c1, path, rcs):
    '''
    Creates a straight segment between two locations using Bresenham algorithm.
    
    Parameters
    ----------
    
    r0,c0: ints
        row and column of origin
    
    r1,c1: ints
        row and column of end
    
    path: 2D numpy array
        current path
    
    rcs: list
        rows and columns of current path        
    
    Notes
    -----
    This is an internal function. First location is included in the segment but not the last one.
    
    '''
  
    # find row and column differences
    delta_c = abs(c1 - c0)
    delta_r = abs(r1 - r0)

    # adjust signs depending on target's quadrant
    if c0 < c1:
        sign_c = 1
    else:
        sign_c = -1

    if r0 < r1:
        sign_r = 1
    else:
        sign_r = -1
    e = delta_c - delta_r

    r, c = r0, c0

    while (r != r1) or (c != c1):
        path[r, c] = 1
        rcs = rcs + [[r,c]]

        e2 = 2 * e
        if e2 > -delta_r:
            e -= delta_r
            c += sign_c
        if e2 < delta_c:
            e += delta_c
            r += sign_r
            
    return path, rcs


def create_paths(blx, bly, origin, destinations, start_path=0):
    '''
    Creates a path to each destination.
    
    Parameters
    ----------
    
    blx: 2D numpy array
        horizontal backlink 
    
    bly: 2D numpy array
        vertical backlink
    
    origin: list
        list of origins [[row, colum],...]    
    
    destinations: list
        list of destinations [[row, colum],...]

    start_path: int
        path identifier, optional
    
    Returns
    -------
    paths: 2D numpy array
        array that results from adding all paths
    
    path_lst: list
        a list of paths. Each path is represented by a dictionary containing three entries:
        
        - 'destination': [row,col] of destination
        - 'origin': [row, col] of origin
        - 'track' : list containing two lists: [rows], [cols] for each cell making up the path 
    
    Notes
    -----
    
    Depending on the size of the chamfer window used, backlink arrays may contain jumps that are greater than one cell.
        
    '''

    # array to store path/s
    paths = np.zeros_like(blx, dtype=np.int16)
    
    # path number
    path_num= start_path
    
    # several destinations?
    if len(destinations) > 1:

        # initialize network paths dictionary
        path_lst = []

        for destination in destinations:
            
            # array to store current path
            pth = np.zeros_like(blx, dtype=np.uint16)
            
            # row & columns of current path
            rcs = []

            # initialize first cell location
            r0, c0 = destination

            while (blx[r0, c0] != 0) or (bly[r0, c0] != 0):

                # update to new path location
                r1 = r0 + blx[r0, c0]
                c1 = c0 + bly[r0, c0]

                # add segment to new path location
                pth, rcs = __segment(r0, c0, r1, c1, pth, rcs)

                # update current location
                r0, c0 = r1, c1

            # add very last location
            pth[r0, c0] = 1
            rcs = rcs + [[r0, c0]]
            
            # store path information
            path_num += 1
            path_lst.append({'id': path_num,
                            'origin': origin[0],
                            'destination': destination,
                            'track': np.array(rcs).T})
            
            # add new path
            paths += pth
        
        return paths, path_lst

    else: # only one destination

        # array to store current path
        pth = np.zeros_like(blx, dtype=np.uint16)
        
        # row & columns of current path
        rcs = []

        # initialize first cell location
        r0, c0 = destinations[0]

        while (blx[r0, c0] != 0) or (bly[r0, c0] != 0):

            # update to new path location
            r1 = r0 + blx[r0, c0]
            c1 = c0 + bly[r0, c0]

            # add segment to new path location
            pth, rcs = __segment(r0, c0, r1, c1, pth, rcs)

            # update current location
            r0, c0 = r1, c1

        # add very last location
        pth[r0, c0] = 1
        rcs = rcs + [[r0, c0]]
        
        # store path information
        path_num += 1
        path_dict = {'id': path_num,
                        'origin': origin[0],
                        'destination': destinations[0],
                        'track': np.array(rcs).T}
        
        # add new path
        paths += pth

        return paths, path_dict


def path_stats(df_paths, ras, df, fun_dic={'fun':np.sum, 'name':'sum'}):
    '''
    Applies ``<function>`` to values in *ras* along each path in *df_paths* dataframe
    
    Parameters
    ----------
    
    df_paths: dataframe
        contains information for various paths
    
    ras: 2D numpy array
        raster from where values are going to be extracted
    
    df: dataframe
        original dataframe with location information
    
    fun_dic: dictionary 
       a dictionary with two entries:

       - **'fun'**:  a numpy function to compute values along a path track
       - **'name'**: function name
       
    Returns
    -------
    
    df_paths: dataframe
        updated version of *df_paths* with an additional *name* column containing the results obtained 
        after applying *<function>* on *ras* values along each path.
    
    Notes
    -----
    
    **df_paths** dataframe must contain a column with the path track (a 2D numpy array with the row and columns
    that make up a path)
    
    '''
    # unpack fun
    f= fun_dic['fun']
    name = fun_dic['name']
    
    # initialize variables
    path_ids = []
    path_stats = []  
    i=0   
    
    for _,pth in df_paths.iterrows():
        
        # extract current path values
        path_values = ras[pth['track'][0], pth['track'][1]]
        
        # find origin and destination ids
        sel = (df['r'] == pth['origin'][0]) & (df['c'] == pth['origin'][1])
        o = df.loc[sel]['id'].values[0]
        sel = (df['r'] == pth['destination'][0]) & (df['c'] == pth['destination'][1])
        d = df.loc[sel]['id'].values[0] 
        
        # generate statistic
        path_ids += [(o,d)]
        path_stats += [f(path_values)]
    
    # Update with new information
    df_paths['path_ids'] = path_ids
    df_paths[name] = path_stats
    
    return df_paths[['id', 'path_ids', 'origin', 'destination', 'track', name]]