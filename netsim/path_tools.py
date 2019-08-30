'''
This module contains functions used to generate paths from backlink files
'''

import numpy as np

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
        row, column of origin    
    
    destinations: list of destinations
        list of destinations [[row, colum]]
    
    Returns
    -------
    paths: 2D numpy array
        array that results from adding all paths
    
    path_lst: list of dictionaries
        a list of paths each represented by a dictionary containing destination, origin and track 
    
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