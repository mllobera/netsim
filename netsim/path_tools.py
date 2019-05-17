'''
Path_tools module with functions to derive paths
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


def create_paths(BLX, BLY, origin, destinations, start_path=0):
    '''
    Creates a path to each destination.
    
    Parameters
    ----------
    
    BLX: 2D numpy array
        horizontal backlink 
    
    BLY: 2D numpy array
        vertical backlink
    
    origin: list
        row, column of origin    
    
    destinations: list of destinations
        list of destinations [[row, colum]]
    
    Returns
    -------
    paths: 2D numpy array
        paths to all destination added
    
    paths_dict: dictionary
        contains destination, origin and track 
    
    Notes
    -----
    
    Depending on the size of the chamfer window used, backlink arrays may contain jumps that are greater than one cell.
    Hence we need to use an auxiliary function, segment(), to connect these.
        
    '''

    # array to store path/s
    paths = np.zeros_like(BLX, dtype=np.int16)
    
    # path number
    path_num= start_path
    
    # initialize network paths dictionary
    path_lst = [] #paths_dict = {} #OrderedDict()
    
    for destination in destinations:
        
        # array to store current path
        pth = np.zeros_like(BLX, dtype=np.uint16)
        
        # row & columns of current path
        rcs = []

        # initialize first cell location
        r0, c0 = destination

        while (BLX[r0, c0] != 0) or (BLY[r0, c0] != 0):

            # update to new path location
            r1 = r0 + BLX[r0, c0]
            c1 = c0 + BLY[r0, c0]

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