name = "netgen"
"""
this module offers all functions needed to generate networks for netsim  

"""

import numpy as np
import pandas as pd
import geopandas as gpd
from collections import OrderedDict
from itertools import permutations, product
from math import factorial


# constants used in netgen 

NSAMPLES = 100
#  maximum number of samples to be drawn from any group.

MAX_PERMUTATION_NUM = 7
# maximum number of points below which permulations are possible.

MIN_NUM_SAMPLE = 3
# minimum number of points above which samples can be drawn (default = max_permutation, minimum = 3, maximum = max_permutation)

MAX_ITERATIONS = 5000
# maximum total number of iterations allowed in the simulation


def setup(df):
    
    '''
    Checks and corrects input geo/pandas dataframe.
    
    Parameters
    ----------
    
    df: dataframe
        dataframe containing locations used for the simulation
        
    
    Returns
    -------
    c_df: dataframe
        correcte dataframe ready for network simulation
    
    Notes
    -----
    
    This function checks for the existence columns needed in the **netsim** simulation. Columns must have the appropriate header (as shown below).
    Columns that do not exist will be generated and populated with default values. The validity of values in existing columns is also checked. Minor
    errors and corrections are notificated. Major errors raise an exception error.
    
    The following columns need to be present:
       - *id*: exclusive identifier for each location. 
       - *group*: identifies a location as being part of a specific group. Groups can be of any size. Groups of size 1 will automatically mixed with the following
       group. A column with a single group affiliation will be created in case this column does not exist (default value is 1).
       - *seq*: identifies rank/ordering of a location within a group. There are two possible scenarios:
          - *No ordering/ ranking (default)*. Within each group a **single** value is used for all sites in that group. Depending on the number of locations in the
          group the simulation with either generate all possible permutations or *num_samples* of randomized samples (with repetition).
          - *Ordering /ranking*. Identify all sites in a group with a increasing monotonic sequence of numbers (no repetitions).
    
    '''
    
    # copy original dataframe
    c_df = df.copy()
    
    # collect basic information from dataframe
    nrows, ncols = c_df.shape
    colnames = list(c_df)
    colnames= [cn.lower() for cn in colnames]
        
    # initialize message 
    msg = []
    error_flag = False

    # id
    if 'id' not in colnames:
        # create a id column
        c_df['id'] = pd.Series(np.arange(nrows, dtype= np.int16), index=c_df.index)
    
    elif 'id' in colnames:
        # any row with the same id?
        if nrows != len(c_df['id'].unique()):
            msg.append('\n ERROR: id column - ids are not unique !!!')
            error_flag = True

    # group
    if 'group' not in colnames:
        # create group column with default (1 for single group)
        c_df['group'] = pd.Series(np.full(nrows, int(1)), index=c_df.index)
        msg.append('group column - created group column with single group !')
    
    groups= c_df['group'].unique()
    ngroup = len(groups)
        
    # seq
    if 'seq' not in colnames:
        # create a seq column with default value
        c_df['seq'] = pd.Series(np.full(nrows, int(1)), index=c_df.index)

    elif 'seq' in colnames:
        # check each group
        for g in groups:
            # Only one group?
            all_one = (c_df.loc[c_df['group'] == g, 'seq'] == 1).all()
            sequence = c_df.loc[c_df['group'] == g, 'seq'].unique()
            unique_sequence = len(c_df.loc[c_df['group'] == g]) == len(sequence)
            if not (all_one or unique_sequence): 
                error_flag = True
                msg.append('\n ERROR: seq column - sequence for group '+str(g)+' is not 1 or sequential!')
    
    # print messages
    if msg != []:
        for m in msg:
            print('\n'+m)
        
        # raise exception if any errors        
        if error_flag:
            raise('\nCheck errors !!! Network simulation ABORTED!! ')    
    else:
        print('\n No corrections or errors !! ')

        
    return c_df



def __shuffle(arr, nsamples):
    '''
    Generator function used to create randomized samples.
    
    Parameters
    ----------
    arr: 1D numpy array
        numpy array to be randomized
    
    nsamples: int
        number of randomized samples we want the generator to produce
        
    Yields
    ------
    
    arr: 1D numpy array
        randomized array

    Notes
    -----
    Internal function called by create_network_generator()

    '''
    
    for i in range(nsamples):
        np.random.shuffle(arr)
        yield tuple(arr)


def create_network_generator(df):
    '''
    Generates network generator.
    
    Parameters
    ---------
    df: dataframe
         list of locations, group membership and parameters used to network simulation
         
    Returns
    -------
    
    netgenetor: list
        network generator
    
    df_net_info: dataframe
        information about each network
    
    total_num_iter: int
        total number of iterations
    
    Notes
    -----

    This function returns three different outputs:
       - It primarily returns a network generator that results from the cartesian product of separate group generators
        (one per set of locations). A *group* generator can be of three different types (single, sample and permutation)
            - *single*: returns always the same combination of locations.
            - *sample*: returns a shuffled version of the locations. *n.b.* repetition can occur.
            - *permutation*: returns permutation of the locations (no repetition).
       - Generates a dataframe, *df_net_info*, containing generator information for each group (number of locations, total
        number of iterations, generator types.
       - total number of iterations that results from the combination of all group generators.
       
    '''

    groups = df['group'].unique()                                                           # number of groups in df?

    netgentor = list(range(len(groups)))                                                    # initialize list to store
                                                                                            # generator functions.

    net_info = {'group':[], 'num_loc':[], 'num_iter':[], 'iter_type':[]}                    # dictionary with generator
                                                                                            # summary.
    total_num_iter = 1

    for i, grp in enumerate(groups):

        df_grp = df.loc[df['group'] == grp]                                                 # select pts with the same group id.
        indxs = list(df_grp['id'])                                                          # generate a list of with pt ids in group.
        num_pts_in_grp = len(df_grp)                                                        # number of points in group?
        num_unique_pts_grp = len(df_grp['seq'].unique())                                    # number of points with distinct pt ids in group?
          
        
        if num_pts_in_grp == num_unique_pts_grp:                                            # is there a single sequence?
            net_info['group'].append(grp)                                                   # single sequence block.
            net_info['num_loc'].append(num_pts_in_grp)
            num_iter_grp = 1
            net_info['num_iter'].append(num_iter_grp)                                                 
            net_info['iter_type'].append('single')
            netgentor[i] = [tuple(indxs)]                                                   # return single seq generator.
 
        else:
            if num_pts_in_grp > MIN_NUM_SAMPLE:                                             # can it be sampled?
                if num_pts_in_grp > MAX_PERMUTATION_NUM:                                    # can it be permutated?
                    net_info['group'].append(grp)                                           # sample sequence block.
                    net_info['num_loc'].append(num_pts_in_grp)
                    num_iter_grp = NSAMPLES
                    net_info['num_iter'].append(num_iter_grp)                        
                    net_info['iter_type'].append('sample')
                    netgentor[i] = __shuffle(indxs, NSAMPLES)                               # return random seq generator.

                else:
                    if factorial(num_pts_in_grp) < NSAMPLES:                                # num_samples > permutations?
                        # permute
                        net_info['group'].append(grp)                                       # permutation sequence block.
                        net_info['num_loc'].append(num_pts_in_grp)
                        num_iter_grp = factorial(num_pts_in_grp)
                        net_info['num_iter'].append(num_iter_grp)            
                        net_info['iter_type'].append('permutation')
                        netgentor[i] = permutations(indxs)                                  # return permutation seq generator.

                    else:
                        net_info['group'].append(grp)                                       # sample sequence block.
                        net_info['num_loc'].append(num_pts_in_grp)
                        num_iter_grp = NSAMPLES
                        net_info['num_iter'].append(num_iter_grp)                        
                        net_info['iter_type'].append('sample')
                        netgentor[i] = __shuffle(indxs, NSAMPLES)                           # return random seq generator.
            else:                                                                           # too small to sample so permute
                net_info['group'].append(grp)                                               # permutation sequence block.
                net_info['num_loc'].append(num_pts_in_grp)
                num_iter_grp = factorial(num_pts_in_grp)
                net_info['num_iter'].append(num_iter_grp)            
                net_info['iter_type'].append('permutation')
                netgentor[i] = permutations(indxs)                                          # return permutation seq generator.

        total_num_iter = total_num_iter * num_iter_grp                                      # update total number of iterations.
    
    df_net_info = pd.DataFrame(net_info)
    print('\n iteration broken per group....\n')
    print(df_net_info)
    print('\n total number of iterations....',total_num_iter)
                                                                                            # return generator for network
                                                                                            # iteration, iteration info and 
    return product(*netgentor), df_net_info, total_num_iter                                 # total number of iterations

def network_layout(df, iteration, iter_num, df_net=None, twoway= False, opt='close'):
    '''
        Creates a dataframe with information for each path in network.
        
    Parameters
    ----------
    
    df: dataframe
         list of locations, group membership used to generate network
        
    net_iteration: tuple of tuples
        a tuple containing one or several tuples, one per group, representing a single network iteration
    
    iter_num: int
        iteration identifier
    
    df_net: dataframe, optional
        dataframe containing the identifiers of the origin and destination of each path plus the iteration identifier
    
    twoway: boolean
        if True two-way paths are generated for each pair of locations in network. Default: False.
    
    opt: string
        Type of network to generate. The options are as follows:
        - *close*: This option defines an independent close network of paths for each group. In this network, the locations are connected
        in order so that the first location is connected to the second one and so one until the last location is conected to the first.
        No network is defined if the first group is made of a single location.
        - *central*: This option defines a network of that consists of centralized set of paths from the locations of the first group to 
        all of the locations in the remaining groups.
        - *decentral*: This option defines a network of paths so that the locations of each group are connected to the locations of the following
        (lower level) group.
        - *distributed*: This option defines a network of paths similar to *decentral*, where the locations of each group are connected to
        the locations of the following (lower level) group, and in addition, locations within each group are interconnected.
        - *all*: This option defines a network of paths from amongst all locations.
    
    Returns
    -------
    
    df_net: dataframe, optional
      dataframe containing the identifiers of the origin and destination of each path plus the iteration identifier
        
    Notes
    -----
    
    This function takes a dataframe with locations, a list of lists containing an ordering of these locations obtained after running 
    ```create_network_generator()``` function and a iteration identifier number. It generates, or updates, a dataframe with the 
    identifiers of the origin and destination of each path that make up the path network for this specific iteration. 

    '''
    
    # option valid?
    options = ['close', 'central', 'decentral','distributed', 'all']
    if opt not in options:
        raise Exception("'opt' not valid!!! Choose from {}".format(options))
    
    # df_net?
    if df_net is None:
        net = {'origin': 'int32', 'destination': 'int32', 'iteration': 'int32'}
        df_net = pd.DataFrame(columns=list(net.keys())).astype(net)      
        
    # distinct groups
    groups = df['group'].unique()
    n_groups = len(groups)
    
    if n_groups == 1:
        if opt == 'close':
            pass
        else:
            # by default do all
            pass
    else:
        # sizes of each group
        grp_sizes = df.groupby(['group']).count().id.values
        
        if opt == 'close':
                       
            for i in range(n_groups):
                if grp_sizes[i] > 1:
                    # generate paths connecting all locations within the same group
                    origins = df.loc[df['group']== groups[i]]['id'].tolist()
                    destinations = origins[1:]+[origins[0]]
                    for o, d in zip(origins, destinations):
                        df_net.loc[len(df_net.index)] = [o,d,iter_num]
                        if twoway & (grp_sizes[i]>2):
                            df_net.loc[len(df_net.index)] = [d,o,iter_num]
                    
                else:
                    print('\nWARNING: No path for first group 1 calculated (size = 1)')                    
        
        elif opt == 'central':
            
            # first group
            if grp_sizes[0] > 1:
                # more than one location. define paths connecting all locations in group
                origins = df.loc[df['group']== groups[0]]['id'].tolist()
                destinations = origins[1:]+[origins[0]]
                for o, d in zip(origins, destinations):
                    df_net.loc[len(df_net.index)] = [o,d,iter_num]
                    if twoway:
                        df_net.loc[len(df_net.index)] = [d,o,iter_num]
            
            # set destinations to all locations in first group 
            destinations_up = df.loc[df['group']== groups[0]]['id'].tolist()
            
            # loop thru lower levels
            for i in range(1,n_groups):
                # define paths from lower level locations to all first group locations
                origins = df.loc[df['group']== groups[i]]['id'].tolist()
                for o in origins:
                    for d in destinations_up:
                        df_net.loc[len(df_net.index)] = [o,d,iter_num]
                        if twoway:
                            df_net.loc[len(df_net.index)] = [d,o,iter_num]

        elif opt == 'decentral':

            # first group
            if grp_sizes[0] > 1:
                # more than one location, define paths connecting all locations in group
                origins = df.loc[df['group']== groups[0]]['id'].tolist()
                destinations = origins[1:]+[origins[0]]
                for o, d in zip(origins, destinations):
                    df_net.loc[len(df_net.index)] = [o,d,iter_num]
                    if twoway:
                        df_net.loc[len(df_net.index)] = [d,o,iter_num]
            
            # set destinations to all locations in first group 
            destinations_up = df.loc[df['group']== groups[0]]['id'].tolist()
            
            # loop thru levels
            for i in range(1,n_groups):
                
                # define paths to from higher level locations to lower level locations
                origins = df.loc[df['group']== groups[i]]['id'].tolist()
                for o in origins:
                    for d in destinations_up:
                        df_net.loc[len(df_net.index)] = [o,d,iter_num]
                        if twoway:
                            df_net.loc[len(df_net.index)] = [d,o,iter_num]
                
                # set destinations to all locations in the previous (higher level) group
                destinations_up = origins

        elif opt == 'distributed':
        
            # first group
            if grp_sizes[0] > 1:
                # more than one location
                origins = df.loc[df['group']== groups[0]]['id'].tolist()
                destinations = origins[1:]+[origins[0]]
                for o, d in zip(origins, destinations):
                    df_net.loc[len(df_net.index)] = [o,d,iter_num]
                    if twoway:
                        df_net.loc[len(df_net.index)] = [d,o,iter_num]

            # set destinations to all locations in first group
            destinations_up = df.loc[df['group']== groups[0]]['id'].tolist()

            # loop thru levels
            for i in range(1,n_groups):

                # define paths to from higher level locations to lower level locations
                origins = df.loc[df['group']== groups[i]]['id'].tolist()
                for o in origins:
                    for d in destinations_up:
                        df_net.loc[len(df_net.index)] = [o,d,iter_num]
                        if twoway:
                            df_net.loc[len(df_net.index)] = [d,o,iter_num]

                # paths between locations of the same level
                indx = range(len(origins))
                for i in indx[:-1]:
                    for j in indx[i+1:]:
                        df_net.loc[len(df_net.index)] = [origins[i], origins[j], iter_num]
                        if twoway:
                            df_net.loc[len(df_net.index)] = [origins[j], origins[i], iter_num]

                # update lower_destinations with previous origins
                destinations_up = origins
            
        else: #opt == 'all'#
            
            # define paths from all locations to all locations
            origins = df['id'].tolist()
            indx = range(len(origins))
            for i in indx[:-1]:
                for j in indx[i+1:]:
                    df_net.loc[len(df_net.index)] = [origins[i], origins[j], iter_num]
                    if twoway:
                        df_net.loc[len(df_net.index)] = [origins[j], origins[i], iter_num]                   

    return df_net