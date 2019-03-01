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
       - *mix*: indicates whether a location from a group can mix with a member of the previous group. This value (0 or 1) must be the same for all members in a group
       (default value is 1). *n.b.* the mix parameter for all locations in the first group is always set to 0
       - *residuality index (ri)*. This correspongs to a fraction $\ge$ 0 and $\le$ 1.0 (default is $\frac{1}{group\_ size}$ for all sites in a group). *n.b.* it 
       must be the same for all members in a group. 
    
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
            msg.append('ERROR: id column - ids are not unique !!!')
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
                msg.append('ERROR: seq column - sequence for group '+str(g)+' is not 1 or sequential!')
    
    # mix
    if 'mix' not in colnames:
        # create mix column with default value (1 for mixing)
        c_df['mix'] = pd.Series(np.full(nrows, int(1)), index=c_df.index)

        # set mix of first group to 0
        c_df.loc[c_df['group'] == groups[0], 'mix'] = 0

    elif 'mix' in colnames:
        # verify for valid codes
        valid = c_df['mix'].isin([0,1]).all()    

        if not valid:
            error_flag = True
            msg.append('ERROR: mix column - mix codes are not 0 or 1 !!')
        else:
           # group_mix for the first group always set to 0 
            valid = (c_df.loc[c_df['group'] == groups[0], 'mix'] == 0).all()
            if not valid:
                c_df.loc[c_df['group'] == groups[0], 'mix'] = 0
                msg.append('mix column - resetting first group mix to 0')

            # set group_mix to the majority mix value in group
            for g in groups:
                unique_mix = c_df.loc[c_df['group'] == g, 'mix'].unique()
                if len(unique_mix) > 1:
                    sel = c_df['group'] == g
                    majority_grp_mix = c_df[sel]['mix'].value_counts().idxmax()
                    c_df.loc[sel,'mix'] = majority_grp_mix
                    msg.append('mix column - resetting mix in group '+str(g)+
                               ' to majority value '+str(majority_grp_mix))
    
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

        
def __paths_table(ids, iteration, df_out, option='cummulative'):
    '''
    
    Generates each path for in a network.
    
    Parameters
    ----------
    
    ids: list
        an iteration of location ids for a group 
   
    iteration: int
        identifier of current iteration
    
    df_out: pandas dataframe
        dataframe used to store the origin, destination and iteration id for each path in a network associated
        with a particular iteration of ids. The dataframe may be empty or not.
    
    option: string
        defines the type of network that will be generated. Currently the options are as follows:
           - *close*: this network is generated by defining a path between each location id and the next location
           id in the sequence that defined by the current ids sequence. A final path is defined between the last 
           location id and the first location id in ids.
           - *cummulative* (*default*) : this network is generated by defining a first path between the first two
           location ids in ids. The remaining location ids are used as origins of the subsequent paths such that 
           destinations are *all* location ids *previous* to the current origin id.
           - *all*: this network is generated by considering every location id in ids as an origin and the remaining
           ids as destinations. 
    
    Returns
    -------
    
    df_out: pandas dataframe
        updated version of dataframe used to store the origin, destination and iteration id for each path in a network
        
    Notes
    -----
        
        'Private' function called by network layout.
        
    '''    
    
    if option == 'cummulative':
        index = range(len(ids)) 
        # first path
        df_out.loc[len(df_out.index)] = [ids[index[0]], ids[index[1]], iteration]

        # loop thru remaining group ids
        for o in index[2:]:
            for d in index[:o]:
                df_out.loc[len(df_out.index)] = [ids[index[o]], ids[index[d]], iteration]
    
    elif option == 'all':
        for o_indx in range(len(ids)):
            d_indxs = ids[:o_indx]+ids[o_indx+1:]
            for d_indx in d_indxs:
                df_out.loc[len(df_out.index)] = [ids[o_indx], ids[d_indx], iteration]        
    
    elif option == 'close':
        for o, d in zip(ids[:-1], ids[1:]):
            df_out.loc[len(df_out.index)] = [o, d, iteration]
        
        # append last path
        df_out.loc[len(df_out.index)] = [ids[-1], ids[0], iteration]
        
    return df_out

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
       - It primarily returns a network generator that results from the (cartesian) product of separate group generators
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


def network_layout(df, net_iteration, iteration, df_out_net=None, option='expanded'):
    '''
    Creates a table with information for each path in  network.
        
    Parameters
    ----------
    
    df: dataframe
         list of locations, group membership and parameters used to generate network
        
    net_iteration: list of list
        a list containing one or several lists, one for each group network that resulted from an single network iteration
    
    iteration_id: int
        iteration identifier
    
    df_out_net: dataframe, optional
        table containing the identifiers of the origin and destination of each path plus the iteration identifier
    
    Returns
    -------
    
    df_out_net: dataframe, optional
       
    Notes
    -----
    
    This function takes a dataframe with locations, a list of lists containing an ordering of these locations obtained after *running 
    create_network_generator()* function and a iteration identifier number. It generates, or updates, a dataframe with the 
    identifiers of the origin and destination of each path that make up the path network for this specific iteration. 
    
    *Mixing*. This refers to whether only paths between locations within the same group are allowed or not. If mixing occurs,
    a path will be generated that connects locations that belong to two distinct groups. More specifically, between a the first
    location in a group and the last location in the *previous* group. The following observations need to be kept in mind with
    regards to mixing:
    - All pts in a group must have the same 'mix' code.
    - The value of the 'mix' column for all members of the first group is always set to 0 (i.e. no mixing).
    - We identify two degenerative cases if the number of points in a group is less than 2:
      - if the group is the first one to be processed the point will be set as an origin and will be mixed with the points of the
      following group.
      - if not the first group, it will be set as a destination and will be mixed with the points of the previous group.
      
    '''

    if df_out_net is None:
        df_out_net = pd.DataFrame(columns=['origin', 'destination', 'iteration'])
        df_out_net['iteration'] = pd.to_numeric(df_out_net['iteration'], downcast='integer')
        
    prev_pt_id = 0                                                                      # indentifies last point used in
                                                                                        # previous group. Only used for
                                                                                        # path between points in different
                                                                                        # groups.

    groups = df['group'].unique()                                                       # distinct groups in df

    for pos, g in enumerate(groups):                                                    # Loop thru groups

        ids = list(net_iteration[pos])                                                  # create list with pts id in group.
        num_pts_in_grp = len(ids)                                                       # number of pts in group.

        # degenerative cases
        if (num_pts_in_grp == 1) and (pos == 0):                                        # one point in first group only.

            prev_pt_id = ids[0]                                                         # set point to origin of next path.
            sel_next = df['group'] == groups[1]                                         # force second group into mixing
            df.at[sel_next, 'mix'] = 1

        elif (num_pts_in_grp == 1) and (pos > 0):                                       # group with only one point.
            
            prev_pt_id = ids[0]                                                         # set pt to next previous_pt_id.
            ids.insert(0, prev_pt_id)                                                   # previous pt id is path origin.
            df_out_net = __paths_table(df, ids, iteration, df_out_net, option)
            
        # non-degenerative cases
        else:

            if df.loc[df['group'] == g, 'mix'].iloc[0] == 0:                            # no mixing.
                
                prev_pt_id = ids[-1]
                df_out_net = __paths_table(df, ids, iteration, df_out_net, option)

            else:                                                                       # mixing.
                ids.insert(0, prev_pt_id)
                prev_pt_id = ids[-1]
                df_out_net = __paths_table(df, ids, iteration, df_out_net, option)

    return df_out_net