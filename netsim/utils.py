'''
This module contains various i/o and utility functions
'''

from __future__ import with_statement
import rasterio as ro
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
import geopandas as gpd
import networkx as nx


def read_raster(fn):
    '''
    Reads raster into a 2D numpy array.
    
    Parameters
    ----------
    
    fn: string
        path to raster image (assume geotiff)
    
    Returns
    -------
    
    ras: 2D numpy array
        raster
    
    profile: dictionary
        raster geospatial information
    
    '''
    
    try:       
        with ro.open(fn) as src:
            profile   = src.profile
            ras    = src.read(1)
            # change nodata value 
            ras[ras == profile['nodata']] = -9999
            
            # make raster c-contiguous
            ras = np.ascontiguousarray(ras, dtype=np.float64)
            profile['dtype']= 'float64'
            profile['bounds'] = src.bounds
            return ras, profile
    
    except EnvironmentError:
        print('Oops! Could not find file')
        
        
def add_polyline(track, gdf):
    '''
    Updates a geopandas dataframe with a polyline.
    
    Parameters
    ----------
    
    track: list of tuples
        list of coordinates representing a polyline
    
    gdf: geodataframe
    
    Returns
    _______
    
    gdf: geodataframe
        Updated geodataframe
    
    '''
    
    from shapely.geometry import LineString
   
    gdf.loc[len(gdf.index), 'geometry'] = LineString(coordinates)
    
    return gdf

def rc2pt(rc, profile):
    '''
    Returns the center x, y coordinates of cells at (row, column) locations.
    
    Parameters
    ----------
    
    rc: list
        list of rows and column pairs
    
    profile: dictionary
        raster geospatial information
    
    Returns
    -------
    pts: list
        list of tuples representing (x,y) coordinates
  
    '''
    
    r, c = rc[0,:], rc[1,:]
    
    # collecting information
    xul, yll, _ , _ = profile['bounds']
    nrows = profile['height']
    cellsize = profile['transform'].a
         
    # calculating x and y
    ys = yll + ((nrows - 1) - r) * cellsize + cellsize / 2
    xs = xul + c * cellsize + cellsize / 2
    
    pts = [(x, y) for x,y in zip(xs, ys)]
    
    return pts

def pt2rc (pts, profile):
    '''
    Returns row and column.
    
    Parameters
    ----------
    
    pts: *shapely* point
        points
    
    profile: dictionary
        raster geospatial information
        
    Returns
    -------
    
    r,c: tuple of numpy arrays
        rows and columns
    
    Notes
    -----
    Expects the point column of a geopandas dataframe and returns the rows and columns for each point. No checks on
    whether points are within bounding box.
    '''
    
    # collecting information
    xul, _, _, yul = profile['bounds']
    nrows, ncols = profile['height'], profile['width']
    cellsize = profile['transform'].a
    
    # convert easting and northings to rows and columns
    r = (yul - pts.y.values) // cellsize
    r[r == nrows] -= 1
    c = (pts.x.values - xul) // cellsize
    c[c == ncols] -= 1
    
    return r.astype('int32'), c.astype('int32')

def plot_map(raster, loc= None, title= None, figsize= (5,5), cmap= 'viridis', cbar= False, save= None, **kwother):
    '''
    Basic raster plot.
    
    Parameters
    ----------
    
    raster: dictionary
        dictionary must have at least two entries:
        
        - **'ras'**: 2D numpy array representing a raster
        - **'profile': raster information (as derived from *rasterio*)
        - *Optional* entries are:
          
          - *'bground'*: a 2D numpy array representing a background raster. Typically a hillshade.
          - *'paths'*:  dataframe representing a network of paths as obtained from ``calculate_paths()`` 
    
    loc: dictionary or geodataframe, optional 
        used to identify point locations. if *dictionary* then it must have at least two entries:
        
        - **'df'**:  geopandas framework holding point data.
        - **'label'**: name of the column in 'df' used to label points.
    
    title: string, optional
        if not empty then title to be used when displaying ras
    
    figsize: tuple
        Size of figure. *Default*: (5,5)
    
    cmap: string
        name of the matplotlib colormap. *Default:* 'viridis'
    
    cbar: boolean
        if True colorbar is displayed. *Default:* False
    
    save: string, optional
        if not empty then name of the output image. *Default: None*
         
    '''
   
    # set colormap
    if cmap:
       cmap = mpl.cm.get_cmap(cmap)
    
    if cbar: # colorbar?       
        wcbar = 0.03
        gridspec_kw={'width_ratios':[1,wcbar], 'wspace': 0.05}
        fig, (ax, cax) = plt.subplots(ncols = 2, figsize= (figsize[1] + figsize[1]*wcbar, figsize[0]), gridspec_kw=gridspec_kw)
    else:
        fig, ax = plt.subplots(ncols = 1, figsize= figsize)

    # raster
    if isinstance(raster, dict):
        if set(['ras', 'profile']) <= set(raster.keys()):
            
            # default
            img = raster['ras']
            profile = raster['profile']
            alpha = 1.0

            # calculate extension
            bounds = raster['profile']['bounds']
            extent = [bounds.left, bounds.right, bounds.bottom, bounds.top]
            
            # background image?
            if 'bground' in raster.keys():
                alpha = 0.5
                im0 = ax.imshow(raster['bground'], extent= extent, origin='upper', cmap= mpl.cm.get_cmap('Greys'),**kwother)

            # nodata?
            if np.any(img == profile['nodata']):
                nodata_mask = img == profile['nodata']
                img = np.ma.array(img, mask = nodata_mask)
                cmap.set_bad('white',0.0)
                
            # paths overlay?    
            if 'paths' in raster.keys():
                
                # assuming img is a dem
                cmap.set_bad('white',0.0)
                ax.contourf(img, 6, extent=extent, origin= 'upper', cmap= cmap, alpha=0.5)
                
                # path_mask = no paths
                paths_mask = raster['paths'] == 0.0

                # combine with existing mask?
                if np.ma.is_masked(img):
                    combined_mask = nodata_mask * paths_mask
                else:
                    combined_mask = paths_mask
                
                # apply mask to path 
                img_paths = np.ma.array(raster['paths'], mask= combined_mask)
                
                # plot paths
                cmap_path = mpl.cm.get_cmap('Oranges_r')
                cmap_path.set_bad('white',0.0)
                im1 = ax.imshow(img_paths, origin= 'upper', cmap= cmap_path, extent= extent, alpha= alpha, **kwother)
                
            else: # regular plot
                im1 = ax.imshow(img, origin= 'upper', cmap= cmap, extent= extent, alpha= alpha, **kwother)
            
            # colorbar?
            if cbar:
                fig.colorbar(im1, cax= cax)
        else:
            raise Exception('raster is missing  ''"ras"'' and/or ''"profile"'' keys! ')
    else:
        raise Exception('raster must be a dictionary')

    
    if title:
        ax.set_title(title)
        
    # locations?
    if not loc is None:
        # is it a dictionary?
        if isinstance(loc, dict):
            if set(['df', 'label']) <= set(loc.keys()):
                xs, ys = loc['df']['geometry'].x.values, loc['df']['geometry'].y.values
                labels = loc['df'][loc['label']].values
                data = zip(labels, xs, ys )  
                ax.scatter(xs,ys, color='darkgray')        
                for id, x, y in data:
                    ax.annotate(str(id), xy=(x,y), color='darkgray', xytext= (2.5,2.5), textcoords='offset points')                
            else:
                 raise Exception('Loc does not have the right keys!')
        else:
             if isinstance(loc, gpd.geodataframe.GeoDataFrame):
                    xs, ys = loc['geometry'].x.values, loc['geometry'].y.values
                    data = zip(xs, ys)
                    ax.scatter(xs,ys, color='darkgray')        
             else:
                raise Exception('Not a dictionary or dataframe!!')

    # saving output?
    if save:
        output= save+'.png'      
        plt.savefig(fname= output, dpi= 300)

    plt.show()

def calculate_hillshade(img, az= 135, elev_angle= 40):
    '''
    Calculates hillshade
    
    Parameters
    ----------
    
    img: 2D numpy array
        array with elevation values
    
    az: float
        Horizontal direction of the source of light (degrees). *Default:* 135 degrees
    
    elev_angle: float
        elevation angle of the source of light (degrees). *Default:* 40 degrees
    
    '''
    az = 360.0 - az    
    x, y = np.gradient(img)
    slope = np.pi/2.0 - np.arctan(np.sqrt(x*x + y*y))
    aspect = np.arctan2(-x, y)
    azrad = np.radians(az)
    altituderad = np.radians(elev_angle)
 
    shaded = np.sin(altituderad)*np.sin(slope) + np.cos(altituderad)*np.cos(slope)*np.cos((azrad - np.pi/2.) - aspect)
    
    return 255*(shaded + 1)/2

def plot_network(df, save = None):
    '''
    Plots network
    
    Parameters
    ----------

    df : dataframe
         contains the ids of origin and destination of each path in the network
    
    save: string
        filename ('.png' added). *Default: None*

    Notes
    -----
        The dataframe above is the output from running ``generate.network_layout()``. To save output user
        must supply *filename*. An image file with '.png' extension will be saved to local directory.

    '''

    # create graph
    G = nx.Graph()
    
    # find and add nodes
    tmp = set(df.loc[:,'origin'].unique().tolist()) | set(df.loc[:,'destination'].unique().tolist())
    G.add_nodes_from(list(tmp))
    
    # add edges
    for indx, row in df.iterrows():
        G.add_edge(row['origin'],row['destination'])
    
    # draw network
    nx.draw(G, with_labels=True)

    # saving output?
    if save:
        output= save+'.png'      
        plt.savefig(fname= output, dpi= 300)
