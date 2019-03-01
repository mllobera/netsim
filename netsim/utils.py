'''

'''

from __future__ import with_statement
import rasterio
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
import geopandas as gpd


def read_raster(fn):
    '''
    Reads raster
    
    Parameters
    ----------
    
    fn: string
        path to raster image (assume geotiff)
    
    Returns
    -------
    
    ras: 2D numpy array
        raster
    
    meta: dictionary
        raster geospatial information
    
    '''
    
    try:       
        with rasterio.open(fn) as src:
            meta   = src.meta
            ras    = src.read(1)
            ras = np.ascontiguousarray(ras, dtype=np.float64)
            meta['bounds'] = src.bounds
            return ras, meta
    
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

def rc2pt(rc, meta):
    '''
    Returns the center x, y coordinates of cells at row, column locations.
    
    Parameters
    ----------
    
    rc: list
        list of rows and column pairs
    
    meta: dictionary
        raster geospatial information
    
    Returns
    -------
    pts: list
        list of tuples representing (x,y) coordinates
    
    No checks.
    '''
    
    r, c = rc[0,:], rc[1,:]
    
    # collecting information
    xul, yll, _ , _ = meta['bounds']
    nrows = meta['height']
    cellsize = meta['affine'].a
         
    # calculating x and y
    ys = yll + ((nrows - 1) - r) * cellsize + cellsize / 2
    xs = xul + c * cellsize + cellsize / 2
    
    pts = [(x, y) for x,y in zip(xs, ys)]
    
    return pts

def pt2rc (pts, meta):
    '''
    Returns row and column
    
    Parameters
    ----------
    
    pts: shapely point
        points
    
    meta: dictionary
        raster geospatial information
        
    Returns
    -------
    
    r,c: numpy arrays
        rows and columns
    
    Notes
    -----
    Expects the point column of a geopandas dataframe and returns the rows and columns for each point. No checks on
    whether points are within bounding box.
    '''
    
    # collecting information
    xul, _, _, yul = meta['bounds']
    nrows, ncols = meta['width'], meta['height']
    cellsize = meta['affine'].a
    
    # convert easting and northings to rows and columns
    r = (yul - pts.y.values) // cellsize
    r[r == nrows] -= 1
    c = (pts.x.values - xul) // cellsize
    c[c == ncols] -= 1
    
    return r.astype('int32'), c.astype('int32')

def plot_map(ras, loc= None, save= False, **kwsplot):
    '''
    Plots a raster.
    
    Parameters
    ----------
    
    ras: 2D numpy array
        raster to be displayed
    
    loc: geodataframe
        location information
    
    save: boolean
        if True save output (default is False)
    
    kwsplot: dictionary
        dictionary that contains other dictionaries (kwsras, kwslbl, kwsother)
    
    Returns
    -------
    None
    
    Notes
    -----
    talk about the internal dictionaries
    
    '''
    
    # separate dictionaries
    if kwsplot:
        if 'kwsras' in kwsplot.keys():
            kws_ras = kwsplot['kwsras']
        else:
            kws_ras ={}

        if 'kwslbl' in kwsplot.keys():
            kws_lbl = kwsplot['kwslbl']
        else:
            kws_lbl ={}
            
        if 'kwsother' in kwsplot.keys():
            kws_other = kwsplot['kwsother']
        else:
            kws_other={}

    if 'label'  in kws_other.keys():
        label = kws_other['label']
    else:
        label = 'id'
        
    if 'cbar_pct' in kws_other.keys():    
        cbar = True
        cbar_pct = kws_other['cbar_pct']
    else:
        cbar= False   

    if 'figsize' in kws_other.keys():
        figsize = kws_other['figsize']
    else:
        figsize = (6,6)

    if 'title' in kws_other.keys():
        title = kws_other['title']
    else:
        title = ''

    if 'output' in kws_other.keys():
        output = kws_other['output']
    else:
        output= './map.png'
        

    if 'cmap' in kws_ras:
        kws_ras['cmap'] = mpl.cm.get_cmap(kws_ras['cmap'])

    # colorbar?
    if cbar:        
        wcbar = cbar_pct/100.0
        gridspec_kw={'width_ratios':[1,wcbar], 'wspace': 0.05}
        fig, (ax, cax) = plt.subplots(ncols = 2, figsize= (figsize[1] + figsize[1]*wcbar, figsize[0]), gridspec_kw=gridspec_kw)
    else:
        fig, ax = plt.subplots(ncols = 1, figsize= figsize)

    if 'bground' in kws_other.keys():
        im2 = ax.imshow(kws_other['bground'], extent=kws_ras['extent'], origin='upper', cmap= mpl.cm.get_cmap('Greys'))

        # plot map
    im = ax.imshow(ras, origin= 'upper', **kws_ras)
    
        
    ax.set_title(title) 

    # coloram
    if cbar:
        fig.colorbar(im, cax= cax)

    # locations?
    if loc is not None:
        xs, ys = loc['geometry'].x.values, loc['geometry'].y.values
        data = zip(loc[label].values, xs, ys )
        ax.scatter(xs,ys, color=kws_lbl['color'])        
        for id, x, y in data:
            ax.annotate(str(id), xy=(x,y), **kws_lbl)

    # save?
    if save:
        plt.savefig(fname= output, dpi= 200)

    plt.show()

def calculate_hillshade(img ,azimuth,angle_altitude):
    '''
    Calculates hillshade
    
    Parameters
    ----------
    
    img: 2D numpy array
        heightmap
    
    azimuth: float
        Horizontal direction of the source of light (degrees)
    
    angle_altitude: float
        elevation angle of the source of light (degrees)
    
    '''
    azimuth = 360.0 - azimuth    
    x, y = np.gradient(img)
    slope = np.pi/2. - np.arctan(np.sqrt(x*x + y*y))
    aspect = np.arctan2(-x, y)
    azimuthrad = np.radians(azimuth)
    altituderad = np.radians(angle_altitude)
 
    shaded = np.sin(altituderad)*np.sin(slope) + np.cos(altituderad)*np.cos(slope)*np.cos((azimuthrad - np.pi/2.) - aspect)
    
    return 255*(shaded + 1)/2