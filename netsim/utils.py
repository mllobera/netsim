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

def plot_map(raster, loc= None, title= None, figsize= (5,5), cmap= 'viridis', cbar= False, save= None, **kwother):
    '''
    Basic raster plot.
    
    Parameters
    ----------
    
    ras: dictionary
        dictionary must have at least two entries: 'ras', 2D numpy array representing a raster; 'meta', meta
        information (from *rasterio*) about the raster. Optional entries are 'bground', a 2D numpy array
        representing a background raster and 'paths', representing a network of paths (it assumed that 'ras'
        represents a DEM).
    
    loc: dictionary or geodataframe
        used to identify point locations. if dictionary then it must have at leat two entries: 'df', identifying
        geopandas framework holding point data, 'label', name of the column in 'df' used to labelling points.
    
    title: string
        if not empty then title to be used when displaying ras
    
    figsize: tuple
        Size of figure. *Default*: (5,5)
    
    cmap: string
        name of the matplotlib colormap. *Default*: 'viridis'
    
    cbar: boolean
        if True colorbar is displayed. *Default*: False
    
    save: string
        if not empty then name of the output image (default is None)
        
    Returns
    -------
    None
        
    '''
   

    # set figure
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
        if set(['ras', 'meta']) <= set(raster.keys()):
            
            # default
            img = raster['ras']
            alpha = 1.0

            # calculate extension
            bounds = raster['meta']['bounds']
            extent = [bounds.left, bounds.right, bounds.bottom, bounds.top]
            
            # background image?
            if 'bground' in raster.keys():
                alpha = 0.5
                im2 = ax.imshow(raster['bground'], extent= extent, origin='upper', cmap= mpl.cm.get_cmap('Greys'))

            # paths overlay?    
            if 'paths' in raster.keys():
                ax.contourf(img, 6, extent=extent, origin= 'upper', cmap= mpl.cm.get_cmap('Purples'), alpha=0.4)
                img = np.ma.array(raster['paths'], mask= raster['paths'] == 0.0)
                
            # plot main raster
            im = ax.imshow(img, origin= 'upper', cmap= cmap, extent= extent, alpha= alpha)
            
            # colorbar?
            if cbar:
                fig.colorbar(im, cax= cax)
        else:
            raise Exception('raster is missing  ''"ras"'' and/or ''"meta"'' keys! ')
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
                ax.scatter(xs,ys, color='lightgray')        
                for id, x, y in data:
                    ax.annotate(str(id), xy=(x,y), color='lightgray', xytext= (2.5,2.5), textcoords='offset points')                
            else:
                 raise Exception('Loc does not have the right keys!')
        else:
             if isinstance(loc, gpd.geodataframe.GeoDataFrame):
                    xs, ys = loc['geometry'].x.values, loc['geometry'].y.values
                    data = zip(xs, ys)
                    ax.scatter(xs,ys, color='lightgray')        
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
        heightmap
    
    az: float
        Horizontal direction of the source of light (degrees)
    
    elev_angle: float
        elevation angle of the source of light (degrees)
    
    '''
    az = 360.0 - az    
    x, y = np.gradient(img)
    slope = np.pi/2. - np.arctan(np.sqrt(x*x + y*y))
    aspect = np.arctan2(-x, y)
    azrad = np.radians(az)
    altituderad = np.radians(elev_angle)
 
    shaded = np.sin(altituderad)*np.sin(slope) + np.cos(altituderad)*np.cos(slope)*np.cos((azrad - np.pi/2.) - aspect)
    
    return 255*(shaded + 1)/2