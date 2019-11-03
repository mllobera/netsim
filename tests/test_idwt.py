import geopandas as gpd
import numpy as np
import netsim.utils as utils
import netsim.path_tools as ptools
from netsim.cost import calculate_dt, calculate_iwdt
from pathlib import Path

data_path = Path.cwd().parent / "netsim" / "data"

fn_dem = data_path / "sample" / "sampleDEM.tif"

dem, profile = utils.read_raster(fn_dem)

cellsize = profile['transform'].a

fn_shp = data_path / "sample" / "sample5.shp"

df_temp = gpd.read_file(fn_shp)

# make a copy
df = df_temp.copy(deep=True)
df['r'], df['c'] = utils.pt2rc(df['geometry'], profile)

vftfn = data_path / "iwdt" / "grad2cost.csv"
vft = np.genfromtxt(vftfn, delimiter=',')

coef = np.polyfit(np.tan(np.radians(vft[: , 0])), vft[:,1], deg=4)

iwdt = np.full_like(dem, 999999.0)

# retrieve the row and column of point 0
sel = df['id'] == 0
row_0 = df.loc[sel, 'r'].values[0]
col_0 = df.loc[sel, 'c'].values[0]

# set to 0.9
iwdt[row_0, col_0]= 0.0

cost_dict={
    'dem': dem,
    'netcost': np.zeros_like(dem), # no previous netcost
    'cellsize': cellsize,
    'weight': 0.0,                 # it does not matter what we put here as netcost is made out of 0s
    'coef': coef
}

def test_iwdt_runs():
    calculate_iwdt(iwdt, cost_dict)
    assert True



