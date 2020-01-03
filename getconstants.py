# -*- coding: utf-8 -*-
"""
Created on Tue Jun 21 12:38:45 2016

@author: Ent00002
"""
#%%
import numpy as np
from netCDF4 import Dataset
from math import radians, cos, sin, acos, sqrt
#%%
# Function to estimate the distance of two points, unit: m
def great_circle(lon1, lat1, lon2, lat2):
    lon1, lat1, lon2, lat2 = map(radians, [lon1, lat1, lon2, lat2])
    Erad = 6.371e6 # [m] Earth radius
    return Erad * (
        acos(sin(lat1) * sin(lat2) + cos(lat1) * cos(lat2) * cos(lon1 - lon2))
    )

def getconstants(latnrs,lonnrs,invariant_data): # def getconstants in Python is the same as function in MATLAB. 
    
    # load the latitude and longitude from the invariants file
    latitude = Dataset(invariant_data, mode = 'r').variables['latitude'][latnrs] # [degrees north]
    longitude = Dataset(invariant_data, mode = 'r').variables['longitude'][lonnrs] # [degrees east]

    # Create land-sea-mask (in this model lakes are considered part of the land)
     # 0 = sea, 1 = land
    lsm=np.loadtxt('Landseamask.csv',delimiter=',')    
    lsm[0,:] = 0 # the northern boundary is always oceanic = 0
    lsm[-1,:] = 0 # the southern boundary is always oceanic = 0
    
    # Constants 
    g = 9.80665 # [m/s2] from ERA-interim archive
    density_water = 1000 # [kg/m3]
    dg = 111089.56 # [m] length of 1 degree latitude
    timestep = 6*3600 # [s] timestep in the ERA-interim archive (watch out! P & E have 3 hour timestep)
         
    # Semiconstants
    gridcell = np.abs(longitude[1] - longitude[0]) # [degrees] grid cell size

    # Area size calculation 
    A_gridcell = np.vstack(np.zeros((len(latitude))))
    lon_w=1
    lon_e=lon_w+gridcell
    l_ew=gridcell * dg
    for i in range(len(latitude)):   
       if (latitude[i] != 90 and latitude[i] != -90):
           lat_n=latitude[i]+gridcell / 2.0
           lat_s=latitude[i]-gridcell / 2.0
           l_n=great_circle(lon_w, lat_n, lon_e, lat_n) # [m] length northern boundary of a cell
           l_s=great_circle(lon_w, lat_s, lon_w, lat_s) # [m] length southern boundary of a cell
           l_diagonal=great_circle(lon_w, lat_s, lon_e, lat_n)
           # estimate the area of the gridcell as two triangles, and the area of 
           # triangle is estimated based on Heron's formula   
           p1=0.5*(l_n+l_ew+l_diagonal)
           p2=0.5*(l_s+l_ew+l_diagonal)
           # [m2] area size of grid cell
           A_gridcell[i] = sqrt(p1*(p1-l_n)*(p1-l_ew)*(p1-l_diagonal))+sqrt(p2*(p2-l_s)*(p2-l_ew)*(p2-l_diagonal))            
       elif latitude[i] == 90:
           lat_s=latitude[i]-gridcell / 2.0
           l_s=great_circle(lon_w, lat_s, lon_w, lat_s)
           p1=0.5*(l_s+2*l_ew)
           A_gridcell[i] = sqrt(p1*(p1-l_s)*(p1-l_ew)*(p1-l_ew))
       elif latitude[i] == -90:
           lat_n=latitude[i]+gridcell / 2.0
           l_n=great_circle(lon_w, lat_n, lon_w, lat_n)
           p1=0.5*(l_n+2*l_ew)
           A_gridcell[i] = sqrt(p1*(p1-l_n)*(p1-l_ew)*(p1-l_ew))
    
    lat_ns=latitude + gridcell / 2.0
    lat_ss=latitude - gridcell / 2.0
    L_N_gridcell = great_circle(lon_w,lat_ns, lon_e, lat_ns) # [m] length northern boundary of a cell
    L_S_gridcell = great_circle(lon_w,lat_ss, lon_e, lat_ss) # [m] length southern boundary of a cell
    L_EW_gridcell = gridcell * dg # [m] length eastern/western boundary of a cell 
        
    return latitude , longitude , lsm , g , density_water , timestep , A_gridcell , L_N_gridcell , L_S_gridcell , L_EW_gridcell , gridcell
