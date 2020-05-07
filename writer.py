import pandas as pd
import numpy as np
from fact.io import to_h5py
import os

def add_tel_location(telescope_events,site_location,positions):

    telescope_location_list = []
    
    site_lat = site_location[0]
    site_lon = site_location[1]
    site_alt = site_location[2]
    r_earth = 6371000.0
    pi = np.pi
    
    for i in range(len(telescope_events.index)):
    
        tel_id = telescope_events['telescope_event_id'][i]
    
        tel_north = positions[tel_id][0].value
        tel_east = - positions[tel_id][1].value
        tel_height = positions[tel_id][2].value

        tel_lat = site_lat  + (tel_north / r_earth) * (180 / pi)
        tel_lon = site_lon + (tel_east / r_earth) * (180 / pi) / np.cos(site_lat * pi/180)
        tel_alt = site_alt + tel_height

        dict_temp = {
            'tel_latitude': tel_lat,
            'tel_longitude': tel_lon,
            'tel_altitude': tel_alt}
        
        telescope_location_list.append(dict_temp)

    telescope_events = telescope_events.join(pd.DataFrame(telescope_location_list))
    
    return telescope_events
    
def add_impact_distance(telescope_events,array_events,positions):

    impact_distance_list = []

    for i in range(len(telescope_events.index)):
        impact_distance = np.nan

        tel_id = telescope_events['telescope_event_id'][i]
        tel_arr_id = telescope_events['array_event_id'][i]
        arr_arr_id = array_events['array_event_id']
        arr_index = np.where(arr_arr_id == tel_arr_id)[0][0]

        if not np.isnan(array_events['mc_core_x'][arr_index]):
            x1 = array_events['mc_core_x'][arr_index]
            y1 = array_events['mc_core_y'][arr_index]
            x2 = positions[tel_id][0].value
            y2 = positions[tel_id][1].value

            v = [x2-x1,y2-y1]
            impact_distance = np.linalg.norm(v)

        dict_temp = {'impact_distance': impact_distance}
        impact_distance_list.append(dict_temp)

    telescope_events = telescope_events.join(pd.DataFrame(impact_distance_list))
    
    return telescope_events

def write(typename,output_path,site_location,array_extras,telescope_extras,runs_all,positions,stereo,id_no):
    
    print('Writing ' + typename + " data...")

    telescope_events = pd.DataFrame(telescope_extras)
    array_events = pd.DataFrame(array_extras)
    runs = pd.DataFrame(runs_all)
    
    # Renaming monte carlo parameters
    
    array_events.rename(columns = {"alt": "true_source_alt", 
                     "az": "true_source_az",
                     "core_x": "mc_core_x",
                     "core_y": "mc_core_y",
                     "energy": "mc_energy",
                     "h_first_int": "mc_h_first_int",
                     "shower_primary_id": "mc_shower_primary_id",
                     "x_max": "mc_x_max",},inplace = True)  
        
    # Calculate and add telescope location to telescope_events
    
    telescope_events = add_tel_location(telescope_events,site_location,positions)
    
    # Calculate and add impact distance to telescope_events
    
    if stereo:
        telescope_events = add_impact_distance(telescope_events,array_events,positions)
    
    if typename == 'gamma-diffuse':
        output_file = output_path + 'gammas-diffuse' + str(id_no) + '.hdf5' 
    else:
        output_file = output_path + typename + 's' + str(id_no) + '.hdf5'

    # Save to hdf5 file
    
    to_h5py(telescope_events, output_file, key='telescope_events', mode='w')
    to_h5py(array_events, output_file, key='array_events', mode='a')
    to_h5py(runs, output_file, key='runs', mode='a')

