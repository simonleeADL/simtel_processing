import pandas as pd
import numpy as np
from fact.io import to_h5py
import os
from datetime import datetime

# Add locations of each telescope based on given lat/long
def add_tel_location(telescope_events, site_location, positions):

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

        tel_lat = site_lat + (tel_north / r_earth) * (180 / pi)
        tel_lon = site_lon + (tel_east / r_earth) * \
            (180 / pi) / np.cos(site_lat * pi / 180)
        tel_alt = site_alt + tel_height

        dict_temp = {
            'tel_x': tel_north,
            'tel_y': -tel_east,
            'tel_latitude': tel_lat,
            'tel_longitude': tel_lon,
            'tel_altitude': tel_alt}

        telescope_location_list.append(dict_temp)

    telescope_events = telescope_events.join(
        pd.DataFrame(telescope_location_list))

    return telescope_events

def write(
        typename,
        output_path,
        site_location,
        array_events_data,
        telescope_events_data,
        runs_all,
        positions,
        stereo,
        id_no):

    print('Writing ' + typename + " data...",datetime.now().time().strftime("%H:%M:%S"))

    telescope_events = pd.DataFrame(telescope_events_data)
    array_events = pd.DataFrame(array_events_data)
    runs = pd.DataFrame(runs_all)
    
    # Calculate and add telescope location to telescope_events
    telescope_events = add_tel_location(
        telescope_events, site_location, positions)

    if typename == 'gamma-diffuse':
        output_file = output_path + 'gammas-diffuse' + str(id_no) + '.hdf5'
    else:
        output_file = output_path + typename + 's' + str(id_no) + '.hdf5'

    # Save to hdf5 file
    to_h5py(telescope_events, output_file, key='telescope_events', mode='w')
    to_h5py(array_events, output_file, key='array_events', mode='a')
    to_h5py(runs, output_file, key='runs', mode='a')
