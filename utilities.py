import pandas as pd
import numpy as np
import astropy.units as u
from astropy.coordinates import SkyCoord, AltAz
from ctapipe.io import event_source
from ctapipe.calib import CameraCalibrator
from ctapipe.image.cleaning import tailcuts_clean, number_of_islands, apply_time_delta_cleaning
from ctapipe.image import hillas_parameters, leakage, concentration
from ctapipe.image.timing_parameters import timing_parameters
from ctapipe.reco import HillasReconstructor
#from ctapipe import utils
from datetime import datetime
import glob
import traceback

def obtain_cleaning_mask(geom, image, time):

    # Settings
    picture_thresh = 10
    boundary_thresh = 5
    min_number_picture_neighbors = 3
    time_limit = 5
    min_number_neighbors = 2

    # Select picture pixels
    pixels_above_picture = image >= picture_thresh

    # Require at least min_number_picture_neighbors.
    number_of_neighbors_above_picture = geom.neighbor_matrix_sparse.dot(
        pixels_above_picture.view(np.byte))
    pixels_in_picture = pixels_above_picture & (
        number_of_neighbors_above_picture >= min_number_picture_neighbors
        )

    # Select all boundary pixels (including picture pixels)
    pixels_above_boundary = image >= boundary_thresh

    # Remove boundary pixels not arrived in a given time frame
    pixels_to_keep = apply_time_delta_cleaning(
         geom, pixels_above_boundary, time, min_number_neighbors, time_limit
         ).astype(np.bool)
    pixels_with_picture_neighbors = geom.neighbor_matrix_sparse.dot(
        pixels_in_picture)
    mask = pixels_with_picture_neighbors & pixels_to_keep

    # remove isolated pixels (pixels with no neighbors)
    number_of_neighbors = geom.neighbor_matrix_sparse.dot(
        mask.view(np.byte))
    mask = mask & (number_of_neighbors >= 1)

    return mask

def process_telescope(tel_id,dl1,event,telescopes,subarray,stereo):
    
    # Skipping if this telescope isn't one that's been asked for
    
    if telescopes != None:
        if (tel_id not in telescopes):
            return

    # Cleaning, then saving Hillas paramaters, and skipping inadequate events

    tel = subarray.tels[tel_id]
    geom = tel.camera.geometry
    image = dl1.image / 0.15552960672840146
    peak_time = dl1.peak_time
    
    clean = obtain_cleaning_mask(geom, image, peak_time) 

    if clean.sum() == 0:
        return

    if stereo and clean.sum() < 5:
        return

    hillas_c = hillas_parameters(geom[clean], image[clean])                    

    if hillas_c.width == 0 or np.isnan(hillas_c.width.value):
        return

    # Get leakage and islands

    leakage_c = leakage(geom, image, clean)
    n_islands, island_ids = number_of_islands(geom, clean)
    
    global event_tels
    event_tels.append(tel.type)

    # Calculating mc impact distance

    x1 = event.mc.core_x.value
    y1 = event.mc.core_y.value
    x2 = subarray.positions[tel_id][0].value
    y2 = subarray.positions[tel_id][1].value

    v = [x2-x1,y2-y1]
    mc_impact_distance = np.linalg.norm(v)

    # Get time gradient
    tgrad = np.nan

    try:
        timing_c = timing_parameters(geom, image, peak_time, hillas_c, clean)
        tgrad = timing_c.slope.value
    except:
        print("Timing parameters didn't work. clean.sum() = " + str(clean.sum()),"\n")
    
    # Saving extra info for telescope_events

    dict_temp = {
        'array_event_id': event.index.event_id,
        'run_id': event.index.obs_id,
        'telescope_id': tel_id,
        'telescope_event_id': tel_id,
        'nislands': n_islands,
        'telescope_type': tel.type,
        'camera_type' : tel.camera.camera_name,
        'focal_length': tel.optics.equivalent_focal_length.value,
        'mc_impact_distance': mc_impact_distance,
        'n_survived_pixels': clean.sum(),
        'tgradient': tgrad,
        
        'x': hillas_c.x.value,
        'y': hillas_c.y.value,
        'r': hillas_c.r.value,
        'phi': hillas_c.phi.value,
        'intensity': hillas_c.intensity,
        'length': hillas_c.length.value,
        'width': hillas_c.width.value,
        'psi': hillas_c.psi.value,
        'skewness': hillas_c.skewness,
        'kurtosis': hillas_c.kurtosis,
        
        'pixels_width_1': leakage_c.pixels_width_1,
        'pixels_width_2': leakage_c.pixels_width_2,
        'intensity_width_1': leakage_c.intensity_width_1,
        'intensity_width_2': leakage_c.intensity_width_2}
    
    global telescope_extras
    telescope_extras.append(dict_temp)

    # Doing geometric stereo reconstruction
    
    global hillas_containers, telescope_pointings
    
    if stereo:
        hillas_containers[tel_id] = hillas_c

        telescope_pointings[tel_id] = SkyCoord(
            alt=event.mc.tel[tel_id].altitude_raw * u.rad,
            az=event.mc.tel[tel_id].azimuth_raw * u.rad,
            frame=horizon_frame)
        

def process_event(event,site_altitude,telescopes,subarray,stereo):

    calib = CameraCalibrator(subarray)

    # Skip file if it's not at the given altitude

    if site_altitude != None:
        if event.mcheader.prod_site_alt.value != site_altitude:
            raise ValueError

    # Exit if a telescope ID has been given that doesn't exist in the array

    if telescopes != None:
        for i in telescopes:
            if (i not in subarray.tel_ids):
                sys.exit('Error: tel_id given that is not in array')

    global hillas_containers, telescope_pointings
    
    hillas_containers = {}
    telescope_pointings = {}
    
    global horizon_frame
    
    horizon_frame = AltAz()
    reco = HillasReconstructor()
    
    # Calibrate the event
    
    calib(event)
    
    # Setting up a container to count what type of telescopes were triggered in this event
    
    global event_tels
    event_tels = []

    for tel_id, dl1 in event.dl1.tel.items():        
        process_telescope(tel_id,dl1,event,telescopes,subarray,stereo)

    if len(event_tels) == 0:
        return

    dict_temp = {}

    if stereo and len(event_tels) > 1:
        array_pointing = SkyCoord(
        az=event.mcheader.run_array_direction[0],
        alt=event.mcheader.run_array_direction[1],
        frame=horizon_frame)

        reconst = reco.predict(
            hillas_containers, subarray, array_pointing, telescope_pointings)

        # Grabbing the extra stereo info for array_events

        dict_temp.update({
            'alt': reconst.alt.value,
            'alt_uncert': reconst.alt_uncert.value,
            'average_intensity': reconst.average_intensity,
            'az': reconst.az.value,
            'az_uncert': reconst.az_uncert,
            'core_uncert': reconst.core_uncert,
            'core_x': reconst.core_x.value,
            'core_y': reconst.core_y.value,
            'goodness_of_fit': reconst.goodness_of_fit,
            'h_max': reconst.h_max.value,
            'h_max_uncert': reconst.h_max_uncert,
            'is_valid': reconst.is_valid,
            'stereo_flag': True})

    if stereo and len(event_tels) == 1:

        # If only one telescope is triggered in a stereo array, replace all stereo features with NaN

        dict_temp.update({
            'alt': np.nan,
            'alt_uncert': np.nan,
            'average_intensity': np.nan,
            'az': np.nan,
            'az_uncert': np.nan,
            'core_uncert': np.nan,
            'core_x': np.nan,
            'core_y': np.nan,
            'goodness_of_fit': np.nan,
            'h_max': np.nan,
            'h_max_uncert': np.nan,
            'is_valid': np.nan,
            'stereo_flag': False})

    # Grabbing the general extra info for array_events

    dict_temp.update({
        'array_event_id': event.index.event_id,
        'run_id': event.index.obs_id,
        'azimuth_raw': event.mc.tel[tel_id].azimuth_raw,
        'altitude_raw': event.mc.tel[tel_id].altitude_raw,
        'azimuth_cor': event.mc.tel[tel_id].azimuth_cor,
        'altitude_cor': event.mc.tel[tel_id].altitude_cor,
        'num_triggered_lst': event_tels.count('LST'),
        'num_triggered_mst': event_tels.count('MST'),
        'num_triggered_sst': event_tels.count('SST'),
        
        'energy': event.mc.energy.value,
        'alt': event.mc.alt.value,
        'az': event.mc.az.value,
        'core_x': event.mc.core_x.value,
        'core_y': event.mc.core_y.value,
        'h_first_int': event.mc.h_first_int.value,
        'x_max': event.mc.x_max.value,
        'shower_primary_id': event.mc.shower_primary_id})
    
    global array_extras    
    array_extras.append(dict_temp)

def process_file(filename,max_events,site_altitude,telescopes,stereo):

    try:
        source = event_source(filename, max_events=max_events, back_seekable=True)
    except:
        print("Error: file " + filename + " could not be read","\n")
        return

    try:
        for event in source:
            try:
                process_event(event,site_altitude,telescopes,source.subarray,stereo)
            except ValueError:
                pass     
    except:
        print("Error: Something unanticipated went wrong with processing an event in file ",filename)
        print(traceback.format_exc(),"\n")
        return
    
    if site_altitude != None:
        if event.mcheader.prod_site_alt.value != site_altitude:
            return

    # Grabbing all the info for runs
    try:
        test = event.mcheader.atmosphere
    except:
        print("For some reason it won't event in file",filename)
        return

    dict_temp = {
        'atmosphere': event.mcheader.atmosphere,
        'core_pos_mode': event.mcheader.core_pos_mode,
        'corsika_bunchsize': event.mcheader.corsika_bunchsize,
        'corsika_high_E_detail': event.mcheader.corsika_high_E_detail,
        'corsika_high_E_model': event.mcheader.corsika_high_E_model,
        'corsika_iact_options': event.mcheader.corsika_iact_options,
        'corsika_low_E_detail': event.mcheader.corsika_low_E_detail,
        'corsika_low_E_model': event.mcheader.corsika_low_E_model,
        'corsika_version': event.mcheader.corsika_version,
        'corsika_wlen_max': event.mcheader.corsika_wlen_max.value,
        'corsika_wlen_min': event.mcheader.corsika_wlen_min.value,
        'detector_prog_id': event.mcheader.detector_prog_id,
        'detector_prog_start': event.mcheader.detector_prog_start,
        'diffuse': event.mcheader.diffuse,
        'energy_range_max': event.mcheader.energy_range_max.value,
        'energy_range_min': event.mcheader.energy_range_min.value,
        'injection_height': event.mcheader.injection_height.value,
        'max_alt': event.mcheader.max_alt.value,
        'max_az': event.mcheader.max_az.value,
        'max_scatter_range': event.mcheader.max_scatter_range.value,
        'max_viewcone_radius': event.mcheader.max_viewcone_radius.value,
        'min_alt': event.mcheader.min_alt.value,
        'min_az': event.mcheader.min_az.value,
        'min_scatter_range': event.mcheader.min_scatter_range.value,
        'min_viewcone_radius': event.mcheader.min_viewcone_radius.value,
        'num_showers': event.mcheader.num_showers,
        'prod_site_alt': event.mcheader.prod_site_alt.value,
        'prod_site_array': event.mcheader.prod_site_array,
        'prod_site_B_declination': event.mcheader.prod_site_B_declination.value,
        'prod_site_B_inclination': event.mcheader.prod_site_B_inclination.value,
        'prod_site_B_total': event.mcheader.prod_site_B_total.value,
        'prod_site_coord': event.mcheader.prod_site_coord,
        'prod_site_subarray': event.mcheader.prod_site_subarray,
        'run_id': event.index.obs_id,
        'shower_prog_id': event.mcheader.shower_prog_id,
        'shower_prog_start': event.mcheader.shower_prog_start,
        'shower_reuse': event.mcheader.shower_reuse,
        'simtel_version': event.mcheader.simtel_version,
        'spectral_index': event.mcheader.spectral_index,
        }
    global runs_all
    runs_all.append(dict_temp)

def process_type(files,max_files,max_events,site_altitude,telescopes,chop):

    subarray = event_source(files[0], max_events=1).subarray         
    positions = subarray.positions
    
    stereo = False
    
    if telescopes != None and len(telescopes) > 1:
        stereo = True
    else:
        stereo = subarray.num_tels > 1
    
    global telescope_extras, array_extras, runs_all
    
    telescope_extras = []
    array_extras = []
    runs_all = []
    
    # Paring down to the relevant files if "chop" is used
    
    if chop != None:
        first_file = choppoints[0] - 1
        after_last_file = choppoints[1]
        files = files[first_file:after_last_file]
    
    files_available = len(files)
    files_to_process = files_available
    
    if max_files != None and files_available > max_files:
        files_to_process = max_files
         
    file_no = 0
        
    for filename in files:
    
        file_no += 1

        if file_no > files_to_process:
            break
            
        print("File",file_no,"of",files_to_process,datetime.now().time().strftime("%H:%M:%S"))
        
        process_file(filename,max_events,site_altitude,telescopes,stereo)
    
    empty = False
    if len(telescope_extras) == 0:
        empty = True
        
    return empty, telescope_extras, array_extras, runs_all, stereo, positions
