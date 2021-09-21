import pandas as pd
import numpy as np
import astropy.units as u
from astropy.coordinates import SkyCoord, AltAz
from ctapipe.io import event_source
from ctapipe.calib import CameraCalibrator
from ctapipe.image.cleaning import apply_time_delta_cleaning
from ctapipe.image import hillas_parameters, leakage, number_of_islands, timing_parameters
from ctapipe.reco import HillasReconstructor
from ctapipe.image.extractor import NeighborPeakWindowSum
from datetime import datetime
import glob
import traceback
import os
import sys

def obtain_cleaning_mask(geom, image, time, camera_name):
	
# Cleaning levels taken from github.com/tudo-astroparticlephysics/cta_preprocessing

    cleaning_level = {
        'FlashCam': (10, 5),
        'CHEC': (3, 1.5),
    }

    # Settings
    picture_thresh, boundary_thresh = cleaning_level[camera_name]
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


def process_telescope(tel, dl1, stereo):

    geom = tel.camera.geometry
    camera_name = tel.camera.camera_name 
    image = dl1.image
    peak_time = dl1.peak_time

    # Cleaning using CHEC method
    clean = obtain_cleaning_mask(geom, image, peak_time, camera_name)

    # Skipping inadequate events
    if clean.sum() == 0:
        return None, None, None

    if stereo and clean.sum() < 5:
        return None, None, None
    
    # Get hillas parameters
    hillas_c = hillas_parameters(geom[clean], image[clean])

    if hillas_c.width == 0 or np.isnan(hillas_c.width.value):
        return None, None, None
    # Get leakage and islands
    leakage_c = leakage(geom, image, clean)
    n_islands, island_ids = number_of_islands(geom, clean)

    # Get time gradient
    tgrad = np.nan

    try:
        timing_c = timing_parameters(geom, image, peak_time, hillas_c, clean)
        tgrad = timing_c.slope.value
    except BaseException:
        print("Timing parameters didn't work. clean.sum() = " +
              str(clean.sum()), "\n")

    # Grab info for telescope_events
    tel_data_dict = {
        'nislands': n_islands,
        'telescope_type': tel.type,
        'camera_type': tel.camera.camera_name,
        'focal_length': tel.optics.equivalent_focal_length.value,
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

    return tel.type, tel_data_dict, hillas_c


def process_event(
        event,
        telescopes,
        subarray,
        stereo,
        calib,
        quality_cuts=[None,None,None,None]):
    
    if stereo:
        hillas_containers = {}
        telescope_pointings = {}
    
    horizon_frame = AltAz()
    reco = HillasReconstructor()
    
    # Calibrate the event
    calib(event)
    
    # Container to count what telescopes types were triggered in this event
    event_tels = []
    tel_data = []

    for tel_id, dl1 in event.dl1.tel.items():

        # If specific telescopes were specified, skip those that weren't
        if telescopes is not None:
            if (tel_id not in telescopes):
                continue

        tel = subarray.tels[tel_id]

        # Process the telescope event
        tel_type, tel_data_dict, hillas_c = process_telescope(tel, dl1, stereo)

        # Add the telescope event data if it wasn't skipped
        if tel_type is not None:

            if quality_cuts != [None,None,None,None]:

                intensity_cut, nislands_cut, n_survived_pixels_cut, intensity_width_1_cut = quality_cuts

                if intensity_cut != None and tel_data_dict['intensity'] <= intensity_cut:
                        continue
                if nislands_cut != None and tel_data_dict['nislands'] >= nislands_cut:
                        continue
                if n_survived_pixels_cut != None and tel_data_dict['n_survived_pixels'] <= n_survived_pixels_cut:
                        continue
                if intensity_width_1_cut != None and tel_data_dict['intensity_width_1'] >= intensity_width_1_cut:
                        continue

            # Calculate mc impact distance
            x1 = event.mc.core_x.value
            y1 = event.mc.core_y.value
            x2 = subarray.positions[tel_id][0].value
            y2 = subarray.positions[tel_id][1].value

            v = [x2 - x1, y2 - y1]
            mc_impact_distance = np.linalg.norm(v)


            # Adding a few extras to the telescope event data
            tel_data_dict.update({
                'array_event_id': event.index.event_id,
                'run_id': event.index.obs_id,
                'telescope_id': tel_id,
                'telescope_event_id': tel_id,
                'mc_impact_distance': mc_impact_distance})

            # Add data to the tables used for geometric reconstruction
            if stereo:
                
                tel_data_dict.update({'impact_distance': np.nan})

                hillas_containers[tel_id] = hillas_c

                telescope_pointings[tel_id] = SkyCoord(
                    alt=event.mc.tel[tel_id].altitude_raw * u.rad,
                    az=event.mc.tel[tel_id].azimuth_raw * u.rad,
                    frame=horizon_frame)

            event_tels.append(tel_type)
            tel_data.append(tel_data_dict)

    # Skip event if no telescopes were processed
    if len(event_tels) == 0:
        return None, None

    arr_data = {}

    if stereo and len(event_tels) > 1:
        array_pointing = SkyCoord(
            az=event.mcheader.run_array_direction[0],
            alt=event.mcheader.run_array_direction[1],
            frame=horizon_frame)

        # Do geometric direction reconstruction
        reconst = reco.predict(
            hillas_containers,
            subarray,
            array_pointing,
            telescope_pointings)
        
        # Calculate impact distance from reconstructed core position
        for i in range(len(tel_data)):
            x1 = reconst.core_x.value
            y1 = reconst.core_y.value
            x2 = subarray.positions[tel_data[i]['telescope_id']][0].value
            y2 = subarray.positions[tel_data[i]['telescope_id']][1].value

            v = [x2 - x1, y2 - y1]
            impact_distance = np.linalg.norm(v)
            
            tel_data[i]['impact_distance'] = impact_distance
        
        # Grab the extra stereo info for array_events

        reconst_az = reconst.az.value
        if reconst_az < -np.pi:
            reconst_az += 2*np.pi
        if reconst_az > np.pi:
            reconst_az -= 2*np.pi
        
        arr_data.update({
            'alt': reconst.alt.value,
            'alt_uncert': reconst.alt_uncert.value,
            'average_intensity': reconst.average_intensity,
            'az': reconst_az,
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

        # If only one telescope is triggered in a stereo array, replace
        # stereo features with NaN
        arr_data.update({
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

    # Grab info for array event data
    arr_data.update({
        'array_event_id': event.index.event_id,
        'run_id': event.index.obs_id,
        'azimuth_raw': event.mc.tel[tel_id].azimuth_raw,
        'altitude_raw': event.mc.tel[tel_id].altitude_raw,
        'azimuth_cor': event.mc.tel[tel_id].azimuth_cor,
        'altitude_cor': event.mc.tel[tel_id].altitude_cor,
        'num_triggered_lst': event_tels.count('LST'),
        'num_triggered_mst': event_tels.count('MST'),
        'num_triggered_sst': event_tels.count('SST'),

        'mc_energy': event.mc.energy.value,
        'true_source_alt': event.mc.alt.value,
        'true_source_az': event.mc.az.value,
        'mc_core_x': event.mc.core_x.value,
        'mc_core_y': event.mc.core_y.value,
        'mc_h_first_int': event.mc.h_first_int.value,
        'mc_x_max': event.mc.x_max.value,
        'mc_shower_primary_id': event.mc.shower_primary_id})

    return tel_data, arr_data

def process_file(
        filename,
        max_events,
        site_altitude,
        telescopes,
        stereo,
        quality_cuts=[None,None,None,None]):
    
    # Read source file
    try:
        source = event_source(
            filename,
            max_events=max_events,
            #back_seekable=True
            )
    except BaseException:
        print("Error: file " + filename + " could not be read", "\n")
        print(traceback.format_exc(), "\n")
        return None, None, None
    
    subarray = source.subarray

    # Exit if a telescope ID has been given that doesn't exist in the array
    if telescopes is not None:
        telescopes = [int(i) for i in telescopes.split(',')]
        for i in telescopes:
            if (i not in subarray.tel_ids):
                sys.exit('Error: tel_id ' + str(i) + ' is not in array ' + str(subarray.tel_ids))

    camera_name = source.subarray.tel[1].camera.camera_name

    if camera_name == 'CHEC':
        width = 6
        shift = 3
    elif camera_name == 'FlashCam':
        width = 4
        shift = 1
    else:
        sys.exit('Error: ' + camera_name + " camera not supported")
    
    image_extractor = NeighborPeakWindowSum(subarray=subarray, window_width=width, window_shift=shift)
    calib = CameraCalibrator(subarray=source.subarray, image_extractor=image_extractor)
    
    tested_altitude = False

    file_tel_data = []
    file_arr_data = []

    try:
        for event in source:
            if not tested_altitude:
                # Skip file if it's not at the given altitude
                tested_altitude = True
                if site_altitude is not None:
                    if event.mcheader.prod_site_alt.value != site_altitude:
                        raise ValueError

            # Process event
            tel_data, arr_data = process_event(
                event,
                telescopes,
                subarray,
                stereo,
                calib,
                quality_cuts=quality_cuts)

            if tel_data is not None:
                file_tel_data += tel_data
                file_arr_data.append(arr_data)
    # Skip file if it's at a different altitude than specified
    except ValueError:
        head, tail = os.path.split(filename)
        print(
            tail, "sims are not at an altitude of", str(site_altitude) + "m")
        return None, None, None

    # Skip file if it throws an unanticipated error on its own
    except BaseException:
        print(
            "Error: Something unanticipated went wrong\
            with processing an event in file ",
            filename)
        print(traceback.format_exc(), "\n")
        return None, None, None

    if len(file_tel_data) == 0:
        head, tail = os.path.split(filename)
        print(
            tail, "did not have any events output (maybe too low-energy?)")
        return None, None, None

    # Grab info for run data
    run_data = {
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
        #'prod_site_array': event.mcheader.prod_site_array,
        'prod_site_B_declination': event.mcheader.prod_site_B_declination.value,
        'prod_site_B_inclination': event.mcheader.prod_site_B_inclination.value,
        'prod_site_B_total': event.mcheader.prod_site_B_total.value,
        #'prod_site_coord': event.mcheader.prod_site_coord,
        #'prod_site_subarray': event.mcheader.prod_site_subarray,
        'run_id': event.index.obs_id,
        'shower_prog_id': event.mcheader.shower_prog_id,
        'shower_prog_start': event.mcheader.shower_prog_start,
        'shower_reuse': event.mcheader.shower_reuse,
        'simtel_version': event.mcheader.simtel_version,
        'spectral_index': event.mcheader.spectral_index,
    }

    return file_tel_data, file_arr_data, run_data


def process_type(
        files,
        max_files,
        max_events,
        site_altitude,
        telescopes,
        choppoints,
        quality_cuts=[None,None,None,None]):

    subarray = event_source(files[0], max_events=1).subarray
    positions = subarray.positions
    
    stereo = False
    
    # Flag array as stereo if need be
    if telescopes is not None and len(telescopes) > 1:
        stereo = True
    else:
        stereo = subarray.num_tels > 1

    telescope_events_data = []
    array_events_data = []
    runs_all = []

    # Paring down to the relevant files if "chop" is used
    if choppoints is not None:
        first_file = choppoints[0] - 1
        after_last_file = choppoints[1]
        files = files[first_file:after_last_file]

    # Define how many files to process
    files_available = len(files)
    files_to_process = files_available
    
    if max_files is not None and files_available > max_files:
        files_to_process = max_files

    file_no = 0

    for filename in files:

        file_no += 1

        if file_no > files_to_process:
            break

        print("File", file_no, "of", files_to_process,
              datetime.now().time().strftime("%H:%M:%S"),flush=True)
              
        # Process file
        file_tel_data, file_arr_data, run_data = process_file(
            filename, max_events, site_altitude, telescopes, stereo, quality_cuts)
        
        # Append run data if something was processed
        if run_data is not None:
            telescope_events_data += file_tel_data
            array_events_data += file_arr_data
            runs_all.append(run_data)

    if len(telescope_events_data) == 0:
        return None, None, None, None, None

    return telescope_events_data, array_events_data, runs_all, stereo, positions
