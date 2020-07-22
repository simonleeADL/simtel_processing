import pandas as pd
import numpy as np
import astropy.units as u
from astropy.coordinates import SkyCoord, AltAz
from ctapipe.io import event_source
from ctapipe.calib import CameraCalibrator
from ctapipe.image.cleaning import number_of_islands, apply_time_delta_cleaning
from ctapipe.image import hillas_parameters, leakage
from ctapipe.image.timing_parameters import timing_parameters
from ctapipe.reco import HillasReconstructor
from datetime import datetime
import glob
import traceback
import os


def obtain_cleaning_mask(geom, image, time, camera_name):
	
# Cleaning levels taken from github.com/tudo-astroparticlephysics/cta_preprocessing

    cleaning_level = {
        'ASTRICam': (5, 7),  # (5, 10)?
        'FlashCam': (12, 15),
        'LSTCam': (3.5, 7.5),
        'NectarCam': (3, 5.5),
        'FlashCam': (10, 5),  # there is some scaling missing?
        'DigiCam': (2, 4.5),
        'CHEC': (10, 5),
        'SCTCam': (1.5, 3)
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
    
    # This is unfinished because I don't know what other
    # telescopes use to convert to photoelectrons
    pe_conversion_factor = 1
    if camera_name == 'CHEC':
        pe_conversion_factor = 0.15552960672840146
    image = dl1.image / pe_conversion_factor

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
    tel_data = {
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

    return tel.type, tel_data, hillas_c


def process_event(
        event,
        telescopes,
        subarray,
        stereo,
        telescope_events_data,
        calib):

    if stereo:
        hillas_containers = {}
        telescope_pointings = {}

    horizon_frame = AltAz()
    reco = HillasReconstructor()
    # Calibrate the event
    calib(event)

    # Container to count what telescopes types were triggered in this event
    event_tels = []
    for tel_id, dl1 in event.dl1.tel.items():

        # If specific telescopes were specified, skip those that weren't
        if telescopes is not None:
            if (tel_id not in telescopes):
                print("There's a telescope that wasn't specified")
                continue

        tel = subarray.tels[tel_id]

        # Process the telescope event
        tel_type, tel_data, hillas_c = process_telescope(tel, dl1, stereo)
        # Add the telescope event data if it wasn't skipped
        if tel_type is not None:

            # Calculate mc impact distance
            x1 = event.mc.core_x.value
            y1 = event.mc.core_y.value
            x2 = subarray.positions[tel_id][0].value
            y2 = subarray.positions[tel_id][1].value

            v = [x2 - x1, y2 - y1]
            mc_impact_distance = np.linalg.norm(v)

            # Adding a few extras to the telescope event data
            tel_data.update({
                'array_event_id': event.index.event_id,
                'run_id': event.index.obs_id,
                'telescope_id': tel_id,
                'telescope_event_id': tel_id,
                'mc_impact_distance': mc_impact_distance})

            # Update telescope events table
            telescope_events_data.append(tel_data)
            event_tels.append(tel_type)

            # Add data to the tables used for geometric reconstruction
            if stereo:
                hillas_containers[tel_id] = hillas_c

                telescope_pointings[tel_id] = SkyCoord(
                    alt=event.mc.tel[tel_id].altitude_raw * u.rad,
                    az=event.mc.tel[tel_id].azimuth_raw * u.rad,
                    frame=horizon_frame)

    # Skip event if no telescopes were processed
    if len(event_tels) == 0:
        return telescope_events_data, None

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

        # Grab the extra stereo info for array_events
        arr_data.update({
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

        'energy': event.mc.energy.value,
        'alt': event.mc.alt.value,
        'az': event.mc.az.value,
        'core_x': event.mc.core_x.value,
        'core_y': event.mc.core_y.value,
        'h_first_int': event.mc.h_first_int.value,
        'x_max': event.mc.x_max.value,
        'shower_primary_id': event.mc.shower_primary_id})

    return telescope_events_data, arr_data


def process_file(
        filename,
        max_events,
        site_altitude,
        telescopes,
        stereo,
        telescope_events_data,
        array_events_data):
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
        return telescope_events_data, array_events_data, None
    subarray = source.subarray
    calib = CameraCalibrator(subarray)

    # Exit if a telescope ID has been given that doesn't exist in the array
    if telescopes is not None:
        for i in telescopes:
            if (i not in subarray.tel_ids):
                sys.exit('Error: tel_id given that is not in array')

    tested_altitude = False
    appended_anything = False

    try:
        for event in source:
            if not tested_altitude:
                # Skip file if it's not at the given altitude
                tested_altitude = True
                if site_altitude is not None:
                    if event.mcheader.prod_site_alt.value != site_altitude:
                        raise ValueError

            # Process event
            telescope_events_data, arr_data = process_event(
                event,
                telescopes,
                subarray,
                stereo,
                telescope_events_data,
                calib)

            if arr_data is not None:
                array_events_data.append(arr_data)
                appended_anything = True

    # Skip file if it's at a different altitude than specified
    except ValueError:
        head, tail = os.path.split(filename)
        print(
            tail, "sims are not at an altitude of", str(site_altitude) + "m")
        return telescope_events_data, array_events_data, None

    # Skip file if it throws an unanticipated error on its own
    except BaseException:
        print(
            "Error: Something unanticipated went wrong\
            with processing an event in file ",
            filename)
        print(traceback.format_exc(), "\n")
        return telescope_events_data, array_events_data, None

    if not appended_anything:
        head, tail = os.path.split(filename)
        print(
            tail, "did not have any events output (maybe too low-energy?)")
        return telescope_events_data, array_events_data, None

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

    return telescope_events_data, array_events_data, run_data


def process_type(
        files,
        max_files,
        max_events,
        site_altitude,
        telescopes,
        choppoints):

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
        telescope_events_data, array_events_data, run_data = process_file(
            filename, max_events, site_altitude, telescopes, stereo, telescope_events_data, array_events_data)
        
        # Append run data if something was processed
        if run_data is not None:
            runs_all.append(run_data)

    if len(telescope_events_data) == 0:
        return None, None, None, None, None

    return telescope_events_data, array_events_data, runs_all, stereo, positions
