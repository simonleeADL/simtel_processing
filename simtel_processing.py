import click

@click.command()

@click.option('-l','--location','location', type=str, required=True, help='latitude, longitude, and altitude (metres) of the site in the form lat,lon,alt e.g.: -34.9195,138.6030,52')
@click.option('-t','--type','source_type', type=click.Choice(['g', 'p', 'd']), help='Process only gamma (g), proton (p), or diffuse gammas (d)')
@click.option('-i','--input','input_path', default='./', type=click.Path(exists=True), help='Input directory, must include "gamma", "proton", and "gamma-diffuse" folders full of simtel.gz files. Defaults to current folder.')
@click.option('-o', '--output','output_path', type=click.Path(), help='Output path, default will create Data in input path')
@click.option('--max_files', type=int, help='Maximum number of files to process per type, defaults to all files')
@click.option('--max_events', type=int, default=99999999, help='Maximum number of events to be processed per file, defaults to 99999999')
@click.option('--tels','telescopes', type=str, help='Limit processing to only these telescope IDs (minimum one). List multiple with commas in between, e.g.: 1,3,12')
@click.option('--altitude','site_altitude', type=int, help='Only process runs of this altitude (useful when there are multiple runs mixed up)')
@click.option('--chop','chop', type=str, help='Process only the nth through mth files (starting at 1) with optional output ID. Can only be used when processing specific type (g/p/d) and not with max_files. Formatted with commas between first and last file no., plus an optional ID no. as the third, e.g.: 30,49 or 0,19,3')

def main(input_path,output_path,max_files,max_events,site_altitude,source_type,telescopes,location,chop):

    import os
    import sys
    import glob

    # Making sure all the inputs are sensible and won't throw errors

    if output_path == None:
        output_path = input_path + 'data/'

    if not output_path.endswith('/'):
        output_path = output_path + '/'
    
    if not os.path.exists(output_path):
        os.mkdir(output_path)

    if source_type == None:
        types=['gamma','proton','gamma-diffuse']
        if chop != None:
            sys.exit("Error: Cannot use 'chop' when processing 'all' types")
    if source_type == 'g':
        types=['gamma']
    if source_type == 'p':
        types=['proton']
    if source_type == 'd':
        types=['gamma-diffuse']

    for typename in types:
        if not os.path.isdir(input_path + typename):
            sys.exit('There is no ' + typename + ' directory')

    if telescopes != None:
        try:
            telescopes = [int(i) for i in telescopes.split(',')]
        except:
            sys.exit("Error: Some of those telescope IDs aren't integers")
            
    id_no = ""

    if chop != None:
        if max_files != None:
            sys.exit("Error: Can't chop while also specifying max_files")
        try:
            choppoints = [int(i) for i in chop.split(',')]
        except:
            sys.exit("Error: Invalid chop (should be two integers separated by commas, indicating nth and mth files to process, plus an optional third ID no.)")
        if len(choppoints) < 2:
            sys.exit("Error: Not enough chop points!")
        if len(choppoints) > 3:
            sys.exit("Error: Too many chop points!")
        if choppoints[0] < 1 or choppoints[1] > len(glob.glob(input_path + typename + '/*.simtel.gz')):
            sys.exit("Error: Chop out of file range!")
        if len(choppoints) == 3:
            id_no = choppoints[2]
            
    try:
        site_location = [float(i) for i in location.split(',')]
    except:
        sys.exit("Error: Invalid lat/lon/alt input (should be three numbers separated by commas)")
    if len(site_location) < 3:
        sys.exit("Error: Not enough location arguments")
    if len(site_location) > 3:
        sys.exit("Error: Too many location arguments")
    if not -90 <= site_location[0] <= 90:
        sys.exit("Error: Invalid latitude")
    if not -180 <= site_location[1] <= 180:
        sys.exit("Error: Invalid longitude")
    if not -450 <= site_location[2] <= 9000:
        sys.exit("Error: Invalid height (no longer on land)")

    print("Importing (ctapipe.calib takes a while)...")

    from timeit import default_timer as timer
    import pandas as pd
    import numpy as np
    import astropy.units as u
    from astropy.coordinates import SkyCoord, AltAz
    from ctapipe.io import event_source
    from ctapipe.utils.datasets import get_dataset_path
    from ctapipe.calib import CameraCalibrator
    from ctapipe.image.cleaning import tailcuts_clean, number_of_islands
    from ctapipe.image import hillas_parameters, leakage, concentration
    from ctapipe.image.timing_parameters import timing_parameters
    from ctapipe.reco import HillasReconstructor
    from ctapipe.io import HDF5TableWriter
    from ctapipe import utils
    from fact.io import to_h5py
    import datetime
    import warnings
    import random
    import traceback
    
    from CHECOnsky.calib import obtain_cleaning_mask

    print('...done importing!')

    # I don't like warnings

    #warnings.filterwarnings("ignore")

    # Cleaning levels taken from github.com/tudo-astroparticlephysics/cta_preprocessing

    cleaning_level = {
        'ASTRICam': (5, 7, 2),  # (5, 10)?
        'FlashCam': (12, 15, 2),
        'LSTCam': (3.5, 7.5, 2),
        'NectarCam': (3, 5.5, 2),
        'FlashCam': (4, 8, 2),  # there is some scaling missing?
        'DigiCam': (2, 4.5, 2),
        'CHEC': (2, 4, 2),
        'SCTCam': (1.5, 3, 2)
    }

    # Filename for a temporary file that the HDF5 writer can writer to
    
    process_temp = str(datetime.datetime.now().strftime("%Y-%m-%d_%H:%M:%S:%f")) + str(random.randrange(100000))

    stereo = False
    checked_stereo = False
    secondsleft = 0.0
    nothingprocessed = True
    
    # If telescopes are being specified, check if it's mono or stereo

    if telescopes != None:
        checked_stereo = True
        if len(telescopes) > 1:
            stereo = True
    
    horizon_frame = AltAz()
    reco = HillasReconstructor()

    for typename in types:
        
        # Setting up containers for the extra data that isn't generated in the DL1s by default

        telescope_extras = []
        array_extras = []
        runs_all = []
        
        # Finding the files and making a list
        
        files = glob.glob(input_path + typename + '/*.simtel.gz')
        
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
        
        with HDF5TableWriter(filename=process_temp, overwrite=True,  mode='a', group_name='events') as writer:
        
            for filename in files:

                file_no += 1

                if file_no > files_to_process:
                    break

                # Print out what file we're up to and estimate when it will finish

                finishtime = datetime.datetime.now() + datetime.timedelta(seconds = secondsleft)
                processing_text = typename + ' file no. ' + str(file_no) + ' of ' + str(files_to_process)
                if secondsleft != 0.0:
                    processing_text = processing_text + '. Finish: ' + str(finishtime.strftime("%H:%M")) + ' (' + str(round(secondsleft/60,2)) + ' minutes)'
                print(processing_text, end='\r', flush=True)

                if file_no == 1:
                    start = timer()

                try:
                    source = event_source(filename, max_events=max_events, back_seekable=True)
                except:
                    print("Error: file " + filename + " could not be read","\n")
                    continue

                try:
                    for event in source:

                        calib = CameraCalibrator(event.inst.subarray)

                        # Skip file if it's not at the given altitude

                        if site_altitude != None:
                            if event.mcheader.prod_site_alt.value != site_altitude:
                                break

                        # Turn on the stereo flag if there's more than one telescope

                        if not checked_stereo:
                            stereo = event.inst.subarray.num_tels > 1
                            checked_stero = True

                        # Exit if a telescope ID has been given that doesn't exist in the array

                        if telescopes != None:
                            for i in telescopes:
                                if (i not in event.inst.subarray.tel_ids):
                                    sys.exit('Error: tel_id given that is not in array')

                        hillas_containers = {}
                        telescope_pointings = {}

                        
                        # Calibrate the event
                        
                        calib(event)
                        
                        # Setting up a container to count what type of telescopes were triggered in this event
                        
                        event_tels = []

                        for tel_id, dl1 in event.dl1.tel.items():

                            # Skipping if this telescope isn't one that's been asked for

                            if telescopes != None:
                                if (tel_id not in telescopes):
                                    continue

                            # Cleaning, then saving Hillas paramaters, and skipping inadequate events

                            tel = event.inst.subarray.tels[tel_id]
                            geom = tel.camera.geometry
                            image = dl1.image / 0.15552960672840146
                            pulse_time = dl1.pulse_time

                            boundary, picture, min_neighbors = cleaning_level[tel.camera.camera_name]
                            
                            clean = obtain_cleaning_mask(geom, image, pulse_time)

                            #clean = tailcuts_clean(camera, image, boundary_thresh=boundary, picture_thresh=picture, min_number_picture_neighbors=min_neighbors)   

                            if clean.sum() == 0:
                                continue

                            if stereo and clean.sum() < 5:
                                continue

                            hillas_c = hillas_parameters(geom[clean], image[clean])                    

                            if hillas_c.width == 0 or np.isnan(hillas_c.width.value):
                                continue

                            # Get leakage and islands

                            leakage_c = leakage(geom, image, clean)
                            n_islands, island_ids = number_of_islands(geom, clean)

                            event_tels.append(tel.type)

                            # Calculating mc impact distance

                            x1 = event.mc.core_x.value
                            y1 = event.mc.core_y.value
                            x2 = event.inst.subarray.positions[tel_id][0].value
                            y2 = event.inst.subarray.positions[tel_id][1].value

                            v = [x2-x1,y2-y1]
                            mc_impact_distance = np.linalg.norm(v)

                            # Get time gradient
                            tgrad = np.nan

                            try:
                                timing_c = timing_parameters(geom, image, pulse_time, hillas_c, clean)
                                tgrad = timing_c.slope.value
                            except:
                                print("Timing parameters didn't work. clean.sum() = " + str(clean.sum()),"\n")
                            
                            # Saving extra info for telescope_events

                            dict_temp = {
                                'array_event_id': event.r0.event_id,
                                'run_id': event.r0.obs_id,
                                'telescope_id': tel_id,
                                'telescope_event_id': tel_id,
                                'nislands': n_islands,
                                'telescope_type': tel.type,
                                'camera_type' : tel.camera.camera_name,
                                'focal_length': tel.optics.equivalent_focal_length.value,
                                'mc_impact_distance': mc_impact_distance,
                                'n_survived_pixels': clean.sum(),
                                'tgradient': tgrad}
                            telescope_extras.append(dict_temp)

                            writer.write("telescope_events", hillas_c)
                            writer.write("leakage", leakage_c)
                            if nothingprocessed:
                                nothingprocessed = False

                            # Doing geometric stereo reconstruction

                            if stereo:

                                hillas_containers[tel_id] = hillas_c

                                telescope_pointings[tel_id] = SkyCoord(
                                    alt=event.mc.tel[tel_id].altitude_raw * u.rad,
                                    az=event.mc.tel[tel_id].azimuth_raw * u.rad,
                                    frame=horizon_frame)

                        if len(event_tels) == 0:
                            continue

                        dict_temp = {}

                        if stereo and len(event_tels) > 1:
                            array_pointing = SkyCoord(
                            az=event.mcheader.run_array_direction[0],
                            alt=event.mcheader.run_array_direction[1],
                            frame=horizon_frame)

                            reconst = reco.predict(
                                hillas_containers, event.inst, array_pointing, telescope_pointings)

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

                        writer.write("array_events", event.mc)

                        # Grabbing the general extra info for array_events

                        dict_temp.update({
                            'array_event_id': event.r0.event_id,
                            'run_id': event.r0.obs_id,
                            'azimuth_raw': event.mc.tel[tel_id].azimuth_raw,
                            'altitude_raw': event.mc.tel[tel_id].altitude_raw,
                            'azimuth_cor': event.mc.tel[tel_id].azimuth_cor,
                            'altitude_cor': event.mc.tel[tel_id].altitude_cor,
                            'num_triggered_lst': event_tels.count('LST'),
                            'num_triggered_mst': event_tels.count('MST'),
                            'num_triggered_sst': event_tels.count('SST')})
                            
                        array_extras.append(dict_temp)
                except Exception as e: #except:
                    print("Error: Something unanticipated went wrong with processing an event in file ",filename)
                    print(traceback.format_exc(),"\n")
                    continue
                
                if site_altitude != None:
                    if event.mcheader.prod_site_alt.value != site_altitude:
                        break

                # Grabbing all the info for runs
                try:
                    test = event.mcheader.atmosphere
                except:
                    print("For some reason it won't event in file",filename)
                    continue

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
                    'run_id': event.r0.obs_id,
                    'shower_prog_id': event.mcheader.shower_prog_id,
                    'shower_prog_start': event.mcheader.shower_prog_start,
                    'shower_reuse': event.mcheader.shower_reuse,
                    'simtel_version': event.mcheader.simtel_version,
                    'spectral_index': event.mcheader.spectral_index,
                    }
                runs_all.append(dict_temp)

                end = timer()
                secondsleft = round((((end - start)/file_no)*(files_to_process-file_no)),2)
        
        sys.stdout.write('\n')

        if nothingprocessed:
            print("Nothing was succesfuly processed for",typename)
            break

        print('Writing ' + typename + " file")

        # Read tables from processing file
        
        telescope_events = pd.read_hdf(process_temp, key='events/telescope_events')
        array_events = pd.read_hdf(process_temp, key='events/array_events')
        leakage_df = pd.read_hdf(process_temp, key='events/leakage')
        
        # Renaming monte carlo parameters
        
        array_events.rename(columns = {"alt": "true_source_alt", 
                         "az": "true_source_az",
                         "core_x": "mc_core_x",
                         "core_y": "mc_core_y",
                         "energy": "mc_energy",
                         "h_first_int": "mc_h_first_int",
                         "shower_primary_id": "mc_shower_primary_id",
                         "x_max": "mc_x_max",},inplace = True)

        # Some necessary unit conversion because the HDF5 writer saves angles as degrees
        
        array_events[['true_source_az','true_source_alt']] = array_events[['true_source_az','true_source_alt']].apply(np.radians)
        telescope_events[['phi','psi']] = telescope_events[['phi','psi']].apply(np.radians)
        
        # Add in general extra data
        
        telescope_events = telescope_events.join(leakage_df)
        array_events = array_events.join(pd.DataFrame(array_extras))
        telescope_events = telescope_events.join(pd.DataFrame(telescope_extras))
        runs = pd.DataFrame(runs_all)
        
        # Calculate and add telescope location to telescope_events
        
        telescope_location_list = []
        
        site_lat = site_location[0]
        site_lon = site_location[1]
        site_alt = site_location[2]
        r_earth = 6371000.0
        pi = np.pi
        
        for i in range(len(telescope_events.index)):
        
            tel_id = telescope_events['telescope_event_id'][i]
        
            tel_north = event.inst.subarray.positions[tel_id][0].value
            tel_east = - event.inst.subarray.positions[tel_id][1].value
            tel_height = event.inst.subarray.positions[tel_id][2].value

            tel_lat = site_lat  + (tel_north / r_earth) * (180 / pi)
            tel_lon = site_lon + (tel_east / r_earth) * (180 / pi) / np.cos(site_lat * pi/180)
            tel_alt = site_alt + tel_height

            dict_temp = {
                'tel_latitude': tel_lat,
                'tel_longitude': tel_lon,
                'tel_altitude': tel_alt}
            
            telescope_location_list.append(dict_temp)

        telescope_events = telescope_events.join(pd.DataFrame(telescope_location_list))
        
        # Calculate and add impact distance to telescope_events

        if stereo:
            impact_distance_list = []

            for i in range(len(telescope_events.index)):
                impact_distance = np.nan

                tel_id = telescope_events['telescope_event_id'][i]
                tel_arr_id = telescope_events['array_event_id'][i]
                arr_arr_id = array_events['array_event_id']
                arr_index = np.where(arr_arr_id == tel_arr_id)[0][0]

                if not np.isnan(array_events['core_x'][arr_index]):
                    x1 = array_events['core_x'][arr_index]
                    y1 = array_events['core_y'][arr_index]
                    x2 = event.inst.subarray.positions[tel_id][0].value
                    y2 = event.inst.subarray.positions[tel_id][1].value

                    v = [x2-x1,y2-y1]
                    impact_distance = np.linalg.norm(v)

                dict_temp = {'impact_distance': impact_distance}
                impact_distance_list.append(dict_temp)

            telescope_events = telescope_events.join(pd.DataFrame(impact_distance_list))
        
        if typename == 'gamma-diffuse':
            output_file = output_path + 'gammas-diffuse' + str(id_no) + '.hdf5' 
        else:
            output_file = output_path + typename + 's' + str(id_no) + '.hdf5'

        # Save to hdf5 file
        
        to_h5py(telescope_events, output_file, key='telescope_events', mode='w')
        to_h5py(array_events, output_file, key='array_events', mode='a')
        to_h5py(runs, output_file, key='runs', mode='a')

        os.remove(process_temp)

    print('Finished')

if __name__ == '__main__':
    main()
