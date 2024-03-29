from ctapipe.io import event_source
import glob
import writer
from utilities import process_type


def process(
        input_path,
        output_path,
        max_files,
        max_events,
        site_altitude,
        types,
        telescopes,
        site_location,
        choppoints,
        id_no,
        quality_cuts):

    for typename in types:

        print("Processing", typename)
        
        # Get a list of the files for this source type
        files = glob.glob(input_path + typename + '/*.simtel.zst')
        files = files + glob.glob(input_path + typename + '/*.simtel.gz')

        if len(files)==0:
            print("No ",typename," files found")
            continue
        
        # Process the files
        telescope_events_data, array_events_data, runs_all, stereo, positions = process_type(
            files, max_files, max_events, site_altitude, telescopes, choppoints, quality_cuts)
        
        # Skip writing if nothing was processed
        if telescope_events_data is None:
            print(
            typename, "did not have any events output (maybe too low-energy?)")
            continue
        
        site_location.append(runs_all[0]['prod_site_alt'])

        writer.write(
            typename,
            output_path,
            site_location,
            array_events_data,
            telescope_events_data,
            runs_all,
            positions,
            stereo,
            id_no)
