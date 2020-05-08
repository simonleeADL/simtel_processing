from ctapipe.io import event_source
import glob
import writer
from utilities import obtain_cleaning_mask, process_telescope, process_event, process_file, process_type

def process(input_path,output_path,max_files,
max_events,site_altitude,types,
telescopes,site_location,chop,id_no):
    
    for typename in types:
        
        print("Processing",typename)
        
        files = glob.glob(input_path + typename + '/*.simtel.gz')
        
        empty, telescope_extras, array_extras, runs_all, stereo, positions = process_type(
        files,max_files,max_events,site_altitude,telescopes,chop)
        
        if empty:
            continue
            
        writer.write(typename,output_path,site_location,
        array_extras,telescope_extras,runs_all,positions,stereo,id_no)
