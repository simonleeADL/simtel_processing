import sys
import processor
import os
import glob
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
    
    input_path, output_path, types, site_location, choppoints, id_no = process_inputs(input_path,output_path,source_type, location, chop)

    input_validity = validate(input_path,output_path,max_files,
    max_events,site_altitude,source_type,telescopes,location,
    chop,types)
    
    if input_validity != "Valid":
        sys.exit(input_validity)
            
    processor.process(input_path,output_path,max_files,
    max_events,site_altitude,
    source_type,telescopes,
    site_location,chop,id_no)
    
    print("Wow we actually made it, cool")
   
def process_inputs(input_path,output_path,source_type, location, chop):

    id_no = ""
    
    if not input_path.endswith('/'):
        input_path = input_path + '/'

    if output_path == None:
        output_path = input_path + 'data/'

    if not output_path.endswith('/'):
        output_path = output_path + '/'
    
    if not os.path.exists(output_path):
        os.mkdir(output_path)

    if source_type == None:
        types=['gamma','proton','gamma-diffuse']
    if source_type == 'g':
        types=['gamma']
    if source_type == 'p':
        types=['proton']
    if source_type == 'd':
        types=['gamma-diffuse'] 
        
    try:
        if chop != None:
            choppoints = [int(i) for i in chop.split(',')]
            if len(choppoints) == 3:
                id_no = choppoints[2]
        else:
            choppoints = None
    except:
        choppoints = None
        
    try:
        site_location = [float(i) for i in location.split(',')]
    except:
        site_location = None
           
    return input_path, output_path, types, site_location, choppoints, id_no

def validate(input_path,output_path,max_files,max_events,
site_altitude,source_type,telescopes,location,chop,types):

    # Making sure all the inputs are sensible and won't throw errors
    
    if source_type == None and chop != None:
            return "Error: Cannot use 'chop' when processing 'all' types"

    for typename in types:
        if not os.path.isdir(input_path + typename):
            print(input_path + typename)
            return 'There is no ' + typename + ' directory'

    if telescopes != None:
        try:
            telescopes = [int(i) for i in telescopes.split(',')]
        except:
            return "Error: Some of those telescope IDs aren't integers"
            
    id_no = ""

    if chop != None:
        if max_files != None:
            return "Error: Can't chop while also specifying max_files"
        try:
            choppoints = [int(i) for i in chop.split(',')]
        except:
            return "Error: Invalid chop (should be two integers separated by commas, indicating nth and mth files to process, plus an optional third ID no.)"
        if len(choppoints) < 2:
            return "Error: Not enough chop points!"
        if len(choppoints) > 3:
            return "Error: Too many chop points!"
        if choppoints[0] < 1 or choppoints[1] > len(glob.glob(input_path + typename + '/*.simtel.gz')):
            return "Error: Chop out of file range!"
            
    try:
        site_location = [float(i) for i in location.split(',')]
    except:
        return "Error: Invalid lat/lon/alt input (should be three numbers separated by commas)"
    if len(site_location) < 3:
        return "Error: Not enough location arguments"
    if len(site_location) > 3:
        return "Error: Too many location arguments"
    if not -90 <= site_location[0] <= 90:
        return "Error: Invalid latitude"
    if not -180 <= site_location[1] <= 180:
        return "Error: Invalid longitude"
    if not -450 <= site_location[2] <= 9000:
        return "Error: Invalid height (no longer on land)"
        
    return "Valid"

if __name__ == '__main__':
    main()
