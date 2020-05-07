def validate(input_path,output_path,max_files,max_events,site_altitude,source_type,telescopes,location,chop):
    import os
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
            return "Error: Cannot use 'chop' when processing 'all' types")
    if source_type == 'g':
        types=['gamma']
    if source_type == 'p':
        types=['proton']
    if source_type == 'd':
        types=['gamma-diffuse']

    for typename in types:
        if not os.path.isdir(input_path + typename):
            return 'There is no ' + typename + ' directory')

    if telescopes != None:
        try:
            telescopes = [int(i) for i in telescopes.split(',')]
        except:
            return "Error: Some of those telescope IDs aren't integers"
            
    id_no = ""

    if chop != None:
        if max_files != None:
            return "Error: Can't chop while also specifying max_files")
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
        if len(choppoints) == 3:
            id_no = choppoints[2]
            
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
        return "Error: Invalid longitude")
    if not -450 <= site_location[2] <= 9000:
        return "Error: Invalid height (no longer on land)"
        
    return "Valid"
