#!/g/data/hh5/public/apps/miniconda3/envs/analysis3-22.10/bin/python3

import os
import mule
from netCDF4 import Dataset
import numpy as np
import cftime
import datetime

stashm = mule.STASHmaster.from_version('7.3')

for filename in ( "n48.veg.func_igbp.shiftedAusNZ.nc", "n48.veg.func_seas.shiftedAusNZ.nc" ):
### Load netcdf
    dataset = Dataset(filename)
    missing_data_value = -1073741824.0

    ### Account for differences within the file
    if filename == "n48.veg.func_seas.shiftedAusNZ.nc":
        z_var = 'surface'
        time_type = 2 ### periodic time series
    elif filename == "n48.veg.func_igbp.shiftedAusNZ.nc":
        z_var = 'pseudo'
        time_type = 0 ### single valid time

    ### Use its contents to generate a template
    template = {}
    template['integer_constants'] = {}
    template['integer_constants']['num_times'] = dataset.dimensions['t'].size
    template['integer_constants']['num_levels'] = dataset.dimensions[z_var].size
    template['integer_constants']['num_cols'] = dataset.dimensions['longitude'].size
    template['integer_constants']['num_rows'] = dataset.dimensions['latitude'].size

    template['real_constants'] = {}
    template['real_constants']['col_spacing'] = 360.0 / dataset.dimensions['longitude'].size
    template['real_constants']['row_spacing'] = 180.0 / (dataset.dimensions['latitude'].size-1)
    template['real_constants']['start_lat'] = -90.0
    template['real_constants']['start_lon'] = 0.0
    template['real_constants']['north_pole_lat'] = 90.0
    template['real_constants']['north_pole_lon'] = 0.0


    template['fixed_length_header'] = {}
    ### Valid times for seasonal file
    if filename == "n48.veg.func_seas.shiftedAusNZ.nc":

        start_day = 16
        start_time =  cftime.datetime.strptime(dataset.variables['t'].time_origin,'%d-%b-%Y:%H:%M:%S',calendar=dataset.variables['t'].calendar).replace(day=start_day)

        template['fixed_length_header']['t1_year'] = start_time.year
        template['fixed_length_header']['t1_month'] = start_time.month
        template['fixed_length_header']['t1_day'] = start_time.day
        template['fixed_length_header']['t1_hour'] = start_time.hour
        template['fixed_length_header']['t1_minute'] = start_time.minute
        template['fixed_length_header']['t1_second'] = 0
        template['fixed_length_header']['t1_year_day_number'] = start_day


        end_time = (start_time + datetime.timedelta(days=dataset['t'][-1].data.item())).replace(day=start_day)
        template['fixed_length_header']['t2_year'] = end_time.year
        template['fixed_length_header']['t2_month'] = end_time.month
        template['fixed_length_header']['t2_day'] = end_time.day
        template['fixed_length_header']['t2_hour'] = end_time.hour
        template['fixed_length_header']['t2_minute'] = end_time.minute
        template['fixed_length_header']['t2_second'] = 0
        template['fixed_length_header']['t2_year_day_number'] = (end_time - start_time).days + start_day

        ### Interval between valid times
        template['fixed_length_header']['t3_year'] = 0
        template['fixed_length_header']['t3_month'] = 1
        template['fixed_length_header']['t3_day'] = 0
        template['fixed_length_header']['t3_hour'] = 0
        template['fixed_length_header']['t3_minute'] = 0
        template['fixed_length_header']['t3_second'] = 0
        template['fixed_length_header']['t3_year_day_number'] = 0

        ### Other important values for the seasonal file
        template['fixed_length_header']['total_prognostic_fields'] = 180 ### 3*60
        template['integer_constants']['num_field_types'] = 15

    ### Set all of these to the same value for the other file
    elif filename == "n48.veg.func_igbp.shiftedAusNZ.nc":

        start_time =  cftime.datetime.strptime(dataset.variables['t'].time_origin,'%d-%b-%Y:%H:%M:%S',calendar=dataset.variables['t'].calendar).replace(day=1)

        template['fixed_length_header']['t1_year'] = start_time.year
        template['fixed_length_header']['t1_month'] = start_time.month
        template['fixed_length_header']['t1_day'] = start_time.day
        template['fixed_length_header']['t1_hour'] = start_time.hour
        template['fixed_length_header']['t1_minute'] = start_time.minute
        template['fixed_length_header']['t1_second'] = 0
        template['fixed_length_header']['t1_year_day_number'] = 0

        template['fixed_length_header']['t2_year'] = start_time.year
        template['fixed_length_header']['t2_month'] = start_time.month
        template['fixed_length_header']['t2_day'] = start_time.day
        template['fixed_length_header']['t2_hour'] = start_time.hour
        template['fixed_length_header']['t2_minute'] = start_time.minute
        template['fixed_length_header']['t2_second'] = 0
        template['fixed_length_header']['t2_year_day_number'] = 0

        template['fixed_length_header']['t3_year'] = 0
        template['fixed_length_header']['t3_month'] = 0
        template['fixed_length_header']['t3_day'] = 0
        template['fixed_length_header']['t3_hour'] = 0
        template['fixed_length_header']['t3_minute'] = 0
        template['fixed_length_header']['t3_second'] = 0
        template['fixed_length_header']['t3_year_day_number'] = 0

        ### Interval between valid times
        template['fixed_length_header']['t3_year'] = 0
        template['fixed_length_header']['t3_month'] = 1
        template['fixed_length_header']['t3_day'] = 0
        template['fixed_length_header']['t3_hour'] = 0
        template['fixed_length_header']['t3_minute'] = 0
        template['fixed_length_header']['t3_second'] = 0
        template['fixed_length_header']['t3_year_day_number'] = 0

        ### Other important values for the seasonal file
        template['fixed_length_header']['total_prognostic_fields'] = 180 ### 3*60
        template['integer_constants']['num_field_types'] = 15

    ### Other important settings
    template['fixed_length_header']['data_set_format_version'] = 15 ### MASS storage
    template['fixed_length_header']['sub_model'] = 1 ### Atmosphere
    template['fixed_length_header']['vert_coord_type'] = 1 ### Matches 
    template['fixed_length_header']['dataset_type'] = 4 ### Ancillary dataset
    template['fixed_length_header']['grid_staggering'] = 3 ### ENDGame
    template['fixed_length_header']['time_type'] = time_type
    template['fixed_length_header']['model_version'] = 703 ### Based on UM7.3 STASHMaster
    template['fixed_length_header']['horiz_grid_type'] = 0 ### Global grid
    template['fixed_length_header']['calendar'] = 2 ### 360 day calendar

    ### create the ancil file
    new_ancil = mule.AncilFile.from_template(template)

    coord_vars=('t',z_var,'longitude','latitude','t_1')

    ### Copy netcdf data to new fields
    for var in dataset.variables:
        if var in coord_vars: continue
        print(var)
        stash_id = 0
        for k,v in stashm.items():
            if dataset[var].title == v.name:
                stash_id=int(k)
                break

        
        #for t in range(len(dataset[var])):
        for t in range(len(dataset['t'])):
            try:
                field_3d = dataset[var][t]
            except IndexError:
                ### This must be field1384 - this field needs to repeat 5 times
                field_3d=dataset[var][t%12]
            for pli in range(len(field_3d)):

                field_2d = field_3d[pli]
                out_field_2d = np.where(field_2d.mask,missing_data_value,field_2d)

                new_field = mule.Field2.empty()

                ### Programmatically derive as much metadata as we can
                new_field.lbuser4 = stash_id
                new_field.lbnpt = dataset.dimensions['longitude'].size
                new_field.bzx = -template['real_constants']['col_spacing']
                new_field.bdx = template['real_constants']['col_spacing']
                new_field.lbrow = dataset.dimensions['latitude'].size
                new_field.bzy = -90.0 - template['real_constants']['row_spacing']
                new_field.bdy = template['real_constants']['row_spacing']
                new_field.bmks = 1.0

                new_field.lbfc = int(var[5:])

                ### Derive time bounds
                if isinstance(dataset['t'][t],np.float32):
                    field_start_time = start_time + datetime.timedelta(days=dataset['t'][t].item())
                elif isinstance(dataset['t'][t],np.ma.core.MaskedArray):
                    field_start_time = start_time + datetime.timedelta(days=dataset['t'][t].data.item())
                else:
                    raise TypeError
            
                field_start_time = field_start_time.replace(day=1)
                field_end_time = field_start_time.replace(day=30,hour=23,minute=59)

                print(field_end_time)

                ### Start time for this field
                new_field.lbyr = field_start_time.year
                new_field.lbmon = field_start_time.month
                new_field.lbdat = field_start_time.day
                new_field.lbhr = field_start_time.hour
                new_field.lbmin = field_start_time.minute
                new_field.lbday = 0

                ### End time for this field (only matters for seasonal file)
                new_field.lbyrd = field_end_time.year
                new_field.lbmond = field_end_time.month
                new_field.lbdatd = field_end_time.day
                new_field.lbhrd = field_end_time.hour
                new_field.lbmind = field_end_time.minute
                new_field.lbdayd = 0

                ### Constants derived from existing ancils and F03
                new_field.lbrel = 2 ### UM <vn8.1
                new_field.lbcode = 1 ### Regular lat/lon grid
                new_field.lbhem = 0 ### Global
                new_field.lbpack = 0 ### No packing
                new_field.lbuser1 = 1 ### Real field
                new_field.lbuser3 = 0 ### Boundary datasets - not relevant
                new_field.lbuser5 = int(dataset[z_var][pli].data.item()) ### Pseudo-level for this field
                new_field.lbuser6 = 0 ### Unused
                new_field.lbuser7 = 1 #### Atmosphere
                new_field.lbtim = time_type ### Time indicator
                new_field.lblev = 8888 ### level code - 8888
                new_field.lbproc = 0 ### "None of the above"
                new_field.lbvc = 129 ### Vertical coord type - 129 = surface
                new_field.bmdi = missing_data_value ### Fill value
                new_field.bplat = 90.0 ### Pole latitude
                new_field.lbext = 0 ### No extra data
                new_field.lbft = 0 ### No forecast period
                new_field.lbproj = 900 ### fields file projection number
                new_field.lbuser2 = 1 ### data start
                new_field.lbsrce = 1111 ### If YYYY=1111, Words 39-45 = LBUSER area for UM use.


                ### Various things that need to be 0
                new_field.lbrvc = 0 ### Vertical co-ordinate type for reference level. Always zero.
                new_field.lbcfc = 0 ### LBCFC Second field code (=0 in UM)
                new_field.lbexp = 0 ### Experiment number (optional parameter for users reference only).
                new_field.lbtyp = 0 ### CF Field code
                new_field.lbrsvd1 = 0
                new_field.lbrsvd2 = 0
                new_field.lbrsvd3 = 0
                new_field.lbrsvd4 = 0

                array_provider = mule.ArrayDataProvider(out_field_2d)
                new_field.set_data_provider(array_provider)
                new_ancil.fields.append(new_field)
                
    ### write file
    #new_ancil.to_file(f"{os.path.splitext(filename)[0]}.ancil")
    new_ancil.to_file(f"{filename.split('.')[2]}.ancil")
