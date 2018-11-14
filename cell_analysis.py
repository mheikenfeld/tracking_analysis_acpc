import logging

def extract_cell_cubes_subset(cubelist_in,mask,track,cell,
                              z_coord='model_level_number', height_levels=None):
    
    from iris.analysis import SUM
    from iris import Constraint
    from iris.cube import CubeList
    from iris.coords import AuxCoord
    import numpy as np
    from tobac import mask_cell,mask_cell_surface,get_bounding_box
    from copy import deepcopy
    
    track_i=track[track['cell']==cell]
    
    cubelist_cell_integrated_out=CubeList()
    cubelist_cell_sum=CubeList()

    for time_i in track_i['time'].values:    
        
        logging.debug('start extracting cubes for cell '+str(cell) + ' and time '+str(time_i))

        constraint_time = Constraint(time=time_i)
        mask_i=mask.extract(constraint_time)
        mask_cell_i=mask_cell(mask_i,cell,track_i,masked=False)
        mask_cell_surface_i=mask_cell_surface(mask_i,cell,track_i,masked=False,z_coord=z_coord)

        x_dim=mask_cell_surface_i.coord_dims('projection_x_coordinate')[0]
        y_dim=mask_cell_surface_i.coord_dims('projection_y_coordinate')[0]
        x_coord=mask_cell_surface_i.coord('projection_x_coordinate')
        y_coord=mask_cell_surface_i.coord('projection_y_coordinate')
    
        if (mask_cell_surface_i.core_data()>0).any():
            box_mask_i=get_bounding_box(mask_cell_surface_i.core_data(),buffer=1)
    
            box_mask=[[x_coord.points[box_mask_i[x_dim][0]],x_coord.points[box_mask_i[x_dim][1]]],
                     [y_coord.points[box_mask_i[y_dim][0]],y_coord.points[box_mask_i[y_dim][1]]]]
        else:
            box_mask=[[np.nan,np.nan],[np.nan,np.nan]]
    
        width=20
        dx=500
        x=track_i[track_i['time'].values==time_i]['projection_x_coordinate'].values[0]
        y=track_i[track_i['time'].values==time_i]['projection_y_coordinate'].values[0]

        n_add_width=2
        
        box_slice=[[x-(width+n_add_width)*dx,x+(width+n_add_width)*dx],[y-(width+n_add_width)*dx,y+(width+n_add_width)*dx]]
               
        x_min=np.nanmin([box_mask[0][0],box_slice[0][0]])
        x_max=np.nanmax([box_mask[0][1],box_slice[0][1]])
        y_min=np.nanmin([box_mask[1][0],box_slice[1][0]])
        y_max=np.nanmax([box_mask[1][1],box_slice[1][1]])

        constraint_x=Constraint(projection_x_coordinate=lambda cell: int(x_min) < cell < int(x_max))
        constraint_y=Constraint(projection_y_coordinate=lambda cell: int(y_min) < cell < int(y_max))

        constraint=constraint_time & constraint_x & constraint_y
    
        mask_cell_i=mask_cell_i.extract(constraint_x & constraint_y)
        mask_cell_surface_i=mask_cell_surface_i.extract(constraint_x & constraint_y)
           
        cubelist_i=cubelist_in.extract(constraint)
       
        cubelist_cell_sum.extend(sum_profile_mask(cubelist_i,height_levels,mask_cell_i))
    cubelist_cell_sum_out=cubelist_cell_sum.merge()        
    for cube in cubelist_cell_sum_out:
        cell_time_coord=AuxCoord(track_i['time_cell'].dt.total_seconds().values, units='s',long_name='time_cell')
        cube.add_aux_coord(cell_time_coord,cube.coord_dims('time')[0])

    for cube in cubelist_cell_sum_out:
        cubelist_cell_integrated_out.append(cube.collapsed(('geopotential_height'),SUM))
            
    track_cell_integrated=deepcopy(track_i)
#
    for cube in cubelist_cell_integrated_out:
        track_cell_integrated[cube.name()]=cube.core_data()
        
    return cubelist_cell_sum_out, cubelist_cell_integrated_out,track_cell_integrated


def collapse_profile_mask(variable_cube,height_levels_borders,mask_cell,coordinate='geopotential_height',method=None,fillnan=None):
    from iris.coords import DimCoord
    from iris.cube import Cube
    from iris import Constraint
    from numpy import array,transpose
    import numpy as np
    from tobac import mask_cube

    # Set midpoints of height bounds as coordinate values
    height_levels=height_levels_borders[0:-1]+0.5*np.diff(height_levels_borders)
    # Create array of height bounds
    bounds=transpose(array([height_levels_borders[0:-1],height_levels_borders[1:]]))
    #Create DimCoord of height levels and add coordinates
    dim_coordinate=DimCoord(height_levels,standard_name=coordinate,units=variable_cube.coord(coordinate).units,bounds=bounds)
    variable_cube_out=Cube(height_levels,units=variable_cube.units,long_name=variable_cube.name(),dim_coords_and_dims=[(dim_coordinate,0)])
    variable_cube_out.add_aux_coord(variable_cube.coord('time'))    
    for i_height,height in enumerate(height_levels_borders[:-1]):
        constraint_height=Constraint(geopotential_height= lambda cell: height_levels_borders[i_height]<= cell <height_levels_borders[i_height+1])
        mask=mask_cell
        variable_cube_i=mask_cube(variable_cube,mask.core_data()).extract(constraint_height)        
        variable_cube_i=variable_cube_i.collapsed(('model_level_number','x','y'),method)
        variable_cube_out_i=variable_cube_i.core_data()
        if fillnan is not None:
            variable_cube_out_i=fillnan
        variable_cube_out.data[i_height]=variable_cube_out_i
    return variable_cube_out
    
def sum_profile_mask(cubelist_in,height_levels_borders,mask_cell):    
    from iris.cube import CubeList
    from iris.analysis import SUM
    from dask.array.ma import masked_invalid
    cubelist_out = CubeList()
    for variable in cubelist_in:
        #Sum up values for height slices and mask for all cubes in cubelist_in:
        cube_cell_profile=collapse_profile_mask(variable_cube=variable,height_levels_borders=height_levels_borders,mask_cell=mask_cell,
                                             coordinate='geopotential_height',method=SUM)
        cube_cell_profile.rename(variable.name())
        cube_cell_profile.data=masked_invalid(cube_cell_profile.core_data())
        cubelist_out.append(cube_cell_profile)
    return cubelist_out

