def add_geopotential_height(cube):    
    from acpc_intercomparison.make_geopotential_height_coord import geopotential_height_coord,geopotential_height_coord_stag
    coord_names=[coord.name() for coord in cube.coords()]
    if ('model_level_number' in coord_names) and ('geopotential_height' not in coord_names):
        if cube.coord('model_level_number').shape[0]==95:
            cube.add_aux_coord(geopotential_height_coord_stag,data_dims=cube.coord_dims('model_level_number'))
        if cube.coord('model_level_number').shape[0]==94:
            cube.add_aux_coord(geopotential_height_coord,data_dims=cube.coord_dims('model_level_number'))
    return cube
