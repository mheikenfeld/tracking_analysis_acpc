import iris
import numpy as np

geopotential_height_array=np.array([-24.0999279 ,    24.6950779 ,    75.92983246,   129.72631836,
         186.21264648,   245.52326965,   307.79940796,   373.18936157,
         441.84881592,   513.94122314,   589.63830566,   669.12017822,
         752.57617188,   840.20489502,   932.21508789,  1028.82580566,
        1130.26696777,  1236.78015137,  1348.61901855,  1466.04980469,
        1589.35205078,  1718.81933594,  1854.76000977,  1997.49768066,
        2147.37231445,  2304.74047852,  2469.97705078,  2643.47558594,
        2825.64892578,  3016.93066406,  3217.77685547,  3428.66503906,
        3650.09765625,  3882.60180664,  4126.73144531,  4383.06689453,
        4652.21972656,  4935.07373047,  5230.72851562,  5531.390625  ,
        5831.390625  ,  6131.390625  ,  6431.390625  ,  6731.390625  ,
        7031.390625  ,  7331.390625  ,  7631.390625  ,  7931.390625  ,
        8231.390625  ,  8531.390625  ,  8831.390625  ,  9131.390625  ,
        9431.390625  ,  9731.390625  , 10031.390625  , 10331.390625  ,
       10631.390625  , 10931.390625  , 11231.390625  , 11531.390625  ,
       11831.390625  , 12131.390625  , 12431.390625  , 12731.390625  ,
       13031.390625  , 13331.390625  , 13631.390625  , 13931.390625  ,
       14231.390625  , 14531.390625  , 14831.390625  , 15131.390625  ,
       15431.390625  , 15731.390625  , 16031.390625  , 16331.390625  ,
       16631.390625  , 16931.390625  , 17231.390625  , 17531.390625  ,
       17831.390625  , 18131.390625  , 18431.390625  , 18731.390625  ,
       19031.390625  , 19331.390625  , 19631.390625  , 19931.390625  ,
       20231.390625  , 20531.390625  , 20831.390625  , 21131.390625  ,
       21431.390625  , 21731.390625  , 22031.390625 ])

geopotential_height_coord=iris.coords.AuxCoord(geopotential_height_array[:-1],standard_name='geopotential_height',units='m')
geopotential_height_coord_stag=iris.coords.AuxCoord(geopotential_height_array,standard_name='geopotential_height',units='m')


def add_geopotential_height(cube):
    coord_names=[coord.name() for coord in cube.coords()]
    if ('model_level_number' in coord_names) and ('geopotential_height' not in coord_names):
        if cube.coord('model_level_number').shape[0]==95:
            cube.add_aux_coord(geopotential_height_coord_stag,data_dims=cube.coord_dims('model_level_number'))
        if cube.coord('model_level_number').shape[0]==94:
            cube.add_aux_coord(geopotential_height_coord,data_dims=cube.coord_dims('model_level_number'))
    return cube

