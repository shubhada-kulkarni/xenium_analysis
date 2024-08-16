#!/usr/bin/env python

# python script to separate the main morphology image of Xenium into separate levels
# so that cellpose can process this

import tifffile
import sys, os

tiff = sys.argv[1]
outpath = os.path.dirname(tiff)

# Variable 'LEVEL' determines the level to extract. It ranges from 0 (highest
# resolution) to 6 (lowest resolution) for morphology.ome.tif
# LEVEL = 6   # highest resolution 0 did not work, so gave 1
level_array = [1, 2, 3, 4, 5, 6]


for LEVEL in level_array:
    print(LEVEL)

    with tifffile.TiffFile(tiff) as tif:
        image = tif.series[0].levels[LEVEL].asarray()
    
    tifffile.imwrite(outpath+'/cellpose_level_'+str(LEVEL)+'_morphology.ome.tif',
        image,
        photometric='minisblack',
        dtype='uint16',
        tile=(1024, 1024),
        compression='JPEG_2000_LOSSY',
        metadata={'axes': 'ZYX'},
    )


