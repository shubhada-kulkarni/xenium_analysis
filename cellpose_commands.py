# python script for cell segmentation analysis

import numpy as np
import time, os, sys
from urllib.parse import urlparse
import matplotlib
import matplotlib.pyplot as plt
import matplotlib as mpl
#%matplotlib inline
mpl.rcParams['figure.dpi'] = 300
import tifffile
from cellpose import models, io
from cellpose import plot

## read the image with cell boundary marker channel
#file = "/prj/XeniumProbeDesign/heart_human_29072024/output-XETG00046__0018072__Region_1__20240725__112631/morphology_focus/morphology_focus_0001.ome.tif"
## img = cellpose.io.imread(file)
#data = tifffile.imread(file)
#print(data.shape)

# DEFINE CELLPOSE MODEL
# model_type='cyto3' or model_type='nuclei'
model = models.Cellpose(gpu=False, model_type='cyto3')

channels = [[2,3], [0,0], [0,0]]
# or
chan=[2,3]

filename = "/prj/XeniumProbeDesign/heart_human_29072024/output-XETG00046__0018072__Region_1__20240725__112631/morphology_focus/morphology_focus_0001.ome.tif"
img = io.imread(filename)

# running cellpose model evaluation
masks, flows, styles, diams = model.eval(img, diameter=None, channels=chan)      # run cellpose size model, mask model and get masks
io.masks_flows_to_seg(img, masks, flows, filename, channels=chan, diams=diams)   # saves the npy file
io.save_to_png(img, masks, flows, filename)

# plotting the results
fig = plt.figure(figsize=(12,12))
plot.show_segmentation(fig, img, masks, flows[0], channels=chan)
plt.tight_layout()
plt.savefig("segmented.png")