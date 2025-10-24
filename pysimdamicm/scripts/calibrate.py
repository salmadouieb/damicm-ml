"""
Take raw fits image, pedestal subtract and fit dark current.

Necessary inputs:
1. image fits file path as 'path'
2. specify json file name (must be in json folder) as "json"
3. whether or not you want to print the calibration curves

"""


# loading damicm module
import pysimdamicm as ccd
from pysimdamicm import utils, io, processes

# loading other packages
import numpy as np
from matplotlib import pyplot as plt
import ROOT
from astropy.io import fits

file_path = path
# open the file as a hdul
json_name = json

# fresh load of the packaged JSON
cfg_file = f"{ccd.__path__[0]}/json/" + json_name
cfg = ccd.utils.config.Config(cfg_file, simulations=False)

# load the rawdata object
rdata = io.rawdata.BuilderRawData(file_path, cfg.configuration["input"])

rdata.prepare_data()

# compress skips (necessary even though there are no skips)
comp = ccd.processes.skipper_analysis.CompressSkipperProcess()
comp.func_to_compress = ['mean']
comp.id_skip_start = 2
comp.execute_process(rdata)

# initialise pedestal subtraction process
ped = ccd.processes.skipper_analysis.PedestalSubtractionProcess()
