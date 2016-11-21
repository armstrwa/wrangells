'''
Script to analyze best fully automated results, now that I've done a lot of cherry picking..
20 nov 2016
william armstrong
'''

import json
import matplotlib.pyplot as plt
import numpy as np
import wrangellsGeospatialTools as wgt 

jsonFn = '/Users/anderson/Desktop/ARMSTRONG/wrangells/corr/select/filtered/swath/kennCL_swathVelocitySampling_2016-10-03.json'

data = wgt.readJsonFile(jsonFn)