# -*- coding: utf-8 -*-
"""
Created on Wed Nov  4 16:07:28 2020

@author: Julius Kricheldorff MSc.

@affilation: University of Oldenburg

@contact: julius@kricheldorff.de
"""
import random
import pandas as pd
import numpy as np
from time import time
from math import floor

# Create a dataseries with the seeds we are using for the randomization
PilotSeeds = pd.DataFrame({'Seed': [None]*30, 
                               'Used': ['No']*30})
for i in range(30):
    PilotSeeds.Seed[i]= random.randrange(2000)
# Save the data frame as csv
PilotSeeds.to_csv(r'C:\Users\juliu\OneDrive\Dokumente\PhD-Thesis\Teaching\Adaptive Control\Experiment\AdaptiveControl\OS_Experiment\Lists\PilotSeeds.csv', header=True)    

 
for i in range(5,30):
    random.seed(PilotSeeds.Seed[i])
    # Workspace
    block = BlockListRandomizer() 
    LIn1, LIn2 = block.create_LWPC_sets('MI')
    LCon1, LCon2 = block.create_LWPC_sets('MC')
    I1, I2, I3, I4 = block.create_IWPC_Blocks('size')
    LCon1['Block']= 3
    LCon2['Block']= 4
    I1['Block']= 5
    I2['Block']= 6
    I3['Block']= 7
    I4['Block']= 8
    
    # import repair functions
    repair = List_repair()
    
    # Find the minimum index per list
    indices_LIn1 = repair.constrain_checks(LIn1, 0)
    indices_LIn2 = repair.constrain_checks(LIn2, 0)
    indices_LCon1 = repair.constrain_checks(LCon1, 0)
    indices_LCon2 = repair.constrain_checks(LCon2, 0)
    indices_I1 = repair.constrain_checks(I1, 0)
    indices_I2 = repair.constrain_checks(I2, 0)
    indices_I3 = repair.constrain_checks(I3, 0)
    indices_I4 = repair.constrain_checks(I4, 0)
    
    
    # repair the data frame
    
    fina_LIn1 = repair.repair(indices_LIn1,  LIn1)
    fina_LIn2 = repair.repair(indices_LIn2,  LIn2)
    fina_LCon1 = repair.repair(indices_LCon1, LCon1)
    fina_LCon2 = repair.repair(indices_LCon2, LCon2)
    fina_I1 = repair.repair(indices_I1,  I1)
    fina_I2 = repair.repair(indices_I2,  I2)
    fina_I3 = repair.repair(indices_I3,  I3)
    fina_I4 = repair.repair(indices_I4,  I4)
    
    frames = [fina_LIn1, fina_LIn2, fina_LCon1, fina_LCon2, fina_I1, fina_I2, fina_I3, fina_I4]
    df = pd.concat(frames)
    
    # save file
    filename = r'C:\Users\juliu\OneDrive\Dokumente\PhD-Thesis\Teaching\Adaptive Control\Experiment\AdaptiveControl\OS_Experiment\Lists\PilotList' + str(i) +'.csv'
    df.to_csv(filename, header=True)
    
data = pd.read_csv(r'C:\Users\juliu\OneDrive\Dokumente\PhD-Thesis\Teaching\Adaptive Control\Experiment\AdaptiveControl\OS_Experiment\Lists\PilotSeeds.csv')
data.Used[0]=='Yes'