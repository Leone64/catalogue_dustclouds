import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.io import fits
import time

from dustmaps.edenhofer2023 import Edenhofer2023Query
import dustmaps.edenhofer2023
from dustmaps.config import config

import json


import FuncDef_Simon as fd

from tqdm import tqdm

#Download dustmaps data and where to store it
config["data_dir"] = "/home/leone/Documents/Uni/Bachelor/DustmapsData"
# dustmaps.edenhofer2023.fetch()

working_folder = "/home/leone/Documents/Uni/Bachelor/TestScripts/Checkpoint/"

#load cloud data
data = fd.load_data(working_folder + "cloud-data.csv", sep=";")

#setup empty dataframe to store results
result = pd.DataFrame(columns=["i", "mean", "std", "staterr"])

start_time_total = time.time()
total_clouds = len(data.index)

for idx, i in enumerate(tqdm(data.index), start=4):
    if i == 5:
        break
    cloud_start = time.time()

    #this is the long one
    masses = fd.massCalc(*data.loc[i], var=True, folder=working_folder, cloud_id=i)

    # Total ETA calculation
    elapsed_total = time.time() - start_time_total
    completed_clouds = idx + 1
    avg_time = elapsed_total / completed_clouds
    remaining_clouds = total_clouds - completed_clouds
    total_eta = avg_time * remaining_clouds

    # Append total ETA to progress.txt
    with open(working_folder + "progress.txt", "a") as f:
        f.write(
            f"\n[Total] Completed: {completed_clouds}/{total_clouds} | "
            f"Elapsed: {elapsed_total/60:.1f}m | ETA: {total_eta/60:.1f}m\n\n"
        )

    #calc mean, std, and statistical error
    mean = np.nanmean(list(masses.values()))
    std = np.nanstd(list(masses.values()))
    cloud = fd.Cloud(*data.loc[i])
    staterr = cloud.getMassStatError()

    #save/checkpoint results
    with open(working_folder + "masses_values.txt", "a") as f:
        f.write(list(masses.values()))
    result.loc[i] = [i, mean, std, staterr]
    result.to_csv(working_folder + "results.csv")

    with open(working_folder + "results.txt", "a") as f:
        f.write(f"{i}: {mean} +/- {std} ({staterr})\n")

#save final results 
result.to_csv(working_folder + "results.csv")
