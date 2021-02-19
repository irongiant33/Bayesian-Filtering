import numpy as np
import math
import os

import plot
import trueTrack
import observations
import error
import kalmanFilter


class Parameters:
    def __init__(self):
        self.num_steps = 200
        self.scan_time = 1
        self.sigma_driving_noise = .1
        self.sigma_measurement_noise = 1
        self.start_state = [0, 0, 0, 0]
        self.prior_covariance = list(np.diagflat([100, 100, 1, 1]))
        random_covariance = math.sqrt(self.prior_covariance) * \
            list(np.random.rand(4,1))
        self.priorMean = self.start_state + random_covariance

def ask_for_parameters():
    path_sym = "/"
    if os.name == "nt":
        path_sym = "\\"
    f = []
    mypath = f"..{path_sym}parameters{path_sym}"
    for (_, _, filenames) in os.walk(mypath):
        parameter_files = [files if files[-3:] == "txt" for files in filenames]
        f.extend(parameter_files)
        break
    if len(f) != 0:
        

def main():
    parameters = ask_for_parameters()
    true_tracks = get_true_track(parameters)
    observations = get_observations(true_tracks, parameters)
    observed_rmse = get_error(trueTracks, observations)
    estimated_tracks = kalman_filter(observations, parameters)
    estimated_rmse = get_error(trueTracks, estimatedTracks)
    plot(observations, observed_rmse, estimated_tracks, estimated_rmse)

if __name__ == "__main__":
    main()