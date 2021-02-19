import numpy as np
import math
import os
import copy
import json

import plot
import trueTrack
import observations
import error
import kalmanFilter


class Parameters:
    def __init__(self):
        self._num_steps = 200
        self._scan_time = 1
        self._sigma_driving_noise = .1
        self._sigma_measurement_noise = 1
        self._start_state = np.array([0, 0, 0, 0])
        self._prior_covariance = np.diagflat([100, 100, 1, 1])
        random_covariance = np.sqrt(self._prior_covariance) * \
            np.random.rand(len(self._start_state), 1)
        self._priorMean = self._start_state + random_covariance
        self.param_dict = {
            "Number of Steps":self._num_steps,
            "Scan Time (s)":self._scan_time,
            "Driving Noise Standard Deviation":self._sigma_driving_noise,
            "Measurement Noise Standard Deviation":self._sigma_measurement_noise,
            "Start State (m)":self._start_state,
            "Prior Covariance":self._prior_covariance
        }

    def set(self, field, value):
        if field == "Number of Steps":
            try:
                value = int(value)
                self.num_steps = value
                return True
            except:
                return False
        elif field == "Scan Time (s)":
            try:
                value = float(value)
                self.scan_time = value
                return True
            except:
                return False
        elif field == "Driving Noise Standard Deviation":
            try:
                value = float(value)
                self.sigma_driving_noise = value
                return True
            except:
                return False
        elif field == "Measurement Noise Standard Deviation":
            try:
                value = float(value)
                self.sigma_measurement_noise = value
                return True
            except:
                return False
        elif field == "Start State (m)":
            try:
                value = value.split(',')
                value = [float(val) for val in value]
                self.prior_covariance = value
                return True
            except:
                return False
        elif field == "Prior Covariance":
            try:
                value = value.split(',')
                value = [float(val) for val in value]
                self.prior_covariance = value
                return True
            except:
                return False
        else:
            return False
        

    @staticmethod
    def _is_number(value):
        return (isinstance(value, int) or isinstance(value, float))

    @property
    def num_steps(self):
        return self._num_steps

    @num_steps.setter
    def num_steps(self, value):
        if isinstance(value, int):
            if value <= 0:
                raise ValueError
            self._num_steps = value
        raise ValueError

    @property
    def scan_time(self):
        return self._scan_time

    @scan_time.setter
    def scan_time(self, value):
        if self._is_number(value):
            if value <= 0:
                raise ValueError
            self._scan_time = value
        else:
            raise ValueError

    @property
    def sigma_driving_noise(self):
        return self._sigma_driving_noise

    @sigma_driving_noise.setter
    def sigma_driving_noise(self, value):
        if self._is_number(value):
            if value <= 0:
                raise ValueError
            self._sigma_driving_noise = value
        else:
            raise ValueError

    @property
    def sigma_measurement_noise(self):
        return self._sigma_measurement_noise

    @sigma_measurement_noise.setter
    def sigma_measurement_noise(self, value):
        if self._is_number(value):
            if value <= 0:
                raise ValueError
            self._sigma_measurement_noise = value
        else:
            raise ValueError

    @property
    def start_state(self):
        return self._start_state

    @start_state.setter
    def start_state(self, value):
        if isinstance(value, list):
            for val in value:
                if not self._is_number(val):
                    raise ValueError
            self._start_state = copy.deepcopy(np.array(value))
        else:
            raise ValueError

    @property
    def prior_covariance(self):
        return self._prior_covariance

    @prior_covariance.setter
    def prior_covariance(self, value):
        if((isinstance(value, list)) and (len(value)==len(self._start_state))):
            for val in value:
                if not self._is_number(val):
                    raise ValueError
            cov_matrix = np.diagflat([100, 100, 1, 1])
            self._start_state = copy.deepcopy(cov_matrix)
        else:
            return False

    @property
    def prior_mean(self):
        random_covariance = np.sqrt(self._prior_covariance) * \
            np.random.rand(len(self._start_state), 1)
        return self._start_state + random_covariance 

def load_saved_parameters(filename):
    saved_params = Parameters()
    with open(filename, "r") as f:
        params = json.load(f)
        for param in params:
            valid = saved_params.set(param, params[param])
            if not valid:
                return False
    return saved_params

def ask_for_parameters():
    path_sym = "/"
    if os.name == "nt":
        path_sym = "\\"
    f = []
    mypath = f"..{path_sym}parameters{path_sym}"
    for (_, _, filenames) in os.walk(mypath):
        parameter_files = []
        [parameter_files.append(i) if i[-5:] == ".json" else 0 for i in filenames]
        f.extend(parameter_files)
        for k in range(1, 1 + len(parameter_files)):
            print(f"({k}) {parameter_files[k - 1]}")
        break
    if len(f) != 0:
        user_in = input("If you would like to load a previous configuration,"+\
            "enter \n the number here. Otherwise, enter any other character: ")
        if(user_in.isnumeric() and int(user_in) <= len(parameter_files)):
            file_num = int(user_in) - 1
            valid = load_saved_parameters(mypath + parameter_files[file_num])
            if isinstance(valid, bool):
                print("Bad file parameters.")
            else:
                return valid
    new_params = Parameters()
    print("Enter your desired parameters in the following fields.\n")
    print("Ensure values are appropriate and the start state and prior \n" + \
        "covariance are comma separated values with the same dimensionality.")
    for field in new_params.param_dict:
        valid = False
        while(not valid):
            user_in = input(f"{field}: ")
            try:
                valid = new_params.set(field, user_in)
            except:
                print("Enter a valid value")
    save_name = input("Enter a filename to save the parameters: ")
    with open(f"{mypath}{save_name}.json", "w") as f:
        json.dump(new_params.param_dict, f)
    return new_params


def main():
    parameters = ask_for_parameters()
    #true_tracks = get_true_track(parameters)
    #observations = get_observations(true_tracks, parameters)
    #observed_rmse = get_error(trueTracks, observations)
    #estimated_tracks = kalman_filter(observations, parameters)
    #estimated_rmse = get_error(trueTracks, estimatedTracks)
    #plot(observations, observed_rmse, estimated_tracks, estimated_rmse)

if __name__ == "__main__":
    main()