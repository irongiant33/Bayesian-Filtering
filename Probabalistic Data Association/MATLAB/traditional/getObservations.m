function observations = getObservations(trueTracks,parameters)
numSteps = parameters.numSteps;
sigmaMeasurementNoise = parameters.sigmaMeasurementNoise;

observations = trueTracks + sigmaMeasurementNoise*randn(4,numSteps);

end

