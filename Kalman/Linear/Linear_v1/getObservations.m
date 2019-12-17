function observations = getObservations(trueTracks,parameters)
numSteps = parameters.numSteps;
sigmaMeasurementNoise = parameters.sigmaMeasurementNoise;

observations = trueTracks(1:2,:) + sigmaMeasurementNoise*randn(2,numSteps);

end

