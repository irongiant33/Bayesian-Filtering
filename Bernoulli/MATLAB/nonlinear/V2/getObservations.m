function observations = getObservations(trueTracks,parameters)
numSteps = parameters.numSteps;
sigmaMeasurementNoiseRange = parameters.sigmaMeasurementNoiseRange;
sigmaMeasurementNoiseBearing = parameters.sigmaMeasurementNoiseBearing;
sensorPosition = parameters.sensorPosition;

%we only observe range and bearing but we need to estimate xPos and yPos.
observations = zeros(2,numSteps);
observations(1,:) = sqrt((trueTracks(1,:) - sensorPosition(1)).^2+(trueTracks(2,:) - sensorPosition(2)).^2) + sigmaMeasurementNoiseRange*randn(1,numSteps);
observations(2,:) =  atan2d(trueTracks(1,:) - sensorPosition(1), trueTracks(2,:) - sensorPosition(2)) + sigmaMeasurementNoiseBearing*randn(1,numSteps);
end
