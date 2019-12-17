function observations = getObservations(trueTracks,parameters)
numDevices = parameters.numDevices;
numSteps = parameters.numSteps*numDevices;
sigmaMeasurementNoiseRange = parameters.sigmaMeasurementNoiseRange;
sigmaMeasurementNoiseBearing = parameters.sigmaMeasurementNoiseBearing;
sensorPosition = parameters.sensorPosition;

%we only observe range and bearing but we need to estimate xPos and yPos.
observations = zeros(2,numSteps);
for i = 1:numDevices
    observations(1,i:numDevices:numSteps) = sqrt((trueTracks(1,:) - sensorPosition(1)).^2+(trueTracks(2,:) - sensorPosition(2)).^2) + sigmaMeasurementNoiseRange*randn(1,numSteps/numDevices);
    observations(2,i:numDevices:numSteps) =  atan2d(trueTracks(1,:) - sensorPosition(1), trueTracks(2,:) - sensorPosition(2)) + sigmaMeasurementNoiseBearing*randn(1,numSteps/numDevices);
end
end
