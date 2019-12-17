function observations = getObservations(trueTracks,parameters)
numDevices = parameters.numDevices;
numSteps = parameters.numSteps*numDevices;
sigmaMeasurementNoiseRange = parameters.sigmaMeasurementNoiseRange;
sigmaMeasurementNoiseVelocity = parameters.sigmaMeasurementNoiseVelocity;
sigmaMeasurementNoiseBearing = parameters.sigmaMeasurementNoiseBearing;
sensorPosition = parameters.sensorPosition; 

%we only observe range and bearing but we need to estimate xPos and yPos.
observations = zeros(3,numSteps);

%only the monostatic case gets a velocity measurement
observations(1,1:numDevices:numSteps) = sqrt((trueTracks(1,:) - sensorPosition(1,1)).^2+(trueTracks(2,:) - sensorPosition(1,2)).^2) + sigmaMeasurementNoiseRange*randn(1,numSteps/numDevices);
observations(3,1:numDevices:numSteps) = 0.5*(((trueTracks(1,:) - sensorPosition(1,1)).^2+(trueTracks(2,:) - sensorPosition(1,2)).^2).^(-0.5)).*(2*(trueTracks(1,:) - sensorPosition(1,1)).*trueTracks(3,:)+2*(trueTracks(2,:) - sensorPosition(1,2)).*trueTracks(4,:)) + sigmaMeasurementNoiseVelocity*randn(1,numSteps/numDevices);
observations(2,1:numDevices:numSteps) =  atan2d(trueTracks(1,:) - sensorPosition(1,1), trueTracks(2,:) - sensorPosition(1,2)) + sigmaMeasurementNoiseBearing*randn(1,numSteps/numDevices);

for i = 2:numDevices
    observations(1,i:numDevices:numSteps) = sqrt((sensorPosition(i,1)-trueTracks(1,:)).^2 + (sensorPosition(i,2)-trueTracks(2,:)).^2) + sqrt((sensorPosition(1,1)-trueTracks(1,:)).^2 + (sensorPosition(1,2)-trueTracks(2,:)).^2) + sigmaMeasurementNoiseRange*randn(1,numSteps/numDevices);
    observations(2,i:numDevices:numSteps) =  atan2d(trueTracks(1,:) - sensorPosition(i,1), trueTracks(2,:) - sensorPosition(i,2)) + sigmaMeasurementNoiseBearing*randn(1,numSteps/numDevices);
end
end
