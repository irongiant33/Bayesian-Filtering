function observations = getObservations(trueTracks,parameters)
numSteps = parameters.numSteps;
sigmaMeasurementNoiseRange = parameters.sigmaMeasurementNoiseRange;
sigmaMeasurementNoiseBearing = parameters.sigmaMeasurementNoiseBearing;
sigmaMeasurementNoiseVelocity = parameters.sigmaMeasurementNoiseVelocity;
sensorPosition = parameters.sensorPosition;

%we only observe range and bearing but we need to estimate xPos and yPos.
observations = zeros(3,numSteps);
observations(1,:) = sqrt((trueTracks(1,:) - sensorPosition(1)).^2+(trueTracks(2,:) - sensorPosition(2)).^2) + sigmaMeasurementNoiseRange*randn(1,numSteps);
observations(2,:) =  atan2d(trueTracks(1,:) - sensorPosition(1), trueTracks(2,:) - sensorPosition(2)) + sigmaMeasurementNoiseBearing*randn(1,numSteps);
observations(3,:) = 0.5*(((trueTracks(1,:) - sensorPosition(1)).^2+(trueTracks(2,:) - sensorPosition(2)).^2).^(-0.5)).*(2*(trueTracks(1,:) - sensorPosition(1)).*trueTracks(3,:)+2*(trueTracks(2,:) - sensorPosition(2)).*trueTracks(4,:)) + sigmaMeasurementNoiseVelocity*randn(1,numSteps);
end

