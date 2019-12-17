function observations = getObservations(trueTracks,parameters)
numDevices = parameters.numDevices;
networkMode = parameters.networkMode;
numSteps = parameters.numSteps*numDevices;
sigmaMeasurementNoiseRange = parameters.sigmaMeasurementNoiseRange;
sigmaMeasurementNoiseVelocity = parameters.sigmaMeasurementNoiseVelocity;
sigmaMeasurementNoiseBearing = parameters.sigmaMeasurementNoiseBearing;
sensorPosition = parameters.sensorPosition;
shareMode = 0;
if(numDevices > 1)
    shareMode = 1;
end

%we only observe range and bearing but we need to estimate xPos and yPos.
if(networkMode > 0)
    for i = 1:numDevices
        observations(:,:,i) = zeros(3,numSteps);
    end
else
    observations = zeros(3,numSteps);
end

%only the monostatic case gets a velocity measurement
for k = 1:(numDevices^networkMode)
    observations(1,k:numDevices^shareMode:numSteps,k) = sqrt((trueTracks(1,k:numDevices^shareMode:numSteps) - sensorPosition(k,1)).^2+(trueTracks(2,k:numDevices^shareMode:numSteps) - sensorPosition(k,2)).^2) + sigmaMeasurementNoiseRange*randn(1,numSteps/(numDevices^shareMode));
    observations(2,k:numDevices^shareMode:numSteps,k) = atan2d(trueTracks(1,k:numDevices^shareMode:numSteps) - sensorPosition(k,1), trueTracks(2,k:numDevices^shareMode:numSteps) - sensorPosition(k,2)) + sigmaMeasurementNoiseBearing*randn(1,numSteps/(numDevices^shareMode));
    observations(3,k:numDevices^shareMode:numSteps,k) = 0.5*(((trueTracks(1,k:numDevices^shareMode:numSteps) - sensorPosition(k,1)).^2+(trueTracks(2,k:numDevices^shareMode:numSteps) - sensorPosition(k,2)).^2).^(-0.5)).*(2*(trueTracks(1,k:numDevices^shareMode:numSteps) - sensorPosition(k,1)).*trueTracks(3,k:numDevices^shareMode:numSteps)+2*(trueTracks(2,k:numDevices^shareMode:numSteps) - sensorPosition(k,2)).*trueTracks(4,k:numDevices^shareMode:numSteps)) + sigmaMeasurementNoiseVelocity*randn(1,numSteps/(numDevices^shareMode));

    for i = 1:numDevices
        if(i ~= k)
            observations(1,i:numDevices^shareMode:numSteps,k) = sqrt((sensorPosition(i,1)-trueTracks(1,i:numDevices^shareMode:numSteps)).^2 + (sensorPosition(i,2)-trueTracks(2,i:numDevices^shareMode:numSteps)).^2) + sqrt((sensorPosition(k,1)-trueTracks(1,i:numDevices^shareMode:numSteps)).^2 + (sensorPosition(k,2)-trueTracks(2,i:numDevices^shareMode:numSteps)).^2) + sigmaMeasurementNoiseRange*randn(1,numSteps/(numDevices^shareMode));
            observations(2,i:numDevices^shareMode:numSteps,k) =  atan2d(trueTracks(1,i:numDevices^shareMode:numSteps) - sensorPosition(k,1), trueTracks(2,i:numDevices^shareMode:numSteps) - sensorPosition(k,2)) + sigmaMeasurementNoiseBearing*randn(1,numSteps/(numDevices^shareMode));
        end
    end
end
end
