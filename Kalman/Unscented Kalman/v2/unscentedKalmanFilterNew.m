function estimatedTracks = unscentedKalmanFilterNew(observations,parameters)
    %load parameters
    priorMean = parameters.priorMean;
    numSteps = parameters.numSteps;
    scanTime = parameters.scanTime;
    sensorPosition = parameters.sensorPosition;
    rangeVar = parameters.sigmaMeasurementNoiseRange^2;
    angleVar = parameters.sigmaMeasurementNoiseBearing^2;
    drivingVar = parameters.sigmaDrivingNoise^2;
    priorCovariance = parameters.priorCovariance;

    
    %declare variables
    stateSize = length(priorMean);
    measurementCovariance = diag([rangeVar;angleVar]); %we're not measuring velocity
    [stateTransition,noiseGain] = getTransitionMatrices(scanTime); %4x4, 4x2 respectively
    drivingCovariance = noiseGain*diag([drivingVar;drivingVar])*noiseGain'; %4x4 matrix
    estimatedState = zeros(stateSize,numSteps+1);
    estimatedState(:,1) = priorMean;
    estimatedCovariance = priorCovariance;
    
    for k=1:numSteps

        %predict
        predictedMean = stateTransition * estimatedState(:,k);
        predictedCovariance = stateTransition * estimatedCovariance * stateTransition' + drivingCovariance;
        
        %generate sigma points
        [pointsState, weights] = getSigmaPoints(predictedMean, predictedCovariance);
        numPoints = size(pointsState,2);

        %propagate sigma points through observation model
        pointsMeasurements = zeros(2,numPoints); %2x5 matrix
        pointsMeasurements(1,:) = sqrt((pointsState(1,:) - sensorPosition(1)).^2+(pointsState(2,:) - sensorPosition(2)).^2);
        pointsMeasurements(2,:) =  atan2d(pointsState(1,:) - sensorPosition(1),pointsState(2,:) - sensorPosition(2));

        %calculate statistics of transformed state and prediction state
        meanMeasurements = zeros(2,1);
        for i=1:numPoints
            meanMeasurements(1) = meanMeasurements(1) + weights(i) * pointsMeasurements(1,i);
        end
        meanMeasurements(2) = angleFromPoints(pointsMeasurements(2,:)', weights);
        %a major problem in the whole implementation is that 2pi ambiguities have not been considered
        

        covPrediction = zeros(2,2);
        diffMeasurements = zeros(2,numPoints);
        for i=1:numPoints
            diffMeasurements(:,i) = (pointsMeasurements(:,i) - meanMeasurements);
            diffMeasurements(2,i) = wrapTo180(diffMeasurements(2,i));
            covPrediction = covPrediction + weights(i) * (diffMeasurements(:,i) * diffMeasurements(:,i)');
        end
        covPrediction = covPrediction + measurementCovariance; 
        %there was a bug here (measurement covariance it only added after the covarince matrix from sigma points has been fully calculated)

        %calculate cross covariance matrix
        covCross = zeros(stateSize,2);
        for i=1:numPoints
            covCross = covCross + weights(i) * (pointsState(:,i) - predictedMean) * diffMeasurements(:,i)';
        end

        %update with Kalman Filter equations
        kalmanGain = covCross / covPrediction;
        diff = observations(:,k) - meanMeasurements;
        diff(2,1) = wrapTo180(diff(2,1));
        estimatedState(:,k+1) = predictedMean + kalmanGain * wrapTo180(diff);
        estimatedCovariance = predictedCovariance - kalmanGain * covPrediction * kalmanGain';
    end
    estimatedTracks = estimatedState(:,2:numSteps+1);
end

function [A,W] = getTransitionMatrices(scanTime)
%this is state transistion matrix
A = diag(ones(4,1));
A(1,3) = scanTime;
A(2,4) = scanTime;

%W is noise gain
W = zeros(4,2);
W(1,1) = 0.5*scanTime^2;
W(2,2) = 0.5*scanTime^2;
W(3,1) = scanTime;
W(4,2) = scanTime;
end


function [points, weights] = getSigmaPoints(mean, covariance)
[lengthState, ~] =  size(mean);
numSigmaPoints = 2*lengthState;

%determine sigma point weights
weights = 1/numSigmaPoints * ones(numSigmaPoints,1);

%calculate sigma points
sqrtmCovariance = sqrtm(covariance); %matrix square root
points = zeros(lengthState,2*lengthState);
points(:,1:lengthState) = repmat(mean,1,lengthState) + sqrt(lengthState)*sqrtmCovariance;
points(:,lengthState+1:2*lengthState) = repmat(mean,1,lengthState) - sqrt(lengthState)*sqrtmCovariance;
end

% this function calculates the mean angle from a set of weighted points
function [angle] = angleFromPoints(points, weights)
numPoints = size(points,2);

signs = sign(points);
if(sum(signs*signs(1)) ~= numPoints)
    totalWeightPositive = sum(weights(signs>=0));
    positiveMean = sum(points(signs>=0).*weights(signs>=0)/totalWeightPositive);
    totalWeightNegative = sum(weights(signs<0));
    negativeMean = sum(points(signs<0).*weights(signs<0)/totalWeightNegative);
    
    if((positiveMean - negativeMean) >= 180)
        negativeMean = negativeMean+360;
    end
    angle = positiveMean*totalWeightPositive+negativeMean*totalWeightNegative;
else
    angle = sum(points.*weights);
end

angle = wrapTo180(angle);

end