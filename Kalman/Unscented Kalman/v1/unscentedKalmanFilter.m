function estimatedTracks = unscentedKalmanFilter(observations,parameters)
    %load parameters
    priorMean = parameters.priorMean;
    dim = parameters.dimensionality;
    numPoints = 2 * numel(priorMean);
    numSteps = parameters.numSteps;
    numMeasurements = parameters.numMeasurements;
    dt = parameters.scanTime;
    sensorPosition = parameters.sensorPosition;
    rangeVar = parameters.sigmaMeasurementNoiseRange^2;
    angleVar = parameters.sigmaMeasurementNoiseBearing^2;
    drivingVar = parameters.sigmaDrivingNoise^2;
    priorCovariance = parameters.priorCovariance;
    alpha = parameters.alpha;
    kappa = parameters.kappa;
    beta = parameters.beta;
    
    %declare variables
    stateSize = length(priorMean);
    measurementCovariance = diag([rangeVar;angleVar]); %we're not measuring velocity
    [stateTransition,noiseGain] = getTransitionMatrices(dt); %4x4, 4x2 respectively
    drivingCovariance = noiseGain*diag([drivingVar;drivingVar])*noiseGain'; %4x4 matrix
    estimatedState = zeros(stateSize,numSteps+1);
    estimatedState(:,1) = priorMean;
    estimatedCov = priorCovariance;
    
    for k=1:numSteps
        %predict (normal KF predict equations)
        predictedMean = stateTransition * estimatedState(:,k);
        predictedCovariance = stateTransition * estimatedCov * stateTransition' + drivingCovariance;
        
        %generate weights
        Wm=[0.5/stateSize+zeros(1,2*stateSize)];
        Wc=Wm;

        %generate sigma points
        sigmaPoints = getSigmaPoints(predictedMean,predictedCovariance);
        
        %transform sigma points through process model (not necessary for
        %linear model)
        %transformedState = stateTransition * sigmaPoints; %4x5 matrix

        %propagate sigma points through observation model
        predictionState = zeros(numMeasurements,numPoints); %2x5 matrix
        predictionState(1,:) = sqrt((sigmaPoints(1,:) - sensorPosition(1)).^2+(sigmaPoints(2,:) - sensorPosition(2)).^2);
        predictionState(2,:) =  atan2d(sigmaPoints(1,:) - sensorPosition(1),sigmaPoints(2,:) - sensorPosition(2));
        
        %we don't have states 3 and 4 because we're not measuring velocity

        %calculate statistics of transformed state and prediction state
        meanPrediction = zeros(numMeasurements,1);
        for i=1:numPoints
            meanPrediction(1) = meanPrediction(1) + Wm(i) * predictionState(1,i);
        end
        meanPrediction(2) = anglesFromPoints(predictionState(2,:)',Wm);

        covPrediction = zeros(numMeasurements,numMeasurements);
        for i=1:numPoints
            diffMeasurements = predictionState(:,i) - meanPrediction;
            covPrediction = covPrediction + Wc(i) * (diffMeasurements * diffMeasurements');
        end
        covPrediction = covPrediction + measurementCovariance;

        %calculate cross covariance matrix
        covCross = zeros(stateSize,numMeasurements);
        for i=1:numPoints
            covCross = covCross + Wc(i) * (sigmaPoints(:,i)- predictedMean) * (predictionState(:,i)-meanPrediction)';
        end

        %update with Kalman Filter equations
        kalmanGain = covCross / covPrediction;
        estimatedState(:,k+1) = predictedMean + kalmanGain * wrapTo180(observations(:,k) - meanPrediction);
        temp = kalmanGain * covPrediction * kalmanGain';
        estimatedCov = predictedCovariance - temp;
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

function sigmaPoints=getSigmaPoints(mean, covariance)
    lengthState = length(mean);
    sigmaPoints = zeros(lengthState,2*lengthState);
    sigmaPoints(:,1:lengthState)=repmat(mean,1,lengthState)+sqrt(lengthState)*sqrtm(covariance);
    sigmaPoints(:,lengthState+1:lengthState*2)=repmat(mean,1,lengthState)-sqrt(lengthState)*sqrtm(covariance);
end

function angles=anglesFromPoints(points,weights)
    numPoints = size(points,2);

    signs = sign(points);
    if(sum(signs*signs(1)) ~= numPoints)
        totalWeightPositive = sum(weights(signs>=0));
        positiveMean = sum(points(signs>=0).*weights(signs>=0)'/totalWeightPositive);
        totalWeightNegative = sum(weights(signs<0));
        negativeMean = sum(points(signs<0).*weights(signs<0)'/totalWeightNegative);

        if((positiveMean - negativeMean) >= 180)
            negativeMean = negativeMean+360;
        end
        angle = positiveMean*totalWeightPositive+negativeMean*totalWeightNegative;
    else
        angle = sum(points.*weights);
    end

    angles = wrapTo180(angle);
end