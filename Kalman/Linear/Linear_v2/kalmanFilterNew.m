function estimatedTracks = kalmanFilterNew(observations,parameters)
    % This code implements the Kalman filter
    dt=parameters.scanTime;
    priorMean = parameters.priorMean;%xhat
    priorCovariance = parameters.priorCovariance;%P matrix
    %stateTransition = [1 0 dt 0;0 1 0 dt;0 0 1 0;0 0 0 1]; %phi matrix, assume constant acceleration.
    stateTransform = diag([1;1;1;1]); %H matrix
    mVar=parameters.sigmaMeasurementNoise^2;
    dVar=parameters.sigmaDrivingNoise^2;
    measurementCovariance = diag([mVar;mVar;mVar;mVar]);%R matrix
    %drivingCovariance = diag([dVar;dVar;dVar;dVar]);%Q matrix
    numSteps = parameters.numSteps;
    
    [stateTransition,noiseMatrix] = getTransitionMatrices(dt);
    drivingCovariance = noiseMatrix*diag([dVar;dVar])*noiseMatrix'; 
    
    
    estimatedTracks = zeros(size(priorMean,1),numSteps+1);
    
    estimatedTracks(:,1) = priorMean;
    estimatedCovariance = priorCovariance;
    for k=1:numSteps
        %predict
        priorMean = stateTransition*estimatedTracks(:,k);
        priorCovariance = stateTransition*estimatedCovariance*stateTransition'+drivingCovariance;
        predictedMeasurement = stateTransform*priorMean;
        measurementPredictionCovariance = stateTransform * priorCovariance*stateTransform' + measurementCovariance;
        
        %update
        kalmanGain = (priorCovariance * stateTransform') / measurementPredictionCovariance;
        estimatedTracks(:,k+1) = priorMean + kalmanGain * (observations(:,k) - predictedMeasurement);
        estimatedCovariance = priorCovariance - kalmanGain * stateTransform * priorCovariance; 
    end    
    estimatedTracks = estimatedTracks(:,2:numSteps+1); 
    
end


function [A,W] = getTransitionMatrices(scanTime)
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
