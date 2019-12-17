% Sam Miller, 2019 (reviewed by F. Meyer)

function [stateEstimates,existences] = bernoulliFilter(observations,parameters)
numSteps = parameters.numSteps;
numDevices = parameters.numDevices;
scanTime = parameters.scanTime;
ps = parameters.survivalProb;
pb = parameters.birthProb;
dVar = parameters.sigmaDrivingNoise^2;
rangeVar = parameters.sigmaMeasurementNoiseRange^2;
velocityVar = parameters.sigmaMeasurementNoiseVelocity^2;
angleVar = parameters.sigmaMeasurementNoiseBearing^2;
measurementRange = parameters.measurementRange;
velocityRange = parameters.velocityRange;
meanClutter = parameters.meanClutter;
pd = parameters.detectionProbability;
detectionThreshold = parameters.detectionThreshold;
maxComponents = parameters.maxComponents;
sensorPosition = parameters.sensorPosition;

lambdaCmono = meanClutter /  (measurementRange * 2*velocityRange * 360);
lambdaCbi = meanClutter / (measurementRange * 360);
measurementCovMono = diag([rangeVar;angleVar;velocityVar]);
measurementCovBi = diag([rangeVar;angleVar]);

[stateTransition,noiseMatrix] = getTransitionMatrices(scanTime);
drivingCov = noiseMatrix*diag([dVar;dVar])*noiseMatrix';

%initialize relevant variables
birthMeans = parameters.birthMeans;
birthCovariances = parameters.birthCovariances;
birthWeights = parameters.birthWeights;
numBirthComponents = size(birthMeans,2);

spatialMeans = zeros(4,0);
spatialCovariances = zeros(4,4,0);
spatialWeights = zeros(0,1);
numComponents = 0;

stateEstimates = nan(4,numSteps);
existences = zeros(1,numSteps);
currentExistence = 0;

for step = 1:numSteps
    step
    
    % ---prediction---
    % predict existence probability
    predictedExistence = pb * (1 - currentExistence) + ps * currentExistence;
    
    % predict mean and covariance of spatial components
    for n = 1:numComponents
        spatialMeans(:,n) = stateTransition * spatialMeans(:,n);
        spatialCovariances(:,:,n) = stateTransition * spatialCovariances(:,:,n) * stateTransition' + drivingCov;
    end
    
    for device = 1:numDevices
        numMeasurements = 2;
        if(device == 1)
            numMeasurements = 3;
        end
        % fill in birth mean/covariance
        spatialMeans = cat(2,spatialMeans,birthMeans);
        spatialCovariances = cat(3,spatialCovariances,birthCovariances);
        
        % ---pre-update---
        currentObs = observations{device + numDevices*(step-1)};
        numObs = size(currentObs,2);

        % get predicted measurements and evaluate likelihood for all component-observation combinations
        measurementMeans = zeros(numMeasurements,numComponents + numBirthComponents);
        measurementCovariances = zeros(numMeasurements,numMeasurements,numComponents + numBirthComponents);
        crossCovariances= zeros(4,numMeasurements,numComponents + numBirthComponents);
        measurementLikelihood = zeros(numComponents + numBirthComponents,numObs);
        for n = 1:(numComponents + numBirthComponents)

            [sigmas,sigmaWeights] = getSigmaPoints(spatialMeans(:,n),spatialCovariances(:,:,n));
            numPoints = size(sigmas,2);

            %propagate sigma points through observation model
            pointsMeasurements = zeros(numMeasurements,numPoints);
            if(device == 1)
                pointsMeasurements(1,:) = sqrt((sigmas(1,:) - sensorPosition(1,1)).^2+(sigmas(2,:) - sensorPosition(1,2)).^2);
                pointsMeasurements(2,:) =  atan2d(sigmas(1,:) - sensorPosition(1,1),sigmas(2,:) - sensorPosition(1,2));
                pointsMeasurements(3,:) = 0.5*(((sigmas(1,:) - sensorPosition(1,1)).^2+(sigmas(2,:) - sensorPosition(1,2)).^2).^(-0.5)).*(2*(sigmas(1,:) - sensorPosition(1,1)).*sigmas(3,:)+2*(sigmas(2,:) - sensorPosition(1,2)).*sigmas(4,:));
            else
                pointsMeasurements(1,:) = sqrt((sensorPosition(device,1)-sigmas(1,:)).^2 + (sensorPosition(device,2)-sigmas(2,:)).^2) + sqrt((sensorPosition(1,1)-sigmas(1,:)).^2 + (sensorPosition(1,2)-sigmas(2,:)).^2);
                pointsMeasurements(2,:) = atan2d(sigmas(1,:) - sensorPosition(device,1),sigmas(2,:) - sensorPosition(device,2));
            end

            %calculate statistics of prediction measurements
            for k = 1:numPoints
                measurementMeans(1,n) = measurementMeans(1,n) + sigmaWeights(k) * pointsMeasurements(1,k);
                if(device == 1)
                    measurementMeans(3,n) = measurementMeans(3,n) + sigmaWeights(k) * pointsMeasurements(3,k);
                end
            end
            measurementMeans(2,n) = angleFromPoints(pointsMeasurements(2,:)', sigmaWeights);

            diffMeasurements = zeros(numMeasurements,numPoints);
            for k = 1:numPoints
                diffMeasurements(:,k) = (pointsMeasurements(:,k) - measurementMeans(:,n));
                diffMeasurements(2,k) = wrapTo180(diffMeasurements(2,k));
                measurementCovariances(:,:,n) = measurementCovariances(:,:,n) + sigmaWeights(k) * (diffMeasurements(:,k) * diffMeasurements(:,k)');
            end
            if(device == 1)
                measurementCovariances(:,:,n) = measurementCovariances(:,:,n) + measurementCovMono;
            else
                measurementCovariances(:,:,n) = measurementCovariances(:,:,n) + measurementCovBi;
            end

            %calculate cross covariance matrixes
            for k = 1:numPoints
                crossCovariances(:,:,n) = crossCovariances(:,:,n) + sigmaWeights(k) * (sigmas(:,k) - spatialMeans(:,n)) * diffMeasurements(:,k)';
            end        

            for m = 1:numObs
                measurementLikelihood(n,m) = gaussian(currentObs(1:numMeasurements,m),measurementMeans(:,n),measurementCovariances(:,:,n));
            end
        end

        % calculate deltaK
        deltaK = 0;
        for n = 1:numComponents
            for m = 1:numObs
                if(device == 1)
                    deltaK = deltaK + (spatialWeights(n) * measurementLikelihood(n,m))/lambdaCmono;
                else
                    deltaK = deltaK + (spatialWeights(n) * measurementLikelihood(n,m))/lambdaCbi;
                end
            end
        end
        deltaK = pd * (1 - deltaK);

        % update birth/spatial weights
        birthComp = (pb * (1 - currentExistence))/predictedExistence;
        survivalComp = (ps * currentExistence)/predictedExistence;
        spatialWeights = cat(1,survivalComp*spatialWeights,birthComp*birthWeights);
        numComponents = numComponents + numBirthComponents;


        % ---update---
        % preallocate variables for update step
        spatialMeans = repmat(spatialMeans,[1,1,numObs+1]);
        spatialCovariances = repmat(spatialCovariances,[1,1,1,numObs+1]);
        spatialWeights = repmat(spatialWeights,[1,numObs+1]);

        % update weights for missed detection case (component means and variances remain unchanged)
        spatialWeights(:,1) = ( (1 - pd)/(1 - deltaK) ) * spatialWeights(:,1);

        % update existence
        currentExistence = ((1 - deltaK)/(1 - deltaK * predictedExistence)) * predictedExistence;

        % this is the detection component of the updated spatial pdf
        updateConstant = pd/(1 - deltaK);
        for n = 1:numComponents
            for m = 1:numObs
                if(device == 1)
                    spatialWeights(n,m+1) = updateConstant * ( (spatialWeights(n,m+1) * measurementLikelihood(n,m)) / lambdaCmono );
                else
                    spatialWeights(n,m+1) = updateConstant * ( (spatialWeights(n,m+1) * measurementLikelihood(n,m)) / lambdaCbi );
                end
                kalmanGain = crossCovariances(:,:,n)/measurementCovariances(:,:,n);

                innovation = (currentObs(1:numMeasurements,m) - measurementMeans(:,n));
                innovation(2) = wrapTo180(innovation(2));            
                spatialMeans(:,n,m+1) = spatialMeans(:,n,m+1) + kalmanGain * innovation;
                spatialCovariances(:,:,n,m+1) = spatialCovariances(:,:,n,m+1) - kalmanGain * measurementCovariances(:,:,n) * kalmanGain';
            end
        end
        % update total number of components in the updated spatial pdf
        numComponents = numComponents + numObs * numComponents;

        % reshape components vectors
        spatialWeights = reshape(spatialWeights,[numComponents, 1]);
        spatialMeans = reshape(spatialMeans,[4,numComponents]);
        spatialCovariances = reshape(spatialCovariances,[4,4,numComponents]);
        
        % ---pruning---
        % keep only the maxComponents components with the heighest weight
        [~,indexes] = sort(spatialWeights,1,'descend');
        indexes = indexes(1:maxComponents);
        spatialWeights = spatialWeights(indexes);
        spatialMeans = spatialMeans(:,indexes);
        spatialCovariances = spatialCovariances(:,:,indexes);
        numComponents = maxComponents;

        % renormalize Weights
        spatialWeights = spatialWeights/sum(spatialWeights,1);
    end

    % ---detection and estimation---
    existences(step) = currentExistence;
    if(existences(step) > detectionThreshold)
        stateEstimates(:,step) = zeros(4,1);
        for n=1:numComponents
            stateEstimates(:,step) = stateEstimates(:,step) + spatialWeights(n) * spatialMeans(:,n);
        end
    end
end
end

function [A,W] = getTransitionMatrices(scanTime)
A = diag([1;1;1;1]);
A(1,3) = scanTime;
A(2,4) = scanTime;

%W is noise gain
W = zeros(4,2);
W(1,1) = 0.5*scanTime^2;
W(2,2) = 0.5*scanTime^2;
W(3,1) = scanTime;
W(4,2) = scanTime;
end

% Important to note this gaussian wraps around for angle
function likelihood = gaussian(x,mean,covariance)
[dim,~] = size(covariance);
detInverse = (det(covariance))^-0.5;
coeff = 1/((2*pi)^(dim/2)) * detInverse;
diff = x - mean;
diff(2) = wrapTo180(diff(2));
likelihood = coeff * exp(-0.5 * diff' * covariance^-1 * diff);
end

% Produce 8 sigma points with the same mean and covariance as the parameters.
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