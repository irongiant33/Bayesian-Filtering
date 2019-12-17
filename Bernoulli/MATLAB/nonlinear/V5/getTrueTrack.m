function trueTracks = getTrueTrack(parameters,apperanceFromTo)
numSteps = parameters.numSteps;
scanTime = parameters.scanTime;
sigmaDrivingNoise = parameters.sigmaDrivingNoise;
startState = parameters.startState;
[stateTransition,noiseGain] = getTransitionMatrices(scanTime);

trueTracks = nan(4,numSteps+1);
currentState = startState;
for step = 1:numSteps
    
    currentState = stateTransition*currentState + noiseGain*sigmaDrivingNoise*randn(2,1);
    
    if (step >= apperanceFromTo(1) && step <= apperanceFromTo(2))
        %next state is previous state (distance plus scanTime*velocity)
        %plus some process noise drawn from the normal distribution. The noise
        %is affected by the scanTime according to the noiseGain
        trueTracks(:,step+1) = currentState;
    end
    
end
trueTracks = trueTracks(:,2:numSteps+1);

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