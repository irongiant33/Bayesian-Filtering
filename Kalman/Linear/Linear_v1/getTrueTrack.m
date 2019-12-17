function trueTracks = getTrueTrack(parameters)
numSteps = parameters.numSteps;
scanTime = parameters.scanTime;
sigmaDrivingNoise = parameters.sigmaDrivingNoise;
startState = parameters.startState;
[A,W] = getTransitionMatrices(scanTime);

trueTracks = zeros(4,numSteps+1);
trueTracks(:,1) = startState;
for step = 1:numSteps
    %next state is previous state plus some process noise drawn from the normal distribution.
    trueTracks(:,step+1) = A*trueTracks(:,step) + W*sigmaDrivingNoise*randn(2,1);
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