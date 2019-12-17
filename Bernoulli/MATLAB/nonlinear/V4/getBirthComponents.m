function [parameters] = getBirthComponents(parameters,rangeStd,velStd)
    %make this a grid of births such that 1 standard deviation of
    %each birth mean covers the entire area.
    limits = parameters.limits;
    [numMeasurements,~] = size(parameters.startState);

    birthMeans = zeros(numMeasurements,0);
    birthCovariances = zeros(numMeasurements,numMeasurements,0);
    for x=(limits(1,1)+rangeStd):2*rangeStd:limits(1,2)
        for y = (limits(2,1)+rangeStd):2*rangeStd:limits(2,2)
            for xVel=(limits(3,1)+velStd):2*velStd:limits(3,2)
                for yVel=(limits(4,1)+velStd):2*velStd:limits(4,2)
                    birthMeans = cat(2,birthMeans,[x;y;xVel;yVel]);
                    birthCovariances = cat(3,birthCovariances,diag([rangeStd^2;rangeStd^2;velStd^2;velStd^2]));
                end
            end
        end
    end
    
    parameters.birthMeans = birthMeans;
    parameters.birthCovariances = birthCovariances;
    
    numComponents = size(birthMeans,2);
    parameters.birthWeights = (1/numComponents)*ones(numComponents,1);
end