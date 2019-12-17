%maybe improve this so that the estimated tracks' size of the scatter point
%is related to the estimated covariance at that step to get an idea of 
%how the estimate improves over time. Would need to save covariance over
%time in the filter code though.
function plotOverTime(trueTracks,observations,estimatedTracks)
    [~,numSteps]=size(trueTracks);
    for step=1:numSteps
        currentObs = cell2mat(observations(step));
        [~,numObs] = size(currentObs);
        p = plot(estimatedTracks(1,1:step),estimatedTracks(2,1:step));
        hold on
        axis([-100 100 -100 100])
        q = plot(trueTracks(1,1:step),trueTracks(2,1:step));
        estColor = get(p,'color');
        trueColor = get(q,'color');
        scatter(estimatedTracks(1,step),estimatedTracks(2,step),30,estColor);
        scatter(trueTracks(1,step),trueTracks(2,step),30,trueColor);
        for obs=1:numObs
            scatter(currentObs(1,obs),currentObs(2,obs),30,'r');
        end
        pause(0.1)
        hold off
    end
end