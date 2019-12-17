function rmse = getError(trueTracks,estimatedTracks)
    rmse=zeros(1,length(trueTracks));
    for k=1:length(trueTracks)
        rmse(k)=sqrt((trueTracks(1,k)-estimatedTracks(1,k))^2+(trueTracks(2,k)-estimatedTracks(2,k))^2);
    end
end