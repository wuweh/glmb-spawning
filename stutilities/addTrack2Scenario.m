function scenario = addTrack2Scenario( track, scenario, model )

scenario.trackCount = scenario.trackCount + 1;
scenario.tBirth = [scenario.tBirth,track.t0];

Xsplit = cell(model.K,1);

numEpochs = track.tf - track.t0; 
X = zeros(model.x_dim, numEpochs + 1);
X(:,1) = track.X0;
for k = 2:numEpochs + 1
    X(:,k) = model.A * X(:,k-1);
end

scenario.truthPlot.X = [scenario.truthPlot.X;X];
scenario.truthPlot.type = [scenario.truthPlot.type,track.type];

count = 0;
for k = track.t0:track.tf
    count = count + 1;
    Xsplit{k} = X(:,count);
    scenario.X{k} = [scenario.X{k},X(:,count)];
    scenario.L{k}  = concat_labels( scenario.L{k},track.L );
    scenario.track_list{k} = [scenario.track_list{k},scenario.trackCount];
    scenario.N_true(k) = scenario.N_true(k) + 1;
end

scenario.Xsplit{scenario.trackCount} = Xsplit;

end

