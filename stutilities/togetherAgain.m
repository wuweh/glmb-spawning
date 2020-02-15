function [ logWTVtemp, tracksTemp ] = togetherAgain( incoming, model )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Let's do the prediction and update, ... together, ... again. With Gaussian Mixtures (GM).
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hypsUpdate = incoming.hyps_update;
tracksUpdate = incoming.tracks_update;
tracksTemp = incoming.tracks_temp;
npredtracks = incoming.npredtracks;
logWTVupdate = incoming.log_wtv_update;
logCardUpdate = incoming.log_cdn_update;
logWTVtemp = cell(size(logWTVupdate));
hbescell = incoming.hbescell;
rankCostsArray = incoming.rankCostsArray;
assnmtArray = incoming.assnmtArray;
comptposArray = incoming.comptposArray;
clutterIntensity = model.clutterIntensity;
lambda_c = model.lambda_c;
Z_k = incoming.Z;
numMeas = size(Z_k,2); 
spawn = model.spawn;
xdim = model.x_dim;

P_S = model.P_S;    Q_S = 1-P_S;
P_D = model.P_D;    Q_D = 1-P_D;
P_T = model.P_T;    Q_T = 1-P_T;
Nspawn = length(model.spawn.w);
N_max = model.N_max;
for n= 0:N_max %loop over updated cardinality
    nidx= n+1;
    numcmp= length(hypsUpdate{nidx});
    if n==0
        numcmp=1; %trick to force loop entry for 0 cardinality prediction and  update using same code
    end
    
    for cidx= 1:numcmp %loop over all components
        hbes= hbescell{nidx}(cidx); %number of h-best to generate (use to allocate proportionally)
        if hbes ~= 0
            %cost matrix
            nbirthtracks= length(model.bar_q);
            nspawntracks = Nspawn;
            nexisttracks= n;
            ntotaltracks= nbirthtracks+nexisttracks*(1+nspawntracks);
            rankcosts = rankCostsArray{nidx}{cidx};
            A = assnmtArray{nidx}{cidx};
            assnmt = A;
            assnmt=assnmt-ntotaltracks; assnmt(assnmt<=0)= assnmt(assnmt<=0)-1; %set "not exist" states to negative assignment
            C = logCardUpdate(nidx)+sum(logWTVupdate{nidx}(cidx))+(-lambda_c) + numMeas*log(clutterIntensity);
            for hidx=1:min(hbes,length(rankcosts)) % loop over assignment array rows
                nupdatetracks= sum(assnmt(hidx,:)>0);
                if nupdatetracks == 0
                    continue;
                end
                if nupdatetracks <= N_max % if <= N_max
                    nuidx = nupdatetracks+1;
                    comptpos = comptposArray{nidx}{cidx}{hidx};
                    w_hyp = 1; % Initialize hypothesis weight
                    
                    % Birth First : All we need here are w_hyp values, no need to re-compute track info since birth is independent
                    % of survival and spawning
                    for tidx = 1:nbirthtracks % loop over birth tracks
                        asstmp = assnmt(hidx,tidx);
                        if asstmp > ntotaltracks
                            zidx = asstmp-ntotaltracks;
                            tempidx = npredtracks*zidx + tidx;
                            w_hyp = w_hyp * model.bar_q(tidx) * (tracksTemp{tempidx}.qz)/clutterIntensity; % P_D already included in .qz
                        elseif asstmp > 0
                            w_hyp = w_hyp * model.bar_q(tidx) * Q_D;
                        else
                            w_hyp = w_hyp * (1-model.bar_q(tidx));
                        end
                    end % loop over birth tracks
                    
                    % Surv and spawn
                    for eidx = 1:nexisttracks % loop over surv tracks
                        tidx = tidx + 1;
                        asstmp = assnmt(hidx,tidx:tidx+nspawntracks);
                        updtabidx = hypsUpdate{nidx}{cidx}{eidx};
                        pretabidx = ((Nspawn+1)*updtabidx-Nspawn) + nbirthtracks;
                        if asstmp(1) > ntotaltracks && isempty(find(asstmp(2:end)>ntotaltracks,1))
                            % If the parent track is detected, but all its spawn tracks are missed, there is no need to recompute 
                            % track values. We can just use the values already in tracksTemp, but, we still need w_hyp values 
                            % just in case.
                            zidx = asstmp(1)-ntotaltracks;
                            tempidx = npredtracks*zidx(end) + pretabidx;
                            w_hyp = w_hyp * P_S * (tracksTemp{tempidx}.qz)/clutterIntensity; % P_D already included in .qz
                            for ess = 1:nspawntracks
                                if asstmp(ess+1) > 0
                                    w_hyp = w_hyp * P_T * Q_D;
                                else
                                    w_hyp = w_hyp * Q_T;
                                end
                            end
                            tidx = tidx + nspawntracks;
                            continue
                        elseif asstmp(1) <= ntotaltracks && ismember(length(find(asstmp(2:end)>ntotaltracks,1)),[0,1])
                            % If the parent is not detected and at most only one spawn track is detected, still no need to 
                            % recompute track values. But again, we still need w_hyp values just in case.
                            if asstmp(1) > 0
                                w_hyp = w_hyp * P_S * Q_D;
                            else
                                w_hyp = w_hyp * Q_S;
                            end
                            for ess = 1:nspawntracks
                                pretabidx = pretabidx + 1;
                                if asstmp(ess+1) > ntotaltracks
                                    zidx = asstmp(ess+1)-ntotaltracks;
                                    tempidx = npredtracks*zidx(end) + pretabidx;
                                    w_hyp = w_hyp * P_T * (tracksTemp{tempidx}.qz)/clutterIntensity; % P_D already included in .qz
                                elseif asstmp(ess+1) > 0
                                    w_hyp = w_hyp * P_T * Q_D;
                                else
                                    w_hyp = w_hyp * Q_T;
                                end
                            end
                            tidx = tidx + nspawntracks;
                            continue
                        end
                        
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        % If we make it to this point in the code, then at a minimum there are two tracks (e.g., parent & spawn, 
                        % or spawn & spawn) that need to be measurement updated together. Not sure if spawn tracks need to be
                        % predicted \textit{AND} updated together if the parent is missed, but doing so requires less code, I
                        % think.
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        
                        
                        % We need to parse the tracks according to whether they are measurement updated or not. If they are not,
                        % there's no need to process them, except to get the values necessary to compute w_hyp. If they are, then
                        % we group them and their associated parameters together.
                        
                        tracks = cell(1+Nspawn,1);
                        parentDetectFlag = 0;
                        zcount = 0;
                        tempidx = [];
                        
                        % Grab all observations and form a [xdim*Nz x 1] array, where Nz is the number of associated observations
                        zidx = asstmp(asstmp>ntotaltracks)-ntotaltracks;
                        Z = Z_k(:,zidx); Z = Z(:);
                        
                        % Parent track
                        tracks{1} = tracksUpdate{updtabidx}; % this is the one parent track we are currently concerned with
                        % Gather survival deviation vectors and process noise covariances
                        numGMcomp = length(tracks{1}.w);
                        tracks{1}.d = zeros(xdim,numGMcomp);
                        tracks{1}.Q = zeros(xdim,xdim,numGMcomp);
                        for gidx = 1:length(tracks{1}.w)
                            tracks{1}.d(:,gidx) = model.d_s;
                            tracks{1}.Q(:,:,gidx) = model.Q_s;
                        end
                        
                        if asstmp(1) > ntotaltracks
                            w_hyp = w_hyp*P_S*P_D/clutterIntensity; % likelihood comes later
                            zcount = zcount + 1;
                            tempidx = [tempidx; npredtracks*zidx(zcount) + pretabidx]; %#ok<AGROW>
                            parentDetectFlag = 1;
                        elseif asstmp(1) > 0
                            w_hyp = w_hyp*P_S*Q_D;
                            % No need to update tracksTemp or hypsTemp
                        else
                            w_hyp = w_hyp*Q_S;
                        end
                        JT = zeros(1,nspawntracks);
                        for ess = 1:nspawntracks % Loop over spawn tracks
                            pretabidx = pretabidx + 1;
                            if asstmp(ess+1) > ntotaltracks
                                tracks{ess+1}.idx = ess;
                                tracks{ess+1}.w = spawn.w{ess};
                                JT(ess) = length(spawn.w{ess});
                                tracks{ess+1}.d = spawn.posDev{ess};
                                tracks{ess+1}.v = spawn.velUnitVec{ess};
                                tracks{ess+1}.Qpos = spawn.Qpos{ess};
                                tracks{ess+1}.Qvel = spawn.Qvel{ess};
                                w_hyp = w_hyp*P_T*P_D/clutterIntensity; % likelihood comes later
                                zcount = zcount + 1;
                                tempidx = [tempidx; npredtracks*zidx(zcount) + pretabidx]; %#ok<AGROW>
                            elseif asstmp(ess+1) > 0
                                w_hyp = w_hyp*P_T*Q_D;
                                % No need to update tracksTemp
                            else
                                w_hyp = w_hyp*Q_T;
                            end
                        end % Loop over spawn tracks
                        JT(JT==0) = [];
                        tracks = tracks(~cellfun('isempty',tracks));
                        
                        [ tracks_out,QzForHyp ] = GMjointPreUpd( tracks,parentDetectFlag,Z,JT,model );
                        w_hyp = w_hyp*QzForHyp;
                        if length(tracks_out) ~= length(tempidx)
                            keyboard; %error('Number of GM tracks and number of tracks_temp indices are not equal.');
                        end
                        tidx = tidx + nspawntracks;
                    end % loop over surv tracks
                    logWTVtemp{nuidx}(comptpos)= C + log(w_hyp); % C is computed at top of assignment row loop
                    if isinf(logWTVtemp{nuidx}(comptpos)) || isnan(logWTVtemp{nuidx}(comptpos))
                        keyboard;
                    end
                end % if <= N_max
            end % loop over assignment array rows
        end % if hbes~=0
    end % loop over all components
end % loop over updated cardinality

end % of the .m function you're in




