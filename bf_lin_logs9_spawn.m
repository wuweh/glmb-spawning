%==== M-Best- code base from CPHD_LIN_F3 


%==== Run parameters

X = scenario.X;          % True multi-target states
L = scenario.L;          % "True" track labels
N_true = scenario.N_true;    % True cardinality


marg_flag= 0; %0/1 on/off for MDGLMB approximation
elim_threshold= 1e-5; GMthresh.elim_threshold = elim_threshold; %for pruning of Gaussians inside tracks
merge_threshold= 4; GMthresh.merge_threshold = merge_threshold; %for merging of Gaussians inside tracks
cap_threshold= 5; GMthresh.cap_threshold =cap_threshold; %for capping of Gaussians inside tracks


Hbes= 1000; %cap number of components

stdvpts= [0.1 0.2 1.1]; %cardinality std dev cut off points for internal requested number of components
Cmin= 0; %min components per cardinality for pruning
chop_threshold = 1e-19; %pruning threshold for component weights

gate_flag= 1;               %0/1 off/on
% P_G= 0.999999999;             %gate probability (for analog and digital gating)  
P_G= 0.9999999;
gamma= chi2inv(P_G,z_dim);  %for analog gating only - inv chi^2 dn gamma value

Z_orig = Z;                 %copy origina l measurements - Z{k} will be replaced with gated measurements
origclutrate= lambda_c; log_origclutrate= log_lambda_c; %original clutter rate
origclutpdf= clutterpdf; log_origclutpdf= log_clutterpdf; %original clutter pdf


%==== Start filtering

hat_N= zeros(K,1);
hat_X= cell(K,1);
hat_P= cell(K,1);
hat_T= cell(K,1);

ospa = zeros(K,1);
eloc = zeros(K,1);
ecar = zeros(K,1);

% ospaLab = ospa;


cdn_update_stack  = zeros(N_max+1,K);
cdn_update_mean   = zeros(K,1);
cdn_update_mode   = zeros(K,1);
cdn_update_var    = zeros(K,1);

emm=0;
for k=1:K
    emm= max(emm,size(Z{k},2));
end
tmp_N_max= max(emm,N_max)+3;

%---precalculate constants used prediction and update calculations
nvector  = [0:N_max]'; %#ok<*NBRAK>
logPSpow = [0:tmp_N_max]'*log(P_S);
logQSpow = [0:tmp_N_max]'*log(1-P_S);
logPDpow = [0:tmp_N_max]'*log(P_D);
logQDpow = [0:tmp_N_max]'*log(1-P_D);

% PULL THESE LINES OUTSIDE OF THIS SCRIPT FOR MC RUNS
% THESE ONLY DEPEND ON N_MAX AND DO NOT CHANGE FOR SIMULATION RUNS
% precalculate values of P(n,j) and C(n,j)
logCcoef = zeros(tmp_N_max+1,tmp_N_max+1);
logPcoef = zeros(tmp_N_max+1,tmp_N_max+1);
logfactorial = zeros(tmp_N_max+1,1);
logfactorial(1) = 0;
for n=1:tmp_N_max
    logfactorial(n+1) = log(n)+logfactorial(n);
end

for ell=0:tmp_N_max
    for j=0:ell
        logPcoef(ell+1,j+1) = logfactorial(ell+1)-logfactorial(ell-j+1);
        logCcoef(ell+1,j+1) = logPcoef(ell+1,j+1)-logfactorial(j+1);
    end
end
% end calculations for P(n,j) and C(n,j)

if ~exist('run_flag','var')
    run_flag = 'disp';
end

%--- i'll need to record these stuff for performance analysis
gaus_size_lg= zeros(K,1);
cpu_time_lg= zeros(K,1);

log_cdn_update = -realmax*ones(N_max+1,1); log_cdn_update(1)= 0;
cdn_update= zeros(N_max+1,1); cdn_update(1)= 1;
log_wtv_update= cell(N_max+1,1);  log_wtv_update{1}= 0;
wtv_update= cell(N_max+1,1);  wtv_update{1}= 1;
tracks_update= cell(0,1); 
hyps_update= cell(N_max+1,1);


tstart = 1;


for k=1:K
    time_start= cputime;
    
    %---stochastic component selection/allocation
    prev_cvar= ([0:N_max].^2*cdn_update(:)) - ([0:N_max]*cdn_update(:))^2;
    prev_stdv= sqrt(prev_cvar);
    Hbesuse= Hbesreq(find(cumsum(stdvpts)>prev_stdv,1)); 
    hbescell= cell(N_max+1,1);
    sampledn= resample(cdn_update,Hbesuse);
    for n=0:N_max
        nidx= n+1;
        hbescell{nidx}= zeros(length(wtv_update{nidx}),1);
        sampledc= resample(wtv_update{nidx},nnz(sampledn==nidx));
        for cidx=1:length(wtv_update{nidx})
            hbescell{nidx}(cidx)= nnz(sampledc==cidx);
        end
    end
    
    %---update
    %init params
    log_cdn_temp = -realmax*ones(N_max+1,1); 
    log_wtv_temp= cell(N_max+1,1); 
    tracks_predict= cell(length(model.bar_q)+length(tracks_update)*(1+model.nSpawn),1);  
    hyps_temp= cell(N_max+1,1); 
    
    %no of measurements
    m = size(Z{k},2);

    %h-best update
        %create birth tracks
        for tabbidx=1:length(model.bar_q)
            tracks_predict{tabbidx}.m = model.bar_x{tabbidx};
            tracks_predict{tabbidx}.P = model.bar_Q{tabbidx};
            tracks_predict{tabbidx}.w = model.lambda_b{tabbidx}(:);
            tracks_predict{tabbidx}.l = [k;tabbidx];
            tracks_predict{tabbidx}.ah = [];
            tracks_predict{tabbidx}.avps = model.bar_q(tabbidx);
            tracks_predict{tabbidx}.avqs = 1-model.bar_q(tabbidx);
            tracks_predict{tabbidx}.avps_tempd = model.bar_q_tempd(tabbidx);
            tracks_predict{tabbidx}.avqs_tempd = 1-model.bar_q_tempd(tabbidx);
            tracks_predict{tabbidx}.velHist = [];
            tracks_predict{tabbidx}.bearing.hist = [];
            tracks_predict{tabbidx}.bearing.est = [];
            tracks_predict{tabbidx}.stateHist = [];
        end
        
        %create predicted surviving and spawn tracks
        tmpPreIdx = [];
        tmpParIdxPre = [];
        pretabidx = length(model.bar_q);
        for tabsidx=1:length(tracks_update) % loop over surviving tracks
            pretabidx = pretabidx + 1;
            
            [wtemp_predict,mtemp_predict,Ptemp_predict]= kalman_predict_sum(1,model.A,model.Q,tracks_update{tabsidx}.w,tracks_update{tabsidx}.m,tracks_update{tabsidx}.P);
            tracks_predict{pretabidx}.velHist = tracks_update{tabsidx}.velHist;
            tracks_predict{pretabidx}.stateHist = tracks_update{tabsidx}.stateHist;
            tracks_predict{pretabidx}.bearing = tracks_update{tabsidx}.bearing;
            tracks_predict{pretabidx}.m = mtemp_predict;
            tracks_predict{pretabidx}.P = Ptemp_predict;
            tracks_predict{pretabidx}.w = wtemp_predict;
            tracks_predict{pretabidx}.l = tracks_update{tabsidx}.l;
            tracks_predict{pretabidx}.ah = tracks_update{tabsidx}.ah;
            tracks_predict{pretabidx}.avps = P_S;
            tracks_predict{pretabidx}.avqs = Q_S;
            tracks_predict{pretabidx}.avps_tempd = P_S_tempd;
            tracks_predict{pretabidx}.avqs_tempd = Q_S_tempd;

            for ess = 1:model.nSpawn % loop over spawn tracks
                pretabidx = pretabidx + 1;
                
                
                [ wStemp_predict,mStemp_predict,PStemp_predict ] = get_spawnState( tracks_update{tabsidx},ess,model );
                
                tracks_predict{pretabidx}.w = wStemp_predict;
                tracks_predict{pretabidx}.m = mStemp_predict;
                tracks_predict{pretabidx}.P = PStemp_predict;
                tracks_predict{pretabidx}.l = [tracks_update{tabsidx}.l;[k;ess]];
                tracks_predict{pretabidx}.ah = tracks_update{tabsidx}.ah;
                tracks_predict{pretabidx}.avps = P_T;
                tracks_predict{pretabidx}.avqs = Q_T;
                tracks_predict{pretabidx}.avps_tempd = P_T_tempd;
                tracks_predict{pretabidx}.avqs_tempd = Q_T_tempd;
                tracks_predict{pretabidx}.velHist = [];
                tracks_predict{pretabidx}.stateHist = [];
                tracks_predict{pretabidx}.bearing.hist = [];
                tracks_predict{pretabidx}.bearing.est  = [];
            end % loop over spawn tracks

        end % loop over surviving tracks
       
        %track level gating
        if gate_flag          
            if m~=0
                %find gated measurements
                %(operates track by track, not by predicted PHD)
                valid_idx= [];
                for tabidx=1:length(tracks_predict)
                    tracks_predict{tabidx}.gate_meas= [];
                    for j=1:length(tracks_predict{tabidx}.w)
                        Sj= model.R + model.C_posn*tracks_predict{tabidx}.P(:,:,j)*model.C_posn';
                        Vs= chol(Sj); det_Sj= prod(diag(Vs))^2; inv_sqrt_Sj= inv(Vs);
                        iSj= inv_sqrt_Sj*inv_sqrt_Sj';
                        nu= Z{k}- model.C_posn*repmat(tracks_predict{tabidx}.m(:,j),[1 m]);
                        dist= sum((inv_sqrt_Sj'*nu).^2);
                        valid_idx= union(valid_idx,find( dist < gamma ));
                        tracks_predict{tabidx}.gate_meas= union(tracks_predict{tabidx}.gate_meas,find( dist < gamma ));
                    end
                end
                
                %delete non gated measurements and reindex gated measurements
                %(optional but eliminates globally unused measurements)
                gate_meas_newidx= zeros(1,size(Z{k},2));
                gate_meas_newidx(valid_idx)= 1:length(valid_idx);
                
                for tabidx=1:length(tracks_predict)
                    tracks_predict{tabidx}.gate_meas= gate_meas_newidx(tracks_predict{tabidx}.gate_meas);    
                end
                
                Z{k} = Z{k}(:,valid_idx); %#ok<SAGROW>
                m = size(Z{k},2);
                
            end
        end
        
        tracks_temp= cell(length(tracks_predict)*(1+size(Z{k},2)),1);
        %create temporary updated update tracks (legacy ones first)
        for tabidx= 1:length(tracks_predict)
            tracks_temp{tabidx}= tracks_predict{tabidx};
            tracks_temp{tabidx}.qz= NaN;
            tracks_temp{tabidx}.ah= [tracks_predict{tabidx}.ah(:); 0];
            tracks_temp{tabidx}.avpd= P_D;
            tracks_temp{tabidx}.avqd= Q_D;
            tracks_temp{tabidx}.avpd_tempd= P_D_tempd;
            tracks_temp{tabidx}.avqd_tempd= Q_D_tempd;
            tracks_temp{tabidx}.used= 0;
            [~,iw] = max(tracks_temp{tabidx}.w);
            tracks_temp{tabidx}.stateHist = [tracks_temp{tabidx}.stateHist,tracks_temp{tabidx}.m(:,iw)];
        end
        
        %create temporary updated update tracks (now measurement updated ones, organized in blocks of predicted tracks, one for each received measurement)
        for emm= 1:m
            for tabidx= 1:length(tracks_predict)
                if ~gate_flag || any(emm == tracks_predict{tabidx}.gate_meas) %if gating is off do all updates automatically, or if meas is validated do single target update
                    stoidx= length(tracks_predict)*emm + tabidx; %index of predicted track i updated with measurement j is (number_predicted+tracks*j + i)
                    [wtemp_update,mtemp_update,Ptemp_update] = kalman_update_sum(Z{k}(:,emm),1,model.C_posn,zeros(z_dim,1),model.R,tracks_predict{tabidx}.w,tracks_predict{tabidx}.m,tracks_predict{tabidx}.P);
                    tracks_temp{stoidx}= tracks_predict{tabidx};
                    tracks_temp{stoidx}.m = mtemp_update;
                    tracks_temp{stoidx}.P = Ptemp_update;
                    tracks_temp{stoidx}.qz = P_D*sum(wtemp_update); % bayes evidence for construction of cost matrix in data association step
                    tracks_temp{stoidx}.qz_tempd = P_D_tempd*sum(wtemp_update); % tempered bayes evidence for construction of cost matrix in data association step
                    tracks_temp{stoidx}.w = wtemp_update/sum(wtemp_update);
                    [~,iw] = max(tracks_temp{stoidx}.w);
                    tracks_temp{stoidx}.stateHist = [tracks_temp{stoidx}.stateHist,tracks_temp{stoidx}.m(:,iw)];
                    tracks_temp{stoidx}.ah= [tracks_predict{tabidx}.ah(:); emm];
                    tracks_temp{stoidx}.avpd_tempd= P_D_tempd;
                    tracks_temp{stoidx}.avqd_tempd= Q_D_tempd;
                    tracks_temp{stoidx}.used= 0;
                    
                    tracks_temp{stoidx}.w= tracks_temp{stoidx}.w/sum(tracks_temp{stoidx}.w);
                    [tracks_temp{stoidx}.w,tracks_temp{stoidx}.m,tracks_temp{stoidx}.P]= gaus_prune(tracks_temp{stoidx}.w,tracks_temp{stoidx}.m,tracks_temp{stoidx}.P,elim_threshold);
                    [tracks_temp{stoidx}.w,tracks_temp{stoidx}.m,tracks_temp{stoidx}.P]= gaus_merge(tracks_temp{stoidx}.w,tracks_temp{stoidx}.m,tracks_temp{stoidx}.P,merge_threshold);
                    [tracks_temp{stoidx}.w,tracks_temp{stoidx}.m,tracks_temp{stoidx}.P]= gaus_cap(tracks_temp{stoidx}.w,tracks_temp{stoidx}.m,tracks_temp{stoidx}.P,cap_threshold);
                        
                    [~, widx] = max( tracks_temp{stoidx}.w );
                    vel = model.C_vel*tracks_temp{stoidx}.m(:,widx);
                    tracks_temp{stoidx}.velHist = [tracks_temp{stoidx}.velHist,vel];
                    bear = get_velBearing( vel(1),vel(2) );
                    if ~isnan(bear) && ~isinf(bear)
                        tracks_temp{stoidx}.bearing.hist = [tracks_temp{stoidx}.bearing.hist,bear];
                    end
                    if ~isempty(tracks_temp{stoidx}.bearing.hist)
                        tracks_temp{stoidx}.bearing.est = mean(tracks_temp{stoidx}.bearing.hist);
                    end
                    
                    
                else %meas is not validated, don't bother updating, set evidence to zero
                    stoidx= length(tracks_predict)*emm + tabidx; %index of predicted track i updated with measurement j is (number_predicted+tracks*j + i)
                    tracks_temp{stoidx}.qz= 0; %bayes evidence is identically zero
                    tracks_temp{stoidx}.used= 0;
                    tracks_temp{stoidx}.gate_meas= [];
                    tracks_temp{stoidx}.velHist = [];
                    tracks_temp{stoidx}.stateHist = [];
                    tracks_temp{stoidx}.bearing.hist = [];
                    tracks_temp{stoidx}.bearing.est = [];
                end
            end
        end

       
        %update hypotheses
        assnmtArray = cell(size(hyps_update));
        rankCostsArray = cell(size(assnmtArray));
        comptposArray = cell(size(assnmtArray));
        for n= 0:N_max %loop over updated cardinality
                nidx= n+1;
                numcmp= length(hyps_update{nidx});
                if n==0 
                   numcmp=1; %trick to force loop entry for 0 cardinality prediction and  update using same code 
                end
                    
                for cidx= 1:numcmp %loop over all components
                    hbes= hbescell{nidx}(cidx); %number of h-best to generate (use to allocate proportionally)
                    if hbes ~= 0
                        %cost matrix
                        nbirthtracks= length(model.bar_q);
                        nspawntracks = model.nSpawn;
                        nexisttracks= n;
                        ntotaltracks= nbirthtracks+nexisttracks*(1+nspawntracks);
                        %--- [Begin] Make association table for given update component
                        % row 1: tracks_update table index, 0 for birth, parent and child have same index
                        % row 2: prediction table index
                        assocTable = zeros(2,ntotaltracks);
                        %--- [End] Make association table for given update component
                        PSvec= zeros(ntotaltracks,1);
                        PDvec= zeros(ntotaltracks,1);
                        QSvec= zeros(ntotaltracks,1);
                        QDvec= zeros(ntotaltracks,1);
                        costm= zeros(ntotaltracks,m);
                        %tempered cost matrix
                        PSvec_tempd= zeros(ntotaltracks,1);
                        PDvec_tempd= zeros(ntotaltracks,1);
                        QSvec_tempd= zeros(ntotaltracks,1);
                        QDvec_tempd= zeros(ntotaltracks,1);
                        costm_tempd= zeros(ntotaltracks,m);
                        %calculate values for birth tracks
                        for bidx= 1:nbirthtracks
                            assocTable(:,bidx) = [0;bidx];
                            PSvec(bidx)= tracks_temp{bidx}.avps;
                            PDvec(bidx)= tracks_temp{bidx}.avpd;
                            QSvec(bidx)= tracks_temp{bidx}.avqs;
                            QDvec(bidx)= tracks_temp{bidx}.avqd;
                            PSvec_tempd(bidx)= tracks_temp{bidx}.avps_tempd;
                            PDvec_tempd(bidx)= tracks_temp{bidx}.avpd_tempd;
                            QSvec_tempd(bidx)= tracks_temp{bidx}.avqs_tempd;
                            QDvec_tempd(bidx)= tracks_temp{bidx}.avqd_tempd;
                            for emm= 1:m
                                linidx= length(tracks_predict)*emm+bidx;
                                if tracks_temp{linidx}.qz    %i.e. must be non-zero                             
                                    costm(bidx,emm)= PSvec(bidx)/QSvec(bidx)*tracks_temp{linidx}.qz/(lambda_c*clutterpdf*QDvec(bidx));
                                    costm_tempd(bidx,emm)= PSvec_tempd(bidx)/QSvec_tempd(bidx)*tracks_temp{linidx}.qz_tempd/(lambda_c*clutterpdf*QDvec_tempd(bidx));
                                end
                            end 
                        end
                        %calculate values for existing and spawn tracks
                        cost_row = nbirthtracks;
                        for tee= 1:nexisttracks % loop over existing tracks
                            cost_row = cost_row + 1;
                            updtabidx = hyps_update{nidx}{cidx}{tee};
                            pretabidx = ( (nspawntracks+1)*updtabidx - nspawntracks ) + nbirthtracks; 
                            assocTable(:,cost_row) = [updtabidx,pretabidx]';
                            PSvec(cost_row)= tracks_temp{pretabidx}.avps;
                            PDvec(cost_row)= tracks_temp{pretabidx}.avpd;
                            QSvec(cost_row)= tracks_temp{pretabidx}.avqs;
                            QDvec(cost_row)= tracks_temp{pretabidx}.avqd;
                            PSvec_tempd(cost_row)= tracks_temp{pretabidx}.avps_tempd;
                            PDvec_tempd(cost_row)= tracks_temp{pretabidx}.avpd_tempd;
                            QSvec_tempd(cost_row)= tracks_temp{pretabidx}.avqs_tempd;
                            QDvec_tempd(cost_row)= tracks_temp{pretabidx}.avqd_tempd;
                            for emm= 1:m
                                linidx= length(tracks_predict)*emm+pretabidx;
                                if tracks_temp{linidx}.qz    %i.e. must be non-zero  
                                    costm(cost_row,emm)= PSvec(cost_row)/QSvec(cost_row)*tracks_temp{linidx}.qz/(lambda_c*clutterpdf*QDvec(cost_row));
                                    costm_tempd(cost_row,emm)= PSvec_tempd(cost_row)/QSvec_tempd(cost_row)*tracks_temp{linidx}.qz_tempd/(lambda_c*clutterpdf*QDvec_tempd(cost_row));
                                end
                            end
                            for ess = 1:nspawntracks % loop over spawn tracks
                                cost_row = cost_row + 1;
                                pretabidx = pretabidx + 1;
                                assocTable(:,cost_row) = [updtabidx,pretabidx]';
                                PSvec(cost_row)= tracks_temp{pretabidx}.avps;
                                PDvec(cost_row)= tracks_temp{pretabidx}.avpd;
                                QSvec(cost_row)= tracks_temp{pretabidx}.avqs;
                                QDvec(cost_row)= tracks_temp{pretabidx}.avqd;
                                PSvec_tempd(cost_row)= tracks_temp{pretabidx}.avps_tempd;
                                PDvec_tempd(cost_row)= tracks_temp{pretabidx}.avpd_tempd;
                                QSvec_tempd(cost_row)= tracks_temp{pretabidx}.avqs_tempd;
                                QDvec_tempd(cost_row)= tracks_temp{pretabidx}.avqd_tempd;
                                for emm= 1:m
                                    linidx= length(tracks_predict)*emm+pretabidx;
                                    if tracks_temp{linidx}.qz    %i.e. must be non-zero
                                        costm(cost_row,emm)= PSvec(cost_row)/QSvec(cost_row)*tracks_temp{linidx}.qz/(lambda_c*clutterpdf*QDvec(cost_row));
                                        costm_tempd(cost_row,emm)= PSvec_tempd(cost_row)/QSvec_tempd(cost_row)*tracks_temp{linidx}.qz_tempd/(lambda_c*clutterpdf*QDvec_tempd(cost_row));
                                    end
                                end
                            end % loop over spawn tracks
                        end % loop over existing tracks
                        costm= [diag(1./QDvec) diag(PSvec./QSvec) costm]; %#ok<AGROW> %[notsurvived survived_but_not_detected survived_detected_and_generated_measurement]
                        costm_tempd= [diag(1./QDvec_tempd) diag(PSvec_tempd./QSvec_tempd) costm_tempd]; %#ok<AGROW> %[notsurvived survived_but_not_detected survived_detected_and_generated_measurement]
                        
                        neglogcostm= -log(costm); %DON'T transpose to leave tracks on rows instead of measurements (track to measurement assignment)
                        neglogcostm_tempd= -log(costm_tempd); %DON'T transpose to leave tracks on rows instead of measurements (track to measurement assignment)
                         
                        switch model.cost_method
                            case 'Gibbs'
                                 %USE EITHER: gibbs sampling trick (one target per measurement, and one measurement per measurement)
                                [assnmt,nlcost]= mbestwrap_updt_gibbsamp_tempered(neglogcostm,neglogcostm_tempd,hbes); rankscosts= exp(-nlcost);
                            case 'Murtys'
                                %OR: optimal assignment trick (one target per measurement, and one measurement per measurement)
                                [assnmt,nlcost]= mbestwrap_updt_joint(neglogcostm,hbes); rankscosts= exp(-nlcost);
                        end
                        
                        assnmtArray{nidx}{cidx} = assnmt;
                        rankCostsArray{nidx}{cidx} = rankscosts;
                        
                        assnmt=assnmt-ntotaltracks; assnmt(assnmt<=0)= assnmt(assnmt<=0)-1; %set not born/not survived states to negative assignment
                        %meas update/clutter update for tracks
                        for hidx=1:min(hbes,length(rankscosts))    
                            nupdatetracks= sum(assnmt(hidx,:)>0);
                            if nupdatetracks <= N_max
                                nuidx= nupdatetracks+1;
                                comptpos= length(hyps_temp{nuidx})+1;
                                comptposArray{nidx}{cidx}{hidx} = comptpos;
                                trackcounter=1;
                                for tidx=1:ntotaltracks
                                    asstmp= assnmt(hidx,tidx);
                                    if asstmp > 0
                                        
                                        %offset of current track is position in predicted track table
                                        newoffset= assocTable(2,tidx);
                                        
                                        %index of corresponding updated track
                                        if asstmp > ntotaltracks    %measurement assignment
                                            linidx= length(tracks_predict)*(asstmp-ntotaltracks)+newoffset;
                                        elseif asstmp > 0           %missed detection
                                            linidx= newoffset;
                                        end
                                        hyps_temp{nuidx}{comptpos}{trackcounter}= linidx; trackcounter=trackcounter+1;
                                        tracks_temp{linidx}.used = 1;
                                    end
                                end
                            end

                            if nupdatetracks == 0
                                log_wtv_temp{nuidx}= logsumexp([log_wtv_temp{nuidx},log_cdn_update(nidx)+sum(log_wtv_update{nidx}(cidx))+(-lambda_c)+sum(log(QSvec))+sum(log(QDvec))+m*(log_lambda_c+log_clutterpdf)+(-nlcost(hidx))]);
                            elseif nupdatetracks <= N_max
                                log_wtv_temp{nuidx}= cat(1,log_wtv_temp{nuidx},log_cdn_update(nidx)+sum(log_wtv_update{nidx}(cidx))+(-lambda_c)+sum(log(QSvec))+sum(log(QDvec))+m*(log_lambda_c+log_clutterpdf)+(-nlcost(hidx)));
                            end
                        end
                    end
                end               
        end

      
        % [Begin: Re-Predict/Update] --------------------------------------------------
        clutterIntensity = lambda_c*clutterpdf;
        if do_joint
            
            incoming.hyps_update = hyps_update;
            incoming.tracks_update = tracks_update;
            incoming.tracks_temp = tracks_temp;
            incoming.npredtracks = length(tracks_predict);
            incoming.log_wtv_update = log_wtv_update;
            incoming.log_cdn_update = log_cdn_update;
            incoming.hbescell = hbescell;
            incoming.rankCostsArray = rankCostsArray;
            incoming.assnmtArray = assnmtArray;
            incoming.comptposArray = comptposArray;
            incoming.Z = Z{k};
            
            [ log_wtv_temp, tracks_temp ] = togetherAgain( incoming, model );
        end
%         wRE = exp(log_wtv_RE);
        
        
        % [End: Re-Predict/Update] ----------------------------------------------------
    
    %calc cardinality distribution and normalize component weights
    for n= 0:N_max
        nidx= n+1;
        log_wtv_temp{nidx}= log_wtv_temp{nidx}+eps(0);
        if isempty(log_wtv_temp{nidx}), log_cdn_temp(nidx)=-inf; else log_cdn_temp(nidx)= logsumexp(log_wtv_temp{nidx}); end
        log_wtv_temp{nidx}= log_wtv_temp{nidx}-logsumexp(log_wtv_temp{nidx});
    end
    log_wtv_temp{1}=0;
    wtv_temp= cell(N_max+1,1); for n=0:N_max, nidx= n+1; wtv_temp{nidx}= exp(log_wtv_temp{nidx}); end
    
    log_cdn_temp= log_cdn_temp - logsumexp(log_cdn_temp); cdn_temp= exp(log_cdn_temp); 
    cdn_temp= cdn_temp/sum(cdn_temp);
    
    cdn_update= cdn_temp; log_cdn_update= log_cdn_temp;
 
    tracks_update= tracks_temp; hyps_update= hyps_temp; log_wtv_update= log_wtv_temp; wtv_update= wtv_temp; %#ok<NASGU>
    
    %--- D-GLMB marginalization step (optional) [new marginalized tracks are appended to the track table, redundant tracks are deleted in the next step]
    if marg_flag
        tracks_mdglmb= tracks_update; hyps_mdglmb= cell(N_max+1,1); log_wtv_mdglmb= cell(N_max+1,1); wtv_mdglmb= cell(N_max+1,1); %#ok<UNRCH>
        
        hyps_mdglmb{1}= hyps_update{1};
        log_wtv_mdglmb{1}= 0;
        wtv_mdglmb{1}= 1;
        
        trackcount= length(tracks_mdglmb)+1;
        for n=1:N_max
            nidx= n+1;
            hypcount= 1;
            markermerged= zeros(length(log_wtv_update{nidx}),1);
            tmp_wtv_mdglmb= [];
            tmp_wtv_update= exp(log_wtv_update{nidx});
            for cidx=1:length(hyps_update{nidx})
                
                if markermerged(cidx) == 0
                    %init vector of component indices for which will be merged
                    tmp_merge_cidxs= [cidx];
                    
                    %grab label set for current component
                    vectorlabels= zeros(2,n);
                    for tidx=1:n
                        vectorlabels(:,tidx)= tracks_update{hyps_update{nidx}{cidx}{tidx}}.l;
                    end
                    
                    %scan for other label sets matching current component and build vector of component indices for merging
                    for vidx=cidx+1:length(hyps_update{nidx})
                        %if the component hasn't already been merged check it
                        if markermerged(vidx) == 0
                            %grab label set for next component
                            checkvectorlabels= zeros(2,n);
                            for tidx=1:n
                                checkvectorlabels(:,tidx)= tracks_update{hyps_update{nidx}{vidx}{tidx}}.l;
                            end
                            
                            %test for match with label sets and add index for merging later
                            if isempty(setdiff(vectorlabels',checkvectorlabels','rows')) %use set diff to compare unordered sets of labels for match
                                tmp_merge_cidxs= cat(1,tmp_merge_cidxs,vidx);
                                markermerged(vidx)= 1;
                            end
                            
                        end
                    end
                    
                    %if the hypothesis/component is not merged with another, simply copy the existing track out (saves a lot of memory compared to creating new redundant tracks everytime)
                    if length(tmp_merge_cidxs)==1
                        hyps_mdglmb{nidx}{hypcount}= hyps_update{nidx}{cidx};
                        tmp_wtv_mdglmb(hypcount)= tmp_wtv_update(cidx);
                        for tidx=1:n
                            tracks_mdglmb{hyps_update{nidx}{cidx}{tidx}}.ah= [];
                        end
                        %else do the merging indicated by the vector of indices
                    else
                        %create merged component
                        for tidx=1:n
                            hyps_mdglmb{nidx}{hypcount}{tidx}= trackcount;
                            tracks_mdglmb{trackcount}.m= [];
                            tracks_mdglmb{trackcount}.P= [];
                            tracks_mdglmb{trackcount}.w= [];
                            tracks_mdglmb{trackcount}.l= tracks_update{hyps_update{nidx}{cidx}{tidx}}.l;
                            tracks_mdglmb{trackcount}.ah= [];
                            tracks_mdglmb{trackcount}.used= 1;
                            tracks_mdglmb{trackcount}.newidx= [];
                            tracks_mdglmb{trackcount}.qz= [];
                            trackcount= trackcount+1;
                        end
                        %loop over hypotheses/components and create merged tracks
                        for midx= 1:length(tmp_merge_cidxs)
                            vidx= tmp_merge_cidxs(midx);
                            for tidx=1:n
                                mergetrackidx= trackcount-(n+1)+tidx;
                                tracks_mdglmb{mergetrackidx}.m= cat(2,tracks_mdglmb{mergetrackidx}.m,tracks_update{hyps_update{nidx}{vidx}{tidx}}.m);
                                tracks_mdglmb{mergetrackidx}.P= cat(3,tracks_mdglmb{mergetrackidx}.P,tracks_update{hyps_update{nidx}{vidx}{tidx}}.P);
                                tracks_mdglmb{mergetrackidx}.w= cat(1,tracks_mdglmb{mergetrackidx}.w,tmp_wtv_update(vidx)*tracks_update{hyps_update{nidx}{vidx}{tidx}}.w);
                            end
                        end
                        tmp_wtv_mdglmb(hypcount)= sum(tmp_wtv_update(tmp_merge_cidxs));
                        
                        %normalize track weights for newly merged tracks then perform GM mixture reduction by pruning, merging and truncation on indiviudal tracks
                        for tidx=1:n
                            mergetrackidx= trackcount-(n+1)+tidx;
                            tracks_mdglmb{mergetrackidx}.w= tracks_mdglmb{mergetrackidx}.w/sum(tracks_mdglmb{mergetrackidx}.w);
                            [tracks_mdglmb{mergetrackidx}.w,tracks_mdglmb{mergetrackidx}.m,tracks_mdglmb{mergetrackidx}.P]= gaus_prune(tracks_mdglmb{mergetrackidx}.w,tracks_mdglmb{mergetrackidx}.m,tracks_mdglmb{mergetrackidx}.P,elim_threshold);
                            [tracks_mdglmb{mergetrackidx}.w,tracks_mdglmb{mergetrackidx}.m,tracks_mdglmb{mergetrackidx}.P]= gaus_merge(tracks_mdglmb{mergetrackidx}.w,tracks_mdglmb{mergetrackidx}.m,tracks_mdglmb{mergetrackidx}.P,merge_threshold);
                            [tracks_mdglmb{mergetrackidx}.w,tracks_mdglmb{mergetrackidx}.m,tracks_mdglmb{mergetrackidx}.P]= gaus_cap(tracks_mdglmb{mergetrackidx}.w,tracks_mdglmb{mergetrackidx}.m,tracks_mdglmb{mergetrackidx}.P,cap_threshold);
                        end
                    end
                    %increment hypothesis counter
                    hypcount= hypcount+1;
                end
            end

            %store merged (normalized) hypothesis/component weights
            if ~isempty(tmp_wtv_mdglmb)
                log_wtv_mdglmb{nidx}= log(tmp_wtv_mdglmb(:)); wtv_mdglmb{nidx}= tmp_wtv_mdglmb(:);
            end

        end
        hyps_update= hyps_mdglmb; tracks_update= tracks_mdglmb; log_wtv_update= log_wtv_mdglmb; wtv_update= wtv_mdglmb; 
    end
    
    
    %--- Prune/merge/cap Gaussian Mixtures of used tracks
    for updidx = 1:length(tracks_update)
        if tracks_update{updidx}.used
            [tracks_update{updidx}.w,tracks_update{updidx}.m,tracks_update{updidx}.P]= gaus_prune(tracks_update{updidx}.w,tracks_update{updidx}.m,tracks_update{updidx}.P,GMthresh.elim_threshold);
            try
                [tracks_update{updidx}.w,tracks_update{updidx}.m,tracks_update{updidx}.P]= gaus_merge(tracks_update{updidx}.w,tracks_update{updidx}.m,tracks_update{updidx}.P,GMthresh.merge_threshold);
            catch
                keyboard;
            end
            [tracks_update{updidx}.w,tracks_update{updidx}.m,tracks_update{updidx}.P]= gaus_cap(tracks_update{updidx}.w,tracks_update{updidx}.m,tracks_update{updidx}.P,GMthresh.cap_threshold);
        end
    end
    
    
    %--- merge duplicate hypotheses from posterior
    %flatten posterior into component table
    hcell= cell(0,0);
    tcell= cell(0,0);
    nvect= zeros(0,0);
    wvect= zeros(0,0);
    hidx=1;
    for n=0:N_max
        nidx=n+1;
        if n==0
            hcell{hidx}= sprintf('%i*',[]);
            tcell{hidx}= {};
            nvect(hidx)= n;
            wvect(hidx)= log_wtv_update{nidx};
            hidx= hidx+1;
        else
            for cidx= 1:length(hyps_update{nidx})
                trackpointerstemp= zeros(n,1);
                for tidx=1:n
                    trackpointerstemp(tidx)= hyps_update{nidx}{cidx}{tidx};
                end
                hcell{hidx}= sprintf('%i*',sort(trackpointerstemp(:)')); %#ok<TRSRT>
                tcell{hidx}= hyps_update{nidx}{cidx};
                nvect(hidx)= n;
                wvect(hidx)= log_wtv_update{nidx}(cidx);
                hidx= hidx+1;
            end
        end
    end
    
    %find unique components and preallocate memory
    [nc,ia,ic]= unique(hcell);
    hyps_temp= cell(N_max+1,1); log_wtv_temp= cell(N_max+1,1); wtv_temp= cell(N_max+1,1);
    for n=0:N_max
        nidx= n+1;
        ncomps= sum(nvect(ia)==n);
        if ncomps ~=0 
            log_wtv_temp{nidx}= NaN*zeros(ncomps,1); wtv_temp{nidx}= NaN*zeros(ncomps,1);
        end
    end
    %write new components
    for hidx= 1:length(nc)
        write_idx_n= length(tcell{ia(hidx)})+1;
        write_idx_c= find(isnan(log_wtv_temp{write_idx_n}),1);
        if ~isempty(tcell{ia(hidx)}), hyps_temp{write_idx_n}{write_idx_c}= tcell{ia(hidx)}; end
        log_wtv_temp{write_idx_n}(write_idx_c,1)= logsumexp(wvect(ic==hidx));
        wtv_temp{write_idx_n}(write_idx_c,1)= exp(log_wtv_temp{write_idx_n}(write_idx_c,1));
    end
    hyps_update= hyps_temp; log_wtv_update= log_wtv_temp; wtv_update= wtv_temp;
    
    %--- compact track table and reindex hypthesis cell
    tracks_temp= cell(0,1); hyps_temp= cell(N_max+1,1);
    trackcount= 0;
    for tabidx= 1:length(tracks_update)
        if tracks_update{tabidx}.used == 1
            trackcount= trackcount+1;
            tracks_update{tabidx}.newidx= trackcount;
            tracks_temp{trackcount,1}= tracks_update{tabidx};
        end
    end
    
    for n=0:N_max
        nidx= n+1;
        for cidx=1:length(hyps_update{nidx})
            for tidx= 1:n
                hyps_temp{nidx}{cidx}{tidx}= tracks_update{hyps_update{nidx}{cidx}{tidx}}.newidx;
            end
        end
    end
      
    for tabidx= 1:length(tracks_temp)
        tracks_temp{tabidx}= rmfield(tracks_temp{tabidx},{'used','newidx','qz'});
    end
    
    tracks_update= tracks_temp; hyps_update= hyps_temp;

    %--- chopping and truncating
 
    %rank hypothesis by weight and store cardinality and birth labels
    comparstack= [];
    weightstack= [];
    nlabelstack= [];
    clabelstack= [];
    for n= 1:N_max
        nidx= n+1;
        numcmp= length(hyps_update{nidx});
        comparstack= cat(1, comparstack, cdn_update(nidx)*wtv_update{nidx});
        weightstack= cat(1, weightstack, wtv_update{nidx});
        nlabelstack= cat(1, nlabelstack, nidx*ones(numcmp,1));
        clabelstack= cat(1, clabelstack, [1:numcmp]');
    end
    totcompraw= length(clabelstack);
    
    [idxkeep]= find(comparstack > chop_threshold);
    comparstack= comparstack(idxkeep);
    weightstack= weightstack(idxkeep); 
    nlabelstack= nlabelstack(idxkeep); 
    clabelstack= clabelstack(idxkeep);
    
    [~,idxsort]= sort(-comparstack); totcomp= length(comparstack); idxsort= idxsort(1:min(totcomp,Hbes));
    comparstack= comparstack(idxsort); weightstack= weightstack(idxsort); nlabelstack= nlabelstack(idxsort); clabelstack= clabelstack(idxsort);
    
    wtv_temp= cell(N_max+1,1); hyps_temp= cell(N_max+1,1);
    
    for u=1:min(totcomp,Hbes) %copy best hypotheses into temp variables
        nidx= nlabelstack(u);
        cidx= clabelstack(u);
        
        wtv_temp{nidx}= cat(1,wtv_temp{nidx},wtv_update{nidx}(cidx)); 
        hyps_temp{nidx}{end+1}= hyps_update{nidx}{cidx};
    end 
    wtv_temp{1}= wtv_update{1}; hyps_temp{1}= hyps_update{1}; %copy zero cardinality
    
    for n= 1:N_max %enforce minimum number per cardinality
        nidx= n+1;
        if length(wtv_temp{nidx}) < Cmin && length(wtv_update{nidx}) ~= 0 %#ok<ISMT>
            [~,idxcomp]= sort(-wtv_update{nidx});
            idxcomp= idxcomp(1:min(Cmin,length(wtv_update{nidx})));
            for cidx=1:length(idxcomp)
                wtv_temp{nidx}(cidx)= wtv_update{nidx}(idxcomp(cidx));
                hyps_temp{nidx}{cidx}= hyps_update{nidx}{idxcomp(cidx)};
            end
        end
        wtv_temp{nidx}= wtv_temp{nidx}/sum(wtv_temp{nidx});
    end
    wtv_update= wtv_temp; hyps_update= hyps_temp;
    for n=0:N_max
        nidx=n+1; 
        log_wtv_update{nidx}= log(wtv_update{nidx}); 
        cdn_update(nidx)= cdn_update(nidx)*sum(wtv_update{nidx}); %resets the cdn to zero if all components in this cardinality were truncated 
        log_cdn_update(nidx)= log(cdn_update(nidx));
    end
    
    
    %--- state extraction

    %use MAP estimate for cardinality dn and pick highest weights
    [~,mode] = max(cdn_update);
    hat_N_MAP(k) = mode-1; %#ok<SAGROW>
    hat_N(k) = hat_N_MAP(k);
    [~,idx1]= max(wtv_update{mode});
    hat_X{k}= []; 
    hat_P{k}= [];
    hat_T{k}= [];
    tmpidx = [];
    
%     chosenIdx = [];
    
    
    if exist('cardFig','var')
        figure(cardFig);
        plot(k,hat_N(k),'k.','MarkerSize',8);
    end
    for n=1:hat_N(k)
        [~,idx2]= max(tracks_update{hyps_update{mode}{idx1}{n}}.w);
        
        updIDX = hyps_update{mode}{idx1}{n};
        
%         chosenIdx = [chosenIdx;updIDX];
        est_m= tracks_update{hyps_update{mode}{idx1}{n}}.m(:,idx2);
        
        est_P= tracks_update{hyps_update{mode}{idx1}{n}}.P(:,:,idx2);
        est_l= tracks_update{hyps_update{mode}{idx1}{n}}.l;
        
        est_l = [est_l;zeros(2*k-length(est_l),1)]; %#ok<AGROW> %----------------------------------- Change for spawn: pad with zeros to account for varying label lengths
        
        hat_X{k} = [hat_X{k}, est_m];
        hat_P{k} = [hat_P{k}, est_P];
        hat_T{k} = [hat_T{k}, est_l];
    end
    hat_T{k}(sum(hat_T{k},2)==0,:) = []; %---------------------------------------------- Change for spawn: remove excess padding

    if exist('plotOverTime','var') && plotOverTime == 1 && exist('noPlots','var') && ~noPlots
        figure(XYfig);
        zk = Z_tracks{k};
        for zidx = 1:size(zk,2)
            hz = plot(zk(1,zidx),zk(2,zidx),'kx','MarkerSize',6);
            
        end
        if exist('cardFig','var')
            figure(cardFig);
            plot(k,hat_N(k),'k.','MarkerSize',12);
        end
        
        if exist('cardFig','var') || exist('XYfig','var')
            pause(0.01);
        end
    end
    
    %compute diagnoistics for cardinality distribution at update
    [~,mode] = max(cdn_update);
    cdn_update_mean(k) = sum(nvector .* cdn_update);
    cdn_update_mode(k) = mode-1;
    cdn_update_var(k)  = sum(nvector.^2 .* cdn_update) - cdn_update_mean(k)^2;

    
    %--- stats
    cpu_time_lg(k)= cputime-time_start;  
    
    %--- store
    cdn_update_stack(:,k) = cdn_update; 
    
    Xtrue = [X{k}(1,:);X{k}(3,:)];
    if isempty(hat_X{k})
        Xest = [];
    else
        Xest = [hat_X{k}(1,:);hat_X{k}(3,:)];
    end
    [ospa(k),eloc(k),ecar(k)] = new_dist(Xtrue,Xest,100,1);
    
    if exist('plotOverTime','var') && plotOverTime == 1 && exist('noPlots','var') && ~noPlots
        figure(ospaFig);
        plot(1:k,ospa(1:k),'m--');
    end
    
    labels =  fliplr(hat_T{k});
 
    fprintf('time= %4d, EN(k)  = %4.3f, VARN(k)= %4.5f, #comp=      %4d, #trak= %4d, #comp merg= %4d\n\n',k,cdn_update_mean(k),cdn_update_var(k),totcompraw,length(tracks_update),min(totcomp,Hbes));
    fprintf('N(k)= %4d, MODN(k)= %5d,    OSPA= %4.5f\n\n',N_true(k),cdn_update_mode(k),ospa(k));


    fprintf('label estimates =\n');
    [nlab,mlab] = size(labels);
    for nidx = 1:nlab
        for midx = 1:mlab-1
            fprintf('  %2d',labels(nidx,midx));
        end
        fprintf('  %2d\n',labels(nidx,end));
    end
    fprintf('\n\n');

end
    
%put the measurements back once we're done
Z_gate = Z; Z = Z_orig;
lambda_c= origclutrate; log_lambda_c= log_origclutrate; clutterpdf= origclutpdf; log_clutterpdf= log_origclutpdf;

