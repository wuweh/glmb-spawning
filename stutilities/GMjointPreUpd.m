function [ tracks_out,QzForHyp ] = GMjointPreUpd( tracks,parentDetectFlag,Z,JT,model )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                                                      
%
%
% Input:    tracks = {[parent track], [spawn track]^{(i)} : i = 1,...,N_{\beta,upd}} where N_{\beta,upd} is the number of tracks 
%                     spawned from [parent_track] that are to be measurement updated
% 
%           parentDetectFlag - 1 if parent track detected, 0 otherwise 
%
%           Z - measurement vector, already organized w.r.t. associated tracks
%
%           JT = [JT^{(1)},...,JT^{(N_{\beta,upd})}] where each JT^{(j)} is the number of GM components of the jth spawn track
%
%
%
% Assumptions/Notes:
%       [1] All input spawn tracks are measurement updated
%       [2] The parent track may or may not be updated, but is required to predict the spawn tracks
%
%
% Programmer: Daniel 'Stu' Bryant
%       Date: 19-May-2017
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xdim = model.x_dim;
zdim = model.z_dim;

idxStart = @(n) xdim*n-(xdim-1);
idxStop  = @(n) xdim*n;

numTracks = length(tracks);
if parentDetectFlag; numUpdTracks = numTracks; else numUpdTracks = numTracks - 1; end %#ok<SEPEX>
JS = length(tracks{1}.w);

idxTable = makeGMindexTable(JS,JT);
numComps = size(idxTable,2);

%% Prediction
w = ones(numComps,1);
m = zeros(xdim*numTracks,numComps);
P = zeros(xdim*numTracks,xdim*numTracks,numComps);
F = repmat(model.A,numTracks,1);
for cidx = 1:numComps % loop over GM components (i.e., columns of idxTable)
    zeroFlag = 0;
    mPrior = tracks{1}.m(:,idxTable(1,cidx));
    vPrior = model.J*mPrior;
    bear = get_velBearing( vPrior(1), vPrior(2) );
    if isnan(bear) || isinf(bear)
        zeroFlag = 1;
    else
        R = [cos(bear),-sin(bear);sin(bear),cos(bear)];
    end
    Pprior = tracks{1}.P(:,:,idxTable(1,cidx));
    d = zeros(xdim*numTracks,1);
    Q = zeros(xdim*numTracks,xdim*numTracks);
    vMagPre = norm(vPrior)*model.spawn.velMult;
    for tidx = 1:numTracks % loop over rows of idxTable : tidx = 1 is the parent, all others are its spawn
        span = idxStart(tidx):idxStop(tidx);
        tabidx = idxTable(tidx,cidx);
        w(cidx) = w(cidx)*tracks{tidx}.w(tabidx);
        if tidx == 1
            d(span,1) = tracks{tidx}.d(:,(idxTable(tidx,cidx)));
            Q(span,span) = tracks{tidx}.Q(:,:,tabidx);
        else
            gidx = idxTable(tidx,cidx);
            if zeroFlag
                d(span,1) = [1e8,1e8,1e8,1e8]';
                Q(span,span) = 100*eye(model.x_dim);
            else
                dpos = R*tracks{tidx}.d(:,gidx);
                velDir = R*tracks{tidx}.v(:,gidx);
                dvec = [dpos(1),velDir(1)*vMagPre-1*mPrior(2),dpos(2),velDir(2)*vMagPre-1*mPrior(4)]';
                Qpos = R*tracks{tidx}.Qpos(:,:,tabidx)*R';
                Qvel = R*tracks{tidx}.Qvel(:,:,tabidx)*R';
                
                p = [Qpos(:,1),zeros(2,1),Qpos(:,2),zeros(2,1)];
                p = [p(1,:);zeros(1,4);p(2,:);zeros(1,4)];
                
                v = [zeros(2,1),Qvel(:,1),zeros(2,1),Qvel(:,2)];
                v = [zeros(1,4);v(1,:);zeros(1,4);v(2,:)];
                
                d(span,1) = dvec;
                Q(span,span) = p+v;
            end
        end
        
    end % loop over rows of idxTable
    
    m(:,cidx) = F*mPrior + d;
    P(:,:,cidx) = Q + F*Pprior*F';
end % loop over GM components

%% Update
if parentDetectFlag
    mBar = m;
    Pbar = P;
else
    mBar = zeros(xdim*numUpdTracks,numComps);
    Pbar = zeros(xdim*numUpdTracks,xdim*numUpdTracks,numComps);
    for cidx = 1:numComps % loop over GM components (i.e., columns of idxTable)
        mBar(:,cidx) = m(xdim+1:end,cidx);
        Pbar(:,:,cidx) = P(xdim+1:end,xdim+1:end,cidx);
    end % loop over GM components (i.e., columns of idxTable)
end

R = kron(eye(numUpdTracks),model.R); % block diagonal obs. error matrix
H = kron(eye(numUpdTracks),model.C_posn); % block diagonal state to obs. matrix

wTilde = zeros(size(w));
mTilde = zeros(size(mBar));
Ptilde = zeros(size(Pbar));
for cidx = 1:numComps
    mu = H*mBar(:,cidx);
    S  = R+H*Pbar(:,:,cidx)*H'; Vs= chol(S); det_S= prod(diag(Vs))^2; inv_sqrt_S= inv(Vs); iS= inv_sqrt_S*inv_sqrt_S';
    K  = Pbar(:,:,cidx)*H'*iS;
    logqz = -0.5*zdim*log(2*pi) - 0.5*log(det_S) - 0.5*(Z-mu)'*iS*(Z-mu);
    wTilde(cidx) = w(cidx)*exp(logqz);
    mTilde(:,cidx) = mBar(:,cidx) + K*(Z-mu);
    Ptilde(:,:,cidx) = (eye(size(Pbar(:,:,cidx))) - K*H)*Pbar(:,:,cidx);
end

QzForHyp = sum(wTilde); % weight that contributes to current hypothesis weight
if isnan(QzForHyp); error('Whoa whoa whoa! What''s this NaN doing here?'); end
wTilde = wTilde./sum(wTilde);

%% Break each Gaussian in the mixture down according to which track it belongs to
tracks_out = cell(numUpdTracks,1);
for tidx = 1:numUpdTracks
    span = idxStart(tidx):idxStop(tidx);
    tracks_out{tidx}.w = wTilde;
    tracks_out{tidx}.m = mTilde(span,:);
    tracks_out{tidx}.P = Ptilde(span,span,:);
    tracks_out{tidx}.qz = QzForHyp;
end


%% Functions(s)

    function [ tab ] = makeGMindexTable(parent_GMidx,spawn_GMidx) %#ok<INUSL>
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % This function creates an index table used to form "stacked" GM components. Stacked means that, for a given component,
        % its mean and covariance contain elements for multiple tracks. E.g., if a parent and spawn track are predicted together,
        % one GM component mean will have the parent mean "stacked" on top of the spawn mean. Covariances are not stacked, per se,
        % but are formed in a diagonal fashion, however, note that they are not considered block diagonal. 
        %
        % Example: Consider a case where we have a parent track and 2 spawn tracks such that
        %   parent track has:    2 GM components,
        %   spawn track (1) has: 3 GM components,
        %   spawn track (2) has: 2 GM components.
        %
        %   Then, the following index table is formed. Note that 2*3*2 yields 12 new GM components, therefore the number of 
        %   columns = number of new GM components. The number of rows = the number of tracks being predicted together.
        % 
        %         1     1     1     1     1     1     2     2     2     2     2     2
        %         1     1     2     2     3     3     1     1     2     2     3     3
        %         1     2     1     2     1     2     1     2     1     2     1     2
        %
        %   The top row pertains to the incoming parent track GM, the second row pertains to the first incoming spawn track GM,
        %   and so on. So, if we look at the 4th column, it tells us that to form the 4th new GM component, we need
        %       * the parent track's 1st GM comp.
        %       * the first spawn tracks's 2nd GM comp.
        %       * the second spawn track's 2nd GM comp.
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        str = 'tab = combvec([1:parent_GMidx]';
        for i = 1:length(spawn_GMidx)
            str = [str,',[1:spawn_GMidx(',num2str(i),')]'];
        end
        str = [str,');'];   %#ok<*AGROW>
        eval(str)
        tab0 = tab; %#ok<NODEF>
        tab = unique(tab','rows')'; 
        if size(tab0)~=size(tab)
            error('Well, you''ve gone and done it Stu. This function doesn''t work the way you hoped it would. If you are not Stu, please contact him and tell him he sucks.');
        end
    end



end % of the .m function you're in