function [ ws,ms,Ps ] = get_spawnState( ptrack,ess,model )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Adjusted so that a "fair" comparison with the CPHD filter can be made
%   - process noise added to parent cov
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

w_p = ptrack.w;
m_p = ptrack.m;
P_p = ptrack.P;
ng_p = length(w_p);

spawn = model.spawn;
w_t = spawn.w{ess};
d_t = spawn.posDev{ess};
v_t = spawn.velUnitVec{ess};
Q_t_pos = spawn.Qpos{ess};
Q_t_vel = spawn.Qvel{ess};
ng_t = length(w_t);

numComps = ng_p*ng_t;
ws = zeros(numComps,1);
ms = zeros(model.x_dim,numComps);
Ps = zeros(model.x_dim,model.x_dim,numComps);

F = model.A;
count = 0;
for pgidx = 1:ng_p % loop over parent GM comps.
    zeroFlag = 0;
    mPrior = m_p(:,pgidx);
    Pprior = P_p(:,:,pgidx);
    vMagPrior = norm(mPrior([2,4]));
    vMagPre = vMagPrior*spawn.velMult;
    
    vPrior = mPrior([2,4]);
    bear = get_velBearing( vPrior(1), vPrior(2) );
    if isnan(bear) || isinf(bear)
        zeroFlag = 1;
    else
        R = [cos(bear),-sin(bear);sin(bear),cos(bear)];
    end
    for tgidx = 1:ng_t % loop over spawn GM comps.
        count = count + 1;
        
        if zeroFlag
            ws(count) = 1;
            ms(:,count) = [1e8,1e8,1e8,1e8]';
            Ps(:,:,count) = 100*eye(model.x_dim);
            continue;
        end
        
        dpos = R*d_t(:,tgidx);
        velDir = R*v_t(:,tgidx);
        ws(count) = w_p(pgidx)*w_t(tgidx);
        dvec = [dpos(1),velDir(1)*vMagPre-1*mPrior(2),dpos(2),velDir(2)*vMagPre-1*mPrior(4)]';
        ms(:,count) = F*mPrior + dvec;
        Qpos = R*Q_t_pos(:,:,tgidx)*R';
        Qvel = R*Q_t_vel(:,:,tgidx)*R';
        
        p = [Qpos(:,1),zeros(2,1),Qpos(:,2),zeros(2,1)];
        p = [p(1,:);zeros(1,4);p(2,:);zeros(1,4)];
        
        v = [zeros(2,1),Qvel(:,1),zeros(2,1),Qvel(:,2)];
        v = [zeros(1,4);v(1,:);zeros(1,4);v(2,:)];
        
        Ps(:,:,count) = p+v + F*Pprior*F';

    end % loop over spawn GM comps.
end % loop over parent GM comps.

