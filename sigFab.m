declare_problem;

scenario.X = cell(K,1);          % True multi-target states
scenario.Xsplit = {};
scenario.L = cell(K,1);          % "True" track labels
scenario.N_true = zeros(K,1);    % True cardinality
scenario.track_list = cell(K,1); % Track number list at time k
scenario.truthPlot.X = {};       % True multi-target states parsed for plotting
scenario.truthPlot.type = {};    % Target type, e.g., 'birth' or 'spawn', for plotting
scenario.track_list = cell(K,1); % Index location of each track in scenario.X{k} for each time k
scenario.tBirth = [];
scenario.trackCount = 0;

vmags = [];

%% Generate scenario birth and spawn tracks
bSpeed = 10;

oneParentOneSpawn = 0;
connect = {};
connCount = 0;

% Birth 1
track{1}.t0 = 1;
track{1}.tf = 35;
track{1}.region = 1;
track{1}.speed = bSpeed;
track{1}.bearing = deg2rad(300+15);
pos0 = model.C_posn*model.bar_x{track{1}.region}(:,1);
vel0 = track{1}.speed*[cos(track{1}.bearing);sin(track{1}.bearing)];
track{1}.X0 = [pos0(1),vel0(1),pos0(2),vel0(2)]';
track{1}.L = [track{1}.t0;track{1}.region];
track{1}.type = 'birth';

scenario = addTrack2Scenario( track{1}, scenario, model );

vmags = [vmags; norm(track{1}.X0([2,4]))];

if ~oneParentOneSpawn
    
    % Birth 2
    track{2}.t0 = 2;
    track{2}.tf = 37;
    track{2}.region = 2;
    track{2}.speed = bSpeed;
    track{2}.bearing = deg2rad(180+15);
    pos0 = model.C_posn*model.bar_x{track{2}.region}(:,1);
    vel0 = track{2}.speed*[cos(track{2}.bearing);sin(track{2}.bearing)];
    track{2}.X0 = [pos0(1),vel0(1),pos0(2),vel0(2)]';
    track{2}.L = [track{2}.t0;track{2}.region];
    track{2}.type = 'birth';
    
    scenario = addTrack2Scenario( track{2}, scenario, model );
    
    vmags = [vmags; norm(track{2}.X0([2,4]))];
    
    
    % Birth 3
    track{3}.t0 = 3;
    track{3}.tf = 40;
    track{3}.region = 3;
    track{3}.speed = bSpeed;
    track{3}.bearing = deg2rad(60+15);
    pos0 = model.C_posn*model.bar_x{track{3}.region}(:,1);
    vel0 = track{3}.speed*[cos(track{3}.bearing);sin(track{3}.bearing)];
    track{3}.X0 = [pos0(1),vel0(1),pos0(2),vel0(2)]';
    track{3}.L = [track{3}.t0;track{3}.region];
    track{3}.type = 'birth';
    
    scenario = addTrack2Scenario( track{3}, scenario, model );
    
    vmags = [vmags; norm(track{3}.X0([2,4]))];
end

t_intsect = 45;
p_intsect = [0;0];

% Spawn 1
track{4}.t0 = 10;
track{4}.tf = K;
track{4}.parent = track{1};
track{4}.L = [track{4}.parent.L;track{4}.t0;1];
xPar = model.F(track{4}.t0-track{4}.parent.t0)*track{4}.parent.X0;
vPar = xPar([2,4]);
vBear = get_velBearing( vPar(1),vPar(2) );
theta = vBear - pi/2;
pos0 =  spawn.devMag.mean*[cos(theta),sin(theta)]' + model.C_posn*xPar;
[ phi, vSpawn ] = get_posBearing( pos0, p_intsect,t_intsect-track{4}.t0 );
vel0 = vSpawn*[cos(phi);sin(phi)];
track{4}.X0 = [pos0(1),vel0(1),pos0(2),vel0(2)]';

dBear = rad2deg(phi - vBear);
while dBear > 360; dBear = dBear - 360;end
while dBear < 0; dBear = dBear + 360; end
track{4}.dBear = dBear;

track{4}.type = 'spawn';

scenario = addTrack2Scenario( track{4}, scenario, model );

vmags = [vmags; norm(track{4}.X0([2,4]))];

connCount = connCount + 1;
connect{connCount} = [(model.C_posn*xPar)';pos0'];





if ~oneParentOneSpawn
    % Spawn 2
    track{5}.t0 = 11;
    track{5}.tf = K;
    track{5}.parent = track{2};
    track{5}.L = [track{5}.parent.L;track{5}.t0;1];
    xPar = model.F(track{5}.t0-track{5}.parent.t0)*track{5}.parent.X0;
    vPar = xPar([2,4]);
    vBear = get_velBearing( vPar(1),vPar(2) );
    theta = vBear - pi/2;
    pos0 =  spawn.devMag.mean*[cos(theta),sin(theta)]' + model.C_posn*xPar;
    [ phi, vSpawn ] = get_posBearing( pos0, p_intsect,t_intsect-track{5}.t0 );
    vel0 = vSpawn*[cos(phi);sin(phi)];
    track{5}.X0 = [pos0(1),vel0(1),pos0(2),vel0(2)]';
    
    dBear = rad2deg(phi - vBear);
    while dBear > 360; dBear = dBear - 360;end
    while dBear < 0; dBear = dBear + 360; end
    track{5}.dBear = dBear;
    
    track{5}.type = 'spawn';
    
    scenario = addTrack2Scenario( track{5}, scenario, model );
    
    vmags = [vmags; norm(track{5}.X0([2,4]))];
    
    connCount = connCount + 1;
    connect{connCount} = [(model.C_posn*xPar)';pos0'];


    
    % Spawn 3
    track{6}.t0 = 12;
    track{6}.tf = K;
    track{6}.parent = track{3};
    track{6}.L = [track{6}.parent.L;track{6}.t0;1];
    xPar = model.F(track{6}.t0-track{6}.parent.t0)*track{6}.parent.X0;
    vPar = xPar([2,4]);
    vBear = get_velBearing( vPar(1),vPar(2) );
    theta = vBear - pi/2;
    pos0 =  spawn.devMag.mean*[cos(theta),sin(theta)]' + model.C_posn*xPar;
    [ phi, vSpawn ] = get_posBearing( pos0, p_intsect,t_intsect-track{6}.t0 );
    vel0 = vSpawn*[cos(phi);sin(phi)];
    track{6}.X0 = [pos0(1),vel0(1),pos0(2),vel0(2)]';
    
    dBear = rad2deg(phi - vBear);
    while dBear > 360; dBear = dBear - 360;end
    while dBear < 0; dBear = dBear + 360; end
    track{6}.dBear = dBear;
    
    track{6}.type = 'spawn';
    
    scenario = addTrack2Scenario( track{6}, scenario, model );
    
    vmags = [vmags; norm(track{6}.X0([2,4]))];
    
    connCount = connCount + 1;
    connect{connCount} = [(model.C_posn*xPar)';pos0'];

end




if ~oneParentOneSpawn
    
    % Birth 4
    track{7}.t0 = 57;
    track{7}.tf = K;
    track{7}.region = 1;
    track{7}.speed = 10;
    track{7}.bearing = deg2rad(210-15);
    pos0 = model.C_posn*model.bar_x{track{7}.region}(:,1);
    vel0 = track{7}.speed*[cos(track{7}.bearing);sin(track{7}.bearing)];
    track{7}.X0 = [pos0(1),vel0(1),pos0(2),vel0(2)]';
    track{7}.L = [track{7}.t0;track{7}.region];
    track{7}.type = 'birth';
    
    scenario = addTrack2Scenario( track{7}, scenario, model );
    
    vmags = [vmags; norm(track{7}.X0([2,4]))];
    
    
    
    
    % Birth 5
    track{8}.t0 = 59;
    track{8}.tf = K;
    track{8}.region = 2;
    track{8}.speed = 10.6;
    track{8}.bearing = deg2rad(90-15);
    pos0 = model.C_posn*model.bar_x{track{8}.region}(:,1);
    vel0 = track{8}.speed*[cos(track{8}.bearing);sin(track{8}.bearing)];
    track{8}.X0 = [pos0(1),vel0(1),pos0(2),vel0(2)]';
    track{8}.L = [track{8}.t0;track{8}.region];
    track{8}.type = 'birth';
    
    scenario = addTrack2Scenario( track{8}, scenario, model );
    
    vmags = [vmags; norm(track{8}.X0([2,4]))];
    
    
    % Birth 6
    track{9}.t0 = 55;
    track{9}.tf = K;
    track{9}.region = 3;
    track{9}.speed = 9.6;
    track{9}.bearing = deg2rad(330-15);
    pos0 = model.C_posn*model.bar_x{track{9}.region}(:,1);
    vel0 = track{9}.speed*[cos(track{9}.bearing);sin(track{9}.bearing)];
    track{9}.X0 = [pos0(1),vel0(1),pos0(2),vel0(2)]';
    track{9}.L = [track{9}.t0;track{9}.region];
    track{9}.type = 'birth';
    
    scenario = addTrack2Scenario( track{9}, scenario, model );
    
    vmags = [vmags; norm(track{9}.X0([2,4]))];
    
    
    
    % Spawn 4
    track{10}.t0 = 56;
    track{10}.tf = K;
    track{10}.parent = track{4};
    track{10}.L = [track{10}.parent.L;track{10}.t0;1];
    intTrackId = 9;
    t_intsect = 82;
    p_intsect = model.C_posn*model.F(t_intsect-track{intTrackId}.t0)*track{intTrackId}.X0;
    xPar = model.F(track{10}.t0-track{10}.parent.t0)*track{10}.parent.X0;
    vPar = xPar([2,4]);
    vBear = get_velBearing( vPar(1),vPar(2) );
    theta = vBear - pi/2;
    pos0 =  spawn.devMag.mean*[cos(theta),sin(theta)]' + model.C_posn*xPar;
    [ phi, vSpawn ] = get_posBearing( pos0, p_intsect,t_intsect-track{10}.t0 );
    vel0 = vSpawn*[cos(phi);sin(phi)];
    track{10}.X0 = [pos0(1),vel0(1),pos0(2),vel0(2)]';
    
    dBear = rad2deg(phi - vBear);
    while dBear > 360; dBear = dBear - 360;end
    while dBear < 0; dBear = dBear + 360; end
    track{10}.dBear = dBear;
    
    track{10}.type = 'spawn';
    
    scenario = addTrack2Scenario( track{10}, scenario, model );
    
    vmags = [vmags; norm(track{10}.X0([2,4]))];
    
    connCount = connCount + 1;
    connect{connCount} = [(model.C_posn*xPar)';pos0'];
    
    
    % Spawn 5
    track{11}.t0 = 58;
    track{11}.tf = K;
    track{11}.parent = track{5};
    track{11}.L = [track{11}.parent.L;track{11}.t0;1];
    intTrackId = 7;
    t_intsect = 84;
    p_intsect = model.C_posn*model.F(t_intsect-track{intTrackId}.t0)*track{intTrackId}.X0;
    xPar = model.F(track{11}.t0-track{11}.parent.t0)*track{11}.parent.X0;
    vPar = xPar([2,4]);
    vBear = get_velBearing( vPar(1),vPar(2) );
    theta = vBear - pi/2;
    pos0 =  spawn.devMag.mean*[cos(theta),sin(theta)]' + model.C_posn*xPar;
    [ phi, vSpawn ] = get_posBearing( pos0, p_intsect,t_intsect-track{11}.t0 );
    vel0 = vSpawn*[cos(phi);sin(phi)];
    track{11}.X0 = [pos0(1),vel0(1),pos0(2),vel0(2)]';
    
    dBear = rad2deg(phi - vBear);
    while dBear > 360; dBear = dBear - 360;end
    while dBear < 0; dBear = dBear + 360; end
    track{11}.dBear = dBear;
    
    track{11}.type = 'spawn';
    
    scenario = addTrack2Scenario( track{11}, scenario, model );
    
    vmags = [vmags; norm(track{11}.X0([2,4]))];
    
    connCount = connCount + 1;
    connect{connCount} = [(model.C_posn*xPar)';pos0'];
    
    
    
    
    
    % Spawn 6
    track{12}.t0 = 60;
    track{12}.tf = K;
    track{12}.parent = track{6};
    track{12}.L = [track{12}.parent.L;track{12}.t0;1];
    intTrackId = 8;
    t_intsect = 86;
    p_intsect = model.C_posn*model.F(t_intsect-track{intTrackId}.t0)*track{intTrackId}.X0;
    xPar = model.F(track{12}.t0-track{12}.parent.t0)*track{12}.parent.X0;
    vPar = xPar([2,4]);
    vBear = get_velBearing( vPar(1),vPar(2) );
    theta = vBear - pi/2;
    pos0 =  spawn.devMag.mean*[cos(theta),sin(theta)]' + model.C_posn*xPar;
    [ phi, vSpawn ] = get_posBearing( pos0, p_intsect,t_intsect-track{12}.t0 );
    vel0 = vSpawn*[cos(phi);sin(phi)];
    track{12}.X0 = [pos0(1),vel0(1),pos0(2),vel0(2)]';
    
    dBear = rad2deg(phi - vBear);
    while dBear > 360; dBear = dBear - 360;end
    while dBear < 0; dBear = dBear + 360; end
    track{12}.dBear = dBear;
    
    track{12}.type = 'spawn';
    
    scenario = addTrack2Scenario( track{12}, scenario, model );
    
    vmags = [vmags; norm(track{12}.X0([2,4]))];
    
    connCount = connCount + 1;
    connect{connCount} = [(model.C_posn*xPar)';pos0'];
    
end




% % Velocities
% vmags
% 
% for tidx = 1:length(track)
%     if isfield(track{tidx},'dBear')
%         fprintf('Track %d, dBear = %f deg\n',tidx,360-track{tidx}.dBear);
%     end
% end




