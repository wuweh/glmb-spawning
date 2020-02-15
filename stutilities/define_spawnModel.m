spawn.devMag.mean = 70;
spawn.devMag.sigma = 10;
spawn.Az.mean = deg2rad( -90 );
spawn.Az.sigma = deg2rad( 10 ); 


% What it was with pretty good results
xSigma = 5;
ySigma = 5;


p = [xSigma^2,0;0,ySigma^2];
dTheta = acos(1-((2*ySigma)^2/(2*(spawn.devMag.mean^2))));
degSpan = 2*spawn.Az.sigma + 1;

devMag = [spawn.devMag.mean-spawn.devMag.sigma,...
          spawn.devMag.mean,...
          spawn.devMag.mean+spawn.devMag.sigma];

numAngles = ceil(degSpan/dTheta);
numDist   = length(devMag);  
numComps = numAngles*numDist + numDist;

posDev = zeros(x_dim/2,numComps);
velUnitVec = zeros(x_dim/2,numComps);
Qpos = zeros(x_dim/2,x_dim/2,numComps);
Qvel = Qpos;
     
count = 0;
for rngidx = 1:length(devMag)
    dev = devMag(rngidx);
    
    for idx = 1:numAngles/2
        count = count + 1;
        theta = spawn.Az.mean + idx*dTheta;
        R = [cos(theta),-sin(theta);sin(theta),cos(theta)];
        posDev(:,count) = R*[dev,0]';
        velUnitVec(:,count) = R*[1,0]';
        q = R*p*R';
        Qpos(:,:,count) = q;
        Qvel(:,:,count) = q;
    end
    for idx = 1:numAngles/2
        count = count + 1;
        theta = spawn.Az.mean - idx*dTheta;
        R = [cos(theta),-sin(theta);sin(theta),cos(theta)];
        posDev(:,count) = R*[dev,0]';
        velUnitVec(:,count) = R*[1,0]';
        q = R*p*R';
        Qpos(:,:,count) = q;
        Qvel(:,:,count) = q;
    end
    
    count = count + 1;
    theta = spawn.Az.mean;
    R = [cos(theta),-sin(theta);sin(theta),cos(theta)];
    posDev(:,count) = R*[dev,0]';
    velUnitVec(:,count) = R*[1,0]';
    q = R*p*R';
    Qpos(:,:,count) = q;
    Qvel(:,:,count) = q;
end

numComps = size(posDev,2);

spawn.velMult = 0;%1;
spawn.w{1} = ones(1,numComps)/numComps;
spawn.posDev{1} = posDev;
spawn.Qpos{1} = Qpos;
spawn.Qvel{1} = Qvel;
spawn.velUnitVec{1} = velUnitVec;


