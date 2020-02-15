
%==== This M.file declares the tracking problem structure ===
% here we declare a 2-D constant velocity model with coordinate
% measurements
%=== the following parameters is directly assessed by 
% siggen.m and phdfilter.m
K= 100;  %data length
model.cost_method = 'Gibbs';
do_joint = 1;
x_dim= 4;   %dimension of state vector
z_dim= 2;   %dimension of observation vector
model.x_dim = x_dim;
M_max= 4;      %max. number of targets
P_S= .99;  %probability of target survival
P_D= .88;  %probability of detection in measurements
% P_D= .90;  %probability of detection in measurements
miss_on = 1;
Q_S= 1-P_S; %probability of death
Q_D = 1-P_D;   %probability of misdetection
P_S_tempd= 0.90*P_S;  %tempered probability of target survival
P_D_tempd= 0.90*P_D;  %tempered probability of detection in measurements
Q_S_tempd= 1-P_S_tempd; %tempered probability of death
Q_D_tempd= 1-P_D_tempd;   %tempered probability of misdetection
N_max = 10; %max number of terms to calculate the cardinality distribution
%WARNING: DO NOT set N_max  <= max |Z_k| (out of bounds ezrror will occur)
%NOTE:    Advisable to set N_max >> max |Z_k| 
%         and              N_max >> max |X_k| 
%This will ensure capturing the tail of the clutter and state cardinality!!
lambda_c = 66; %make sure lambda_c << N_max !!!
% lambda_c = 1; %make sure lambda_c << N_max !!!
log_lambda_c = log(lambda_c);
cdn_clutter = poisspdf([0:N_max]',lambda_c); %cardinality distribution of clutter - must be Poisson
log_cdn_clutter = log(cdn_clutter+eps(0));
range_c= [ -1000 1000; -1000 1000 ];    %clutter intervals
clutterpdf = 1/prod(range_c(:,2)-range_c(:,1));
log_clutterpdf = log(clutterpdf);
run_flag = 'disp';

% Hbesreq= [5000 5000 5000];  %min num generated internally based on cardinality std dev
Hbesreq= [1000 1000 1000];  %min num generated internally based on cardinality std dev
% Hbesreq= [500 500 500];  %min num generated internally based on cardinality std dev
% Hbesreq= [300,300,300];  %min num generated internally based on cardinality std dev

%=== the following parameters are more problem dependent
% 'siggen', 'phdfilter', etc. do not assess these parameters directly.
% These parameters are used by problem dependent functions such as
%
% gen_newstate, gen_birthstate, gen_observation, compute_likelihood
% 
% we create a structure array 'model' to store these specific parameters

%===here we set up the state spc eqn. x= Ax_old + Bv
T= 1;   %sampling period
A0= [ 1 T; 0 1 ];                       
model.A= [ A0 zeros(2,2); zeros(2,2) A0 ];
F0 = @(dt) [ 1 dt; 0 1 ];   
model.F = @(dt) [ F0(dt) zeros(2,2); zeros(2,2) F0(dt) ];


B0= [ (T^2)/2; T ];
model.B= [ B0 zeros(2,1); zeros(2,1) B0 ];
% model.sigma_v= 5;
model.sigma_v= 1;
model.Q= (model.sigma_v)^2* model.B*model.B';
model.Q_s = model.Q;
model.d_s = zeros(x_dim,1);

%=== parameters for the observation 
model.C_posn= [ 1 0 0 0 ; 0 0 1 0 ]; 
model.C_vely= [ 0 1 0 0 ; 0 0 0 1 ]; 
model.H = model.C_posn;
model.C_vel = [0,1,0,0; 0,0,0,1];
model.J = model.C_vel;
model.D= diag([ 10; 10 ]);   %std for angle and range noise
% model.D= diag([ 5e-0; 5e-0 ]);   %std for angle and range noise
model.R= model.D*model.D';

%===MB BIRTH --- here is the parameter for the birth target states
L_b= 3;         %no. of MB birth terms

estNumBirth = 6;
estNumSpawn = 6;
estNumParentPerTimeStep = 6;
numBirthPerTimestep = 3;
numSpawnPerTimeStepPerParent = 1;
wBirth = estNumBirth/K*(1/numBirthPerTimestep);
wSpawn = estNumSpawn/K*(1/(numSpawnPerTimeStepPerParent*estNumParentPerTimeStep));

% wSpawn = 0.015;


wBirth_tempd = 1.0*wBirth;
wSpawn_tempd = 1.0*wSpawn;


pBirth = diag([10,10,10,10]);

N_birth= zeros(L_b,1);
model.bar_q= zeros(L_b,1);
model.bar_q_tempd= zeros(L_b,1);
model.lambda_b= cell(L_b,1);
model.bar_x= cell(L_b,1);
model.bar_B= cell(L_b,1); 
model.bar_Q= cell(L_b,1);

x0 = 0;
y0 = 0;

d_top = 500;
d =2*d_top*sind(60);
d_bot = d*sind(30)/(2*sind(60));

N_birth(1)= 1;                                     %no. of Gaussians in birth term 1
model.bar_q(1)=wBirth;                            %prob of existence for term 1
model.bar_q_tempd(1)=wBirth_tempd;                      %tempered prob of existence for term 1
model.lambda_b{1}(1)= 1;                        %weight of Gaussians
model.bar_x{1}(:,1)= [ x0; 0; y0+d_top ;0 ];           %mean of Gaussians
% model.bar_x{1}(:,1) = [-400, 0, 0, 0]';
model.bar_B{1}(:,:,1)= pBirth; %std of Gaussians

N_birth(2)= 1;          %no. of Gaussians in birth term 2
model.bar_q(2)=wBirth;                            %prob of existence for term 2
model.bar_q_tempd(2)=wBirth_tempd;                      %tempered prob of existence for term 2
model.lambda_b{2}(1)= 1;                        %weight of Gaussians
model.bar_x{2}(:,1)= [ x0+d/2;0;y0-d_bot;0 ];     %mean of Gaussians
% model.bar_x{2}(:,1)= [ 400,0,0,0 ];
model.bar_B{2}(:,:,1)= pBirth; %std of Gaussians

N_birth(3)= 1;          %no. of Gaussians in birth term 3
model.bar_q(3)=wBirth;                            %prob of existence for term 3
model.bar_q_tempd(3)=wBirth_tempd;                      %tempered prob of existence for term 3
model.lambda_b{3}(1)= 1;                        %weight of Gaussians
model.bar_x{3}(:,1)= [ x0-d/2;0;y0-d_bot;0 ];     %mean of Gaussians
model.bar_B{3}(:,:,1)= pBirth;%std of Gaussians

for i=1:L_b
    for g=1:N_birth(i)
    model.bar_Q{i}(:,:,g)= model.bar_B{i}(:,:,g)*model.bar_B{i}(:,:,g)'; %cov of Gaussians
    end
end

%% ===MB SPAWN --- here is the parameter for the spawn target states
P_T = wSpawn;
Q_T = 1-P_T;
P_T_tempd = wSpawn_tempd;
Q_T_tempd = 1-P_T_tempd;

define_spawnModel;

%% Variables for 'togetherAgain.m'
model.P_S = P_S;
model.P_D = P_D;
model.P_T = P_T;
model.N_max = N_max;
model.z_dim = z_dim;
model.spawn = spawn;
model.lambda_c = lambda_c;
model.clutterIntensity = lambda_c*clutterpdf;
model.K = K;

model.nBirth = 3;
model.nSpawn = 1;
