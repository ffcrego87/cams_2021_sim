%% Example problem

%% Type of simulation
Nb_test = 0; % 0 - one simulation 1 - several simulations
tr_loc_state = 1; % Transmit local state: 1 - exchange state estimates 0 -exchange measurements
debug_observer = 0; % 0 - use estimates on the controller 1 - use true state on the controller

%% Nb test parameters
Nb_min = 3; % minimum number of bits
Nb_max = 9; % maximum number of bits
n_runs = 25; % number of runs

%% System

dt = 1; % sampling period
dim = 4; % dimension of agent state

A  = kron([1 dt;
           0 1],eye(2));          % Dynamics matrix
       
B  = kron([0;
           dt],eye(2));          % Input matrix
       
%% Noise

% process noise observer
wstd_o = 0.001; % process noise standard deviation

% measurement noise observer
vstd_o = 0.05; % measurement noise standard deviation x-0.01 y-0.05
v_distance_flag_o = 1; % is the noise dependent on distance
dvlstd_o = 0.01; % dvl measurement noise standard deviation

% initial condition observer
xstd_o = 0.1; % Standard deviation of initial position
sstd_o = 0.05; % Standard deviation of initial speed

% process noise
wstd  = 0.001;       % process noise standard deviation

% measurement noise
vstd  = 0.01; % measurement noise standard deviation per distance
v_distance_flag = 1; % is the noise dependent on distance
dvlstd  = 0.01; % dvl measurement noise standard deviation

% initial condition
xstd = 0.1; % Standard deviation of initial position
sstd = 0.05; % Standard deviation of initial speed

%% Nominal positions
p_nominal = [ 0  0  0 -5  5;
             50 55 45 50 50];
p_beacon = [0;
            0];
N_dvl = 1;

%% Trajectory
% vel = [1;0]; % formation velocity
% Freq = 0.02; % Perturbation frequancy
% Amp_pert = 0.2; % Amplitude of frequency perturbation as a fraction of speed; 
% 
% dimp = dim/2;
% R = eye(dimp);
% R(1:2,1:2) = [0 1;-1 0];
% Pert_basis = Amp_pert*R*vel*Freq^-1;
% p0 = p_nominal(:,1);
% 
% pd = @(t) p0+vel*t+Pert_basis*sin(2*pi*Freq*t);
% dpd = @(t) vel+2*pi*Freq*Pert_basis*cos(2*pi*Freq*t);

p0 = p_nominal(:,1);
max_speed = 0.5;
Size_pd = 50;

pd = @(t) p0+Size_pd*[cos((max_speed/Size_pd)*t+pi/2);sin(2*(max_speed/Size_pd)*t+pi)/2];
dpd = @(t) max_speed*[-sin((max_speed/Size_pd)*t+pi/2);cos(2*(max_speed/Size_pd)*t+pi)];

%% Quantizer parameters
alpha=2;
beta=0.5;
Lambda_zero=10;
Nb = 2;
Ndelay = 0;

%% Control gains
Kp=0.2;
Ks=1;

%% Final iteration
ifinal = ceil((Size_pd/max_speed)*4*pi);

%% Formation plot 
T_plot_form = (Size_pd/max_speed)*pi/4;

%% plot ranges for distance objactives
dreps = 0.1;
Distance_ranges = [5-dreps 5*sqrt(2)-dreps 10-dreps; 5+dreps 5*sqrt(2)+dreps 10+dreps];
