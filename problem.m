%% Example problem

%% System

dt= 1; % sampling period
speed = [1;0]; % formation speed
dim = 4; % dimension of agent state
Freq = 0.02; % Perturbation frequancy
Amp_pert = 0.2; % Amplitude of frequency perturbation as a fraction of speed; 

A  = kron([1 dt;
           0 1],eye(2));          % Dynamics matrix
       
B  = kron([0;
           dt],eye(2));          % Input matrix

debug_observer = 0;
       
%% Noise

% process noise observer
wstd_o  = 0.05;       % process noise standard deviation

% measurement noise observer
vstd_o  = 0.01; % measurement noise standard deviation
v_distance_flag_o = 1; % is the noise dependent on distance
dvlstd_o  = 0.05; % dvl measurement noise standard deviation

% initial condition observer
xstd_o = 0.01; % Standard deviation of initial position
sstd_o = 0.001; % Standard deviation of initial speed

% process noise
wstd  = 0.01;       % process noise standard deviation

% measurement noise
vstd  = 0.0005; % measurement noise standard deviation per distance
v_distance_flag = 1; % is the noise dependent on distance
dvlstd  = 0.05; % dvl measurement noise standard deviation

% initial condition
xstd = 0.5; % Standard deviation of initial position
sstd = 0.1; % Standard deviation of initial speed

tr_loc_state = 0; % Transmit local state

%% Nominal positions
p_nominal = [ -100 -100 -100 -105 -95;
                50   55   45   50  50];
p_beacon = [0 400;
            0   0];
N_dvl = 1;

%% Quantizer parameters
alpha=2;
beta=0.5;
Lambda_zero=10;
Nb = 3;

%% Control gains
Kp=0.5;
Ks=1;

%% Final iteration
ifinal = 201;

%% Formation plot 
T_plot_form = 50;

%% plot ranges for distance objactives
dreps = 0.1;
Distance_ranges = [5-dreps 5*sqrt(2)-dreps 10-dreps; 5+dreps 5*sqrt(2)+dreps 10+dreps];
