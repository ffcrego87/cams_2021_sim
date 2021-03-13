%% Example problem

%% System

dt= 1; % sampling period
speed = [1;0]; % formation speed
dim = 4; % dimension of agent state
Freq = 0.02; % Perturbation frequancy
Amp_pert = 0.5; % Amplitude of frequency perturbation as a fraction of speed; 

A  = kron([1 dt;
           0 1],eye(2));          % Dynamics matrix
       
B  = kron([0;
           dt],eye(2));          % Input matrix

debug_observer = 0;
       
%% Noise

% process noise observer
wstd_o  = 0.1;       % process noise standard deviation

% measurement noise observer
vstd_o  = 0.05; % measurement noise standard deviation
dvlstd_o  = 0.005; % dvl measurement noise standard deviation

% process noise
wstd  = 0.01;       % process noise standard deviation

% measurement noise
vstd  = 0.1; % measurement noise standard deviation per distance
dvlstd  = 0.05; % dvl measurement noise standard deviation

% initial condition
xstd = 2; % Standard deviation of initial position
sstd = 0.1; % Standard deviation of initial speed


%% Nominal positions
p_nominal = [ -20 -20 -20 -25 -15;
               50  55  45  50  50];
p_beacon = [0;
            0];
N_dvl = 1;

%% Control gains
Kp=0.5;
Ks=1;

%% Final iteration
ifinal = 201;

%% Formation plot 
T_plot_form = 50;
