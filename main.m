% Main file for UKF state estimation
clear
close all
clc

%% main file for UKF state and parameter estimation

% Model parameters
Num_StateVar = 2;      %number of state variables
Num_ObsVar = 1;        %number of observation variables
Num_sigma = 2*Num_StateVar + 1;      % Number of sigma points
Num_Param = 1;          % number of parameters to be estimated



% Tuning parameters for state estimation
alpha = 0.01;       % Tuning parameter to define the spread of sigma points
beta = 2;           % Tuning parameter for distribution of state 
kappa = 1;          % Tuning parameter for scaling of sigma points            
lambda = ((alpha^2)*(Num_StateVar + kappa)) - Num_StateVar; % Composite scaling parameter
gamma = sqrt(Num_StateVar + lambda);       % scaling factor


% Tuning parameters for parameter estimation
alpha_q = 0.1;       % Tuning parameter to define the spread of sigma points
beta_q = 2;           % Tuning parameter for distribution of state
kappa_q = 0;          % Tuning parameter for scaling of sigma points          
lambda_q = ((alpha_q^2)*(Num_StateVar + kappa_q)) - Num_StateVar; % Composite scaling parameter
gamma_q = sqrt(Num_StateVar + lambda_q);       % scaling factor



% Load the data from open-loop model
load('Temperature.mat');        % Load temperature data (state variable 1)
load('Stress.mat');             % Load stress data (state variable 2)
load('Resistance.mat');         % Load resistance data (observation variable)
load('Voltage.mat');            % Load voltage data (input variable)
load('Zeta_dot.mat');           % Load phase transformation data
load('Zeta.mat');               % Load phase transformation
load('Displacement.mat');       % Load phase transformation
load('Temp_rate.mat');          % Load temperature gradient transformation
load('Stress_rate.mat');        % Load stress gradient transformation
time = T.time;                  % Extract time array
t_sample = time(end-10) - time(end-11);     % Sampling time
T = T.Data;                     % Extract temperature array
S = S.Data;                     % Extract stress array
Res = Res.Data;                 % Extract resistance array
V = V.Data;                     % Extract voltage array
Zeta_dot = Zeta_dot.Data;       % Extract phase transformation rate array
Zeta = Zeta.Data;               % Extract phase fraction array
Dis = Dis.Data;                 % Extract actuator displacement
temp_dot = T_dot.Data;                 % Extract actuator displacement
stress_dot = S_dot.Data;                 % Extract actuator displacement



% Process and Measurement noise characteristics of states
w1 = 1e-5;
w2 = 10000000;
Qv = [w1*w1 0;0 w2*w2];  % covariance of process noise
Rv = 1e-20;        % covariance of measurement noise

W1 = 0*rand(length(time),1);
W2 = 0*rand(length(time),1);

noise_param = 0.00001;
w_param = -noise_param + (noise_param - (-noise_param))*rand(length(time),1);


% Process and Measurement noise characteristics of parameters
Qq = 1e-6;      
Rq = 5e-7;



% Initial conditions of states
T0 = T(1,1);         % Initial temperature (C)
S0 = S(1,1);         % Initial stress (Pa)
xp_km = [T0;S0];     % initial state vector
Pp_km = eye(Num_StateVar);     % initial state covraiance



% Initial conditions of parameters
Ks_in = 95;
qp_km = Ks_in;
Pqp_km = eye(Num_Param);



% Weighting terms for states
Wm = zeros(Num_sigma,1);    % Weights on mean
Wc = zeros(Num_sigma,1);    % Weights on covariance

Wm(1,1) = lambda/(Num_StateVar+lambda);        % weight for i=0
Wm(2:end,1) = 1/(2*(Num_StateVar+lambda));     % weight for i = 1 to 2*Num_sigma

Wc(1,1) = (lambda/(Num_StateVar+lambda))+(1 - alpha^2 +beta); % weight for i=0 (1 - alpha^2 +beta)
Wc(2:end,1) = 1/(2*(Num_StateVar+lambda));     % weight for 1 = 1 to 2*Num_sigma


% Weighting terms for parameters (Num_sigma by 1) or (5 vs 1)
Wm_q = zeros(Num_sigma,1);    % Weights on mean
Wc_q = zeros(Num_sigma,1);    % Weights on covariance

Wm_q(1,1) = lambda_q/(Num_StateVar+lambda_q);        % weight for i=0
Wm_q(2:end,1) = 1/(2*(Num_StateVar+lambda_q));     % weight for i = 1 to 2*Num_sigma

Wc_q(1,1) = (lambda_q/(Num_StateVar+lambda_q))+(1 - alpha_q^2 +beta_q); % weight for i=0 (1 - alpha^2 +beta)
Wc_q(2:end,1) = 1/(2*(Num_StateVar+lambda_q));     % weight for 1 = 1 to 2*Num_sigma



% Storing matrix for state variables
state_estimate = zeros(Num_StateVar,length(time));
state_estimate(:,1) = xp_km;    % Enforce the initial conditions


% Storing matrix for estimated parameter
param_estimate = zeros(Num_Param,length(time));
param_estimate(:,1) = qp_km;    % Enforce the initial condition


% Storing matrix for observation variable
obs_estimate = zeros(Num_ObsVar,length(time));
obs_estimate(1,1) = Res(1,1);   % Enforce the initial measurement

obs_estimate_param = zeros(Num_Param,length(time));
obs_estimate_param(1,1) = Res(1,1);   % Enforce the initial measurement



% Storing matrix for erate variables of sigma points
X_dot_km = zeros(Num_StateVar,Num_sigma);
X_dot_km(1,:) = temp_dot(1,1);
X_dot_km(2,:) = stress_dot(1,1);


% Storing matrix for erate variables of sigma points
x_dot_km = zeros(Num_StateVar,1);
x_dot_km(1,1) = temp_dot(1,1);
x_dot_km(2,1) = stress_dot(1,1);




% Main UKF Solver
for i = 2:length(time)
    
    
    %% State estimation
    % Calculate sigma point at k-1 time after measurement update
    Xp_km = sigma_point(xp_km,Pp_km,gamma);

    
    %%%%%%%%%%%%%%%%%%%%%
    % Time update 
    %%%%%%%%%%%%%%%%%%%%%
    
    % Unscented transform the sigma at k time before measurement update
    u_km = V(:,:,i-1);                % voltage input at k-1 time
    Res_km = Res(i-1,1);               % Actual measurement (taken from open-loop model)
    W1_k = W1(i,1);
    W2_k = W2(i,1);
    [Xm_k,X_dot_km] = ut_transform(Xp_km,u_km,t_sample,Res_km,X_dot_km,Num_sigma,W1_k,W2_k,qp_km);
    
    
    % A priori state estimate at k time
    xm_k = timeupdate_state(Wm,Xm_k,Num_StateVar,Num_sigma);
    
    
    % A priori error covariance estimate at k time
    Pm_k = timeupdate_cov(Wc,Xm_k,xm_k,Qv,Num_sigma,Num_StateVar);
    
    
    % Unscented transform the sigma points using observation function
    zeta = Zeta(i-1,1);     % phase fraction at k-1 time
    Ym_k = ut_obs_transform(Xm_k,zeta,Num_ObsVar,qp_km);
    
    
    % A priori measurement estimate at k time
    ym_k = timeupdate_obs(Wm,Ym_k,Num_ObsVar,Num_sigma);
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Measurement update
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Covariance of the estimated measurement at k time
    Pyy_k = meaupdate_Pyy(Wc,Ym_k,ym_k,Rv,Num_sigma,Num_ObsVar);
    
    % Cross covariance between a piori state estimate and measurement
    % estimate at k time
    Pxy_k = meaupdate_Pxy(Wc,Xm_k,xm_k,Ym_k,ym_k,Num_sigma,Num_ObsVar,Num_StateVar);
    
    % Calculate the Kalman gain
    K_k = Pxy_k/Pyy_k;
    
    % State update after measurement
    Res_k = Res(i,1);               % Actual measurement (taken from open-loop model)
    xp_k = real(xm_k + K_k*(Res_k - ym_k));
    
    % Covariance update after measurement
    Pp_k = real(Pm_k - (K_k*Pyy_k*K_k'));
    
    
    % Store the estimated states
    state_estimate(:,i) = xp_k;
    obs_estimate(:,i) = ym_k;
    
    
    
    %% Parameter estimation
    
    %%%%%%%%%%%%%%%%%%%
    % Time update
    %%%%%%%%%%%%%%%%%%%
    
   
    % Calculate sigma point at k-1 time after measurement update
    Qm_k = sigma_point_param(qp_km,Pqp_km,gamma_q);
    
       
    % A priori parameter estimate at k time
    qm_k = qp_km + w_param(i,1);
    Pqm_k = Pqp_km + Qq;
    
    
    % Unscented transform the sigma points using observation function
    zeta = Zeta(i-1,1);     % phase fraction at k-1 time
    Ym_k = ut_transform_obs_param(xp_km,Qm_k,Num_ObsVar,Num_sigma,zeta);
    
    
    % A priori measurement estimate at k time
    ym_k = timeupdate_obs_param(Wm_q,Ym_k,Num_ObsVar,Num_sigma);
    

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Measurement update
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Covariance of the estimated measurement at k time
    Pyy_k = meaupdate_Pyy_param(Wc_q,Ym_k,ym_k,Rq,Num_sigma,Num_ObsVar);
    
    
    % Cross covariance between a piori state estimate and measurement
    % estimate at k time
    Pqy_k = meaupdate_Pqy(Wc_q,Qm_k,qm_k,Ym_k,ym_k,Num_sigma,Num_ObsVar,Num_Param);
    
    
    % Calculate the Kalman gain
    K_k = Pqy_k/Pyy_k;
    
    
    % Parameter update after measurement
    Res_k = Res(i,1);               % Actual measurement (taken from open-loop model)
    qp_k = qm_k + K_k*(Res_k - ym_k);
    
    
    % Covariance update after measurement
    Pqp_k = Pqm_k - (K_k*Pyy_k*K_k');
    
    
    
    %% Post processing
    
    % Store the estimated states
    param_estimate(:,i) = qp_k;
    obs_estimate_param(1,i) = ym_k;
       
    % Update the estimated state derivative
    x_dot_km = (xp_k - xp_km)/t_sample;
       
    % Update the state and covariance matrix for next iteration
    xp_km = real(xp_k);
    Pp_km = real(Pp_k);
          
    % Update the state and covariance matrix for next iteration
    qp_km = qp_k;
    Pqp_km = Pqp_k;
    
      
end


%% Plot    

% T_real = real(state_estimate(1,:));
% T_img = imag(state_estimate(1,:));
% TT = sqrt((T_real.^2) + (T_img.^2));
% figure(1)
% plot(time,TT);
% hold on
% plot(time,T(:,1));
% 
% 
% figure(1)
% plot(time,state_estimate(1,:));
% hold on
% plot(time,T(:,1));
% 
% 
% S_real = real(state_estimate(2,:));
% SS = S_real;
% S_img = imag(state_estimate(2,:));
% SS = sqrt(S_real.^2 + S_img.^2);
% figure(2)
% plot(time,SS);
% hold on
% plot(time,S(:,1));
% 
% 
% figure(2)
% plot(time,state_estimate(2,:));
% hold on
% plot(time,S(:,1));
%     
% jj = state_estimate(1,:);
% oo = real(jj);
% oo = -1*oo+20;
% tt = imag(jj);
% pp = ((oo.^2)+(tt.^2)).^0.5;
% figure(3)
% plot(time,oo);
% hold on
% plot(time,T(:,1));
% axis([0 time(end,1) 0 100])    
% 
% 
% Resistance
% figure(4)
% plot(time,obs_estimate');
% hold on
% plot(time,Res(:,1));
% hold on
% plot(time,obs_estimate_param');
% 
%   
% Displacement
% d = 0.6e-3;        % Diameter of SMA wire (m)
% Ac = (pi*d^2)/4;    % Cross sectional area of SMA wire (m^2)
% Ks = 100;           % Spring constant (N/m)
% Dis_est = (SS*Ac)/Ks;
% 
% figure(5)
% plot(time,Dis_est);
% hold on
% plot(time,Dis(:,1));
% 
% figure(6)
% Dis_diff = abs(Dis_est-Dis(:,1)');
% plot(time,Dis_diff);
% 
% 
% Input Voltage
% figure(7)
% plot(time,V(1,:));
% 
% 
% Input Voltage
% figure(8)
% plot(time,param_estimate(1,:));
