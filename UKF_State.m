function [xp_k,Pp_k,X_dot_km] = UKF_State(xp_km,Pp_km,u,Res,Zeta,X_dot_km,qp_km)

u = real(u);
Res = real(Res);
Zeta = real(Zeta);
X_dot_km = real(X_dot_km);
qp_km = real(qp_km);

% Model parameters
Num_StateVar = real(2);      %number of state variables
Num_ObsVar = real(1);        %number of observation variables
Num_sigma = real(2*Num_StateVar + 1);      % Number of sigma points 

t_sample = real(1e-5);


% Tuning parameters
alpha = real(0.01);       % Tuning parameter to define the spread of sigma points
beta = real(2);           % Tuning parameter for distribution of state
kappa = real(1);          % Tuning parameter for scaling of sigma points           
lambda = real(((alpha^2)*(Num_StateVar + kappa)) - Num_StateVar); % Composite scaling parameter
gamma = real(sqrt(real(Num_StateVar + lambda)));       % scaling factor


% Process and Measurement noise characteristics
w1 = real(1e-10);
w2 = real(1000);
Qv = real([w1*w1 0;0 w2*w2]);  % covariance of process noise (from UKF Gurung [14])
Rv = 0;        % covariance of measurement noise (from UKF Gurung [14])


% Weighting terms (Num_sigma by 1) or (5 vs 1)
Wm = zeros(Num_sigma,1);    % Weights on mean
Wc = zeros(Num_sigma,1);    % Weights on covariance

Wm(1,1) = real(lambda/(Num_StateVar+lambda));        % weight for i=0
Wm(2:end,1) = real(1/(2*(Num_StateVar+lambda)));     % weight for i = 1 to 2*Num_sigma

Wc(1,1) = real((lambda/(Num_StateVar+lambda))+(1 - alpha^2 +beta)); % weight for i=0 (1 - alpha^2 +beta)
Wc(2:end,1) = real(1/(2*(Num_StateVar+lambda)));     % weight for 1 = 1 to 2*Num_sigma

xp_km = real(xp_km);
Pp_km = real(Pp_km);


%% UKF state estimation

% Calculate sigma point at k-1 time after measurement update
Xp_km = real(sigma_point(xp_km,Pp_km,gamma));


%% Time update 

% Unscented transform the sigma at k time before measurement update
u_km = u;                % voltage input at k-1 time
Res_km = Res;               % Actual measurement (taken from open-loop model)
[Xm_k,X_dot_km] = ut_transform(Xp_km,u_km,t_sample,Res_km,X_dot_km,Num_sigma,qp_km);
X_dot_km = real(X_dot_km);

% A priori state estimate at k time
xm_k = real(timeupdate_state(Wm,Xm_k,Num_StateVar,Num_sigma));

% A priori error covariance estimate at k time (2 by 2 matrix)
Pm_k = real(timeupdate_cov(Wc,Xm_k,xm_k,Qv,Num_sigma,Num_StateVar));

% Unscented transform the sigma points using observation function
zeta = real(Zeta);     % phase fraction at k-1 time
Ym_k = real(ut_obs_transform(Xm_k,zeta,Num_ObsVar,qp_km));

% A priori measurement estimate at k time
ym_k = real(timeupdate_obs(Wm,Ym_k,Num_ObsVar,Num_sigma));


%% Measurement update

% Covariance of the estimated measurement at k time
Pyy_k = real(meaupdate_Pyy(Wc,Ym_k,ym_k,Rv,Num_sigma,Num_ObsVar));

% Cross covariance between a piori state estimate and measurement
% estimate at k time
Pxy_k = real(meaupdate_Pxy(Wc,Xm_k,xm_k,Ym_k,ym_k,Num_sigma,Num_ObsVar,Num_StateVar));

% Calculate the Kalman gain
K_k = real(Pxy_k/Pyy_k);

% State update after measurement
Res_k = Res;               % Actual measurement (taken from open-loop model)
xp_k = real(xm_k + K_k*(Res_k - ym_k));
xp_k = real(xp_k);

% Covariance update after measurement
Pp_k = real(Pm_k - (K_k*Pyy_k*K_k'));
Pp_k = real(Pp_k);





