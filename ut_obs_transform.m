% Unscented utransformation of sigma points using observation function
function Ym_k = ut_obs_transform(Xm_k,zeta,Num_ObsVar,qp_km)

% Xm_k = a priori sigma points at k time before measurement update (2 by 5)
% zeta = phase fraction at k-1 time
% Num_ObsVar = number of observation variable (just one for resistance)

% Ym_k = unscented transformed sigma points using observation function at k
% time


% Initialization
Ym_k = zeros(Num_ObsVar,size(Xm_k,2));

%% Resistance observation function update
% Relevant parameters
d = 0.6e-3;         % Diameter of SMA wire (m)
rho_a = 8.37e-7;    % Resistivity of austenite phase (ohm*m)
rho_m = 9.603e-7;   % Resistivity of martensite phase (ohm*m)
L0 = 0.25;          % Initial length of SMA wire (m)
disp_i = 0.02;      % Initial displacement (m) prior to voltage application
eps_i = disp_i/L0;  % Initial strain prior to voltage application
Ac = (pi*d^2)/4;    % Cross sectional area of SMA wire (m^2)
Ks = real(qp_km);           % Spring constant (N/m)


% Resistance unscented transform
for i = 1:size(Xm_k,2)
    
    % Resistance model in terms of state variables
    Ym_k(1,i) = (((rho_m*zeta) + (1-zeta)*rho_a)*L0*(1 + (eps_i - ((Xm_k(2,i)*Ac)/(Ks*L0)))))/Ac;
    
end

Ym_k = real(Ym_k);

