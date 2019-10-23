function Ym_k = ut_transform_obs_param(xp_km,Qm_k,Num_ObsVar,Num_sigma,zeta)


% Initialize
Ym_k = zeros(Num_ObsVar,Num_sigma);
Qm_k = real(Qm_k);


%% Process function for h2 (actual h2 = 0.001)

% Relevant parameters
d = 0.6e-3;         % Diameter of SMA wire (m)
rho_a = 8.37e-7;    % Resistivity of austenite phase (ohm*m)
rho_m = 9.603e-7;   % Resistivity of martensite phase (ohm*m)
L0 = 0.25;          % Initial length of SMA wire (m)
disp_i = 0.02;      % Initial displacement (m) prior to voltage application
eps_i = disp_i/L0;  % Initial strain prior to voltage application
Ac = (pi*d^2)/4;    % Cross sectional area of SMA wire (m^2)



for i = 1:Num_sigma
    
    % Resistance model in terms of state variables
    Ym_k(1,i) = (((rho_m*zeta) + (1-zeta)*rho_a)*L0*(1 + (eps_i - ((xp_km(2,1)*Ac)/(Qm_k(1,i)*L0)))))/Ac;
    
end
 

Ym_k = real(Ym_k);
    
    
    
    
    
    
    
    
    
    
