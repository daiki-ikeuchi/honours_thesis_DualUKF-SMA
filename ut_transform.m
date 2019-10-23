function [Xm_k,X_dot_km] = ut_transform(Xp_km,u_km,t_sample,Res_km,X_dot_km,Num_sigma,qp_km)

% Xm_k = Sigma points at k time before measurement update
% Xp_km = Sigma points at k-1 time (2 (Num_StateVar) by 5 (Num_sigma) matrix)
% u_km = control input or voltage at k-1 time (1 by 1 scalar)
% zeta_dot = phase transformation rate at k-1 time

% Initialize
Xm_k = zeros(size(Xp_km,1),size(Xp_km,2));

% Assign a priori phase transformation rate
Zeta_dot = real(zeta_dot(Xp_km,X_dot_km,Num_sigma));


%% Temperature process function update
% Relevant parameters
d = 0.6e-3;        % Diameter of SMA wire (m)
L0 = 0.25;      % Initial length of SMA wire (m)
Relec = Res_km;     % Electrical resistance per unit length (ohm/m) 
h0 = 120;     % Convectivity (W/(m^2)*degC) 
h2 = 0.001;   % h-temperature dependent factor (W/(m^2)*((degC)^3))
As = pi*d*L0;   % Circumferential surface area for convection (m^2)
T_amb = 20;     % Ambient temperature (degC)
m = 1.14e-4;    % Mass per unit length (kg/m)
Cp = 840;       % Specific heat (J/kg*degC)  Note: check if its Kelvin


% Temperature unscented transform
for i = 1:size(Xp_km,2)
    
    % Calculate partially the temperature gradient
    AA = (((u_km^2)/Relec)-((h0 + h2*Xp_km(1,i)*Xp_km(1,i))*As*(Xp_km(1,i) - T_amb)));
    
    % Calculate the temperature gradient at k-1
    Tdot = AA/(m*Cp);
    
    % Record the a priori temperature rate of sigma points
    X_dot_km(1,i) = Tdot;
    
    % Unscented transform the temperature process function
    Xm_k(1,i) = Xp_km(1,i) + Tdot*t_sample;
    

end
   


%% Stress process function update

% Relevant parameters
theta_T = 550000; % Thermal coefficient of expansion of SMA (Pa/degC)
Da = 75e9;      % Modulus of elasticity at austenite (Pa)
Dm = 28e9;      % Modulus of elasticity at martensite (Pa)
D = (Da+Dm)/2;  % Average modulus of elasticity (Pa)
eps_l = 0.067;   % Maximum residual strain    
Ac = (pi*d^2)/4;    % Cross sectional area of SMA wire (m^2)
Ks = qp_km;           % Spring constant (N/m)


% Stress unscented transform
for j = 1:size(Xp_km,2)
    
    % Calculate partially the temperature gradient
    AA = (((u_km^2)/Relec)-((h0 + h2*Xp_km(1,i)*Xp_km(1,i))*As*(Xp_km(1,i) - T_amb)));
    
    % Calculate the temperature gradient at k-1
    Tdot = AA/(m*Cp);
    
    % Calculate the stress gradient at k-1
    Sdot = ((theta_T*Tdot)-(eps_l*D*Zeta_dot(1,i)))/(1 + ((Ac*D)/(Ks*L0)));
    
    % Record the a priori temperature rate of sigma points
    X_dot_km(2,i) = Sdot;
    
    % Unscented transform the stress process function
    Xm_k(2,j) = Xp_km(2,j) + Sdot*t_sample;
    
end
    
Xm_k = real(Xm_k);    
X_dot_km = real(X_dot_km);
    
    
    
    