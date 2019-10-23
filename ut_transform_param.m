function Qm_k = ut_transform_param(Qp_km,u_km,Res_km,x_dot_km,xp_km,Num_sigma,Num_Param)

% Xm_k = Sigma points at k time before measurement update
% Xp_km = Sigma points at k-1 time (2 (Num_StateVar) by 5 (Num_sigma) matrix)
% u_km = control input or voltage at k-1 time (1 by 1 scalar)
% zeta_dot = phase transformation rate at k-1 time

% Initialize
Qm_k = zeros(Num_Param,Num_sigma);



%% Temperature process function update
% Relevant parameters
d = 0.6e-3;        % Diameter of SMA wire (m)
L0 = 0.25;      % Initial length of SMA wire (m)
Relec = Res_km;     % Electrical resistance per unit length (ohm/m) 
h0 = 120;     % Convectivity (W/(m^2)*degC) 
As = pi*d*L0;   % Circumferential surface area for convection (m^2)
T_amb = 20;     % Ambient temperature (degC)
m = 1.14e-4;    % Mass per unit length (kg/m)
Cp = 840;       % Specific heat (J/kg*degC)  Note: check if its Kelvin


% Temperature unscented transform
for i = 1:Num_sigma
    
    % Calculate partially the temperature gradient
    AA = (((u_km^2)/Relec)-((h0 + h2*Xp_km(1,i)*Xp_km(1,i))*As*(Xp_km(1,i) - T_amb)));
    
    % Calculate the temperature gradient at k-1
    Tdot = AA/(m*Cp);
    
    % Record the a priori temperature rate of sigma points
    X_dot_km(1,i) = Tdot;
    
    % Unscented transform the temperature process function
    Xm_k(1,i) = Xp_km(1,i) + Tdot*t_sample + W1_k;
    

end
   

    
    
    
    