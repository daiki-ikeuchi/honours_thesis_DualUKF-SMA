%% SMA wire actuator modelling

clear
close all
clc


% State estimate
xp_km_in = [20;7.0736e6];
Pp_km_in = eye(2);
X_dot_km_in = [0 0 0 0 0;0 0 0 0 0];


% Parameter estimate
Ks_in = 99;
qp_km_in = Ks_in;
Pqp_km_in = eye(1);


% Geometry of SMA wire
d = 0.6e-3;        % Diameter of SMA wire (m)
L0 = 0.25;      % Initial length of SMA wire (m)

% Constitutive model
Da = 75e9;      % Modulus of elasticity at austenite (Pa)
Dm = 28e9;      % Modulus of elasticity at martensite (Pa)
D = (Da+Dm)/2;  % Average modulus of elasticity (Pa)
eps_l = 0.067;   % Maximum residual strain
theta_T = 550000; % Thermal coefficient of expansion of SMA (Pa/degC)


% Heat transfer model
m = 1.14e-4;    % Mass per unit length (kg/m)
Cp = 840;       % Specific heat (J/kg*degC)  Note: check if its Kelvin
R = 0.8491;     % Electrical resistance per unit length (ohm/m)  Note: How to deal with in self-sensing
As = pi*d*L0;   % Circumferential surface area for convection (m^2)
T_amb = 20;     % Ambient temperature (degC)
h0 = 120;     % Convectivity (W/(m^2)*degC) 
h2 = 0.001;   % h-temperature dependent factor (W/(m^2)*((degC)^3))


% Phase transformation model
Ta_s = 68;      % Austenite phase start temperature (degC)
Ta_f = 78;      % Austenite phase finish temperature (degC)
Tm_s = 52;      % Martensite phase start temperature (degC)
Tm_f = 42;      % Martensite phase finish temperature (degC)

aA = pi/(Ta_f-Ta_s);    % Curve fitting parameter (1/degC)
Ca = 12000000;  % Austenite transformation constant (Pa/degC)
bA = -aA/Ca;    % Curve fitting parameter (1/Pa)

aM = pi/(Tm_s-Tm_f);    % Curve fitting parameter (1/degC)
Cm = 10000000;  % Martensite transformation constant (Pa/degC)
bM = -aM/Cm;    % Curve fitting parameter (1/Pa)
zeta0 = 1;      % SMA initial martensite fraction factor


% SMA wire actuator dynamics
Ac = (pi*d^2)/4;    % Cross sectional area of SMA wire (m^2)
Ks = 100;           % Spring constant (N/m)
disp_i = 0.02;      % Initial displacement (m) prior to voltage application
eps_i = disp_i/L0;  % Initial strain prior to voltage application


% Resistance dynamics
rho_a = 8.37e-7;    % Resistivity of austenite phase (ohm*m)
rho_m = 9.603e-7;   % Resistivity of martensite phase (ohm*m)







