function zeta_dot = zeta_dot(Xp_km,X_dot_km,Num_sigma)

% Initialization
zeta_dot = zeros(1,Num_sigma);

% Define the parameters
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

for i = 1:Num_sigma
    
    % Heating
    if X_dot_km(1,i) > 0 && Xp_km(1,i) > Ta_s && Xp_km(1,i) < Ta_f
        
        zeta_dot(1,i) = (-0.5*zeta0)*(sin(aA*(Xp_km(1,i)-Ta_s)+bA*Xp_km(2,i))*(aA*X_dot_km(1,i) + bA*X_dot_km(2,i)));
    
    % Cooling    
        if X_dot_km(1,i) < 0 && Xp_km(1,i) > Tm_f && Xp_km(1,i) < Tm_s
        zeta_dot(1,i) = (-0.5*zeta0)*(sin(aM*(Xp_km(1,i)-Tm_f)+bM*Xp_km(2,i))*(aM*X_dot_km(1,i) + bM*X_dot_km(2,i)));  
        end
        
    else
        zeta_dot(1,i) = 0;
    end
end
    
        
        
        