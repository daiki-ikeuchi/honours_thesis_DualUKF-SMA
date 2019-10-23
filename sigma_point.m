% calculate sigma points from posterior estimates at k-1 time
function Xp_km = sigma_point(xp_km,Pp_km,gamma)

% Initialization
Xp_km = zeros(size(xp_km,1),5);
Pp_km = real(Pp_km);

% Calculate gamma
A = (gamma*sqrt(Pp_km));
Y = xp_km(:,ones(1,numel(xp_km)));

% Correspond to "Num_StateVar" row vs "Num_sigma" column (2 by 5 matrix)
% Xp_km = [xp_km Y+A Y-A]; 
Xp_km(:,1) = xp_km;
Xp_km(:,2:3) = Y+A;
Xp_km(:,4:5) = Y-A;

Xp_km = real(Xp_km);
