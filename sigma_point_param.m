% calculate sigma points from posterior estimates at k-1 time

function Qm_k = sigma_point_param(qp_km,Pqp_km,gamma_q)

% Initialization
Qm_k = zeros(size(qp_km,1),5);
qp_km = real(qp_km);
Pqp_km = real(Pqp_km);
gamma_q = real(gamma_q);

% Calculate gamma
A = real((gamma_q*sqrt(Pqp_km)));
Y = qp_km(:,ones(1,numel(qp_km)));

% Correspond to "Num_Param" row vs "Num_sigma" column (1 by 5 matrix)
% Qp_km = [qp_km Y+A Y-A]; 
Qm_k(:,1) = qp_km;
Qm_k(:,2:3) = Y+A;
Qm_k(:,4:5) = Y-A;

Qm_k = real(Qm_k);