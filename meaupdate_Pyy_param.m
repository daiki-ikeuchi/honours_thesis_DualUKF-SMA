% Measurement update: Covariance of the estimated measurement at k time
function Pyy_k = meaupdate_Pyy_param(Wc_q,Ym_k,ym_k,Rq,Num_sigma,Num_ObsVar)

% Wc = weights for covariance (5 by 1)
% Ym_k = a priori sigma points at k time (1 by 5) where row = resistance
% ym_k = a priori estimate of observed variable (resistance) (1 by 1)
% Rv = measurement noise (1 by 1)

% Initialize
Pyy_k = zeros(Num_ObsVar);

% Calculate the a priori error covariance matrix (Pyy)
for i = 1:Num_sigma
    
    % Sum for the number of sigma points
    Pyy_k = Pyy_k + Wc_q(i,1)*(Ym_k(:, i) - ym_k)*(Ym_k(:, i) - ym_k)';
    
end

% Add the process noise
Pyy_k = real(Pyy_k + Rq);


