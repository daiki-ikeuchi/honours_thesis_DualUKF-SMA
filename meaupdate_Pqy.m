% Cross covariance between a piori state estimate and measurement estimate at k time
function Pqy_k = meaupdate_Pqy(Wc_q,Qm_k,qm_k,Ym_k,ym_k,Num_sigma,Num_ObsVar,Num_Param)

% Wc = weights for covariance (5 by 1)
% Ym_k = a priori sigma points at k time (1 by 5) where row = resistance
% ym_k = a priori estimate of observed variable (resistance) (1 by 1)
% Rv = measurement noise (1 by 1)

% Initialize
Pqy_k = zeros(Num_Param,Num_ObsVar);

% Calculate the a priori error covariance matrix (Pyy)
for i = 1:Num_sigma
    
    % Sum for the number of sigma points
    Pqy_k = Pqy_k + Wc_q(i,1)*(Qm_k(:, i) - qm_k)*(Ym_k(:, i) - ym_k)';
    
end

Pqy_k = real(Pqy_k);
