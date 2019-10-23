% Cross covariance between a piori state estimate and measurement estimate at k time
function Pxy_k = meaupdate_Pxy(Wc,Xm_k,xm_k,Ym_k,ym_k,Num_sigma,Num_ObsVar,Num_StateVar)

% Wc = weights for covariance (5 by 1)
% Ym_k = a priori sigma points at k time (1 by 5) where row = resistance
% ym_k = a priori estimate of observed variable (resistance) (1 by 1)
% Rv = measurement noise (1 by 1)

% Initialize
Pxy_k = zeros(Num_StateVar,Num_ObsVar);

% Calculate the a priori error covariance matrix (Pyy)
for i = 1:Num_sigma
    
    % Sum for the number of sigma points
    Pxy_k = Pxy_k + Wc(i,1)*(Xm_k(:, i) - xm_k)*(Ym_k(:, i) - ym_k)';
    
end


Pxy_k = real(Pxy_k);