% A priori error covariance estimate at k time
function Pm_k = timeupdate_cov(Wc,Xm_k,xm_k,Qv,Num_sigma,Num_StateVar)

% Wc = weights for covariance (5 by 1)
% Xm_k = a priori sigma points at k time (2 by 5) where row 1 = temp and
% row 2 - stress

% xm_k = a priori estimate of state variables (2 by 1) where row 1 = temp and
% row 2 - stress

% Qv = process noise (2 by2)

% Initialize
Pm_k = zeros(Num_StateVar);

% Calculate the a priori error covariance matrix
for i = 1:Num_sigma
    
    % Sum for the number of sigma points
    Pm_k = Pm_k + Wc(i,1)*(Xm_k(:, i) - xm_k)*(Xm_k(:, i) - xm_k)';
    
end

% Add the process noise
Pm_k = real(Pm_k + Qv);




