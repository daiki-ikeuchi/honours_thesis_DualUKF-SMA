% A priori measurement estimate at k time
function ym_k = timeupdate_obs_param(Wm_q,Ym_k,Num_ObsVar,Num_sigma)

% Wm = weights for states (5 by 1)
% Ym_k = unscented transformed sigma points using observation function (1
% by 5)
% Num_ObsVar = number of observed variable (just 1 for resistance)
% Num_sigma = number of sigma points (5 in this case)


% Initialize
ym_k = zeros(Num_ObsVar, 1);

% Calculate a priori estimates of measurement (1 by 1)
% where (1,1) = resistance
for i = 1:Num_sigma    
   
    ym_k = ym_k + Wm_q(i,1)*Ym_k(:, i);    
end  

ym_k = real(ym_k);