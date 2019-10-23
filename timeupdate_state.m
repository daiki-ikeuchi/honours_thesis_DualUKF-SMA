% A priori estimate for state variables at k time 
function xm_k = timeupdate_state(Wm,Xm_k,Num_StateVar,Num_sigma)

% Wm = weights for mean state (5 by 1 matrix(
% Xm_k = sigma points at k time before measurement update (2 by 5 matrix)
% Num_StateVar = 2
% Num_sigma = 5;

% Initialize
xm_k = zeros(Num_StateVar, 1);

% Calculate a priori estimates of state variables (2 by 1)
% where (1,1) = tempearure and (2,1) = stress
for i = 1:Num_sigma    
   
    xm_k = xm_k + Wm(i,1)*Xm_k(:, i);    
end  

xm_k = real(xm_k);