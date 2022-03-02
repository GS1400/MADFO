%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% updateInfo %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% updateInfo updates s and M
function [s,M]=updateInfo(OffspringPop,Parent,wi,M,s,ranks,sqrt_s,...
                                          sqrt_cm,cm,c_s,c_1,c_mu,tune)
s = (1-c_s)*s + sqrt_s*Parent.z;
if ~tune.ver
    M = (1 - 0.5*c_1 - 0.5*c_mu) * M + (0.5*c_1)*(M*s)*s';
    for m = 1:tune.mu
      M = M + ((0.5*c_mu*wi(m))*OffspringPop(ranks(m)).d)*...
           OffspringPop(ranks(m)).z';
    end
else
     for i=1:tune.mmax
        M(:, i) = (1-cm(i))*M(:, i) + sqrt_cm(i)*Parent.z;
     end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

