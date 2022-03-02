%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% funf.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute function value at x
%
% input
%   x:     current point
%   noise: data struture for noise
%
% output
%   ftrue  true function value
%   f      noisy function value
function [f,ftrue] = funf(x,fun,noise)
xs = shiftPoint(length(x));
x = x+xs; % a shifted point is done
f=fun(x);
disp(['true function value ',num2str(f)])       
% add noise to f
ftrue=f;
if noise.noisefun
     switch noise.distr
     case 1, noise.epsilon = (2*rand-1)*noise.level;
     case 2, noise.epsilon = randn*noise.level;
         otherwise
        disp('error: "distr" should be either "uniform" or "Gauss" ')
     end
     switch noise.type
     case 1, f=f+noise.epsilon; 
     case 2, f=(1+noise.epsilon)*f;
     otherwise 
          disp('error: noisy.type should be either ab or rel')
     end
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%