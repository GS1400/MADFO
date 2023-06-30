
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% driver.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% driver for MADFO
%
clear
clc
rng('default')
% standard initialization
cas=2;  
switch cas
    case 1, st.m=2; st.n=3;   % problem size, [m,n]=size(x)
           x=[-1 -1 -1; -2,-2,-2];  % starting point 
           xx = [x(1,2);x(2,1)];
    case 2, st.m=1; st.n=2; x=[-1 -1]; xx = x;
end                      
% function handle
fun=@(xx)(xx(1)-1)^2+100*(xx(2)-xx(1)^2)^2; % Rosenbrock function
tune=[];
st.prt=0;      %   printlevel (default: -1)
%              %   -1: nothing, 0: litte, >=1: more and more
st.nfmax=5000; %   stop if nf>=nfmax   (default: inf)
st.secmax=inf; %   stop if sec>=secmax (default: inf)
st.ftarget=0;  % target function value (to compute qf)
               %  qf=(ftarget-f)/(finit-f)
st.accf=1e-10; % target accuracy for stopping test qf<=accf
               % (default: 1e-4)
% MADFO repeatedly calls MADFOstep
[xbest,fbest,info] = MADFO(fun,x,st,tune)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
