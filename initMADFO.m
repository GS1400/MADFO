%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% initTune.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [wi,sigma,iter,d_s,echi,c_s,sqrt_s,nhp,RecomInfo,Fnm,... 
%                      next, c_1,c_mu,s,M,Xt,Ft]=initMADFO(dim,tune)
%
% initialize M and sigma and related parameters used for updating them
%
% tune         % on input, tune contains those tuning parameters
%              % that are explicitly set; the others take default values
% dim          % problem dimension

function  [wi,sigma,iter,d_s,echi,c_s,sqrt_s,nhp,RecomInfo,Fnm,... 
                     next, c_1,c_mu,s,M,Xt,Ft]=initMADFO(dim,tune)

wi_raw   = log(tune.lambda/2 + 0.5) - log((1:tune.mu));
wi       = wi_raw/sum(wi_raw);
mu_eff   = 1/sum(wi .^2); sigma  = tune.sigmainit; iter = 0; 
c_s    = min(1.9999,(mu_eff + 2)/(dim + mu_eff + 5));
sqrt_s   = sqrt(c_s*(2-c_s)*mu_eff);
d_s      = 1 + c_s + 2*max( 0, sqrt((mu_eff-1)/(dim+1)) - 1 );
echi     = sqrt(dim)*(1 - 1/dim/4 - 1/dim/dim/21);
c_1  = 2/((dim+1.3)^2 + mu_eff);
c_mu = min(1-c_1,2*(mu_eff-2 + 1/mu_eff)/((dim + 2)^2 + mu_eff));
s        = zeros(dim, 1);
M = eye(dim); Xt=[]; Ft=[]; RecomInfo=[];  Fnm=[];
nhp=0; next=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%