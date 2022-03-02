%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MABBO.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [x,f,info] = MABBO(fun,x,st,tune);
%
% solve unconstrained noisy black box optimization problem 
%    min f(x) 
%  
% fun      % function handle for f(.)
% x        % starting point (must be specified)
% st       % structure with stop and print criteria
%          % (indefinite run if no stopping criterion is given)
%  .secmax       %   stop if sec>=secmax (default: inf)
%  .nfmax        %   stop if nf>=nfmax   (default: inf)
%  .qf           %   stop if qf<=accf    (default: 1e-4)
%  .fbest        %   optimal function value (default: 0)
%  .prt          %   printlevel (default: -1)
%                %   -2: nothing, -1: litte, >=0: more and more
% tune     % optional structure containing tuning parameters
%          %   for details see below
%
% x        % best point found 
% f        % function value at best point found 
% info     % structure with performance info
%          %   for details see below
% 
function [x, f,info] = MABBO(fun,x,st,tune)

info.error='';
% check function handle
if isempty(fun)
  message = 'subTR needs the function handle fun to be defined';
  info.error= message; 
  return
elseif ~isa(fun,'function_handle')
  message = 'fun should be a function handle';
  info.error= message;
  return
end
% starting point
if isempty(x)
  message = 'starting point must be defined';
  info.error= message;
  return      
elseif ~isa(x,'numeric')
  message = 'x should be a numeric vector'; 
  info.error= message;
  return       
end

n=length(x);

if ~exist('tune'), tune=[]; end;
 
if ~isfield(tune,'alg'), tune.alg = 0; end
 
if ~isfield(tune,'mem'), tune.mem = 10; end

if ~isfield(tune,'lambda'), tune.lambda = 4 + floor(3*log(n)); end
  
    
lambda = tune.lambda;

if ~isfield(tune,'mu'), tune.mu = floor(tune.lambda/2); end

mu=tune.mu;

if ~isfield(tune,'gamma'), tune.gamma = tune.lambda; end

gamma=tune.gamma;


if ~isfield(tune,'gammaE'), tune.gammaE = 3; end

if ~isfield(tune,'E'), tune.E = inf; end


if ~isfield(tune,'sigma_init'), tune.sigma_init = 1; end

sigma_init=tune.sigma_init;

if ~isfield(tune,'zeta'), tune.zeta = 1e-12; end  

if ~isfield(tune,'sigmamin'), tune.sigmamin = 1e-12; end

if ~isfield(tune,'sigmabar'), tune.sigmabar = 1e10; end

if ~isfield(tune,'sigmamax'), tune.sigmamax = 1e4; end

if ~isfield(tune,'sigmalbar'), tune.sigmalbar = 0.99; end

if ~isfield(tune,'q'), tune.q = 5; end
if ~isfield(tune,'epsa'), tune.epsa = 1e-2; end
if ~isfield(tune,'epsb'), tune.epsb = 0.85; end


if ~exist('st'), st=[]; end;
if isfield(st,'prt'), info.prt = st.prt; 
else,  info.prt = -1; 
end;
% stopping criteria
if isfield(st,'secmax'), info.secmax=st.secmax;
else, info.secmax=2000*n;
end;
if isfield(st,'nfmax'), info.nfmax=st.nfmax;
else, info.nfmax=2000*n;
end;
if isfield(st,'fbest'), info.fbest=st.fbest;
else, info.fbest=-inf; 
end;
if isfield(st,'accf'), info.accf=st.accf;
else, info.accf=-inf;
end;
if isfield(st,'finit'), info.finit=st.finit;
else, info.finit=fun(x);
end;
if isfield(st,'reallife'), info.reallife=st.reallife;
else, info.reallife=0;
end;
info.initTime=cputime;  info.nf=0;
dim      = length(x);
wi_raw   = log(lambda/2 + 0.5) - log((1:mu));
wi       = wi_raw/sum(wi_raw);
mu_eff   = 1/sum(wi .^2);
ordering = 'ascend'; sigma  = sigma_init; g = 0; 
xbest = x; fbest = fun(xbest);
c_s    = min(1.9999,(mu_eff + 2)/(dim + mu_eff + 5));
sqrt_s   = sqrt(c_s*(2-c_s)*mu_eff);
d_s      = 1 + c_s + 2*max( 0, sqrt((mu_eff-1)/(dim+1)) - 1 );
echi     = sqrt(dim)*(1 - 1/dim/4 - 1/dim/dim/21);
if ~tune.alg
    c_1  = 2/((dim+1.3)^2 + mu_eff);
    c_mu = min(1-c_1,2*(mu_eff-2 + 1/mu_eff)/((dim + 2)^2 + mu_eff));
    cp = []; cd = []; sqrt_cp=[];
else
    cp       =  min(1.9999,lambda/dim*4.^-(0:gamma-1));
    cd       = 1.5.^-(0:gamma-1) / dim;  
    sqrt_cp  = sqrt(cp .* (2-cp) * mu_eff);
    c_1=[]; c_mu=[];
end
s        = zeros(dim, 1);
if ~tune.alg, M = eye(dim);
else, M = zeros(dim,tune.gamma);  
end
ext=1; xss=[]; fss=[]; Parent=[];
while(1)
     info.g=g;
     if g>1, Parent.d0=Parent.d; end
     [OffspringPop,ranks,Parent,fs,info]= ...
          getMS(fun,Parent,cd,M,dim,xbest,wi,sigma,ordering,tune,info);
     if info.done, break; end
     [s,M] = updateInfo(OffspringPop,Parent,wi,M,s,ranks,sqrt_s,...
                                           sqrt_cp,cp,c_s,c_1,c_mu,tune);                                          
     [sigma,Parent] = getStepSize(xbest,ext,s,c_s,d_s,sigma,echi,...
                                                            Parent,tune);
     [Parent]=subspaceDir(g,Parent,tune);   
     Parent.y = xbest + sigma*Parent.d  ;
     [Parent.f]= fun(Parent.y);
     info.nf=info.nf+1;
     % check stopping test
     sec       = (cputime-info.initTime);
     info.done = (sec>info.secmax)|(info.nf>=info.nfmax);
     %info.qf   = abs((ftrue-info.fbest)/(info.finit-info.fbest));
     %info.done = (info.done|info.qf<=info.accf);
     info.sec  = sec;
     if info.done, break; end
     [ext,xbest,fbest,info]=extStepTri(fun,fs,Parent,...
                         sigma,xbest,fbest,tune,info);
     if info.done, break; end
     [xbest,fbest,xss,fss,info]= ...
         heuristicPoint(fun,ext,xss,fss,Parent,xbest,fbest,info,tune);
     g = g+1;          
end % of MABBO

x=xbest; f=fbest;

% status detection 

% if info.qf<=info.accf
%     info.status = 'accuracy reached';
% elseif info.nf>=info.nfmax
%      info.status = 'nfmax reached';
% elseif info.sec>=info.secmax
%     info.status = 'secmax reached';
% else
%     info.status = 'unknown';
% end


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

