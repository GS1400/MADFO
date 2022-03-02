%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MADFO.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [x,f,info] = MADFO(fun,x,st,tune);
%
% solve unconstrained noisy derivative free optimization problem 
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
%                %   -1: nothing, 0: litte, >=0: more and more
% tune     % optional structure containing tuning parameters
%          %   for details see below
%
% x        % best point found 
% f        % function value at best point found 
% info     % structure with performance info
%          %   for details see below
% 
function [x, f,info] = MADFO(fun,x,st,tune)
info.error='';
% check function handle
if isempty(fun)
  message = 'MADFO needs the function handle fun to be defined';
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
% tune % structure containing all tuning parameters 
%      % all parameters have a default that can be overwritten 
%      % by specifying it as an input
if ~exist('tune'), tune=[]; end
% type of algorithm: 0 (fast version), 1 (limited version)
if ~isfield(tune,'ver'), tune.ver = 0; end
%  memory for non-monotone term 
if ~isfield(tune,'mem'), tune.mem = 10; end
%  sample size 
if ~isfield(tune,'lambda'), tune.lambda = 4 + floor(3*log(n)); end
lambda = tune.lambda;
% number of selected search points in the populations
if ~isfield(tune,'mu'), tune.mu = floor(tune.lambda/2); end
mu=tune.mu;
% memory for evaluation matrx
if ~isfield(tune,'mmax'), tune.mmax = tune.lambda; end
mmax=tune.mmax;
% parameters for expanding step size
if ~isfield(tune,'gammaE'), tune.gammaE = 2; end
% maximum iterations of extrapolation step
if ~isfield(tune,'E'), tune.E = inf; end
% initial step sizee
if ~isfield(tune,'sigma_init'), tune.sigmainit = 1; end
sigmainit=tune.sigmainit;
% parameter for line search condition
if ~isfield(tune,'gamma'), tune.gamma = 1e-12; end  
% minimum threshold for sigma
if ~isfield(tune,'sigmamin'), tune.sigmamin = 1e-12; end
% uppper bound on the heuristic step size
if ~isfield(tune,'sigmabar'), tune.sigmabar = 1e10; end
% maximum value for sigma
if ~isfield(tune,'sigmamax'), tune.sigmamax = 1e4; end
%  lower bound on the heuristic step size
if ~isfield(tune,'sigmalbar'), tune.sigmalbar = 0.99; end
% parameter for updating heuristic step size in mutation
if ~isfield(tune,'q'), tune.q = 5; end
% parameters for updating step sizes for subspace directions 
if ~isfield(tune,'epsa'), tune.epsa = 1e-2; end
if ~isfield(tune,'epsb'), tune.epsb = 0.85; end
if ~exist('st'), st=[]; end
if isfield(st,'prt'), info.prt = st.prt; 
else,  info.prt = -1; 
end
% stopping criteria
if isfield(st,'secmax'), info.secmax=st.secmax;
else, info.secmax=180;
end
if isfield(st,'nfmax'), info.nfmax=st.nfmax;
else, info.nfmax=500*n;
end
if isfield(st,'fbest'), info.fbest=st.fbest;
else, info.fbest=0; 
end
if isfield(st,'accf'), info.accf=st.accf;
else, info.accf=1e-4;
end
if isfield(st,'finit'), info.finit=st.finit;
else, info.finit=fun(x);
end
info.initTime=cputime;  info.nf=0;
dim      = length(x);
wi_raw   = log(lambda/2 + 0.5) - log((1:mu));
wi       = wi_raw/sum(wi_raw);
mu_eff   = 1/sum(wi .^2);
ordering = 'ascend'; sigma  = sigmainit; g = 0; 
xbest = x; fbest = fun(xbest);
c_s    = min(1.9999,(mu_eff + 2)/(dim + mu_eff + 5));
sqrt_s   = sqrt(c_s*(2-c_s)*mu_eff);
d_s      = 1 + c_s + 2*max( 0, sqrt((mu_eff-1)/(dim+1)) - 1 );
echi     = sqrt(dim)*(1 - 1/dim/4 - 1/dim/dim/21);
if ~tune.ver
    c_1  = 2/((dim+1.3)^2 + mu_eff);
    c_mu = min(1-c_1,2*(mu_eff-2 + 1/mu_eff)/((dim + 2)^2 + mu_eff));
    cm = []; cd = []; sqrt_cm=[];
else
    cm       =  min(1.9999,lambda/dim*4.^-(0:mmax-1));
    cd       = 1.5.^-(0:mmax-1) / dim;  
    sqrt_cm  = sqrt(cm .* (2-cm) * mu_eff);
    c_1=[]; c_mu=[];
end
s        = zeros(dim, 1);
if ~tune.ver, M = eye(dim);
else, M = zeros(dim,tune.mmax);  
end
ext=1; xss=[]; fss=[]; Parent=[];
while(1)
     info.g=g;
     if info.prt>=0
        disp('=======================================================')
        disp('=======================================================')
        disp('=======================================================')
        disp(['MADFO at iter = ', num2str(g)])
        disp(['number of func. eval. used so far = ',num2str(info.nf)])
        disp('=======================================================')
     end
     if g>1, Parent.d0=Parent.d; end
     if info.prt>=0
        disp(' ')
        disp('start of getMS')
        disp(' ')
     end
     [OffspringPop,ranks,Parent,fs,info]= ...
          getMS(fun,Parent,cd,M,dim,xbest,wi,sigma,ordering,tune,info);
      if info.prt>=0
        disp(' ') 
        disp('end of getMS')
     end 
     if info.done, break; end
     [s,M] = updateInfo(OffspringPop,Parent,wi,M,s,ranks,sqrt_s,...
                                       sqrt_cm,cm,c_s,c_1,c_mu,tune);
     if info.prt>=0
        disp(' ')
        disp('updateInfo was done')
     end
     [sigma,Parent] = getStepSize(xbest,ext,s,c_s,d_s,sigma,echi,...
                                                        Parent,tune);
     if info.prt>=0
        disp(' ') 
        disp('getStepSize was done')
        disp(' ')
     end
     [Parent]=subspaceDir(g,Parent,tune);   
      if info.prt>=0
        disp('subspaceDir was done')
     end
     if info.prt>=0
        disp(' ') 
        disp('start of extStepTri')
        disp(' ')
     end
     Parent.y = xbest + sigma*Parent.d  ;
     [Parent.f,info.ftrue]= fun(Parent.y);
     info.nf=info.nf+1;
     % check stopping test
     sec       = (cputime-info.initTime);
     info.done = (sec>info.secmax)|(info.nf>=info.nfmax);
     info.qf   = abs((info.ftrue-info.fbest)/(info.finit-info.fbest));
     info.done = (info.done|info.qf<=info.accf);
     info.sec  = sec;
     if info.done, break; end
     [ext,xbest,fbest,info]=extStepTri(fun,fs,Parent,...
                         sigma,xbest,fbest,tune,info);
     if info.prt>=0
        if info.prt>=1&&~ext
            disp('extrapolation could not be performed')
        elseif info.prt>=1
           disp('extrapolation has been performed successfully')
        end
        disp(' ')
        disp('end of extStepTri')
     end                
     if info.done, break; end
     if info.prt>=0
        disp(' ')
        disp('start of heuristicPoint')
        disp(' ')
     end
     [xbest,fbest,xss,fss,info]= ...
         heuristicPoint(fun,ext,xss,fss,Parent,xbest,fbest,info,tune);
     if info.prt>=0
        disp(' ')
        disp('end of heuristicPoint')
        disp('=======================================================')
        disp('=======================================================')
        disp('=======================================================')
     end
     g = g+1;          
end % of MADFO
x=xbest; f=fbest; 
% status detection 
if info.qf<=info.accf, info.status = 'accuracy reached';
elseif info.nf>=info.nfmax, info.status = 'nfmax reached';
elseif info.sec>=info.secmax,  info.status = 'secmax reached';
else, info.status = 'unknown';
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

