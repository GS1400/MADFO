
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MADFO.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [x,f,info] = MADFO(fun,x,st,tune);
%
% solves the unconstrained noisy derivative free optimization problem 
%    min f(x) 
%  
% fun          % function handle for objective function
% x            % starting point (must be specified)
% st           % structure with stop and print criteria
%              %   (indefinite run if no stopping criterion is given)
%  .secmax     %   stop if sec>=secmax (default: inf)
%  .nfmax      %   stop if nf>=nfmax   (default: inf)
%  .qf         %   stop if qf<=acc    (default: 1e-4)
%  .ftarget    %   function value accepted as optimal (default: 0)
%  .prt        %   printlevel (default: -1)
%              %   -1: nothing, 0: litte, >=1: more and more
% tune         % optional structure specifying tuning parameters
%              %   for details see initTune.m
%
% x            % best point found 
% f            % function value at best point found 
% info         % performance information for MADFO
%  .finit      %   initial function value
%  .ftarget    %   target function value (to compute qf)
%  .qf         %   (ftarget-f)/(finit-f)
%  .initTime   %   inital cputime
%  .done       %   done with the search?
%  .acc        %   stop when qf<=acc
%  .secmax     %   stop if sec>=secmax 
%  .nfmax      %   stop if nf>=nfmax 
%  .finit      %   the initial f
%  .prt        %   printlevel 
%              %     -1: nothing, 0: litte, >=1: more and more
% 
function [xbest,fbest,info] = MADFO(fun,x,st,tune);
persistent m n dim 
%%%%%%%%%%%%%%%%%%%%%%%%
%%%% initial checks %%%%
%%%%%%%%%%%%%%%%%%%%%%%%
prt=st.prt;
if prt>=0,
  disp(' ')
  disp('==============================================================')
  disp('start MADFO')
  disp('==============================================================')
end;
% check function handle
if isempty(fun)
  message = 'MADFO needs the function handle fun to be defined';
  disp(message)  
  return
elseif ~isa(fun,'function_handle')
  message = 'fun should be a function handle';
  disp(message)
  return
end;
% starting point
if isempty(x)
  message = 'starting point must be defined';
  disp(message) 
  return      
elseif ~isa(x,'numeric')
  message = 'x should be a numeric vector'; 
  disp(message ) 
  return       
end
m=st.m; n=st.n; dim=m*n;
[mx,nx]=size(x);
if mx~=m && nx~=n,
  sizex=[mx,nx], sizeNeeded=[m,n] 
  error('dimension mismatch');
end
% add info on stopping criteria
if isfield(st,'secmax'), info.secmax=st.secmax;
else, info.secmax=inf;
end;
if isfield(st,'nfmax'), info.nfmax=st.nfmax;
else, info.nfmax=inf;
end;
if isfield(st,'ftarget'), info.ftarget=st.ftarget;
else, info.ftarget=0;
end;
if isfield(st,'accf'), info.accf=st.accf;
else, info.acc=-inf;
end;
info.prt = prt; % print level
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% initialize solver environment %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MADFOstep(tune,fun,prt,dim);
initTime=cputime;
%%%%%%%%%%%%%%%%%%%
%%%% main loop %%%%
%%%%%%%%%%%%%%%%%%%
nf=0;
while 1
  % get function value
  f  = fun(x); nf = nf+1;
  % get new point for function evaluation
  x         = MADFOstep(x(:),f,prt);
  x         = reshape(x,m,n);
  % restore original format of x
  sec       = (cputime-initTime);
  info.done = (sec>st.secmax)|(nf>=st.nfmax);
  if nf>1
     qf        = (f-st.farget)/(finit-st.farget);
     info.qf   = qf;
     info.done = (info.done|info.qf<=st.acc);
     if f<fb
         fb=f;
          if prt>=0
             disp(['improved function value f=',num2str(fb),...
               ' at ',num2str(nf)])
          end
     end
  else
      finit = f; fb=f;
      if prt>=0
         disp(['improved function value f=',num2str(fb),...
               ' at ',num2str(nf)])
       end
  end
  
  info.sec  = sec;
  if info.done, break; end
end
%%%%%%%%%%%%%%%%%%%%%
%%%% return info %%%%
%%%%%%%%%%%%%%%%%%%%%
[xbest,fbest,info]= MADFOstep;
% update info
info.qs_convergeSpeed  = qf;
info.time_In_Second    = sec;
%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% solution status  %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
if info.qs_convergeSpeed<=st.acc
  info.status_of_converge = 'accuracy reached';
elseif nf>=st.nfmax
   info.status_of_converge = 'nfmax reached';
elseif sec>=st.secmax
  info.status_of_converge = 'secmax reached';
else
  info.status_of_converge = 'unknown';
end;
if prt>=0,
  disp('==============================================================')
  disp('end MADFO')
  disp('==============================================================')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

