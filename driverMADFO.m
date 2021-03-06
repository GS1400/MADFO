%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% driverVRBBON.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:3
    for j=1:100
        fprintf('=')
    end
    fprintf('\n')
end
fprintf(['MADFO solves a unconstrained noisy derivative free',...
    ' optimization of a not necessarily smooth function of \n',...
    ' many continuous arguments. No derivatives are needed.',...
    ' A limited amount of noise is tolerated.\n\n']);
disp('===============================================================')
clear;
solverPath = input(['Insert the MADFO path \n',...
    '>> solverPath='],'s');
if ~exist(solverPath, 'dir')
    disp('the directory does not exist')
    return
end
eval(['addpath ',solverPath,'/MADFO'])
disp('===============================================================')
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create noise
fprintf(['noise.level: 0.0001/0.01/0.1 \n',...
         'noise.type:  1 (absolute) or  2 (relative)\n',...
         'noise.distr: 1 (uniform)  or  2 (Gauss)\n'])
noise = struct('noisefun',1,'level',0.01,'type',2,'distr',2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create problem definition

% define problem parameters (to be adapted to your problem)
n=10; % dimension
p=2;  % Norm in objective function
e=1;  % Exponent in objective function function

% create random matrix and right hand side for objective function
% (specific to the model problem; replace by whatever data you
% need to provide to your objectiv function)
A=rand(n)-0.5; 
b=-sum(A,2);

% create objective function f(x)=||Ax-b||_p^e
fun0  = @(x) norm(A*x-b,p).^e; 
fun  = @(x) funf(x,fun0,noise);  

% To solve your own problem, simply replace in the previous line 
% the expression after @(x) by your expression or function call. 
% Parameters in this expression are passed by value, hence are 
% fixed during minimization.

% start and stop info
x      = 2*rand(n,1);  % starting point

% problem definition complete
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('===============================================================')
fullinfo=1;  Tuning=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% solve the problem with MADFO
tic; % set clock
if fullinfo    % define stop and print criteria
               % (indefinite run if no stopping criterion is given)
    % pass stop and print criteria
    % (indefinite run if no stopping criterion is given)
    st = struct('secmax',180,'nfmax',10000*n,'finit',fun(x),...
         'fbest',0.001*fun(x),'accf',0.001,'prt',1)
else
    st = []; % budgets are chosen inside MADFO
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  solve the problem with MADFO
if Tuning % self-tuning and info
    % given are the defaults
    % only the deviating parameters need to be set!
     fprintf([...
        '====================================================\n',...
        'type of algorithm: 0(fast version), 1(limited version): ',...
        'tune.ver = 0 \n',...
        '====================================================\n',...
        'memory for non-monotone term : ',...
        'tune.mem = 10;  \n',...
        '====================================================\n',...
        'maximal number of extrapolation step: ',...
        'tune.E= Inf \n',...
        '====================================================\n',...
        '%  sample size: ',...
        'tune.lambda = 4 + floor(3*log(n)); \n',...
        '====================================================\n',...
        'number of selected search points in the populations:  \n',...
        'tune.mu = floor(tune.lambda/2); \n',...
        '====================================================\n',...
        'memory for evaluation matrx: ',...
        'tune.mmax = tune.lambda; \n',...
        '====================================================\n',...
        'factor for extrapolation test: tune.gammaE = 2 \n',...
        '====================================================\n',...
        'initial step size: tune.sigmainit = 1; \n',...
        '====================================================\n',...
        'minimum threshold for sigma: ',...
        'tune.sigmamin = 1e-12; \n',...
        '====================================================\n',...
        'uppper bound on the heuristic step size: ',...
        'tune.sigmabar = 1e10',...
        '====================================================\n',...
        'maximum value for sigma: ',...
        'tune.sigmamax = 1e4; \n',...
        '====================================================\n',...
        'lower bound on the heuristic step size: ',...
        'tune.sigmalbar = 0.99; \n',...
        '====================================================\n',...
        'parameter for updating heuristic step size in mutation: ',...
        'tune.q = 5;  \n',...
        '====================================================\n',...
        'parameter for non-monotone line search condition: ',...
        'tune.gamma = 1e-12;  \n',...
        '====================================================\n',...
        'parameters for updating step sizes for subspace',...
        ' directions: ',...
        ' tune.epsa = 1e-2; \n',...
        ' tune.epsb = 0.85; \n',...
        '====================================================\n']);     
        tune = struct('ver',0,'mem',10,'E',inf,...
            'lambda',4 + floor(3*log(n)),...
            'mu',2+floor(1.5*log(n)),...
            'mmax',4 + floor(3*log(n)),'gammaE',2,...
            'sigmainit',1,'sigmamin',1e-12,'sigmabar',1e10,...
            'sigmamax',1e4,'sigmalbar',0.99,'q',5,...
            'epsa',1e-2,'epsb',0.85,'gamma',1e-12);
else
   tune = []; % full tuning inside MADFO is used
end
[x,f,info] = MADFO(fun,x,st,tune);
% the problem solved possibly by MADFO
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(' '), disp(' '), disp(' '), disp(' ')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
% display output
if ~isempty(info.error), error = info.error
else
if info.prt>=-1
disp('==============================================================');
disp('MADFO completed silently');
disp('==============================================================');
disp('display output:');

format long
info               % progress report by MADFO
disp('------------------------------------------------------------');
nfused=info.nf;               % number of function evaluations used
disp(['the number of function evaluations used: ',num2str(nfused)])
disp('------------------------------------------------------------');
secused=cputime-info.initTime;   % time used up to now
disp(['time used up to now: ',num2str(secused)])
disp('-------------------------------------------------------------');
if noise.noisefun
  ftrue = info.ftrue;
  disp(['the noisy function value at xbest (f): ',num2str(f)])
else
  disp(['the function value at xbest (f): ',num2str(f)])
end
disp('-------------------------------------------------------------');
qf = info.qf;
disp(['qf:=(ftrue-fbest)/(finit-fbest): ',num2str(qf)])
disp('-------------------------------------------------------------');
end
end
for i=1:3
    for j=1:100
        fprintf('=')
    end
    fprintf('\n')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

