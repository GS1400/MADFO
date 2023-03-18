

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% initTune.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function tune=initTune(tune,n);
% initialize tuning parameters
%
% tune         % on input, tune contains those tuning parameters
%              % that are explicitly set; the others take default values
%
% n            % problem dimension

function tune=initTune(tune,n)

% tune % structure containing all tuning parameters 
%      % all parameters have a default that can be overwritten 
%      % by specifying it as an input
if isempty('tune')||nargin==1, tune=[]; end
if nargin==0, error('initTune needs the dimension n as input'); end
% type of algorithm: 
%      0: defaults
%      1: no good step sizes in the mutation and recombination phases
%      2: no recombination subspace direction
%      3: no nonmonotone line search
%      4: no heuristic point
if ~isfield(tune,'ver'), tune.ver = 0; end
%  sample size 
if ~isfield(tune,'lambda'), tune.lambda = 4 + floor(3*log(n)); end
% number of selected search points in the populations
if ~isfield(tune,'mu'), tune.mu = floor(tune.lambda/2); end
% parameters for expanding step size
if ~isfield(tune,'gammaE'), tune.gammaE = 4; end
% maximum iterations of extrapolation step
if ~isfield(tune,'E'), tune.E = inf; end
% initial step sizee
if ~isfield(tune,'sigmainit'), tune.sigmainit = 1; end
% parameter for line search condition
if ~isfield(tune,'gamma'), tune.gamma = 1e-12; end  
% minimum threshold for sigma
if ~isfield(tune,'sigmamin'), tune.sigmamin = 1e-2; end
% maximum value for sigma
if ~isfield(tune,'sigmamax'), tune.sigmamax = 0.5; end
% factor for the memory of the non-monotone term 
if ~isfield(tune,'kappa'), tune.kappa = 5; end
% memory for the list of function values at best points
% to choose three vertices of the first triangle
if ~isfield(tune,'N'), tune.N = 10; end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 