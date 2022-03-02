%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% stochasticNM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% stochasticNM computes stochastic non-monotone term
function R=stochasticNM(fs,fbest,fnew,tune)
n = length(fs);
if n>tune.mem, I=randperm(n,tune.mem);
else, I=1:n;
end
fs=(fs(I)); fbest =min(fbest,min(fs)); fmax = max(fs); fmed = median(fs);
eta1 = (fmed-fbest)/(fmax-fbest);
eta2 = (fmax-fmed)/(fmax-fbest);
if eta1==0 && eta2==0, eta = rand;
elseif eta1==0, eta = eta2;
elseif eta2==0, eta = eta1;
else, eta = min(eta1,eta2);
end
if isnan(eta), eta = rand; end
eta = eta/(2+rand);
if fnew >= fmax, R  = (1-eta)*fmax+eta*fmed;
 elseif fnew>=fmed, R = (1-eta)*fmed+eta*fmax;
elseif fnew>=fbest, R = (1-eta)*fmed+eta*fbest;
else, R = (1-eta)*fbest+eta*fmed;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

