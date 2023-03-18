%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% randNM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% randNM computes randomized non-monotone term
%
function R=randNM(fs,fbest,fnew,tune)
n = length(fs);
if n>tune.lambda, I=randperm(n,tune.lambda);
else, I=1:n;
end
fs=(fs(I)); fmin =min(fbest,min(fs)); fmax = max(fs); fmed = median(fs);
eta1 = (fmed-fmin)/(fmax-fmin);
eta2 = (fmax-fmed)/(fmax-fmin);
if eta1==0 && eta2==0, eta = rand;
elseif eta1==0, eta = eta2;
elseif eta2==0, eta = eta1;
else, eta = min(eta1,eta2);
end
if isnan(eta), eta = rand; end
eta = max(0.5+rand,eta);
if fnew>=fmax, R = eta*fmed+(1-eta)*fmax;
else, R = (1-eta)*fmed+eta*fmin;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

