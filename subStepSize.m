%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% subStepSize %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% subStepSize generates step sizes greater than one
function beta = subStepSize(x,d,g,tune)
if any(x~=0)
    alp     = abs(x)./abs(d);
    II      = find(isnan(alp)|alp>=tune.sigmabar);
    alp(II) = [];
 else
    alp     = 1./abs(d);
    II      = find(isnan(alp)|alp>=tune.sigmabar);
    alp(II) = [];    
end
if ~isempty(alp) && any(alp~=0)
     beta = max(1+rand,rand*(tune.epsa/(1+g)^tune.epsb)*max(alp));
else, beta = 1+rand;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

