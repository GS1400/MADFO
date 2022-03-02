%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% subspaceDir %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% subspaceDir generates subspace direction
function [Parent]=subspaceDir(g,Parent,tune)
if g>1
    if any(Parent.d~=0), alp = abs(Parent.d)./abs(Parent.d0);
    else, alp = 1./abs(Parent.d0); 
    end
    II = find(isnan(alp)|alp>=tune.sigmabar); alp(II) = [];
    if ~isempty(alp) 
       beta = (tune.epsa/(1+g)^tune.epsb)*rand*max(alp);
        if beta<1,  Parent.d = Parent.d + beta*Parent.d0; end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

