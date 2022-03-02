%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% getStepSize %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% getStepSize updates step size sigma
function [sigma,Parent]=getStepSize(xbest,ext,s,c_s,d_s,sigma,...
                                                      echi,Parent,tune)
tt = (c_s/d_s)*((norm(s)/echi)-1);
if tt>0 && ~ext,tt=-tt; end
sigma = min(tune.sigmamax,sigma*exp(tt));
if sigma <=tune.sigmamin
    if any(xbest~=0), alp = abs(xbest)./abs(Parent.d);
    else, alp   = 1./abs(Parent.d);  
    end
    II = find(isnan(alp)|alp>=tune.sigmabar); alp(II) = [];
   if ~isempty(alp), 
       sigma=max(tune.sigmalbar,min(tune.sigmamax,exp(tt)*max(alp))); 
   end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




