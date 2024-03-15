 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% goodStepSize %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% goodStepSize compute mutation and recombination step sizes
%
function sigma=goodStepSize(Vstate,xbest,Mutation,RecomInfo,...
                                              sigma0,sigma,tune,nf,prt)
switch Vstate
    case 3 % in mutation phase
         alp = abs(xbest)./abs(Mutation.d);
         II  = find(isnan(alp)|alp==0|alp>=tune.sigmamax|...
                                               alp<=tune.sigmamin);
         alp(II) = []; 
         if ~isempty(alp)
             % good mutation step size is computed
             sigma= max(tune.sigmamin,min(tune.sigmamax,sigma0));
             sigma = sqrt(sigma*median(alp));
         end 
    case 6 % in recombination phase
            alp = abs(xbest)./abs(RecomInfo.d);
            II  = find(isnan(alp)|alp==0|alp>=tune.sigmamax|...
                                              alp<=tune.sigmamin); 
            alp(II) = [];
            if ~isempty(alp)
               % good recombination step size is computed
               sigma = max(tune.sigmamin,min(tune.sigmamax,sigma0));
               sigma = sqrt(sigma*median(alp));
               if prt>=1
                 disp(' ')
                 disp(['good recombination step size',...
                       ' sigma = ',num2str(sigma),...
                       ' at nf=',num2str(nf)])
               end 
            end
    otherwise
        error('Vstate should be 3 or 6' )
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%