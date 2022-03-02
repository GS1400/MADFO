%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% extStepDone %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% extStepDone performs an extrapolation along d
function [xbest,fbest,info]=extStepDone(fun,fs,...
                                  Parent,sigma,xbest,fbest,tune,info)
F=Parent.f; alpha=sigma; ii=0; 
while 1
    sigma = tune.gammaE*sigma; alpha=[alpha,sigma];
    Parent.y = xbest + sigma*Parent.d;
    [Parent.f,ftrue] = fun(Parent.y); F=[F,Parent.f]; 
    info.nf=info.nf+1;
    % check stopping test
    sec       = (cputime-info.initTime);
    info.done = (sec>info.secmax)|(info.nf>=info.nfmax);
    info.qf   = abs((ftrue-info.fbest)/(info.finit-info.fbest));
    info.done = (info.done|info.qf<=info.accf);
    info.sec  = sec;
    ii        = ii+1;
    fnm       = stochasticNM([F,fs],fbest,Parent.f,tune);
    if fnm<Parent.f+tune.gamma*sigma^2 || ii==tune.E
        [~,ii] = min(F); ii=ii(1); sigma=alpha(ii); fbest=F(ii);
        xbest=xbest+sigma*Parent.d;   
        break;
    end
    if info.done, break; end
end    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


