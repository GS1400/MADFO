%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% extStepTri %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% extStepTri tries to do extrapolation along +d or -d
function [ext,xbest,fbest,info]=extStepTri(fun,fs,...
                                Parent,sigma,xbest,fbest,tune,info)        
ext=0;
fnm=stochasticNM(fs,fbest,Parent.f,tune);
if fnm>Parent.f+tune.gamma*sigma^2
    ext=1;
    [xbest,fbest,info]=extStepDone(fun,fs,...
                               Parent,sigma,xbest,fbest,tune,info);                             
    if info.done, return; end
end
if ext==0
   Parent.d = -Parent.d;
   Parent.y = xbest +sigma*Parent.d;
   [Parent.f,info.ftrue]= fun(Parent.y);
   info.nf=info.nf+1;
   % check stopping test
   sec       = (cputime-info.initTime);
   info.done = (sec>info.secmax)|(info.nf>=info.nfmax);
   info.qf   = abs((info.ftrue-info.fbest)/(info.finit-info.fbest));
   info.done = (info.done|info.qf<=info.accf);
   info.sec  = sec;
   fnm=stochasticNM(fs,fbest,Parent.f,tune);
   if fnm>Parent.f+tune.gamma*sigma^2
      ext=1;
      [xbest,fbest,info]=extStepDone(fun,fs,Parent,sigma,xbest,...
                                                    fbest,tune,info);
      if info.done, return; end
   end
   if info.done, return; end
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

