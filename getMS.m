%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% getMS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% getMS do mutation and selection
function [OffspringPop,ranks,Parent,fs,info]=getMS(fun,Parent,cd,M,...
                                  dim,xbest,wi,sigma,ordering,tune,info)
ranks=[]; fs=[];
for ll=1:tune.lambda
      Offspring.z = randn(dim, 1);
      if ~tune.ver, Offspring.d = M*Offspring.z;
      else
          d = Offspring.z;
          for j = 1:min(info.g,tune.mmax)
             d = (1-cd(j))*d + (cd(j)*(M(:, j)'*d)) * M(:, j);
          end
          Offspring.d = d;  
      end
     if any(xbest~=0), alp = abs(xbest)./abs(Offspring.d);
     else,  alp = 1./abs(Offspring.d);
     end
     II  = find(isnan(alp)|alp<2*sigma); alp(II) = []; 
     if ~isempty(alp) && any(alp~=0)
         alpmax = max(sigma,(min(alp)*sigma)^(1/tune.q));
     else, alpmax = sigma;
     end
      xs  = xbest + alpmax*Offspring.d;
      [Offspring.f,info.ftrue] = fun(xs);
      OffspringPop(ll) = Offspring;
      info.nf=info.nf+1;
      % check stopping test
      sec       = (cputime-info.initTime);
      info.done = (sec>info.secmax)|(info.nf>=info.nfmax);
      info.qf   = abs((info.ftrue-info.fbest)/(info.finit-info.fbest));
      info.done = (info.done|info.qf<=info.accf);
      info.sec  = sec;
      if info.done, break; end
end
if ~info.done
    [ranks,fs] = RankPop(OffspringPop, ordering);
    sum_z = zeros(dim, 1);
    sum_d = zeros(dim, 1);
    for m = 1:tune.mu
      sum_z = sum_z + wi(m)*OffspringPop(ranks(m)).z;
      sum_d = sum_d + wi(m)*OffspringPop(ranks(m)).d;
    end
    Parent.z = sum_z;
    Parent.d = sum_d;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

