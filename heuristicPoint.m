%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% heuristicPoint %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% heuristicPoint generates at most five heuristic points, one of
% which is accepted
function  [xbest,fbest,xs,fs,info]=heuristicPoint(fun,ext,xs,fs,...
                                        Parent,xbest,fbest,info,tune)                                  
g = info.g;  
fnm = stochasticNM(fs,fbest,Parent.f,tune);                                    
if size(xs,2)<=3
  xs=[xs,Parent.y]; fs=[fs,Parent.f];
else
  if Parent.f < fnm
     [~,iw] = max(fs); xs(:,iw)=Parent.y; fs(iw)=Parent.f;
     [fs,Is] = sort(fs); xs=xs(:,Is); 
  else
      [~,iw] = max(fs); xs(:,iw) = xbest; fs(iw) = fbest;
      [fs,Is] = sort(fs); xs=xs(:,Is); 
  end
end
if info.prt>=1
   disp(' ')
   disp('vertices of triangle were updated')
end
if ext==0
     if Parent.f < fnm
        fbest= Parent.f; xbest = Parent.y; 
        if info.prt>=1
           disp(' ') 
           disp('heuristicPoint updates a better point')
        end
     else    
         if size(xs,2)>=3
             x1 = xs(:,1); x2 = xs(:,2); x3 = xs(:,3); 
             x23 = 0.5*(x2+x3); x13 = 0.5*(x1+x3); x12 = 0.5*(x1+x2);
             for i=1:5
                 switch i
                     case 1
                          d    = x1-x23;  d1 = d;
                          alp  = subStepSize(x23,d,g,tune);
                          xt   = x23+ alp*d; xx(:,1)=xt;
                     case 2
                           d = x12-x23; 
                           Parent.d=d; Parent.d0=d1;
                           Parent = subspaceDir(g,Parent,tune);
                           d      =  Parent.d;
                           alp    = subStepSize(x23,d,g,tune);
                           xt     = x2+ alp*d; xx(:,2)=xt;
                     case 3
                           d = x13-x23; 
                           Parent.d=d; Parent.d0=d1;
                           Parent = subspaceDir(g,Parent,tune);
                           d      = Parent.d;
                           alp    = subStepSize(x23,d,g,tune);
                           xt     = x23+ alp*d; xx(:,3)=xt;
                     case 4
                         xt = randsubPoint(x1,x13,x12); xx(:,4)=xt;
                     case 5
                         xt = randsubPoint(x23,x13,x12); xx(:,5)=xt;
                 end              
               [ft,info.ftrue] = fun(xt); ff(:,i) = ft;
               info.nf= info.nf+1;
               % check stopping test
              sec       = (cputime-info.initTime);
              info.done = (sec>info.secmax)|(info.nf>=info.nfmax);
              info.qf   = abs((info.ftrue-info.fbest)/...
                          (info.finit-info.fbest));
              info.done = (info.done|info.qf<=info.accf);
              info.sec  = sec;
              fnm=stochasticNM(fs,fbest,ft,tune);
              if ft<fnm, xbest=xt; fbest =ft; 
                 if info.prt>=1
                     disp(' ')
                     disp('heuristicPoint found a better point')
                 end 
                 break;
              elseif i==5
                  [~,ib] = min(ff); fbest = ff(ib); xbest=xx(:,ib);  
                  if info.prt>=1
                     disp(' ')
                     disp('heuristicPoint cannot found a better point')
                     disp(['a point from 5 heuristic points',...
                          ' with lowest inexact f is accepted'])
                  end 
                  break;
              end
              if info.done, return; end  
             end     
         end
     end    
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


