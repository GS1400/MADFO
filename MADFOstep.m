%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% MADFOstep.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% performs one step of MADFO
% 
%
% function MADFOstep(tune,fun,prt,n); % initialize solver
%
% tune         % structure contains those tuning parameters that are 
%              %   explicitly specified; the others take default values
%              %   (for choices and defaults see initTune.m)
% fun          % function handle (empty if used in reentrant mode)
% prt          % print level
% n            % problem dimension
%
% function x=MADFOstep(x,f,prt); % suggest new point
%
% x            % point evaluated (input) or to be evaluated (output)
% f            % function value at input x
% prt          % change print level (optional)
%
% function [xbest,fbest,info]=MADFOstep();  % read results
%
% xbest        % current best point
% fbest        % function value found at xbest
% info         % performance information for MADFO
%  .finit      %   initial function value
%  .ftarget    %   target function value (to compute qf)
%  .qf         %   (ftarget-f)/(finit-f)
%  .initTime   %   inital cputime
%  .done       %   done with the search?
% 
function [x,f,info1] = MADFOstep(x,f,p,n)
persistent Vstate tune info  iter dim  nf prt
% persistent variables for computing s, sigma, and M
persistent  wi d_s echi c_s sqrt_s c_1 c_mu s sigma M  
% persistent variables for the recombination phase
persistent RecomInfo sigma0
% persistent variables for the mutation phase
persistent Mutation MutationInfo lmut
% persistent variables for the selection phase
persistent Fmut ranks
% persistent variables for computing non-monotone term
persistent Fnm  fnm
% persistent variables for extrapolation
persistent F alpha sigmae netxrapol ie
% persistent variables for heuristic
persistent Xt Ft x13 x12 x1 nhp
% persistent variables for best point and its function value
persistent xbest fbest  xfinal ffinal
if nargin==0
  %%%%%%%%%%%%%%%%%%%%%
  %%%% return info %%%%
  %%%%%%%%%%%%%%%%%%%%%
  if ffinal<fbest, f=ffinal; x=xfinal;
  else, f=fbest; x=xbest;
  end
  info.nfuncEvaluations=nf;
  info.nextrapolations=netxrapol;
  info.nheuristicPoint=nhp;
  info1=info;
  return;
end;
if nargout==0
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%% initialize solver environment %%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  info = []; tune = initTune(x,n); 
  prt = p; info1=info; dim = n; Vstate=1; 
  return;
end
if Vstate==1 % the initial function value was computed
    [wi,sigma,iter,d_s,echi,c_s,sqrt_s,nhp,RecomInfo,Fnm,... 
           netxrapol, c_1,c_mu,s,M,Xt,Ft]=initMADFO(dim,tune);
    xbest = x; f= max(-1e50,min(f,1e50)); fbest = f; nf = 1;
    Vstate=2; % go to start the first mutation phase
end
% save the evaluated point and its function value
xfinal = x; ffinal = f;
while 1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%% Mutation phase %%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      if Vstate==2 % initialization for the mutation phase;
                   % start of the mutation phase
         if iter>1, RecomInfo.d0=RecomInfo.d; end
         MutationInfo=[]; Mutation=[];  ranks=[]; lmut = 1; Vstate=3; 
          if prt>=1
              disp(' ')
              disp('+++++++++++++++++++++++++++++++++++++++++++++')
              disp(['start of ',num2str(iter),'th mutation phase'])
              disp('+++++++++++++++++++++++++++++++++++++++++++++')
              disp(' ')
          end
      end
      while Vstate==3 || Vstate==4  % loop for the mutation phase
          if Vstate==3
             % distribution direction is computed
             Mutation.z = randn(dim, 1); 
             % mutation direction is computed
             Mutation.d = M*Mutation.z;
             sigma0 = sigma;
             if tune.ver~=1 % good mutation step size 
                sigma=goodStepSize(Vstate,xbest,Mutation,...
                                 RecomInfo,sigma0,sigma,tune,nf,prt);
             end
              % mutation point is computed
              xs  = xbest + sigma*Mutation.d; x = xs;
              Vstate=4; % go to compute f at x
              return;
          end
          if Vstate==4 % f at x was computed
              Mutation.f = max(-1e50,min(f,1e50)); nf=nf+1; 
              MutationInfo = [MutationInfo,Mutation]; lmut=lmut+1;
              if lmut>tune.lambda
                  Vstate=5; % end of the mutation phase
                  if prt>=1
                      disp(' ')
                      disp('+++++++++++++++++++++++++++++++++++++++++++')
                      disp(['end of ',num2str(iter),'th mutation phase'])
                      disp('+++++++++++++++++++++++++++++++++++++++++++')
                      disp(' ')
                  end
                  break;
              else
                  Vstate=3; % go back to continue the mutation phase
              end
          end
      end
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%%%%%%%%%%%%%%%%%%%%%%%%%%% Selection phase %%%%%%%%%%%%%%%%%%%%%
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      if Vstate==5
          if prt>=1
              disp(' ')
              disp('+++++++++++++++++++++++++++++++++++++++++++++++')
              disp(['start of ',num2str(iter),'th selection phase'])
              disp('+++++++++++++++++++++++++++++++++++++++++++++++')
              disp(' ')
          end
          % sort lambda function values at the mutation point
          Fmut=[];
          for l=1:tune.lambda, Fmut(l) = MutationInfo(l).f; end
          [~, ranks] = sort(Fmut,'ascend');
          % select mu distribution and mutation directions
          % compute recombination mutation and distribution directions
           sum_z = zeros(dim,1); sum_d = zeros(dim,1);
           for m = 1:tune.mu
               sum_z = sum_z + wi(m)*MutationInfo(ranks(m)).z;
               sum_d = sum_d + wi(m)*MutationInfo(ranks(m)).d;
           end
           RecomInfo.z = sum_z; % recd direction was computed
           RecomInfo.d = sum_d; % recm direction was computed
           Vstate=6;
           if prt>=1
              disp(' ')
              disp('+++++++++++++++++++++++++++++++++++++++++++++')
              disp(['end of ',num2str(iter),'th selection phase'])
              disp('+++++++++++++++++++++++++++++++++++++++++++++')
              disp(' ')
          end
      end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%% Recombination phase %%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % update Fnm and compute the recm subspace direction  
    if Vstate==6 
        if prt>=1
           disp(' ')
           disp('++++++++++++++++++++++++++++++++++++++++++++++++++')
           disp(['start of ',num2str(iter),'th recombination phase'])
           disp('++++++++++++++++++++++++++++++++++++++++++++++++++')
           disp(' ')
        end
        if tune.ver~=3 % update the list Fnm required to compute fnm
                       % this is done by updateNM
           Fnm=updateNM(Fnm,Fmut,iter,tune);
        end
        if tune.ver==2 % no subspace direction
           if prt>=2
              disp(' ')
              disp('no recm subspace direction')
          end 
        else %  computation of recm subspace direction
             if iter>1
             RecomInfo=recomSubDir(Xt,Ft,RecomInfo,dim,iter);
                if prt>=2
                  disp(' ')
                  disp('recm subspace direction was computed')
                end 
            end
        end
         if tune.ver~=1 % good recombination step size
            sigma = goodStepSize(Vstate,xbest,Mutation,RecomInfo,...
                                 sigma0,sigma,tune,nf,prt);
         end
        sigmae = sigma;
        % recombinstion point is computed
        RecomInfo.y = xbest + sigmae*RecomInfo.d; x = RecomInfo.y;
        Vstate=7;
        return % go to compute f at x
    end
    %%%%%%%%%%%%%%%%%%%%%%%
    % start of randomNMLS %
    %%%%%%%%%%%%%%%%%%%%%%%
    if Vstate==7 && tune.ver==3 % no non-monotone line search
        RecomInfo.f=  max(-1e50,min(f,1e50)); nf=nf+1;
        if RecomInfo.f<fbest % decrease in f found
            Vstate=14;  % go to update s, sigma and M
            xbest = RecomInfo.y; fbest = RecomInfo.f;
        else
            Vstate=12; % go to perform heuPoint
        end
    end
    if Vstate==7 && tune.ver~=3 % non-monotone line search
        RecomInfo.f= max(-1e50,min(f,1e50)); nf=nf+1;
        fnm = randNM(Fnm,fbest,RecomInfo.f,tune);
        if fnm>RecomInfo.f+tune.gamma*sigmae^2 % a decrease in f is found
            F=RecomInfo.f; alpha=sigmae; ie=0; Vstate=8; 
            % go to do extrapolation
        else
            Vstate=11; % opposite direction is tried
        end
    end
     while Vstate==8 || Vstate==9 % extrapolation loop
           if Vstate==8
                sigmae = tune.gammaE*sigmae; alpha=[alpha,sigmae];
                RecomInfo.y=xbest+sigmae*RecomInfo.d; x=RecomInfo.y;
                Vstate=9;
                return; % go to compute f at x
           end
           if Vstate==9 % f at x was computed
             RecomInfo.f = max(-1e50,min(f,1e50)); 
             F=[F,RecomInfo.f]; nf=nf+1; ie = ie+1;
             nodec = fnm<RecomInfo.f+tune.gamma*sigmae^2;
             if nodec || ie==tune.E % end of extrapolation
                                    % update the best point and traingle
                netxrapol = netxrapol+1;
                [~,ib] = min(F); ib=ib(1); sigmae=alpha(ib);
                fbest=F(ib); xbest=xbest+sigmae*RecomInfo.d;
                % update the vertices of the first triangle
                [Xt,Ft]=triUpdate1(Xt,Ft,xbest,fbest,tune);
                Vstate=14; % go to update s, sigma and M
                F=[]; alpha=[];
                 if prt>=2
                   disp(' ')
                   disp(['extrapolation was terminated,',...
                         ' resulting in a new best point'])
                 end 
                if prt>=1
                   disp(' ')
                   disp('+++++++++++++++++++++++++++++++++++++++')
                   disp(['end of ',num2str(iter),...
                         'th recombination phase'])
                   disp('+++++++++++++++++++++++++++++++++++++++')
                   disp(' ')
                end
                break;
             else
                Vstate=8; % go back to go on extrapolation
             end
           end
     end
    % opposite recom,bination direction is tried
    if Vstate==11 && tune.ver~=3 
        RecomInfo.d = -RecomInfo.d; 
        RecomInfo.y = xbest +sigmae*RecomInfo.d; x = RecomInfo.y;
        Vstate=88;
         if prt>=2
            disp(' ')
            disp('opposite recm subspace direction was tried')
          end 
        return; % go to compute f at x
    end
    if Vstate==88
         RecomInfo.f= max(-1e50,min(f,1e50)); nf=nf+1;
         if fnm>RecomInfo.f+tune.gamma*sigmae^2
             F=RecomInfo.f; alpha=sigmae; ie=0;  
             Vstate=888; % a decrease in f is found
         else
             Vstate=12; % no decrease in f 
         end
    end
    while Vstate==888 || Vstate==999 % extrapolation loop
          if Vstate==888 % expand extrapolation step size
              sigmae = tune.gammaE*sigmae; alpha=[alpha,sigmae];
              RecomInfo.y = xbest + sigmae*RecomInfo.d;
              x = RecomInfo.y;
              Vstate=999;
              return; % go to compute f at x
           end
           if Vstate==999 % f at x was computed
             RecomInfo.f = max(-1e50,min(f,1e50)); 
             F=[F,RecomInfo.f]; nf=nf+1; ie = ie+1;
             nodec = fnm<RecomInfo.f+tune.gamma*sigmae^2;
             if nodec || ie==tune.E % end of extrapolation
                netxrapol = netxrapol+1;
                [~,ib] = min(F); ib=ib(1); sigmae=alpha(ib);
                fbest=F(ib); xbest=xbest+sigmae*RecomInfo.d;
                % update the vertices of the first triangle
                [Xt,Ft]=triUpdate1(Xt,Ft,xbest,fbest,tune);
                Vstate=14; % go to update s, sigma and M
                F=[]; alpha=[];
                if prt>=2
                   disp(' ')
                   disp(['extrapolation was terminated,', ...
                         ' resulting in a new best point'])
                end 
                if prt>=1
                   disp(' ')
                   disp('+++++++++++++++++++++++++++++++++++++')
                   disp(['end of ',num2str(iter),...
                         'th recombination phase'])
                   disp('++++++++++++++++++++++++++++++++++++')
                   disp(' ')
                end
                break;
             else
                Vstate=888; % go back to go on extrapolation
             end
           end
    end
    %%%%%%%%%%%%%%%%%%%%%
    % end of randomNMLS %
    %%%%%%%%%%%%%%%%%%%%%
   if Vstate==12 % no decrease in f found;
        if prt>=1
           disp(' ')
           disp('++++++++++++++++++++++++++++++++++++++++++++++++++')
           disp(['end of ',num2str(iter),'th recombination phase'])
           disp('++++++++++++++++++++++++++++++++++++++++++++++++++')
           disp(' ')
        end
       % try to find one heuristic point
        if RecomInfo.f < fbest % a decrease in f was found
            % update best point and its function value
            fbest= RecomInfo.f; xbest = RecomInfo.y; 
            % update the vertices of the first triangle
            [Xt,Ft]=triUpdate1(Xt,Ft,xbest,fbest,tune);
            Vstate=14; % go to update sigma and M 
             if prt>=2
               disp(' ')
               disp('a decrease in f was found')
             end 
        else
            if size(Xt,2)>=3 && tune.ver~=4 
                % vertices of the second triangle are obtained
                [x1,x13,x12,Xt,Ft]=triUpdate2(Xt,Ft);
                Vstate=15;
           else % no heuPoint; because the first triangle 
                % has not been formed yet
               Vstate=14;
                if prt>=2
                  disp(' ')
                  disp('no random triangle subspace point is tried')
                end 
            end
        end
   end
   if Vstate==15 % triangle subsapce subspace point inside Delta2
                 % Delta2 is the second triangle
       xt = TriSubPoint(x1,x13,x12); x=xt; nhp = nhp+1;
       if prt>=2
            disp(' ')
            disp('random triangle subspace point was computed')
       end 
       Vstate=16; 
       return; % go to compute f at x
   end 
   if Vstate==16 % f at x was computed; check whether or not
                 % x can be a new best point
       f = max(-1e50,min(f,1e50)); nf= nf+1;
       % the second decrease condition in f is checked
       if f<fbest, xbest=x; fbest =f; % best point is updated
           if prt>=2
               disp(' ')
               disp('random triangle subspace point was accpeted')
           end 
       end
       Vstate=14; % go to update s, sigma and M
   end
   if Vstate==14 % update s, sigma and M
        iter = iter+1;  
        s = (1-c_s)*s + sqrt_s*RecomInfo.z;
        M = (1 - 0.5*c_1 - 0.5*c_mu) * M + (0.5*c_1)*(M*s)*s';
        for m = 1:tune.mu
          M = M + ((0.5*c_mu*wi(m))*MutationInfo(ranks(m)).d)*...
          MutationInfo(ranks(m)).z';
        end
         pow = (c_s/d_s)*((norm(s)/echi)-1); sigma = sigma*exp(pow);
         if prt>=2
            disp(' ')
            disp('s, sigma and M were updated')
         end 
         Vstate=2; % go back a new mutation phase
   end
end % of MADFO
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

