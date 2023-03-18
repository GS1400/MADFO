 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% updateNM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% updateNM Fnm that is used to compute fnm
%
function Fnm=updateNM(Fnm,Fmut,iter,tune)
if iter==1, Fnm=Fmut; % initialization was done
else % update of Fnm is done
    lenFnm = length(Fnm);
    if lenFnm<tune.kappa*tune.lambda, Fnm=[Fnm Fmut];
    else % pick randomly lambda indices and update Fnm
       Inm = randperm(lenFnm,tune.lambda); Fnm(Inm) = Fmut;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%