 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% recomSubDir %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% recomSubDir compute recm subspce direction
%
function RecomInfo=recomSubDir(Xt,Ft,RecomInfo,dim,iter)
if length(Ft)==3 % compute scaling vector
    ib=1; II = setdiff(1:3,ib); XX = Xt(:,II);
    for i=1:2, dX(:,i)  = XX(:,i)-Xt(:,ib);end
    sc = max(dX')'; sc(sc==0) = 1;
else
   sc  = ones(dim,1);
end
randfac=rand;
if randfac>0.5 
    beta1=randfac; beta2 = sqrt(1-randfac^2)*exp(-iter); 
else
 beta1=sqrt(1-randfac^2); beta2 = randfac*exp(-iter);
end
d = beta1*RecomInfo.d + beta2*RecomInfo.d0; 
RecomInfo.d = sc.*d; % recm subspace direction was computed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%