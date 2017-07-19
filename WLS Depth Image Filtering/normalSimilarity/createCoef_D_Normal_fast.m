function [ out  ] = createCoef_D_Normal_fast(D, Landa, X)
out = zeros (size(X,1), size(X,1));

for i=1: size(X,1)
   
   
   
   for k=1: size(X,1)
       sumK=0;
       
       sumJ =sum(D(:,i).*D(:,k));
       if (i==k ) sumK = Landa*sumJ +1; else sumK= Landa *sumJ; end 
       out(i,k) = sumK;
   end
   
end