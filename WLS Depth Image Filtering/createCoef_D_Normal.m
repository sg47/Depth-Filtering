function [ out  ] = createCoef_D_Normal(D, Landa, X)
out = zeros (size(X,1), size(X,1));

for i=1: size(X,1)
   
   sumK=0;
   
   for k=1: size(X,1)
       sumJ =0;
       for j=1:size(X,1)
          sumJ = sumJ + D(j,i)*D(j,k); 
       end
   end
   
   if (i==k ) sumK = Landa*sumJ +1; else sumK= Landa *sumJ; end 
   out(i,k) = sumK;
   i
end