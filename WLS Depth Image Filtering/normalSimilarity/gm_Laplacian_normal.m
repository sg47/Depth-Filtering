function [ out  ] = gm_Laplacian_normal( rangePatch,normalPatch )

[m,n]=size(rangePatch);
rangePatch=double(rangePatch);
normalPatch=double(normalPatch);
Nx=normalPatch(:,:,1); Ny=normalPatch(:,:,2); Nz=normalPatch(:,:,3);
normalPatchVector= [Nx(:),Ny(:),Nz(:)];
out = zeros(size(rangePatch));
outNormal = zeros(size(rangePatch));

classNum = 2;
guassian_obj = gmdistribution.fit(normalPatchVector,classNum);

[idx2,c]=kmeans(normalPatchVector,classNum);
[mu,sigma]=GMM_parameter(normalPatchVector,idx2,classNum);

%%
[idx,nlogl] = cluster(guassian_obj,normalPatchVector);
%idx = idx';
ind_patch = reshape(idx2,size(rangePatch));
%imagesc(ind_patch);
%imshow(normalPatch);
% rangePatch = zeros(5);
% [m,n] = size(rangePatch);

rangePatch_padded = padarray(rangePatch,[1,1],0);
rangePatch_padded = rangePatch;
normalPatch_padded = padarray(normalPatch,[1,1],0);
w=3;


T_normal_weight= ones(size(rangePatch_padded,1)*size(rangePatch_padded,2),m*n);
T_normal_weight_col = ones(size(rangePatch_padded,1)*size(rangePatch_padded,2),1);
B = ones(size(rangePatch_padded,1)*size(rangePatch_padded,2),1);

%%
D=[0 -1 0;-1 4 -1;0 -1 0];
D3=[-1 -1 -1;-1 8 -1;-1 -1 -1];
D2=[0 0 -1 0 0; 0 -1 -1 -1 0; -1 -1 12 -1 -1; 0 -1 -1 -1 0; 0 0 -1 0 0];
Landa=0.1;
T=convmtx2(D,[m,n]);
T1=T'*T;
I1=speye(size(T1));
B = I1+Landa*T1;  % B initializing with previous DX

%%
T_normal_weight= T1;
T_normal_weight_col = ones(size(rangePatch_padded,1)*size(rangePatch_padded,2),1);

col_index= m+2;
for i=floor(w/2)+1:size(rangePatch,1)-floor(w/2)-1
    for j=floor(w/2)+1:size(rangePatch,2)-floor(w/2)-1
        outPoint=zeros(double(1));
        outPoint=zeros(double(1));
        classCenter = ind_patch(i,j);
        sumSimilarity=zeros(double(1));
        
        weight=zeros(m,n);
        for ii=i-floor(w/2):i+floor(w/2)
            for jj=j-floor(w/2):j+floor(w/2)
              
                similarity = gauss_distribution( normalPatch(ii,jj,:),mu(classCenter,:),sigma(:,:,classCenter));
               
              % similarity =1;
               weight(ii,jj) =1;
               outPoint = outPoint+ rangePatch(ii,jj)*similarity;
               sumSimilarity = sumSimilarity + similarity;
            end
        end

       weight = -1 * weight;
       weight(i,j) =  sum(-1 * weight(:))-1;
       weight_paded = weight;
       T_normal_weight_col = weight_paded(:);
       T_normal_weight(:,col_index) = T_normal_weight_col;
       
       if (mod(col_index+1,(m))==0)  col_index=col_index+2; end
       
       col_index=col_index+1;
       
      
       % imagesc(weight);  %if(i==13 && j==20)
       outNormal(i,j) = outPoint/sumSimilarity;
       % out(i,j) = rangePatch(i,j);
       
    end
end

%% Laplacian Fitler
x1=rangePatch(:);
%B2 = createCoef_D_Normal_fast(T_normal_weight, Landa, x1);
B2 = createCoef_D_Normal_fast2(T1, Landa, x1);

x1 = double(rangePatch);
%C = T_normal_weight'*T_normal_weight * x1(:);
%C =  x1(:);
C=outNormal(:);
Res=inv(B2)*C;
out=reshape(Res,[m,n]);
%out = outNormal;
end

