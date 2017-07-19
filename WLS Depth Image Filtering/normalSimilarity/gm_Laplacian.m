function [ out  ] = gm_Laplacian( rangePatch,normalPatch )

[m,n]=size(rangePatch);
rangePatch=double(rangePatch);
normalPatch=double(normalPatch);
Nx=normalPatch(:,:,1); Ny=normalPatch(:,:,2); Nz=normalPatch(:,:,3);
normalPatchVector= [Nx(:),Ny(:),Nz(:)];
out = zeros(size(rangePatch));

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
normalPatch_padded = padarray(normalPatch,[1,1],0);
w=3;


T_normal_weight= ones(size(rangePatch_padded,1)*size(rangePatch_padded,2),m*n);
T_normal_weight_col = ones(size(rangePatch_padded,1)*size(rangePatch_padded,2),1);




col_index= m+2;
for i=floor(w/2)+1:size(rangePatch,1)-floor(w/2)
    for j=floor(w/2)+1:size(rangePatch,2)-floor(w/2)
        outPoint=zeros(double(1));
        outPoint=zeros(double(1));
        classCenter = ind_patch(i,j);
        sumSimilarity=zeros(double(1));
        
        weight=zeros(m,n);
        for ii=i-floor(w/2):i+floor(w/2)
            for jj=j-floor(w/2):j+floor(w/2)
               % similarity = exp(-(normalPatch(ii,jj,:)-guassian_obj.mu(classCenter,:)).^2/(2*guassian_obj.Sigma(:,:,classCenter)^2));
               %similarity = gauss_distribution( normalPatch(ii,jj,:),guassian_obj.mu(classCenter,:),guassian_obj.Sigma(:,:,classCenter));
               similarity = gauss_distribution( normalPatch(ii,jj,:),mu(classCenter,:),sigma(:,:,classCenter));
               
              % similarity =1;
               %outPoint = outPoint+ rangePatch(ii,jj)*similarity;
               weight(ii,jj) =similarity;
               if (ii==i && jj==j) weight(ii,jj)=0; end
               %sumSimilarity = sumSimilarity + similarity;
            end
        end
        
       weight = -1 * weight;
       weight(i,j) =  sum(-1 * weight(:));
       
       weight_paded = padarray(weight,[1,1],0)';
       T_normal_weight_col = weight_paded(:);
       T_normal_weight(:,col_index) = T_normal_weight_col;
       
       if (mod(col_index+1,(m))==0)  col_index=col_index+2; end
       
       col_index=col_index+1;
       
       %weight=im2bw(weight);
       %weight = fspecial('gaussian',5,25);
       %if (max(max(weight))==0) weight = ones(w); end
       % imagesc(weight);  %if(i==13 && j==20)
       % out(i,j) = outPoint/sumSimilarity;
       % out(i,j) = rangePatch(i,j);
       
    end
end

%% Laplacian Fitler
Landa=0.1;

x1=rangePatch(:);
D=[0 -1 0;-1 4 -1;0 -1 0];
D3=[-1 -1 -1;-1 8 -1;-1 -1 -1];

D2=[0 0 -1 0 0; 0 -1 -1 -1 0; -1 -1 12 -1 -1; 0 -1 -1 -1 0; 0 0 -1 0 0];

T=convmtx2(D,[m,n]);
%T = T .* (T_normal_weight);
T =T_normal_weight;
T1=T'*T;
I1=speye(size(T1));
J = double(rangePatch);
Res=inv(I1+Landa*T1)*J(:);
%ResNorm=(Res-min(Res))/(max(Res)-min(Res));
out=reshape(Res,[m,n]);


end

