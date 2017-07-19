function [ out  ] = gmdistFilter( rangePatch,normalPatch, filterstatus )


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
w=5;
for i=floor(w/2)+1:size(rangePatch,1)-floor(w/2)-1
    for j=floor(w/2)+1:size(rangePatch,2)-floor(w/2)-1
         outPoint=zeros(double(1));
        outPoint=zeros(double(1));
        classCenter = ind_patch(i,j);
        sumSimilarity=zeros(double(1));
        
        weight=zeros(w);
        for ii=i-floor(w/2):i+floor(w/2)
            for jj=j-floor(w/2):j+floor(w/2)
               % similarity = exp(-(normalPatch(ii,jj,:)-guassian_obj.mu(classCenter,:)).^2/(2*guassian_obj.Sigma(:,:,classCenter)^2));
               %similarity = gauss_distribution( normalPatch(ii,jj,:),guassian_obj.mu(classCenter,:),guassian_obj.Sigma(:,:,classCenter));
               similarity = gauss_distribution( normalPatch(ii,jj,:),mu(classCenter,:),sigma(:,:,classCenter));
               
               if (filterstatus) similarity =1; end
               
               outPoint = outPoint+ rangePatch(ii,jj)*similarity;
               weight(floor(w/2)+(ii-i)+1,floor(w/2)+(jj-j)+1) =similarity;
               sumSimilarity = sumSimilarity + similarity;
            end
        end

        
       %weight=im2bw(weight);
       %weight = fspecial('gaussian',5,25);
       %if (max(max(weight))==0) weight = ones(w); end
       %if (filterstatus) weight = ones(5); end
       out(i,j) = sum(sum(rangePatch(i-floor(w/2):i+floor(w/2),j-floor(w/2):j+floor(w/2)).*weight))/sum(sum(weight));
       % imagesc(weight);  %if(i==13 && j==20)
        out(i,j) = outPoint/sumSimilarity;
       % out(i,j) = rangePatch(i,j);
    end
end


end

