 function [Normal] =CalNormal(Points)
 N= [1,0;
     0,1;
     -1,0;
     0,-1;];
Normal=Points;
 for i = 1: size(Points,1)
     for j = 1:size(Points,2)
        Vec=zeros(4,3);
        for k =1 :4

            if(i+N(k,1)>0&&j+N(k,2)>0&&i+N(k,1)<=size(Points,1)&&j+N(k,2)<=size(Points,2))
                Vec(k,:) = reshape(Points(i+N(k,1),j+N(k,2),:)-Points(i,j,:),1,3);
            end
        end
            M(1,:)= cross(Vec(1,:),Vec(2,:));
            
            M(2,:)= cross(Vec(2,:),Vec(3,:));
            
            M(3,:)= cross(Vec(3,:),Vec(4,:));

            M(4,:)= cross(Vec(1,:),Vec(4,:));
            W = find(M(:,1)~=0&M(:,2)~=0&M(:,3)~=0);
            M=M(W,:);
            M=M./repmat(sqrt(sum(M.^2,2)),1,3);
            Normal(i,j,:) = mean(M);
     end
 end
Normal = Normal./repmat(sqrt(sum(Normal.^2,3)),1,1,3);