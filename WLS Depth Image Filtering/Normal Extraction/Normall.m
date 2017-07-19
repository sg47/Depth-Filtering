function [Normal] = Normall(Points)
if (size(Points,2) == 3)
    Normal = points2normals(Points);
    Normal = Normal';
else
      NormData1(:,1:3) = reshape(Points,size(Points,1)*size(Points,2),3);
      NormData = points2normals(NormData1);
      NormData = NormData';
      Normal(:,:,1) = reshape(NormData(:,1),size(Points,1),size(Points,2));
      Normal(:,:,2)= reshape(NormData(:,2),size(Points,1),size(Points,2));
      Normal(:,:,3)= reshape(NormData(:,3),size(Points,1),size(Points,2));
end