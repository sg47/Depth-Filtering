function n_C = snn_win(x,y,z,winSize)
% function normC = snn(x,y,z) calculate surface normal norm of depth map
[m,n] = size(x);
U1 = cat(3,diff(x,1,2),diff(y,1,2),diff(z,1,2));
V1 = cat(3,diff(x,1,1),diff(y,1,1),diff(z,1,1));
U1(end,:,:)=[];    V1(:,end,:)=[];

h = fspecial('average',[winSize,winSize]);
U1 = imfilter(U1,h,'replicate');
V1 = imfilter(V1,h,'replicate');

U2 = cat(3,-diff(x,1,2),-diff(y,1,2),-diff(z,1,2));
V2 = cat(3,-diff(x,1,1),-diff(y,1,1),-diff(z,1,1));
U2(1,:,:)=[];    V2(:,1,:)=[];

U2 = imfilter(U2,h,'replicate');
V2 = imfilter(V2,h,'replicate');

C1 = cross(U1,V1);
C2 = cross(U2,V2);

C = ones(m,n,3);

C(1:m-1,1:n-1,:) = C1;

C(2:m,end,:) = C2(:,end,:);
C(end,2:n,:) = C2(end,:,:);

% C(2:m-1,2:n-1,:) = (C1(2:end, 2:end,:) + C2(1:end-1, 1:end-1,:))/2;

winSize = sqrt(C(:,:,1).^2 + C(:,:,2).^2 + C(:,:,3).^2);

n_C = C./cat(3,winSize,winSize,winSize);
