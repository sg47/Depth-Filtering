 function [Normal] =CalNormal1(Points)
x=Points(:,:,1);
y=Points(:,:,2);
z=Points(:,:,3);
[m,n] = size(x);

U1 = cat(3,diff(x,1,2),diff(y,1,2),diff(z,1,2));
V1 = cat(3,diff(x,1,1),diff(y,1,1),diff(z,1,1));
%U1(end,:,:)=[];   % V1(:,end,:)=[];
U2 = cat(3,-diff(x,1,2),-diff(y,1,2),-diff(z,1,2));
V2 = cat(3,-diff(x,1,1),-diff(y,1,1),-diff(z,1,1));

%U2(1,:,:)=[];    V2(:,1,:)=[];
A=U1; B=V1;
A(end,:,:)=[];  B(:,end,:)=[]; 
C1 = cross(B,A);
A=U1; B=V2;
A(1,:,:)=[];  B(:,end,:)=[]; 
C2 = cross(A,B);
A=U2; B=V2;
A(1,:,:)=[];  B(:,1,:)=[]; 
C3 = cross(B,A);
A=U2; B=V1;
A(end,:,:)=[];  B(:,1,:)=[]; 
C4 = cross(A,B);
Normal= zeros(m,n,3);
find(C1(sum(C1.^2,3)==0))
C1=C1./repmat(sqrt(sum(C1.^2,3)),1,1,3);
C2=C2./repmat(sqrt(sum(C2.^2,3)),1,1,3);
C3=C3./repmat(sqrt(sum(C3.^2,3)),1,1,3);
C4=C4./repmat(sqrt(sum(C4.^2,3)),1,1,3);

Normal(1:end-1,1:end-1,:) = C1;
Normal(2:end,1:end-1,:) = Normal(2:end,1:end-1,:)+C2;
Normal(2:end,2:end,:) = Normal(2:end,2:end,:)+C3;
Normal(1:end-1,2:end,:) = Normal(1:end-1,2:end,:)+C4;

Normal(2:end-1,2:end-1,:) = Normal(2:end-1,2:end-1,:)/4;

Normal (1,2:end-1,:)= Normal(1,2:end-1,:)/3;
Normal(end,2:end-1,:) = Normal(end,2:end-1,:)/3;
Normal(2:end-1,1,:) = Normal(2:end-1,1,:)/3;
Normal(2:end-1,end,:) = Normal(2:end-1,end,:)/3;
Normal(1,1,:) = Normal(1,1,:)/2;
Normal(end,end,:) = Normal(end,end,:)/2;
Normal(1,end,:) = Normal(1,end,:)/2;
Normal(end,1,:) = Normal(end,1,:)/2;

Normal = Normal./repmat(sqrt(sum(Normal.^2,3)),1,1,3);