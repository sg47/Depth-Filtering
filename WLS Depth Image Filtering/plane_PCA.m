function [normal,varargout] = plane_PCA(varargin);

% [normal,{basis_vect_1,basis_vect_2,center,std_dev}] = ...
%   plane_PCA({mat3D}|{Xe,Ye,Ze}|{XYZe})
%
% Best fitting plane computed with principal component analysis.
%
% Ref.: getBestFittingPlanePCA.m (jan.weingarten@epfl.ch)
%
% 2004.8.19 stefan.gaechter@epfl.ch

switch nargin,
  case 1,
    if size( varargin{1}, 2) == 3,
      XYZe = varargin{1};
      Xe = XYZe(:,1);
      Ye = XYZe(:,2);
      Ze = XYZe(:,3);
    else,
      Xe = varargin{1}.X;
      Ye = varargin{1}.Y;
      Ze = varargin{1}.Z;
    end; %of if
  case 3,
    Xe = varargin{1};
    Ye = varargin{2};
    Ze = varargin{3};
  otherwise,
    error(['ERROR: False number of input arguments: ',num2str(nargin),'.']);
end; %of switch

Xe = Xe(:);
Ye = Ye(:);
Ze = Ze(:);

[m,n]=size(Xe);
for i=1:m
    for j=1:n
        if (isnan(Xe(i,j)))
         Xe(i,j)=0;
            
        end
        if (isnan(Ye(i,j)))
         Ye(i,j)=0;
            
        end
        if (isnan(Ze(i,j)))
         Ze(i,j)=0;
            
        end
    end
end

% Mean Vector (Center of Gravitiy)
xc = mean(Xe);
yc = mean(Ye);
zc = mean(Ze);

% Covariance Matrix
SS = cov( [Xe-xc, Ye-yc, Ze-zc] );

% Eigenvalues and Eigenvectors
[eigen_vects, eigen_vals_mat] = eig(SS);
[eigen_vals, indices] = sort( diag( eigen_vals_mat ) );
eigen_vects = eigen_vects( :, indices );

% Surface Normal
%  Eigenvector corresponding to smallest eigenvalue.
normal = eigen_vects(:,1);
basis_vect_2 = eigen_vects(:,2);
basis_vect_1 = eigen_vects(:,3);

% Normalization
norm_factor = norm( normal );
if norm_factor > 0,
  normal = normal / norm_factor;
end; %of if
norm_factor = norm( basis_vect_1 );
if norm_factor > 0,
  basis_vect_1 = basis_vect_1 / norm_factor;
end; %of if
norm_factor = norm(basis_vect_2);
if norm_factor > 0,
  basis_vect_2 = basis_vect_2/norm_factor;
end; %of if

% Right handed coord. system.
if normal'*cross( basis_vect_1, basis_vect_2 ) < 0,
  normal = -normal;
end; %of if

switch nargout,
  case 0,
    % NOP
  case 1,
    % NOP
  case 2,
    varargout{1} = [xc, yc, zc]';
  case 3,
    varargout{1} = basis_vect_1;
    varargout{2} = basis_vect_2;
  case 4,
    varargout{1} = basis_vect_1;
    varargout{2} = basis_vect_2;
    varargout{3} = [xc, yc, zc]';
  case 5,
    varargout{1} = basis_vect_1;
    varargout{2} = basis_vect_2;
    varargout{3} = [xc, yc, zc]';
    varargout{4} = sqrt( [eigen_vals(1), eigen_vals(3), eigen_vals(2)]' );
  otherwise,
    error(['ERROR: False number of output arguments: ',num2str(nargout),'.']);
end; %of switch
