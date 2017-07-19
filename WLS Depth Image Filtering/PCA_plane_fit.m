% [normal center_point convex_edge_xyz <cov_mat> <projected_points>] =  PCA_plane_fit(x, y, z, boundary, <outlier_thr>)
function [normal, center_point, conv_edge_xyz, conv_edge_xy, axis_xy, conv_area, inlier_ratio, varargout] = ...
         PCA_plane_fit(x, y, z, edge_idx, varargin);

global min_inlier_ratio sensor_offset

data = [x(:), y(:), z(:)];
mean_data = mean(data,1);
N_point = size(data, 1);
if nargin>4
    outlier_thr = varargin{1};
    [eigen_vects, score, eigen_vals] = princomp(data); 
    outlier_idx = find(abs((data-repmat(mean_data,N_point,1))*eigen_vects(:,3))>outlier_thr);
    data(outlier_idx,:) = [];
    inlier_ratio = size(data, 1)/N_point;
    if inlier_ratio < min_inlier_ratio        % not a plane probably!
        normal = []; center_point = []; conv_edge_xyz = []; conv_edge_xy = []; axis_xy = []; conv_area = [];
        for i=1:nargout-1, varargout{i} = []; end
        return;
    end
    edge_flag = zeros(N_point,1); edge_flag(edge_idx) = 1; edge_flag(outlier_idx) = [];
    edge_idx = find(edge_flag == 1);
    mean_data = mean(data,1);
    N_point = size(data, 1);
else
    inlier_ratio = 1;
end
[eigen_vects, score, eigen_vals] = princomp(data);
normal = eigen_vects(:,3);                                      % Coresponding to the smallest eigen value
axis_xy = eigen_vects(:,1:2);

x = score(edge_idx,1); y = score(edge_idx,2);
[conv_edge_idx, conv_area] = convhull(x, y);

% Projecting on the first two bases (onto the plane)
conv_edge_xyz = (repmat(mean_data,length(conv_edge_idx),1) + [x(conv_edge_idx), y(conv_edge_idx)]*axis_xy')';  
conv_edge_xy = [x(conv_edge_idx), y(conv_edge_idx)]';

if (mean_data-sensor_offset')*normal>0       % Point toward origin, now 3D frame is not right handed any more,
    normal = -normal;       % but the 2D frame in the plane is right handed with respect to robot frame
end

if normal'*cross(axis_xy(:,1), axis_xy(:,2))<0              % not a right handed coordinate system
    axis_xy(:,1) = -axis_xy(:,1);                           % scores should also be changed if used later
end

center_point = mean_data';

switch nargout
    case 8
        varargout{1} = eigen_vals(3);
    case 9
        varargout{1} = eigen_vals(3);
        varargout{2} = repmat(mean_data,N_point,1) + score(:,1:2)*eigen_vects(:,1:2)';  % Projecting on the first two bases
    otherwise
        error(['ERROR: False number of output arguments: ',num2str(nargout),'.']);
end %of switch
