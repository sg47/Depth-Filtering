function a = iRDT(s, d, phi, varargin)

% TBD: take care of update mode (weighted update by sample distance)
%      array input instead of cell input
narg = 3;
ext_mode = 'periodic';
shift = uint32(0);
switch (nargin)
    case narg
        EPS = 0;
    case narg +1
        EPS = varargin{1};
    case narg + 2
        EPS = varargin{1};
        ext_mode = varargin{2};
    case narg + 3
        EPS = varargin{1};
        ext_mode = varargin{2};
        shift = varargin{3};
    otherwise
        warning('iRDT: Error in number of arguments');
end

global SULD LD_Thr;

%s = num2cell(s);
%d = num2cell(d);
%EPS = 0.0005;   % Enlarge to use linear wavelet on details of previous transform

sig_len = length(s{end});
level = length(d);
phi = phi*2^(level-1);
cos_phi = cos(phi);
l = length(d{1});   % l is the length of details in each scale

if strcmp(ext_mode, 'periodic')
    for j=1:level
        % d{j} = d{j}*sqrt(2); s{j} = s{j}*sqrt(2);
        for i=1:l
            n = mod(i, l)+1;      % next index in details
            if SULD == 1
                dp = sign(d{j}(i))*min(abs(d{j}(i)), LD_Thr);
                dn = sign(d{j}(n))*min(abs(d{j}(n)), LD_Thr);
            else
                dp = d{j}(i);
                dn = d{j}(n);
            end
            s{j+1}(2*i) = s{j}(i) - (dp+dn)/4;      
        end
        for i=1:l
            p = mod(2*i-3, 2*l)+1; n = 2*i; % previous and next indices to i
            if abs(s{j+1}(p)+s{j+1}(n)) > EPS   % to avoid devide by zero or very small denominator
                s{j+1}(mod(2*i-2, 2*l)+1) = d{j}(i) + 2*s{j+1}(p)*s{j+1}(n)*cos_phi/(s{j+1}(p)+s{j+1}(n));
            else   % in tiny triangles, use zero prediction or normal average
                s{j+1}(mod(2*i-2, 2*l)+1) = d{j}(i);% d{j}(i) + (s{j+1}(p)+s{j+1}(n))/2;%
            end
        end
        if bitget(shift,level-j+1)  % odd samples should be pushed to even positions
            s{j+1} = s{j+1}([2:end,1]);
        end
        phi = phi/2; cos_phi = cos(phi);
        l = ceil(sig_len/(2^(level-j)));
        s{j+1} = s{j+1}(1:l);      % removed temprarily for array input
     %   s{j+1}(find((isnan(s{j+1})==1)))=0;
    end
else
    for j=1:level
        l = length(d{j});   % l is the length of details in each scale
        % d{j} = d{j}*sqrt(2); s{j} = s{j}*sqrt(2);
        for i=2:l
            s{j+1}(2*i) = s{j}(i-1) - (d{j}(i-1)+d{j}(i))/4;
        end
        for i=2:l-1 %%% loop variable start change by mahdi aghaee 1->2
            if abs(s{j+1}(2*i)+s{j+1}(2*i+2)) > EPS   % to avoid devide by zero or very small denominator
                s{j+1}(2*i+1) = d{j}(i) + 2*s{j+1}(2*i)*s{j+1}(2*i+2)*cos_phi/(s{j+1}(2*i)+s{j+1}(2*i+2));
            else   % in tiny triangles, use zero prediction or normal average
                s{j+1}(2*i+1) = d{j}(i);% d{j}(i) + (s{j+1}(p)+s{j+1}(n))/2;%
            end
        end
        phi = phi/2; cos_phi = cos(phi);
        %l = ceil(sig_len/(2^(level-j)))+2;  
        s(2*l+1:end) = []; %%% insert by mahdi aghaee
        if bitget(shift,level-j+1)  % odd samples should be pushed to even positions
            s{j+1}(1:3) = [];
        else                % odd samples should remain in odd positions
            s{j+1}(1:4) = [];
        end
    end
    %%% s{level+1} = s{level+1}(1:sig_len); comment by mahdi aghaee
end
a = s{level+1};

return