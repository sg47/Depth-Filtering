function pad = signalPadding(signal,method, varargin)
if strcmp(method, 'linear')
    % At the moment for non-periodic extension, linear padding is used
    if nargin>2
        cos_phi = varargin{1};
    else
        disp('Error: 3 argument used in linear Padding');
    end
    t3 = signal(1)*signal(2)/(2*signal(2)*cos_phi-signal(1));
    t2 = t3*signal(1)/(2*signal(1)*cos_phi-t3);
    t1 = t2*t3/(2*t3*cos_phi-t2);
    t4 = signal(end)*signal(end-1)/(2*signal(end-1)*cos_phi-signal(end));
    t5 = t4*signal(end)/(2*signal(end)*cos_phi-t4);
    t6 = t5*t4/(2*t4*cos_phi-t5);
    t7 = t6*t5/(2*t5*cos_phi-t6);
    pad = [t1; t2; t3; t4; t5; t6; t7];
elseif strcmp(method, 'mirror')
    % At the moment for non-periodic extension, mirror padding is used
    pad = [signal(4:-1:2); signal(end-1:-1:end-4)];
elseif strcmp(method, 'const')
    % At the moment for non-periodic extension, constant padding is used
    tBegin = mean(signal(1:3));
    tEnd = mean(signal(end-2:end));
    pad = [ones(3,1)*tBegin ; ones(4,1)*tEnd];
end