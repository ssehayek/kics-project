%
function out = ICSFitGen(params,xiLags,etaLags,varargin)

errBool = 0;
includeZeroLag = 0;

for i = 1:length(varargin)
    if any(strcmpi(varargin{i},{'error','err','residual','res'}))
        errBool = 1;
        corrData = varargin{i+1};
    elseif any(strcmpi(varargin{i},{'includeZeroLag','zeroLag'}))
        includeZeroLag = 1;
    end
end

if errBool && ~includeZeroLag
    zeroLagInd = find(xiLags==0 & etaLags==0);
    xiLags(zeroLagInd) = [];
    etaLags(zeroLagInd) = [];
   
    corrData(zeroLagInd) = [];    
end

s = struct('amp',params(1),'offset',params(2),'w1',params(3),'w2',params(4),'theta',params(5));

if errBool 
    nan_inds = find(isnan(corrData));
    
    xiLags(nan_inds) = [];
    etaLags(nan_inds) = [];
    corrData(nan_inds) = [];
end

% fields = fieldnames(s);
%
% for n = 4:length(fields)
%     if n <= length(params)   
%         fieldname = fields{n};
%         s.(fieldname) = params(n);
%     end
% end

% K = k_on + k_off;
    
a = cos(s.theta).^2/(s.w1.^2)+sin(s.theta).^2/(s.w2.^2);
b = -sin(2*s.theta)/(2*s.w1.^2)+sin(2*s.theta)/(2*s.w2^2);
c = sin(s.theta).^2/(s.w1.^2)+cos(s.theta).^2/(s.w2.^2);

F = s.amp*exp(-a*xiLags.^2+2*b*xiLags.*etaLags-c*etaLags.^2)+s.offset;

% K/k_on*1/(s.density*pi*min(s.w1,s.w2)^2)

if errBool
    err = norm(corrData - F,2);
end

if ~errBool
    out = F;
else
    out = err;
end
