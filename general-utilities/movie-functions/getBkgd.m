function [bkgd_noise] = getBkgd(J,varargin)

bkgd_type = 'getMean';
bkgd_rgn = [];

for ii = 1:length(varargin)
    if any(strcmpi(varargin{ii},{'type','bkgdType','method','routine'}))
        if any(strcmpi(varargin{ii+1},{'getMean','mean','meanOnly'}))
            % already default
        elseif any(strcmpi(varargin{ii+1},{'includeBleach','useBleach'...
                'bleach','bleaching'}))
            bkgd_type = 'includeBleach';
        elseif any(strcmpi(varargin{ii+1},{'kICSNoise','kICS'}))
            bkgd_type = 'kICS';
        else
            warning(['Unknown option for ''',varargin{ii},...
                ''', using default options.'])
        end
    elseif any(strcmpi(varargin{ii},{'bkgdRgn','noiseRgn'}))
        if isnumeric(varargin{ii+1}) && isequal(size(varargin{ii+1}),[2,2])
            bkgd_rgn = varargin{ii+1};
        else
            warning(['Unknown option for ''',varargin{ii},...
                ''', using default options; input as [xmin,xmax;ymin,ymax].'])
        end
    end
end

if isempty(bkgd_rgn)
    J_bkgd = J;
else
    x_vec = bkgd_rgn(1,1):bkgd_rgn(1,2);
    y_vec = bkgd_rgn(2,1):bkgd_rgn(2,2);
    
    J_bkgd = J(y_vec,x_vec,:);
end

switch bkgd_type
    case 'getMean'
        bkgd_noise = mean(J_bkgd(:));
        %         J_sub = J - bkgd_noise;
    case 'includeBleach'
        size_y = size(J,1);
        size_x = size(J,2);
        n_frames = size(J,3);
        
        t = 0:n_frames-1;
        p = bleachFit(J_bkgd,varargin{:});
        
        if length(p) == 2 % fit was done w/o offset
            bleach_profile = @(p,t) p(1)*exp(-p(2)*t);
        else
            bleach_profile = @(p,t) p(1)*exp(-p(2)*t) + p(3);
        end
        
        bkgd_noise = repmat(reshape(bleach_profile(p,t),[1,1,n_frames]),...
            [size_y,size_x,1]);
        %         J_sub = J - bleach_fit;
    case 'kICS'
        bkgd_noise = numel(J(:,:,1))*var(J_bkgd(:));
end

if strcmp(bkgd_type,'includeBleach')
    figure()
    hold on
    
    int_trace = squeeze(mean(mean(J_bkgd,1),2));
    plot(t,int_trace)
    plot(t,bleach_profile(p,t))
end
