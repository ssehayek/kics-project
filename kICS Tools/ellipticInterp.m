% Written by: Xiyu Yi, Simon Sehayek
%
% Email: chinahelenyxy@gmail.com, simon.sehayek@gmail.com
%
% this code is meant to Fourier interpolate the kICS function along
% concentric circles in |k| centered on the lowest frequency value.
%
% INPUT
%
% A: normalized kICS function at single lag value, which has lowest
% fequency at center of matrix
%
% ksq_vec: |k|^2 query points
%
% n_theta: number of theta points to sample for each value of |k| in the
% interval (0,2*pi]. "n_theta" has same dimensions as "ksq_vec"
function [B] = ellipticInterp(A,ksq_vec,n_theta,varargin)

[ydim,xdim]=size(A); % note inverted order definition of x and y

O = [xdim,ydim]/2; % coordinates of correlation center

% rectangle defining boundaries of image
bound_x = [1,1,xdim,xdim];
bound_y = [1,ydim,ydim,1];
image_bound = @(xq,yq) inpolygon(xq,yq,bound_x,bound_y);
%

FA = fft2(A); 

%% now make sure inverse Fourier transform works:
% take inverse-Fourier transform of A, manually. don't use fast fourier
% transform.
Fbx=zeros(xdim,1);Fbx(2)=1;Fx=fft(Fbx);
Fby=zeros(1,ydim);Fby(2)=1;Fy=fft(Fby);

cFy=conj(Fy);
cFx=conj(Fx);
%

B=cell(1,length(ksq_vec));
for i0=1:length(ksq_vec)
    kmag=sqrt(ksq_vec(i0));
    
    % scale factors used to define ellipse in pixel space
    a = xdim*kmag/(2*pi); 
    b = ydim*kmag/(2*pi);
    
    % for i0th rho value, take theataS1 intervals between 0 to 2pi for
    % theta
    thetaSl=[2*pi/n_theta(i0):2*pi/n_theta(i0):2*pi];
    x=zeros(1,length(thetaSl));
    for i1=1:length(thetaSl)
        thetaS=thetaSl(i1);
        
        % convert |k| coordinate to pixels
        % for fixed |k|, the corresponding pixel as a function of theta
        % lies along an ellipse 
        n = a*b/sqrt((b*cos(thetaS))^2+(a*sin(thetaS))^2);
        ifxp=n*cos(thetaS)+O(1); % pixel x-coord
        ifyp=n*sin(thetaS)+O(2); % y-coord
        %
        % check if (ifxp,ifyp) is within confines of image
        if image_bound(ifxp,ifyp)
            x(i1)=(cFy.^ifyp)*(FA)*(cFx.^ifxp)./ydim./xdim;
        else
            x(i1) = NaN;
        end
        %
    end
    x(isnan(x)) = [];
    B{i0}=x;
end
