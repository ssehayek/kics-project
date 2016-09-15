%% Input Variables

clear all

addpath(genpath('C:\Users\SimonS\Dropbox (Personal)\Research\MSc\SOFI'));

% Load Data
load_data = 1;
T = 2048; % num of frames to load in

% Fourier filtering
fft_filt_bool = 1;

% Multiplicative Window
window_bool = 0;
win_x_div = 8;
win_y_div = win_x_div;

% Normalization
normByLag = 0;

% kICS Corr Plot
plotTauLags = [1,5,10,20]; % actual lag values i.e. tau=0 is 0th lag
kSqMin = 0.01;
kSqMax = 5;

%% Load Data

if load_data
    suggestedPath = 'C:\Users\SimonS\SharePoint\Simon Sehayek, Mr\PaulWiseman\Saved Sessions';
    [filename,filepath,~] = uigetfile('C:\Users\SimonS\SharePoint\Simon Sehayek, Mr\PaulWiseman\Saved Sessions\*.mat');
    cd(filepath)
    load(filename)
end

%% Load in Movie

% tRange = session_info.tRange;
rect_pos = session_info.rectPos;
ROI = [rect_pos(1),rect_pos(1)+rect_pos(3);rect_pos(2),rect_pos(2)+rect_pos(4)];

if ~exist('loadedMovie','var')
    loadedMovie = im3read(session_info.fileparts.moviepath,1:T,'ROI',ROI);
end

%% Fourier Filtering

if fft_filt_bool
    J = loadedMovie;
    J_w = fft(loadedMovie,[],3);
    J_w(:,:,1) = 0; % remove temporal DC component
    J_t = ifft(J_w,[],3);
    loadedMovie = J_t;
end
% if fft_filt_bool
%     loadedMovie = loadedMovie - repmat(mean(loadedMovie(:,:,1:T),3),[1,1,T]);
% end
 
%% Run kICS

r_k_norm = kICS3(win_movie,'normByLag',normByLag); % temporal correlation function

size_y = size(r_k_norm,1);
size_x = size(r_k_norm,2);

min_size = min(size_x,size_y);
[r_k_circ,kSqVector] = circular(r_k_norm,floor(min_size/2)-2);

r_k_abs = abs(r_k_circ);

%% kICS Corr Plot

% close all

kSqMinIndex = find(kSqVector >= kSqMin,1,'first')
kSqMaxIndex = find(kSqVector <= kSqMax,1,'last')
kSqVectorSubset = kSqVector(kSqMinIndex:kSqMaxIndex);
kSqSubsetInd = kSqMinIndex:kSqMaxIndex;

plotLegend = cell(1,length(plotTauLags));

figure()
hold on

color = lines(length(plotTauLags));

for tauInd = 1:length(plotTauLags)
    plot(kSqVectorSubset,abs(squeeze(r_k_abs(kSqSubsetInd,plotTauLags(tauInd)+1))),'.','Color',color(tauInd,:),'MarkerSize',16)
    plotLegend{tauInd} = ['$\tau = $~' num2str(plotTauLags(tauInd))];
end

hold off

title('kICS correlation function','interpreter','latex')
xlabel('$|{\bf k}|^2$~(pixels\textsuperscript{-2})','interpreter','latex')
ylabel('$R(|{\bf k}|^2,\tau)$','interpreter','latex')
legend(plotLegend,'interpreter','latex')

tightfig(gcf)