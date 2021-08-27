function organizeKICSSims(basedir,targetfile)

targetfile_m = matfile(targetfile,'Writable',true);

% directory level 1
D = targetfile_m.D;
k_on = targetfile_m.k_on;
k_off = targetfile_m.k_off;
k_p = targetfile_m.k_p;

% directory level 2
prob_place = targetfile_m.prob_place;
mean_agg_num = targetfile_m.mean_agg_num;
std_agg_dist = targetfile_m.std_agg_dist;
num_filaments = targetfile_m.num_filaments;

% filename
N = targetfile_m.N;
N_diff = targetfile_m.N_diff;
if exist('targetfile_m.frac_diff','var')
    frac_diff = targetfile_m.frac_diff;
else
    frac_diff = [];
    targetfile_m.frac_diff = frac_diff;
end
w0 = targetfile_m.w0;
sz = targetfile_m.sz;
T = targetfile_m.T;
nsub = targetfile_m.n_sub_frames;
noise_type = targetfile_m.noise_type;

% this is meant to better organize files & speed up simulations
if num_filaments == 0 || prob_place == 0
    % no filaments means no immobile particles, or aggregates
    num_filaments = 0; prob_place = 0;
    mean_agg_num = 0; std_agg_dist = 0;
    
    targetfile_m.num_filaments = 0;
    targetfile_m.prob_place = 0;
    targetfile_m.mean_agg_num = [];
    targetfile_m.std_agg_dist = [];
end
if mean_agg_num == 0
    % no aggregates
    std_agg_dist = 0;
    
    targetfile_m.std_agg_dist = [];
end
if N_diff == 0
    % no diffusers
    D = 0; frac_diff = 0;
    
    targetfile_m.D = [];
end
if k_off == 0
    % no blinking
    k_on = 1;
    
    targetfile_m.k_on = 1;
end
%

dir_1 = ['D_',num2str(D),'_kon_',num2str(k_on),'_koff_',...
    num2str(k_off),'_kp_',num2str(k_p)];

dir_2 = ['n_fils_',num2str(num_filaments),'_prob_place_',...
    num2str(prob_place),'_mean_agg_num_',num2str(mean_agg_num),...
    '_std_agg_dist_',num2str(std_agg_dist)];

if isempty(frac_diff)
    filename = ['N_',num2str(N),'_N_diff_',num2str(N_diff),'_w0_',...
        num2str(w0),'_sz_',num2str(sz),'_T_',num2str(T),'_nsub_',...
        num2str(nsub),'_noise_type_',num2str(noise_type),'.mat'];
else
    filename = ['N_',num2str(N),'_frac_diff_',num2str(frac_diff),'_w0_',...
        num2str(w0),'_sz_',num2str(sz),'_T_',num2str(T),'_nsub_',...
        num2str(nsub),'_noise_type_',num2str(noise_type),'.mat'];
end

K = k_on + k_off;
targetfile_m.true_params = [D,k_on/K,K,N_diff/N];

filedir = [basedir,filesep,dir_1,filesep,dir_2];
filepath = [filedir,filesep,filename];

if ~exist(filedir,'dir')
    mkdir(filedir)
end

if ~strcmp(targetfile_m.Properties.Source,filepath)
    if exist(filepath,'file')
        filepath = iterateFilename(filepath);
    end
    movefile(targetfile_m.Properties.Source,filepath);
end