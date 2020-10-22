%% preliminary code

% reseed the random number generator; otherwise same random
% numbers are generated each time MATLAB restarts
seed = rng('shuffle');

run('analysisInput.m') % gets all parameters from "analysisInput.m"

addpath(genpath(code_path))

%% create save path for simulation

% this is meant to better organize files & speed up simulations
%
if (num_filaments == 0 || prob_place == 0) && (use_mask == 0)
    % no filaments means no immobile particles, or aggregates
    num_filaments = 0; prob_place = 0;
    mean_agg_num = 0; std_agg_dist = 0;
end
if use_mask
    num_filaments = 0; prob_place = 0;
end
if mean_agg_num == 0
    % no aggregates
    std_agg_dist = 0;
end
if N_diff == 0
    % no diffusers
    D = 0; frac_diff = 0;
end
if k_off == 0
    % no blinking
    k_on = 1;
end
%

dir_1 = ['D_',num2str(D),'_kon_',num2str(k_on),'_koff_',...
    num2str(k_off),'_kp_',num2str(k_p)];

if use_mask
    dir_2 = ['N_imm_',num2str(N_imm),'_mean_agg_num_',num2str(mean_agg_num),...
        '_std_agg_dist_',num2str(std_agg_dist),'_mask'];
else
    dir_2 = ['n_fils_',num2str(num_filaments),'_prob_place_',...
        num2str(prob_place),'_mean_agg_num_',num2str(mean_agg_num),...
        '_std_agg_dist_',num2str(std_agg_dist)];
end

if isempty(frac_diff)
    dir_3 = ['N_diff_',num2str(N_diff),'_w0_',...
        num2str(w0),'_sz_',num2str(sz),'_T_',num2str(T),'_nsub_',...
        num2str(nsub),'_noise_type_',num2str(noise_type)];
else
    dir_3 = ['frac_diff_',num2str(frac_diff),'_w0_',...
        num2str(w0),'_sz_',num2str(sz),'_T_',num2str(T),'_nsub_',...
        num2str(n_sub_frames),'_noise_type_',num2str(noise_type)];
end

if isempty(sim_tag)
    sim_tag = 'simulation';
end

save_simpath = iterateFilename([base_dir,filesep,'simulations',filesep...
    dir_1,filesep,dir_2,filesep,dir_3,filesep],'prepStr','_rep_');
save_simfile = [save_simpath,filesep,sim_tag,'.mat'];

if ~exist(save_simpath,'dir')
    mkdir(save_simpath)
end

disp(['creating simulation ',save_simpath,'.'])

%% create run directory

cd(save_simpath)

% move all files from temp_dir to save_simpath
movefile([tmp_dir,filesep,'*'],save_simpath);
% delete unnecessary files
delete submitMatlabJob condor_script
% rename analysisInput.m -> parameters.m
movefile('analysisInput.m','parameters.m');

%% create simulations

% struct for mask option
mask_struct.use_mask = use_mask;
mask_struct.mask_filepath = mask_filepath;
mask_struct.N_imm = N_imm;

[J,true_params] = kicsSim(sz,T,w0,N_diff,D,k_on,k_off,k_p,...
    mean_agg_num,std_agg_dist,num_filaments,prob_place,'parallel',sim_parallel,...
    'subFrames',n_sub_frames,'savepath',save_simfile,'noiseType',...
    noise_type,'kerVar',kernel_varargin,'laserVar',laser_varargin,...
    'noiseVar',noise_varargin,'fracDiff',frac_diff,'mask',mask_struct);