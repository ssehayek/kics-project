contents = dir('*.mat');
mat_files = {contents.name};

for i = 1:length(mat_files)
    try
        file_i = matfile(mat_files{i});
        
        % directory level 1
        D_i = file_i.D;
        k_on_i = file_i.k_on;
        k_off_i = file_i.k_off;
        k_p_i = file_i.k_p;
        
        dir_1 = ['D_',num2str(D_i),'_kon_',num2str(k_on_i),'_koff_',...
            num2str(k_off_i),'_kp_',num2str(k_p_i)];
        
        % directory level 2
        mean_agg_num_i = file_i.mean_agg_num;
        std_agg_dist_i = file_i.std_agg_dist;
        
        dir_2 = ['mean_agg_num_',num2str(mean_agg_num_i),'_std_agg_dist_',...
            num2str(std_agg_dist_i)];
        
        % filename
        N_i = file_i.N;
        N_diff_i = file_i.N_diff;
        w0_i = file_i.w0;
        sz_i = file_i.sz;
        T_i = file_i.T;
        nsub_i = file_i.n_sub_frames;
        noise_type_i = file_i.noise_type;
        
        filename = ['N_',num2str(N_i),'_N_diff_',num2str(N_diff_i),'_w0_',...
            num2str(w0_i),'_sz_',num2str(sz_i),'_T_',num2str(T_i),'_nsub_',...
            num2str(nsub_i),'_noise_type_',num2str(noise_type_i),'.mat'];
        
        filedir = [dir_1,filesep,dir_2];
        filepath = [filedir,filesep,filename];
        
        mkdir(filedir)
        movefile(file_i.Properties.Source,filepath);
    end
end