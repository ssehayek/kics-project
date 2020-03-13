function [r_k_0_sub] = kICSSubNoise(r_k,ksq_min,ksq_max)

[~,~,noise_inds] = getKSqVector(r_k,'kSqMin',ksq_min,'kSqMax',ksq_max);
r_k_0 = r_k(:,:,1);
r_k_0_sub = r_k(:,:,1) - mean(r_k_0(noise_inds));