function [aggregates,dye_positions,N_agg_tot] = generateAggPositions(...
    agg_positions,mean_agg_num,std_agg_dist,varargin)

% number of aggregates
N = size(agg_positions,1);
% number of dyes per aggregate (not including initial dye)
n_dyes = poissrnd(mean_agg_num,[1,N]);
% cumulative dye count over indexed vesicles (not including initial dye)
dye_csum = [0,cumsum(n_dyes)];
% number of aggregated dyes (not including initial dye)
N_agg = dye_csum(end);
% relative dye positions to respective aggregate centers
dye_rel_pos = std_agg_dist*randn([N_agg,2]);

% absolute dye positions
aggregates.position = cell(1,N);
for m = 1:N
    m1 = dye_csum(m) + 1;
    m2 = dye_csum(m+1);
    %
    aggregates.position{m} = [agg_positions(m,:);agg_positions(m,:) + ...
        dye_rel_pos(m1:m2,:)];
end

% concatenated positions
dye_positions = cell2mat(aggregates.position');
% total number of dyes per aggregate
aggregates.nDyes = n_dyes + 1;
% total number of aggregated dyes
N_agg_tot = N_agg + N;