function output = getFittingScatterSlide_xy(sim_info,kSq,tauVector,kICSCorrSubset,itNum,ub,lb,ScatterN)
% scatter vectors in parameter space
Pool = rand(ScatterN,6);

for i0 = 1:6
Pool(:,i0) = Pool(:,i0).*(ub(i0)-lb(i0)) + ones(ScatterN,1).*lb(i0);
end
% calculate err for each parameter vector
for i0 = 1:ScatterN
    disp(['sliding parameter vector #',num2str(ScatterN-i0),'/',num2str(ScatterN),';'])
    params_xy.D=Pool(i0,1); % diffusion constant
    params_xy.k_on=Pool(i0,2); % k_on;
    params_xy.k_off=Pool(i0,3); % k_off;
    params_xy.w0=Pool(i0,4); % w0 of PSF;
    params_xy.fd=Pool(i0,5); % diffusion fraction;
    params_xy.etaP=Pool(i0,6); % noise factor (\eta' in doc);

    % slide this parameter vector.
    params0 = params_xy;
    
    for i1 = 1:itNum
        [params, Err]= xy_Slide(params0, kSq, tauVector, kICSCorrSubset, 100, lb, ub);
        params0 = params;
    end
    
    output.ErrSlide(i0)=Err;
    PoolSlid(i0,:)=[params.D, params.k_on, params.k_off, params.w0, params.fd, params.etaP];
end

output.GTerr = sum(sum(abs(kICSNorm_xy(sim_info,kSq,tauVector)-kICSCorrSubset)));
Ind = find(output.ErrSlide == min(output.ErrSlide));
Ind = Ind(1);

params.D=PoolSlid(Ind,1); % diffusion constant
params.k_on=PoolSlid(Ind,2); % k_on;
params.k_off=PoolSlid(Ind,3); % k_off;
params.w0=PoolSlid(Ind,4); % w0 of PSF;
params.fd=PoolSlid(Ind,5); % diffusion fraction;
params.etaP=PoolSlid(Ind,6); % noise factor (\eta' in doc);

output.params = params;
output.paramsTrue = [sim_info.D, sim_info.k_on, sim_info.k_off, sim_info.w0,sim_info.fd,sim_info.etaP];
output.paramsFit  = [params.D, params.k_on, params.k_off, params.w0,params.fd,params.etaP];
output.Pool = PoolSlid;
output.Err = output.ErrSlide;

output.CalcErr = min(output.ErrSlide);
end