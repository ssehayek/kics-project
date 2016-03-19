function [params, Err]= xy_Slide(params0, kSq, tauVector, kICSCorrSubset, sN, lb, ub)
%% slide along Var dimension
         Var = {'w0','etaP','D','k_on','k_off','fd'};% order matters! this order is determined from parameter vector space testing result.
        VarI = [4,    6,     1,  2,     3,     5];

for i0 = 1:6
        eval(['Varlb = lb(',num2str(VarI(i0)),');'])
        eval(['Varub = ub(',num2str(VarI(i0)),');'])
        
    
        
        ErrWalk = zeros(1, sN);
        VarSearch = linspace(Varlb,Varub, sN);

        for i1 = 1 : sN;
            paramsWalk = params0;
            eval(['paramsWalk.', Var{i0}, ' = VarSearch(i1);'])
            ErrWalk(i1) = sum(sum(abs(kICSNorm_xy(paramsWalk, kSq, tauVector) - kICSCorrSubset)));
        end
        
        Ind = find(ErrWalk == min(ErrWalk));
        eval(['params0.',Var{i0},' = VarSearch(Ind(1));'])
        % suggested lb and ub for next iteration.
        a = VarSearch(min(Ind+3, sN));
        b = VarSearch(max(Ind-3, 1));
        minX = min(a, b);
        maxX = max(a, b);
        
        eval(['lb(',num2str(VarI(i0)),') = minX(1);']);
        eval(['ub(',num2str(VarI(i0)),') = maxX(1);']);
        
end

    params = params0;
    Err = ErrWalk(Ind(1));
end