parpool;
num_datasets = 9600;

outdir = 'C:\\Users\\TSB-MTO\\Desktop\\Tra-Ap2025\\revision2\\IGSCA';
%warningfile = fullfile(outdir, 'warning_iterations.mat');

%warning_iters = nan(num_datasets,1);  % preallocate (NaN = no warning)

parfor i = 8001:num_datasets
        fname = sprintf('C:\\Users\\TSB-MTO\\Desktop\\Tra-Ap2025\\revision2\\matlab\\DATA\\data%d.mat', i);
        S = load(fname);

        X = S.list.X;
        nfactors = size(S.list.Ptrue, 2);

        J = size(X,2);
        P = nfactors;

        W0 = ones(J,P)*99;
        C0 = W0';
        B0 = zeros(P,P);
        dimtype = ones(1,P);

        Data = S.list.X;

        % Regularization settings
        lnorm = 1;
        num_folder = 5;
        LC = [0,10.^(-3:0.2:0)];
        LB = 0;
        nbt = 100;

    % Reset warnings
    lastwarn('');
    flagged = 0;

    % Run function
    R = igsca_reg_v3(Data,W0,C0,B0,dimtype,lnorm,1,num_folder,LC,LB,nbt);

    % Check if warning occurred
    [warnMsg, ~] = lastwarn;
    if ~isempty(warnMsg)
        flagged = 1;
    end

    % Save using parsave
    outname = fullfile(outdir, sprintf('IGSCA_%d.mat', i));
    parsave(outname, R, flagged);

end
