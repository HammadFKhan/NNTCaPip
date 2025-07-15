function [C_dec,F_detrended] = convolveFluroescence(F,Fneu,conFlag)
p = 2;

[N,T] = size(F);
C_dec = zeros(N,T);         % deconvolved DF/F traces
S_dec = zeros(N,T);         % deconvolved neural activity
bl = zeros(N,1);            % baseline for each trace (should be close to zero since traces are DF/F)
neuron_sn = zeros(N,1);     % noise level at each trace
g = cell(N,1);              % discrete time constants for each trace
if p == 1; model_ar = 'ar1'; elseif p == 2; model_ar = 'ar2'; else; error('This order of dynamics is not supported'); end
options.df_prctile = 5;
options.detrend_window = 2;
options.fr = 30;

dF = [];
F_detrended = detrend_calcium(F,Fneu,options);
options = CNMFSetParms;
options.spk_SNR = 2;
if conFlag
    for i = 1:N
        spkmin = options.spk_SNR*GetSn(F_detrended(i,:));
        lam = choose_lambda(exp(-1/(options.fr*options.decay_time/2)),GetSn(F_detrended(i,:)),options.lam_pr);
        [cc,spk,opts_oasis] = deconvolveCa(F_detrended(i,:),model_ar,'method','foopsi','optimize_pars',true,'maxIter',20,...
            'window',100,'lambda',lam,'smin',spkmin);
        bl(i) = opts_oasis.b;
        C_dec(i,:) = cc(:)' + bl(i);
        S_dec(i,:) = spk(:);
        neuron_sn(i) = opts_oasis.sn;
        g{i} = opts_oasis.pars(:)';
        disp(['Performing deconvolution. Trace ',num2str(i),' out of ',num2str(N),' finished processing.'])
    end

    [corr_coeff, p_value] = corrcoef(F_detrended, C_dec);
    fprintf('r = %.2f\np = %.2e \n', corr_coeff(2), p_value(2));
end

end
function F_dff = detrend_calcium(F, f, options)
    % Default parameters
    if ~isfield(options, 'df_prctile'), options.df_prctile = 20; end  % 20th percentile baseline
    if ~isfield(options, 'detrend_window'), options.detrend_window = 100; end % 60s window (assuming 10Hz)
    
    % Step 1: Compute ΔF/F using percentile baseline
    Fd = prctile(F, options.df_prctile, 2);
    background_est = f;  % Estimate background fluorescence
    F0_bg = prctile(background_est, options.df_prctile, 2);
    F0 = repmat(F0_bg + Fd, 1, size(F, 2));
    F_dff = (F - repmat(Fd, 1, size(F, 2))) ./ F0;
end
