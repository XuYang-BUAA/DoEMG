function [decomp_result] = PD3(sig, mu_len_t, isSilence)
% =========================================================================
%    Decompose sEMG signals using PD-III Algorithm.                       *
%                                                                         *
%    Luca D , C. J . Decomposition of Surface EMG Signals[J].             *
%    Journal of Neurophysiology, 2006, 96(3):1646-1657.                   *
%                                                                         *
% INPUT:                                                                  *
%    sig             -- filtered multi-channel sEMG data(generally 4)     *
%    mu_len_t        -- length of MUAP(window for data progressing)       *
%                       (approximately 15ms)                              *
%    isSilence       -- index of whether it contains related information  *
%                                                                         *
% OUTPUT:                                                                 *
%    decomp_result   -- cell for the decompsition result                  *
%                                                                         *
%  WARNINGS:   none                                                       *
%                                                                         *
%  HISTORY:                                                               *
%                                                                         *
% =========================================================================
% ========== Data information ========================================-----
    data = sig.data;
    % fdata = doFilter10KHz(data);
    fdata = data;
    t0 = sig.t0; dt = sig.dt;
    [data_len, ch_num] = size(fdata);
    % Time series
    t = t0 + dt * (0:data_len-1);
    % Prompt
    if (~isSilence)
        fprintf('[Info]: Numbers of Channels: %d\n', ch_num);
        fprintf('[Info]: Sampling Time: %.2f s\n', data_len * dt);
        fprintf('[Info]: Sampling rate: %d Hz\n', 1/dt);
        fprintf('[Info]: Data Length: %d\n', data_len);
    end
% ========== End of Data information ======================================


% ========== Settings =====================================================
    max_ipi_t = 0.35;
    max_ipi = round(max_ipi_t/sig.dt);
    min_firing_times = 10;
% ========== End of Settings ==============================================


% ========== sEMG decomposition ===========================================
% Peak detection
    % Get active threshold & baseline noise sd for each channel.
    [thresh, noise_sigma] = getthreshold(fdata(1:15000,:), 4);
    % Get appropriate spike width.
    spike_width = goodwidth(mu_len_t/sig.dt);
    % Detect peaks whose amplitudes exceed the threshold.
    [spikes, spike_times] = getspikes(fdata, thresh, spike_width,'Plot', 1);
% Generate MU templates

    mu_tmplts = genMUTemplates(spikes, spike_times, max_ipi, noise_sigma);
    %% test template
    mu_len = mu_tmplts.mu_len
    for i=1:25
        subplot(5,5,i)
        for j=1:4
            plot(mu_tmplts.shape(((j-1)*mu_len+1):j*mu_len,i));
            hold on
        end
    end
    if (mu_tmplts.mu_num == 0)
        decomp_result = [];
        return;
    end
    mu_tmplts = mergeMU(mu_tmplts, 0.7);
    
% Generate firing trains
    [spikes, spike_times] = getspikes(fdata,thresh,spike_width,'MergeMode','All');
    [result, resdata] = getFirings(fdata, spikes, spike_times, mu_tmplts);
% ========== End of sEMG decomposition ====================================

    decomp_result.firings = result;
    decomp_result.residue = resdata;
    decomp_result.tmplts  = mu_tmplts;
    
% =========================================================================


end

% ==============================EOF========================================
