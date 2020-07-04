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
    %thresh = 0.22;
    [spikes, spike_times] = getspikes(fdata, thresh, spike_width,'Plot', 0);
% Generate MU templates
    mu_tmplts = genMUTemplates(spikes, spike_times, max_ipi, noise_sigma);
    %% test template
%     %%test template mem
%     figure()
%     for i=1:mu_tmplts.mu_num
%         subplot(5,5,i)
%         for j=1:length(mu_tmplts.tmplt_mem)
%             if(mu_tmplts.mu_id==mu_tmplts.tmplt_mem(1,j))
%                 plot(mu_tmplts.tmplt_mem(3:77,j));
%                 hold on
%             end
%         end
%     end
            
    if (mu_tmplts.mu_num == 0)
        decomp_result = [];
        return;
    end
    mu_tmplts = mergeMU(mu_tmplts, 0.7);
    %test
    mu_len = mu_tmplts.mu_len;
    t_mu = (t0:mu_len-1)' * dt;
    figure
    for i=1:mu_tmplts.mu_num
        subplot(5,5,i)
        
        for j=1:ch_num
            plot(t_mu,mu_tmplts.shape(((j-1)*mu_len+1):j*mu_len,i));%+(j-1)/2);
            hold on;
            axis([t_mu(1),t_mu(mu_len),-1,1])
        end
    end
    figure()
    selecet_mu = 4;
    test_index = find(mu_tmplts.tmplt_mem(1,:)==mu_tmplts.mu_id(selecet_mu));
    for j=1:4
        subplot(2,2,j)
        for i=1:mu_tmplts.fire_times(selecet_mu)
        %for i =1:5
            plot(mu_tmplts.tmplt_mem(3+mu_len*(j-1):mu_len*(j-1)+mu_len+2,test_index(i)));
            hold on
            plot(mu_tmplts.shape(1+mu_len*(j-1):mu_len*j,selecet_mu),'*')
        end
    end
%     figure()
%     w=hanning(75);
%     for j=1:4
%         subplot(2,2,j)
%         for i=1:10
%             plot(w.*mu_tmplts.tmplt_mem(3+75*(j-1):75*(j-1)+77,test_index(3*i)));
%             hold on
%         end
%     end
%     figure
%     for i=1:4
%         subplot(2,2,i)
%         plot(t_mu,mu_tmplts.shape(1+75*(i-1):75*i,16))
%         hold on
%     end
% Generate firing trains
%     [spikes, spike_times] = getspikes(fdata,thresh,spike_width,'MergeMode','All');
%     [result, resdata] = getFiringsNew(fdata, spikes, spike_times, mu_tmplts);
    [result, resdata] = getFiringsNew(fdata, mu_tmplts,thresh,spike_width);
% ========== End of sEMG decomposition ====================================

    decomp_result.firings = result;
    decomp_result.residue = resdata;
    decomp_result.tmplts  = mu_tmplts;
    
% =========================================================================


end

% ==============================EOF========================================
