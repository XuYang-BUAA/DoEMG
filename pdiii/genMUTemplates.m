function [ mu_tmplts ] = genMUTemplates(spikes, spike_times, max_ipi, noise_sigma)
% =========================================================================
%              Generate MUAP templates from raw sEMG data.                *
%                                                                         *
%                                                                         *
%  INPUT:                                                                 *
%    spikes          -- *
%    spike_times     -- *
%    max_ipi         -- *
%    noise_sigma     -- *
%                                                                         *
%  OUTPUT:                                                                *
%    mu_tmplts       -- *
%                                                                         *
%  WARNINGS:   none                                                       *
%                                                                         *
%  HISTORY:                                                               *
%    7/1/2020 : XuY update                                                *
% =========================================================================
%========== Settings and Parameters =======================================
    % Settings
    op_chn_merge = 'Series';
    % Parameters
    p_accept_region.k1sq = 0.5*0.5;
    p_accept_region.k2sq = 0.5;
    p_accept_region.k3sq = 1;
%     drawacceptregion(p_accept_region);
    min_firing_times = 10;
    [spike_width, spike_num, ch_num] = size(spikes);
%==========================================================================


%========== MU template Extraction ========================================
    mu_num = 0;
    template_s = zeros(spike_width * ch_num, mu_num);
    tmplt_mem = zeros(spike_width * ch_num, spike_num);
    tmplt_muidx = zeros(1, spike_num);
    tmplt_cnt = 0;
    lastfirings = zeros(1,mu_num);
    firecnt = zeros(1,mu_num);
    ipi_mius = zeros(1,mu_num);
    ipi_sigmas = zeros(1,mu_num);
    ipi_nums = zeros(1,mu_num);
% Process of finding MU Tamplates.
    for spike_ind = 1:spike_num
        candi = permute(spikes(:,spike_ind,:),[1 3 2]);
        tn = spike_times(spike_ind);
         
        % Check if candi is a new template.
        switch op_chn_merge
            case 'Series'
                % Connect the muti-channels end-to-end.
                candi_s = candi(:);
                 [k1sq, k2sq, k3sq] = caldiff(candi_s, template_s);
                flag_mu_match = all([k1sq < p_accept_region.k1sq; 
                                     k2sq < p_accept_region.k2sq;
                                     k3sq < p_accept_region.k3sq]);
            case 'Logic'
                % Use OR logic across channels.
                % candi would be a new MU, if there was at least one 
                % channel that satisfied the criteria of a new template.
                k1sq = zeros(ch_num, mu_num);
                k2sq = k1sq; k3sq = k1sq; alpha = k1sq;
                for i = 1:mu_num
                    template_p = reshape(template_s(:,i), [spike_width, ch_num]);
                    [k1sqrx, k2sqrx, k3sqrx, alphax] = caldiff(candi, template_p);
                    k1sq(:,i) = diag(k1sqrx);
                    k2sq(:,i) = diag(k2sqrx);
                    k3sq(:,i) = diag(k3sqrx); 
                    alpha(:,i) = diag(alphax);
                end
                flag_mu_match = all([all(k1sq < p_accept_region.k1sq);
                                     all(k2sq < p_accept_region.k2sq);
                                     all(k3sq < p_accept_region.k3sq)]);
            otherwise
                error('Options for "ChannelMerge": "Series"; "Logic".');
        end

        % Generate new template / Assign candi to a particular MU.
        if any(flag_mu_match)
            % candi is not a new template.
            % MU assignment
            % Calculate MAP
            Pui = caloccurprob(tn - lastfirings, ipi_mius, ipi_sigmas);
            Pui = Pui/sum(Pui);
            Pcandi = sum((repmat(candi_s, 1, mu_num) - template_s).^2);
            %Among these possible mu templates, the one for which the cost
            %has the least value is detected as firing MU.
            cost = Pcandi - 2 * mean(noise_sigma)^2 * log(Pui);
            mu_ind = find(cost == min(cost(flag_mu_match)));
            % Update MU
            template_s(:, mu_ind) = (5 * template_s(:, mu_ind) + candi_s) / 6;
            % Store the template data.
            tmplt_cnt = tmplt_cnt + 1;
            tmplt_mem(:,tmplt_cnt) = template_s(:,mu_ind);
            tmplt_muidx(tmplt_cnt) = mu_ind;
            % IPI calculation
            ipi = tn - lastfirings(mu_ind);
            lastfirings(mu_ind) = tn;
            firecnt(mu_ind) = firecnt(mu_ind) + 1;
            if  ipi <= max_ipi
                miu = ipi_mius(mu_ind);
                sigma = ipi_sigmas(mu_ind);
                [ipi_mius(mu_ind), ipi_sigmas(mu_ind)] = ...
                    calestimateIPI(miu, sigma, ipi);
                ipi_nums(mu_ind) = ipi_nums(mu_ind) + 1;
            end
        else
            % candi is a new template.
            mu_num = mu_num + 1;
            template_s(:,mu_num) = candi_s;
            lastfirings(mu_num) = tn;
            firecnt(mu_num) = 1;
            ipi_mius(mu_num) = max_ipi / 2;
            ipi_sigmas(mu_num) = ipi_mius(mu_num) / 3;
            ipi_nums(mu_num) = 0;
            % Store the template data.
            tmplt_cnt = tmplt_cnt + 1;
            tmplt_mem(:,tmplt_cnt) = template_s(:,mu_num);
            tmplt_muidx(tmplt_cnt) = mu_num;
        end
        
    end
%==========================================================================


%========== Identify decomposable mu templates ============================
    valid_flag = firecnt > min_firing_times;
    valid_idx = find(valid_flag);
    tmp = [tmplt_muidx; spike_times; tmplt_mem];
%   Validation
%     for i = valid_idx
%         figure;
%         subplot(2,1,1);
%         plot(tmp(3:end,tmp(1,:)==i));
%         subplot(2,1,2);
%         spk = permute(spikes(:,tmp(1,:)==i,:),[1 3 2]);
%         spk = reshape(spk, spike_width*ch_num,[]);
%         plot(spk);
%     end
    mu_num = sum(valid_flag);
    template_s = template_s(:,valid_flag);

    tmp_flag = sum(eq(tmplt_muidx, valid_idx'),1) > 0;
    tmp = tmp(:,tmp_flag);
    
%     lastfirings = lastfirings(:,valid_flag);
    firecnt = firecnt(:,valid_flag);
%     ipi_mius = ipi_mius(:,valid_flag);
%     ipi_sigmas = ipi_sigmas(:,valid_flag);
%     ipi_nums = ipi_nums(:,valid_flag);

%==========================================================================


%========== Construct the MUPool structure ================================
    mu_tmplts = struct('shape', template_s, 'ch_num', ch_num, 'tmplt_mem', tmp,...
        'mu_num', mu_num, 'mu_len', spike_width, 'fire_times',firecnt, 'mu_id', valid_idx);

end
%==============================EOF=========================================
