function [ firings, resdata ] = getFirings(rdata, spikes, spike_train, mus)
% =========================================================================
%                          Get MU firing trains.                          *
%                                                                         *
%  INPUT:                                                                 *
%    rdata           -- sEMG data                                         *
%    spikes          -- firing segments                                   *
%    spike_train     -- timings of firing segments                        *
%    mus             -- mu templates                                      *
%                    -- members:                                          *
%                    -- shape      :template signal segments              *
%                    -- ch_num     :channel number                        *
%                    -- tmplt_mem  :memory of all templates               *
%                    -- mu_num     :MU number                             *
%                    -- mu_len     :length of MU segment                  *
%                    -- fire_times :firing times of each MU               *
%                    -- mu_id      :id index of MU templates              *
%                                                                         *
%  OUTPUT:                                                                *
%    firings         -- merged MU templates                               *
%    resdata         -- residual data                                     *
%                                                                         *
%  WARNINGS:   none                                                       *
%                                                                         *
%  HISTORY:                                                               *
%    7/2/2020 : XuY update                                                *
% =========================================================================
% =========================================================================
% =========== Basic information & Settings=================================
% =========================================================================
    data_len = length(rdata);
    resdata = rdata;
    
    mu_num = mus.mu_num;
    mu_len = mus.mu_len;
    mu_lenH = floor(mu_len/2);
    ch_num = mus.ch_num;
    mu_shapes = mus.shape;
    mu_ids = mus.mu_id;
    tmplt_mem = mus.tmplt_mem;
    
    asr_match_thres = 0.2;
    
    [~, spk_num, ~] = size(spikes);
    spikes = reshape(permute(spikes,[1 3 2]), [], spk_num);
    
% =========================================================================
% =========== Calculate time-varied templates & ASR =======================
% =========================================================================
    asr = zeros(data_len, mu_num);
    tmplts = mu_shapes;
    for spk_idx = 1:spk_num
        % Calculate time-varied templates 
        t = spike_train(spk_idx);
        mempointer = tmplt_mem(1,tmplt_mem(2,:) == t);
        if mempointer
            mu_idx = find(mu_ids==mempointer);
            tmplts(:,mu_idx) = (5*tmplts(:,mu_idx)+spikes(:,spk_idx))/6;
            tmplt_mem(3:end,spk_idx) = tmplts(:,mu_idx);
        end

        % Calculate ASR
        tmp = calASR(tmplts, spikes(:,spk_idx));
        asr(t,:) = tmp';
    end
    
    % Find local maxima in asr.
%     for mu_idx = 1:mu_num
%         [peaks, timings] = findpeaks(asr(:,mu_idx), ...
%             'MinPeakHeight',asr_match_thres,'MinPeakDistance', mu_len);
%         asr(:,mu_idx) = 0;
%         asr(timings, mu_idx) = peaks;
%     end  

% =========================================================================
% =========== Aliasing rejection analysis =================================
% =========================================================================
%{
    tmplts = mu_shapes;
    lobe_loc = zeros(5, ch_num, mu_num);
    lobe_amp = zeros(5, ch_num, mu_num);
    for idx = 1:mu_num
        [lobe_loc(:,:,idx),lobe_amp(:,:,idx)] = getlobe(tmplts(:,idx), ch_num, 5);
    end
    for spk_idx = 1:spk_num
        % Calculate time-varied templates.
        t = spike_train(spk_idx);
        mempointer = tmplt_mem(1,tmplt_mem(2,:)==t);
        if mempointer
            mu_idx = find(mu_ids==mempointer);
            
            tmplts(:,mu_idx) = (5*tmplts(:,mu_idx)+spikes(:,spk_idx))/6;
            [lobe_loc(:,:,mu_idx),lobe_amp(:,:,mu_idx)] = getlobe(tmplts(:,mu_idx), ch_num, 5);
            
            tmplt_mem(3:end,spk_idx) = tmplts(:,mu_idx);
        end
        
        % if no actived mu, go to next time.
        if all(asr(t,:) == 0)
            continue;
        end
        
        % Active mu index
        mu_idx_act = find(asr(t,:) > 0);
        
        % Analysis 'main lobe'
        p1 = 0.5; p2 = 0.1;
        [~, mu_num_act] = size(mu_idx_act);

        % if there are MORE than 1 actived mu, do 'main lobe' analysis.
        if mu_num_act > 1
            % Check all combinations
            ac_mu_idx = [];
            ac_mu_asr = 0;
            for i = 1:mu_num_act
                cmb = nchoosek(mu_idx_act, i);
                [cmb_num,~] = size(cmb);
                for cmb_idx = 1:cmb_num
                    chk_mu_idxs = cmb(cmb_idx,:);
                    mlsum = sum(lobe_amp(3,:,chk_mu_idxs), 3);
                    D = spikes(:,spk_idx);
                    lowlim = (1 - p2) ./ p1 .* D;
                    uplim = (1 + p2) ./ p1 .* D;
                    if any(and(lowlim<=mlsum, mlsum<=uplim))
                    % if peak height sum is acceptable.
                        if sum(asr(t,chk_mu_idxs)) > ac_mu_asr
                        % if asr sum is the biggest, record it.
                            ac_mu_idx = chk_mu_idxs;
                            ac_mu_asr = sum(asr(t,chk_mu_idxs));
                        end
                    end
                end
            end
            
            if ac_mu_asr == 0
            % if no acceptable combinations, select the one having the largest
            % asr value.
                [~,ac_mu_idx] = max(asr(t,:));
            end
            % End check all combinations
            
            % Reject the aliasings
            temp = asr(t,ac_mu_idx);
            asr(t,:) = 0; asr(t,ac_mu_idx) = temp;
        end
        
    end
%}

% =========================================================================
% =========== Probability Assignment ======================================
% =========================================================================
%
    max_hypo_num = 2;
    hypos = zeros(spk_num, mu_num, max_hypo_num);
    ida_probs = zeros(spk_num, mu_num);
    tmplts = mu_shapes;
%     mu_pthres = calPHat(tmplts, tmplts);
%     mu_pthres(1:mu_num+1:end) = 0;
%     mu_pthres = max(mu_pthres, [], 2);
    for spk_idx = 1:spk_num
        % Calculate time-varied templates.
        t = spike_train(spk_idx);
        mempointer = tmplt_mem(1,tmplt_mem(2,:)==t);
        if mempointer
            mu_idx = find(mu_ids==mempointer);
            tmplts(:,mu_idx) = (5*tmplts(:,mu_idx)+spikes(:,spk_idx))/6;
            tmplt_mem(3:end,spk_idx) = tmplts(:,mu_idx);
        end
        
        res = spikes(:,spk_idx);
        asr_l = asr(t,:);
        
        [probs, pwrs] = reIDA(tmplts, asr_l>asr_match_thres, res, 0);
        [hypo_num, ~] = size(probs);
        for hypo_idx = 1:min(hypo_num, max_hypo_num)
            hypos(spk_idx, :, hypo_idx) = probs(hypo_idx,:);
        end
        if (hypo_num > 0)
            ida_probs(spk_idx, :) = max(probs, [],1);
        end
    end
%}

% =========================================================================
% =========== Decision Making =============================================
% =========================================================================
    firings = zeros(data_len, mu_num);
    for mu_idx = 1:mu_num
        spk_idxs = ida_probs(:,mu_idx) > 0;
        f_tm = getMaxUti(spike_train(spk_idxs), ida_probs(spk_idxs,mu_idx), mu_len);
        firings(f_tm, mu_idx) = 1;
    end

% ABANDON CODE SEGMENT.
% Method:
% for every non-overlapping data segment, extend the templates to the
% segment length and using least square to minimize the RMSE.
%{
% =========================================================================
% =========== Decision Making =============================================
% =========================================================================
    dis = diff([spike_train(1)-mu_len,spike_train,spike_train(end)+mu_len]);
    dis = dis >= mu_len;
    s = spike_train(dis(1:end-1)) - mu_lenH;
    e = spike_train(dis(2:end)) + mu_lenH;
    dataseg = [s',e'];

    tmplts = mu_shapes;
    [seg_num,~] = size(dataseg);
    alphas = cell(seg_num,1);
    for sig_idx = 1:seg_num
        st = dataseg(sig_idx,1); et = dataseg(sig_idx,2); seg_len = et-st+1;
        d = rdata(st:et,:);
        local_asr = asr(st:et,:);
        
        while(true)
            gdts = find(local_asr>0);
            t = zeros(seg_len*ch_num, length(gdts));  
            tmplts = reshape(tmplts, mu_len, ch_num, mu_num);
            for gdt_idx = 1:length(gdts)
                gdt = gdts(gdt_idx);
                mu_idx = ceil(gdt / seg_len);
                dt = mod(gdt, seg_len);
                train = zeros(seg_len,1); train(dt) = 1;
                tmp = zeros(seg_len, ch_num);
                for chn_idx = 1:ch_num
                    tmp(:,chn_idx) = conv(train, tmplts(:,chn_idx,mu_idx), 'same');
                end
                t(:,gdt_idx) = tmp(:);
            end
            alpha = t \ d(:);
            if all(alpha>=0)
                break;
            end
            local_asr(gdts(alpha<0)) = 0;
        end
    end
%}

% ABANDON CODE SEGMENT.
% Method:
% Try all possible combinations to find which comb can constitute the power
% of the data.
%{
% =========================================================================
% =========== Generate hypothesis =========================================
% =========================================================================
    max_hypo_num = 3;
    hypo = zeros(spk_num, mu_num, max_hypo_num);
    tmplts = mu_shapes;
    for spk_idx = 1:spk_num
        % Calculate time-varied templates
        mempointer = tmplt_mem(1,tmplt_mem(2,:)==spike_train(spk_idx));
        if mempointer
            mu_idx = find(mu_ids==mempointer);
            tmplts(:,mu_idx) = (5*tmplts(:,mu_idx)+spikes(:,spk_idx))/6;
            tmplt_mem(3:end,spk_idx) = tmplts(:,mu_idx);
        end
        
        spk_pwr = sum(spikes(:,spk_idx).^2);
        ac_mu_idx = find(asr(spk_idx,:)>0);
        ac_mu_num = length(ac_mu_idx);
        hypo_idx = 0;
        for i = 1:ac_mu_num
            cmb = nchoosek(ac_mu_idx, i);
            [cmb_num,~] = size(cmb);
            for cmb_idx = 1:cmb_num
                res = spikes(:,spk_idx) - sum(tmplts(:,cmb(cmb_idx,:)),2);
                if sum(res.^2)/spk_pwr < 0.5 && hypo_idx < max_hypo_num
                    hypo_idx = hypo_idx + 1;
                    hypo(spk_idx,cmb(cmb_idx,:),hypo_idx) = 1;
                end
            end
        end
    end
%}

% =========================================================================
% =========== Data Segmentation ===========================================
% =========================================================================
% Data Segment: Get the data segment where MU fires.
    tmp = find(sum(firings,2)');
    dis = diff([tmp(1)-mu_len,tmp,tmp(end)+mu_len]);
    dis = dis >= mu_len;
    s = tmp(dis(1:end-1)) - mu_lenH;
    e = tmp(dis(2:end)) + mu_lenH;
    s(1) = max([s(1), 1]);
    e(end) = min([e(end), data_len]);
    dataseg = [s',e'];
% =========================================================================
% =========== Results Refinement ==========================================
% =========================================================================
    tmplts = mu_shapes;
    [seg_num,~] = size(dataseg);
    for sig_idx = 1:seg_num
        st = dataseg(sig_idx,1); et = dataseg(sig_idx,2); seg_len = et-st+1;
        d = resdata(st:et,:);
        local_f = firings(st:et,:);
        
        gdts = find(local_f>0);
        t = zeros(seg_len*ch_num, length(gdts));
        
        tmplts = reshape(tmplts, mu_len, ch_num, mu_num);
        
        for gdt_idx = 1:length(gdts)
            gdt = gdts(gdt_idx);
            mu_idx = ceil(gdt / seg_len);
            dt = mod(gdt, seg_len);
            train = zeros(seg_len,1); train(dt) = 1;
            tmp = zeros(seg_len, ch_num);
            for chn_idx = 1:ch_num
                tmp(:,chn_idx) = conv(train, tmplts(:,chn_idx,mu_idx), 'same');
            end
            t(:,gdt_idx) = tmp(:);
        end
        alpha = t \ d(:);
        resdata(st:et,:) = d - reshape(t * alpha,[],ch_num);
    end
end

%========================================================================--
function [Probs, Amps] = calPHat(t, d)
% Calculate the estimated probability.
% If p(Mp x Np) or d(Md x Nd) is matrix, calculate on column, 
% and the result is Np x Nd.

    [~,~,~,Amps] = caldiff(d, t);Amps = Amps';
    B = Amps;
    B(Amps > 1) = 1./Amps(Amps > 1); B(Amps < 0) = 0; B = B .* abs(Amps);
    Probs = B.*sqrt(sum(t.^2)' * (1./sum(d.^2)));

end

%========================================================================--
function [probs, pwrs] = reIDA(t, mask, d, ord)
% Iterative MUAP Discrimination Analysis.(Recursion)
    [mu_len, mu_num] = size(t);    
    probs  = zeros(0, mu_num);
    pwrs   = zeros(0, 0);
    if all(~mask)
        return;
    end
    
    [p,a] = calPHat(t, d);
    p(~mask') = 0; a(~mask') = 0;
    
    [~, mu_idx_asr] = sort(p, 'descend');
    for mu_idx = mu_idx_asr'
        if ~mask(mu_idx)
            continue;
        end
        mask(mu_idx) = false;
        
        res = d - a(mu_idx) * t(:,mu_idx);
        pwr = sum(res.^2) / sum(d.^2);
        if pwr <= 0.25
            subp = zeros(1, mu_num);
            subp(mu_idx) = (0.5)^(ord)*p(mu_idx);
            probs = [probs;subp];
            pwrs  = [pwrs;  pwr];
            return;
        else
            [subp,subpwr] = reIDA(t, mask, res, ord+1);
            subp(:,mu_idx) = (0.5)^(ord)*p(mu_idx);
            probs  = [probs; subp];
            pwrs   = [pwrs;subpwr];
        end
        
        mask(mu_idx) = true;
    end
    
    if (ord == 0) && (isempty(pwrs))
        mu_idx = mu_idx_asr(1);
        subp = zeros(1, mu_num);
        subp(mu_idx) = (0.5)^(ord)*p(mu_idx);
        
        res = d - a(mu_idx) * t(:,mu_idx);
        pwr = sum(res.^2) / sum(d.^2);
        
        probs = [probs;subp];
        pwrs  = [pwrs;  pwr];
    end

end

%==========================================================================
% function [amp, chn] = getMUAmp(mu_shape)
% % Calculate MU Amplitude.
%     amp_of_chns = max(abs(mu_shape));
%     [amp, chn] = max(amp_of_chns);
% end
% 
% function [ firings, resdata ] = getFirings( raw, spikes, spike_train, mus )
% %Get MU firing trains.
% 
% %=========-- Basic information & Settings=================================-
%     data_len = length(raw);
%     resdata = raw;
%     
%     mu_num = mus.mu_num;
%     mu_len = mus.mu_len;
%     mu_lenH = floor(mu_len/2);
%     ch_num = mus.ch_num;
%     mu_shapes = mus.shape;
%     mu_ids = mus.mu_id;
%     tmplt_mem = mus.tmplt_mem;
%     
%     [~, spk_num, ~] = size(spikes);
%     spikes = reshape(permute(spikes,[1 3 2]), [], spk_num);
% %=========-- End of Basic information & Settings===========================
% 
% %=========-- Calculate time-varied templates & ASR ========================
%     mu_amp = max(abs(mu_shapes));
%     [~, muIdxofAmp] = sort(mu_amp, 'descend');
%     mu_amp_thresh = mu_amp / max(mu_amp);
%     asr = zeros(spk_num, mu_num);
%     tmplts = mu_shapes;
%     for spk_idx = 1:spk_num
%         % Calculate time-varied templates
%         mempointer = tmplt_mem(1,tmplt_mem(2,:)==spike_train(spk_idx));
%         if mempointer
%             mu_idx = find(mu_ids==mempointer);
%             tmplts(:,mu_idx) = (5*tmplts(:,mu_idx)+spikes(:,spk_idx))/6;
%             tmplt_mem(3:end,spk_idx) = tmplts(:,mu_idx);
%         end
%         
%         % Calculate ASR
%         norm2t = sum(tmplts.^2);
%         norm2d = sum(spikes(:,spk_idx).^2);
%         dotmulti = sum(repmat(spikes(:,spk_idx),1,mu_num) .* tmplts);
%         A = dotmulti ./ norm2t;
%         coeff = A;
%         coeff(abs(A-1)<mu_amp_thresh) = 1;
%         coeff(A>(1+mu_amp_thresh)) = 1./coeff(A>(1+mu_amp_thresh));
%         coeff(A<0) = 0;
%         asr(spk_idx,:) = dotmulti ./ sqrt(norm2d.*norm2t) .* coeff;
%     end
% %=========-- End of Calculate time-varied templates & ASR ===============--
%     
% %=========-- Generate hypothesis ==========================================
%     max_hypo_num = 3;
%     hypo = zeros(spk_num, mu_num, max_hypo_num);
%     tmplts = mu_shapes;
%     for spk_idx = 1:spk_num
%         % Calculate time-varied templates
%         mempointer = tmplt_mem(1,tmplt_mem(2,:)==spike_train(spk_idx));
%         if mempointer
%             mu_idx = find(mu_ids==mempointer);
%             tmplts(:,mu_idx) = (5*tmplts(:,mu_idx)+spikes(:,spk_idx))/6;
%             tmplt_mem(3:end,spk_idx) = tmplts(:,mu_idx);
%         end
%         
%         spk_pwr = sum(spikes(:,spk_idx).^2);
%         ac_mu_idx = find(asr(spk_idx,:)>0);
%         ac_mu_num = length(ac_mu_idx);
%         hypo_idx = 0;
%         for i = 1:ac_mu_num
%             cmb = nchoosek(ac_mu_idx, i);
%             [cmb_num,~] = size(cmb);
%             for cmb_idx = 1:cmb_num
%                 res = spikes(:,spk_idx) - sum(tmplts(:,cmb(cmb_idx,:)),2);
%                 if sum(res.^2)/spk_pwr < 0.25 && hypo_idx < max_hypo_num
%                     hypo_idx = hypo_idx + 1;
%                     hypo(spk_idx,cmb(cmb_idx,:),hypo_idx) = 1;
%                 end
%             end
%         end
%     end
% %=========-- End of Generate hypothesis =================================--
%     
% end

%==============================EOF=======================================--
