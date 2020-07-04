function [mu_firings , res_data] = checkMaxMu(data,mu_templates,thresh,spike_width)
% =========================================================================
%      determine locations where the mu_template may be present           *
%                                                                         *
%                                                                         *
%  INPUT:                                                                 *
%    data            -- sEMG data                                         *
%    spike_segment   -- spike segments                                    *
%    spike_train     -- timings of spile segments                         *
%    mu_templates    -- max mu shape                                      *
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
%    mu_firings      -- merged MU templates                               *
%    res_data        -- residual data                                     *
%                                                                         *
%  WARNINGS:   none                                                       *
%                                                                         *
%  HISTORY:                                                               *
%    7/2/2020 : XuY created                                               *
% =========================================================================
% =========== Basic information & Settings=================================
    mu_shapes = mu_templates.shape;
    mu_len = mu_templates.mu_len;
    mu_len_h = floor(mu_len/2);
    mu_ch_num = mu_templates.ch_num;
    mu_shapes_res = mu_shapes;
    mu_num = mu_templates.mu_num;
    
    [~,iteration_times] = size(mu_shapes_res); 
    
    
    
    data_len = length(data);
    res_data = data;
    mu_firings = zeros(data_len, mu_num);
% =========================================================================
    for i = 1:iteration_times
        [spikes, spike_train] = getspikes(res_data,thresh,spike_width,'MergeMode','Vary','Plot', 0);
        spikes = hanning(mu_len).*spikes;
        [~, spk_num, ~] = size(spikes);
        spikes = reshape(permute(spikes,[1 3 2]), [], spk_num);
        
        max_mu_index = findMaxMu(mu_shapes_res);
        max_mu_shape = mu_shapes_res(:,max_mu_index);
        %mu_shapes_res(:,max_mu_index) = [];
        mu_shapes_res(:,max_mu_index) = 0;
        P = calProb(max_mu_shape,mu_shapes_res,mu_len,mu_len_h);
        prob_threshold = max(P);
        
        P_spk = calProb(max_mu_shape,spikes,mu_len,mu_len_h);
        [~,firing_train] = findpeaks(P_spk,'MinPeakDistance', mu_len_h,'MinPeakHeight',prob_threshold);   %......................................
        %[~,firing_train] = findpeaks(P_spk,'MinPeakHeight',0.7);
        max_mu_firings = spike_train(firing_train);
        mu_firings(max_mu_firings, max_mu_index) = 1; 
        figure
        for k=1:4
            %subplot(4,1,k)
            plot(res_data(:,k)+1.5*(k-1),'b')
            hold on
            line([0,data_len],[1.5*(k-1),1.5*(k-1)],'linestyle','-')
        end
        res_data = subtractMu(res_data,max_mu_shape,max_mu_firings,mu_len_h);
        %test
        
        for k=1:4
            %subplot(4,1,k)
%             plot(data(:,k))
%             hold on
            plot(res_data(:,k)+1.5*(k-1),'r')
        end
        for n=1:length(max_mu_firings)
            line([max_mu_firings(n),max_mu_firings(n)],[-1,6],'linestyle',':')
        end
        figure
        for k=1:4
            plot(max_mu_shape(1+mu_len*(k-1):mu_len*k)+1.5*(k-1));
            hold on;
        end
%         figure
%         for k=1:mu_ch_num
%             subplot(2,2,k)
%             plot(data(max_mu_firings(10)-37:max_mu_firings(10)+37,k))
%             hold on
%             plot(res_data(max_mu_firings(10)-37:max_mu_firings(10)+37,k),'-')
%             plot(max_mu_shape(75*(k-1)+1:75*k),'*')
%         end
        hold on;
    end
    
end %end of checkMaxMu
% find mu with max amplitude
function max_mu_index = findMaxMu(mu_shapes)
    [~,max_mu_index] = max(sum(abs(mu_shapes)));
end % end of findMaxMu

%estimate of the probability that MU p occurred at segment d
function P = calProb(p , d, mu_len,mu_len_h)
    A = p\d;
    B = A;
    B(A>1) = 1./A(A>1);
    B(A<=1) = A(A<=1);
    B(A<0) = 0;
    norm2d = sum(d.^2);
    orthogonal = d-A.*p;
    P = B.*sqrt(1-sum(orthogonal.^2)./norm2d);
end % end of calProb

%subtract mu from emg data. return residual data for next iteration
function res_data = subtractMu(res_data,max_mu_shape,max_mu_firings,mu_len_h)
    for i=1:length(max_mu_firings)
        timing = max_mu_firings(i);
        A = max_mu_shape\reshape(res_data(timing-mu_len_h:timing+mu_len_h,:),[],1);
        res_data(timing-mu_len_h:timing+mu_len_h,:) = res_data(timing-mu_len_h:timing+mu_len_h,:) - A*reshape(max_mu_shape,[],4);
    end
end % end of subtractMu












