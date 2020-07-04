function [spikes, spike_times] = getspikes(sig, thresh, width, varargin)
% =========================================================================
%                  Find spike shapes and spike times.                     *
%                                                                         *
%  INPUT:                                                                 *
%    sig             -- filtered multi-channel sEMG data                  *
%    thresh          -- trigger threshold                                 *
%    width           -- wave width (typically 5-15ms)                     *
%    Name,Value      -- 
%                                                                         *
%  OUTPUT:                                                                *
%    spikes          -- spikes waves with input width                     *
%    spike_times     -- times of peaks                                    *
%                                                                         *
%  WARNINGS:                                                              *
%    None                                                                 *
%                                                                         *
%  HISTORY:                                                               *
%    06/29/2020 XY : Update.                                              *
% ========================================================================= 
    if nargin < 3
        error('Not enough arguments!');
    else
        [op_merge_mode, op_flag_plot] = parse(varargin);
    end
    
    if isstruct(sig)
		sig = sig.data;
    end
    
    [sig_len, chn_num] = size(sig);
    
    % Detect spike timing in each channel.
    spike_times = cell(chn_num, 1);
    spike_peaks = cell(chn_num, 1);
    for chn = 1:chn_num
        [peaks, timings] = findpeaks(abs(sig(:,chn)), ...
            'MinPeakHeight', thresh(chn));
%             'MinPeakDistance', floor(width/2));
        spike_times{chn} = timings;
        spike_peaks{chn} = peaks;      %peaks 
    end
    
    % Merge the spike times.
    switch op_merge_mode
        case 'None'
            mergedspktms = cell2mat(spike_times);
        case 'All'
            mergedspktms = cell2mat(spike_times);
            mergedspktms = sort(mergedspktms);
        case 'MaxPeak'
            peaks = zeros(sig_len,1);
            for chn = 1:chn_num
                peaks(spike_times{chn}) = max(peaks(spike_times{chn}), spike_peaks{chn});
            end
            [~, mergedspktms] = findpeaks(peaks, 'MinPeakDistance', floor(width/2));
        case 'Vary'
            peaks = zeros(sig_len,1);
            for chn = 1:chn_num
                peaks(spike_times{chn}) = max(peaks(spike_times{chn}), spike_peaks{chn});
            end
            distance = floor(width/2);
            vary_range = floor(distance/2);
            [~, mergedspktms] = findpeaks(peaks, 'MinPeakDistance', distance);
            spik_range = repmat((-vary_range:vary_range),length(mergedspktms),1);
            
            mergedspktms = repmat(mergedspktms,1,2*vary_range+1)+spik_range;
            mergedspktms = sort(reshape(mergedspktms,[],1));
        otherwise
            error('Unkown value for "MergeMode"!');
    end
    
    % Mark the peaks and plot.
    if op_flag_plot
        figure();
        h = markpeaks(sig, mergedspktms);
        %h = markpeaks(peaks,mergedspktms);
    end
    
    % Extract signal segments.
    spikes = sigsegment(sig, mergedspktms, width);
    spike_times = mergedspktms';
    
end

%--------------------------------------------------------------------------
function [op_merge_mode, op_flag_plot] = parse(opt)
% Parse options of "getspikes".
    [~,n] = size(opt);
    % Default arguments.
    op_merge_mode = 'MaxPeak'; op_flag_plot = false;
    
    i = 1;
    while i <= n
        switch opt{i}
            case 'MergeMode'
                % How to merge the spike times
                op_merge_mode = opt{i+1};
            case 'Plot'
                % Flag if there need a plot.
                op_flag_plot = opt{i+1};
            otherwise
                error('Unknown option!');
        end
        i = i+2;
    end
    
end

