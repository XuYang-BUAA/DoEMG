function [h] = markpeaks( sig, spike_time, varargin )
% Mark peaks on the signal.
    if nargin == 3
        h = varargin{1};
    else
        h = axes();
    end

    t0 = 1;
    dt = 1;
    if isstruct(sig)
        t0 = sig.t0;
        dt = sig.dt;
        sig = sig.data;
    end
    
    [sig_len, chn_num] = size(sig);
    
    sig = 0.5 * sig ./ max((max(abs(sig)))) + ...
        repmat(0:chn_num-1, sig_len, 1);
    t = t0 + dt * (0:sig_len-1)';
    
%     if ~iscell(spike_time)
%         spike_time = {spike_time,spike_time,spike_time,spike_time};
%     end
    
    %l = plot(h, t, sig, '--');
    l = plot(h, sig);

    h.NextPlot = 'add';
    for chn = 1:chn_num
        %spike_index = round((spike_time - t0)/dt) + 1;
        hold on
        plot(h, spike_time, sig(spike_time,chn),'o');
    end
    h.NextPlot = 'replacechildren';
end

