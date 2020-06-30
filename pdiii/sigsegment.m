function [ sigout ] = sigsegment(sig, t, seglen, opt)
% Extract segments from a signal.

    if nargin==3 || lower(opt) == 'c'
        % t is the center of the segment.
		shift = -1;
	elseif lower(opt) == 'l'
        % t is the left edge of the segment.
		shift = 0;
    elseif lower(opt) == 'r'
        % t is the right edge of the segment.
        shift = -2;
    else
		error ('Unknown option.');
    end
    
    if isnan(t); t=[]; end
    
    t0 = 1;
    dt = 1;
    if isstruct(sig)
        sig = sig.data;
    end
    
    [sig_len, chn_num] = size(sig);
    
    seg_num = length(t);
	samp_num = round (seglen/dt);
    
    sig = [zeros(samp_num, chn_num); sig; zeros(samp_num, chn_num)];
    
    ibeg = round((t - t0)/dt) + shift * floor(samp_num/2) + 1 + samp_num;
    iend = ibeg + samp_num - 1;
    
    sigout = zeros(samp_num, seg_num, chn_num);
    for i = 1:seg_num
        sigout(:,i,:) = permute(sig(ibeg(i):iend(i),:), [1 3 2]); 
    end
end

