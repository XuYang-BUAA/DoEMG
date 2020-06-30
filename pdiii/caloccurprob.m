function [z] = caloccurprob(t, miu, sigma, varargin)
% Calculate the estimate of MU occurrence probability ('detection' hazrd
% function).
 
    [~,n] = size(varargin);
    % Default parameters.
    op_miss_prob = 0.1;
    op_max_miss_number = 10;
    op_ignore = 3.6;
    % Parse parameters.
    i = 1;
    while i <= n
        switch varargin{i}
            case 'MissProb'
                % The probability of miss detection of the last firing.
                op_miss_prob = varargin{i+1};
            case 'MaxMissNumber'
                % The max number of miss detection.
                op_max_miss_number = varargin{i+1};
            case 'IgnoreThresh'
                % Ignore the result when dt > ignore*miu, and return 1/miu.
                op_ignore = varargin{i+1};
            otherwise
                error('Unknown option!');
        end
        i = i+2;
    end
    
    % Calculate hazard function.
    z = zeros(size(t));
    t = t(:); miu = miu(:); sigma = sigma(:);
    n = (1:op_max_miss_number)'; p_pown = op_miss_prob.^(n-1);
    nsigma2 = n * sigma'.^2; nt = (repmat(t',size(n))-n*miu')./sqrt(nsigma2);
    Q = 1 - normcdf(nt);
    E = exp(-nt.^2./2);
    p = (p_pown' * diag(1./sqrt(2*pi*n)) * E)./ sigma' ./ (p_pown' * Q);
    z(:) = p(:);
    z(t-op_ignore*miu>0) = 1./(miu(t-op_ignore*miu>0));
    
end

