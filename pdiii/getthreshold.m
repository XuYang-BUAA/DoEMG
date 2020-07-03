function [thresholds, sds] = getthreshold (sig, COEF)
% =========================================================================
% Get detection threshold and baseline noise sd based on baseline noise.  *
% If sig/sig.data is a Matrix, threshold and sd will be calculated for    *
% every COLUMN.                                                           *
%                                                                         *
% INPUT:                                                                  *
%    sig             -- filtered multi-channel sEMG data                  *
%    COEF            -- typically 4                                       *
%                                                                         *
% OUTPUT:                                                                 *
%    thresholds      -- thresholds for peak detection algorithm           *
%    sds             -- standard deviation of baseline noise              *
%                                                                         *
%  WARNINGS:                                                              *
%    None                                                                 *
%                                                                         *
%  HISTORY:                                                               *
%    06/29/2020 XY : Update.                                              *
% =========================================================================
    if nargin<2
        COEF = 4; 
    end
    if isstruct(sig)
		sig = sig.data;
    end
    [l, n] = size(sig);
    % sorts the elements of A in ascending order.
    % If A is a matrix, then sort(A) treats the columns of A as vectors and sorts each column
	s = sort(abs(sig(1:l,:)));
	sd = sqrt(cumsum(s.^2)./repmat((1:l)', 1, n));
    sds = zeros(1, n);
    for chn = 1:n
        i = find(COEF*sd(:,chn)>s(:,chn), 1, 'last');
        if isempty (i)
            fprintf(sprintf('[Warning]:No physiological baseline.'));
            thresholds = s(round(.99*l), chn);  % non-physiological baseline
            sds(chn) = sd(round(.99*l), chn);
        else
            thresholds(chn) = 2*s(i,chn);
            sds(chn) = sd(i,chn);
        end
    end
end
