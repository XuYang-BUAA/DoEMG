function [asr] = calASR(t, d)
% =========================================================================
%        estimate of the probability that MUAP occured                    *
%                                                                         *
%  INPUT:                                                                 *
%    t               -- mu templates shap                                         *
%    d               -- spike timing                                  *
%                                                                         *
%  OUTPUT:                                                                *
%    asr             --                           *
%                                                                         *
%  WARNINGS:   none                                                       *
%                                                                         *
%  HISTORY:                                                               *
%    7/2/2020 : XuY update                                                *
% =========================================================================
    [t_len,t_num] = size(t);
    [d_len,d_num] = size(d);
    assert(t_len == d_len, 'Matrix t and d should have same number of raws!');
    
    amp = max(abs(t));
    amp_thresh = amp / max(amp);
    amp_thresh = repmat(amp_thresh',1,d_num);
    norm2t = sum(t.^2); % 1 x t
    norm2d = sum(d.^2); % 1 x d
    dotmulti = t' * d; % t x d
    A = dotmulti ./ repmat(norm2t', 1, d_num); % t x d
    coeff = A; % t x d
    coeff(abs(A-1)<amp_thresh) = 1;
    coeff(A>(1+amp_thresh)) = 1./coeff(A>(1+amp_thresh));
    coeff(A<0) = 0;
    asr = dotmulti ./ sqrt(norm2t' * norm2d) .* coeff;
end

