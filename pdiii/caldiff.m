function [K1sqr, K2sqr, K3sqr, A] = caldiff(candi, template)
% =========================================================================
% Calculate the degree of difference between candi & template.            *
% If candi(Mc x Nc) or template(Mt x Nt) is matrix, calculate on column,  *
% and the result is Nc x Nt.                                              *
%                                                                         *
%                                                                         *
%  INPUT:                                                                 *
%    candi           -- candidate firing segment for classification       *
%    template        -- MU templates                                      *
%                                                                         *
%  OUTPUT:                                                                *
%  Denote candi, template as vector x, y.                                 *
%    K1sqr           -- K1^2 , K1 = |A - 1|                               *
%    K2sqr           -- K2^2 , K2 = sqrt(|x|^2 / |y|^2 - A^2);            *
%    K3sqr           -- K3^2 , K3 = sqrt(K2^2 + (A - 1)^2).               *
%    A               -- A = x' * y/(y' * y);                              *
%                                                                         *
%  WARNINGS:   none                                                       *
%                                                                         *
%  HISTORY:                                                               *
%    7/1/2020 : XuY update                                                *
% =========================================================================
    [cnd_len, cnd_num] = size(candi);
    [tmp_len, ~] = size(template);
    %hanning window
%     hanning_coef = hanning(tmp_len);
%     candi = candi.*hanning_coef;
%     template = template.*hanning_coef;
    
    if (cnd_len ~= tmp_len)
        error('Vectors should have the same length!');
    end
    tmp2 = sum(template.^2);
    cnd2 = sum(candi.^2);
    A = candi' * template ./ repmat(tmp2, cnd_num, 1);
    K1sqr = abs(A - 1).^2;
    K2sqr = cnd2' * (1./tmp2) - A.^2;
    K3sqr = K2sqr + K1sqr;
    
end

%------------------------------EOF-----------------------------------------

