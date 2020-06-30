%--------------------------------------------------------------------------
function [new_miu, new_sigma] = calestimateIPI(miu, sigma, new_x, varargin)
% New estimate of interpulse interval using statistics or Kalman filter.
% If there are 4 parameters in, use statistics or use Kalman filter.

    if nargin == 3
        N = 10;
        new_miu = (1 - 1/N)*miu + 1/N*new_x;
        new_sigma = sqrt((1 - 1/N)*sigma^2+1/(N-1)*(new_x - new_miu)^2);
    else
        n = varargin{1};
        [new_miu, new_sigma] = calmiusigma(miu, sigma, n, new_x);
    end
end


%--------------------------------------------------------------------------
function [new_miu, new_sigma] = calmiusigma(miu, sigma, n, new_x)
% Use new data to updata the estimated mean value & standard deviation.
    new_miu = (miu*n + new_x)/(n+1);
    dmiu = new_miu - miu; dx = new_x - new_miu;
    new_sigma = sqrt(((n-1).*sigma.^2 + dx.^2 + n*dmiu.^2)/n);
    new_sigma = max(0, new_sigma);
end

