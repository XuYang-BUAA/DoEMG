function w = goodwidth (width, dt)
% GOODWIDTH		Find highly composite odd number
%	GOODWIDTH (N) returns a highly composite odd number
%	that is greater than or equal to N. This number
%	provides an efficient length for trigonometric polynomials.
%
%	GOODWIDTH (W, DT) returns a length in seconds that 
%	corresponds to an efficient odd number of samples when
%	the sampling interval is DT.

% Copyright (c) 2006-2009. Kevin C. McGill and others.
% Part of EMGlab version 1.0.
% This work is licensed under the Aladdin free public license.
% For copying permissions see license.txt.
% email: emglab@emglab.net

	if nargin==1
		n = ceil(width);
		w = inf;
        for n7 = 0: ceil(log(n)/log(7))
			r7 = n/7^n7;
            for n5 = 0: ceil(log(r7)/log(5))
				r5 = r7/5^n5;
				n3 = max(ceil(log(r5)/log(3)),0);
				w = min (w, 3^n3*5^n5*7^n7);
            end
        end
	else
		w = goodwidth (ceil(width/dt)) * dt;
	end
end
    
