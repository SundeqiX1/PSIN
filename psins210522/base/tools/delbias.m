function x = delbias(x, b)
% Delete bias.
%
% Prototype:  out = delbias(in, b)
% Inputs: in - input date with bias
%         b - bias
% Output: out - output date with no bias
%
% See also  imudeldrift, adddt, sumn, meann.

% Copyright(c) 2009-2016, by Gongmin Yan, All rights reserved.
% Northwestern Polytechnical University, Xi An, P.R.China
% 23/07/2016, 03/07/2018 
    for k=1:size(x,2)
        if exist('b','var')
            if k<=length(b), x(:,k) = x(:,k)-b(k); end
        else
            x(:,k) = x(:,k)-mean(x(:,k));
        end
    end