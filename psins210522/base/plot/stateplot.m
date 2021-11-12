function sstate = stateplot(state, maxst)
% Runing-state plot.
%
% Prototype: sstate = stateplot(state, maxst)
% Inputs: state - runing state
%         maxst - max state
%          
% See also  kffile, imuplot, insplot, inserrplot, kfplot, rvpplot, gpsplot.

% Copyright(c) 2009-2017, by Gongmin Yan, All rights reserved.
% Northwestern Polytechnical University, Xi An, P.R.China
% 18/03/2017
    [n, m] = size(state);
    if m==1,  t = (1:n)';
    else      t = state(:,end); state = state(:,end-1); end
    sstate = [];
    for k=0:31
        tmp = bitand(state,2^k);
        mtmp = max(tmp);
        if mtmp>0, sstate = [sstate, tmp/mtmp*(k+1)]; end
    end
    myfigure, plot(t, sstate, '*'); xygo('Running-state');
    maxst0 = max(max(sstate));
    if nargin<2, maxst=0; end;
    maxst = max(maxst0,maxst)+0.2;
    ylim([0,maxst]);
