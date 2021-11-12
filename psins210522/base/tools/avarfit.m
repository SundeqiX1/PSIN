function avarfit(sigma, tau)
% Gyro Allan variance fit tool.
%
% Prototype: avarfit(sigma, tau)
% Inputs: sigma - gyro Allan variance array, in deg/hur
%         tau - Allan variance correlated time array, in sec
%
% Example: 
%     y = avarsimu([0.001,0.05,0.01,.01], [1, 1], 0.01, 1000000, 1);
%     [sigma,tau] = avar(y, 0.01);  avarfit(sigma, tau);
%
% See also  avar, avarsimu.

% Copyright(c) 2009-2021, by Gongmin Yan, All rights reserved.
% Northwestern Polytechnical University, Xi An, P.R.China
% 18/05/2021
global gtau hls ys mkva1 mkvb1 mkva2 mkvb2
    if nargin<2,
        tau=0.01*2.^(0:length(sigma)-1)'; % tau=0.01*[1,2,4,8,...]';
        if size(sigma,2)==2, tau=sigma(:,end); sigma=sigma(:,1); end
    end
    gtau = tau;
    ys = [tau.^-1, tau.^-0.5, tau.^0, tau.^0.5, tau.^1, ararmkv(1,1), ararmkv(1,1)];
    ys(:,1) = ys(:,1)*(sigma(1)/ys(1,1))*0.7;
    ys(:,2) = ys(:,2)*(sigma(1)/ys(1,2))*0.7;
    [m,idx] = min(sigma); mkva1=m/1000; mkvb1=m/5000; mkva2=m/5000; mkvb2=m/1000;
    ys(:,3) = ys(:,3)*(sigma(idx)/ys(idx,3))*0.7;
    ys(:,4) = ys(:,4)*(sigma(end)/ys(end,4))*0.7;
    ys(:,5) = ys(:,5)*(sigma(end)/ys(end,5))*0.7;
    ys(:,6) = ararmkv(mkva1,mkvb1);
    ys(:,7) = ararmkv(mkva2,mkvb2);
    ys(:,8) = normv(ys(:,1:7),2);
    myfig; hls(8)=loglog(tau, ys(:,8), '-g', 'linewidth',2); xygo('\tau / s', '\sigma_A / ( \circ / h )');
    hls(9)=loglog(tau, sigma, '-m*', 'linewidth',2);  legend('fitted AVAR', 'original AVAR');
    for k=1:7
        hls(k) = plot(tau, ys(:,k), '--');
        set(hls(k), 'ButtonDownFcn', {@bdfcn,k});
    end
    dispQNBKR();

    function bdfcn(~,~,hk)
        global prep hls
        set(hls(hk),'LineWidth',2,'Selected','on');
        prep = get(gca,'CurrentPoint');
        set(gcf,'WindowButtonMotionFcn',{@wbmfcn,hk});
        set(gcf,'WindowButtonUpFcn',{@wbufcn,hk});

    function wbmfcn(~,~,hk)
        global gtau prep hls ys mkva1 mkvb1 mkva2 mkvb2
        curp = get(gca,'CurrentPoint');
        if hk<6
            ys(:,hk) = ys(:,hk)*(curp(1,2)/prep(1,2));
        elseif hk==6
            [~,idx] = max(ys(:,6));
            if idx>2 && curp(1,1)<gtau(idx)
                mkva1 = mkva1*curp(1,2)/prep(1,2);
            elseif idx<length(gtau)-1 && gtau(idx)<curp(1,1)
                mkvb1 = mkvb1*curp(1,2)/prep(1,2);
            end
            ys(:,6) = ararmkv(mkva1, mkvb1);
        elseif hk==7
            [~,idx] = max(ys(:,7));
            if idx>2 && curp(1,1)<gtau(idx)
                mkva2 = mkva2*curp(1,2)/prep(1,2);
            elseif idx<length(gtau)-1 && gtau(idx)<curp(1,1)
                mkvb2 = mkvb2*curp(1,2)/prep(1,2);
            end
            ys(:,7) = ararmkv(mkva2, mkvb2);
        end
        prep = curp;
        set(hls(hk),'YData',ys(:,hk));
        set(hls(8),'YData',normv(ys(:,1:7),2));
        dispQNBKR();

    function wbufcn(~,~,hk)
        global hls
        set(gcf,'WindowButtonMotionFcn','');
        set(gcf,'WindowButtonUpFcn','');
        set(hls(hk),'LineWidth',1,'Selected','off');

    function dispQNBKR()
        global gtau ys mkva1 mkvb1 mkva2 mkvb2
        y1(5) = 0;
        for k=1:5, y1(k) = interp1(gtau, ys(:,k), 1.0);  end
        [y1(6),idx1] = max(ys(:,6));  [y1(6), tau1] = maxmkv(mkva1, mkvb1, gtau(idx1-1), gtau(idx1+1));
        [y1(7),idx2] = max(ys(:,7));  [y1(7), tau2] = maxmkv(mkva2, mkvb2, gtau(idx2-1), gtau(idx2+1));
        str = sprintf('Q=%.2f (\\prime\\prime),  N=%.6f (\\circ/h^{1/2}),  B=%.4f (\\circ/h),  K=%.6f (\\circ/h^{3/2}),  R=%.4f (\\circ/h^2),  M1=%.4f (\\circ/h)@%.3f (s),  M2=%.4f (\\circ/h)@%.3f (s)', ...
                y1(1)/sqrt(3), y1(2)/60, y1(3), y1(4)*sqrt(3)/60, y1(5)*sqrt(2)*3600, y1(6)/0.62, tau1/1.9, y1(7)/0.62, tau2/1.9);
        title(str);
        
	function s = ararmkv(a, b)
        global gtau
        s = sqrt(1./(1/a*gtau.^-1+1/b*gtau.^1));
        
	function [m, tau] = maxmkv(a, b, t1, t2)
        t = linspace(t1,t2);
        s = sqrt(1./(1/a*t.^-1+1/b*t.^1));
        [m, idx] = max(s);
        tau = t(idx);