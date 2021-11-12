function imu = imuadderr(imu, imuerr, ts)
% SIMU outputs adding errors simulation, denoted as:
%    imu = K*imu + drift error.
%
% Prototype: imu = imuadderr(imu, imuerr, ts)
% Inputs: imu - raw SIMU data
%         imuerr - SIMU error struture array
%         ts - SIMU sample interval
% Output: imu - output SIMU data added errors
%
% See also  imuerrset, imuclbt, imudeldrift, imugsensi, avperrset, trjsimu, insupdate.

% Copyright(c) 2009-2014, by Gongmin Yan, All rights reserved.
% Northwestern Polytechnical University, Xi An, P.R.China
% 11/09/2013, 06/03/2014, 21/07/2015, 17/08/2016, 24/02/2017, 24/07/2020
    if nargin<3,  ts = imu(2,7)-imu(1,7);  end  % the last column implies sampling interval
    if isfield(imuerr, 'rx')  % inner-lever-arm effect for accelerometers
        wb = imu(:,1:3)/ts;
        dotwf = imudot(imu, 5.0);
        fL = [ (-dotwf(:,2).^2-dotwf(:,3).^2+0)*imuerr.rx(1) + (dotwf(:,1).*dotwf(:,2)-wb(:,3))*imuerr.rx(2) + (dotwf(:,1).*dotwf(:,3)+wb(:,2))*imuerr.rx(3), ... 
               (dotwf(:,1).*dotwf(:,2)+wb(:,3))*imuerr.ry(1) + (-dotwf(:,1).^2-dotwf(:,3).^2+0)*imuerr.ry(2) + (dotwf(:,2).*dotwf(:,3)-wb(:,1))*imuerr.ry(3), ... 
               (dotwf(:,1).*dotwf(:,3)-wb(:,2))*imuerr.rz(1) + (dotwf(:,2).*dotwf(:,3)+wb(:,1))*imuerr.rz(2) + (-dotwf(:,1).^2-dotwf(:,2).^2+0)*imuerr.rz(3)  ]; 
        imu(:,4:6) = imu(:,4:6) + fL*ts;
    end
    if isfield(imuerr, 'dtGA') % time-asynchrony between gyro & acc
        acc = [ [imu(1,4:6),imu(1,7)-ts]; imu(:,4:7); [repmat(imu(end,4:6),10,1),imu(end,7)+(1:10)'*ts] ];
        acc(:,1:3) = cumsum(acc(:,1:3),1);
        t = [imu(1,end)-ts; imu(:,end)] + imuerr.dtGA;
        for k=1:3
            acc1 = interp1(acc(:,end), acc(:,k), t, 'linear');
            imu(:,3+k) = diff(acc1,1);
        end
    end
    [m, n] = size(imu); sts = sqrt(ts);
    drift = [ ts*imuerr.eb(1) + sts*imuerr.web(1)*randn(m,1), ...
              ts*imuerr.eb(2) + sts*imuerr.web(2)*randn(m,1), ...
              ts*imuerr.eb(3) + sts*imuerr.web(3)*randn(m,1), ...
              ts*imuerr.db(1) + sts*imuerr.wdb(1)*randn(m,1), ...
              ts*imuerr.db(2) + sts*imuerr.wdb(2)*randn(m,1), ...
              ts*imuerr.db(3) + sts*imuerr.wdb(3)*randn(m,1) ];
    if min(abs(imuerr.sqg))>0
        mvg = markov1(imuerr.sqg.*sqrt(imuerr.taug/2), imuerr.taug, ts, m);   % q = 2*sigma.^2.*beta
        drift(:,1:3) = drift(:,1:3) + mvg*ts;
    end
    if min(abs(imuerr.sqa))>0
        mva = markov1(imuerr.sqa.*sqrt(imuerr.taua/2), imuerr.taua, ts, m);
        drift(:,4:6) = drift(:,4:6) + mva*ts;
    end
    if min(abs(imuerr.KA2))>0
        imu(:,4:6) = [ imu(:,4)+imuerr.KA2(1)/ts*imu(:,4).^2, ...
                       imu(:,5)+imuerr.KA2(2)/ts*imu(:,5).^2, ...
                       imu(:,6)+imuerr.KA2(3)/ts*imu(:,6).^2 ];
    end
    if isfield(imuerr, 'dKg')
        Kg = eye(3)+imuerr.dKg; Ka = eye(3)+imuerr.dKa;
        imu(:,1:6) = [imu(:,1:3)*Kg', imu(:,4:6)*Ka'];
    end
    imu(:,1:6) = imu(:,1:6) + drift;
