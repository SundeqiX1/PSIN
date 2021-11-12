% 19-rotation-postion systematic calibration simulation.
% Ref:XieBo,Multiposition calibratiuon method of laser gyro SINS,JCIT,2011.
% Copyright(c) 2009-2019, by Gongmin Yan, All rights reserved.
% Northwestern Polytechnical University, Xi An, P.R.China
% 17/04/2019
glvs
ts = 0.01;
att0 = [0; -90; -90]*glv.deg;  pos0 = posset(34,0,0);
% 19-position setting
paras = [
    1    0,1,0, 90, 9
    2    0,1,0, 90, 9
    3    0,1,0, 90, 9
    4    0,1,0, -90, 9
    5    0,1,0, -90, 9
    6    0,1,0, -90, 9
    7    0,0,1, 90, 9
    8    1,0,0, 90, 9
    9    1,0,0, 90, 9
    10   1,0,0, 90, 9
    11   -1,0,0, 90, 9
    12   -1,0,0, 90, 9
    13   -1,0,0, 90, 9
    14    0,0,1, 90, 9
    15    0,0,1, 90, 9
    16    0,0,-1, 90, 9
    17    0,0,-1, 90, 9
    18    0,0,-1, 90, 9
];  paras(:,5) = paras(:,5)*glv.deg;
att = attrottt(att0, paras, ts);
% att1=att; att1(:,end)=att1(:,end)/2; att3ddemo(att1(1:10:end,:),1);
% IMU simulation
len = length(att);
avp = [att(:,1:3),zeros(len,3),repmat(pos0',len,1),att(:,4)];
imu = avp2imu(avp);
imuplot(imu,1);
% systematic calibration
imu1 = imuclbt(imu);
[clbt, av] = sysclbt(imu1, pos0);
clbtfile(clbt);
% clbtfile(clbt, 'clbt19.bin');
% binfile('imu19.bin', [imu1,imu1(:,end)]);

