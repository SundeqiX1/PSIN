#include "KFApp.h"

/***************************  class CKFApp  *********************************/
CKFApp::CKFApp(double ts):CSINSGNSS(19, 6, ts)
{
//state: 0-2 phi; 3-5 dvn; 6-8 dpos; 9-11 eb; 12-14 db; 15-17 lever; 18 dt
//meas:  0-2 dvn; 3-5 dpos
}

void CKFApp::Init(const CSINS &sins0, int grade)
{
	CSINSGNSS::Init(sins0);
	Pmax.Set2(fPHI(600,600),  fXXX(500),  fPOS(1e6),  fDPH3(5000),  fMG3(10), fXXX(10),  0.1);
	Pmin.Set2(fPHI(0.1,1.0),  fXXX(0.001),  fPOS(0.1),  fDPH3(0.1),  fUG3(10), fXXX(0.01),  0.0001);
	Pk.SetDiag2(fPHI(60,600),  fXXX(1.0),  fPOS(10.0),  fDPH3(100),  fMG3(3.0), fXXX(1.0),  0.01);
	Qt.Set2(fDPSH3(0.1),  fUGPSHZ3(1.0),  fOOO,  fOO6,	fOOO, 0.0);
	Rt.Set2(fXXZ(0.5,1.0),   fLLH(10.0,30.0));
	Rmax = Rt*100;  Rmin = Rt*0.01;  Rb = 0.6;
	FBTau.Set(fXX9(0.1),  fXX6(1.0),  fINF3, INF);
}

void CKFApp::SetMeas(const CVect3 &pgps, const CVect3 &vgps)
{
	if(!IsZero(pgps))
	{
		*(CVect3*)&Zk.dd[3] = sins.pos - pgps;
		SetMeasFlag(000070);
	}
	if(!IsZero(vgps))
	{
		*(CVect3*)&Zk.dd[0] = sins.vn - vgps;
		SetMeasFlag(000007);
	}
}

int CKFApp::Update(const CVect3 *pwm, const CVect3 *pvm, int nSamples, double ts)
{
	int res=TDUpdate(pwm, pvm, nSamples, ts, 5);
	return res;
}

