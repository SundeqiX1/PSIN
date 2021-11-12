#include ".\PSINSCore\kfapp.h"
#include ".\PSINSCore\psins_demo.h"

void main(void)
{
	if(PSINSDemo>=0) {	psinsdemo(); exit(0); }
	// else (PSINSDemo==-1) ...
	CFileRdWt::Dir("D:\\ygm2020\\PSINS��վ\\�ߵ�����\\", "D:\\psins210207\\VC60\\Data\\");
	CFileRdWt fins("ins.bin"), fkf("kf.bin");
	CFileRdSr fimu("mimuattgps.bin");  // download from: http://www.psins.org.cn/newsinfo/958984.html
	DataSensor *pDS=(DataSensor*)fimu.buff, *pDS0=&fimu.DS0;

	CKFApp kf(TS);
	kf.Init(CSINS(pDS0->att, pDS0->gpsvn, pDS0->gpspos, pDS0->t));

	for(int i=0; i<5000*FRQ; i++)
	{
		if(!fimu.load(1)) break;
		if(pDS->gpspos.i>0.1 && !hit3(pDS->t,500,600,900,1000,2000,2100))
		{
			kf.SetMeas(pDS->gpspos, pDS->gpsvn);
		}
		kf.Update(&pDS->wm, &pDS->vm, 1, TS);

		if(i%5==0||pDS->gpspos.i>0.1)
		{
			fins << kf.sins << pDS->att;
			fkf << kf;
		}

		disp(i, FRQ, 100);
	}
}
