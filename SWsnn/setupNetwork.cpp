
#include "mysnn.h"

#include <sys/stat.h> // mkdir

#include <sys/time.h>    // for gettimeofday() chenged!!

#define COMPACTION_ALIGNMENT_PRE  16
#define COMPACTION_ALIGNMENT_POST 0

inline int isnan(double x) { return x != x; }
inline int isinf(double x) { return !isnan(x) && isnan(x - x); }
inline int iserr(double x) {return isnan(x)||isinf(x);}


static float getWeights(grpInfo_t *gInfo,int connProp, float initWt, float maxWt, unsigned int nid, int gid);
static void connectFull(snnInfo_t *sInfo,connInfo_t* cInfo,grpInfo_t *gInfo); 
static void connectPart(snnInfo_t *sInfo,connInfo_t* cInfo,grpInfo_t *gInfo); 
static void resetNeuron(neurInfo_t *nInfo,int nid,grpInfo_t *gInfo,int gid);
static void resetNeuron2(snnInfo_t *sInfo,int nid,grpInfo_t *gInfo,int gid);
static void resetSynapticConnections();
static void updateSpikeGeneratorsInit(grpInfo_t* grpInfo, int numGrp); 

void setConductances(snnInfo_t *sInfo,bool isSet,int tdAMPA,int trNMDA,int tdNMDA,int tdGABAa,int trGABAb,int tdGABAb) {
	if (isSet) {
		assert(tdAMPA>0); assert(tdNMDA>0); assert(tdGABAa>0); assert(tdGABAb>0);
		assert(trNMDA>=0); assert(trGABAb>=0); // 0 to disable rise times
		assert(trNMDA!=tdNMDA); assert(trGABAb!=tdGABAb); // singularity
	}

	// set conductances globally for all connections
	sInfo->sim_with_conductances  |= isSet;
	sInfo->dAMPA  = 1.0-1.0/tdAMPA;
	sInfo->dNMDA  = 1.0-1.0/tdNMDA;
	sInfo->dGABAa = 1.0-1.0/tdGABAa;
	sInfo->dGABAb = 1.0-1.0/tdGABAb;
#if 0
	if (trNMDA>0) {
		sInfo->sim_with_NMDA_rise = 1;
		sInfo->rNMDA = 1.0-1.0/trNMDA;
		double tmax = (-tdNMDA*trNMDA*log(1.0*trNMDA/tdNMDA))/(tdNMDA-trNMDA); // t at which cond will be max
		sInfo->sNMDA = 1.0/(exp(-tmax/tdNMDA)-exp(-tmax/trNMDA)); // scaling factor, 1 over max amplitude
		assert(!isinf(tmax) && !isnan(tmax) && tmax>=0);
		assert(!isinf(sInfo->sNMDA)&&!isnan(sInfo->sNMDA)&&sInfo->sNMDA>0);
	}

	if (trGABAb>0) {
		sInfo->sim_with_GABAb_rise = 1;
		sInfo->rGABAb = 1.0-1.0/trGABAb;
		double tmax=(-tdGABAb*trGABAb*log(1.0*trGABAb/tdGABAb))/(tdGABAb-trGABAb); // t at which cond will be max
		sInfo->sGABAb=1.0/(exp(-tmax/tdGABAb)-exp(-tmax/trGABAb)); // scaling factor, 1 over max amplitude
		assert(!isinf(tmax) && !isnan(tmax));
		assert(!isinf(sInfo->sGABAb)&&!isnan(sInfo->sGABAb)&&sInfo->sGABAb>0);
	}
#endif
	//if (sInfo->sim_with_conductances) {
	//	//printf("AMPA decay time= %5d ms\n", tdAMPA);
	//} else {
	//}
}

// set Izhikevich parameters for group
void setNeuronParameters(grpInfo_t *gInfo,int gid,float izh_a,float izh_a_sd,float izh_b,float izh_b_sd,float izh_c,float izh_c_sd,float izh_d,float izh_d_sd)
{
	assert(gid>=0);
	assert(izh_a_sd>=0);assert(izh_b_sd>=0);
	assert(izh_c_sd>=0);assert(izh_d_sd>=0);

	gInfo[gid].Izh_a = izh_a;
	gInfo[gid].Izh_a_sd  =  izh_a_sd;
	gInfo[gid].Izh_b = izh_b;
	gInfo[gid].Izh_b_sd = izh_b_sd;
	gInfo[gid].Izh_c = izh_c;
	gInfo[gid].Izh_c_sd = izh_c_sd;
	gInfo[gid].Izh_d = izh_d;
	gInfo[gid].Izh_d_sd = izh_d_sd;
}

static void connectFull(snnInfo_t *sInfo,connInfo_t* cInfo,grpInfo_t *gInfo) {
}

//#define PROPAGATED_BUFFER_SIZE  (1023)

void initNetwork(snnInfo_t *sInfo, grpInfo_t *gInfo,connInfo_t *cInfo,
		int numGrp,int numConn,int randSeed){//????
	int i,g;

	srand48(randSeed); sInfo->randSeed=randSeed;

	sInfo->simTimeMs = 0; sInfo->simTimeSec = 0; sInfo->simTimeMs = 0;

	sInfo->maxDelay = 1; //sInfo->maxDelay_= 1;

	sInfo->numGrp = numGrp;	sInfo->numConn= numConn;

	sInfo->numN = 0; sInfo->numNReg = 0; 
	sInfo->numNExcReg = 0; sInfo->numNInhReg = 0; sInfo->numNPois = 0;

	sInfo->maxSpikesD1 = 0; sInfo->maxSpikesD2 = 0;

	sInfo->simTimeRunStart = 0; sInfo->simTimeRunStop = 0;
	
	sInfo->secD1fireCntHost = 0; sInfo->secD2fireCntHost = 0;
	sInfo->spikeCountAll1secHost = 0;
	
	sInfo->MaxFiringRate = 60;//60Hz ????
	sInfo->sim_with_homeostasis = 0; sInfo->sim_with_conductances = 1;
	sInfo->sim_with_NMDA_rise = 0; sInfo->sim_with_GABAb_rise = 0;

	//sInfo->pbuf = (Pbuf_t*)malloc(sizeof(Pbuf_t));
	//pbufInit(sInfo->pbuf,0,PROPAGATED_BUFFER_SIZE);

	// some default decay and rise times
	sInfo->dAMPA = 1.0-1.0/5.0;
	sInfo->rNMDA = 1.0-1.0/10.0; sInfo->dNMDA = 1.0-1.0/150.0;
        sInfo->sNMDA = 1.0;
        sInfo->dGABAa = 1.0-1.0/6.0;
        sInfo->rGABAb = 1.0-1.0/100.0; sInfo->dGABAb = 1.0-1.0/150.0;
        sInfo->sGABAb = 1.0;

	sInfo->timeTableD1 = NULL; sInfo->timeTableD2 = NULL;
	sInfo->firingTableD1 = NULL; sInfo->timeTableD2 = NULL;
	//sInfo->grpIds = NULL; //maybe unused!!
	sInfo->nSpikeCnt = NULL; sInfo->lastSpikeTime = NULL;
	
	for (i=0;i<numGrp;i++){
		gInfo[i].Type = 0;
		gInfo[i].StartN = -1;gInfo[i].EndN = -1;gInfo[i].SizeN = -1;
		
		gInfo[i].MaxDelay = sInfo->maxDelay;
		gInfo[i].FiringCount1sec = 0;

		gInfo[i].isSpikeGenerator= 0;gInfo[i].NewTimeSlice = 0;
		gInfo[i].CurrTimeSlice = 0;gInfo[i].SliceUpdateTime = 0;
		gInfo[i].RefractPeriod = 0;gInfo[i].RatePtr = NULL;

		gInfo[i].WithSTP = 0;gInfo[i].WithHomeostasis = 0;
		
		gInfo[i].numSynPre = 0;	gInfo[i].numSynPost = 0;

		gInfo[i].Izh_a = 0.; gInfo[i].Izh_a_sd = 0.;
		gInfo[i].Izh_b = 0.; gInfo[i].Izh_b_sd = 0.;
		gInfo[i].Izh_c = 0.; gInfo[i].Izh_c_sd = 0.;
		gInfo[i].Izh_d = 0.; gInfo[i].Izh_d_sd = 0.;

	}
	
	for(i=0;i<numConn;i++){
		cInfo[i].connId = i;
		cInfo[i].maxDelay = sInfo->maxDelay;
		cInfo[i].grpSrc = -1;cInfo[i].grpDest = -1;
		cInfo[i].initWt = 0.;cInfo[i].maxWt = 0.;
		cInfo[i].mulSynFast = 1.;cInfo[i].mulSynSlow = 1.;
		cInfo[i].maxSynPost = -1;cInfo[i].maxSynPre = -1;
		cInfo[i].numSyn = -1;cInfo[i].connProp = 0;cInfo[i].p = 0.;
	}
	
	return;
}

void freeNetwork(snnInfo_t *sInfo){

	MPI_Group_free(&sInfo->lcgroup);
	MPI_Comm_free(&sInfo->lccomm);
	return;
}

static void connectPart(snnInfo_t *sInfo,connInfo_t* cInfo,grpInfo_t *gInfo) {
    	int rank, nproc;
    	MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	//int gSrc=cInfo->grpSrc;
	//int gDest=cInfo->grpDest;
	uint8_t minD=1;
	uint8_t maxD=cInfo->maxDelay;
	assert(maxD-1==sInfo->Ndelay);
	bool noDirect = 0;
	float synWt;
	int preN=sInfo->preN;
	int Nsyn=sInfo->Nsyn;
	int Ndelay=sInfo->Ndelay;
	int MaxN=sInfo->MaxN;
	int SizeN=sInfo->swInfo[0].SizeN;
	assert(sInfo->numNReg%sInfo->Nsyn==0);
//	printf("connectPart::preN=%d Ndelay=%d NTh=%d MaxN=%d SizeN=%d\n",preN,Ndelay,NTh,MaxN,SizeN);
	/******init*****/
	long is;
	long len=sInfo->preN*Ndelay*NTh*MaxN;
	//printf("connectPart::sInfoHost %d\n",sInfo->sInfoHost);
	//printf("connectPart::%ld %d %d %d\n",len,preN,NTh,MaxN);
	for(is=0;is<len;is++){
//printf("connectPart::is=%d\n,size=%d\n",is,sizeof(synInfo_t));
		sInfo->sInfoHost[is].postId=0xffffffff;
		///sInfo->sInfoHost[is].dl=0;
		//sInfo->sInfoHost[is].wt=0./0.;
	}
	/***************/
	
		assert(sInfo->Ndelay<=20);
	int MIND=minD<<sInfo->Nop;
	assert(MIND==sInfo->Ndt);
	int MAXD=maxD<<sInfo->Nop;
	assert((MAXD-MIND)%sInfo->Ndt==0);
		int NNN; //select 1,2,3,4
		NNN = sInfo->numNReg/sInfo->Nsyn;
		assert(sInfo->numNReg%sInfo->Nsyn==0);
		assert(NNN>0);
	int i,j;
	int Ndma[NTh];
	for(i=0;i<NTh;i++)Ndma[i]=0;
	for(i=0;i<sInfo->preN;i++)  {
		int offset=i*Ndelay*NTh*MaxN;
		if(i==sInfo->numNReg) cInfo++;

		//float fac=1.0;
		//if(i>=gInfo[1].StartN && i<sInfo->numNReg) fac=-1.0;

		int iD=0,iTh=0,iN[20];
		for(iD=0;iD<Ndelay;iD++)iN[iD]=0;
		//for(j=0;j<sInfo->numNReg;j++) {
		for(j=i%NNN;j<sInfo->numNReg;j+=NNN) {
			if(j<sInfo->gStart||j>=sInfo->gStart+sInfo->gSize) continue;
			if(j>=sInfo->swInfo[iTh].StartN+sInfo->swInfo[iTh].SizeN){
				iTh++;
				for(iD=0;iD<Ndelay;iD++)iN[iD]=0;
			}
			if((noDirect)&&i==j) continue;
			uint8_t dVal;
			for(;;){
				//dVal=MIND+rand()%(MAXD-MIND);//for delay!!
				//dVal=(dVal>>sInfo->Nop)<<sInfo->Nop;//for test
				dVal=((i+j/NNN)%sInfo->Ndelay+1)<<sInfo->Nop;//for test
				assert((dVal>=MIND)&&(dVal<MAXD));
				iD = (dVal-MIND)>>sInfo->Nop;
				assert(iD<sInfo->Ndelay && iD>=0);
				if(iN[iD]<MaxN) break;
				else {
//printf("pre:%d post:%d iD:%d iN:%d %d %d %d %d %d %d %d %d %d\n",i,j,iD,iN[0],iN[1],iN[2],iN[3],iN[4],iN[5],iN[6],iN[7],iN[8],iN[9]);
					//assert(0);
				}
			}
			synWt=cInfo->initWt;
			//synWt=cInfo->maxWt*(float)drand48();
			assert(synWt>=0.);
			//synWt=synWt*fac;//for Inh neuron
			int addr=offset+iD*NTh*MaxN+iTh*MaxN+iN[iD];
			sInfo->sInfoHost[addr].postId=j;
			sInfo->sInfoHost[addr].wt=synWt;
/*****************************************************************/
			///sInfo->sInfoHost[addr].dl=dVal;
/*********************************************************************/
			iN[iD]++;
			if(iN[iD]>Ndma[iTh]) Ndma[iTh]=iN[iD];
			/******************/
		}
//printf("pre:%d post:%d iD:%d iN:%d %d %d %d %d %d %d %d %d %d\n",i,j,iD,iN[0],iN[1],iN[2],iN[3],iN[4],iN[5],iN[6],iN[7],iN[8],iN[9]);
	}
	//for(i=0;i<NTh;i++)sInfo->swInfo[i].Ndma = Ndma[i];
	for(i=0;i<NTh;i++)sInfo->swInfo[i].Ndma = MAX(Ndma[i],1);
	if(rank==0)printf("connectPart::Ndma=%d Ndma2=%d MaxN=%d\n",Ndma[0],Ndma[NTh-1],MaxN);
	return;
}//connectPart

void buildModel(snnInfo_t *sInfo,grpInfo_t *gInfo,connInfo_t *cInfo){//????
	int i,j;
    	int rank, nproc;
    	MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    	int lcrank, lcnproc;
    	MPI_Comm_size(sInfo->lccomm, &lcnproc);
    	MPI_Comm_rank(sInfo->lccomm, &lcrank);

	//allocate neurons in blocks
	int Size=sInfo->BaSize/lcnproc;
	int least=sInfo->BaSize%lcnproc;
	if (lcrank<least){
		sInfo->gStart = sInfo->BaStart+lcrank*(Size+1);
		sInfo->gSize  = Size+1;
	}else{
		sInfo->gStart = sInfo->BaStart+lcrank*Size+least;
                sInfo->gSize  = Size;
	}

	if (sInfo->Nsyn==0) sInfo->Nsyn=1000;
	//if (sInfo->Nsyn>sInfo->numNReg) sInfo->Nsyn=sInfo->numNReg;
	assert(sInfo->Nsyn<=sInfo->numNReg);
fprintf(stderr,"rank=%d nproc=%d gSize=%d gStart=%d numNReg=%d Nsyn=%d\n",rank,nproc,sInfo->gSize,sInfo->gStart,sInfo->numNReg,sInfo->Nsyn);fflush(NULL);

	sInfo->Ndelay=sInfo->maxDelay-1;
	if(sInfo->Ndelay==0) assert(0);	
	sInfo->Ndt= 8;sInfo->dt=1.0/sInfo->Ndt;sInfo->Nop= 3;

	sInfo->simTime = 0;
	sInfo->preN = sInfo->numN;
	//sInfo->MaxN = (sInfo->Nsyn/NTh/sInfo->Ndelay/nproc/4+1)*4;
	sInfo->MaxN = (sInfo->Nsyn/NTh/sInfo->Ndelay/nproc+1);
	if(rank==0)printf("preN=%d Nsyn=%d Ndelay=%d MaxN=%d\n",sInfo->preN,sInfo->Nsyn,sInfo->Ndelay,sInfo->MaxN);

	for(i=0;i<NTh;i++){
		int SizeN=(sInfo->gSize)/NTh;
		int least=sInfo->gSize%NTh;
		if(i<least){
			sInfo->swInfo[i].StartN = sInfo->gStart+i*(SizeN+1);
			sInfo->swInfo[i].SizeN = SizeN+1;
		}
		else{
			sInfo->swInfo[i].StartN = sInfo->gStart+i*(SizeN)+least;
			sInfo->swInfo[i].SizeN = SizeN;
		}
		sInfo->swInfo[i].MaxN = sInfo->MaxN;
		sInfo->swInfo[i].Ndt = sInfo->Ndt;
		sInfo->swInfo[i].Nop = sInfo->Nop;
		sInfo->swInfo[i].dt = sInfo->dt;
		sInfo->swInfo[i].sliceTime = 0;
		sInfo->swInfo[i].simTime = 0;
	}

	///allocate host memory
	//neuron data arrays
	sInfo->swPara.Izh_a = (float*)malloc(sInfo->gSize*sizeof(float));
	sInfo->swPara.Izh_b = (float*)malloc(sInfo->gSize*sizeof(float));
	sInfo->swPara.Izh_c = (float*)malloc(sInfo->gSize*sizeof(float));
	sInfo->swPara.Izh_d = (float*)malloc(sInfo->gSize*sizeof(float));

	sInfo->swPara.voltage = (float*)malloc(sInfo->gSize*sizeof(float));
	sInfo->swPara.recovery = (float*)malloc(sInfo->gSize*sizeof(float));
	sInfo->swPara.gAMPA = (float*)malloc(sInfo->gSize*sizeof(float));
	sInfo->swPara.gNMDA_d = (float*)malloc(sInfo->gSize*sizeof(float));
	sInfo->swPara.gGABAa = (float*)malloc(sInfo->gSize*sizeof(float));
	sInfo->swPara.gGABAb_d = (float*)malloc(sInfo->gSize*sizeof(float));

	//synapse data arrays
	/*sInfo->swPara.wt = (float*)malloc((long)sInfo->preN*sInfo->Ndelay*sInfo->MaxN*NTh*sizeof(float));
	sInfo->swPara.postId = (unsigned int*)malloc((long)sInfo->preN*sInfo->Ndelay*sInfo->MaxN*NTh*sizeof(unsigned int));*/

if(rank==0)fprintf(stderr,"rank=%d SizeN=%d StartN=%d  SizeN1=%d StartN1=%d \n",rank,sInfo->swInfo[0].SizeN,sInfo->swInfo[0].StartN,sInfo->swInfo[NTh-1].SizeN,sInfo->swInfo[NTh-1].StartN);
	//mpi mem alloc
	sInfo->NNgroup = sInfo->gSize;
	//sInfo->NN = sInfo->NNgroup*nproc;
	//sInfo->NN = sInfo->numNReg;
	sInfo->NN = sInfo->BaSize+sInfo->BaSizeconn[0]+sInfo->BaSizeconn[1];
	sInfo->ND = sInfo->Ndelay;
	sInfo->firingTableHost=(spikeTime_t*)malloc(sInfo->NNgroup*sizeof(spikeTime_t));
	sInfo->firingTableAll=(spikeTime_t*)malloc(sInfo->NN*sInfo->ND*sizeof(spikeTime_t));

	sInfo->displs = (int *) malloc(sizeof(int) * lcnproc); // displs array
    	sInfo->recvCount = (int *) malloc(sizeof(int) * lcnproc); // send num array
//+++++++++run mpi++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#if 1
{
    	spikeTime_t *sendBuf, *recvBuf;
    	int *displs, *recvCount;
    	int root = 0;
	int senddatanum,nsall=0;

       	sendBuf = &(sInfo->firingTableHost[0]);
       	senddatanum= sInfo->gSize;

    	displs = sInfo->displs; // displs
    	recvCount = sInfo->recvCount; // send num

 	MPI_Gather(&senddatanum,1,MPI_INT,recvCount,1, MPI_INT,root,sInfo->lccomm);

    	if(!lcrank){
		displs[0] = 0;
		nsall = recvCount[0];
        	for(i=1;i<lcnproc;i++){
            		//printf("%d \n",recvCount[i]);
        		displs[i] = displs[i-1]+recvCount[i-1];
			nsall += recvCount[i];
		}
           	//printf("mpi::%d %d %d %d NSall=%d\n",recvCount[0],recvCount[1],recvCount[2],recvCount[3],NSall);
    	}

	recvBuf = &(sInfo->firingTableAll[0]);
	
	assert(sizeof(spikeTime_t)==sizeof(long));
	MPI_Gatherv(sendBuf, 
		senddatanum, 
		MPI_LONG, 
		recvBuf, 
		recvCount, 
		displs, 
		MPI_LONG,
		root, 
		sInfo->lccomm);
#if 1	
 	MPI_Bcast(&nsall,1,MPI_INT,root,sInfo->lccomm);
	MPI_Bcast(recvBuf,nsall,MPI_LONG,root,sInfo->lccomm);
#endif 
}
#endif
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	/*******set syn and neur*******/
	int gid;
	assert(gInfo[0].StartN==0);
	for(gid=0;gid<1;gid++){
		int S = MAX(gInfo[gid].StartN, sInfo->gStart);
		int E = MIN(gInfo[gid].EndN, sInfo->gStart+sInfo->gSize-1);
		for(i=S;i<=E;i++)resetNeuron2(sInfo,i-sInfo->gStart,gInfo,gid);//nInfo
	}
	sInfo->synWt = cInfo[0].initWt;
	/******************************/
	return;
}

//void createNetwork(snnInfo_t *sInfo,grpInfo_t *gInfo,connInfo_t *cInfo){//????
void setswInfo(snnInfo_t *sInfo,grpInfo_t *gInfo,connInfo_t *cInfo){//????
	int i;
	/*****set SW*****/
	for(i=0;i<NTh;i++){
		if(sInfo->swInfo[i].StartN>=sInfo->numNReg){assert(0);}//???

		sInfo->swInfo[i].dAMPA = sInfo->dAMPA;
		sInfo->swInfo[i].dNMDA = sInfo->dNMDA;
		sInfo->swInfo[i].rNMDA = sInfo->rNMDA;
		sInfo->swInfo[i].sNMDA = sInfo->sNMDA;

		sInfo->swInfo[i].dGABAa = sInfo->dGABAa;
		sInfo->swInfo[i].dGABAb = sInfo->dGABAb;
		sInfo->swInfo[i].rGABAb = sInfo->rGABAb;
		sInfo->swInfo[i].sGABAb = sInfo->sGABAb;

		sInfo->swInfo[i].sim_with_conductances = sInfo->sim_with_conductances;
		sInfo->swInfo[i].WithSTDP = 0;

	}
	/******************************/
	return;
}

void setupNetwork(snnInfo_t *sInfo,grpInfo_t *gInfo) {
}

static void resetNeuron(neurInfo_t *nInfo,int n,grpInfo_t *gInfo,int g) {
	nInfo[n].Izh_a=gInfo[g].Izh_a+gInfo[g].Izh_a_sd*(float)drand48();
	nInfo[n].Izh_b=gInfo[g].Izh_b+gInfo[g].Izh_b_sd*(float)drand48();
	nInfo[n].Izh_c=gInfo[g].Izh_c+gInfo[g].Izh_c_sd*(float)drand48();
	nInfo[n].Izh_d=gInfo[g].Izh_d+gInfo[g].Izh_d_sd*(float)drand48();

	nInfo[n].voltage=nInfo[n].Izh_c;// initial values for new_v
	nInfo[n].recovery=nInfo[n].Izh_b*nInfo[n].voltage;//for new_u

	/*************/
	nInfo[n].gAMPA = 0.;
	nInfo[n].gNMDA_d = 0.;
	//nInfo[n].gNMDA_r = 0.;
	nInfo[n].gGABAa = 0.;
	nInfo[n].gGABAb_d = 0.;
	//nInfo[n].gGABAb_r = 0.;
	/*************/
}
static void resetNeuron2(snnInfo_t *sInfo,int n,grpInfo_t *gInfo,int g) {
	sInfo->swPara.Izh_a[n]=gInfo[g].Izh_a+gInfo[g].Izh_a_sd*(float)drand48();
	sInfo->swPara.Izh_b[n]=gInfo[g].Izh_b+gInfo[g].Izh_b_sd*(float)drand48();
	sInfo->swPara.Izh_c[n]=gInfo[g].Izh_c+gInfo[g].Izh_c_sd*(float)drand48();
	sInfo->swPara.Izh_d[n]=gInfo[g].Izh_d+gInfo[g].Izh_d_sd*(float)drand48();

	sInfo->swPara.voltage[n]=sInfo->swPara.Izh_c[n];// initial values for new_v
	//sInfo->swPara.recovery[n]=sInfo->swPara.Izh_b[n]*sInfo->swPara.voltage[n];//for new_u
	sInfo->swPara.recovery[n]=50;//for new_u

	/*************/
	sInfo->swPara.gAMPA[n] = 0.;
	sInfo->swPara.gNMDA_d[n] = 0.;
	sInfo->swPara.gGABAa[n] = 0.;
	sInfo->swPara.gGABAb_d[n] = 0.;
	/*************/
}

//resets nSpikeCnt[]
//static void resetSpikeCnt(snnInfo_t* sInfo,grpInfo_t* gInfo,int gid) {
//}
//static int updateSpikeTables(snnInfo_t *s) {//generate firingTable
//}
