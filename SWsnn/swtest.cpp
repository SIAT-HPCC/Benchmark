#include "mysnn.h"
#include <sys/time.h>    // for gettimeofday() chenged!!

static void setBrainArea(snnInfo_t *sInfo, int Nneuron, int BaNum, int Nsyn, int exNsyn);

int main(int argc, char *argv[]) {
//========================================================================
    	int rank, nproc;
    	MPI_Init(&argc, &argv);
    	MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//=======================================================================
	// create a network on CPU or GPU
	int numGrp=1, numConn=2;
	int randSeed = 42;
	int NE=10240000; int NI=0; int NP=0;
	int NReg=NE+NI; int NInput=NP; int NPre=NReg+NInput;
	//assert(NReg<=120000); //maxmum neurons
	snnInfo_t snnInfo;
	grpInfo_t *grpInfo=(grpInfo_t*)malloc(numGrp*sizeof(grpInfo_t));
	connInfo_t *connInfo=(connInfo_t*)malloc(numConn*sizeof(connInfo_t));
	///init set grpInfo connInfo snnInfo
	initNetwork(&snnInfo,grpInfo,connInfo,numGrp,numConn,randSeed);

	setBrainArea(&snnInfo, NReg, 16, 1000, 500);

/*************************/
struct timeval t0,t1,t2,t3,t4,t5,t6,t7,t8,t9;
struct timeval *time1,*time2;
gettimeofday( &t0, NULL );
/*************************/
gettimeofday( &t1, NULL );
	// configure the network
	float a,b,c,d,a_sd,b_sd,c_sd,d_sd;
	
	/*******grp para*****/
	int i=0;
	grpInfo[i].SizeN = NE;
	grpInfo[i].StartN= 0; 
	grpInfo[i].EndN = grpInfo[i].StartN+grpInfo[i].SizeN-1;
	a   =0.02; b   =0.2; c   =-65; d   =8.0;
	a_sd=0.0 ; b_sd=0.05; c_sd=0.0;  d_sd=0.0;
	setNeuronParameters(grpInfo,i,a,a_sd,b,b_sd,c,c_sd,d,d_sd);
	/*********************/
	///*******set snn para******/
	snnInfo.numGrp=i+1;
	snnInfo.numN=grpInfo[i].EndN+1;
	snnInfo.numNExcReg=grpInfo[0].SizeN;
	snnInfo.numNInhReg=0;
	snnInfo.numNReg=grpInfo[0].SizeN;
	assert(snnInfo.numNReg==NReg);
	assert(snnInfo.numNPois==NInput);
	assert(snnInfo.numNExcReg==NE);
	assert(snnInfo.numNInhReg==NI);
	assert(snnInfo.numN==NPre);
	//snnInfo.Nsyn=1000;//connection number+++++++++++++++++++++++++++
	assert(snnInfo.numNReg%snnInfo.Nsyn==0);///+++++++++++++++++++++++
	//assert(snnInfo.numNReg%nproc==0);
	
	/*********************/
	int minD,maxD; double minW,maxW,initW;
	int gSrc,gDest;
	/*******conn para******/
	//NE->NE
	int DL=2;assert(DL>=2);
	minD=1; maxD=DL; minW=0.0;initW=0.000001;maxW=0.000001;
	gSrc=0;gDest=0;
	int j=0;
	connInfo[j].connId=j;
	connInfo[j].grpSrc=gSrc;
	connInfo[j].grpDest=gDest;
	
	connInfo[j].maxDelay=maxD;
	connInfo[j].initWt=initW;
	connInfo[j].maxWt=maxW;

	connInfo[j].p=1.;
	connInfo[j].maxSynPost=grpInfo[gDest].SizeN-1;
	connInfo[j].maxSynPre=grpInfo[gSrc].SizeN-1;

	connInfo[j].numSyn=grpInfo[gSrc].SizeN*(grpInfo[gSrc].SizeN-1);
	
	/*******conn para******/
	//NP->NE
	minD=1; maxD=DL; minW=0.000;initW=0.0005;maxW=0.0005;
	gSrc=1;gDest=0;
	j++;
	connInfo[j].connId=j;
	connInfo[j].grpSrc=gSrc;
	connInfo[j].grpDest=gDest;
	
	connInfo[j].maxDelay=maxD;
	connInfo[j].initWt=initW;
	connInfo[j].maxWt=maxW;

	connInfo[j].p=1.;
	connInfo[j].maxSynPost=grpInfo[gDest].SizeN;
	connInfo[j].maxSynPre=grpInfo[gSrc].SizeN;

	connInfo[j].numSyn=grpInfo[gSrc].SizeN*grpInfo[gDest].SizeN;
	/*******************************/
	snnInfo.maxDelay=maxD;//????
gettimeofday( &t2, NULL );
gettimeofday( &t3, NULL );
gettimeofday( &t4, NULL );
	//COBA or CUBA (true or false)
	int tdAMPA=5,tdNMDA=150,tdGABAa=6,tdGABAb=150;
	setConductances(&snnInfo,1,tdAMPA,0,tdNMDA,tdGABAa,0,tdGABAb);

gettimeofday( &t5, NULL );
	// build the network
	buildModel(&snnInfo,grpInfo,connInfo);///malloc, block, and connectPart
	setswInfo(&snnInfo,grpInfo,connInfo);///set SW parameters
	
	MPI_Barrier(MPI_COMM_WORLD);//mpi syn
gettimeofday( &t6, NULL );
gettimeofday( &t7, NULL );
{gettimeofday( &t8, NULL );}
	int nmsec = 1000;
	runNetwork(&snnInfo,grpInfo,connInfo,snnInfo.nInfoHost,
               snnInfo.sInfoHost,snnInfo.swInfo,nmsec,0);
{gettimeofday( &t9, NULL );}
	freeNetwork(&snnInfo);
    	MPI_Finalize();
//==========================================================================
/****time cost*******/
#if 0
#endif
if(rank==0){
time1=&t8; time2=&t9;
printf("main::t9-t8: %f\n", (double)(time2->tv_sec-time1->tv_sec)+(double)(time2->tv_usec-time1->tv_usec)*1e-6);

time1=&t0; time2=&t9;
printf("main::all  : %f\n", (double)(time2->tv_sec-time1->tv_sec)+(double)(time2->tv_usec-time1->tv_usec)*1e-6);
}
	return 0;
}

static void setBrainArea(snnInfo_t *sInfo, int Nneuron, int BaNum,
			 int Nsyn, int exNsyn){
    	int rank, nproc;
    	MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Group world_group,local_group;
     	MPI_Comm_group(MPI_COMM_WORLD, &world_group);

        int BaId;
        int *BaStart, *BaSize, *BaEnd;
        int *BaStartP, *BaSizeP, *BaEndP;
	
	BaStart = (int*)malloc(BaNum*sizeof(int));
	BaSize = (int*)malloc(BaNum*sizeof(int));
	BaStartP = (int*)malloc(BaNum*sizeof(int));
	BaSizeP = (int*)malloc(BaNum*sizeof(int));

	
	
	int Size=nproc/BaNum;
	int least=nproc%BaNum;
	assert(BaNum<=nproc);
	for(int ib=0;ib<BaNum;ib++){
	if (ib<least){
		BaStartP[ib] = ib*(Size+1);
		BaSizeP[ib] = Size+1;
	}else{
		BaStartP[ib] = ib*Size+least;
                BaSizeP[ib] = Size;
	}
	}
	
	for(int ib=rank/(Size+1);;ib++){
		if (rank>=BaStartP[ib]&&rank<BaStartP[ib]+BaSizeP[ib])
		{ 	BaId=ib; break;}
	}
	
	Size=Nneuron/BaNum;
	least=Nneuron%BaNum;
	for(int ib=0;ib<BaNum;ib++){
	if (ib<least){
		BaStart[ib] = ib*(Size+1);
		BaSize[ib] = Size+1;
	}else{
		BaStart[ib] = ib*Size+least;
                BaSize[ib] = Size;
	}
	}

	int nrank = BaSizeP[BaId];
	int *ranks = (int*)malloc(nrank*sizeof(int));
	for(int ir=0;ir<nrank;ir++){ranks[ir]=BaStartP[BaId]+ir;}

	MPI_Group_incl(world_group, nrank, ranks, &local_group);
	MPI_Comm local_comm;
	MPI_Comm_create(MPI_COMM_WORLD, local_group, &local_comm);

	
	sInfo->BaId = BaId;
	sInfo->BaNum = BaNum;
	sInfo->BaStart = BaStart[BaId];
	sInfo->BaSize = BaSize[BaId];
	sInfo->BaEnd = BaStart[BaId]+BaSize[BaId]-1;
	sInfo->BaStartP = BaStartP[BaId];
	sInfo->BaSizeP = BaSizeP[BaId];
	sInfo->BaEndP = BaStartP[BaId]+BaSizeP[BaId]-1;
	
	sInfo->Nsyn=Nsyn;
	sInfo->BaNconn = 2;
	if(BaId==0||BaId==BaNum-1)sInfo->BaNconn = 1;
	
	sInfo->BaIdconn[0] = BaId-1;
	sInfo->BaIdconn[1] = BaId+1;
	if(BaId==0)sInfo->BaIdconn[0] = 1;
	if(BaId==BaNum-1)sInfo->BaIdconn[0] = BaId-1;

	sInfo->BaSizeconn[0]=0;
	sInfo->BaSizeconn[1]=0;
	
	for(int ic=0;ic<sInfo->BaNconn;ic++){
		int ib = sInfo->BaIdconn[ic];
		sInfo->BaProot[ic]=BaStartP[ib];
		sInfo->BaStartconn[ic]=BaStart[ib];
		sInfo->BaSizeconn[ic]=BaSize[ib];
		sInfo->BaEndconn[ic]=BaStart[ib]+BaSize[ib]-1;
		sInfo->exNsyn[ic] = exNsyn;
	}
	int root0;
	MPI_Group_translate_ranks(world_group,1,&sInfo->BaStartP,local_group,&root0);
	assert(root0==0);

	sInfo->lcgroup = local_group;
	sInfo->lccomm = local_comm;

	fprintf(stderr,"Proot::%d %d %d\n",rank,sInfo->BaProot[0],sInfo->BaProot[1]);
	
	/*
    	int lcrank=-1, lcnproc=-1;
    	MPI_Comm_size(sInfo->lccomm, &lcnproc);
    	MPI_Comm_rank(sInfo->lccomm, &lcrank);
	fprintf(stderr,"rank: %d %d np: %d %d\n",rank,lcrank,nproc,lcnproc);
	int buffer[4];
	buffer[0]=rank;buffer[1]=rank;buffer[2]=rank;buffer[3]=rank;
	MPI_Bcast(buffer,4,MPI_INT,1,sInfo->lccomm);
	fprintf(stderr,"lc_bcast::%d %d %d %d %d\n",rank,buffer[0],buffer[1],buffer[2],buffer[3]);*/
	
	free(BaStart);
	free(BaSize);
	free(BaStartP);
	free(BaSizeP);
	free(ranks);
	return;
}
