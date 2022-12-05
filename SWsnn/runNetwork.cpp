
#include "mysnn.h"
#include "mygpu.h"

#define true 1
#define false 0
#define PROPAGATED_BUFFER_SIZE  (1023)
#define ALL -1
#define MAX_SIMULATION_TIME     ((uint32_t)(0x7fffffff))

#define UNKNOWN_NEURON  (0)
#define POISSON_NEURON  (1 << 0)
#define TARGET_AMPA     (1 << 1)
#define TARGET_NMDA     (1 << 2)
#define TARGET_GABAa    (1 << 3)
#define TARGET_GABAb    (1 << 4)

unsigned int dma[64],spike[64],NS_group,NSall[3],numSpike[64],nspikeall;
/******************************/
static bool updateTime(snnInfo_t *snnInfo);
static void doSnnSim(snnInfo_t*,grpInfo_t*,connInfo_t*,neurInfo_t*,synInfo_t*, swInfo_t*);
/********global**********/
double DO_T1,DO_T2,DO_T3,DO_T4,DO_T5,DO_T6,DO_T7,DO_T8,DO_TA;
long CC0,CC1;
/*************************/
int runNetwork(snnInfo_t *sInfo, grpInfo_t *gInfo,
	       connInfo_t *cInfo, neurInfo_t *nInfo,
	       synInfo_t *synInfo, swInfo_t *swInfo,
	       int _nmsec, bool printRun) {
/*************************/
struct timeval t0,t1,t2,t3,t4,t5,t6,t7,t8,t9;
struct timeval *time1,*time2;
/*************************/
	//mpi function
    	int rank, nproc;
    	MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	assert(_nmsec >= 0);
	int runDurMs = _nmsec;
	sInfo->simTimeRunStart = sInfo->simTime;
	sInfo->simTimeRunStop  = sInfo->simTime+runDurMs;
/******************************/
DO_T1=0.0; DO_T2=0.0; DO_T3=0.0; DO_T4=0.0;
DO_T5=0.0; DO_T6=0.0; DO_T7=0.0; DO_T8=0.0; DO_TA=0.0;CC0=0;CC1=0;

sInfo->DE_T0=0.0;
sInfo->DE_T1=0.0;
sInfo->DE_T2=0.0;
sInfo->DE_T3=0.0;
sInfo->DE_T4=0.0;
/********************************/
{gettimeofday( &t0, NULL );}
	initSW(sInfo,swInfo);  //GPU init**************************************
{gettimeofday( &t1, NULL );}
	int i;
//printf("0000\n");
	nspikeall=0;

if(rank==0)printf("**************start simulation***************\n");	
	for(i=0; i<runDurMs; i++) {
if(rank==0&&(i&0x0000001f)==0)printf("->time %d ms\n",i);	
		doSnnSim(sInfo,gInfo,cInfo,nInfo,synInfo,swInfo);
		updateTime(sInfo);
	}
//printf("0001\n");
	freeSW(sInfo); //GPU free************************************************
/**********************/
unsigned long cdma=0,cspike=0;
for(i=0;i<63;i++){
	cdma=MAX(dma[i],cdma);
	cspike=MAX(spike[i],cspike);
}
int nSpikeAll=0;
for(i=0;i<64;i++){
	nSpikeAll+=numSpike[i];
}
printf("rank=%d nspikeall=%d nSpikeAll=%d cdma=%ld,cspike=%ld\n",rank,nspikeall,nSpikeAll,cdma,cspike);
/**********************/
if(rank==0){
printf("doSnnSim::T5: %f s [StateUpdate]\n", DO_T5);
printf("doSnnSim::T6: %f s [mpi communication]\n", DO_T6);
printf("doSnnSim::T7: %f s [SpikeDeliver]\n", DO_T7);

printf("doSnnSim::T1: %f s [StateUpdate-Kernel]\n", sInfo->DE_T1);
printf("doSnnSim::T2: %f s [StateUpdate-DtoH]\n", sInfo->DE_T2);
printf("doSnnSim::T3: %f s [SpikeDeliver-HtoD]\n", sInfo->DE_T3);
printf("doSnnSim::T4: %f s [SpikeDeliver-Kernel]\n", sInfo->DE_T4);
}
/************************************************/
if(rank==0){
time1=&t0; time2=&t1;
printf("initSW::t1-t0: %f s [initSW cost]\n", (double)(time2->tv_sec-time1->tv_sec)+(double)(time2->tv_usec-time1->tv_usec)*1e-6);
printf("snnSim::T5+T6+T7: %f s [Simulate cost]\n", DO_T5+DO_T6+DO_T7);
}

//SNN模拟总耗时、神经元更新耗时、MPI通信耗时、脉冲传送耗时、HtoD数据拷贝耗时、DtoH数据拷贝耗时
if(rank==0){
FILE *fout = fopen("b.out","w");
fprintf(fout,"SNN模拟总耗时: %f s\n",DO_T5+DO_T6+DO_T7);
fprintf(fout,"神经元更新耗时(Kernel): %f s\n",sInfo->DE_T1);
fprintf(fout,"脉冲传送耗时(Kernel): %f s\n",sInfo->DE_T4);
fprintf(fout,"MPI通信耗时: %f s\n",DO_T6);
fprintf(fout,"数据拷贝耗时(HtoD): %f s\n",sInfo->DE_T3);
fprintf(fout,"数据拷贝耗时(DtoH): %f s\n",sInfo->DE_T2);
fclose(fout);
}
	return 0;
}

static void doSnnSim(snnInfo_t *sInfo, grpInfo_t *gInfo,
               connInfo_t *cInfo, neurInfo_t *nInfo,
               synInfo_t *synInfo, swInfo_t *swInfo) {
/*************************/
struct timeval t0,t1,t2,t3,t4,t5,t6,t7,t8,t9;
struct timeval *tm1,*tm2;
gettimeofday( &t0, NULL );
/*************************/
long st,st0,ed;
//printf("d0000\n");
{{gettimeofday( &t4, NULL );}
	StateUpdate(sInfo); //GPU running***************************************
{gettimeofday( &t5, NULL );}}
	sInfo->simTime++;
//printf("d0001\n");
	//mpi function?????
    	int lcrank, lcnproc;
    	MPI_Comm_size(sInfo->lccomm, &lcnproc);
    	MPI_Comm_rank(sInfo->lccomm, &lcrank);
    	int rank, nproc;
    	MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	int i;
	int iND = sInfo->simTime%sInfo->Ndelay;
#if 1
//printf("d0002\n");
    	spikeTime_t *sendBuf, *recvBuf;
    	int *displs, *recvCount;
    	int root = 0;
	int senddatanum;

       	//sendBuf = &(sInfo->firingTableHost[0]);
       	sendBuf = &(sInfo->swParaHost.firingTableHost[0]);
       	//senddatanum= NS_group;
       	senddatanum= (sInfo->gSize-1)/100+1;///???constant for mpi
       	senddatanum= rand()%(senddatanum*2);///???constant for mpi

    	displs = sInfo->displs; // displs
    	recvCount = sInfo->recvCount; // send num

//printf("d0003\n");
 	MPI_Gather(&senddatanum,1,MPI_INT,recvCount,1, MPI_INT,root,sInfo->lccomm);

//printf("d0004\n");
	NSall[0]=0;NSall[1]=0;NSall[2]=0;
    	if(!lcrank){
		displs[0] = 0;
		NSall[0] = recvCount[0];
        	for(i=1;i<lcnproc;i++){
            		//printf("%d \n",recvCount[i]);
        		displs[i] = displs[i-1]+recvCount[i-1];
			NSall[0] += recvCount[i];
		}
    	}
//printf("d0005\n");
	//recvBuf = &(sInfo->firingTableAll[iND*sInfo->NN]);
	recvBuf = &(sInfo->swParaHost.firingTableAllHost[iND*sInfo->NN]);
	MPI_Gatherv(sendBuf, 
		senddatanum, 
		MPI_LONG, 
		recvBuf, 
		recvCount, 
		displs, 
		MPI_LONG,
		root, 
		sInfo->lccomm);




//++++++++brain area comm+++++++++++++++++++++++++++++++++++++++++
if(lcrank==0){
	MPI_Request request0[2],request1[2];
	MPI_Request request2[2],request3[2];
	MPI_Status status2[2],status3[2];
	MPI_Status status0[2],status1[2];
	for(int ic=0;ic<sInfo->BaNconn;ic++){
		int dest = sInfo->BaProot[ic];
		MPI_Isend(&NSall[0],1,MPI_INT,dest,dest,MPI_COMM_WORLD,&request0[ic]);
		MPI_Irecv(&NSall[ic+1],1,MPI_INT,dest,rank,MPI_COMM_WORLD,&request1[ic]);
	}
	int offset=0;
	for(int ic=0;ic<sInfo->BaNconn;ic++){
		MPI_Wait(&request0[ic],&status0[ic]);
		MPI_Wait(&request1[ic],&status1[ic]);
		int dest = sInfo->BaProot[ic];
		MPI_Isend(&(sInfo->swParaHost.firingTableAllHost[iND*sInfo->NN]),
			NSall[0],MPI_LONG,dest,dest,MPI_COMM_WORLD,&request2[ic]);
		offset += NSall[ic];
		MPI_Irecv(&(sInfo->swParaHost.firingTableAllHost[iND*sInfo->NN+offset]),
			NSall[ic+1],MPI_LONG,dest,rank,MPI_COMM_WORLD,&request3[ic]);
	}
	for(int ic=0;ic<sInfo->BaNconn;ic++){
		MPI_Wait(&request2[ic],&status2[ic]);
		MPI_Wait(&request3[ic],&status3[ic]);
	}
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


//printf("d0006\n");
#if 1	
 	MPI_Bcast(&NSall[0],3,MPI_INT,root,sInfo->lccomm);
	MPI_Bcast(recvBuf,NSall[0]+NSall[1]+NSall[2],MPI_LONG,root,sInfo->lccomm);
	nspikeall += NSall[0];
	nspikeall += NSall[1];
	nspikeall += NSall[2];
//printf("d0007\n");
#endif 
#endif
gettimeofday( &t6, NULL );
{	
	SpikeDeliver(sInfo);
{gettimeofday( &t7, NULL );
}}
tm1=&t4; tm2=&t5;
DO_T5+=((double)(tm2->tv_sec-tm1->tv_sec)+(double)(tm2->tv_usec-tm1->tv_usec)*1e-6);
tm1=&t5; tm2=&t6;
DO_T6+=((double)(tm2->tv_sec-tm1->tv_sec)+(double)(tm2->tv_usec-tm1->tv_usec)*1e-6);
tm1=&t6; tm2=&t7;
DO_T7+=((double)(tm2->tv_sec-tm1->tv_sec)+(double)(tm2->tv_usec-tm1->tv_usec)*1e-6);
        return;
}
static bool updateTime(snnInfo_t *snnInfo) {
	bool finishedOneSec = false;
	// update relevant parameters...now
	if(++snnInfo->simTimeMs == 1000) {
		snnInfo->simTimeMs = 0;
		snnInfo->simTimeSec++;
		finishedOneSec = true;
	}
	//snnInfo->simTime++;

	return finishedOneSec;
}

