#ifndef _MY_SNN_H_
#define _MY_SNN_H_
#include <stdint.h>
#include <stdlib.h> 
#include <stdio.h>
//#include <athread.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <assert.h>
#include <sys/time.h>    // for gettimeofday() chenged!!
#include "mpi.h"

#include "mypbuf.h"

// typedef int bool;
 
#include "swstruct.h"
#define NTh 64

/*
static inline unsigned long rpcc()
{
        unsigned long time;
        asm("rtc %0": "=r" (time) : );
        return time;
}
*/

typedef struct PoissonRate_s{
	int numPois;
	float *h_rates_;

}PoissonRate_t;
typedef struct snnInfo_s{

	unsigned int	simTimeMs;
	uint64_t        simTimeSec;
	unsigned int	simTime;
	uint8_t                                  maxDelay;
	//int				maxDelay_;
	unsigned int    numGrp,numConn;	
	int	        	numN,numNReg,numNExcReg,numNInhReg,numNPois;
	unsigned int		maxSpikesD1;
	unsigned int		maxSpikesD2;
	unsigned int    simTimeRunStart;
	unsigned int    simTimeRunStop;

	double dAMPA,dNMDA,rNMDA,sNMDA;
	double dGABAa,dGABAb,rGABAb,sGABAb;

	unsigned int		*timeTableD1;
	unsigned int		*timeTableD2;
	unsigned int		*firingTableD1;
	unsigned int		*firingTableD2;
	unsigned int    MaxFiringRate;
	unsigned int	secD1fireCntHost;
	unsigned int	secD2fireCntHost;
 	unsigned int	spikeCountAll1secHost;
	Pbuf_t 		*pbuf;
	//short int *grpIds;
	bool sim_with_homeostasis;
	bool sim_with_conductances;
	//bool sim_with_homeostasis;
	bool sim_with_NMDA_rise;
	bool sim_with_GABAb_rise;
	int         	*nSpikeCnt;
	uint32_t    	*lastSpikeTime;

	/*****for SW******/
	swInfo_t swInfo[NTh];
	swInfo_t *swInfoDevice;
	
	swPara_t swParaHost;
	swPara_t *swParaDevice;
	swPara_t swPara;
	//swInfo_t *swInfo;
	int preN,MaxN;
	synInfo_t* sInfoHost;
	neurInfo_t* nInfoHost;
	synInfo_t* sInfoDevice;
	neurInfo_t* nInfoDevice;
	/*****************/
	int Nsyn;
	int Ndelay;
	//float DT;//mindelay
	float dt;
	int Ndt, Nop;
	/*********/
	int randSeed;

	//mpi firingTable
	spikeTime_t *firingTableHost;
	spikeTime_t *firingTableAll;
	spikeTime_t *firingTableDevice;
	spikeTime_t *firingTableAllDevice;
	int NNgroup,NN,ND;
	int gStart,gSize;
	int *displs, *recvCount;


	float synWt;

	//////timer
	double DE_T0,DE_T1,DE_T2,DE_T3,DE_T4;
	
	float *extrvol;

	//////Brain area
	int BaId,BaNum;
        int BaStart, BaSize,BaEnd;
        int BaStartP, BaSizeP, BaEndP;
        int BaNconn;
        int BaIdconn[2];
        int BaProot[2];
        int BaStartconn[2], BaSizeconn[2],BaEndconn[2];
        int exNsyn[2];

	MPI_Group lcgroup;
	MPI_Comm lccomm;

}snnInfo_t;

typedef struct grpInfo_s{
	int StartN;
	int EndN;
	int SizeN;
	int8_t		MaxDelay;
	int 		FiringCount1sec;

	bool	isSpikeGenerator;
	int			NewTimeSlice;
	int			CurrTimeSlice;
	uint32_t 	SliceUpdateTime;
	float   	RefractPeriod;
	PoissonRate_t*	RatePtr;

	bool	WithSTP;
	bool 		WithHomeostasis;
	unsigned int	Type;

	unsigned int numSynPre,numSynPost;

	float Izh_a,Izh_a_sd;
	float Izh_b,Izh_b_sd;
	float Izh_c,Izh_c_sd;
	float Izh_d,Izh_d_sd;

}grpInfo_t;

typedef struct connInfo_s{
	uint8_t connId;
	uint8_t maxDelay;
	int grpSrc,grpDest;
	float initWt,maxWt;
	float mulSynFast,mulSynSlow;
	int maxSynPost,maxSynPre;
	int numSyn;
	//conType_t type;
	uint32_t connProp;//bit info & op
	float p;//conn probability
}connInfo_t;

#if 1
int runNetwork(snnInfo_t *snnInfo, grpInfo_t *grpInfo,
               connInfo_t *connInfo, neurInfo_t *neurInfo,
               synInfo_t *synInfo_t, swInfo_t *swInfo,
               int _nmsec, bool printRun);

void initNetwork(snnInfo_t *sInfo, grpInfo_t *gInfo,connInfo_t *cInfo,
                int numGrp,int numConn,int randSeed);
void setNeuronParameters(grpInfo_t *gInfo,int gid,float izh_a,
		float izh_a_sd,float izh_b,float izh_b_sd,float izh_c,
		float izh_c_sd,float izh_d,float izh_d_sd);
void setConductances(snnInfo_t *sInfo,bool isSet,int tdAMPA,int trNMDA,
		int tdNMDA,int tdGABAa,int trGABAb,int tdGABAb);
void createNetwork(snnInfo_t *sInfo,grpInfo_t *gInfo,connInfo_t *cInfo);
void setswInfo(snnInfo_t *sInfo,grpInfo_t *gInfo,connInfo_t *cInfo);
void buildModel(snnInfo_t *sInfo,grpInfo_t *gInfo,connInfo_t *cInfo);
void setupNetwork(snnInfo_t *sInfo,grpInfo_t *gInfo);
void setRates(PoissonRate_t *ratePtr,int numPois,float rate);
void setSpikeRate(grpInfo_t *gInfo,int gid,PoissonRate_t* rPtr,int refPeriod);
void freeNetwork(snnInfo_t *sInfo);
#endif


#if 0
typedef struct synInfo_s{
	unsigned int postId;//post neuron id
	float wt;
	//float wtChange;
}synInfo_t;

typedef struct neurInfo_s{
	
	
}neurInfo_t;

typedef struct swInfo_s{
	int StartN;
	int SizeN;
	//int EndN;
	//int GrpId;


	//properties
	int WithSTDP;
	
	//connectInfo
	int preN; //num of pre neurons (3rd dim)
	//int NTh=64; //num of threads (2nd dim)
	int maxConnPT; //max num of connections in a Th for a preN (1st dim)

	//related pointers
	synInfo_t *synInfo;
	int *firingTableD1;
	
}swInfo_t;
#endif
#endif


