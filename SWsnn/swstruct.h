#ifndef _SWSTRUCT_H_
#define _SWSTRUCT_H_

#include <assert.h>

#define MIN(a,b) ((a)<(b)?(a):(b))
#define MAX(a,b) ((a)>(b)?(a):(b))

typedef struct spikeTime_s{
	unsigned int time;//unit: resolution interval
	unsigned int nid;
}spikeTime_t;

typedef struct synInfo_s{
	unsigned int postId;//post neuron id
	//unsigned short dl;//unit: resolution interval
	float wt;
}synInfo_t;

typedef struct neurInfo_s{
	//int nSpikeCnt;
	
	float Izh_b,Izh_a,Izh_c,Izh_d;
	
	float voltage;
	float recovery;

	float gAMPA;
	float gNMDA_d;
	float gGABAa;
	float gGABAb_d;
	int nSpikeCnt;
	
}neurInfo_t;

typedef struct swInfo_s{
	int StartN;
	int SizeN;
	float dAMPA,dNMDA,rNMDA,sNMDA;
	float dGABAa,dGABAb,rGABAb,sGABAb;

	//properties
	int sim_with_conductances;
	int WithSTDP;
	//int WithSTP;
	//int WithHomeostasis;
	//int sim_with_homeostasis;
	//int sim_with_NMDA_rise;
	//int sim_with_GABAb_rise;
	//int MaxDelay;
	
	//connectInfo
	//int preN; //num of pre neurons (3rd dim)
	//int Ndelay;
	int MaxN; //max num of connections in a Th for a preN (1st dim)
	int Ndma;

	int Nop,Ndt;
	float dt;
	int sliceTime, simTime;
	//int nSpikePoisAll;

	//int fireCnt;

	//host address
	//synInfo_t *sInfoHost;
	//neurInfo_t *nInfoHost;
	//int *firingTableHost;
	
	//mpi firingTable
	//spikeTime_t *firingTableHost;
	//spikeTime_t *firingTableAll;
	//int NNgroup,NN,ND;

	//int rank,nproc;
	//int gStart,gSize;
	
}swInfo_t;

typedef struct swPara_s{
	//int nSpikeCnt;
	// neuron data arrays
	float *Izh_b;
	float *Izh_a;
	float *Izh_c;
	float *Izh_d;
	
	float *voltage;
	float *recovery;

	float *gAMPA;
	float *gNMDA_d;
	float *gGABAa;
	float *gGABAb_d;

	// synapse data arrays
	float *wt;
	unsigned int *postId;//post neuron id
	//unsigned short *dl;//unit: resolution interval
	float synWt;	

	// Brain area
	float *exwt[2];
        int   *expostId[2];

	// firing table
	spikeTime_t *firingTable;
	spikeTime_t *firingTableAll;
	spikeTime_t *firingTableHost;
	spikeTime_t *firingTableAllHost;
	
	// ring buffer
	float *ringBuffer;

	//global variables
	int gStart;
	int gSize;
	int offset;//ringBuffer offset
	int lenRB; //ringBuffer length

	int endST; //firingTable location
	
	int preN;
	int MaxDelay;
	int MaxN;
	int Ndelay;
	int nSpikeCnt;

	int Dconn;

	float *extrvol;	

	int BaSize,BaStart;
		
}swPara_t;








#endif
