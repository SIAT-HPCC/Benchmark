#include "hip/hip_runtime.h"
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include "swstruct.h"
#include "mygpu.h"

// 检查显卡错误
void checkError() {
	hipError_t err = hipGetLastError();
    	if (err != hipSuccess) {
        	std::cout << hipGetErrorString(err) << std::endl;
    	}
}

__global__ void kernel_connectPart(swInfo_t *swInfo,swPara_t *sw0);

void initSW(snnInfo_t *sInfo, swInfo_t *swInfo){//????
    	int rank, nproc;
    	MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	///set sInfo
	
	
	///allocate host memory
	//neuron data arrays
	/*sInfo->swPara.Izh_a = (float*)malloc(sInfo->gSize*sizeof(float));
	sInfo->swPara.Izh_b = (float*)malloc(sInfo->gSize*sizeof(float));
	sInfo->swPara.Izh_c = (float*)malloc(sInfo->gSize*sizeof(float));
	sInfo->swPara.Izh_d = (float*)malloc(sInfo->gSize*sizeof(float));

	sInfo->swPara.voltage = (float*)malloc(sInfo->gSize*sizeof(float));
	sInfo->swPara.recovery = (float*)malloc(sInfo->gSize*sizeof(float));
	sInfo->swPara.gAMPA = (float*)malloc(sInfo->gSize*sizeof(float));
	sInfo->swPara.gNMDA_d = (float*)malloc(sInfo->gSize*sizeof(float));
	sInfo->swPara.gGABAa = (float*)malloc(sInfo->gSize*sizeof(float));
	sInfo->swPara.gGABAb_d = (float*)malloc(sInfo->gSize*sizeof(float));
	*/
	//synapse data arrays
	/*sInfo->swPara.wt = (float*)malloc((long)sInfo->preN*sInfo->Ndelay*sInfo->MaxN*NTh*sizeof(float));
	sInfo->swPara.postId = (unsigned int*)malloc((long)sInfo->preN*sInfo->Ndelay*sInfo->MaxN*NTh*sizeof(unsigned int));
	*///move to buildModel function
	//firingTable and ringBuffer in device
if(rank==0)printf("initSW::set swPara and swInfo\n");	
	///set swPara and swInfo
	sInfo->swParaHost.gStart = sInfo->gStart;
	sInfo->swParaHost.gSize = sInfo->gSize;
	sInfo->swParaHost.offset = 0;
	sInfo->swParaHost.lenRB = 2*sInfo->Ndt;
	sInfo->swParaHost.endST = 0;
	sInfo->swParaHost.MaxDelay = sInfo->maxDelay;
	sInfo->swParaHost.preN = sInfo->preN;
	sInfo->swParaHost.MaxN = sInfo->MaxN;
	sInfo->swParaHost.Ndelay = sInfo->Ndelay;
	sInfo->swParaHost.nSpikeCnt = 0;
	sInfo->swParaHost.Dconn = (sInfo->numNReg-1)/sInfo->Nsyn+1;
	sInfo->swParaHost.synWt = sInfo->synWt;

	sInfo->swParaHost.BaSize = sInfo->BaSize;
	sInfo->swParaHost.BaStart = sInfo->BaStart;
	
	
if(rank==0)printf("initSW::allocate device memory\n");	
	///allocate device memory
	//neuron data arrays
  	hipMalloc(&(sInfo->swParaHost.Izh_a), sInfo->gSize*sizeof(float));  
  	hipMalloc(&(sInfo->swParaHost.Izh_b), sInfo->gSize*sizeof(float));  
  	hipMalloc(&(sInfo->swParaHost.Izh_c), sInfo->gSize*sizeof(float));  
  	hipMalloc(&(sInfo->swParaHost.Izh_d), sInfo->gSize*sizeof(float));  

  	hipMalloc(&(sInfo->swParaHost.voltage), sInfo->gSize*sizeof(float));  
  	hipMalloc(&(sInfo->swParaHost.recovery), sInfo->gSize*sizeof(float));  
  	hipMalloc(&(sInfo->swParaHost.gAMPA), sInfo->gSize*sizeof(float));  
  	hipMalloc(&(sInfo->swParaHost.gNMDA_d), sInfo->gSize*sizeof(float));  
  	hipMalloc(&(sInfo->swParaHost.gGABAa), sInfo->gSize*sizeof(float));  
  	hipMalloc(&(sInfo->swParaHost.gGABAb_d), sInfo->gSize*sizeof(float));  
//printf("i0000\n");
	/****************extract voltages***************************/
	//if(rank==0){
	int NN = 10;
	sInfo->extrvol = (float*)malloc(NN*8008*sizeof(float));
	hipMalloc(&(sInfo->swParaHost.extrvol), NN*8008*sizeof(float));
	//}
	
if(rank==0)printf("initSW::copy neuron data from host to device\n");	
  	hipMemcpy(sInfo->swParaHost.Izh_a,sInfo->swPara.Izh_a,sInfo->gSize*sizeof(float), hipMemcpyHostToDevice);
  	hipMemcpy(sInfo->swParaHost.Izh_b,sInfo->swPara.Izh_b,sInfo->gSize*sizeof(float), hipMemcpyHostToDevice);
  	hipMemcpy(sInfo->swParaHost.Izh_c,sInfo->swPara.Izh_c,sInfo->gSize*sizeof(float), hipMemcpyHostToDevice);
  	hipMemcpy(sInfo->swParaHost.Izh_d,sInfo->swPara.Izh_d,sInfo->gSize*sizeof(float), hipMemcpyHostToDevice);
  	hipMemcpy(sInfo->swParaHost.voltage,sInfo->swPara.voltage,sInfo->gSize*sizeof(float), hipMemcpyHostToDevice);
  	hipMemcpy(sInfo->swParaHost.recovery,sInfo->swPara.recovery,sInfo->gSize*sizeof(float), hipMemcpyHostToDevice);
  	hipMemcpy(sInfo->swParaHost.gAMPA,sInfo->swPara.gAMPA,sInfo->gSize*sizeof(float), hipMemcpyHostToDevice);
  	hipMemcpy(sInfo->swParaHost.gNMDA_d,sInfo->swPara.gNMDA_d,sInfo->gSize*sizeof(float), hipMemcpyHostToDevice);
  	hipMemcpy(sInfo->swParaHost.gGABAa,sInfo->swPara.gGABAa,sInfo->gSize*sizeof(float), hipMemcpyHostToDevice);
  	hipMemcpy(sInfo->swParaHost.gGABAb_d,sInfo->swPara.gGABAb_d,sInfo->gSize*sizeof(float), hipMemcpyHostToDevice);

//printf("i0001\n");
	//synapse data arrays
  	//hipMalloc(&(sInfo->swParaHost.wt),sInfo->preN*sInfo->Ndelay*sInfo->MaxN*NTh*sizeof(float));
  	//hipMalloc(&(sInfo->swParaHost.postId),sInfo->preN*sInfo->Ndelay*sInfo->MaxN*NTh*sizeof(int));
  	hipMalloc(&(sInfo->swParaHost.wt),sInfo->BaSize*sInfo->Ndelay*sInfo->MaxN*sizeof(float));
  	hipMalloc(&(sInfo->swParaHost.postId),sInfo->BaSize*sInfo->Ndelay*sInfo->MaxN*sizeof(int));

//printf("i0002\n");
  	///set neurInfo in device
  	///set synInfo in device
  	///set firing table in device
  	hipHostMalloc(&sInfo->swParaHost.firingTableHost, sInfo->gSize*sizeof(spikeTime_t)); 
  	hipMalloc(&sInfo->swParaHost.firingTable, sInfo->gSize*sizeof(spikeTime_t));
	
  	///set firing table all in device
  	hipHostMalloc(&sInfo->swParaHost.firingTableAllHost,sInfo->NN*sInfo->Ndelay*sizeof(spikeTime_t));
  	hipMalloc(&sInfo->swParaHost.firingTableAll,sInfo->NN*sInfo->Ndelay*sizeof(spikeTime_t));

	///set ringBuffer
  	hipMalloc(&sInfo->swParaHost.ringBuffer,sInfo->swParaHost.lenRB*sInfo->gSize*sizeof(float));
	

if(rank==0)printf("initSW::set parameters in device\n");	
  	///set parameters in device
  	hipMalloc(&(sInfo->swInfoDevice), NTh*sizeof(swInfo_t)); 
  	hipMemcpy(sInfo->swInfoDevice, sInfo->swInfo,NTh*sizeof(swInfo_t), hipMemcpyHostToDevice);
  	hipMalloc(&(sInfo->swParaDevice), sizeof(swPara_t)); 
  	hipMemcpy(sInfo->swParaDevice, &(sInfo->swParaHost),sizeof(swPara_t), hipMemcpyHostToDevice);
	///build models in device

  	hipLaunchKernelGGL(kernel_connectPart, dim3(NBLOCK), dim3(NTHREAD), 0, 0, sInfo->swInfoDevice, sInfo->swParaDevice);
  	hipDeviceSynchronize();  

  	checkError();  // 检测错误
	return;	
}

__global__ void kernel_connectPart(swInfo_t *swInfo,swPara_t *sw0) {
printf("c000\n");
	int i,j;
	j = threadIdx.x;
	int iB = blockIdx.x;
	//int preN = sw0->preN;
	int preN = sw0->BaSize;
	int MaxN = sw0->MaxN;
	int Ndelay = sw0->Ndelay;
	int Dconn = sw0->Dconn;
	float *wt = sw0->wt;
	unsigned int *postId = sw0->postId;
	int StartN = swInfo->StartN;
	int EndN = swInfo->StartN+swInfo->SizeN-1;
	float synWt = sw0->synWt;
	///init
printf("c001\n");
	for(int iN=0;iN<preN;iN++){
		int DD = iN % Dconn;
		for(int iD=0;iD<Ndelay;iD++){
			//int offset = iN*Ndelay*NTh*MaxN+iD*NTh*MaxN+iB*MaxN;
			int offset = iN*Ndelay*MaxN;
			int S = StartN/Dconn*Dconn+DD;
			if (S<StartN) S+=Dconn;
			for (i=0; i+j<MaxN; i+=blockDim.x){
				postId[offset+i+j] = 0xffffffff;
				if(S+(i+j)*Dconn<=EndN){
					postId[offset+i+j] = S+(i+j)*Dconn;
					wt    [offset+i+j] = synWt;
				}
			}
			
		}
	}

}
void freeSW(snnInfo_t *sInfo){
    	int rank, nproc;
    	MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	if(rank==0){
  	hipMemcpy(sInfo->extrvol, sInfo->swParaHost.extrvol,8000*10*sizeof(float), hipMemcpyDeviceToHost);  
	FILE *fout = fopen("a.out","w");
	for(int i=0;i<80000;i++)fprintf(fout,"%e\n",sInfo->extrvol[i]);
	fclose(fout);}


  	return;
}
__global__ void kernel_StateUpdate(swInfo_t *swInfo,swPara_t *sw0);
__global__ void kernel_SpikeDeliver(swInfo_t *swInfo,swPara_t *sw0);

void StateUpdate(snnInfo_t *sInfo){
	struct timeval t0,t1,t2,t3,t4,t5,t6,t7,t8,t9;
	struct timeval *tm1,*tm2;
//printf("s0000\n");
	{gettimeofday( &t0, NULL );}
  	hipLaunchKernelGGL(kernel_StateUpdate, dim3(NBLOCK), dim3(NTHREAD), 0, 0, sInfo->swInfoDevice, sInfo->swParaDevice);
//printf("s0001\n");
  	hipDeviceSynchronize();  
//printf("s0002\n");
	//DtoH脉冲事件显存到主存传递
	{gettimeofday( &t1, NULL );}
  	hipMemcpy(sInfo->swParaHost.firingTableHost, sInfo->swParaHost.firingTable,
    		((sInfo->gSize-1)/100+1)*sizeof(spikeTime_t), hipMemcpyDeviceToHost);  // 拷贝数据
//printf("s0003\n");

	{gettimeofday( &t2, NULL );}
	tm1=&t0; tm2=&t1;
	sInfo->DE_T1+=((double)(tm2->tv_sec-tm1->tv_sec)+(double)(tm2->tv_usec-tm1->tv_usec)*1e-6);
	tm1=&t1; tm2=&t2;
	sInfo->DE_T2+=((double)(tm2->tv_sec-tm1->tv_sec)+(double)(tm2->tv_usec-tm1->tv_usec)*1e-6);
//printf("s0004\n");
	return;
}
void SpikeDeliver(snnInfo_t *sInfo){
	struct timeval t0,t1,t2,t3,t4,t5,t6,t7,t8,t9;
	struct timeval *tm1,*tm2;

//printf("sd000\n");
	{gettimeofday( &t0, NULL );}
  	hipMemcpy(sInfo->swParaHost.firingTableAll, sInfo->swParaHost.firingTableAllHost,
    	((sInfo->NN-1)/100+1)*sizeof(spikeTime_t), hipMemcpyHostToDevice);  // 拷贝数据
	{gettimeofday( &t1, NULL );}
	//HtoD脉冲事件主存到显存传递
	//kernel 函数启动一次
//printf("sd001\n");
  	hipLaunchKernelGGL(kernel_SpikeDeliver, dim3(NBLOCK), dim3(NTHREAD), 0, 0, sInfo->swInfoDevice, sInfo->swParaDevice);
//printf("sd002\n");
  	hipDeviceSynchronize();  // 等待核函数运行结束
//printf("sd003\n");
	{gettimeofday( &t2, NULL );}
	tm1=&t0; tm2=&t1;
	sInfo->DE_T3+=((double)(tm2->tv_sec-tm1->tv_sec)+(double)(tm2->tv_usec-tm1->tv_usec)*1e-6);
	tm1=&t1; tm2=&t2;
	sInfo->DE_T4+=((double)(tm2->tv_sec-tm1->tv_sec)+(double)(tm2->tv_usec-tm1->tv_usec)*1e-6);
	return;
}

#define dvdtIzh(v,u,tmpI,h) (((0.04*(v)+5.0)*(v)+140.0-(u)+(tmpI))*(h))
#define dudtIzh(v,u,a,b,h) ((a)*((b)*(v)-(u))*(h))

// kernel_StateUpdate
__global__ void kernel_StateUpdate(swInfo_t *swInfo,swPara_t *sw0) {
//printf("g000\n");
	if(blockIdx.x==0&&threadIdx.x==0)sw0->endST = 0;
	for(;sw0->endST!=0;){}
	__shared__ swInfo_t sw;
	//if(threadIdx.x==0) sw = swInfo[blockIdx.x];
	sw = swInfo[blockIdx.x];
	///decayConduct
//printf("g111\n");
	{
	int i,j;
	j = threadIdx.x;
	for (i=0; i+j<sw.SizeN; i+=blockDim.x){
		int offset = sw.StartN - sw0->gStart + i + j;
		//int offset = i + j;
		if(sw.sim_with_conductances){
			sw0->gAMPA[offset] *= sw.dAMPA;
                        sw0->gGABAa[offset] *= sw.dGABAa;

                        sw0->gNMDA_d[offset] *= sw.dNMDA;//instantaneous rise
                        sw0->gGABAb_d[offset] *= sw.dGABAb;//instantaneous rise
		}
		else {
			sw0->gAMPA[offset]=0.0f; //in CUBA current,sum up all wts
		}
		
	}
	}
	//int endST_mpi=0; int topST_mpi=0; int usedST_mpi=0;
	
printf("g222\n");
	for(int it=0;it<sw.Ndt;it++){
		///neuronUpdate
		{
		float tmpiNMDA, tmpI;
		float tmpgNMDA, tmpgGABAb;
		int i=0,j;
		j = threadIdx.x;
		/////for(i=0;i<swInfo.SizeN;i++) 
		for (i=0; i+j<sw.SizeN; i+=blockDim.x){
			int offset = sw.StartN - sw0->gStart + i + j;
			//int offset = i + j;
#if 0
			if (sw.sim_with_conductances) {
				float volt = sw.voltage[offset];
				tmpiNMDA=(volt+80.0)*(volt+80.0)/60.0/60.0;

				tmpgNMDA=sw.gNMDA_d[offset];
				tmpgGABAb=sw.gGABAb_d[offset];

				tmpI=-(sw.gAMPA[offset]*(volt-0)
				 +tmpgNMDA*tmpiNMDA/(1+tmpiNMDA)*(volt-0)
				 +sw.gGABAa[offset]*(volt+70)
				 +tmpgGABAb*(volt+90)   );
			} else {
				tmpI=sw.gAMPA[offset];
			}
#else

			if (sw.sim_with_conductances) {
				float volt = -60.0;
				tmpiNMDA=(volt+80.0)*(volt+80.0)/60.0/60.0;

				tmpgNMDA=sw0->gNMDA_d[offset];
				tmpgGABAb=sw0->gGABAb_d[offset];

				tmpI=-(sw0->gAMPA[offset]*(volt-0)
				 	+tmpgNMDA*tmpiNMDA/(1+tmpiNMDA)*(volt-0)
				 	+sw0->gGABAa[offset]*(volt+70)
				 	+tmpgGABAb*(volt+90)   );
			} else {
				tmpI=sw0->gAMPA[offset];
			}


#endif	


			float v = sw0->voltage[offset];
			float u = sw0->recovery[offset];
			float h = sw.dt;
			float a = sw0->Izh_a[offset];
			float b = sw0->Izh_b[offset];

			float k1 = dvdtIzh(v, u, tmpI, h);
			float l1 = dudtIzh(v, u, a, b, h);

			float k2 = dvdtIzh(v+0.5*k1, u+0.5*l1, tmpI, h);
			float l2 = dudtIzh(v+0.5*k1, u+0.5*l1, a, b, h);

			float k3 = dvdtIzh(v+0.5*k2, u+0.5*l2, tmpI, h);
			float l3 = dudtIzh(v+0.5*k2, u+0.5*l2, a, b, h);

			float k4 = dvdtIzh(v+k3, u+l3, tmpI, h);
			float l4 = dudtIzh(v+k3, u+l3, a, b, h);
			v += (k1+2.0*k2+2.0*k3+k4)/6.0;
			/********extract voltages**********/
			if (offset<1000 && offset%100==0){
				sw0->extrvol[(sw.sliceTime+it)*10+offset/100] = v;
			}
			/*************************************/
			if (v > 30.0) v = 30.0;
			if (v < -90.0) v = -90.0;

			u += (l1+2.0*l2+2.0*l3+l4)/6.0;

			sw0->voltage[offset] = v;
			sw0->recovery[offset]= u;

			/*** findFiring ***/
			if (v >= 30.0) {
				sw0->voltage[offset] = sw0->Izh_c[offset];
				sw0->recovery[offset] += sw0->Izh_d[offset];
				//if(addSpikeToTable_mpi(i)) assert(0);//????
				//addSpikeToTable
				/*int addr = atomicAdd(&(sw0->endST),1);
				sw0->firingTable[addr].nid = offset;
				sw0->firingTable[addr].time = sw.sliceTime+it;*/
			}
			

		} // end SizeN

#if 1
		/// reset ringBuffer
		int dIndex=sw0->offset+it;////????????
		//assert(dIndex<lenRB);
		
		//int i=0,j;
		j = threadIdx.x;
		for (i=0; i+j<sw.SizeN; i+=blockDim.x){
			int offset = sw.StartN - sw0->gStart + i + j;
			//int offset = i + j;
			//int addr = offset*sw0->lenRB + dIndex;
			int addr = offset + dIndex*sw0->gSize;
			if (sw.sim_with_conductances) {
				sw0->gAMPA[offset] += sw0->ringBuffer[addr];
				sw0->gNMDA_d[offset] += sw0->ringBuffer[addr];
			} else {
				///assert(0);
			}
			sw0->ringBuffer[addr]=0.;
		}
#endif
		}
		
	}
	if(threadIdx.x==0){
		sw.sliceTime += sw.Ndt;
		sw.simTime ++;
		///assert(sw.simTime == (sw.sliceTime>>sw.Nop));
		if(blockIdx.x==0){
		sw0->offset += sw.Ndt;//ringBuffer offset
		if(sw0->offset >= sw0->lenRB) sw0->offset -= sw0->lenRB;
		}
		swInfo[blockIdx.x] = sw;
		//SpikeDmaWrite_mpi(ptr);//mpi++++
	}
	  
}


__global__ void kernel_SpikeDeliver(swInfo_t *swInfo,swPara_t *sw0){
	__shared__ swInfo_t sw;
	sw = swInfo[blockIdx.x];
	int NSall = (sw0->preN-1)/100+1;//???
	int i,j;
	j = threadIdx.x;
	int iB = blockIdx.x;
	int preN = sw0->preN;
	int MaxN = sw0->MaxN;
	int Ndelay = sw0->Ndelay;
	int gStart = sw0->gStart;
	int ofs = sw0->offset;
	int lenRB = sw0->lenRB;
	float *wt = sw0->wt;
	unsigned int *postId = sw0->postId;
	float *ringBuffer = sw0->ringBuffer;
	//int StartN = swInfo->StartN;
	//int EndN = swInfo->StartN+swInfo->SizeN-1;
	//float synWt = sw0->synWt;

	for (int is=0; is<NSall; is++){
		//unsigned int time = sw0->firingTableAll[is].time;
		//unsigned int nid  = sw0->firingTableAll[is].nid;
		unsigned int time = 2;
		unsigned int nid  = 0;
		long offset = nid*Ndelay*NTh*MaxN;
		time = time - (time>>swInfo->Nop)<<swInfo->Nop;
		/*for (i=0; i+j<MaxN; i+=blockDim.x){
			int addr = offset + iB*MaxN + i + j;
			unsigned int pId = postId[addr];
			float Wt = wt[addr];
			if(pId!=0xffffffff){
				addr = (pId-gStart)*lenRB + ofs + time;
				ringBuffer[addr] += Wt;
			}
		}*/
	}

	//InputCurrent
	float Wt=0.0005/1.0;
	int   nspike=220;
       	float change = Wt*nspike;
	int dIndex = ofs+sw.Ndt/2;
	for (i=0; i+j<sw.SizeN; i+=blockDim.x){
		int offset = sw.StartN - sw0->gStart + i + j;
        	//if (dIndex>=lenRB)dIndex-=lenRB;
       		//ringBuffer[offset*lenRB+dIndex] += change;
       		ringBuffer[offset+dIndex*sw0->gSize] += change;
	}

}



