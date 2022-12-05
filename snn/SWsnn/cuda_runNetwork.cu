#include <stdio.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include "swstruct.h"
#include "mygpu.h"

//#include"./common/book.h"
//#include <unistd.h>        // for sleep()
//#define NBLOCK 64
//#define NTHREAD 128 

/*swInfo_t *cu_swInfo;
neurInfo_t *cu_nInfo;
synInfo_t *cu_sInfo;*/

/**************for delay*************************/
/*spikeTime_t *cu_firingTable;
float *cu_ringBuffer;*/
/*************************************************/

// 检查显卡错误
void checkError() {
    cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess) {
        std::cout << cudaGetErrorString(err) << std::endl;
    }
}
void initSW(snnInfo_t *sInfo, swInfo_t *swInfo){//????
  cudaMallocHost(&sInfo->nInfoHost, sInfo->gSize*sizeof(neurInfo_t));  // 内存空间
  cudaMalloc(&sInfo->nInfoDevice, sInfo->gSize*sizeof(neurInfo_t));  // 显存空间
  cudaMemcpy(sInfo->nInfoDevice, sInfo->nInfoHost,
    sInfo->gSize*sizeof(neurInfo_t), cudaMemcpyHostToDevice);  // 拷贝数据
  
	//以下两个变量(神经元、突触)可有可无
	//cudaMalloc 神经元数据 
	//显存变量：sInfo->nInfoDevice 对应主存变量：sInfo->nInfoHost (类型neurInfo_t*)
	//
	//空间大小sInfo->gSize*sizeof(neurInfo_t)
	
	//cudaMemcpy 以上主存到显存数据拷贝
	
  cudaMallocHost(&sInfo->sInfoHost,
    sInfo->preN*sInfo->Ndelay*sInfo->MaxN*NTh*sizeof(synInfo_t));  // 内存空间
  cudaMalloc(&sInfo->sInfoDevice,
    sInfo->preN*sInfo->Ndelay*sInfo->MaxN*NTh*sizeof(synInfo_t));  // 显存空间
  cudaMemcpy(sInfo->nInfoDevice, sInfo->nInfoHost,
    sInfo->preN*sInfo->Ndelay*sInfo->MaxN*NTh*sizeof(synInfo_t),
    cudaMemcpyHostToDevice);  // 拷贝数据
	//cudaMalloc 突触数据 
	//显存变量：sInfo->sInfoDevice 对应主存变量：sInfo->sInfoHost (类型synInfo_t*)
	//空间大小(long)sInfo->preN*sInfo->Ndelay*sInfo->MaxN*NTh*sizeof(synInfo_t)
	                
	//cudaMemcpy 以上主存到显存数据拷贝

  cudaMallocHost(&sInfo->firingTableHost, sInfo->gSize*sizeof(spikeTime_t));  // 内存空间
  cudaMalloc(&sInfo->firingTableDevice, sInfo->gSize*sizeof(spikeTime_t));  // 显存空间
  cudaMemcpy(sInfo->firingTableDevice, sInfo->firingTableHost,
    sInfo->gSize*sizeof(spikeTime_t), cudaMemcpyHostToDevice);  // 拷贝数据
  // 没有清零
	//以下两个变量(局部、全局脉冲事件信息)测试HtoD通信，必须建立
	//cudaMalloc 本地脉冲事件表 
	//显存变量：sInfo->firingTableDevice 对应主存变量sInfo->firingTableHost (类型spikeTime_t*)
	//空间大小sInfo->gSize*sizeof(spikeTime_t)
	//cudaMemset 空间清零(可有可无)
	
  cudaMallocHost(&sInfo->firingTableAll,
    sInfo->NN*sInfo->Ndelay*sizeof(spikeTime_t));  // 内存空间
  cudaMalloc(&sInfo->firingTableAllDevice,
    sInfo->NN*sInfo->Ndelay*sizeof(spikeTime_t));  // 显存空间
  cudaMemcpy(sInfo->firingTableAllDevice, sInfo->firingTableAll,
    sInfo->NN*sInfo->Ndelay*sizeof(spikeTime_t), cudaMemcpyHostToDevice);  // 拷贝数据
  // 没有清零
	//cudaMalloc 全局脉冲事件表 
	//显存变量: sInfo->firingTableAllDevice 对应主存变量sInfo->firingTableAll (类型spikeTime_t*)
	//空间大小sInfo->NN*sInfo->Ndelay*sizeof(spikeTime_t)
	//cudaMemset 空间清零(可有可无)



  checkError();  // 检测错误
	return;	
}

void freeSW(){
  // 没有释放
	return;
}
// kernel_empty 空核函数 只是启动
__global__ void kernel_empty() {
  // 没有内容
}
void StateUpdate(snnInfo_t *sInfo){
	//kernel 函数启动一次
  kernel_empty<<<128, 128>>>();
  cudaDeviceSynchronize();  // 等待核函数运行结束
	//DtoH脉冲事件显存到主存传递
  cudaMemcpy(sInfo->firingTableAll, sInfo->firingTableAllDevice,
    ((sInfo->gSize-1)/100+1)*sizeof(spikeTime_t), cudaMemcpyDeviceToHost);  // 拷贝数据
	//显存变量: sInfo->firingTableAllDevice 对应主存变量sInfo->firingTableAll (类型spikeTime_t*)
	//传递大小: ((sInfo->gSize-1)/100+1)*sizeof(spikeTime_t)

	return;
}
void SpikeDeliver(snnInfo_t *sInfo){

  cudaMemcpy(sInfo->firingTableAllDevice, sInfo->firingTableAll,
    ((sInfo->NN-1)/100+1)*sizeof(spikeTime_t), cudaMemcpyHostToDevice);  // 拷贝数据

	//HtoD脉冲事件主存到显存传递
	//显存变量: sInfo->firingTableAllDevice 对应主存变量sInfo->firingTableAll (类型spikeTime_t*)
	//传递大小: ((sInfo->NN-1)/100+1)*sizeof(spikeTime_t)
	
	//kernel 函数启动一次
  kernel_empty<<<128, 128>>>();
  cudaDeviceSynchronize();  // 等待核函数运行结束
	return;
}
//inline float dvdtIzh(float v, float u, float tmpI, float h);
//inline float dudtIzh(float v, float u, float a, float b, float h);
//inline float dvdtIzh(float v, float u, float tmpI, float h) {return (((0.04*v+5.0)*v+140.0-u+tmpI)*h);}
//inline float dudtIzh(float v, float u, float a, float b, float h) {return (a*(b*v-u)*h);}
