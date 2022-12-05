#ifndef _MY_GPU_H_
#define _MY_GPU_H_
#include <stdint.h>
#include <stdlib.h> 
#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <assert.h>
#include <sys/time.h>    // for gettimeofday() chenged!!
#include "mpi.h"
#include "swstruct.h"
#include "mysnn.h"
#include <pthread.h>
#define NTh 64

//#include <unistd.h>        // for sleep()
//#define NBLOCK 64
#define NBLOCK NTh
#define NTHREAD 128 

void initSW(snnInfo_t *sInfo, swInfo_t *swInfo);
void freeSW(snnInfo_t *sInfo);
void StateUpdate(snnInfo_t *sInfo);
void SpikeDeliver(snnInfo_t *sInfo);

#endif
