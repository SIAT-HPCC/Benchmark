#include <stdio.h>
#include <math.h>
#include <string.h>
#include "slave.h"
#include <simd.h>
#include <dma.h>

#include "my_slave.h"

//==============Registor communication===================================
//simd data: intv8,uintv8,int256,uint256,floatv4,doublev4
//simd_load,simd_set_intv8
int dmaWait(int* reply,int value)
{
	while(reply[0]!=value);	
    	return 0;
}
