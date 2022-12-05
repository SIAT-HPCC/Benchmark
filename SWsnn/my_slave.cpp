#include <stdio.h>
#include <math.h>
#include <string.h>
#include "slave.h"
#include <simd.h>
#include <dma.h>

#include "my_slave.h"

static __inline long rpcc()
{
   long a; asm volatile ("memb");
   asm volatile ("rcsr %0,4":"=r"(a));
   return a;
}

//==============Registor communication===================================
//simd data: intv8,uintv8,int256,uint256,floatv4,doublev4
//simd_load,simd_set_intv8
intv8 put_get_intv8(intv8 _v,int srcId)
{
    int dst=8;
    int col=COL(srcId),row=ROW(srcId);
    if(COL(_MYID)==col&&ROW(_MYID)==row){
	REG_PUTR(_v,dst);
	REG_PUTC(_v,dst);
    }
    else if(COL(_MYID)!=col&&ROW(_MYID)==row){
        REG_GETR(_v);
	REG_PUTC(_v,dst);
    }
    else if(COL(_MYID)!=col&&ROW(_MYID)!=row){
        REG_GETC(_v);
    }
    else if(COL(_MYID)==col&&ROW(_MYID)!=row){
        REG_GETC(_v);
    }
    return _v;
}
#if 0
    intv8 _v;
        ((int*)(&_v))[0]=_MYID;
        ((int*)(&_v))[1]=_MYID;
        ((int*)(&_v))[2]=_MYID;
        ((int*)(&_v))[3]=_MYID;
        ((int*)(&_v))[4]=_MYID;
        ((int*)(&_v))[5]=_MYID;
        ((int*)(&_v))[6]=_MYID;
        ((int*)(&_v))[7]=_MYID;
#endif
#if 0
void func()
{
        int i,j;

        my_id = athread_get_id(-1);
        
	get_reply = 0;
	athread_get(PE_MODE,&a[my_id][0],&a_[0],sizeof(int),&get_reply,0,0,0);
	athread_get(PE_MODE,&b[my_id][0],&b_[0],sizeof(int),&get_reply,0,0,0);
	while(get_reply!=2);
        for(i=0;i<I;i++){
                c_[i]=a_[i]+b_[i];
        }
	put_reply=0;
	athread_put(PE_MODE,&c_[0],&c[my_id][0],I*sizeof(int),&put_reply,0,0);
        while(put_reply!=1);	
}
#endif

#if 0
void func3()
{
        int i,j;
        my_id = athread_get_id(-1);
        //=============================
	get_reply = 0;
	dma_desc dma_get;
    if(_MYID==0){
	dma_set_size(&dma_get,1000*sizeof(int));
	dma_set_op(&dma_get,DMA_GET);
	dma_set_reply(&dma_get,&reply);
	dma_set_mode(&dma_get,PE_MODE);
	dma_set_bsize(&dma_get,0);
	dma_set_stepsize(&dma_get,0);
	dma(dma_get,&a[0][0],&a_slave[0]);
	dma_wait(&reply,1);
    }
	
	if(_MYID==0){
		put_reply=0;
		athread_put(PE_MODE,&time[0],&counter[0],13*sizeof(unsigned long),&put_reply,0,0);
        	while(put_reply!=1);	
	}
}
#endif
#if 0
void func5()
{
    unsigned long i=0,j=0;
    intv8 _v;
    //int dst=COL(_MYID)+1+8;
    int dst=8;
    //int dst=COL(_MYID);
    int col=7,row=0;
    if(COL(_MYID)==col&&ROW(_MYID)==row){
        ((int*)(&_v))[0]=_MYID;
        ((int*)(&_v))[1]=_MYID;
        ((int*)(&_v))[2]=_MYID;
        ((int*)(&_v))[3]=_MYID;
        ((int*)(&_v))[4]=_MYID;
        ((int*)(&_v))[5]=_MYID;
        ((int*)(&_v))[6]=_MYID;
        ((int*)(&_v))[7]=_MYID;

    	st = rpcc();
	REG_PUTR(_v,dst);
	REG_PUTC(_v,dst);
    	ed = rpcc();
	//REG_GETR(_v);
    }
    if(COL(_MYID)!=col&&ROW(_MYID)==row){
    	st = rpcc();
        REG_GETR(_v);
    	ed = rpcc();
	REG_PUTC(_v,dst);
    }
    if(COL(_MYID)!=col&&ROW(_MYID)!=row){
    	st = rpcc();
        REG_GETC(_v);
    	ed = rpcc();
    }
    if(COL(_MYID)==col&&ROW(_MYID)!=row){
    	st = rpcc();
        REG_GETC(_v);
    	ed = rpcc();
    }
    /*******************************/
    dst=COL(_MYID)+1;
    if(dst>=8) dst-=8;
    if(COL(_MYID)<=7){
        ((int*)(&_v))[0]=_MYID;
        ((int*)(&_v))[1]=_MYID;
        ((int*)(&_v))[2]=_MYID;
        ((int*)(&_v))[3]=_MYID;
        ((int*)(&_v))[4]=_MYID;
        ((int*)(&_v))[5]=_MYID;
        ((int*)(&_v))[6]=_MYID;
        ((int*)(&_v))[7]=_MYID;

    	REG_PUTR(_v,dst);
    }
    if(COL(_MYID)>=0){
        REG_GETR(_v);
    }
    /*******************************/



	i=(unsigned long)(((int*)(&_v))[2]);
//==============================================       
    result = (ed-st);
       
	put_reply=0;
	athread_put(PE_MODE,&result,&counter[_MYID],sizeof(unsigned long),&put_reply,0,0);
	athread_put(PE_MODE,&i,&counter[_MYID+64],sizeof(unsigned long),&put_reply,0,0);
	athread_put(PE_MODE,&st,&counter[_MYID+128],sizeof(unsigned long),&put_reply,0,0);
	athread_put(PE_MODE,&ed,&counter[_MYID+192],sizeof(unsigned long),&put_reply,0,0);
        while(put_reply!=4);	
//=================================================
}
#endif
