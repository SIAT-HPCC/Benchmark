CC = hipcc
MPI = -I/opt/hpc/software/mpi/openmpi/4.0.4/gcc-7.3.1/include/openmpi -I/opt/hpc/software/mpi/openmpi/4.0.4/gcc-7.3.1/include/openmpi/opal/mca/event/libevent2022/libevent -I/opt/hpc/software/mpi/openmpi/4.0.4/gcc-7.3.1/include/openmpi/opal/mca/event/libevent2022/libevent/include -I/opt/hpc/software/mpi/openmpi/4.0.4/gcc-7.3.1/include -pthread -Wl,-rpath -Wl,/opt/hpc/software/mpi/openmpi/4.0.4/gcc-7.3.1/lib -Wl,--enable-new-dtags -L/opt/hpc/software/mpi/openmpi/4.0.4/gcc-7.3.1/lib -lmpi

all:
	rm -f swtest
	$(CC) $(MPI) -w -c swtest.cpp -O2 
	$(CC) $(MPI) -w -c setupNetwork.cpp -O2 
	$(CC) $(MPI) -w -c runNetwork.cpp -O2 
	$(CC) $(MPI) -w -c cuda_runNetwork.cpp -O2 
	srun -n 1 --gres=dcu:1 $(CC) $(MPI) -w -o swtest swtest.o setupNetwork.o runNetwork.o cuda_runNetwork.o
	rm -f *.o
