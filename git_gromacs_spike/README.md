# Gromacs 新冠病毒REMD模拟测试

## 软件介绍
Gromacs是分子动力模拟学的常用软件之一，是业内CPU并行的标杆软件，支持单节点,多节点，CUDA平台(包括移植到DCU平台)。
Gromacs不仅可以实现常规分子动力学模拟，同时支持多种增强采样算法，增强采样算法可以提高分子动力学模拟的采样效率，使得模拟分子规模和模拟时间得以延长，以副本交换算法为代表的多种增强采样算法已广泛应用在各种生物大分子的模拟计算中。
本次测试使用Gromacs完成新冠病毒S蛋白在不同温度下副本交换模拟，模拟包含的水分子数目为740356个，蛋白分子数44568个。

## 使用介绍

本次测试病毒体系已经生成完毕，只需要指定副本数目，模拟温度选择从290K-436.69K，交换概率为10%
1. 下载gromacs gpu版本并安装。

2. 创建模拟环境，指定副本数为512个, 命令执行完成后会生成副本交换的remd文件512个，为副本交换的路径。

    ``` bash do_build 512 ```

3. 分布式资源配置。编辑作业脚本文件submit.job  设置节点数、每节点任务数、每节点DCU数目、MPI进程数目（进程总数>=副本数目）

4. 在脚本文件submit.job中，指定gromacs运行参数，按照副本交换文件数目，先生成副本初始化文件初始化，然后指定副本交换参数，交换间隔-replex，模拟步长 -nsteps ，模拟路径-multidir 。

```
(for dir in remd{0..511}; do cd $dir; gmx_mpi  grompp -f remd.mdp -c init.gro -p fws_plus.top -o remd.tpr  -maxwarn 5; cd ..; done)
mpirun -np $num_procs gmx_mpi mdrun -multidir remd{0..511} -replex 500 -nsteps 5000 -deffnm remd
```
5. 提交作业，执行以下命令提交作业。

   ``` sbatch submit.job ```

5. 查看输出，每个副本的输出结果在remd文件夹下，包括轨迹文件remd.xtc, 运行日志remd.log等。
