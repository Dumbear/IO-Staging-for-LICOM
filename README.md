# I/O Staging for LICOM 说明

## 目录结构

* `src/LICOM`
    * 修改后的LICOM源代码，目前只有`adios_io.F90`
* `src/test_skel`：skel测试
    * `licom2.xml`没有变化，`licom2_staging.xml`是给staging程序用的，但实际上只是把transport method从DataSpaces改成了MPI
    * `licom_skel.c`是改后的skel程序，用法为：
        * `mpiexec -n <n_proc> ./licom_skel <global_domain_x> <global_domain_y> <global_domain_z> <proc_x> <proc_y> <proc_z>`
        * 总的3维domain即为`global_domain_x * global_domain_y * global_domain_z`，例如低分辨率的LICOM就是`362 * 194 * 30`
        * `proc_x * proc_y * proc_z`必须等于`n_proc`，分别是3个方向划分的进程数，目前`proc_z`必须等于1
    * `staging.cpp`是staging程序，用法为：
        * `mpiexec -n <n_proc> ./staging <n_steps> <timeout> <global_domain_y> <global_domain_z> <proc_x> <proc_y> <proc_z>`
        * `n_steps`目前设为1即可，`timeout`设为一个比较大的数即可，其他参数跟skel意义相同，自由设置
    * `submit.sh`是提交脚本，其大致流程为：
        * 写DataSpaces配置文件`dataspaces.conf`
        * 运行DataSpaces server，并等待其生成配置文件`conf`
        * 运行LICOM skel程序
        * 运行staging程序
        * 等待所有程序结束

## DataSpaces说明

* 一个完整的I/O Staging程序应当是包含一个writer，一个reader，以及介于两者之间的DataSpaces server
* writer
    * 这里就是skel程序，可以替换为真实的LICOM程序（修改`submit.sh`即可）
    * writer可能包含多个进程，由于XML文件里定义的transport method是DataSpaces，所有进程在初始化ADIOS的时候会连接到DataSpaces server
    * writer有自己的data decomposition，每个进程都各自输出属于自己的数据，输出的流程是：
        * writer进程通过RPC请求告诉DataSpaces server要输出数据了
        * DataSpaces server收到RPC请求后响应，从writer进程处通过RDMA直接读取数据
        * writer进程等待server读取完毕后再继续执行
* reader
    * 这里就是staging程序
    * reader同样可能包含多个进程，由于代码里定义的read method是DataSpaces，所有进程在初始化read method的时候会连接到DataSpaces server
    * reader同样有自己的data decomposition，每个进程都各自读取属于自己的数据
* DataSpaces server
    * writer输出的数据在DataSpaces server处形成一个全局的视图，decomposition信息对外界是透明的
