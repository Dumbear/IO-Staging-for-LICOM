/******************************************************************************\
*                         Author:  Dumbear                                     *
*                         Email:   dumbear[#at]163.com                         *
*                         Website: http://dumbear.com                          *
\******************************************************************************/
#include "mpi.h"
#include "adios.h"
#include <cstdio>

using namespace std;

int nx, ny, start_x, start_y, count_x, count_y, data;

int proc_rank, proc_size;
int64_t fp_out;

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &proc_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &proc_size);

    nx = 2;
    ny = 2;
    start_x = proc_rank % 2;
    start_y = proc_rank / 2 % 2;
    count_x = 1;
    count_y = 1;
    data = proc_rank;

    adios_init("config.xml", MPI_COMM_WORLD);
    adios_open(&fp_out, "sample", "sample.bp", "w", MPI_COMM_WORLD);

    uint64_t group_size = sizeof(int) * 6 + sizeof(data), total_size;
    adios_group_size(fp_out, group_size, &total_size);

    adios_write(fp_out, "nx", &nx);
    adios_write(fp_out, "ny", &ny);
    adios_write(fp_out, "start_x", &start_x);
    adios_write(fp_out, "start_y", &start_y);
    adios_write(fp_out, "count_x", &count_x);
    adios_write(fp_out, "count_y", &count_y);
    adios_write(fp_out, "data", &data);

    adios_close(fp_out);
    adios_finalize(proc_rank);

    MPI_Finalize();
    return 0;
}
