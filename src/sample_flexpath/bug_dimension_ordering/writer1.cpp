/******************************************************************************\
*                         Author:  Dumbear                                     *
*                         Email:   dumbear[#at]163.com                         *
*                         Website: http://dumbear.com                          *
\******************************************************************************/
#include "mpi.h"
#include "adios.h"
#include <cstdio>

using namespace std;

int nx, ny, nz, start_x, start_y, start_z, count_x, count_y, count_z, data[2];

int proc_rank, proc_size;
int64_t fp_out;

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &proc_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &proc_size);

    nx = 2;
    ny = 2;
    nz = 2;
    start_x = proc_rank / 2 % 2;
    start_y = proc_rank % 2;
    start_z = 0;
    count_x = 1;
    count_y = 1;
    count_z = 2;
    data[0] = data[1] = proc_rank;

    adios_init("config.xml", MPI_COMM_WORLD);
    adios_open(&fp_out, "sample", "sample.bp", "w", MPI_COMM_WORLD);

    uint64_t group_size = sizeof(int) * 9 + sizeof(data), total_size;
    adios_group_size(fp_out, group_size, &total_size);

    adios_write(fp_out, "nx", &nx);
    adios_write(fp_out, "ny", &ny);
    adios_write(fp_out, "nz", &nz);
    adios_write(fp_out, "start_x", &start_x);
    adios_write(fp_out, "start_y", &start_y);
    adios_write(fp_out, "start_z", &start_z);
    adios_write(fp_out, "count_x", &count_x);
    adios_write(fp_out, "count_y", &count_y);
    adios_write(fp_out, "count_z", &count_z);
    adios_write(fp_out, "data", data);

    adios_close(fp_out);
    adios_finalize(proc_rank);

    MPI_Finalize();
    return 0;
}
