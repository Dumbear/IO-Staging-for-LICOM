/******************************************************************************\
*                         Author:  Dumbear                                     *
*                         Email:   dumbear[#at]163.com                         *
*                         Website: http://dumbear.com                          *
\******************************************************************************/
#include "mpi.h"
#include "adios.h"
#include "adios_read.h"
#include <cstdio>
#include <cstring>

using namespace std;

int proc_rank, proc_size;
ADIOS_FILE *fp_in;

uint64_t dim_start[10], dim_count[10];
ADIOS_SELECTION *sel2d, *sel3d;
int data[10];

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &proc_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &proc_size);

    dim_start[0] = 0; dim_start[1] = 0; dim_start[2] = 0;
    dim_count[0] = 2; dim_count[1] = 2; dim_count[2] = 2;
    sel2d = adios_selection_boundingbox(2, dim_start, dim_count);
    sel3d = adios_selection_boundingbox(3, dim_start, dim_count);

    adios_read_init_method(ADIOS_READ_METHOD_FLEXPATH, MPI_COMM_WORLD, "");
    fp_in = adios_read_open("sample.bp", ADIOS_READ_METHOD_FLEXPATH, MPI_COMM_WORLD, ADIOS_LOCKMODE_NONE, 0.0);

    ADIOS_VARINFO *nx_info = adios_inq_var(fp_in, "nx");
    ADIOS_VARINFO *ny_info = adios_inq_var(fp_in, "ny");
    ADIOS_VARINFO *nz_info = adios_inq_var(fp_in, "nz");
    printf("nx = %d\n", *(int *)nx_info->value);
    printf("ny = %d\n", *(int *)ny_info->value);
    printf("nz = %d\n", *(int *)nz_info->value);

    memset(data, -1, sizeof(data));
    adios_schedule_read(fp_in, sel2d, "data2d", 0, 1, data);
    adios_perform_reads(fp_in, 1);
    for (int i = 0; i < 4; ++i) {
        printf(i == 3 ? "%d\n" : "%d, ", data[i]);
    }

    memset(data, -1, sizeof(data));
    adios_schedule_read(fp_in, sel3d, "data3d", 0, 1, data);
    adios_perform_reads(fp_in, 1);
    for (int i = 0; i < 8; ++i) {
        printf(i == 7 ? "%d\n" : "%d, ", data[i]);
    }

    adios_read_close(fp_in);
    adios_read_finalize_method(ADIOS_READ_METHOD_FLEXPATH);

    MPI_Finalize();
    return 0;
}
