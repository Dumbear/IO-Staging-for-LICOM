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
ADIOS_SELECTION *sel;
int data[10];

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &proc_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &proc_size);

    dim_start[0] = 0; dim_start[1] = 0; dim_start[2] = 0;
    dim_count[0] = 2; dim_count[1] = 2; dim_count[2] = 2;
    sel = adios_selection_boundingbox(3, dim_start, dim_count);
    memset(data, -1, sizeof(data));

    adios_read_init_method(ADIOS_READ_METHOD_FLEXPATH, MPI_COMM_WORLD, "");
    fp_in = adios_read_open("sample.bp", ADIOS_READ_METHOD_FLEXPATH, MPI_COMM_WORLD, ADIOS_LOCKMODE_NONE, 0.0);
    ADIOS_VARINFO *nx_info = adios_inq_var(fp_in, "/global/nx");
    ADIOS_VARINFO *ny_info = adios_inq_var(fp_in, "/global/ny");
    ADIOS_VARINFO *nz_info = adios_inq_var(fp_in, "/global/nz");
    printf("nx = %d\n", *(int *)nx_info->value);
    printf("ny = %d\n", *(int *)ny_info->value);
    printf("nz = %d\n", *(int *)nz_info->value);
    adios_schedule_read(fp_in, sel, "/var/data", 0, 1, data);
    adios_perform_reads(fp_in, 1);

    printf("[");
    for (int i = 0; i < 7; ++i) {
        printf("%d, ", data[i]);
    }
    printf("%d]\n", data[7]);

    adios_read_close(fp_in);
    adios_read_finalize_method(ADIOS_READ_METHOD_FLEXPATH);

    MPI_Finalize();
    return 0;
}
