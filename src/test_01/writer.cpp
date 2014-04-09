/******************************************************************************\
*                         Author:  Dumbear                                     *
*                         Email:   dumbear[#at]163.com                         *
*                         Website: http://dumbear.com                          *
\******************************************************************************/
#include "mpi.h"
#include "adios.h"
extern "C" {
    #include "skel/skel_xml_output.h"
}
#include <algorithm>
#include <bitset>
#include <cctype>
#include <climits>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <functional>
#include <iomanip>
#include <iostream>
#include <list>
#include <map>
#include <numeric>
#include <queue>
#include <set>
#include <sstream>
#include <stack>
#include <string>
#include <typeinfo>
#include <utility>
#include <vector>

using namespace std;

#define output(x) cout << #x << ": " << (x) << endl;

typedef long long LL;
typedef vector<int> VI;
typedef vector<long long> VL;
typedef vector<double> VD;
typedef vector<string> VS;

MPI_Comm io_comm;
int proc_rank, proc_size, proc_x, proc_y, proc_z, rank_x, rank_y, rank_z;

int dim_global_x, dim_global_y, dim_global_z;
int dim_start_x, dim_start_y, dim_start_z;
int dim_count_x, dim_count_y, dim_count_z;
uint64_t group_size, total_size;
float *data;

int64_t fp_out;

double timer[4];

void prepare_writing() {
    int from, to;

    from = (int64_t)dim_global_x * rank_x / proc_x;
    to = (int64_t)dim_global_x * (rank_x + 1) / proc_x;
    dim_start_x = from;
    dim_count_x = to - from;

    from = (int64_t)dim_global_y * rank_y / proc_y;
    to = (int64_t)dim_global_y * (rank_y + 1) / proc_y;
    dim_start_y = from;
    dim_count_y = to - from;

    from = (int64_t)dim_global_z * rank_z / proc_z;
    to = (int64_t)dim_global_z * (rank_z + 1) / proc_z;
    dim_start_z = from;
    dim_count_z = to - from;

    uint64_t dim_count = (uint64_t)dim_count_x * dim_count_y * dim_count_z;
    data = new float[dim_count];
    for (int i = 0; i < dim_count; ++i) {
        data[i] = (float)rand() / RAND_MAX;
    }
    group_size = 0;
    group_size += (sizeof(dim_global_x) + sizeof(dim_start_x) + sizeof(dim_count_x)) * 3;
    group_size += sizeof(*data) * dim_count;
}

void process() {
    prepare_writing();

    adios_init("writer_adios.xml", io_comm);
    adios_open(&fp_out, "test", "writer_test.bp", "w", io_comm);
    adios_group_size(fp_out, group_size, &total_size);

    timer[3] -= MPI_Wtime();

    adios_write(fp_out, "dim/global_x", &dim_global_x);
    adios_write(fp_out, "dim/global_y", &dim_global_y);
    adios_write(fp_out, "dim/global_z", &dim_global_z);
    adios_write(fp_out, "dim/start_x", &dim_start_x);
    adios_write(fp_out, "dim/start_y", &dim_start_y);
    adios_write(fp_out, "dim/start_z", &dim_start_z);
    adios_write(fp_out, "dim/count_x", &dim_count_x);
    adios_write(fp_out, "dim/count_y", &dim_count_y);
    adios_write(fp_out, "dim/count_z", &dim_count_z);

    timer[1] -= MPI_Wtime();
    adios_write(fp_out, "data", data);
    timer[1] += MPI_Wtime();

    timer[3] += MPI_Wtime();

    timer[2] -= MPI_Wtime();
    adios_close(fp_out);
    timer[2] += MPI_Wtime();
    adios_finalize(proc_rank);

    delete[] data;
}

bool parse_arguments(int argc, char **argv) {
    if (argc != 7) {
        return false;
    }

    if (sscanf(argv[1], "%d", &dim_global_x) != 1 || dim_global_x <= 0) {
        return false;
    }
    if (sscanf(argv[2], "%d", &dim_global_y) != 1 || dim_global_y <= 0) {
        return false;
    }
    if (sscanf(argv[3], "%d", &dim_global_z) != 1 || dim_global_z <= 0) {
        return false;
    }

    if (sscanf(argv[4], "%d", &proc_x) != 1 || dim_global_x < proc_x) {
        return false;
    }
    if (sscanf(argv[5], "%d", &proc_y) != 1 || dim_global_y < proc_y) {
        return false;
    }
    if (sscanf(argv[6], "%d", &proc_z) != 1 || dim_global_z < proc_z) {
        return false;
    }
    if (proc_x * proc_y * proc_z != proc_size) {
        return false;
    }

    rank_x = proc_rank % proc_x;
    rank_y = proc_rank / proc_x % proc_y;
    rank_z = proc_rank / proc_x / proc_y % proc_z;

    return true;
}

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);
    io_comm = MPI_COMM_WORLD;
    MPI_Comm_rank(io_comm, &proc_rank);
    MPI_Comm_size(io_comm, &proc_size);

    if (parse_arguments(argc, argv)) {
        process();
        skel_write_coarse_xml_data(timer[0], timer[1], timer[2], timer[3]);
    }

    MPI_Finalize();
    return 0;
}
