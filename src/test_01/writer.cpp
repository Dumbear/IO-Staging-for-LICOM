/******************************************************************************\
*                         Author:  Dumbear                                     *
*                         Email:   dumbear[#at]163.com                         *
*                         Website: http://dumbear.com                          *
\******************************************************************************/
#include "mpi.h"
#include "adios.h"
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
int proc_rank, proc_size;

int dim_global, dim_start, dim_count;
uint64_t group_size, total_size;
float *data;

int64_t fp_out;

void prepare_writing() {
    int from = (int64_t)dim_global * proc_rank / proc_size;
    int to = (int64_t)dim_global * (proc_rank + 1) / proc_size;
    dim_start = from;
    dim_count = to - from;

    data = new float[dim_count];
    for (int i = 0; i < dim_count; ++i) {
        data[i] = (float)rand() / RAND_MAX;
    }
    group_size = 0;
    group_size += sizeof(dim_global) + sizeof(dim_start) + sizeof(dim_count);
    group_size += sizeof(data[0]) * dim_count;
}

void process() {
    prepare_writing();

    adios_init("writer_adios.xml", io_comm);
    adios_open(&fp_out, "test", "writer_test.bp", "w", io_comm);
    adios_group_size(fp_out, group_size, &total_size);
    adios_write(fp_out, "dim/global", &dim_global);
    adios_write(fp_out, "dim/start", &dim_start);
    adios_write(fp_out, "dim/count", &dim_count);
    adios_write(fp_out, "data", data);
    adios_close(fp_out);
    adios_finalize(proc_rank);

    delete[] data;
}

bool parse_arguments(int argc, char **argv) {
    if (argc != 2) {
        return false;
    }

    int total_size;
    if (sscanf(argv[1], "%d", &total_size) != 1) {
        return false;
    }
    if (total_size <= 0 || total_size < proc_size) {
        return false;
    }
    dim_global = total_size * (1 << 20);

    return true;
}

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);
    io_comm = MPI_COMM_WORLD;
    MPI_Comm_rank(io_comm, &proc_rank);
    MPI_Comm_size(io_comm, &proc_size);

    if (parse_arguments(argc, argv)) {
        process();
    }

    MPI_Finalize();
    return 0;
}
