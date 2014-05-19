/******************************************************************************\
*                         Author:  Dumbear                                     *
*                         Email:   dumbear[#at]163.com                         *
*                         Website: http://dumbear.com                          *
\******************************************************************************/
#include "mpi.h"
#include "adios.h"
#include "adios_read.h"
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
int proc_rank, proc_size, n_steps, proc_x, proc_y, proc_z, rank_x, rank_y, rank_z;

ADIOS_READ_METHOD method_read;
string filename_in, filename_out;
ADIOS_FILE *fp_in;
int64_t fp_out;

int ni_global, nj_global, nk_global, ni_offset, nj_offset, nk_offset, ni_local, nj_local, nk_local;
ADIOS_SELECTION *selection_2d, *selection_3d;
map<string, void *> var_map;
char *buf;

bool decompose() {
    ADIOS_VARINFO *ni_info = adios_inq_var(fp_in, "/dimensions/ni_global");
    ADIOS_VARINFO *nj_info = adios_inq_var(fp_in, "/dimensions/nj_global");
    ADIOS_VARINFO *nk_info = adios_inq_var(fp_in, "/dimensions/nk_global");

    ni_global = *(int *)ni_info->value;
    nj_global = *(int *)nj_info->value;
    nk_global = *(int *)nk_info->value;

    if (proc_x > ni_global || proc_y > nj_global || proc_z > nk_global) {
        return false;
    }

    uint64_t from, to;

    from = ni_global * rank_x / proc_x;
    to = ni_global * (rank_x + 1) / proc_x;
    ni_offset = from;
    ni_local = to - from;

    from = nj_global * rank_y / proc_y;
    to = nj_global * (rank_y + 1) / proc_y;
    nj_offset = from;
    nj_local = to - from;

    from = nk_global * rank_z / proc_z;
    to = nk_global * (rank_z + 1) / proc_z;
    nk_offset = from;
    nk_local = to - from;

    uint64_t dim[2][3] = {{nk_offset, nj_offset, ni_offset}, {nk_local, nj_local, ni_local}};
    selection_2d = adios_selection_boundingbox(2, dim[0] + 1, dim[1] + 1);
    selection_3d = adios_selection_boundingbox(3, dim[0], dim[1]);

    var_map.clear();
    var_map["/aux/ni_offset"] = &ni_offset;
    var_map["/aux/nj_offset"] = &nj_offset;
    var_map["/aux/nk_offset"] = &nk_offset;
    var_map["/aux/ni_local"] = &ni_local;
    var_map["/aux/nj_local"] = &nj_local;
    var_map["/aux/nk_local"] = &nk_local;

    return true;
}

void read_variable(int id, void *data, ADIOS_SELECTION *sel = NULL) {
    adios_schedule_read(fp_in, sel, fp_in->var_namelist[id], 0, 1, data);
    adios_perform_reads(fp_in, 1);
}

void write_variable(int id, void *data) {
    adios_write(fp_out, fp_in->var_namelist[id], data);
    printf("Rank[%d]: Variable %d - <%s> written.\n", proc_rank, id, fp_in->var_namelist[id]);
}

void process_step() {
    if (!decompose()) {
        return;
    }

    printf("Rank[%d]: Data decomposed successfully - global(%d, %d, %d), offset(%d, %d, %d), count(%d, %d, %d).\n", proc_rank,
        ni_global, nj_global, nk_global, ni_offset, nj_offset, nk_offset, ni_local, nj_local, nk_local);

    adios_open(&fp_out, "ssaveins", filename_out.c_str(), "w", io_comm);
    printf("Rank[%d]: File %s opened for write.\n", proc_rank, filename_out.c_str());

    uint64_t group_size = 0, total_size;
    group_size += 4 + 4 + 4 * 9;
    group_size += (uint64_t)8 * 5 * ni_local * nj_local * nk_local;
    group_size += (uint64_t)8 * 15 * ni_local * nj_local;
    adios_group_size(fp_out, group_size, &total_size);

    buf = new char[sizeof(double) * ni_local * nj_local * nk_local + 1024];
    for (int i = 0; i < fp_in->nvars; ++i) {
        string name(fp_in->var_namelist[i]);
        ADIOS_VARINFO *v = adios_inq_var_byid(fp_in, i);
        switch (v->ndim) {
            case 0:
                if (var_map.find(name) != var_map.end()) {
                    write_variable(i, var_map[name]);
                } else {
                    ADIOS_VARINFO *info = adios_inq_var_byid(fp_in, i);
                    write_variable(i, info->value);
                }
                break;
            case 2:
                read_variable(i, buf, selection_2d);
                write_variable(i, buf);
                break;
            case 3:
                read_variable(i, buf, selection_3d);
                write_variable(i, buf);
                break;
            default:
                printf("Rank[%d]: Unexpected variable <%s>, dimension = %d!", proc_rank, name.c_str(), v->ndim);
                break;
        }
    }
    delete[] buf;

    adios_release_step(fp_in);
    adios_close(fp_out);
    printf("Rank[%d]: Step processed successfully.\n\n", proc_rank);
}

void process() {
    fp_in = adios_read_open(filename_in.c_str(), method_read, io_comm, ADIOS_LOCKMODE_NONE, 0.0);
    if (fp_in == NULL) {
    } else {
        if (fp_in->current_step == 0) {
            printf("Rank[%d]: File %s opened for read.\n", proc_rank, filename_in.c_str());
            process_step();
        }
        adios_read_close(fp_in);
    }
}

bool parse_arguments(int argc, char **argv) {
    if (argc != 5) {
        return false;
    }

    if (sscanf(argv[1], "%d", &n_steps) != 1 || n_steps <= 0) {
        return false;
    }
    if (sscanf(argv[2], "%d", &proc_x) != 1 || proc_x <= 0) {
        return false;
    }
    if (sscanf(argv[3], "%d", &proc_y) != 1 || proc_y <= 0) {
        return false;
    }
    if (sscanf(argv[4], "%d", &proc_z) != 1 || proc_z <= 0) {
        return false;
    }

    if (proc_x * proc_y * proc_z != proc_size) {
        return false;
    }

    if (proc_rank == 0) {
        printf("<n_steps> = %d\n", n_steps);
        printf("Data decomposed on %d = %d * %d * %d processes\n\n", proc_size, proc_x, proc_y, proc_z);
    }

    rank_x = proc_rank % proc_x;
    rank_y = proc_rank / proc_x % proc_y;
    rank_z = proc_rank / proc_x / proc_y % proc_z;

    method_read = ADIOS_READ_METHOD_FLEXPATH;
    return true;
}

void get_filenames() {
    filename_in = "out_ssaveins.bp";
    filename_out = "out_ssaveins-staged.bp";
}

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);
    io_comm = MPI_COMM_WORLD;
    MPI_Comm_rank(io_comm, &proc_rank);
    MPI_Comm_size(io_comm, &proc_size);

    if (parse_arguments(argc, argv)) {
        adios_read_init_method(method_read, io_comm, "");
        printf("Rank[%d]: FLEXPATH method initialized.\n", proc_rank);
        adios_init("licom2_staging.xml", io_comm);
        printf("Rank[%d]: ADIOS initialized.\n", proc_rank);
        for (int i = 0; i < n_steps; ++i) {
            get_filenames();
            process();
        }
        adios_read_finalize_method(method_read);
        adios_finalize(proc_rank);
    } else {
        if (proc_rank == 0) {
            printf("Usage: %s <n_steps> <proc_x> <proc_y> <proc_z>\n", argv[0]);
        }
    }

    MPI_Finalize();
    printf("Rank[%d]: All done.", proc_rank);
    return 0;
}
