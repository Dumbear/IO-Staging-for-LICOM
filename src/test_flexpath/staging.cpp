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
float timeout_sec;
ADIOS_FILE *fp_in;
int64_t fp_out;
char *buf;

int ni_global, nj_global, nk_global, ni_offset, nj_offset, nk_offset, ni_local, nj_local, nk_local;
uint64_t dim_start[4][8], dim_count[4][8];

template<typename T>
void all_output(T data) {
    cout << "Rank[" << proc_rank << "]: " << data << endl;
}
template<typename T1, typename T2>
void all_output(T1 data1, T2 data2) {
    cout << "Rank[" << proc_rank << "]: " << data1 << ' ' << data2 << endl;
}
template<typename T1, typename T2, typename T3>
void all_output(T1 data1, T2 data2, T3 data3) {
    cout << "Rank[" << proc_rank << "]: " << data1 << ' ' << data2 << ' ' << data3 << endl;
}

template<typename T>
void single_output(T data) {
    if (proc_rank == 0) {
        cout << "Master: " << data << endl;
    }
}
template<typename T1, typename T2>
void single_output(T1 data1, T2 data2) {
    if (proc_rank == 0) {
        cout << "Master: " << data1 << ' ' << data2 << endl;
    }
}
template<typename T1, typename T2, typename T3>
void single_output(T1 data1, T2 data2, T3 data3) {
    if (proc_rank == 0) {
        cout << "Master: " << data1 << ' ' << data2 << ' ' << data3 << endl;
    }
}

void decompose() {
    if (ni_global == -1 || nj_global == -1 || nk_global == -1) {
        return;
    }
    if (proc_x > ni_global || proc_y > nj_global || proc_z > nk_global) {
        return;
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

    dim_start[2][0] = ni_offset; dim_start[2][1] = nj_offset;
    dim_start[3][0] = ni_offset; dim_start[3][1] = nj_offset; dim_start[3][2] = nk_offset;
    dim_count[2][0] = ni_local; dim_count[2][1] = nj_local;
    dim_count[3][0] = ni_local; dim_count[3][1] = nj_local; dim_count[3][2] = nk_local;
    reverse(dim_start[2], dim_start[2] + 2);
    reverse(dim_count[2], dim_count[2] + 2);
    reverse(dim_start[3], dim_start[3] + 3);
    reverse(dim_count[3], dim_count[3] + 3);
}

void read_single(int id, void *data, ADIOS_SELECTION *sel = NULL) {
    all_output("Reading variable", id, string() + "- " + fp_in->var_namelist[id] + "...");
    adios_schedule_read(fp_in, sel, fp_in->var_namelist[id], 0, 1, data);
    all_output("Variable", id, string() + "- " + fp_in->var_namelist[id] + " scheduled.");
    adios_perform_reads(fp_in, 1);
    all_output("Variable", id, string() + "- " + fp_in->var_namelist[id] + " read.");
}

void write_single(int id, void *data) {
    all_output("Writing variable", id, string() + "- " + fp_in->var_namelist[id] + "...");
    adios_write(fp_out, fp_in->var_namelist[id], data);
    all_output("Variable", id, string() + "- " + fp_in->var_namelist[id] + " written.");
}

void process_step() {
    all_output("Processing first step...");
    all_output("Opening write file <" + filename_out + ">...");
    adios_open(&fp_out, "ssaveins", filename_out.c_str(), "w", io_comm);
    all_output("Write file opened.");
    all_output("There are", fp_in->nvars, "variables to write.");

    decompose();
    uint64_t group_size = 4 + 4 + 4 * 9, total_size;
    group_size += 8 * 5 * dim_count[3][0] * dim_count[3][1] * dim_count[3][2];
    group_size += 8 * 15 * dim_count[2][0] * dim_count[2][1];
    adios_group_size(fp_out, group_size, &total_size);

    buf = new char[8 * dim_count[3][0] * dim_count[3][1] * dim_count[3][2] + 1024];
    for (int i = 0; i < fp_in->nvars; ++i) {
        string name(fp_in->var_namelist[i]);
        ADIOS_VARINFO *v = adios_inq_var_byid(fp_in, i);
        if (v->ndim > 0) {
            ADIOS_SELECTION *sel = adios_selection_boundingbox(v->ndim, dim_start[v->ndim], dim_count[v->ndim]);
            read_single(i, buf, sel);
            write_single(i, buf);
        } else if (name == "/dimensions/ni_global") {
            read_single(i, &ni_global);
            write_single(i, &ni_global);
        } else if (name == "/dimensions/nj_global") {
            read_single(i, &nj_global);
            write_single(i, &nj_global);
        } else if (name == "/dimensions/nk_global") {
            read_single(i, &nk_global);
            write_single(i, &nk_global);
        } else if (name == "/aux/ni_offset") {
            write_single(i, &ni_offset);
        } else if (name == "/aux/nj_offset") {
            write_single(i, &nj_offset);
        } else if (name == "/aux/nk_offset") {
            write_single(i, &nk_offset);
        } else if (name == "/aux/ni_local") {
            write_single(i, &ni_local);
        } else if (name == "/aux/nj_local") {
            write_single(i, &nj_local);
        } else if (name == "/aux/nk_local") {
            write_single(i, &nk_local);
        } else {
            read_single(i, buf);
            write_single(i, buf);
        }
    }
    delete[] buf;

    adios_release_step(fp_in);
    adios_close(fp_out);
    all_output("First step processed.");
}

bool process_single() {
    all_output("Opening read file <" + filename_in + ">...");
    fp_in = adios_read_open(filename_in.c_str(), method_read, io_comm, ADIOS_LOCKMODE_ALL, timeout_sec);
    if (fp_in == NULL) {
        if (adios_errno == err_file_not_found) {
            all_output("[ERROR] Read file not found!");
        } else if (adios_errno == err_end_of_stream) {
            all_output("[ERROR] Read file has ended!");
        } else {
            all_output("[ERROR] Read file has some unknown problems!");
        }
        all_output("[ERROR] Read file opend failed!");
        return false;
    } else {
        all_output("Read file opened.");
        if (fp_in->current_step != 0) {
            all_output("[ERROR] First step missing, stop current processing!");
            return false;
        }
        process_step();
    }
    adios_read_close(fp_in);
    return true;
}

void advance_day(int &year, int &month, int &day) {
    int ndays[] = {-1, 31, -1, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
    ndays[2] = (year % 400 == 0 || year % 4 == 0 && year % 100 != 0 ? 29 : 28);
    if (++day > ndays[month]) {
        day = 1;
        if (++month > 12) {
            month = 1;
            ++year;
        }
    }
}

void process() {
    all_output("Initializing read method...");
    adios_read_init_method(method_read, io_comm, "max_chunk_size=100; app_id =32767; \nverbose= 3;poll_interval  =  100;");
    all_output("Read method initialized.");

    all_output("Initializing ADIOS...");
    adios_init("licom2_staging.xml", io_comm);
    all_output("ADIOS initialized.");

    int year = 1, month = 1, day = 3;
    for (int i = 0; i < n_steps; ++i) {
        char buf[128];
        // sprintf(buf, "fort.22.%04d-%02d-%02d", year, month, day);
        sprintf(buf, "out_ssaveins");
        filename_in = string() + buf + ".bp";
        filename_out = string() + buf + "-staged.bp";
        process_single();
        advance_day(year, month, day);
    }

    adios_read_finalize_method(method_read);
    adios_finalize(proc_rank);
}

bool parse_arguments(int argc, char **argv) {
    if (argc != 9) {
        return false;
    }

    if (sscanf(argv[1], "%d", &n_steps) != 1) {
        single_output("[ERROR] <n_steps> is not valid!");
        return false;
    }
    single_output("<n_steps> =", n_steps);

    if (sscanf(argv[2], "%f", &timeout_sec) != 1) {
        single_output("[ERROR] <timeout_sec> is not valid!");
        return false;
    }
    single_output("<timeout_sec> =", timeout_sec);

    if (sscanf(argv[3], "%d", &ni_global) != 1) {
        return false;
    }
    if (sscanf(argv[4], "%d", &nj_global) != 1) {
        return false;
    }
    if (sscanf(argv[5], "%d", &nk_global) != 1) {
        return false;
    }
    if (ni_global <= 0 || nj_global <= 0 || nk_global <= 0) {
        single_output("[ERROR] <ni_global> | <nj_global> | <nk_global> is not valid!");
        return false;
    }

    if (sscanf(argv[6], "%d", &proc_x) != 1) {
        single_output("[ERROR] <proc_x> is not valid!");
        return false;
    }
    single_output("<proc_x> =", proc_x);

    if (sscanf(argv[7], "%d", &proc_y) != 1) {
        single_output("[ERROR] <proc_y> is not valid!");
        return false;
    }
    single_output("<proc_y> =", proc_y);

    if (sscanf(argv[8], "%d", &proc_z) != 1) {
        single_output("[ERROR] <proc_z> is not valid!");
        return false;
    }
    single_output("<proc_z> =", proc_z);

    if (proc_x * proc_y * proc_z != proc_size) {
        single_output("[ERROR] <proc_x> * <proc_y> * <proc_z> not equals <proc_size>!");
        return false;
    }

    rank_x = proc_rank % proc_x;
    rank_y = proc_rank / proc_x % proc_y;
    rank_z = proc_rank / proc_x / proc_y % proc_z;

    method_read = ADIOS_READ_METHOD_FLEXPATH;
    return true;
}

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);
    io_comm = MPI_COMM_WORLD;
    MPI_Comm_rank(io_comm, &proc_rank);
    MPI_Comm_size(io_comm, &proc_size);

    if (parse_arguments(argc, argv)) {
        all_output("Staging application started. #proc =", proc_size);
        process();
        all_output("Staging application ended. #proc =", proc_size);
    }

    MPI_Finalize();
    return 0;
}
