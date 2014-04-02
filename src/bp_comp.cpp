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

ADIOS_FILE *fp1, *fp2;
string filename1, filename2;
char *buf1, *buf2;

bool compare_variable_metadata(ADIOS_VARINFO *v1, ADIOS_VARINFO *v2) {
    if (v1->type != v2->type) {
        cout << "Different types: {" << adios_type_to_string(v1->type) << "}|{" << adios_type_to_string(v2->type) << "}." << endl;
        return false;
    }
    if (v1->nsteps != v2->nsteps) {
        cout << "Different #steps: {" << v1->nsteps << "}|{" << v2->nsteps << "}." << endl;
        return false;
    }
    if (v1->ndim != v2->ndim) {
        cout << "Different #dimensions: {" << v1->ndim << "}|{" << v2->ndim << "}." << endl;
        return false;
    }
    for (int i = 0; i < v1->ndim; ++i) {
        if (v1->dims[i] != v2->dims[i]) {
            cout << "Different dimensions #" << i << ": {" << v1->dims[i] << "}|{" << v2->dims[i] << "}." << endl;
            return false;
        }
    }
    return true;
}

uint64_t read_variable(ADIOS_FILE *fp, ADIOS_VARINFO *v, char **buf) {
    uint64_t size = adios_type_size(v->type, v->value) * v->nsteps;
    for (int i = 0; i < v->ndim; ++i) {
        size *= v->dims[i];
    }
    *buf = new char[size + 1024];
    adios_schedule_read_byid(fp, NULL, v->varid, 0, v->nsteps, *buf);
    adios_perform_reads(fp, 1);
    return size;
}

bool compare_data(char *buf1, char *buf2, uint64_t size) {
    uint64_t cnt = 0;
    for (uint64_t i = 0; i < size; ++i) {
        if (buf1[i] != buf2[i]) {
            ++cnt;
        }
    }
    if (cnt != 0) {
        cout << (double)cnt / size * 100.0 << "% (" << cnt << "/" << size << ") data is different." << endl;
    } else {
        cout << "They are identical!" << endl;
    }
    delete[] buf1;
    delete[] buf2;
    return cnt == 0;
}

void process() {
    adios_read_init_method(ADIOS_READ_METHOD_BP, MPI_COMM_WORLD, "");
    fp1 = adios_read_open_file(filename1.c_str(), ADIOS_READ_METHOD_BP, MPI_COMM_WORLD);
    if (fp1 == NULL) {
        return;
    }
    fp2 = adios_read_open_file(filename2.c_str(), ADIOS_READ_METHOD_BP, MPI_COMM_WORLD);
    if (fp2 == NULL) {
        adios_read_close(fp1);
        return;
    }

    if (fp1->nvars == fp2->nvars) {
        for (int i = 0; i < fp1->nvars; ++i) {
            cout << "Comparing variable #" << i << "..." << endl << '\t';
            if (strcmp(fp1->var_namelist[i], fp2->var_namelist[i]) != 0) {
                cout << "Different names, {" << fp1->var_namelist[i] << "}|{" << fp2->var_namelist[i] << "}." << endl;
                continue;
            }
            cout << "Name: " << fp1->var_namelist[i] << endl << '\t';
            ADIOS_VARINFO *v1 = adios_inq_var_byid(fp1, i);
            ADIOS_VARINFO *v2 = adios_inq_var_byid(fp2, i);
            if (compare_variable_metadata(v1, v2)) {
                uint64_t size1 = read_variable(fp1, v1, &buf1);
                uint64_t size2 = read_variable(fp2, v2, &buf2);
                compare_data(buf1, buf2, min(size1, size2));
            }
            adios_free_varinfo(v1);
            adios_free_varinfo(v2);
        }
    } else {
        cout << "Different #variables, {" << fp1->nvars << "}|{" << fp2->nvars << "}." << endl;
    }

    adios_read_close(fp1);
    adios_read_close(fp2);
    adios_read_finalize_method(ADIOS_READ_METHOD_BP);
}

bool parse_arguments(int argc, char **argv) {
    if (argc != 3) {
        return false;
    }
    filename1 = argv[1];
    filename2 = argv[2];
    return true;
}

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);

    if (parse_arguments(argc, argv)) {
        cout << "Comparing file <" << filename1 << "> to file <" << filename2 << ">..." << endl;
        process();
    }

    MPI_Finalize();
    return 0;
}
