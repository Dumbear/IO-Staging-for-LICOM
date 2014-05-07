#include "adios.h"
#include "mpi.h"
#include "skel/skel_xml_output.h"
#include <stdlib.h>
#include <stdio.h>

int skel_mpi_size, skel_mpi_rank;
int proc_x, proc_y, proc_z, rank_x, rank_y, rank_z;
int ni_global;
int nj_global;
int nk_global;
int ni_offset;
int nj_offset;
int nk_offset;
int ni_local;
int nj_local;
int nk_local;

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
}

int parse_arguments(int argc, char **argv) {
    if (argc != 7) {
        return 0;
    }

    if (sscanf(argv[1], "%d", &ni_global) != 1) {
        return 0;
    }
    if (sscanf(argv[2], "%d", &nj_global) != 1) {
        return 0;
    }
    if (sscanf(argv[3], "%d", &nk_global) != 1) {
        return 0;
    }

    if (ni_global <= 0 || nj_global <= 0 || nk_global <= 0) {
        return 0;
    }

    if (sscanf(argv[4], "%d", &proc_x) != 1) {
        return 0;
    }
    if (sscanf(argv[5], "%d", &proc_y) != 1) {
        return 0;
    }
    if (sscanf(argv[6], "%d", &proc_z) != 1) {
        return 0;
    }

    if (proc_x * proc_y * proc_z != skel_mpi_size || proc_z != 1) {
        return 0;
    }

    rank_x = skel_mpi_rank % proc_x;
    rank_y = skel_mpi_rank / proc_x % proc_y;
    rank_z = skel_mpi_rank / proc_x / proc_y % proc_z;

    return 1;
}

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);

    double skel_init_timer = 0;
    double skel_open_timer = 0;
    double skel_access_timer = 0;
    double skel_close_timer = 0;
    double skel_total_timer = 0;

    // Time the init
    MPI_Barrier(MPI_COMM_WORLD);
    skel_init_timer -= MPI_Wtime();
    adios_init("licom2.xml", MPI_COMM_WORLD);
    skel_init_timer += MPI_Wtime();

    int skel_i;
    uint64_t adios_groupsize;
    MPI_Comm_rank(MPI_COMM_WORLD, &skel_mpi_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &skel_mpi_size);

    if (!parse_arguments(argc, argv)) {
        adios_finalize(skel_mpi_rank);
        MPI_Finalize();
        return 0;
    }

    int64_t adios_handle;
    uint64_t skel_total_size;

    // Scalar declarations
    int number_month;
    int number_day;
    decompose();
    number_month = 1;
    number_day = 1;

    // Array declarations

    double *h0;

    double *u;
    double *v;
    double *at1;
    double *at2;
    double *ws;

    double *su;
    double *sv;
    double *swv;
    double *lwv;
    double *sshf;
    double *lthf;
    double *fresh;

    double *t_cpl;
    double *s_cpl;
    double *u_cpl;
    double *v_cpl;
    double *dhdx;
    double *dhdy;
    double *q;

    h0 = (double*) malloc (ni_local*nj_local * sizeof (double) );
    for (skel_i = 0; skel_i < ni_local*nj_local; skel_i++) 
        h0[skel_i] = (double) skel_mpi_rank;
    u = (double*) malloc (ni_local*nj_local*nk_local * sizeof (double) );
    for (skel_i = 0; skel_i < ni_local*nj_local*nk_local; skel_i++) 
        u[skel_i] = (double) skel_mpi_rank;
    v = (double*) malloc (ni_local*nj_local*nk_local * sizeof (double) );
    for (skel_i = 0; skel_i < ni_local*nj_local*nk_local; skel_i++) 
        v[skel_i] = (double) skel_mpi_rank;
    at1 = (double*) malloc (ni_local*nj_local*nk_local * sizeof (double) );
    for (skel_i = 0; skel_i < ni_local*nj_local*nk_local; skel_i++) 
        at1[skel_i] = (double) skel_mpi_rank;
    at2 = (double*) malloc (ni_local*nj_local*nk_local * sizeof (double) );
    for (skel_i = 0; skel_i < ni_local*nj_local*nk_local; skel_i++) 
        at2[skel_i] = (double) skel_mpi_rank;
    ws = (double*) malloc (ni_local*nj_local*nk_local * sizeof (double) );
    for (skel_i = 0; skel_i < ni_local*nj_local*nk_local; skel_i++) 
        ws[skel_i] = (double) skel_mpi_rank;
    su = (double*) malloc (ni_local*nj_local * sizeof (double) );
    for (skel_i = 0; skel_i < ni_local*nj_local; skel_i++) 
        su[skel_i] = (double) skel_mpi_rank;
    sv = (double*) malloc (ni_local*nj_local * sizeof (double) );
    for (skel_i = 0; skel_i < ni_local*nj_local; skel_i++) 
        sv[skel_i] = (double) skel_mpi_rank;
    swv = (double*) malloc (ni_local*nj_local * sizeof (double) );
    for (skel_i = 0; skel_i < ni_local*nj_local; skel_i++) 
        swv[skel_i] = (double) skel_mpi_rank;
    lwv = (double*) malloc (ni_local*nj_local * sizeof (double) );
    for (skel_i = 0; skel_i < ni_local*nj_local; skel_i++) 
        lwv[skel_i] = (double) skel_mpi_rank;
    sshf = (double*) malloc (ni_local*nj_local * sizeof (double) );
    for (skel_i = 0; skel_i < ni_local*nj_local; skel_i++) 
        sshf[skel_i] = (double) skel_mpi_rank;
    lthf = (double*) malloc (ni_local*nj_local * sizeof (double) );
    for (skel_i = 0; skel_i < ni_local*nj_local; skel_i++) 
        lthf[skel_i] = (double) skel_mpi_rank;
    fresh = (double*) malloc (ni_local*nj_local * sizeof (double) );
    for (skel_i = 0; skel_i < ni_local*nj_local; skel_i++) 
        fresh[skel_i] = (double) skel_mpi_rank;
    t_cpl = (double*) malloc (ni_local*nj_local * sizeof (double) );
    for (skel_i = 0; skel_i < ni_local*nj_local; skel_i++) 
        t_cpl[skel_i] = (double) skel_mpi_rank;
    s_cpl = (double*) malloc (ni_local*nj_local * sizeof (double) );
    for (skel_i = 0; skel_i < ni_local*nj_local; skel_i++) 
        s_cpl[skel_i] = (double) skel_mpi_rank;
    u_cpl = (double*) malloc (ni_local*nj_local * sizeof (double) );
    for (skel_i = 0; skel_i < ni_local*nj_local; skel_i++) 
        u_cpl[skel_i] = (double) skel_mpi_rank;
    v_cpl = (double*) malloc (ni_local*nj_local * sizeof (double) );
    for (skel_i = 0; skel_i < ni_local*nj_local; skel_i++) 
        v_cpl[skel_i] = (double) skel_mpi_rank;
    dhdx = (double*) malloc (ni_local*nj_local * sizeof (double) );
    for (skel_i = 0; skel_i < ni_local*nj_local; skel_i++) 
        dhdx[skel_i] = (double) skel_mpi_rank;
    dhdy = (double*) malloc (ni_local*nj_local * sizeof (double) );
    for (skel_i = 0; skel_i < ni_local*nj_local; skel_i++) 
        dhdy[skel_i] = (double) skel_mpi_rank;
    q = (double*) malloc (ni_local*nj_local * sizeof (double) );
    for (skel_i = 0; skel_i < ni_local*nj_local; skel_i++) 
        q[skel_i] = (double) skel_mpi_rank;

    skel_total_timer -= MPI_Wtime();
    for (skel_i = 0; skel_i < 1; skel_i++) {
        // Time the opens
        MPI_Barrier(MPI_COMM_WORLD);
        skel_open_timer -= MPI_Wtime();
        MPI_Comm comm = MPI_COMM_WORLD;
        adios_open(&adios_handle, "ssaveins", "out_ssaveins.bp", (skel_i == 0 ? "w" : "a"), comm);
        skel_open_timer += MPI_Wtime();

        // Time the writes
        skel_access_timer -= MPI_Wtime();

        // Set the adios group size
        adios_groupsize =
                         4 +
                         4 +
                         4 +
                         4 +
                         4 +
                         4 +
                         4 +
                         4 +
                         4 +
                         4 +
                         4 +
                         (uint64_t)8 * (ni_local) * (nj_local) +
                         (uint64_t)8 * (ni_local) * (nj_local) * (nk_local) +
                         (uint64_t)8 * (ni_local) * (nj_local) * (nk_local) +
                         (uint64_t)8 * (ni_local) * (nj_local) * (nk_local) +
                         (uint64_t)8 * (ni_local) * (nj_local) * (nk_local) +
                         (uint64_t)8 * (ni_local) * (nj_local) * (nk_local) +
                         (uint64_t)8 * (ni_local) * (nj_local) +
                         (uint64_t)8 * (ni_local) * (nj_local) +
                         (uint64_t)8 * (ni_local) * (nj_local) +
                         (uint64_t)8 * (ni_local) * (nj_local) +
                         (uint64_t)8 * (ni_local) * (nj_local) +
                         (uint64_t)8 * (ni_local) * (nj_local) +
                         (uint64_t)8 * (ni_local) * (nj_local) +
                         (uint64_t)8 * (ni_local) * (nj_local) +
                         (uint64_t)8 * (ni_local) * (nj_local) +
                         (uint64_t)8 * (ni_local) * (nj_local) +
                         (uint64_t)8 * (ni_local) * (nj_local) +
                         (uint64_t)8 * (ni_local) * (nj_local) +
                         (uint64_t)8 * (ni_local) * (nj_local) +
                         (uint64_t)8 * (ni_local) * (nj_local);
        adios_group_size(adios_handle, adios_groupsize, &skel_total_size);

        // Write each variable
        adios_write(adios_handle, "/dimensions/ni_global", &ni_global);
        adios_write(adios_handle, "/dimensions/nj_global", &nj_global);
        adios_write(adios_handle, "/dimensions/nk_global", &nk_global);
        adios_write(adios_handle, "/aux/ni_offset", &ni_offset);
        adios_write(adios_handle, "/aux/nj_offset", &nj_offset);
        adios_write(adios_handle, "/aux/nk_offset", &nk_offset);
        adios_write(adios_handle, "/aux/ni_local", &ni_local);
        adios_write(adios_handle, "/aux/nj_local", &nj_local);
        adios_write(adios_handle, "/aux/nk_local", &nk_local);
        adios_write(adios_handle, "/var/number_month", &number_month);
        adios_write(adios_handle, "/var/number_day", &number_day);
        adios_write(adios_handle, "/var/h0", h0);
        adios_write(adios_handle, "/var/u", u);
        adios_write(adios_handle, "/var/v", v);
        adios_write(adios_handle, "/var/at1", at1);
        adios_write(adios_handle, "/var/at2", at2);
        adios_write(adios_handle, "/var/ws", ws);
        adios_write(adios_handle, "/var/su", su);
        adios_write(adios_handle, "/var/sv", sv);
        adios_write(adios_handle, "/var/swv", swv);
        adios_write(adios_handle, "/var/lwv", lwv);
        adios_write(adios_handle, "/var/sshf", sshf);
        adios_write(adios_handle, "/var/lthf", lthf);
        adios_write(adios_handle, "/var/fresh", fresh);
        adios_write(adios_handle, "/var/t_cpl", t_cpl);
        adios_write(adios_handle, "/var/s_cpl", s_cpl);
        adios_write(adios_handle, "/var/u_cpl", u_cpl);
        adios_write(adios_handle, "/var/v_cpl", v_cpl);
        adios_write(adios_handle, "/var/dhdx", dhdx);
        adios_write(adios_handle, "/var/dhdy", dhdy);
        adios_write(adios_handle, "/var/q", q);

        // Stop timing the writes
        skel_access_timer += MPI_Wtime();

        // Time the closes
        skel_close_timer -= MPI_Wtime();
        adios_close(adios_handle);
        skel_close_timer += MPI_Wtime();
    }
    skel_total_timer += MPI_Wtime();

    // Output results
    skel_write_coarse_xml_data(skel_open_timer, skel_access_timer, skel_close_timer, skel_total_timer);
    double skel_total_init, skel_total_open, skel_total_access, skel_total_close, skel_total_total;
    MPI_Reduce(&skel_init_timer, &skel_total_init, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&skel_open_timer, &skel_total_open, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&skel_access_timer, &skel_total_access, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&skel_close_timer, &skel_total_close, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&skel_total_timer, &skel_total_total, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    if (skel_mpi_rank == 0) {
        fprintf(stdout, "\n");
        fprintf(stdout, "\n*************************");
        fprintf(stdout, "\n   Groupsize: %lli", adios_groupsize);
        fprintf(stdout, "\n  Init Time: %f", skel_total_init);
        fprintf(stdout, "\n  Open Time: %f", skel_total_open);
        fprintf(stdout, "\nAccess Time: %f", skel_total_access);
        fprintf(stdout, "\n Close Time: %f", skel_total_close);
        fprintf(stdout, "\n Total Time: %f", skel_total_total);
        fprintf(stdout, "\n*************************");
        fprintf(stdout, "\n");
    }

    // Free the arrays

    free(h0);
    free(u);
    free(v);
    free(at1);
    free(at2);
    free(ws);
    free(su);
    free(sv);
    free(swv);
    free(lwv);
    free(sshf);
    free(lthf);
    free(fresh);
    free(t_cpl);
    free(s_cpl);
    free(u_cpl);
    free(v_cpl);
    free(dhdx);
    free(dhdy);
    free(q);

    // Clean up
    adios_finalize(skel_mpi_rank);
    MPI_Finalize();
    return 0;
}
