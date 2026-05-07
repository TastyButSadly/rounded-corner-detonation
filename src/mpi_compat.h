#ifndef MPI_COMPAT_H
#define MPI_COMPAT_H

#if __has_include(<mpi.h>)
#include <mpi.h>
#else

#include <algorithm>
#include <cstdlib>

using MPI_Comm = int;
using MPI_Datatype = int;
using MPI_Op = int;

struct MPI_Status {
    int dummy = 0;
};

constexpr MPI_Comm MPI_COMM_WORLD = 0;
constexpr MPI_Comm MPI_COMM_NULL = -1;
constexpr int MPI_PROC_NULL = -1;
constexpr MPI_Datatype MPI_DOUBLE = 0;
constexpr MPI_Op MPI_MIN = 0;
inline MPI_Status* MPI_STATUS_IGNORE = nullptr;

inline int MPI_Init(int*, char***) {
    return 0;
}

inline int MPI_Finalize() {
    return 0;
}

inline int MPI_Comm_rank(MPI_Comm, int* rank) {
    *rank = 0;
    return 0;
}

inline int MPI_Comm_size(MPI_Comm, int* size) {
    *size = 1;
    return 0;
}

inline int MPI_Abort(MPI_Comm, int errorcode) {
    std::exit(errorcode);
}

inline int MPI_Dims_create(int nnodes, int ndims, int dims[]) {
    for (int axis = 0; axis < ndims; ++axis) {
        if (dims[axis] == 0) {
            dims[axis] = 1;
        }
    }
    if (ndims >= 1) {
        dims[0] = nnodes;
    }
    if (ndims >= 2) {
        dims[1] = 1;
    }
    return 0;
}

inline int MPI_Cart_create(MPI_Comm, int, const int[], const int[], int, MPI_Comm* newcomm) {
    *newcomm = MPI_COMM_WORLD;
    return 0;
}

inline int MPI_Cart_coords(MPI_Comm, int, int maxdims, int coords[]) {
    for (int axis = 0; axis < maxdims; ++axis) {
        coords[axis] = 0;
    }
    return 0;
}

inline int MPI_Cart_shift(MPI_Comm, int, int, int* source, int* dest) {
    *source = MPI_PROC_NULL;
    *dest = MPI_PROC_NULL;
    return 0;
}

inline int MPI_Barrier(MPI_Comm) {
    return 0;
}

inline int MPI_Allreduce(const double* sendbuf, double* recvbuf, int count,
                         MPI_Datatype, MPI_Op, MPI_Comm) {
    std::copy(sendbuf, sendbuf + count, recvbuf);
    return 0;
}

inline int MPI_Sendrecv(const double*, int, MPI_Datatype, int, int,
                        double*, int, MPI_Datatype, int, int,
                        MPI_Comm, MPI_Status*) {
    return 0;
}

inline int MPI_Send(const double*, int, MPI_Datatype, int, int, MPI_Comm) {
    return 0;
}

inline int MPI_Recv(double*, int, MPI_Datatype, int, int, MPI_Comm, MPI_Status*) {
    return 0;
}

#endif

#endif
