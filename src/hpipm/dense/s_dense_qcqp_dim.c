#include <stdio.h>
#include <stdlib.h>

#include "hpipm/auxiliary/aux_mem.h"
#include "hpipm/auxiliary/aux_string.h"
#include "hpipm/dense/s_dense_qcqp_dim.h"
#include "hpipm/dense/s_dense_qp_dim.h"



#define DENSE_QCQP_DIM s_dense_qcqp_dim
#define DENSE_QP_DIM s_dense_qp_dim
#define DENSE_QP_DIM_CREATE s_dense_qp_dim_create
#define DENSE_QP_DIM_STRSIZE s_dense_qp_dim_strsize
#define DENSE_QP_DIM_MEMSIZE s_dense_qp_dim_memsize
#define DENSE_QP_DIM_SET s_dense_qp_dim_set
#define DENSE_QP_DIM_SET_NV s_dense_qp_dim_set_nv
#define DENSE_QP_DIM_SET_NE s_dense_qp_dim_set_ne
#define DENSE_QP_DIM_SET_NB s_dense_qp_dim_set_nb
#define DENSE_QP_DIM_SET_NG s_dense_qp_dim_set_ng
#define DENSE_QP_DIM_SET_NSB s_dense_qp_dim_set_nsb
#define DENSE_QP_DIM_SET_NSG s_dense_qp_dim_set_nsg
#define DENSE_QP_DIM_SET_NS s_dense_qp_dim_set_ns

#define DENSE_QCQP_DIM_STRSIZE s_dense_qcqp_dim_strsize
#define DENSE_QCQP_DIM_MEMSIZE s_dense_qcqp_dim_memsize
#define DENSE_QCQP_DIM_CREATE s_dense_qcqp_dim_create
#define DENSE_QCQP_DIM_SET s_dense_qcqp_dim_set
#define DENSE_QCQP_DIM_SET_NV s_dense_qcqp_dim_set_nv
#define DENSE_QCQP_DIM_SET_NE s_dense_qcqp_dim_set_ne
#define DENSE_QCQP_DIM_SET_NB s_dense_qcqp_dim_set_nb
#define DENSE_QCQP_DIM_SET_NG s_dense_qcqp_dim_set_ng
#define DENSE_QCQP_DIM_SET_NQ s_dense_qcqp_dim_set_nq
#define DENSE_QCQP_DIM_SET_NSB s_dense_qcqp_dim_set_nsb
#define DENSE_QCQP_DIM_SET_NSG s_dense_qcqp_dim_set_nsg
#define DENSE_QCQP_DIM_SET_NSQ s_dense_qcqp_dim_set_nsq
#define DENSE_QCQP_DIM_SET_NS s_dense_qcqp_dim_set_ns
#define DENSE_QCQP_DIM_GET_NV s_dense_qcqp_dim_get_nv
#define DENSE_QCQP_DIM_GET_NE s_dense_qcqp_dim_get_ne
#define DENSE_QCQP_DIM_GET_NB s_dense_qcqp_dim_get_nb
#define DENSE_QCQP_DIM_GET_NG s_dense_qcqp_dim_get_ng
#define DENSE_QCQP_DIM_GET_NQ s_dense_qcqp_dim_get_nq


#include "x_dense_qcqp_dim.c"



