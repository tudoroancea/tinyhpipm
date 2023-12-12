from ctypes import (
    CDLL,
    POINTER,
    byref,
    c_char_p,
    c_double,
    c_int,
    c_void_p,
    cast,
    create_string_buffer,
)

import numpy as np


class hpipm_ocp_qp_solver:
    def __init__(self, qp_dims, arg):
        # load hpipm shared library
        __hpipm = CDLL("libhpipm.so")
        self.__hpipm = __hpipm

        # set up ipm workspace struct
        sizeof_d_ocp_qp_ipm_workspace = __hpipm.d_ocp_qp_ipm_ws_strsize()
        ipm_ws_struct = cast(
            create_string_buffer(sizeof_d_ocp_qp_ipm_workspace), c_void_p
        )
        self.ipm_ws_struct = ipm_ws_struct

        # allocate memory for ipm workspace
        ipm_size = __hpipm.d_ocp_qp_ipm_ws_memsize(qp_dims.dim_struct, arg.arg_struct)
        ipm_ws_mem = cast(create_string_buffer(ipm_size), c_void_p)
        self.ipm_ws_mem = ipm_ws_mem

        # create C ws
        __hpipm.d_ocp_qp_ipm_ws_create(
            qp_dims.dim_struct, arg.arg_struct, ipm_ws_struct, ipm_ws_mem
        )

        self.arg = arg
        self.dim_struct = qp_dims.dim_struct

        # getter functions for feedback matrices
        self.__getters = {
            "ric_Lr": {
                "n_row": __hpipm.d_ocp_qp_dim_get_nu,
                "n_col": __hpipm.d_ocp_qp_dim_get_nu,
                "var": __hpipm.d_ocp_qp_ipm_get_ric_Lr,
            },
            "ric_Ls": {
                "n_row": __hpipm.d_ocp_qp_dim_get_nx,
                "n_col": __hpipm.d_ocp_qp_dim_get_nu,
                "var": __hpipm.d_ocp_qp_ipm_get_ric_Ls,
            },
            "ric_P": {
                "n_row": __hpipm.d_ocp_qp_dim_get_nx,
                "n_col": __hpipm.d_ocp_qp_dim_get_nx,
                "var": __hpipm.d_ocp_qp_ipm_get_ric_P,
            },
            "ric_lr": {
                "n_row": __hpipm.d_ocp_qp_dim_get_nu,
                "var": __hpipm.d_ocp_qp_ipm_get_ric_lr,
            },
            "ric_p": {
                "n_row": __hpipm.d_ocp_qp_dim_get_nx,
                "var": __hpipm.d_ocp_qp_ipm_get_ric_p,
            },
            "ric_K": {
                "n_row": __hpipm.d_ocp_qp_dim_get_nu,
                "n_col": __hpipm.d_ocp_qp_dim_get_nx,
                "var": __hpipm.d_ocp_qp_ipm_get_ric_K,
            },
            "ric_k": {
                "n_row": __hpipm.d_ocp_qp_dim_get_nu,
                "var": __hpipm.d_ocp_qp_ipm_get_ric_k,
            },
        }

    def solve(self, qp, qp_sol):
        self.__hpipm.d_ocp_qp_ipm_solve(
            qp.qp_struct, qp_sol.qp_sol_struct, self.arg.arg_struct, self.ipm_ws_struct
        )

    def get_feedback(self, qp, field, idx_start, idx_end=None):
        if field not in self.__getters:
            raise NameError("hpipm_ocp_qp_solver.get: wrong field")
        else:
            if idx_end is None:
                idx_end = idx_start
            return self.__get_feedback(self.__getters[field], qp, idx_start, idx_end)

    def __get_feedback(self, getter, qp, idx_start, idx_end):
        n_row = np.zeros(1, dtype=int)
        n_row_ptr = cast(n_row.ctypes.data, POINTER(c_int))
        n_col = np.zeros(1, dtype=int)
        n_col_ptr = cast(n_col.ctypes.data, POINTER(c_int))
        var = []
        for i in range(idx_start, idx_end + 1):
            # get dimension of feedback matrix at stage
            getter["n_row"](self.dim_struct, i, n_row_ptr)
            if "n_col" in getter:
                getter["n_col"](self.dim_struct, i, n_col_ptr)
            else:
                n_col[0] = 1

            var.append(np.zeros((n_row[0], n_col[0]), dtype=float))
            tmp_ptr = cast(var[-1].ctypes.data, POINTER(c_double))
            getter["var"](
                qp.qp_struct, self.arg.arg_struct, self.ipm_ws_struct, i, tmp_ptr
            )

        return var if len(var) > 1 else var[0]

    def get(self, field):
        if field == "stat":
            # get iters
            iters = np.zeros((1, 1), dtype=int)
            tmp = cast(iters.ctypes.data, POINTER(c_int))
            self.__hpipm.d_ocp_qp_ipm_get_iter(self.ipm_ws_struct, tmp)
            # get stat_m
            stat_m = np.zeros((1, 1), dtype=int)
            tmp = cast(stat_m.ctypes.data, POINTER(c_int))
            self.__hpipm.d_ocp_qp_ipm_get_stat_m(self.ipm_ws_struct, tmp)
            # get stat pointer
            res = np.zeros((iters[0][0] + 1, stat_m[0][0]))
            ptr = c_void_p()
            self.__hpipm.d_ocp_qp_ipm_get_stat(self.ipm_ws_struct, byref(ptr))
            tmp = cast(ptr, POINTER(c_double))
            for ii in range(iters[0][0] + 1):
                for jj in range(stat_m[0][0]):
                    res[ii][jj] = tmp[jj + ii * stat_m[0][0]]
            return res
        elif (field == "status") | (field == "iter"):
            res = np.zeros((1, 1), dtype=int)
            tmp = cast(res.ctypes.data, POINTER(c_int))
        elif (
            (field == "max_res_stat")
            | (field == "max_res_eq")
            | (field == "max_res_ineq")
            | (field == "max_res_comp")
        ):
            res = np.zeros((1, 1))
            tmp = cast(res.ctypes.data, POINTER(c_double))
        else:
            raise NameError("hpipm_ocp_qp_solver.get: wrong field")
        field_name_b = field.encode("utf-8")
        self.__hpipm.d_ocp_qp_ipm_get(c_char_p(field_name_b), self.ipm_ws_struct, tmp)
        return res[0][0]
