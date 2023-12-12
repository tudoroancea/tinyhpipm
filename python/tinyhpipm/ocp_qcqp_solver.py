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


class hpipm_ocp_qcqp_solver:
    def __init__(self, qcqp_dims, arg):
        # load hpipm shared library
        __hpipm = CDLL("libhpipm.so")
        self.__hpipm = __hpipm

        # set up ipm workspace struct
        sizeof_d_ocp_qcqp_ipm_workspace = __hpipm.d_ocp_qcqp_ipm_ws_strsize()
        ipm_ws_struct = cast(
            create_string_buffer(sizeof_d_ocp_qcqp_ipm_workspace), c_void_p
        )
        self.ipm_ws_struct = ipm_ws_struct

        # allocate memory for ipm workspace
        ipm_size = __hpipm.d_ocp_qcqp_ipm_ws_memsize(
            qcqp_dims.dim_struct, arg.arg_struct
        )
        ipm_ws_mem = cast(create_string_buffer(ipm_size), c_void_p)
        self.ipm_ws_mem = ipm_ws_mem

        # create C ws
        __hpipm.d_ocp_qcqp_ipm_ws_create(
            qcqp_dims.dim_struct, arg.arg_struct, ipm_ws_struct, ipm_ws_mem
        )

        self.arg = arg
        self.dim_struct = qcqp_dims.dim_struct

    def solve(self, qp, qp_sol):
        self.__hpipm.d_ocp_qcqp_ipm_solve(
            qp.qp_struct, qp_sol.qp_sol_struct, self.arg.arg_struct, self.ipm_ws_struct
        )

    def get(self, field):
        if field == "stat":
            # get iters
            iters = np.zeros((1, 1), dtype=int)
            tmp = cast(iters.ctypes.data, POINTER(c_int))
            self.__hpipm.d_ocp_qcqp_ipm_get_iter(self.ipm_ws_struct, tmp)
            # get stat_m
            stat_m = np.zeros((1, 1), dtype=int)
            tmp = cast(stat_m.ctypes.data, POINTER(c_int))
            self.__hpipm.d_ocp_qcqp_ipm_get_stat_m(self.ipm_ws_struct, tmp)
            # get stat pointer
            res = np.zeros((iters[0][0] + 1, stat_m[0][0]))
            ptr = c_void_p()
            self.__hpipm.d_ocp_qcqp_ipm_get_stat(self.ipm_ws_struct, byref(ptr))
            tmp = cast(ptr, POINTER(c_double))
            for ii in range(iters[0][0] + 1):
                for jj in range(stat_m[0][0]):
                    res[ii][jj] = tmp[jj + ii * stat_m[0][0]]
            return res
        elif field in {"status", "iter"}:
            res = np.zeros((1, 1), dtype=int)
            tmp = cast(res.ctypes.data, POINTER(c_int))
        elif field in {"max_res_stat", "max_res_eq", "max_res_ineq", "max_res_comp"}:
            res = np.zeros((1, 1))
            tmp = cast(res.ctypes.data, POINTER(c_double))
        else:
            raise NameError("hpipm_ocp_qcqp_solver.get: wrong field")
        field_name_b = field.encode("utf-8")
        self.__hpipm.d_ocp_qcqp_ipm_get(c_char_p(field_name_b), self.ipm_ws_struct, tmp)
        return res[0][0]
