from ctypes import (
    CDLL,
    POINTER,
    c_double,
    c_int,
    c_void_p,
    cast,
    create_string_buffer,
)
import numpy as np


class hpipm_ocp_qp_sol:
    def __init__(self, dim):
        # save dim internally
        self.dim = dim

        # load hpipm shared library
        __hpipm = CDLL("libhpipm.so")
        self.__hpipm = __hpipm

        # C qp struct
        qp_sol_struct_size = __hpipm.d_ocp_qp_sol_strsize()
        qp_sol_struct = cast(create_string_buffer(qp_sol_struct_size), c_void_p)
        self.qp_sol_struct = qp_sol_struct

        # C qp internal memory
        qp_sol_mem_size = __hpipm.d_ocp_qp_sol_memsize(dim.dim_struct)
        qp_sol_mem = cast(create_string_buffer(qp_sol_mem_size), c_void_p)
        self.qp_sol_mem = qp_sol_mem

        # create C qp
        __hpipm.d_ocp_qp_sol_create(dim.dim_struct, qp_sol_struct, qp_sol_mem)

        # getter functions for controls, states and slacks
        self.__getters = {
            "u": {
                "n_var": __hpipm.d_ocp_qp_dim_get_nu,
                "var": __hpipm.d_ocp_qp_sol_get_u,
            },
            "x": {
                "n_var": __hpipm.d_ocp_qp_dim_get_nx,
                "var": __hpipm.d_ocp_qp_sol_get_x,
            },
            "sl": {
                "n_var": __hpipm.d_ocp_qp_dim_get_ns,
                "var": __hpipm.d_ocp_qp_sol_get_sl,
            },
            "su": {
                "n_var": __hpipm.d_ocp_qp_dim_get_ns,
                "var": __hpipm.d_ocp_qp_sol_get_su,
            },
        }

    def get(self, field, idx_start, idx_end=None):
        if field not in self.__getters:
            raise NameError("hpipm_ocp_qp_sol.get: wrong field")
        else:
            return self.__get(self.__getters[field], idx_start, idx_end)

    def __get(self, getter, idx_start, idx_end=None):
        # number of variables
        n_var = np.zeros((1, 1), dtype=int)

        if idx_end is None:
            idx_end = idx_start

        var = []
        for i in range(idx_start, idx_end + 1):
            # get number of variables at stage
            tmp_ptr = cast(n_var.ctypes.data, POINTER(c_int))
            getter["n_var"](self.dim.dim_struct, i, tmp_ptr)

            var.append(np.zeros((n_var[0, 0], 1)))
            tmp_ptr = cast(var[-1].ctypes.data, POINTER(c_double))
            getter["var"](i, self.qp_sol_struct, tmp_ptr)

        return var if len(var) > 1 else var[0]

    def print_c_struct(self):
        self.__hpipm.d_ocp_qp_sol_print(self.dim.dim_struct, self.qp_sol_struct)
