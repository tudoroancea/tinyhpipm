from ctypes import (
    CDLL,
    cast,
    create_string_buffer,
    c_void_p,
    POINTER,
    c_int,
    c_double,
)
import numpy as np


class hpipm_dense_qp_sol:
    def __init__(self, dim):
        # save dim internally
        self.dim = dim

        # load hpipm shared library
        __hpipm = CDLL("libhpipm.so")
        self.__hpipm = __hpipm

        # C qp struct
        qp_sol_struct_size = __hpipm.d_dense_qp_sol_strsize()
        qp_sol_struct = cast(create_string_buffer(qp_sol_struct_size), c_void_p)
        self.qp_sol_struct = qp_sol_struct

        # C qp internal memory
        qp_sol_mem_size = __hpipm.d_dense_qp_sol_memsize(dim.dim_struct)
        qp_sol_mem = cast(create_string_buffer(qp_sol_mem_size), c_void_p)
        self.qp_sol_mem = qp_sol_mem

        # create C qp
        __hpipm.d_dense_qp_sol_create(dim.dim_struct, qp_sol_struct, qp_sol_mem)

        # getter functions for primals, duals, and slacks
        self.__getters = {
            "v": {
                "n_var": __hpipm.d_dense_qp_dim_get_nv,
                "var": __hpipm.d_dense_qp_sol_get_v,
            },
            "pi": {
                "n_var": __hpipm.d_dense_qp_dim_get_ne,
                "var": __hpipm.d_dense_qp_sol_get_pi,
            },
            "lam_lb": {
                "n_var": __hpipm.d_dense_qp_dim_get_nb,
                "var": __hpipm.d_dense_qp_sol_get_lam_lb,
            },
            "lam_ub": {
                "n_var": __hpipm.d_dense_qp_dim_get_nb,
                "var": __hpipm.d_dense_qp_sol_get_lam_ub,
            },
            "lam_lg": {
                "n_var": __hpipm.d_dense_qp_dim_get_ng,
                "var": __hpipm.d_dense_qp_sol_get_lam_lg,
            },
            "lam_ug": {
                "n_var": __hpipm.d_dense_qp_dim_get_ng,
                "var": __hpipm.d_dense_qp_sol_get_lam_ug,
            },
        }

    def get(self, field):
        if field not in self.__getters:
            raise NameError("hpipm_dense_qp_sol.get: wrong field")
        else:
            return self.__get(self.__getters[field])

    def __get(self, getter):
        # number of variables
        n_var = np.zeros((1, 1), dtype=int)

        tmp_ptr = cast(n_var.ctypes.data, POINTER(c_int))
        getter["n_var"](self.dim.dim_struct, tmp_ptr)

        var = np.zeros((n_var[0, 0], 1))
        tmp_ptr = cast(var.ctypes.data, POINTER(c_double))
        getter["var"](self.qp_sol_struct, tmp_ptr)

        return var

    def set(self, field, value):
        if field == "v":
            value = np.ascontiguousarray(value, dtype=np.float64)
            tmp = cast(value.ctypes.data, POINTER(c_double))
            self.__hpipm.d_dense_qp_sol_set_v(tmp, self.qp_sol_struct)
        else:
            raise NameError("hpipm_dense_qp_sol.set: wrong field")

    def print_c_struct(self):
        self.__hpipm.d_dense_qp_sol_print(self.dim.dim_struct, self.qp_sol_struct)
