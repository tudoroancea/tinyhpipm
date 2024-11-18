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

from tinyhpipm.common import tinyhpipm_lib_name


class dense_qcqp:
    def __init__(self, dim):
        # save dim internally
        self.dim = dim

        # load tinyhpipm shared library
        __tinyhpipm = CDLL(tinyhpipm_lib_name)
        self.__tinyhpipm = __tinyhpipm

        # C qp struct
        qcqp_struct_size = __tinyhpipm.d_dense_qcqp_strsize()
        qcqp_struct = cast(create_string_buffer(qcqp_struct_size), c_void_p)
        self.qcqp_struct = qcqp_struct

        # C qp internal memory
        qcqp_mem_size = __tinyhpipm.d_dense_qcqp_memsize(dim.dim_struct)
        qcqp_mem = cast(create_string_buffer(qcqp_mem_size), c_void_p)
        self.qcqp_mem = qcqp_mem

        # create C qp
        __tinyhpipm.d_dense_qcqp_create(dim.dim_struct, qcqp_struct, qcqp_mem)

    def set(self, field, value):
        # cast to np array
        if not isinstance(value, np.ndarray) and isinstance(value, (int, float)):
            value_ = value
            value = np.array((1,))
            value[0] = value_
        # convert into column-major
        value_cm = np.ravel(value, "F")
        if field in {"idxb", "idxs", "idxs_rev"}:
            value_cm = np.ascontiguousarray(value_cm, dtype=np.int32)
            tmp = cast(value_cm.ctypes.data, POINTER(c_int))
        else:
            value_cm = np.ascontiguousarray(value_cm, dtype=np.float64)
            tmp = cast(value_cm.ctypes.data, POINTER(c_double))
        field_name_b = field.encode("utf-8")
        self.__tinyhpipm.d_dense_qcqp_set(c_char_p(field_name_b), tmp, self.qcqp_struct)

    def print_c_struct(self):
        self.__tinyhpipm.d_dense_qcqp_print(self.dim.dim_struct, self.qcqp_struct)

    def codegen(self, file_name, mode):
        file_name_b = file_name.encode("utf-8")
        mode_b = mode.encode("utf-8")
        self.__tinyhpipm.d_dense_qcqp_codegen(
            c_char_p(file_name_b),
            c_char_p(mode_b),
            self.dim.dim_struct,
            self.qcqp_struct,
        )


class dense_qcqp_solver:
    def __init__(self, qcqp_dims, arg):
        # load tinyhpipm shared library
        __tinyhpipm = CDLL(tinyhpipm_lib_name)
        self.__tinyhpipm = __tinyhpipm
        # set up ipm workspace struct
        sizeof_d_dense_qcqp_ipm_workspace = __tinyhpipm.d_dense_qcqp_ipm_ws_strsize()
        ipm_ws_struct = cast(
            create_string_buffer(sizeof_d_dense_qcqp_ipm_workspace), c_void_p
        )
        self.ipm_ws_struct = ipm_ws_struct
        # allocate memory for ipm workspace
        ipm_size = __tinyhpipm.d_dense_qcqp_ipm_ws_memsize(
            qcqp_dims.dim_struct, arg.arg_struct
        )
        ipm_ws_mem = cast(create_string_buffer(ipm_size), c_void_p)
        self.ipm_ws_mem = ipm_ws_mem
        # create C ws
        __tinyhpipm.d_dense_qcqp_ipm_ws_create(
            qcqp_dims.dim_struct, arg.arg_struct, ipm_ws_struct, ipm_ws_mem
        )
        self.arg = arg
        self.dim_struct = qcqp_dims.dim_struct

    def solve(self, qcqp, qcqp_sol):
        self.__tinyhpipm.d_dense_qcqp_ipm_solve(
            qcqp.qcqp_struct,
            qcqp_sol.qcqp_sol_struct,
            self.arg.arg_struct,
            self.ipm_ws_struct,
        )

    def get(self, field):
        if field == "stat":
            # get iters
            iters = np.zeros((1, 1), dtype=int)
            tmp = cast(iters.ctypes.data, POINTER(c_int))
            self.__tinyhpipm.d_dense_qcqp_ipm_get_iter(self.ipm_ws_struct, tmp)
            # get stat_m
            stat_m = np.zeros((1, 1), dtype=int)
            tmp = cast(stat_m.ctypes.data, POINTER(c_int))
            self.__tinyhpipm.d_dense_qcqp_ipm_get_stat_m(self.ipm_ws_struct, tmp)
            # get stat pointer
            res = np.zeros((iters[0][0] + 1, stat_m[0][0]))
            ptr = c_void_p()
            self.__tinyhpipm.d_dense_qcqp_ipm_get_stat(self.ipm_ws_struct, byref(ptr))
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
            raise NameError("tinyhpipm_dense_qcqp_solver.get: wrong field")
        field_name_b = field.encode("utf-8")
        self.__tinyhpipm.d_dense_qcqp_ipm_get(
            c_char_p(field_name_b), self.ipm_ws_struct, tmp
        )
        return res[0][0]


class dense_qcqp_solver_arg:
    def __init__(self, dim, mode):
        c_mode = 0
        if mode == "speed_abs":
            c_mode = 0
        elif mode == "speed":
            c_mode = 1
        elif mode == "balance":
            c_mode = 2
        elif mode == "robust":
            c_mode = 3
        else:
            raise NameError("tinyhpipm_dense_qcqp_solver_arg: wrong mode")

        # save dim internally
        self.dim = dim

        # load tinyhpipm shared library
        __tinyhpipm = CDLL(tinyhpipm_lib_name)
        self.__tinyhpipm = __tinyhpipm

        # C qcqp struct
        arg_struct_size = __tinyhpipm.d_dense_qcqp_ipm_arg_strsize()
        arg_struct = cast(create_string_buffer(arg_struct_size), c_void_p)
        self.arg_struct = arg_struct

        # C qcqp internal memory
        arg_mem_size = __tinyhpipm.d_dense_qcqp_ipm_arg_memsize(dim.dim_struct)
        arg_mem = cast(create_string_buffer(arg_mem_size), c_void_p)
        self.arg_mem = arg_mem

        # create C qcqp
        __tinyhpipm.d_dense_qcqp_ipm_arg_create(dim.dim_struct, arg_struct, arg_mem)

        # initialize default arguments
        __tinyhpipm.d_dense_qcqp_ipm_arg_set_default(c_mode, arg_struct)  # mode==SPEED

    def set(self, field, value):
        if field in {
            "mu0",
            "tol_stat",
            "tol_eq",
            "tol_ineq",
            "tol_comp",
            "reg_prim",
            "reg_dual",
        }:
            tmp_in = np.zeros((1, 1))
            tmp_in[0][0] = value
            tmp = cast(tmp_in.ctypes.data, POINTER(c_double))
        elif field in {"iter_max", "pred_corr", "split_step", "warm_start"}:
            tmp_in = np.zeros((1, 1), dtype=int)
            tmp_in[0][0] = value
            tmp = cast(tmp_in.ctypes.data, POINTER(c_int))
        else:
            raise NameError("tinyhpipm_dense_qcqp_solver_arg.set: wrong field")
        field_name_b = field.encode("utf-8")
        self.__tinyhpipm.d_dense_qcqp_ipm_arg_set(
            c_char_p(field_name_b), tmp, self.arg_struct
        )

    def codegen(self, file_name, mode):
        file_name_b = file_name.encode("utf-8")
        mode_b = mode.encode("utf-8")
        self.__tinyhpipm.d_dense_qcqp_ipm_arg_codegen(
            c_char_p(file_name_b),
            c_char_p(mode_b),
            self.dim.dim_struct,
            self.arg_struct,
        )


class dense_qcqp_sol:
    def __init__(self, dim):
        # save dim internally
        self.dim = dim

        # load tinyhpipm shared library
        __tinyhpipm = CDLL(tinyhpipm_lib_name)
        self.__tinyhpipm = __tinyhpipm

        # C qcqp struct
        qcqp_sol_struct_size = __tinyhpipm.d_dense_qcqp_sol_strsize()
        qcqp_sol_struct = cast(create_string_buffer(qcqp_sol_struct_size), c_void_p)
        self.qcqp_sol_struct = qcqp_sol_struct

        # C qp internal memory
        qcqp_sol_mem_size = __tinyhpipm.d_dense_qcqp_sol_memsize(dim.dim_struct)
        qcqp_sol_mem = cast(create_string_buffer(qcqp_sol_mem_size), c_void_p)
        self.qcqp_sol_mem = qcqp_sol_mem

        # create C qcqp
        __tinyhpipm.d_dense_qcqp_sol_create(
            dim.dim_struct, qcqp_sol_struct, qcqp_sol_mem
        )

        # getter functions for primals, duals, and slacks
        self.__getters = {
            "v": {
                "n_var": __tinyhpipm.d_dense_qcqp_dim_get_nv,
                "var": __tinyhpipm.d_dense_qcqp_sol_get_v,
            },
            "pi": {
                "n_var": __tinyhpipm.d_dense_qcqp_dim_get_ne,
                "var": __tinyhpipm.d_dense_qcqp_sol_get_pi,
            },
            "lam_lb": {
                "n_var": __tinyhpipm.d_dense_qcqp_dim_get_nb,
                "var": __tinyhpipm.d_dense_qcqp_sol_get_lam_lb,
            },
            "lam_ub": {
                "n_var": __tinyhpipm.d_dense_qcqp_dim_get_nb,
                "var": __tinyhpipm.d_dense_qcqp_sol_get_lam_ub,
            },
            "lam_lg": {
                "n_var": __tinyhpipm.d_dense_qcqp_dim_get_ng,
                "var": __tinyhpipm.d_dense_qcqp_sol_get_lam_lg,
            },
            "lam_ug": {
                "n_var": __tinyhpipm.d_dense_qcqp_dim_get_ng,
                "var": __tinyhpipm.d_dense_qcqp_sol_get_lam_ug,
            },
            "lam_uq": {
                "n_var": __tinyhpipm.d_dense_qcqp_dim_get_nq,
                "var": __tinyhpipm.d_dense_qcqp_sol_get_lam_uq,
            },
        }

    def get(self, field):
        if field not in self.__getters:
            raise NameError("tinyhpipm_dense_qcqp_sol.get: wrong field")
        else:
            return self.__get(self.__getters[field])

    def __get(self, getter):
        # number of variables
        n_var = np.zeros((1, 1), dtype=int)

        tmp_ptr = cast(n_var.ctypes.data, POINTER(c_int))
        getter["n_var"](self.dim.dim_struct, tmp_ptr)

        var = np.zeros((n_var[0, 0], 1))
        tmp_ptr = cast(var.ctypes.data, POINTER(c_double))
        getter["var"](self.qcqp_sol_struct, tmp_ptr)

        return var

    def set(self, field, value):
        if field == "v":
            value = np.ascontiguousarray(value, dtype=np.float64)
            tmp = cast(value.ctypes.data, POINTER(c_double))
            self.__tinyhpipm.d_dense_qcqp_sol_set_v(tmp, self.qcqp_sol_struct)
        else:
            raise NameError("tinyhpipm_dense_qcqp_sol.set: wrong field")

    def print_c_struct(self):
        self.__tinyhpipm.d_dense_qcqp_sol_print(
            self.dim.dim_struct, self.qcqp_sol_struct
        )


class dense_qcqp_dim:
    def __init__(self):
        # load tinyhpipm shared library
        __tinyhpipm = CDLL(tinyhpipm_lib_name)
        self.__tinyhpipm = __tinyhpipm

        # C dim struct
        dim_struct_size = __tinyhpipm.d_dense_qcqp_dim_strsize()
        dim_struct = cast(create_string_buffer(dim_struct_size), c_void_p)
        self.dim_struct = dim_struct

        # C dim internal memory
        dim_mem_size = __tinyhpipm.d_dense_qcqp_dim_memsize()
        dim_mem = cast(create_string_buffer(dim_mem_size), c_void_p)
        self.dim_mem = dim_mem

        # create C dim
        __tinyhpipm.d_dense_qcqp_dim_create(self.dim_struct, self.dim_mem)

    def set(self, field, value):
        self.__tinyhpipm.d_dense_qcqp_dim_set.argtypes = [c_char_p, c_int, c_void_p]
        field_name_b = field.encode("utf-8")
        self.__tinyhpipm.d_dense_qcqp_dim_set(
            c_char_p(field_name_b), value, self.dim_struct
        )

    def print_c_struct(self):
        self.__tinyhpipm.d_dense_qcqp_dim_print(self.dim_struct)

    def codegen(self, file_name, mode):
        file_name_b = file_name.encode("utf-8")
        mode_b = mode.encode("utf-8")
        self.__tinyhpipm.d_dense_qcqp_dim_codegen(
            c_char_p(file_name_b), c_char_p(mode_b), self.dim_struct
        )
