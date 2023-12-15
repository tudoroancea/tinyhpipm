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
from typing import Union

import numpy as np

from tinyhpipm.common import FloatOrInt, hpipm_lib_name


class hpipm_ocp_qcqp_dim:
    def __init__(self, N: int):
        # load hpipm shared library
        __hpipm = CDLL(hpipm_lib_name)
        self.__hpipm = __hpipm

        # C dim struct
        dim_struct_size = __hpipm.d_ocp_qcqp_dim_strsize()
        dim_struct = cast(create_string_buffer(dim_struct_size), c_void_p)
        self.dim_struct = dim_struct

        # C dim internal memory
        dim_mem_size = __hpipm.d_ocp_qcqp_dim_memsize(N)
        dim_mem = cast(create_string_buffer(dim_mem_size), c_void_p)
        self.dim_mem = dim_mem

        # create C dim
        __hpipm.d_ocp_qcqp_dim_create(N, self.dim_struct, self.dim_mem)

    def set(self, field: str, value: int, idx_start: int, idx_end: int = None):
        self.__hpipm.d_ocp_qcqp_dim_set.argtypes = [c_char_p, c_int, c_int, c_void_p]
        field_name_b = field.encode("utf-8")
        if idx_end is None:
            self.__hpipm.d_ocp_qcqp_dim_set(
                c_char_p(field_name_b), idx_start, value, self.dim_struct
            )
        else:
            for i in range(idx_start, idx_end + 1):
                self.__hpipm.d_ocp_qcqp_dim_set(
                    c_char_p(field_name_b), i, value, self.dim_struct
                )

    def print_c_struct(self):
        self.__hpipm.d_ocp_qcqp_dim_print(self.dim_struct)

    def codegen(self, file_name: str, mode: str):
        file_name_b = file_name.encode("utf-8")
        mode_b = mode.encode("utf-8")
        self.__hpipm.d_ocp_qcqp_dim_codegen(
            c_char_p(file_name_b), c_char_p(mode_b), self.dim_struct
        )


class hpipm_ocp_qcqp:
    def __init__(self, dim: hpipm_ocp_qcqp_dim):
        # save dim internally
        self.dim = dim

        # load hpipm shared library
        __hpipm = CDLL(hpipm_lib_name)
        self.__hpipm = __hpipm

        # C qp struct
        qp_struct_size = __hpipm.d_ocp_qcqp_strsize()
        qp_struct = cast(create_string_buffer(qp_struct_size), c_void_p)
        self.qp_struct = qp_struct

        # C qp internal memory
        qp_mem_size = __hpipm.d_ocp_qcqp_memsize(dim.dim_struct)
        qp_mem = cast(create_string_buffer(qp_mem_size), c_void_p)
        self.qp_mem = qp_mem

        # create C qp
        __hpipm.d_ocp_qcqp_create(dim.dim_struct, qp_struct, qp_mem)

    def set(
        self,
        field: str,
        value: Union[np.ndarray, int, float],
        idx_start: int,
        idx_end: int = None,
    ):
        # cast to np array
        if not isinstance(value, np.ndarray) and isinstance(value, (int, float)):
            value_ = value
            value = np.array((1,))
            value[0] = value_
        # convert into column-major
        value_cm = np.ravel(value, "F")
        # if ( field == "idxbx" or field == "idxbu" or field == "idxb" or field == "idxs" or field == "idxs_rev"):  # or field=='idxe' or field=='idxbue' or field=='idxbxe' or field=='idxge' or field=='idxqe'):
        if field.startswith("idx"):
            value_cm = np.ascontiguousarray(value_cm, dtype=np.int32)
            tmp = cast(value_cm.ctypes.data, POINTER(c_int))
        else:
            value_cm = np.ascontiguousarray(value_cm, dtype=np.float64)
            tmp = cast(value_cm.ctypes.data, POINTER(c_double))
        field_name_b = field.encode("utf-8")
        if idx_end is None:
            self.__hpipm.d_ocp_qcqp_set(
                c_char_p(field_name_b), idx_start, tmp, self.qp_struct
            )
        else:
            for i in range(idx_start, idx_end + 1):
                self.__hpipm.d_ocp_qcqp_set(
                    c_char_p(field_name_b), i, tmp, self.qp_struct
                )

    def print_c_struct(self):
        self.__hpipm.d_ocp_qcqp_print(self.dim.dim_struct, self.qp_struct)

    def codegen(self, file_name: str, mode: str):
        file_name_b = file_name.encode("utf-8")
        mode_b = mode.encode("utf-8")
        self.__hpipm.d_ocp_qcqp_codegen(
            c_char_p(file_name_b), c_char_p(mode_b), self.dim.dim_struct, self.qp_struct
        )


class hpipm_ocp_qcqp_sol:
    def __init__(self, dim: hpipm_ocp_qcqp_dim):
        # save dim internally
        self.dim = dim

        # load hpipm shared library
        __hpipm = CDLL(hpipm_lib_name)
        self.__hpipm = __hpipm

        # C qp struct
        qp_sol_struct_size = __hpipm.d_ocp_qcqp_sol_strsize()
        qp_sol_struct = cast(create_string_buffer(qp_sol_struct_size), c_void_p)
        self.qp_sol_struct = qp_sol_struct

        # C qp internal memory
        qp_sol_mem_size = __hpipm.d_ocp_qcqp_sol_memsize(dim.dim_struct)
        qp_sol_mem = cast(create_string_buffer(qp_sol_mem_size), c_void_p)
        self.qp_sol_mem = qp_sol_mem

        # create C qp
        __hpipm.d_ocp_qcqp_sol_create(dim.dim_struct, qp_sol_struct, qp_sol_mem)

    def get(self, field: str, idx_start: int, idx_end: int = None):
        if field == "u":
            vec = self.__get_u(idx_start, idx_end)
        elif field == "x":
            vec = self.__get_x(idx_start, idx_end)
        else:
            raise NameError("hpipm_ocp_qcqp_sol.get: wrong field")
        return vec

    def __get_u(self, idx_start: int, idx_end: int = None):
        nu = np.zeros((1, 1), dtype=int)
        if idx_end is None:
            # get nu at stage
            tmp_ptr = cast(nu.ctypes.data, POINTER(c_int))
            self.__hpipm.d_ocp_qcqp_dim_get_nu(self.dim.dim_struct, idx_start, tmp_ptr)
            u = np.zeros((nu[0, 0], 1))
            tmp_ptr = cast(u.ctypes.data, POINTER(c_double))
            self.__hpipm.d_ocp_qcqp_sol_get_u(idx_start, self.qp_sol_struct, tmp_ptr)
        else:
            u = []
            for i in range(idx_start, idx_end + 1):
                # get nu at stage
                tmp_ptr = cast(nu.ctypes.data, POINTER(c_int))
                self.__hpipm.d_ocp_qcqp_dim_get_nu(self.dim.dim_struct, i, tmp_ptr)
                u.append(np.zeros((nu[0, 0], 1)))
                tmp_ptr = cast(u[-1].ctypes.data, POINTER(c_double))
                self.__hpipm.d_ocp_qcqp_sol_get_u(i, self.qp_sol_struct, tmp_ptr)
        return u

    def __get_x(self, idx_start: int, idx_end: int = None):
        # nx
        nx = np.zeros((1, 1), dtype=int)
        if idx_end is None:
            # get nx at stage
            tmp_ptr = cast(nx.ctypes.data, POINTER(c_int))
            self.__hpipm.d_ocp_qcqp_dim_get_nx(self.dim.dim_struct, idx_start, tmp_ptr)
            x = np.zeros((nx[0, 0], 1))
            tmp_ptr = cast(x.ctypes.data, POINTER(c_double))
            self.__hpipm.d_ocp_qcqp_sol_get_x(idx_start, self.qp_sol_struct, tmp_ptr)
        else:
            x = []
            for i in range(idx_start, idx_end + 1):
                # get nx at stage
                tmp_ptr = cast(nx.ctypes.data, POINTER(c_int))
                self.__hpipm.d_ocp_qcqp_dim_get_nx(self.dim.dim_struct, i, tmp_ptr)
                x.append(np.zeros((nx[0, 0], 1)))
                tmp_ptr = cast(x[-1].ctypes.data, POINTER(c_double))
                self.__hpipm.d_ocp_qcqp_sol_get_x(i, self.qp_sol_struct, tmp_ptr)
        return x

    def print_c_struct(self):
        self.__hpipm.d_ocp_qcqp_sol_print(self.dim.dim_struct, self.qp_sol_struct)

    def set(
        self,
        field: str,
        value: Union[np.ndarray, float],
        idx_start: int,
        idx_end: int = None,
    ):
        value = np.atleast_1d(value)
        assert field in {
            "u",
            "x",
            "sl",
            "su",
        }, f"hpipm_ocp_qcqp_sol.set: wrong field {field}"
        setter_method = {
            "u": self.__hpipm.d_ocp_qcqp_sol_set_u,
            "x": self.__hpipm.d_ocp_qcqp_sol_set_x,
            "sl": self.__hpipm.d_ocp_qcqp_sol_set_sl,
            "su": self.__hpipm.d_ocp_qcqp_sol_set_su,
        }[field]
        if idx_end is None:
            idx_end = idx_start

        for i in range(idx_start, idx_end + 1):
            value_cm = np.ascontiguousarray(value.ravel("F"), dtype=np.float64)
            tmp = cast(value_cm.ctypes.data, POINTER(c_double))
            setter_method(i, tmp, self.qp_sol_struct)


class hpipm_ocp_qcqp_solver_arg:
    def __init__(self, dim: hpipm_ocp_qcqp_dim, mode: str):
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
            raise NameError("hpipm_ocp_qcqp_solver_arg: wrong mode")

        # save dim internally
        self.dim = dim

        # load hpipm shared library
        __hpipm = CDLL(hpipm_lib_name)
        self.__hpipm = __hpipm

        # C qp struct
        arg_struct_size = __hpipm.d_ocp_qcqp_ipm_arg_strsize()
        arg_struct = cast(create_string_buffer(arg_struct_size), c_void_p)
        self.arg_struct = arg_struct

        # C qp internal memory
        arg_mem_size = __hpipm.d_ocp_qcqp_ipm_arg_memsize(dim.dim_struct)
        arg_mem = cast(create_string_buffer(arg_mem_size), c_void_p)
        self.arg_mem = arg_mem

        # create C qp
        __hpipm.d_ocp_qcqp_ipm_arg_create(dim.dim_struct, arg_struct, arg_mem)

        # initialize default arguments
        __hpipm.d_ocp_qcqp_ipm_arg_set_default(c_mode, arg_struct)  # mode==SPEED

    def set(self, field: str, value: FloatOrInt):
        if field in {
            "mu0",
            "tol_stat",
            "tol_eq",
            "tol_ineq",
            "tol_comp",
            "reg_prim",
        }:
            tmp_in = np.zeros((1, 1))
            tmp_in[0][0] = value
            tmp = cast(tmp_in.ctypes.data, POINTER(c_double))
        elif field in {"iter_max", "warm_start", "pred_corr", "cond_pred_corr"}:
            tmp_in = np.zeros((1, 1), dtype=int)
            tmp_in[0][0] = value
            tmp = cast(tmp_in.ctypes.data, POINTER(c_int))
        else:
            raise NameError("hpipm_ocp_qcqp_solver_arg.set: wrong field")
        field_name_b = field.encode("utf-8")
        self.__hpipm.d_ocp_qcqp_ipm_arg_set(
            c_char_p(field_name_b), tmp, self.arg_struct
        )

    def codegen(self, file_name: str, mode: str):
        file_name_b = file_name.encode("utf-8")
        mode_b = mode.encode("utf-8")
        self.__hpipm.d_ocp_qcqp_ipm_arg_codegen(
            c_char_p(file_name_b),
            c_char_p(mode_b),
            self.dim.dim_struct,
            self.arg_struct,
        )


class hpipm_ocp_qcqp_solver:
    def __init__(self, dim: hpipm_ocp_qcqp_dim, arg: hpipm_ocp_qcqp_solver_arg):
        # load hpipm shared library
        __hpipm = CDLL(hpipm_lib_name)
        self.__hpipm = __hpipm

        # set up ipm workspace struct
        sizeof_d_ocp_qcqp_ipm_workspace = __hpipm.d_ocp_qcqp_ipm_ws_strsize()
        ipm_ws_struct = cast(
            create_string_buffer(sizeof_d_ocp_qcqp_ipm_workspace), c_void_p
        )
        self.ipm_ws_struct = ipm_ws_struct

        # allocate memory for ipm workspace
        ipm_size = __hpipm.d_ocp_qcqp_ipm_ws_memsize(dim.dim_struct, arg.arg_struct)
        ipm_ws_mem = cast(create_string_buffer(ipm_size), c_void_p)
        self.ipm_ws_mem = ipm_ws_mem

        # create C ws
        __hpipm.d_ocp_qcqp_ipm_ws_create(
            dim.dim_struct, arg.arg_struct, ipm_ws_struct, ipm_ws_mem
        )

        self.arg = arg
        self.dim_struct = dim.dim_struct

    def solve(self, qcqp: hpipm_ocp_qcqp, qp_sol: hpipm_ocp_qcqp_sol):
        self.__hpipm.d_ocp_qcqp_ipm_solve(
            qcqp.qp_struct,
            qp_sol.qp_sol_struct,
            self.arg.arg_struct,
            self.ipm_ws_struct,
        )

    def get(self, field: str):
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
