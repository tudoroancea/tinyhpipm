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

from tinyhpipm.common import hpipm_lib_name


class hpipm_ocp_qcqp_solver:
    def __init__(self, qcqp_dims, arg):
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


class hpipm_ocp_qp:
    def __init__(self, dim):
        # save dim internally
        self.dim = dim

        # load hpipm shared library
        __hpipm = CDLL(hpipm_lib_name)
        self.__hpipm = __hpipm

        # C qp struct
        qp_struct_size = __hpipm.d_ocp_qp_strsize()
        qp_struct = cast(create_string_buffer(qp_struct_size), c_void_p)
        self.qp_struct = qp_struct

        # C qp internal memory
        qp_mem_size = __hpipm.d_ocp_qp_memsize(dim.dim_struct)
        qp_mem = cast(create_string_buffer(qp_mem_size), c_void_p)
        self.qp_mem = qp_mem

        # create C qp
        __hpipm.d_ocp_qp_create(dim.dim_struct, qp_struct, qp_mem)

    def set(self, field, value, idx_start, idx_end=None):
        # cast to np array
        if not isinstance(value, np.ndarray) and isinstance(value, (int, float)):
            value_ = value
            value = np.array((1,))
            value[0] = value_
        # convert into column-major
        value_cm = np.ravel(value, "F")
        # 		if(issubclass(value.dtype.type, np.integer)):
        if (
            field == "idxbx"
            or field == "idxbu"
            or field == "idxb"
            or field == "idxs"
            or field == "idxs_rev"
            or field == "idxe"
            or field == "idxbue"
            or field == "idxbxe"
            or field == "idxge"
        ):
            value_cm = np.ascontiguousarray(value_cm, dtype=np.int32)
            tmp = cast(value_cm.ctypes.data, POINTER(c_int))
        else:
            value_cm = np.ascontiguousarray(value_cm, dtype=np.float64)
            tmp = cast(value_cm.ctypes.data, POINTER(c_double))
        field_name_b = field.encode("utf-8")
        if idx_end is None:
            self.__hpipm.d_ocp_qp_set(
                c_char_p(field_name_b), idx_start, tmp, self.qp_struct
            )
        else:
            for i in range(idx_start, idx_end + 1):
                self.__hpipm.d_ocp_qp_set(
                    c_char_p(field_name_b), i, tmp, self.qp_struct
                )

    def print_c_struct(self):
        self.__hpipm.d_ocp_qp_print(self.dim.dim_struct, self.qp_struct)

    def codegen(self, file_name, mode):
        file_name_b = file_name.encode("utf-8")
        mode_b = mode.encode("utf-8")
        self.__hpipm.d_ocp_qp_codegen(
            c_char_p(file_name_b), c_char_p(mode_b), self.dim.dim_struct, self.qp_struct
        )


class hpipm_ocp_qp_sol:
    def __init__(self, dim):
        # save dim internally
        self.dim = dim

        # load hpipm shared library
        __hpipm = CDLL(hpipm_lib_name)
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


class hpipm_ocp_qp_solver_arg:
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
            raise NameError("hpipm_ocp_qp_solver_arg: wrong mode")

        # save dim internally
        self.dim = dim

        # load hpipm shared library
        __hpipm = CDLL(hpipm_lib_name)
        self.__hpipm = __hpipm

        # C qp struct
        arg_struct_size = __hpipm.d_ocp_qp_ipm_arg_strsize()
        arg_struct = cast(create_string_buffer(arg_struct_size), c_void_p)
        self.arg_struct = arg_struct

        # C qp internal memory
        arg_mem_size = __hpipm.d_ocp_qp_ipm_arg_memsize(dim.dim_struct)
        arg_mem = cast(create_string_buffer(arg_mem_size), c_void_p)
        self.arg_mem = arg_mem

        # create C qp
        __hpipm.d_ocp_qp_ipm_arg_create(dim.dim_struct, arg_struct, arg_mem)

        # initialize default arguments
        __hpipm.d_ocp_qp_ipm_arg_set_default(c_mode, arg_struct)  # mode==SPEED

    def set(self, field, value):
        if (
            (field == "mu0")
            | (field == "tol_stat")
            | (field == "tol_eq")
            | (field == "tol_ineq")
            | (field == "tol_comp")
            | (field == "reg_prim")
        ):
            tmp_in = np.zeros((1, 1))
            tmp_in[0][0] = value
            tmp = cast(tmp_in.ctypes.data, POINTER(c_double))
        elif (field == "iter_max") | (field == "pred_corr") | (field == "split_step"):
            tmp_in = np.zeros((1, 1), dtype=int)
            tmp_in[0][0] = value
            tmp = cast(tmp_in.ctypes.data, POINTER(c_int))
        else:
            raise NameError("hpipm_ocp_qp_solver_arg.set: wrong field")

        field_name_b = field.encode("utf-8")
        self.__hpipm.d_ocp_qp_ipm_arg_set(c_char_p(field_name_b), tmp, self.arg_struct)

    def codegen(self, file_name, mode):
        file_name_b = file_name.encode("utf-8")
        mode_b = mode.encode("utf-8")
        self.__hpipm.d_ocp_qp_ipm_arg_codegen(
            c_char_p(file_name_b),
            c_char_p(mode_b),
            self.dim.dim_struct,
            self.arg_struct,
        )


class hpipm_ocp_qp_solver:
    def __init__(self, qp_dims, arg):
        # load hpipm shared library
        __hpipm = CDLL(hpipm_lib_name)
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
