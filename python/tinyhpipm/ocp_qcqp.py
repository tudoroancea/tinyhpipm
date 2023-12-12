from ctypes import (
    CDLL,
    POINTER,
    c_char_p,
    c_double,
    c_int,
    c_void_p,
    cast,
    create_string_buffer,
)
import numpy as np


class hpipm_ocp_qcqp:
    def __init__(self, dim):
        # save dim internally
        self.dim = dim

        # load hpipm shared library
        __hpipm = CDLL("libhpipm.so")
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
        ):  # or field=='idxe' or field=='idxbue' or field=='idxbxe' or field=='idxge' or field=='idxqe'):
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

    def codegen(self, file_name, mode):
        file_name_b = file_name.encode("utf-8")
        mode_b = mode.encode("utf-8")
        self.__hpipm.d_ocp_qcqp_codegen(
            c_char_p(file_name_b), c_char_p(mode_b), self.dim.dim_struct, self.qp_struct
        )
