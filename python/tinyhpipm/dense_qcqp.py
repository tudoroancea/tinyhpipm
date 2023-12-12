from ctypes import (
    CDLL,
    cast,
    create_string_buffer,
    c_void_p,
    POINTER,
    c_int,
    c_double,
    c_char_p,
)
import numpy as np


class hpipm_dense_qcqp:
    def __init__(self, dim):
        # save dim internally
        self.dim = dim

        # load hpipm shared library
        __hpipm = CDLL("libhpipm.so")
        self.__hpipm = __hpipm

        # C qp struct
        qcqp_struct_size = __hpipm.d_dense_qcqp_strsize()
        qcqp_struct = cast(create_string_buffer(qcqp_struct_size), c_void_p)
        self.qcqp_struct = qcqp_struct

        # C qp internal memory
        qcqp_mem_size = __hpipm.d_dense_qcqp_memsize(dim.dim_struct)
        qcqp_mem = cast(create_string_buffer(qcqp_mem_size), c_void_p)
        self.qcqp_mem = qcqp_mem

        # create C qp
        __hpipm.d_dense_qcqp_create(dim.dim_struct, qcqp_struct, qcqp_mem)

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
        self.__hpipm.d_dense_qcqp_set(c_char_p(field_name_b), tmp, self.qcqp_struct)

    def print_c_struct(self):
        self.__hpipm.d_dense_qcqp_print(self.dim.dim_struct, self.qcqp_struct)

    def codegen(self, file_name, mode):
        file_name_b = file_name.encode("utf-8")
        mode_b = mode.encode("utf-8")
        self.__hpipm.d_dense_qcqp_codegen(
            c_char_p(file_name_b),
            c_char_p(mode_b),
            self.dim.dim_struct,
            self.qcqp_struct,
        )
