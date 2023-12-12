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


class hpipm_dense_qp:
    def __init__(self, dim):
        # save dim internally
        self.dim = dim

        # load hpipm shared library
        __hpipm = CDLL("libhpipm.so")
        self.__hpipm = __hpipm

        # C qp struct
        qp_struct_size = __hpipm.d_dense_qp_strsize()
        qp_struct = cast(create_string_buffer(qp_struct_size), c_void_p)
        self.qp_struct = qp_struct

        # C qp internal memory
        qp_mem_size = __hpipm.d_dense_qp_memsize(dim.dim_struct)
        qp_mem = cast(create_string_buffer(qp_mem_size), c_void_p)
        self.qp_mem = qp_mem

        # create C qp
        __hpipm.d_dense_qp_create(dim.dim_struct, qp_struct, qp_mem)

    def set(self, field, value):
        # cast to np array
        if type(value) is not np.ndarray:
            if isinstance(value, int) or isinstance(value, float):
                value_ = value
                value = np.array((1,))
                value[0] = value_

        # convert into column-major
        value_cm = np.ravel(value, "F")
        if field == "idxb" or field == "idxs" or field == "idxs_rev":
            value_cm = np.ascontiguousarray(value_cm, dtype=np.int32)
            tmp = cast(value_cm.ctypes.data, POINTER(c_int))
        else:
            value_cm = np.ascontiguousarray(value_cm, dtype=np.float64)
            tmp = cast(value_cm.ctypes.data, POINTER(c_double))
        field_name_b = field.encode("utf-8")
        self.__hpipm.d_dense_qp_set(c_char_p(field_name_b), tmp, self.qp_struct)

    def print_c_struct(self):
        self.__hpipm.d_dense_qp_print(self.dim.dim_struct, self.qp_struct)

    def codegen(self, file_name, mode):
        file_name_b = file_name.encode("utf-8")
        mode_b = mode.encode("utf-8")
        self.__hpipm.d_dense_qp_codegen(
            c_char_p(file_name_b), c_char_p(mode_b), self.dim.dim_struct, self.qp_struct
        )
