from ctypes import (
    CDLL,
    c_char_p,
    c_int,
    c_void_p,
    cast,
    create_string_buffer,
)


class hpipm_dense_qcqp_dim:
    def __init__(self):
        # load hpipm shared library
        __hpipm = CDLL("libhpipm.so")
        self.__hpipm = __hpipm

        # C dim struct
        dim_struct_size = __hpipm.d_dense_qcqp_dim_strsize()
        dim_struct = cast(create_string_buffer(dim_struct_size), c_void_p)
        self.dim_struct = dim_struct

        # C dim internal memory
        dim_mem_size = __hpipm.d_dense_qcqp_dim_memsize()
        dim_mem = cast(create_string_buffer(dim_mem_size), c_void_p)
        self.dim_mem = dim_mem

        # create C dim
        __hpipm.d_dense_qcqp_dim_create(self.dim_struct, self.dim_mem)

    def set(self, field, value):
        self.__hpipm.d_dense_qcqp_dim_set.argtypes = [c_char_p, c_int, c_void_p]
        field_name_b = field.encode("utf-8")
        self.__hpipm.d_dense_qcqp_dim_set(
            c_char_p(field_name_b), value, self.dim_struct
        )

    def print_c_struct(self):
        self.__hpipm.d_dense_qcqp_dim_print(self.dim_struct)

    def codegen(self, file_name, mode):
        file_name_b = file_name.encode("utf-8")
        mode_b = mode.encode("utf-8")
        self.__hpipm.d_dense_qcqp_dim_codegen(
            c_char_p(file_name_b), c_char_p(mode_b), self.dim_struct
        )
