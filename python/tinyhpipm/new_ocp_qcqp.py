from ctypes import (
    CDLL,
    c_void_p,
    POINTER,
    c_int,
    c_char_p,
    sizeof,
    cast,
    create_string_buffer,
)
from typing import Union, Optional
import numpy as np
from enum import Enum
from tinyhpipm.common import hpipm_lib_name


__all__ = ["Dim", "Solver"]


class Dim:
    __hpipm: CDLL
    __dim_struct: c_void_p
    __dim_mem: c_void_p
    N: int

    def __init__(self, N: int = 1):
        self.__hpipm = CDLL(hpipm_lib_name)
        self.__dim_struct = cast(
            create_string_buffer(self.__hpipm.d_ocp_qcqp_dim_strsize()), c_void_p
        )
        self.__dim_mem = cast(
            create_string_buffer(self.__hpipm.d_ocp_qcqp_dim_memsize()), c_void_p
        )
        assert N >= 1
        self.__hpipm.d_ocp_qcqp_dim_create(N, self.__dim_struct, self.__dim_mem)
        self.N = N
        self.__allowed_fields = {
            "nx",
            "nu",
            "nbx",
            "nbu",
            "ng",
            "nq",
            "ns",
            "nsbx",
            "nsbu",
            "nsg",
            "nsq",
            "nbxe",
            "nbue",
            "nge",
            "nqe",
        }

    def set(
        self,
        field: str,
        value: int,
        idx_start: Optional[int] = None,
        idx_stop: Optional[int] = None,
    ):
        assert field in self.__allowed_fields, f'field "{field}" not allowed'
        assert isinstance(value, int), f"value must be an integer, got {type(value)}"

        self.__hpipm.d_ocp_qcqp_dim_set.argtypes = [c_char_p, c_int, c_int, c_void_p]
        field_name_b = c_char_p(field.encode("utf-8"))
        if idx_start is None:
            idx_start = 0
            idx_stop = self.N
        if idx_stop is None:
            idx_stop = idx_start + 1
        for i in range(idx_start, idx_stop):
            self.__hpipm.d_ocp_qcqp_dim_set(field_name_b, i, value, self.__dim_struct)

    def condense(self, Ncondensed: int):
        """
        Created a new Dim object for partial condensation
        :param Ncondensed: number of stages in condensed problem
        """
        assert Ncondensed <= self.N
        # TODO: implement this

    def codegen(self, file_name: str, open_mode: str = "a"):
        """
        Generates the C code for the problem
        """
        file_name_b = file_name.encode("utf-8")
        mode_b = open_mode.encode("utf-8")
        self.__hpipm.d_ocp_qcqp_dim_codegen(
            c_char_p(file_name_b), c_char_p(mode_b), self.__dim_struct
        )


class Solver:
    class Mode(Enum):
        SPEED = 1
        BALANCE = 2
        ROBUST = 3

    __hpipm: CDLL
    mode: Mode
    __dim: Dim
    __dim_condensed: Optional[Dim]
    __qcqp_struct: c_void_p
    __qcqp_mem: c_void_p
    __ipm_arg_struct: c_void_p
    __ipm_arg_mem: c_void_p
    __ipm_ws_struct: c_void_p
    __ipm_ws_mem: c_void_p
    __sol_struct: c_void_p
    __sol_mem: c_void_p

    def __init__(
        self,
        dim: Dim,
        dim_condensed: Optional[Dim],
        mode: Mode = Mode.ROBUST,
    ) -> None:
        self.__hpipm = CDLL(hpipm_lib_name)
        self.__dim = dim
        self.__dim_condensed = dim_condensed
        self.mode = mode

        # create qcqp
        self.__qcqp_struct = cast(
            create_string_buffer(self.__hpipm.d_ocp_qcqp_strsize()), c_void_p
        )
        self.__qcqp_mem = cast(
            create_string_buffer(self.__hpipm.d_ocp_qcqp_memsize()), c_void_p
        )
        self.__hpipm.d_ocp_qcqp_create(
            self.__dim.__dim_struct, self.__qcqp_struct, self.__qcqp_mem
        )
        # create ipm arg
        self.__ipm_arg_struct = cast(
            create_string_buffer(self.__hpipm.d_ocp_qcqp_ipm_arg_strsize()), c_void_p
        )
        self.__ipm_arg_mem = cast(
            create_string_buffer(self.__hpipm.d_ocp_qcqp_ipm_arg_memsize()), c_void_p
        )
        self.__hpipm.d_ocp_qcqp_ipm_arg_create(
            self.__dim.__dim_struct, self.__ipm_arg_struct, self.__ipm_arg_mem
        )
        self.__hpipm.d_ocp_qcqp_ipm_arg_set_default(mode.value, self.__ipm_arg_struct)

        # create sol
        self.__sol_struct = cast(
            create_string_buffer(self.__hpipm.d_ocp_qcqp_sol_strsize()), c_void_p
        )
        self.__sol_mem = cast(
            create_string_buffer(self.__hpipm.d_ocp_qcqp_sol_memsize()), c_void_p
        )
        self.__hpipm.d_ocp_qcqp_sol_create(
            self.__dim.__dim_struct, self.__sol_struct, self.__sol_mem
        )

        # create pointers for ipm_ws
        self.__ipm_ws_struct = cast(
            create_string_buffer(self.__hpipm.d_ocp_qcqp_ipm_ws_strsize()), c_void_p
        )
        self.__ipm_ws_mem = cast(
            create_string_buffer(self.__hpipm.d_ocp_qcqp_ipm_ws_memsize()), c_void_p
        )

        if self.__dim_condensed is not None:
            # create qcqp condensed
            self.__qcqp_condensed_struct = cast(
                create_string_buffer(self.__hpipm.d_ocp_qcqp_strsize()), c_void_p
            )
            self.__qcqp_condensed_mem = cast(
                create_string_buffer(self.__hpipm.d_ocp_qcqp_memsize()), c_void_p
            )
            self.__hpipm.d_ocp_qcqp_create(
                self.__dim_condensed.__dim_struct,
                self.__qcqp_condensed_struct,
                self.__qcqp_condensed_mem,
            )
            # create sol condensed
            self.__sol_struct = cast(
                create_string_buffer(self.__hpipm.d_ocp_qcqp_sol_strsize()), c_void_p
            )
            self.__sol_mem = cast(
                create_string_buffer(self.__hpipm.d_ocp_qcqp_sol_memsize()), c_void_p
            )
            self.__hpipm.d_ocp_qcqp_sol_create(
                self.__dim.__dim_struct, self.__sol_struct, self.__sol_mem
            )

            # compute the necessary block size
            self.__block_size = cast(
                create_string_buffer((self.__dim.N + 1) * sizeof(c_int)), POINTER(c_int)
            )
            self.__hpipm.d_part_cond_qcqp_compute_block_size(
                self.__dim.N, self.__dim_condensed.N, self.__block_size
            )

            # create part_cond_arg
            self.__part_cond_arg_struct = cast(
                create_string_buffer(self.__hpipm.d_part_cond_qcqp_arg_strsize()),
                c_void_p,
            )
            self.__part_cond_arg_mem = cast(
                create_string_buffer(self.__hpipm.d_part_cond_qcqp_arg_memsize()),
                c_void_p,
            )
            self.__hpipm.d_part_cond_qcqp_arg_create(
                self.__dim_condensed.N,
                self.__part_cond_arg_struct,
                self.__part_cond_arg_mem,
            )
            self.__hpipm.d_part_cond_qcqp_arg_set_default(self.__part_cond_arg_struct)

            # create part_cond_ws
            self.__part_cond_ws_struct = cast(
                create_string_buffer(self.__hpipm.d_part_cond_qcqp_ws_strsize()),
                c_void_p,
            )
            self.__part_cond_ws_mem = cast(
                create_string_buffer(self.__hpipm.d_part_cond_qcqp_ws_memsize()),
                c_void_p,
            )
            self.__hpipm.d_part_cond_qcqp_ws_create(
                self.__dim,
                self.__block_size,
                self.__dim_condensed,
                self.__part_cond_arg_struct,
                self.__part_cond_ws_struct,
                self.__part_cond_ws_mem,
            )

            # create ipm_ws from condensed dim
            self.__hpipm.d_ocp_qcqp_ipm_ws_create(
                self.__dim_condensed.__dim_struct,
                self.__ipm_arg_struct,
                self.__ipm_ws_struct,
                self.__ipm_ws_mem,
            )
        else:
            # create ipm ws from regular dim
            self.__hpipm.d_ocp_qcqp_ipm_ws_create(
                self.__dim.__dim_struct,
                self.__ipm_arg_struct,
                self.__ipm_ws_struct,
                self.__ipm_ws_mem,
            )

        self.__allowed_fields = {
            # qcqp
            # ipm_arg
            # qqcp_sol
            # part_cond_arg
        }

    def set(
        self,
        field: str,
        value: Union[float, np.ndarray],
        idx_start: Optional[int] = None,
        idx_stop: Optional[int] = None,
    ):
        """
        :param field: name of the field to set could be the qcqp, ipm arg,
        """
        pass

    def get(
        self,
        field: str,
        value: Union[float, np.ndarray],
        idx_start: Optional[int] = None,
        idx_stop: Optional[int] = None,
    ):
        """
        :param field: name of the field to set could be the qcqp, ipm arg,
        """
        pass

    def condense(self):
        """
        Performs the partial condensation of the problem. Only should be called if problem data has changed.
        """
        self.__hpipm.d_part_cond_qcqp_cond(
            self.__qcqp_struct,
            self.__qcqp_condensed_struct,
            self.__part_cond_arg_struct,
            self.__part_cond_ws_struct,
        )

    def solve(self) -> tuple[np.ndarray, np.ndarray]:
        """
        Solve the QCQP and returns for convenience the solution
        """
        x = np.empty((0, 3))
        u = np.empty((0, 2))
        return x, u
