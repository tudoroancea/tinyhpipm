from sys import platform
from typing import Union


def dynamic_library_extension():
    if platform.startswith("linux"):
        return ".so"
    elif platform.startswith("darwin"):
        return ".dylib"
    elif platform.startswith("win32"):
        return ".dll"
    else:
        raise NameError("dynamic_library_extension: unknown platform")


tinyhpipm_lib_name = "libtinyhpipm_dynamic" + dynamic_library_extension()


FloatOrInt = Union[float, int]
