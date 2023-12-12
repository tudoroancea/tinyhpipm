from sys import platform


def dynamic_library_extension():
    if platform.startswith("linux"):
        return ".so"
    elif platform.startswith("darwin"):
        return ".dylib"
    elif platform.startswith("win32"):
        return ".dll"
    else:
        raise NameError("dynamic_library_extension: unknown platform")


hpipm_lib_name = "libhpipm_dynamic" + dynamic_library_extension()
