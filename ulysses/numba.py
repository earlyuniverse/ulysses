
try:
    import numba
    from numba import jit
    from numba.typed import List
except ImportError:
    msg="""
    Numba is not available.
    If present, we use its JIT compiler to accelerate computations.
    Try pip install numba --user.
    """
    print(msg)

    # Dummy functions to not break code if numba is not available
    def jit(func):
        return func
    List = list

