## A base CppObject which prevents null pointer access
cdef class Base:
    cdef void* _ptr
    cdef bint _deallocate
    cdef object _owner

    cdef inline void* ptr(self) except NULL:
        if self._ptr == NULL:
            raise MemoryError('Null pointer referenced: perhaps the class is uninitialised.')
        return self._ptr


## Magic for setting pointers

## Use this for setting a pointer that is owned by the object obj
## e.g. if you want deallocation to happen automatically when this object goes
## out of scope
cdef inline set_owned_ptr(Base obj, void* ptr):
    obj._ptr = ptr
    obj._deallocate = True

## Use this for setting a pointer that is *not* owned by the object obj e.g. if
## you were given this object by a class which will deallocate it once it's
## finished
cdef inline set_borrowed_ptr(Base obj, void* ptr, object owner):
    obj._ptr = ptr
    obj._deallocate = False
    obj._owner = owner

## Use this to create a new object of type cls from the pointer ptr. The class is
## one which owns its pointer, and will deallocate when it's done. It's this one
## that you want to use for e.g. loaders or factory methods where the user is
## expected to manage memory. Or where you've explicitly called new.
cdef inline object new_owned_cls(object cls, void* ptr):
    obj = cls.__new__(cls)
    set_owned_ptr(obj, ptr)
    return obj

## Use this to create a thin wrapper around a pointer that will *not* be
## deallocated when it goes out of scope. Useful when you're given a reference
## from a class which is not expecting you to delete its innards!
cdef inline object new_borrowed_cls(object cls, void* ptr, object owner):
    obj = cls.__new__(cls)
    set_borrowed_ptr(obj, ptr, owner)
    return obj


## Translation utility to allow use of 'Pythonic' negative indices
cdef inline size_t pythonic_index(int i, size_t size) except ? 0:
    if i < 0:
        i += size
    if 0 <= i and i < size:
        return i
    else:
        raise IndexError


# ## Utils for handling error conversions to/from std::pair
# from libcpp.pair cimport pair

# cdef inline object read_edge_pair(pair[double, double] es):
#     return EdgePair(es.first, es.second)

# cdef inline object read_error_pair(pair[double, double] es):
#     return ErrorPair(es.first, es.second)

# cdef inline pair[double, double] read_symmetric(object val) except *:
#     try:
#         a, b = val
#     except TypeError:
#         a = b = val
#     return pair[double, double](a, b)
