from cpython.exc cimport PyErr_NewException

# Create public exceptions -- these will be available for inclusion from an
# external C++ module using 'core.h' or similar, and also be available from
# within Cython.

# These are called from translate_yoda_error in errors.cpp

cdef public:
    # Master exception
    object YodaExc_Exception = PyErr_NewException(
        "yoda.Exception", <object> NULL, <object> NULL)

    # Inherited exceptions
    object YodaExc_BinningError = PyErr_NewException(
        "yoda.BinningError", YodaExc_Exception, <object> NULL)
    object YodaExc_RangeError = PyErr_NewException(
        "yoda.RangeError", YodaExc_Exception, <object> NULL)
    object YodaExc_LockError = PyErr_NewException(
        "yoda.LockError", YodaExc_Exception, <object> NULL)
    object YodaExc_GridError = PyErr_NewException(
        "yoda.GridError", YodaExc_Exception, <object> NULL)
    object YodaExc_LogicError = PyErr_NewException(
        "yoda.LogicError", YodaExc_Exception, <object> NULL)
    object YodaExc_LowStatsError = PyErr_NewException(
        "yoda.LowStatsError", YodaExc_Exception, <object> NULL)
    object YodaExc_WeightError = PyErr_NewException(
        "yoda.WeightError", YodaExc_Exception, <object> NULL)
    object YodaExc_AnnotationError = PyErr_NewException(
        "yoda.AnnotationError", YodaExc_Exception, <object> NULL)
    object YodaExc_ReadError = PyErr_NewException(
        "yoda.ReadError", YodaExc_Exception, <object> NULL)
    object YodaExc_UserError = PyErr_NewException(
        "yoda.UserError", YodaExc_Exception, <object> NULL)

# Note that these don't appear in python space due to the cdef. What we will
# have to do is use some kind of magic to make these appear in Python space.
