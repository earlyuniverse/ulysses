dnl check demangling begin
dnl this checks if the demangling code using cxxabi in Error.cc is
dnl supported or not
dnl

dnl Note that according to BOOST (see the code in,
dnl /usr/include/boost/units/detail/utility.hpp or the link
dnl http://marc.info/?l=boost-bugs&m=122765675617968), instead of
dnl going through all this we could simply do
dnl   #if defined(__GLIBCXX__) || defined(__GLIBCPP__)
dnl
AC_DEFUN([ACX_CHECK_DEMANGLE_SUPPORT],
[
    dnl backup compiler/linker flags
    save_LIBS="$LIBS"
    save_LDFLAGS="$LDFLAGS"
    save_CXXFLAGS="$CXXFLAGS"
    save_CPPFLAGS="$CPPFLAGS"

    AC_LANG_PUSH(C++)
    dnl note the use of [[ ]] to allow brackets in the code
    AC_COMPILE_IFELSE([AC_LANG_SOURCE([[
#include <cstdio>     // sscanf
#include <execinfo.h> // for backtrace
#include <cstdlib>
#include <cxxabi.h>   // demangle

#include <iostream>   // cout used only in the test

using namespace std;

std::string demangle(const char* symbol) {
  size_t size;
  int status;
  char temp[128];
  char* demangled;
  //first, try to demangle a c++ name
  if (1 == sscanf(symbol, "%*[^(]%*[^_]%127[^)+]", temp)) {
    if (NULL != (demangled = abi::__cxa_demangle(temp, NULL, &size, &status))) {
      std::string result(demangled);
      free(demangled);
      return result;
    }
  }
  //if that didn't work, try to get a regular c symbol
  if (1 == sscanf(symbol, "%127s", temp)) {
    return temp;
  }
 
  //if all else fails, just return the symbol
  return symbol;
}

int main(void){ 
  void * array[10];
  char ** messages;

  int size = backtrace(array, 10);
  messages = backtrace_symbols(array, size);
  
  cout << "stack:" << endl;
  for (int i = 1; i < size && messages != NULL; ++i){
    cout << "  #" << i << ": " << messages[i] << endl;
    cout << "    " << demangle(messages[i]) << endl;
  }
  free(messages);

  return 0;

}
    ]])], [acx_demangle_support=yes],[acx_demangle_support=no])
    AC_LANG_POP(C++)

    dnl restore the original flags
    LIBS="$save_LIBS"
    LDFLAGS="$save_LDFLAGS"
    CXXFLAGS="$save_CXXFLAGS"
    CPPFLAGS="$save_CPPFLAGS"

    AC_MSG_CHECKING([demangling support])

    if test "$acx_demangle_support" == yes; then 
	AC_MSG_RESULT(yes);
	$1
    else
	AC_MSG_RESULT(no);
        $2
   fi]
)


dnl check demangling end
