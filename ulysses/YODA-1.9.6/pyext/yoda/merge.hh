template<class T> void cython_iadd(T* a, T* b) {
  *a += *b;
}

template<class T> void cython_isub(T* a, T* b) {
  *a -= *b;
}

// template<class T> void cython_imul_dbl(T* a, double b) {
//   *a *= b;
// }

// template<class T> void cython_idiv_dbl(T* a, double b) {
//   *a /= b;
// }

template<class T> T* cython_add(T* a, T* b) {
  return new T(*a + *b);
}

template<class T> T* cython_sub(T* a, T* b) {
  return new T(*a - *b);
}

template<class T> T* cython_div(T* a, T* b) {
  return new T(*a / *b);
}
