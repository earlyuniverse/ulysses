
// actual physical parameters:
//
// coner
// epsilon
// ovlim

extern "C" {
#ifdef WIN32
  void _stdcall PXCONE 
#else
  void          pxcone_
#endif
  (
    const int    &  mode   ,    // 1=>e+e-, 2=>hadron-hadron
    const int    &  ntrak  ,    // Number of particles
    const int    &  itkdm  ,    // First dimension of PTRAK array: 
    const double *  ptrak  ,    // Array of particle 4-momenta (Px,Py,Pz,E)
    const double &  coner  ,    // Cone size (half angle) in radians
    const double &  epslon ,    // Minimum Jet energy (GeV)
    const double &  ovlim  ,    // Maximum fraction of overlap energy in a jet
    const int    &  mxjet  ,    // Maximum possible number of jets
          int    &  njet   ,    // Number of jets found
          double *  pjet,  // 5-vectors of jets
          int    *  ipass,    // Particle k belongs to jet number IPASS(k)-1
                                // IPASS = -1 if not assosciated to a jet
          int    *  ijmul,    // Jet i contains IJMUL[i] particles
          int    &  ierr        // = 0 if all is OK ;   = -1 otherwise
    );
}

#ifdef WIN32
#define pxcone PXCONE
#else
#define pxcone pxcone_
#endif
