#ifndef __FASTJET_CONTRIB_LUNDJSON_HH__
#define __FASTJET_CONTRIB_LUNDJSON_HH__

#include "LundGenerator.hh"

FASTJET_BEGIN_NAMESPACE

namespace contrib{

// declaration of helper function
void lund_elements_to_json(std::ostream & ostr, const LundDeclustering & d);

/// writes json to ostr for an individual declustering
inline std::ostream & lund_to_json(std::ostream & ostr, const LundDeclustering & d) {
  // /// we hardcode 6 digits of precision -- is this the right way to do things?
  // int prec_store = ostr.precision();
  // ostr.precision(6);

  ostr << "{";
  lund_elements_to_json(ostr, d);
  ostr << "}";

  return ostr;
}


/// writes json to ostr for a vector of declusterings
inline std::ostream & lund_to_json(std::ostream & ostr, const std::vector<LundDeclustering> & d) {
  ostr << "[";
  for (std::size_t i = 0; i < d.size(); i++) {
    if (i != 0) ostr << ",";
    lund_to_json(ostr, d[i]);
  } 
  ostr << "]";
  return ostr;
}

// /// writes a full Lund tree recursively to JSON
// inline std::ostream & lund_to_json(std::ostream & ostr, const LundGenerator & generator, const PseudoJet & jet) {
//   // we will use these below in the loop; declare them here to avoid
//   // repeated creation/destruction
//   PseudoJet p1,p2;
//   std::vector<LundDeclustering> d = generator(jet);
//   ostr << "[";
//   for (std::size_t i = 0; i < d.size(); i++) {
//     if (i != 0) ostr << ",";
//     ostr << "{";
//     lund_elements_to_json(ostr, d[i]);
//     // add in information on the leaf if it exists
//     if (d[i].softer().has_parents(p1,p2)) {
//       ostr << ",";
//       ostr << "\"leaf\":";
//       lund_to_json(ostr, generator, d[i].softer());
//     }
//     ostr << "}";
//   } 
//   ostr << "]";
//   return ostr;
// }

// helper function to write individual elements to json
inline void lund_elements_to_json(std::ostream & ostr, const LundDeclustering & d) {
  ostr << "\"p_pt\":" << d.pair()  .pt() << ",";
  ostr << "\"p_m\":"  << d.pair()  .m () << ",";
  ostr << "\"h_pt\":" << d.harder().pt() << ",";
  ostr << "\"s_pt\":" << d.softer().pt() << ",";
  ostr << "\"z\":"      << d.z()     << ",";
  ostr << "\"Delta\":"  << d.Delta() << ",";
  ostr << "\"kt\":"     << d.kt()    << ",";
  ostr << "\"psi\":"    << d.psi()         ;
}


} // namespace contrib

FASTJET_END_NAMESPACE

#endif // __FASTJET_CONTRIB_LUNDJSON_HH__
