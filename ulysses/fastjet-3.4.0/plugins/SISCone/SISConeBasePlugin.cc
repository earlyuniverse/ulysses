// fastjet stuff
#include "fastjet/SISConeBasePlugin.hh"

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

// returns true when the scale associated with jet a is larger than
// the scale associated with jet b
//
// By default this does a simple direct comparison but it can be
// overloaded for higher precision [recommended if possible]
bool SISConeBasePlugin::UserScaleBase::is_larger(const PseudoJet & a, const PseudoJet & b) const{
  return a.structure_of<UserScaleBase>().ordering_var2()
       > b.structure_of<UserScaleBase>().ordering_var2();
}

FASTJET_END_NAMESPACE
