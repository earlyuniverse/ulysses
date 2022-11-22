#include <HepMC3/Attribute.h>
#include <HepMC3/Data/GenEventData.h>
#include <HepMC3/FourVector.h>
#include <HepMC3/GenCrossSection.h>
#include <HepMC3/GenEvent.h>
#include <HepMC3/GenHeavyIon.h>
#include <HepMC3/GenParticle.h>
#include <HepMC3/GenRunInfo.h>
#include <HepMC3/GenVertex.h>
#include <HepMC3/HEPEVT_Wrapper.h>
#include <functional>
#include <ios>
#include <map>
#include <memory>
#include <ostream>
#include <sstream> // __str__
#include <streambuf>
#include <string>
#include <utility>
#include <vector>

#include <pybind11/pybind11.h>
#include <functional>
#include <string>
#include <HepMC3/Version.h>
#include <HepMC3/Reader.h>
#include <HepMC3/Writer.h>
#include <HepMC3/Print.h>
#include <src/stl_binders.hpp>
#include <src/binders.h>


#ifndef BINDER_PYBIND11_TYPE_CASTER
	#define BINDER_PYBIND11_TYPE_CASTER
	PYBIND11_DECLARE_HOLDER_TYPE(T, std::shared_ptr<T>)
	PYBIND11_DECLARE_HOLDER_TYPE(T, T*)
	PYBIND11_MAKE_OPAQUE(std::shared_ptr<void>)
#endif

void bind_pyHepMC3_13(std::function< pybind11::module &(std::string const &namespace_) > &M)
{
	{ // HepMC3::HEPEVT_Wrapper file:HepMC3/HEPEVT_Wrapper.h line:96
		pybind11::class_<HepMC3::HEPEVT_Wrapper, std::shared_ptr<HepMC3::HEPEVT_Wrapper>> cl(M("HepMC3"), "HEPEVT_Wrapper", "");
		cl.def( pybind11::init( [](){ return new HepMC3::HEPEVT_Wrapper(); } ) );
		cl.def_static("zero_everything", (void (*)()) &HepMC3::HEPEVT_Wrapper::zero_everything, "Set all entries in HEPEVT to zero \n\nC++: HepMC3::HEPEVT_Wrapper::zero_everything() --> void");
		cl.def_static("GenEvent_to_HEPEVT", (bool (*)(const class HepMC3::GenEvent *)) &HepMC3::HEPEVT_Wrapper::GenEvent_to_HEPEVT, "Convert GenEvent to HEPEVT\n\nC++: HepMC3::HEPEVT_Wrapper::GenEvent_to_HEPEVT(const class HepMC3::GenEvent *) --> bool", pybind11::arg("evt"));
		cl.def_static("HEPEVT_to_GenEvent", (bool (*)(class HepMC3::GenEvent *)) &HepMC3::HEPEVT_Wrapper::HEPEVT_to_GenEvent, "Convert HEPEVT to GenEvent\n\nC++: HepMC3::HEPEVT_Wrapper::HEPEVT_to_GenEvent(class HepMC3::GenEvent *) --> bool", pybind11::arg("evt"));
		cl.def_static("fix_daughters", (bool (*)()) &HepMC3::HEPEVT_Wrapper::fix_daughters, "Tries to fix list of daughters \n\nC++: HepMC3::HEPEVT_Wrapper::fix_daughters() --> bool");
		cl.def_static("set_hepevt_address", (void (*)(char *)) &HepMC3::HEPEVT_Wrapper::set_hepevt_address, "C++: HepMC3::HEPEVT_Wrapper::set_hepevt_address(char *) --> void", pybind11::arg("c"));
		cl.def_static("max_number_entries", (int (*)()) &HepMC3::HEPEVT_Wrapper::max_number_entries, "C++: HepMC3::HEPEVT_Wrapper::max_number_entries() --> int");
		cl.def_static("event_number", (int (*)()) &HepMC3::HEPEVT_Wrapper::event_number, "C++: HepMC3::HEPEVT_Wrapper::event_number() --> int");
		cl.def_static("number_entries", (int (*)()) &HepMC3::HEPEVT_Wrapper::number_entries, "C++: HepMC3::HEPEVT_Wrapper::number_entries() --> int");
		cl.def_static("status", (int (*)(const int &)) &HepMC3::HEPEVT_Wrapper::status, "C++: HepMC3::HEPEVT_Wrapper::status(const int &) --> int", pybind11::arg("index"));
		cl.def_static("id", (int (*)(const int &)) &HepMC3::HEPEVT_Wrapper::id, "C++: HepMC3::HEPEVT_Wrapper::id(const int &) --> int", pybind11::arg("index"));
		cl.def_static("first_parent", (int (*)(const int &)) &HepMC3::HEPEVT_Wrapper::first_parent, "C++: HepMC3::HEPEVT_Wrapper::first_parent(const int &) --> int", pybind11::arg("index"));
		cl.def_static("last_parent", (int (*)(const int &)) &HepMC3::HEPEVT_Wrapper::last_parent, "C++: HepMC3::HEPEVT_Wrapper::last_parent(const int &) --> int", pybind11::arg("index"));
		cl.def_static("first_child", (int (*)(const int &)) &HepMC3::HEPEVT_Wrapper::first_child, "C++: HepMC3::HEPEVT_Wrapper::first_child(const int &) --> int", pybind11::arg("index"));
		cl.def_static("last_child", (int (*)(const int &)) &HepMC3::HEPEVT_Wrapper::last_child, "C++: HepMC3::HEPEVT_Wrapper::last_child(const int &) --> int", pybind11::arg("index"));
		cl.def_static("px", (double (*)(const int &)) &HepMC3::HEPEVT_Wrapper::px, "C++: HepMC3::HEPEVT_Wrapper::px(const int &) --> double", pybind11::arg("index"));
		cl.def_static("py", (double (*)(const int &)) &HepMC3::HEPEVT_Wrapper::py, "C++: HepMC3::HEPEVT_Wrapper::py(const int &) --> double", pybind11::arg("index"));
		cl.def_static("pz", (double (*)(const int &)) &HepMC3::HEPEVT_Wrapper::pz, "C++: HepMC3::HEPEVT_Wrapper::pz(const int &) --> double", pybind11::arg("index"));
		cl.def_static("e", (double (*)(const int &)) &HepMC3::HEPEVT_Wrapper::e, "C++: HepMC3::HEPEVT_Wrapper::e(const int &) --> double", pybind11::arg("index"));
		cl.def_static("m", (double (*)(const int &)) &HepMC3::HEPEVT_Wrapper::m, "C++: HepMC3::HEPEVT_Wrapper::m(const int &) --> double", pybind11::arg("index"));
		cl.def_static("x", (double (*)(const int &)) &HepMC3::HEPEVT_Wrapper::x, "C++: HepMC3::HEPEVT_Wrapper::x(const int &) --> double", pybind11::arg("index"));
		cl.def_static("y", (double (*)(const int &)) &HepMC3::HEPEVT_Wrapper::y, "C++: HepMC3::HEPEVT_Wrapper::y(const int &) --> double", pybind11::arg("index"));
		cl.def_static("z", (double (*)(const int &)) &HepMC3::HEPEVT_Wrapper::z, "C++: HepMC3::HEPEVT_Wrapper::z(const int &) --> double", pybind11::arg("index"));
		cl.def_static("t", (double (*)(const int &)) &HepMC3::HEPEVT_Wrapper::t, "C++: HepMC3::HEPEVT_Wrapper::t(const int &) --> double", pybind11::arg("index"));
		cl.def_static("number_parents", (int (*)(const int &)) &HepMC3::HEPEVT_Wrapper::number_parents, "C++: HepMC3::HEPEVT_Wrapper::number_parents(const int &) --> int", pybind11::arg("index"));
		cl.def_static("number_children", (int (*)(const int &)) &HepMC3::HEPEVT_Wrapper::number_children, "C++: HepMC3::HEPEVT_Wrapper::number_children(const int &) --> int", pybind11::arg("index"));
		cl.def_static("number_children_exact", (int (*)(const int &)) &HepMC3::HEPEVT_Wrapper::number_children_exact, "C++: HepMC3::HEPEVT_Wrapper::number_children_exact(const int &) --> int", pybind11::arg("index"));
		cl.def_static("set_event_number", (void (*)(const int &)) &HepMC3::HEPEVT_Wrapper::set_event_number, "C++: HepMC3::HEPEVT_Wrapper::set_event_number(const int &) --> void", pybind11::arg("evtno"));
		cl.def_static("set_number_entries", (void (*)(const int &)) &HepMC3::HEPEVT_Wrapper::set_number_entries, "C++: HepMC3::HEPEVT_Wrapper::set_number_entries(const int &) --> void", pybind11::arg("noentries"));
		cl.def_static("set_status", (void (*)(const int &, const int &)) &HepMC3::HEPEVT_Wrapper::set_status, "C++: HepMC3::HEPEVT_Wrapper::set_status(const int &, const int &) --> void", pybind11::arg("index"), pybind11::arg("status"));
		cl.def_static("set_id", (void (*)(const int &, const int &)) &HepMC3::HEPEVT_Wrapper::set_id, "C++: HepMC3::HEPEVT_Wrapper::set_id(const int &, const int &) --> void", pybind11::arg("index"), pybind11::arg("id"));
		cl.def_static("set_parents", (void (*)(const int &, const int &, const int &)) &HepMC3::HEPEVT_Wrapper::set_parents, "C++: HepMC3::HEPEVT_Wrapper::set_parents(const int &, const int &, const int &) --> void", pybind11::arg("index"), pybind11::arg("firstparent"), pybind11::arg("lastparent"));
		cl.def_static("set_children", (void (*)(const int &, const int &, const int &)) &HepMC3::HEPEVT_Wrapper::set_children, "C++: HepMC3::HEPEVT_Wrapper::set_children(const int &, const int &, const int &) --> void", pybind11::arg("index"), pybind11::arg("firstchild"), pybind11::arg("lastchild"));
		cl.def_static("set_momentum", (void (*)(const int &, const double &, const double &, const double &, const double &)) &HepMC3::HEPEVT_Wrapper::set_momentum, "C++: HepMC3::HEPEVT_Wrapper::set_momentum(const int &, const double &, const double &, const double &, const double &) --> void", pybind11::arg("index"), pybind11::arg("px"), pybind11::arg("py"), pybind11::arg("pz"), pybind11::arg("e"));
		cl.def_static("set_mass", (void (*)(const int &, double)) &HepMC3::HEPEVT_Wrapper::set_mass, "C++: HepMC3::HEPEVT_Wrapper::set_mass(const int &, double) --> void", pybind11::arg("index"), pybind11::arg("mass"));
		cl.def_static("set_position", (void (*)(const int &, const double &, const double &, const double &, const double &)) &HepMC3::HEPEVT_Wrapper::set_position, "C++: HepMC3::HEPEVT_Wrapper::set_position(const int &, const double &, const double &, const double &, const double &) --> void", pybind11::arg("index"), pybind11::arg("x"), pybind11::arg("y"), pybind11::arg("z"), pybind11::arg("t"));

		 binder::custom_HEPEVT_Wrapper_binder(cl);
	}
}
