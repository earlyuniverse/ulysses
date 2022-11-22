#include <HepMC3/Attribute.h>
#include <HepMC3/Data/GenEventData.h>
#include <HepMC3/Data/GenParticleData.h>
#include <HepMC3/FourVector.h>
#include <HepMC3/GenCrossSection.h>
#include <HepMC3/GenEvent.h>
#include <HepMC3/GenHeavyIon.h>
#include <HepMC3/GenParticle.h>
#include <HepMC3/GenRunInfo.h>
#include <HepMC3/GenVertex.h>
#include <HepMC3/WriterAscii.h>
#include <functional>
#include <ios>
#include <iterator>
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

// HepMC3::WriterAscii file:HepMC3/WriterAscii.h line:25
struct PyCallBack_HepMC3_WriterAscii : public HepMC3::WriterAscii {
	using HepMC3::WriterAscii::WriterAscii;

	void write_event(const class HepMC3::GenEvent & a0) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::WriterAscii *>(this), "write_event");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return WriterAscii::write_event(a0);
	}
	bool failed() override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::WriterAscii *>(this), "failed");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return WriterAscii::failed();
	}
	void close() override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::WriterAscii *>(this), "close");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return WriterAscii::close();
	}
};

void bind_pyHepMC3_11(std::function< pybind11::module &(std::string const &namespace_) > &M)
{
	{ // HepMC3::WriterAscii file:HepMC3/WriterAscii.h line:25
		pybind11::class_<HepMC3::WriterAscii, std::shared_ptr<HepMC3::WriterAscii>, PyCallBack_HepMC3_WriterAscii, HepMC3::Writer> cl(M("HepMC3"), "WriterAscii", "");
		cl.def( pybind11::init( [](const class std::basic_string<char> & a0){ return new HepMC3::WriterAscii(a0); }, [](const class std::basic_string<char> & a0){ return new PyCallBack_HepMC3_WriterAscii(a0); } ), "doc");
		cl.def( pybind11::init<const std::string &, class std::shared_ptr<class HepMC3::GenRunInfo>>(), pybind11::arg("filename"), pybind11::arg("run") );

		cl.def("write_event", (void (HepMC3::WriterAscii::*)(const class HepMC3::GenEvent &)) &HepMC3::WriterAscii::write_event, "Write event to file\n\n \n Event to be serialized\n\nC++: HepMC3::WriterAscii::write_event(const class HepMC3::GenEvent &) --> void", pybind11::arg("evt"));
		cl.def("write_run_info", (void (HepMC3::WriterAscii::*)()) &HepMC3::WriterAscii::write_run_info, "Write the GenRunInfo object to file.\n\nC++: HepMC3::WriterAscii::write_run_info() --> void");
		cl.def("failed", (bool (HepMC3::WriterAscii::*)()) &HepMC3::WriterAscii::failed, "Return status of the stream\n\nC++: HepMC3::WriterAscii::failed() --> bool");
		cl.def("close", (void (HepMC3::WriterAscii::*)()) &HepMC3::WriterAscii::close, "Close file stream\n\nC++: HepMC3::WriterAscii::close() --> void");
		cl.def("set_precision", (void (HepMC3::WriterAscii::*)(const int &)) &HepMC3::WriterAscii::set_precision, "Set output precision\n\n So far available range is [2,24]. Default is 16.\n\nC++: HepMC3::WriterAscii::set_precision(const int &) --> void", pybind11::arg("prec"));
		cl.def("precision", (int (HepMC3::WriterAscii::*)() const) &HepMC3::WriterAscii::precision, "Return output precision\n\nC++: HepMC3::WriterAscii::precision() const --> int");
	}
}
