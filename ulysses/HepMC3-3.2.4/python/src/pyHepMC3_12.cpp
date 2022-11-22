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
#include <HepMC3/ReaderAscii.h>
#include <HepMC3/ReaderAsciiHepMC2.h>
#include <HepMC3/ReaderHEPEVT.h>
#include <HepMC3/WriterAsciiHepMC2.h>
#include <HepMC3/WriterHEPEVT.h>
#include <functional>
#include <ios>
#include <istream>
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

// HepMC3::ReaderAscii file:HepMC3/ReaderAscii.h line:29
struct PyCallBack_HepMC3_ReaderAscii : public HepMC3::ReaderAscii {
	using HepMC3::ReaderAscii::ReaderAscii;

	bool skip(const int a0) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::ReaderAscii *>(this), "skip");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return ReaderAscii::skip(a0);
	}
	bool read_event(class HepMC3::GenEvent & a0) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::ReaderAscii *>(this), "read_event");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return ReaderAscii::read_event(a0);
	}
	bool failed() override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::ReaderAscii *>(this), "failed");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return ReaderAscii::failed();
	}
	void close() override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::ReaderAscii *>(this), "close");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return ReaderAscii::close();
	}
};

// HepMC3::WriterAsciiHepMC2 file:HepMC3/WriterAsciiHepMC2.h line:26
struct PyCallBack_HepMC3_WriterAsciiHepMC2 : public HepMC3::WriterAsciiHepMC2 {
	using HepMC3::WriterAsciiHepMC2::WriterAsciiHepMC2;

	void write_event(const class HepMC3::GenEvent & a0) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::WriterAsciiHepMC2 *>(this), "write_event");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return WriterAsciiHepMC2::write_event(a0);
	}
	bool failed() override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::WriterAsciiHepMC2 *>(this), "failed");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return WriterAsciiHepMC2::failed();
	}
	void close() override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::WriterAsciiHepMC2 *>(this), "close");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return WriterAsciiHepMC2::close();
	}
};

// HepMC3::ReaderAsciiHepMC2 file:HepMC3/ReaderAsciiHepMC2.h line:30
struct PyCallBack_HepMC3_ReaderAsciiHepMC2 : public HepMC3::ReaderAsciiHepMC2 {
	using HepMC3::ReaderAsciiHepMC2::ReaderAsciiHepMC2;

	bool skip(const int a0) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::ReaderAsciiHepMC2 *>(this), "skip");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return ReaderAsciiHepMC2::skip(a0);
	}
	bool read_event(class HepMC3::GenEvent & a0) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::ReaderAsciiHepMC2 *>(this), "read_event");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return ReaderAsciiHepMC2::read_event(a0);
	}
	bool failed() override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::ReaderAsciiHepMC2 *>(this), "failed");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return ReaderAsciiHepMC2::failed();
	}
	void close() override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::ReaderAsciiHepMC2 *>(this), "close");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return ReaderAsciiHepMC2::close();
	}
};

// HepMC3::WriterHEPEVT file:HepMC3/WriterHEPEVT.h line:27
struct PyCallBack_HepMC3_WriterHEPEVT : public HepMC3::WriterHEPEVT {
	using HepMC3::WriterHEPEVT::WriterHEPEVT;

	void write_hepevt_particle(int a0, bool a1) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::WriterHEPEVT *>(this), "write_hepevt_particle");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return WriterHEPEVT::write_hepevt_particle(a0, a1);
	}
	void write_hepevt_event_header() override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::WriterHEPEVT *>(this), "write_hepevt_event_header");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return WriterHEPEVT::write_hepevt_event_header();
	}
	void write_event(const class HepMC3::GenEvent & a0) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::WriterHEPEVT *>(this), "write_event");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return WriterHEPEVT::write_event(a0);
	}
	void close() override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::WriterHEPEVT *>(this), "close");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return WriterHEPEVT::close();
	}
	bool failed() override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::WriterHEPEVT *>(this), "failed");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return WriterHEPEVT::failed();
	}
};

// HepMC3::ReaderHEPEVT file:HepMC3/ReaderHEPEVT.h line:32
struct PyCallBack_HepMC3_ReaderHEPEVT : public HepMC3::ReaderHEPEVT {
	using HepMC3::ReaderHEPEVT::ReaderHEPEVT;

	bool read_hepevt_event_header() override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::ReaderHEPEVT *>(this), "read_hepevt_event_header");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return ReaderHEPEVT::read_hepevt_event_header();
	}
	bool read_hepevt_particle(int a0) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::ReaderHEPEVT *>(this), "read_hepevt_particle");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return ReaderHEPEVT::read_hepevt_particle(a0);
	}
	bool skip(const int a0) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::ReaderHEPEVT *>(this), "skip");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return ReaderHEPEVT::skip(a0);
	}
	bool read_event(class HepMC3::GenEvent & a0) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::ReaderHEPEVT *>(this), "read_event");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return ReaderHEPEVT::read_event(a0);
	}
	void close() override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::ReaderHEPEVT *>(this), "close");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return ReaderHEPEVT::close();
	}
	bool failed() override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::ReaderHEPEVT *>(this), "failed");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return ReaderHEPEVT::failed();
	}
};

void bind_pyHepMC3_12(std::function< pybind11::module &(std::string const &namespace_) > &M)
{
	{ // HepMC3::ReaderAscii file:HepMC3/ReaderAscii.h line:29
		pybind11::class_<HepMC3::ReaderAscii, std::shared_ptr<HepMC3::ReaderAscii>, PyCallBack_HepMC3_ReaderAscii, HepMC3::Reader> cl(M("HepMC3"), "ReaderAscii", "");
		cl.def( pybind11::init<const std::string &>(), pybind11::arg("filename") );

		cl.def("skip", (bool (HepMC3::ReaderAscii::*)(const int)) &HepMC3::ReaderAscii::skip, "skip events\n\nC++: HepMC3::ReaderAscii::skip(const int) --> bool", pybind11::arg(""));
		cl.def("read_event", (bool (HepMC3::ReaderAscii::*)(class HepMC3::GenEvent &)) &HepMC3::ReaderAscii::read_event, "Load event from file\n\n \n Event to be filled\n\nC++: HepMC3::ReaderAscii::read_event(class HepMC3::GenEvent &) --> bool", pybind11::arg("evt"));
		cl.def("failed", (bool (HepMC3::ReaderAscii::*)()) &HepMC3::ReaderAscii::failed, "Return status of the stream\n\nC++: HepMC3::ReaderAscii::failed() --> bool");
		cl.def("close", (void (HepMC3::ReaderAscii::*)()) &HepMC3::ReaderAscii::close, "Close file stream\n\nC++: HepMC3::ReaderAscii::close() --> void");
	}
	{ // HepMC3::WriterAsciiHepMC2 file:HepMC3/WriterAsciiHepMC2.h line:26
		pybind11::class_<HepMC3::WriterAsciiHepMC2, std::shared_ptr<HepMC3::WriterAsciiHepMC2>, PyCallBack_HepMC3_WriterAsciiHepMC2, HepMC3::Writer> cl(M("HepMC3"), "WriterAsciiHepMC2", "");
		cl.def( pybind11::init( [](const class std::basic_string<char> & a0){ return new HepMC3::WriterAsciiHepMC2(a0); }, [](const class std::basic_string<char> & a0){ return new PyCallBack_HepMC3_WriterAsciiHepMC2(a0); } ), "doc");
		cl.def( pybind11::init<const std::string &, class std::shared_ptr<class HepMC3::GenRunInfo>>(), pybind11::arg("filename"), pybind11::arg("run") );

		cl.def("write_event", (void (HepMC3::WriterAsciiHepMC2::*)(const class HepMC3::GenEvent &)) &HepMC3::WriterAsciiHepMC2::write_event, "Write event to file\n\n \n Event to be serialized\n\nC++: HepMC3::WriterAsciiHepMC2::write_event(const class HepMC3::GenEvent &) --> void", pybind11::arg("evt"));
		cl.def("write_run_info", (void (HepMC3::WriterAsciiHepMC2::*)()) &HepMC3::WriterAsciiHepMC2::write_run_info, "Write the GenRunInfo object to file.\n\nC++: HepMC3::WriterAsciiHepMC2::write_run_info() --> void");
		cl.def("failed", (bool (HepMC3::WriterAsciiHepMC2::*)()) &HepMC3::WriterAsciiHepMC2::failed, "Return status of the stream\n\nC++: HepMC3::WriterAsciiHepMC2::failed() --> bool");
		cl.def("close", (void (HepMC3::WriterAsciiHepMC2::*)()) &HepMC3::WriterAsciiHepMC2::close, "Close file stream\n\nC++: HepMC3::WriterAsciiHepMC2::close() --> void");
		cl.def("set_precision", (void (HepMC3::WriterAsciiHepMC2::*)(const int &)) &HepMC3::WriterAsciiHepMC2::set_precision, "Set output precision\n\n Available range is [2,24]. Default is 16.\n\nC++: HepMC3::WriterAsciiHepMC2::set_precision(const int &) --> void", pybind11::arg("prec"));
		cl.def("precision", (int (HepMC3::WriterAsciiHepMC2::*)() const) &HepMC3::WriterAsciiHepMC2::precision, "Return output precision\n\nC++: HepMC3::WriterAsciiHepMC2::precision() const --> int");
	}
	{ // HepMC3::ReaderAsciiHepMC2 file:HepMC3/ReaderAsciiHepMC2.h line:30
		pybind11::class_<HepMC3::ReaderAsciiHepMC2, std::shared_ptr<HepMC3::ReaderAsciiHepMC2>, PyCallBack_HepMC3_ReaderAsciiHepMC2, HepMC3::Reader> cl(M("HepMC3"), "ReaderAsciiHepMC2", "");
		cl.def( pybind11::init<const std::string &>(), pybind11::arg("filename") );

		cl.def("skip", (bool (HepMC3::ReaderAsciiHepMC2::*)(const int)) &HepMC3::ReaderAsciiHepMC2::skip, "skip events\n\nC++: HepMC3::ReaderAsciiHepMC2::skip(const int) --> bool", pybind11::arg(""));
		cl.def("read_event", (bool (HepMC3::ReaderAsciiHepMC2::*)(class HepMC3::GenEvent &)) &HepMC3::ReaderAsciiHepMC2::read_event, "Implementation of Reader::read_event \n\nC++: HepMC3::ReaderAsciiHepMC2::read_event(class HepMC3::GenEvent &) --> bool", pybind11::arg("evt"));
		cl.def("failed", (bool (HepMC3::ReaderAsciiHepMC2::*)()) &HepMC3::ReaderAsciiHepMC2::failed, "Return status of the stream\n\nC++: HepMC3::ReaderAsciiHepMC2::failed() --> bool");
		cl.def("close", (void (HepMC3::ReaderAsciiHepMC2::*)()) &HepMC3::ReaderAsciiHepMC2::close, "Close file stream\n\nC++: HepMC3::ReaderAsciiHepMC2::close() --> void");
	}
	{ // HepMC3::WriterHEPEVT file:HepMC3/WriterHEPEVT.h line:27
		pybind11::class_<HepMC3::WriterHEPEVT, std::shared_ptr<HepMC3::WriterHEPEVT>, PyCallBack_HepMC3_WriterHEPEVT, HepMC3::Writer> cl(M("HepMC3"), "WriterHEPEVT", "");
		cl.def( pybind11::init<const std::string &>(), pybind11::arg("filename") );

		cl.def("write_hepevt_particle", [](HepMC3::WriterHEPEVT &o, int const & a0) -> void { return o.write_hepevt_particle(a0); }, "", pybind11::arg("index"));
		cl.def("write_hepevt_particle", (void (HepMC3::WriterHEPEVT::*)(int, bool)) &HepMC3::WriterHEPEVT::write_hepevt_particle, "Write particle to file\n\n  \n Particle to be serialized\n  \n\n Format of record\n\nC++: HepMC3::WriterHEPEVT::write_hepevt_particle(int, bool) --> void", pybind11::arg("index"), pybind11::arg("iflong"));
		cl.def("write_hepevt_event_header", (void (HepMC3::WriterHEPEVT::*)()) &HepMC3::WriterHEPEVT::write_hepevt_event_header, "Write event header to file\n\n     \n\nC++: HepMC3::WriterHEPEVT::write_hepevt_event_header() --> void");
		cl.def("write_event", (void (HepMC3::WriterHEPEVT::*)(const class HepMC3::GenEvent &)) &HepMC3::WriterHEPEVT::write_event, "Write event to file\n\n  \n Event to be serialized\n\nC++: HepMC3::WriterHEPEVT::write_event(const class HepMC3::GenEvent &) --> void", pybind11::arg("evt"));
		cl.def("close", (void (HepMC3::WriterHEPEVT::*)()) &HepMC3::WriterHEPEVT::close, "Close file stream \n\nC++: HepMC3::WriterHEPEVT::close() --> void");
		cl.def("failed", (bool (HepMC3::WriterHEPEVT::*)()) &HepMC3::WriterHEPEVT::failed, "Get stream error state flag \n\nC++: HepMC3::WriterHEPEVT::failed() --> bool");
		cl.def("set_vertices_positions_present", (void (HepMC3::WriterHEPEVT::*)(bool)) &HepMC3::WriterHEPEVT::set_vertices_positions_present, "set flag if vertex positions are available. \n  Effectively this adds or removes key \"vertices_positions_are_absent\" \n  to/from the m_options.\n\nC++: HepMC3::WriterHEPEVT::set_vertices_positions_present(bool) --> void", pybind11::arg("iflong"));
		cl.def("get_vertices_positions_present", (bool (HepMC3::WriterHEPEVT::*)() const) &HepMC3::WriterHEPEVT::get_vertices_positions_present, "get flag if vertex positions are available. \n The flag is deduced from m_options. If the m_options have the key \n \"vertices_positions_are_absent\" the result if false. True otherwise. \n\nC++: HepMC3::WriterHEPEVT::get_vertices_positions_present() const --> bool");
	}
	{ // HepMC3::ReaderHEPEVT file:HepMC3/ReaderHEPEVT.h line:32
		pybind11::class_<HepMC3::ReaderHEPEVT, std::shared_ptr<HepMC3::ReaderHEPEVT>, PyCallBack_HepMC3_ReaderHEPEVT, HepMC3::Reader> cl(M("HepMC3"), "ReaderHEPEVT", "");
		cl.def( pybind11::init<const std::string &>(), pybind11::arg("filename") );

		cl.def("read_hepevt_event_header", (bool (HepMC3::ReaderHEPEVT::*)()) &HepMC3::ReaderHEPEVT::read_hepevt_event_header, "Find and read event header line  from file\n\n    \n\nC++: HepMC3::ReaderHEPEVT::read_hepevt_event_header() --> bool");
		cl.def("read_hepevt_particle", (bool (HepMC3::ReaderHEPEVT::*)(int)) &HepMC3::ReaderHEPEVT::read_hepevt_particle, "read particle from file\n\n \n Particle id\n \n\n Event style\n\nC++: HepMC3::ReaderHEPEVT::read_hepevt_particle(int) --> bool", pybind11::arg("i"));
		cl.def("skip", (bool (HepMC3::ReaderHEPEVT::*)(const int)) &HepMC3::ReaderHEPEVT::skip, "skip events\n\nC++: HepMC3::ReaderHEPEVT::skip(const int) --> bool", pybind11::arg(""));
		cl.def("read_event", (bool (HepMC3::ReaderHEPEVT::*)(class HepMC3::GenEvent &)) &HepMC3::ReaderHEPEVT::read_event, "Read event from file\n\nC++: HepMC3::ReaderHEPEVT::read_event(class HepMC3::GenEvent &) --> bool", pybind11::arg("evt"));
		cl.def("close", (void (HepMC3::ReaderHEPEVT::*)()) &HepMC3::ReaderHEPEVT::close, "Close file stream \n\nC++: HepMC3::ReaderHEPEVT::close() --> void");
		cl.def("failed", (bool (HepMC3::ReaderHEPEVT::*)()) &HepMC3::ReaderHEPEVT::failed, "Get stream error state \n\nC++: HepMC3::ReaderHEPEVT::failed() --> bool");
	}
}
