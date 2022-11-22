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
#include <HepMC3/Reader.h>
#include <HepMC3/ReaderLHEF.h>
#include <HepMC3/ReaderPlugin.h>
#include <HepMC3/WriterPlugin.h>
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

// HepMC3::ReaderLHEF file:HepMC3/ReaderLHEF.h line:35
struct PyCallBack_HepMC3_ReaderLHEF : public HepMC3::ReaderLHEF {
	using HepMC3::ReaderLHEF::ReaderLHEF;

	bool skip(const int a0) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::ReaderLHEF *>(this), "skip");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return ReaderLHEF::skip(a0);
	}
	bool read_event(class HepMC3::GenEvent & a0) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::ReaderLHEF *>(this), "read_event");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return ReaderLHEF::read_event(a0);
	}
	void close() override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::ReaderLHEF *>(this), "close");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return ReaderLHEF::close();
	}
	bool failed() override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::ReaderLHEF *>(this), "failed");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return ReaderLHEF::failed();
	}
};

// HepMC3::ReaderPlugin file:HepMC3/ReaderPlugin.h line:23
struct PyCallBack_HepMC3_ReaderPlugin : public HepMC3::ReaderPlugin {
	using HepMC3::ReaderPlugin::ReaderPlugin;

	bool read_event(class HepMC3::GenEvent & a0) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::ReaderPlugin *>(this), "read_event");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return ReaderPlugin::read_event(a0);
	}
	void close() override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::ReaderPlugin *>(this), "close");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return ReaderPlugin::close();
	}
	bool failed() override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::ReaderPlugin *>(this), "failed");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return ReaderPlugin::failed();
	}
	bool skip(const int a0) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::ReaderPlugin *>(this), "skip");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return Reader::skip(a0);
	}
};

// HepMC3::WriterPlugin file:HepMC3/WriterPlugin.h line:23
struct PyCallBack_HepMC3_WriterPlugin : public HepMC3::WriterPlugin {
	using HepMC3::WriterPlugin::WriterPlugin;

	void write_event(const class HepMC3::GenEvent & a0) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::WriterPlugin *>(this), "write_event");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return WriterPlugin::write_event(a0);
	}
	void close() override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::WriterPlugin *>(this), "close");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return WriterPlugin::close();
	}
	bool failed() override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const HepMC3::WriterPlugin *>(this), "failed");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return WriterPlugin::failed();
	}
};

void bind_pyHepMC3_18(std::function< pybind11::module &(std::string const &namespace_) > &M)
{
	{ // HepMC3::ReaderLHEF file:HepMC3/ReaderLHEF.h line:35
		pybind11::class_<HepMC3::ReaderLHEF, std::shared_ptr<HepMC3::ReaderLHEF>, PyCallBack_HepMC3_ReaderLHEF, HepMC3::Reader> cl(M("HepMC3"), "ReaderLHEF", "");
		cl.def( pybind11::init<const std::string &>(), pybind11::arg("filename") );

		cl.def("skip", (bool (HepMC3::ReaderLHEF::*)(const int)) &HepMC3::ReaderLHEF::skip, "skip events\n\nC++: HepMC3::ReaderLHEF::skip(const int) --> bool", pybind11::arg(""));
		cl.def("read_event", (bool (HepMC3::ReaderLHEF::*)(class HepMC3::GenEvent &)) &HepMC3::ReaderLHEF::read_event, "Reading event \n\nC++: HepMC3::ReaderLHEF::read_event(class HepMC3::GenEvent &) --> bool", pybind11::arg("ev"));
		cl.def("close", (void (HepMC3::ReaderLHEF::*)()) &HepMC3::ReaderLHEF::close, "Close \n\nC++: HepMC3::ReaderLHEF::close() --> void");
		cl.def("failed", (bool (HepMC3::ReaderLHEF::*)()) &HepMC3::ReaderLHEF::failed, "State \n\nC++: HepMC3::ReaderLHEF::failed() --> bool");
	}
	{ // HepMC3::ReaderPlugin file:HepMC3/ReaderPlugin.h line:23
		pybind11::class_<HepMC3::ReaderPlugin, std::shared_ptr<HepMC3::ReaderPlugin>, PyCallBack_HepMC3_ReaderPlugin, HepMC3::Reader> cl(M("HepMC3"), "ReaderPlugin", "");
		cl.def( pybind11::init<const std::string &, const std::string &, const std::string &>(), pybind11::arg("filename"), pybind11::arg("libname"), pybind11::arg("newreader") );

		cl.def("read_event", (bool (HepMC3::ReaderPlugin::*)(class HepMC3::GenEvent &)) &HepMC3::ReaderPlugin::read_event, "Reading event \n\nC++: HepMC3::ReaderPlugin::read_event(class HepMC3::GenEvent &) --> bool", pybind11::arg("ev"));
		cl.def("close", (void (HepMC3::ReaderPlugin::*)()) &HepMC3::ReaderPlugin::close, "Close \n\nC++: HepMC3::ReaderPlugin::close() --> void");
		cl.def("failed", (bool (HepMC3::ReaderPlugin::*)()) &HepMC3::ReaderPlugin::failed, "State \n\nC++: HepMC3::ReaderPlugin::failed() --> bool");
	}
	{ // HepMC3::WriterPlugin file:HepMC3/WriterPlugin.h line:23
		pybind11::class_<HepMC3::WriterPlugin, std::shared_ptr<HepMC3::WriterPlugin>, PyCallBack_HepMC3_WriterPlugin, HepMC3::Writer> cl(M("HepMC3"), "WriterPlugin", "");
		cl.def( pybind11::init( [](const class std::basic_string<char> & a0, const class std::basic_string<char> & a1, const class std::basic_string<char> & a2){ return new HepMC3::WriterPlugin(a0, a1, a2); }, [](const class std::basic_string<char> & a0, const class std::basic_string<char> & a1, const class std::basic_string<char> & a2){ return new PyCallBack_HepMC3_WriterPlugin(a0, a1, a2); } ), "doc");
		cl.def( pybind11::init<const std::string &, const std::string &, const std::string &, class std::shared_ptr<class HepMC3::GenRunInfo>>(), pybind11::arg("filename"), pybind11::arg("libname"), pybind11::arg("newwriter"), pybind11::arg("run") );

		cl.def("write_event", (void (HepMC3::WriterPlugin::*)(const class HepMC3::GenEvent &)) &HepMC3::WriterPlugin::write_event, "Reading event \n\nC++: HepMC3::WriterPlugin::write_event(const class HepMC3::GenEvent &) --> void", pybind11::arg("ev"));
		cl.def("close", (void (HepMC3::WriterPlugin::*)()) &HepMC3::WriterPlugin::close, "Close \n\nC++: HepMC3::WriterPlugin::close() --> void");
		cl.def("failed", (bool (HepMC3::WriterPlugin::*)()) &HepMC3::WriterPlugin::failed, "State \n\nC++: HepMC3::WriterPlugin::failed() --> bool");
	}
}
