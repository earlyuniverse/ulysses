#include <HepMC3/Attribute.h>
#include <HepMC3/GenCrossSection.h>
#include <HepMC3/GenHeavyIon.h>
#include <HepMC3/GenParticle.h>
#include <HepMC3/GenVertex.h>
#include <HepMC3/LHEF.h>
#include <functional>
#include <iterator>
#include <map>
#include <memory>
#include <ostream>
#include <set>
#include <sstream> // __str__
#include <string>
#include <utility>

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

void bind_pyHepMC3_2(std::function< pybind11::module &(std::string const &namespace_) > &M)
{
	// std::map file:bits/stl_map.h line:96
	binder::map_binder<std::string,std::shared_ptr<HepMC3::Attribute>,std::less<std::string >,std::allocator<std::pair<const std::string, std::shared_ptr<HepMC3::Attribute> > >>(M("std"), "std_string", "std_shared_ptr_HepMC3_Attribute_t", "std_less_std_string_t", "std_allocator_std_pair_const_std_string_std_shared_ptr_HepMC3_Attribute_t");

	// std::map file:bits/stl_map.h line:96
	binder::map_binder<std::string,std::map<int, std::shared_ptr<HepMC3::Attribute>, std::less<int>, std::allocator<std::pair<const int, std::shared_ptr<HepMC3::Attribute> > > >,std::less<std::string >,std::allocator<std::pair<const std::string, std::map<int, std::shared_ptr<HepMC3::Attribute>, std::less<int>, std::allocator<std::pair<const int, std::shared_ptr<HepMC3::Attribute> > > > > >>(M("std"), "std_string", "std_map_int_std_shared_ptr_HepMC3_Attribute_std_less_int_std_allocator_std_pair_const_int_std_shared_ptr_HepMC3_Attribute_t", "std_less_std_string_t", "std_allocator_std_pair_const_std_string_std_map_int_std_shared_ptr_HepMC3_Attribute_std_less_int_std_allocator_std_pair_const_int_std_shared_ptr_HepMC3_Attribute_t");

	// std::map file:bits/stl_map.h line:96
	binder::map_binder<std::string,std::string,std::less<std::string >,std::allocator<std::pair<const std::string, std::string > >>(M("std"), "std_string", "std_string", "std_less_std_string_t", "std_allocator_std_pair_const_std_string_std_string_t");

	// std::map file:bits/stl_map.h line:96
	binder::map_binder<std::string,std::set<long, std::less<long>, std::allocator<long> >,std::less<std::string >,std::allocator<std::pair<const std::string, std::set<long, std::less<long>, std::allocator<long> > > >>(M("std"), "std_string", "std_set_long_std_less_long_std_allocator_long_t", "std_less_std_string_t", "std_allocator_std_pair_const_std_string_std_set_long_std_less_long_std_allocator_long_t");

}
