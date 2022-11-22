#include "YODA/Histo1D.h"
#include "YODA/WriterYODA.h"
#include <vector>
#include <memory>
#include <type_traits>

using namespace std;
using namespace YODA;

// /// SFINAE struct to check for dereferencing to AnalysisObject at compile time
// template <typename T, typename=void>
// struct XDerefableToAO : std::false_type {};
// //
// template <typename T>
// struct XDerefableToAO<T, typename std::decay< decltype(*std::declval<T>()) >::type> : std::true_type {};

static_assert(std::is_same<AnalysisObject, typename std::decay<const AnalysisObject&>::type>::value, "const AO& doesn't decay to AO?!");
static_assert(std::is_same<AnalysisObject, typename std::decay<decltype( *std::declval<AnalysisObject*>() )>::type>::value, "*(AO*) doesn't decay to AO via decltype?!");
static_assert(std::is_same<AnalysisObject, typename std::decay<decltype( *std::declval<const AnalysisObject*>() )>::type>::value, "*(const AO*) doesn't decay to AO via decltype?!");
static_assert(std::is_same<AnalysisObject, typename std::decay<decltype( *std::declval<shared_ptr<AnalysisObject>>() )>::type>::value, "shared_ptr<AO> doesn't decay to AO via decltype?!");
static_assert(std::is_same<AnalysisObject, typename std::decay<decltype( *std::declval<shared_ptr<const AnalysisObject>>() )>::type>::value, "shared_ptr<const AO> doesn't decay to AO via decltype?!");
static_assert(std::is_same<AnalysisObject, typename std::decay<decltype( *std::declval<const shared_ptr<AnalysisObject>>() )>::type>::value, "const shared_ptr<AO> doesn't decay to AO via decltype?!");
static_assert(std::is_same<AnalysisObject, typename std::decay<decltype( *std::declval<const shared_ptr<const AnalysisObject>>() )>::type>::value, "const shared_ptr<const AO> doesn't decay to AO via decltype?!");


class A {};
class B : A {};
static_assert(std::is_base_of<A, B>::value, "FOO");
static_assert(std::is_base_of<A, const B>::value, "FOO");
static_assert(std::is_base_of<A, std::remove_reference<const B&>::type>::value, "FOO");
static_assert(std::is_base_of<A, std::decay<const B&>::type>::value, "FOO");


static_assert(Iterable< vector<AnalysisObject*> >::value, "Not iterable");

static_assert(Derefable< AnalysisObject* >::value, "AO* not derefable");
static_assert(Derefable< const AnalysisObject* >::value, "const AO* not derefable");
static_assert(Derefable< shared_ptr<AnalysisObject> >::value, "shared_ptr<AO> not derefable");
static_assert(Derefable< shared_ptr<const AnalysisObject> >::value, "shared_ptr<const AO> not derefable");
static_assert(Derefable< shared_ptr<Histo1D> >::value, "shared_ptr<H1> not derefable");
static_assert(Derefable< shared_ptr<const Histo1D> >::value, "shared_ptr<const H1> not derefable");


static_assert(std::is_base_of<AnalysisObject, AnalysisObject>::value, "AO doesn't inherit from AO?!");
static_assert(std::is_base_of<AnalysisObject, const AnalysisObject>::value, "const AO doesn't inherit from AO?!");
static_assert(std::is_base_of<AnalysisObject, typename std::decay< decltype(*std::declval<AnalysisObject*>()) >::type>::value, "decay<*(H1*)> doesn't inherit from AO?!");
static_assert(std::is_base_of<AnalysisObject, Histo1D>::value, "Histo1D doesn't inherit from AO?!");
static_assert(std::is_base_of<AnalysisObject, const Histo1D>::value, "const Histo1D doesn't inherit from AO?!");
static_assert(std::is_base_of<AnalysisObject, typename std::decay< decltype(*std::declval<Histo1D*>()) >::type>::value, "decay<*(H1*)> doesn't inherit from AO?!");
//static_assert(std::is_base_of<AnalysisObject, decltype(*std::declval<Histo1D*>())>::value, "*(H1*) doesn't inherit from AO?!"); // deref gives a reference: decay needed

static_assert(std::is_same<AnalysisObject,
              std::conditional<std::is_base_of<AnalysisObject, typename std::decay< decltype(*std::declval<AnalysisObject*>()) >::type>::value,
              AnalysisObject, void>::type
              >::value, "FOOOO");

// /// SFINAE struct to check for dereferencing to AnalysisObject at compile time
// template <typename T, typename=YODA::AnalysisObject>
// struct XDerefableToAO : std::false_type {};
// //
// template <typename T>
// // struct XDerefableToAO<T, typename std::decay< decltype(*std::declval<T>()) >::type> : std::true_type {};
// struct XDerefableToAO<T, typename std::conditional<std::is_base_of<AnalysisObject, typename std::decay< decltype(*std::declval<T>()) >::type>::value,
//                                                    YODA::AnalysisObject, void>::type> : std::true_type {};

static_assert(DerefableToAO< AnalysisObject* >::value, "AO* not derefable to AO");
static_assert(DerefableToAO< const AnalysisObject* >::value, "const AO* not derefable to AO");
static_assert(DerefableToAO< shared_ptr<AnalysisObject> >::value, "shared_ptr<AO> not derefable to AO");
static_assert(DerefableToAO< shared_ptr<const AnalysisObject> >::value, "shared_ptr<const AO> not derefable to AO");
static_assert(DerefableToAO< Histo1D* >::value, "shared_ptr<H1> not derefable to AO");
static_assert(DerefableToAO< const Histo1D* >::value, "shared_ptr<const H1> not derefable to AO");
static_assert(DerefableToAO< shared_ptr<Histo1D> >::value, "shared_ptr<H1> not derefable to AO");
static_assert(DerefableToAO< shared_ptr<const Histo1D> >::value, "shared_ptr<const H1> not derefable to AO");


int main() {
  return EXIT_SUCCESS;
}
