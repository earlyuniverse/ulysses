#include "WriterHEPEVTZEUS.h"
#include "HepMC3/HEPEVT_Wrapper.h"
namespace HepMC3
{
WriterHEPEVTZEUS::WriterHEPEVTZEUS(const std::string &filename):WriterHEPEVT(filename) {}
void WriterHEPEVTZEUS::write_hepevt_event_header()
{
    char buf[512];//Note: the format is fixed, so no reason for complicatied tratment
    char* cursor = &(buf[0]);
    cursor += sprintf(cursor, " E % 12i% 12i% 12i\n", HEPEVT_Wrapper::event_number(), 0, HEPEVT_Wrapper::number_entries());
    unsigned long length = cursor - &(buf[0]);
    m_stream->write( buf, length );
}
void WriterHEPEVTZEUS::write_hepevt_particle( int index)
{
    char buf[512];//Note: the format is fixed, so no reason for complicatied tratment
    char* cursor = &(buf[0]);
    cursor += sprintf(cursor,"% 12i% 8i", HEPEVT_Wrapper::status(index), HEPEVT_Wrapper::id(index));
    cursor += sprintf(cursor,"% 8i% 8i", HEPEVT_Wrapper::first_parent(index), HEPEVT_Wrapper::last_parent(index));
    cursor += sprintf(cursor,"% 8i% 8i", HEPEVT_Wrapper::first_child(index), HEPEVT_Wrapper::last_child(index));
    cursor += sprintf(cursor,      "% 19.11E% 19.11E% 19.11E% 19.11E% 19.11E\n", HEPEVT_Wrapper::px(index), HEPEVT_Wrapper::py(index), HEPEVT_Wrapper::pz(index), HEPEVT_Wrapper::e(index), HEPEVT_Wrapper::m(index));
    cursor += sprintf(cursor, "%-52s% 19.11E% 19.11E% 19.11E% 19.11E% 19.11E\n", " ", HEPEVT_Wrapper::x(index), HEPEVT_Wrapper::y(index), HEPEVT_Wrapper::z(index), HEPEVT_Wrapper::t(index), 0.0);
    std::ptrdiff_t length = cursor - &(buf[0]);
    m_stream->write( buf, length );
}
}// namespace HepMC3
