// -*- C++ -*-
//
// This file is part of HepMC
// Copyright (C) 2014-2019 The HepMC collaboration (see AUTHORS for details)
//
#ifndef HEPMC3_READERGZ_H
#define HEPMC3_READERGZ_H
///
/// @file  ReaderGZ.h
/// @brief Definition of class \b ReaderGZ
///
/// @class HepMC3::ReaderGZ
/// @brief GenEvent I/O parsing for structured text files compressed with gzip
///
/// @ingroup IO
///
#include <string>
#include <fstream>
#include <istream>
#include <string.h>
#include "HepMC3/Reader.h"
#include "HepMC3/ReaderAscii.h"
#include "HepMC3/ReaderAsciiHepMC2.h"
#include "HepMC3/ReaderHEPEVT.h"
#include "HepMC3/ReaderLHEF.h"
#include "HepMC3/GenEvent.h"
#include "gzstream.h"
namespace HepMC3 {
/** @brief Union to hold first 4 byts of file, i.e. magic bytes */
union magic_t {
    uint8_t bytes[4]; ///< bytes
    uint32_t number;  ///< int
};
class ReaderGZ : public Reader {
public:
    /** @brief Construcor*/
    ReaderGZ(const std::string& filename) : m_gzstream(filename.c_str())
    {
        std::ifstream file(filename);
        if(!file.is_open()) {
            printf("Error in ReaderGZ: could not open file%s\n",filename.c_str());
            return;
        }
        magic_t my_magic = {0x1f, 0x8b, 0x08, 0x08};
        magic_t file_magic;
        file.read((char *) file_magic.bytes, sizeof(file_magic));
        if ( file_magic.number == my_magic.number )
        {
            m_reader=deduce_reader(m_gzstream);
        }
        else
        {
            printf("Error in ReaderGZ: make sure %s is a gziped file!\n",filename.c_str());
            return;
        }
    };

    ~ReaderGZ() {};
    /** @brief Read event */
    bool read_event(GenEvent& evt) {
        return m_reader->read_event(evt);
    };
    /** @brief State */
    bool failed() {
        return m_gzstream.rdstate();
    }
    /** @brief Close  */
    void close()  {
        if (m_reader) m_reader->close();
    };
private:
    igzstream   m_gzstream;  ///< Stream to read
    std::shared_ptr<Reader>      m_reader; ///< Actual reader
};
}
#endif
