// -*- C++ -*-
//
// This file is part of YODA -- Yet more Objects for Data Analysis
// Copyright (C) 2008-2016 The YODA collaboration (see AUTHORS for details)
//
#include "YODA/Exceptions.h"

// !!! DO NOT MOVE THESE DEFINITIONS INTO THE HEADER !!!!

// Exceptions need to be defined in _one_ specific library.
//
// If the defn.s live in the header, any library including them 
// will have its own copy. These copies are _not_ mutually 
// interchangeable, and throw / catch across a library 
// boundary will not work.

// !!! DO NOT MOVE THESE DEFINITIONS INTO THE HEADER !!!!

YODA::Exception::Exception(const std::string& what) 
	: std::runtime_error(what) {}

YODA::BinningError::BinningError(const std::string& what) 
	: YODA::Exception(what) {}

YODA::RangeError::RangeError(const std::string& what) 
	: YODA::Exception(what) {}

YODA::LockError::LockError(const std::string& what) 
	: YODA::Exception(what) {}

YODA::GridError::GridError(const std::string& what) 
	: YODA::Exception(what) {}

YODA::LogicError::LogicError(const std::string& what) 
	: YODA::Exception(what) {}

YODA::WeightError::WeightError(const std::string& what) 
	: YODA::Exception(what) {}

YODA::LowStatsError::LowStatsError(const std::string& what) 
	: YODA::Exception(what) {}

YODA::AnnotationError::AnnotationError(const std::string& what) 
	: YODA::Exception(what) {}

YODA::ReadError::ReadError(const std::string& what) 
	: YODA::Exception(what) {}

YODA::WriteError::WriteError(const std::string& what) 
	: YODA::Exception(what) {}

YODA::UserError::UserError(const std::string& what) 
	: YODA::Exception(what) {}
