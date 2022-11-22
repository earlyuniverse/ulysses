// -*- C++ -*-
//
// This file is part of YODA -- Yet more Objects for Data Analysis
// Copyright (C) 2008-2017 The YODA collaboration (see AUTHORS for details)
//
#ifndef YODA_FORMATTING_H
#define YODA_FORMATTING_H

#include <iostream>
#include <iomanip>
#include <unistd.h>

#define MSG_(msg) do { std::cout << msg; } while (0)

#define MSG(msg) MSG_(msg << std::endl)

#define PAD(n) std::setw(n) << std::left

#define COLOR_(msg, code) \
  (isatty(1) ? code : "") << msg << (isatty(1) ? "\033[0m" : "")


#define RED(msg)        COLOR_(msg, "\033[0;31m")
#define MSG_RED_(x)     MSG_(RED(x))
#define MSG_RED(x)      MSG(RED(x))

#define GREEN(msg)      COLOR_(msg, "\033[0;32m")
#define MSG_GREEN_(x)   MSG_(GREEN(x))
#define MSG_GREEN(x)    MSG(GREEN(x))

#define YELLOW(msg)     COLOR_(msg, "\033[0;33m")
#define MSG_YELLOW_(x)  MSG_(YELLOW(x))
#define MSG_YELLOW(x)   MSG(YELLOW(x))

#define BLUE(msg)       COLOR_(msg, "\033[0;34m")
#define MSG_BLUE_(x)    MSG_(BLUE(x))
#define MSG_BLUE(x)     MSG(BLUE(x))


#endif
