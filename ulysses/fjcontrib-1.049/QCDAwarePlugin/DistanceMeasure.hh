#ifndef __DISTANCEMEASURE_HH__
#define __DISTANCEMEASURE_HH__

#include <float.h>
#include <algorithm>
#include "fastjet/PseudoJet.hh"

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

namespace contrib {
    namespace QCDAwarePlugin {

        /**
         * the DistanceMeasure abstract base class calculates the
         * distance between two PseudoJets. the default implemented
         * DistanceMeasures include kt, C/A, and anti-kt
         */

        class DistanceMeasure {
            public:
                virtual double dij(const fastjet::PseudoJet& pji, const fastjet::PseudoJet& pjj) const = 0;
                virtual double diB(const fastjet::PseudoJet& pji) const = 0;
                virtual double R() const = 0;
                virtual std::string algname() const = 0;

                virtual ~DistanceMeasure() {};
        };


        class KtMeasure : public DistanceMeasure {
            private:
                double _R;

            public:
                KtMeasure(double R) : _R(R) {}

                inline double dij(const fastjet::PseudoJet& pji, const fastjet::PseudoJet& pjj) const {
                    double drbyR2 = pji.squared_distance(pjj) / (_R * _R);
                    return std::min(pji.perp2(), pjj.perp2()) * drbyR2;
                }

                inline double diB(const fastjet::PseudoJet& pji) const {
                    return pji.perp2();
                }

                virtual double R() const {
                    return _R;
                }

                std::string algname() const {
                    return "kt";
                }
        };


        class AntiKtMeasure : public DistanceMeasure {
            private:
                double _R;

            public:
                AntiKtMeasure(double R) : _R(R) {}

                inline double dij(const fastjet::PseudoJet& pji, const fastjet::PseudoJet& pjj) const {
                    double drbyR2 = pji.squared_distance(pjj) / (_R * _R);
                    return 1.0 / std::max(pji.perp2(), pjj.perp2()) * drbyR2;
                }

                inline double diB(const fastjet::PseudoJet& pji) const {
                    return 1.0 / pji.perp2();
                }

                virtual double R() const {
                    return _R;
                }

                std::string algname() const {
                    return "anti-kt";
                }
        };


        class CAMeasure : public DistanceMeasure {
            private:
                double _R;

            public:
                CAMeasure(double R) : _R(R) {}

                inline double dij(const fastjet::PseudoJet& pji, const fastjet::PseudoJet& pjj) const {
                    return pji.squared_distance(pjj) / (_R * _R);
                }


                inline double diB(const fastjet::PseudoJet& pji) const {
                    return 1.0;
                }

                virtual double R() const {
                    return _R;
                }

                std::string algname() const {
                    return "Cambridge-Aachen";
                }
        };

    } // QCDAware
} // contrib

FASTJET_END_NAMESPACE

#endif
