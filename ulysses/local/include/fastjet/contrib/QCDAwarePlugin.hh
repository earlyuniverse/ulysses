//  QCDAware Package
//  Questions/Comments?  abuckley@cern.ch, cpollard@cern.ch
//
//  Copyright (c) 2014
//  Andy Buckley, Chris Pollard
//
// $Id: QCDAwarePlugin.hh 887 2015-10-08 08:27:29Z cspollard $
//
//----------------------------------------------------------------------
// This file is part of FastJet contrib.
//
// It is free software; you can redistribute it and/or modify it under
// the terms of the GNU General Public License as published by the
// Free Software Foundation; either version 2 of the License, or (at
// your option) any later version.
//
// It is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
// or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public
// License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this code. If not, see <http://www.gnu.org/licenses/>.
//----------------------------------------------------------------------

#ifndef __FASTJET_CONTRIB_QCDAWAREPLUGIN_HH__
#define __FASTJET_CONTRIB_QCDAWAREPLUGIN_HH__

#include "fastjet/internal/base.hh"
#include "fastjet/JetDefinition.hh"
#include "fastjet/ClusterSequence.hh"
#include <queue>
#include <string>
#include <vector>
#include "DistanceMeasure.hh"

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh


namespace contrib {
    namespace QCDAwarePlugin {

        /**
         * \brief the PJDist struct: simple storage of distance
         *        information between two pseudojets
         */

        struct PJDist {
            double dist;
            int pj1;
            int pj2;
        };


        /**
         * \brief the QCDAwarePlugin class performs QCD-aware jet
         *        clustering given a particular DistanceMeasure
         *
         *  the QCDAwarePlugin class performs the QCD-aware jet
         *  clustering given a particular DistanceMeasure.
         */
        //------------------------------------------------------------------------
        // @todo
        // this class isn't really QCDAware; the distance measure is...
        class QCDAwarePlugin : public JetDefinition::Plugin {

            public:
                // User still owns the pointer to dm after using this
                // constructor.
                QCDAwarePlugin(const DistanceMeasure *dm)
                    : _dm(dm) {}

                /// default destructor
                // we don't delete _dm here because it is owned by the
                // user.
                virtual ~QCDAwarePlugin() {}

                void run_clustering(fastjet::ClusterSequence& cs) const;

                std::string description() const;

                double R() const;


            private:
                const DistanceMeasure *_dm;

                void insert_pj(ClusterSequence &cs,
                        std::priority_queue<PJDist, std::vector<PJDist>, std::greater<PJDist> >& pjds,
                        unsigned int iJet,
                        std::vector<bool>& ismerged) const;

                void merge_iB(ClusterSequence &cs,
                        const PJDist& dist,
                        std::vector<bool>& ismerged) const;

                void merge_ij(ClusterSequence &cs,
                        std::priority_queue<PJDist, std::vector<PJDist>, std::greater<PJDist> >& pjds,
                        const PJDist& dist,
                        std::vector<bool>& ismerged) const;

                // returns zero if p and q aren't allowed to combine.
                // returns the combined pid otherwise.
                int flavor_sum(const fastjet::PseudoJet& p, const fastjet::PseudoJet& q) const;
        };


    } // QCDAware
} // contrib

FASTJET_END_NAMESPACE

#endif  // __FASTJET_CONTRIB_QCDAWARE_HH__
