#ifndef DataFormats_VertexReco_MtdtimeSoA_h
#define DataFormats_VertexReco_MtdtimeSoA_h

#include "DataFormats/SoATemplate/interface/SoALayout.h"

GENERATE_SOA_LAYOUT(MtdtimeSoALayout,
                    SOA_COLUMN(int32_t, trackAsocMTD),
		    SOA_SCALAR(uint32_t, IndxTrackAsocMTD),
                    SOA_COLUMN(float, time),
                    SOA_COLUMN(float, timeErr),
                    SOA_COLUMN(float, MVAquality),
                    SOA_COLUMN(float, pathLength),
                    SOA_COLUMN(float, beta),
                    SOA_COLUMN(float, posInMTD_x),
                    SOA_COLUMN(float, posInMTD_y),
                    SOA_COLUMN(float, posInMTD_z),
                    SOA_COLUMN(float, momentumWithMTD),
		    SOA_COLUMN(float, btlMatch_chi2),
		    SOA_COLUMN(float, btlMatchTime_chi2),
		    SOA_COLUMN(float, etlMatch_chi2),
                    SOA_COLUMN(float, etlMatchTime_chi2),
		    SOA_COLUMN(float, trackTime_pi),
                    SOA_COLUMN(float, trackTime_k),
		    SOA_COLUMN(float, trackTime_p),
		    SOA_COLUMN(float, track_sigmaTime_pi),
		    SOA_COLUMN(float, track_sigmaTime_k),
		    SOA_COLUMN(float, track_sigmaTime_p))

using MtdtimeSoA = MtdtimeSoALayout<>;
using MtdtimeSoAView = MtdtimeSoA::View;
using MtdtimeSoAConstView = MtdtimeSoA::ConstView;

#endif
