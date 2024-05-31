#ifndef DataFormats_VertexReco_MtdtimeHostCollection_h
#define DataFormats_VertexReco_MtdtimeHostCollection_h

#include "DataFormats/Portable/interface/PortableHostCollection.h"
#include "DataFormats/VertexReco/interface/MtdtimeSoA.h"

// MtdSoA in host memory
using MtdtimeHostCollection = PortableHostCollection<MtdtimeSoA>;
using MtdtimeHostCollectionView = PortableHostCollection<MtdtimeSoA>::View;
using MtdtimeHostCollectionConstView = PortableHostCollection<MtdtimeSoA>::ConstView;

#endif
