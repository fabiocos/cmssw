#include <iostream>
#include <fstream>
#include <string>
#include <utility>
#include <vector>
#include <algorithm>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESTransientHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/Exception.h"

#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "Geometry/Records/interface/DDSpecParRegistryRcd.h"

#include "DetectorDescription/DDCMS/interface/DDDetector.h"
#include "DetectorDescription/DDCMS/interface/DDSolidShapes.h"
#include "DetectorDescription/DDCMS/interface/DDFilteredView.h"
#include "DetectorDescription/DDCMS/interface/DDSpecParRegistry.h"

#include "Geometry/MTDCommonData/interface/MTDBaseNumber.h"
#include "Geometry/MTDCommonData/interface/BTLNumberingScheme.h"
#include "Geometry/MTDCommonData/interface/ETLNumberingScheme.h"
#include "Geometry/MTDCommonData/interface/BTLElectronicsMapping.h"

#include "DataFormats/ForwardDetId/interface/BTLDetId.h"
#include "DataFormats/ForwardDetId/interface/ETLDetId.h"

#include "DataFormats/Math/interface/angle_units.h"
#include "DataFormats/Math/interface/Rounding.h"
#include <DD4hep/DD4hepUnits.h>

using namespace cms;

class Test_BTLElectronicsMapping : public edm::one::EDAnalyzer<> {
public:
  explicit Test_BTLElectronicsMapping(const edm::ParameterSet&);
  ~Test_BTLElectronicsMapping() override = default;

  void beginJob() override {}
  void analyze(edm::Event const&, edm::EventSetup const&) override;
  void endJob() override {}

  void theBaseNumber(cms::DDFilteredView& fv);

private:
  const edm::ESInputTag tag_;
  std::string ddTopNodeName_;

  MTDBaseNumber thisN_;
  BTLNumberingScheme btlNS_;
  ETLNumberingScheme etlNS_;

  BTLElectronicsMapping elMap_;

  edm::ESGetToken<DDDetector, IdealGeometryRecord> dddetToken_;
  edm::ESGetToken<DDSpecParRegistry, DDSpecParRegistryRcd> dspecToken_;
};

using DD3Vector = ROOT::Math::DisplacementVector3D<ROOT::Math::Cartesian3D<double>>;
using angle_units::operators::convertRadToDeg;
using cms_rounding::roundIfNear0;

Test_BTLElectronicsMapping::Test_BTLElectronicsMapping(const edm::ParameterSet& iConfig)
    : tag_(iConfig.getParameter<edm::ESInputTag>("DDDetector")),
      ddTopNodeName_(iConfig.getUntrackedParameter<std::string>("ddTopNodeName", "BarrelTimingLayer")),
      thisN_(),
      btlNS_(),
      etlNS_(),
      elMap_() {
  dddetToken_ = esConsumes<DDDetector, IdealGeometryRecord>(tag_);
  dspecToken_ = esConsumes<DDSpecParRegistry, DDSpecParRegistryRcd>(tag_);
}

void Test_BTLElectronicsMapping::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  auto pDD = iSetup.getTransientHandle(dddetToken_);

  auto pSP = iSetup.getTransientHandle(dspecToken_);

  if (ddTopNodeName_ != "BarrelTimingLayer" && ddTopNodeName_ != "EndcapTimingLayer") {
    edm::LogWarning("DD4hep_BTLElectronicsMapping") << ddTopNodeName_ << "Not valid top MTD volume";
    return;
  }

  if (!pDD.isValid()) {
    edm::LogError("DD4hep_BTLElectronicsMapping") << "ESTransientHandle<DDCompactView> pDD is not valid!";
    return;
  }
  if (pDD.description()) {
    edm::LogInfo("DD4hep_BTLElectronicsMapping") << pDD.description()->type_ << " label: " << pDD.description()->label_;
  } else {
    edm::LogWarning("DD4hep_BTLElectronicsMapping") << "NO label found pDD.description() returned false.";
  }

  if (!pSP.isValid()) {
    edm::LogError("DD4hep_BTLElectronicsMapping") << "ESTransientHandle<DDSpecParRegistry> pSP is not valid!";
    return;
  }

  DDFilteredView fv(pDD.product(), pDD.product()->description()->worldVolume());
  fv.next(0);
  edm::LogInfo("DD4hep_BTLElectronicsMapping") << fv.name();

  DDSpecParRefs specs;
  std::string attribute("ReadOutName"), name;
  if (ddTopNodeName_ == "BarrelTimingLayer") {
    name = "FastTimerHitsBarrel";
  } else if (ddTopNodeName_ == "EndcapTimingLayer") {
    name = "FastTimerHitsEndcap";
  }
  if (name.empty()) {
    edm::LogError("DD4hep_BTLElectronicsMapping") << "No sensitive detector provided, abort";
    return;
  }
  pSP.product()->filter(specs, attribute, name);

  edm::LogVerbatim("Geometry").log([&specs](auto& log) {
    log << "Filtered DD SpecPar Registry size: " << specs.size() << "\n";
    for (const auto& t : specs) {
      log << "\nSpecPar " << t.first << ":\nRegExps { ";
      for (const auto& ki : t.second->paths)
        log << ki << " ";
      log << "};\n ";
      for (const auto& kl : t.second->spars) {
        log << kl.first << " = ";
        for (const auto& kil : kl.second) {
          log << kil << " ";
        }
        log << "\n ";
      }
    }
  });

  bool write = false;
  bool isBarrel = true;
  bool exitLoop = false;
  uint32_t level(0);
  uint32_t count(0);

  do {
    if (dd4hep::dd::noNamespace(fv.name()) == "BarrelTimingLayer") {
      isBarrel = true;
      edm::LogInfo("DD4hep_BTLElectronicsMapping") << "isBarrel = " << isBarrel;
    } else if (dd4hep::dd::noNamespace(fv.name()) == "EndcapTimingLayer") {
      isBarrel = false;
      edm::LogInfo("DD4hep_BTLElectronicsMapping") << "isBarrel = " << isBarrel;
    }

    if (level > 0 && fv.navPos().size() < level) {
      level = 0;
      write = false;
      if (isBarrel) {
        exitLoop = true;
      } else if (!isBarrel && count == 2) {
        exitLoop = true;
      }
    }
    if (dd4hep::dd::noNamespace(fv.name()) == ddTopNodeName_) {
      write = true;
      level = fv.navPos().size();
      count += 1;
    }

#ifdef EDM_ML_DEBUG
    edm::LogVerbatim("Test_BTLElectronicsMapping")
        << "level= " << level << " isBarrel= " << isBarrel << " exitLoop= " << exitLoop << " count= " << count << " "
        << fv.path();
#endif

    // Test only the desired subdetector

    if (exitLoop && isBarrel) {
      break;
    }

    // Actions for MTD volumes: searchg for sensitive detectors

    if (write) {
      std::stringstream ss;

      theBaseNumber(fv);

      auto print_path = [&]() {
        ss << " - OCMS[0]/";
        for (int ii = thisN_.getLevels() - 1; ii-- > 0;) {
          ss << thisN_.getLevelName(ii);
          ss << "[";
          ss << thisN_.getCopyNumber(ii);
          ss << "]/";
        }
      };

      print_path();

      edm::LogInfo("DD4hep_BTLElectronicsMapping") << ss.str();

      bool isSens = false;

      for (auto const& t : specs) {
        for (auto const& it : t.second->paths) {
          if (dd4hep::dd::compareEqual(dd4hep::dd::noNamespace(fv.name()), dd4hep::dd::realTopName(it))) {
            isSens = true;
            break;
          }
        }
      }

      if (isSens) {
        //
        // Test of numbering scheme for sensitive detectors
        //

        std::stringstream sunitt;
        std::stringstream snum;

        if (!isBarrel) {
          continue;
        }
        BTLDetId theId(btlNS_.getUnitID(thisN_));
        sunitt << theId.rawId();
        snum << theId;
        snum << "\n";
        edm::LogInfo("DD4hep_BTLElectronicsMapping") << snum.str();

        BTLElectronicsMapping::SiPMChPair SiPMChs = elMap_.GetSiPMChPair(theId);
        BTLElectronicsMapping::TOFHIRChPair TOFHIRChs = elMap_.GetTOFHIRChPair(theId);

        //
        // Test of positions for sensitive detectors
        //

        std::stringstream spos;

        auto fround = [&](double in) {
          std::stringstream ss;
          ss << std::fixed << std::setw(14) << roundIfNear0(in);
          return ss.str();
        };

        if (!dd4hep::isA<dd4hep::Box>(fv.solid())) {
          throw cms::Exception("DD4hep_BTLElectronicsMapping") << "MTD sensitive element not a DDBox";
          break;
        }
        dd4hep::Box mySens(fv.solid());

        DD3Vector sidePlusLocal(mySens.x(), 0., 0.);
        DD3Vector sideMinusLocal(-mySens.x(), 0., 0.);
        DD3Vector sidePlusGlobal = (fv.rotation())(sidePlusLocal) + fv.translation();
        DD3Vector sideMinusGlobal = (fv.rotation())(sideMinusLocal) + fv.translation();

        spos << "global z = " << fround(fv.translation().z() / dd4hep::mm) << "\n";

        spos << "Side minus local  = " << fround(sideMinusLocal.X() / dd4hep::mm)
             << fround(sideMinusLocal.Y() / dd4hep::mm) << fround(sideMinusLocal.Z() / dd4hep::mm)
             << " global phi = " << fround(convertRadToDeg(sideMinusGlobal.Phi())) << " SiPMChMinus: " << SiPMChs.Minus
             << " TOFHIRChMinus: " << TOFHIRChs.Minus << "\n";

        spos << "Side plus  local  = " << fround(sidePlusLocal.X() / dd4hep::mm)
             << fround(sidePlusLocal.Y() / dd4hep::mm) << fround(sidePlusLocal.Z() / dd4hep::mm)
             << " global phi = " << fround(convertRadToDeg(sidePlusGlobal.Phi())) << " SiPMChPlus:  " << SiPMChs.Plus
             << " TOFHIRChPlus:  " << TOFHIRChs.Plus << "\n";

        if (SiPMChs.Minus != elMap_.SiPMCh(theId, 0) || SiPMChs.Plus != elMap_.SiPMCh(theId, 1)) {
          spos << "DIFFERENCE IN SiPMChs calculation methods \n";
        }
        if (TOFHIRChs.Minus != elMap_.TOFHIRCh(theId, 0) || TOFHIRChs.Plus != elMap_.TOFHIRCh(theId, 1)) {
          spos << "DIFFERENCE IN TOFHIRChs calculation methods \n";
        }

        // edm::LogInfo("DD4hep_BTLElectronicsMapping") << "Xtal from TOFHIR Channel Minus: " << elMap_.THChToXtal(theId.smodule(),  elMap_.TOFHIRCh(theId, 0))
        //                                         << "\nXtal from TOFHIR Channel Plus : " << elMap_.THChToXtal(theId.smodule(),  elMap_.TOFHIRCh(theId, 1));
        // edm::LogInfo("DD4hep_BTLElectronicsMapping") << "Xtal BTLDetId from TOFHIR Channel Minus: " << elMap_.THChToBTLDetId(theId.zside(), theId.mtdRR(), theId.runit(), theId.dmodule(), theId.smodule(),  elMap_.TOFHIRCh(theId, 0))
        //                                         << "\nXtal BTLDetId from TOFHIR Channel Plus: " << elMap_.THChToBTLDetId(theId.zside(), theId.mtdRR(), theId.runit(), theId.dmodule(), theId.smodule(),  elMap_.TOFHIRCh(theId, 1));
        edm::LogInfo("DD4hep_BTLElectronicsMapping") << spos.str();

        sunitt << fround(fv.translation().z() / dd4hep::mm) << " " << fround(sideMinusLocal.X() / dd4hep::mm) << " "
               << fround(convertRadToDeg(sideMinusGlobal.Phi())) << " " << std::setw(2) << SiPMChs.Minus << " "
               << std::setw(2) << TOFHIRChs.Minus << " " << fround(sidePlusLocal.X() / dd4hep::mm) << " "
               << fround(convertRadToDeg(sidePlusGlobal.Phi())) << " " << std::setw(2) << SiPMChs.Plus << " "
               << std::setw(2) << TOFHIRChs.Plus;

        edm::LogVerbatim("MTDUnitTest") << sunitt.str();
      }
    }
  } while (fv.next(0) && !(exitLoop == 1 && count == 2));
}

void Test_BTLElectronicsMapping::theBaseNumber(cms::DDFilteredView& fv) {
  thisN_.reset();
  thisN_.setSize(fv.navPos().size());

  for (uint ii = 0; ii < fv.navPos().size(); ii++) {
    std::string_view name((fv.geoHistory()[ii])->GetName());
    size_t ipos = name.rfind('_');
    thisN_.addLevel(name.substr(0, ipos), fv.copyNos()[ii]);
#ifdef EDM_ML_DEBUG
    edm::LogVerbatim("DD4hep_BTLElectronicsMapping") << ii << " " << name.substr(0, ipos) << " " << fv.copyNos()[ii];
#endif
  }
}

DEFINE_FWK_MODULE(Test_BTLElectronicsMapping);
