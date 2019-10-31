#include "CombineHarvester/TauIDSFMeasurement/interface/HttSystematics_TauIDRun2.h"
#include "CombineHarvester/CombineTools/interface/Process.h"
#include "CombineHarvester/CombineTools/interface/Systematics.h"
#include "CombineHarvester/CombineTools/interface/Utilities.h"
#include <string>
#include <vector>

using namespace std;

namespace ch {

using ch::syst::SystMap;
using ch::syst::SystMapAsymm;
using ch::syst::era;
using ch::syst::channel;
using ch::syst::bin_id;
using ch::syst::process;
using ch::syst::bin;
using ch::JoinStr;

void AddTauIDRun2Systematics(CombineHarvester &cb, bool jetfakes, bool embedding, bool regional_jec, int era) {

  // ##########################################################################
  // Define groups of processes
  // ##########################################################################

  std::vector<std::string> mc_processes;
  // Signal processes
  if (embedding){
      std::vector<std::string> signals = {"EMB"};
      mc_processes =
                  {"TT", "TTL", "TTJ", "W", "ZJ", "ZL", "VV", "VVL", "VVJ", "ST"};
  }
  else {
      std::vector<std::string> signals = {"ZTT"};
      mc_processes =
                  {"ZTT", "TT", "TTT", "TTL", "TTJ", "W", "ZJ", "ZL", "VV", "VVT", "VVL", "VVJ", "ST"};
  }

  // ##########################################################################
  // Uncertainty: Lumi
  // References:
  // - "CMS Luminosity Measurements for the 2016 Data Taking Period"
  //   (PAS, https://cds.cern.ch/record/2257069)
  // Notes:
  // - FIXME: Adapt for fake factor and embedding
  // ##########################################################################

  float lumi_unc = 1.0;
  if (era == 2016) {
      lumi_unc = 1.025;
  } else if (era == 2017) { // TODO: Check if number for 2017 is still valid.
      lumi_unc = 1.023;
  } else if (era == 2018) {
      lumi_unc = 1.025;
  }

  cb.cp()
      .channel({"mt", "mm"})
      .process(mc_processes)
      .AddSyst(cb, "lumi_$ERA", "lnN", SystMap<>::init(lumi_unc));

  // ##########################################################################
  // Uncertainty: Trigger efficiency
  // References:
  // Notes:
  // - FIXME: References?
  // ##########################################################################

  cb.cp()
      .channel({"mt"})
      .process(mc_processes)
      .AddSyst(cb, "CMS_eff_trigger_mt_$ERA", "shape", SystMap<>::init(1.00));

  // 100% uncorrelated for embedded
  cb.cp()
      .channel({"mt"})
      .process({"EMB"})
      .AddSyst(cb, "CMS_eff_trigger_emb_mt_$ERA", "shape", SystMap<>::init(1.00));

  // Muon ID
  cb.cp()
      .channel({"mt"})
      .process(mc_processes)
      .AddSyst(cb, "CMS_eff_mc_m", "lnN", SystMap<>::init(1.014));

  // Muon ID
  cb.cp()
      .channel({"mt"})
      .process({"EMB"})
      .AddSyst(cb, "CMS_eff_emb_m", "lnN", SystMap<>::init(1.014));

  // MC + embedded correlated uncertainty

  // Muon ID
  cb.cp()
      .channel({"mt"})
      .process(JoinStr({mc_processes, {"EMB"}}))
      .AddSyst(cb, "CMS_eff_m", "lnN", SystMap<>::init(1.014));

  // Control region uncertainties
  // TODO: Check how to correlate with embedded if embedded is used
  cb.cp()
      .channel({"mm"})
      .process(mc_processes)
      .AddSyst(cb, "CMS_eff_mc_m", "lnN", SystMap<>::init(1.04));

  // ##########################################################################
  // Uncertainty: b-tag and mistag efficiency
  // References:
  // Notes:
  // - FIXME: References?
  // ##########################################################################

  // cb.cp()
  //     .channel({"mt"})
  //     .process(mc_processes)
  //     .AddSyst(cb, "CMS_htt_eff_b_$ERA", "shape", SystMap<>::init(1.00));

  // cb.cp()
  //     .channel({"mt"})
  //     .process(mc_processes)
  //     .AddSyst(cb, "CMS_htt_mistag_b_$ERA", "shape", SystMap<>::init(1.00));

  // ##########################################################################
  // Uncertainty: Tau energy scale
  // References:
  // Notes:
  // - Tau energy scale is splitted by decay mode.
  // - FIXME: References?
  // ##########################################################################

  //// correlated part between ERAs
  // MC uncorrelated uncertainty

  cb.cp()
      .channel({"mt"})
      .process({{"ZTT", "TTT", "TTL", "VVT", "VVL", "jetFakes"}})
      .AddSyst(cb, "CMS_scale_mc_t_1prong", "shape", SystMap<>::init(0.71));

  cb.cp()
      .channel({"mt"})
      .process({"ZTT", "TTT", "TTL", "VVT", "VVL", "jetFakes"})
      .AddSyst(cb, "CMS_scale_mc_t_1prong1pizero", "shape",
               SystMap<>::init(0.71));

  cb.cp()
      .channel({"mt"})
      .process({"ZTT", "TTT", "TTL", "VVT", "VVL", "jetFakes"})
      .AddSyst(cb, "CMS_scale_mc_t_3prong", "shape", SystMap<>::init(0.71));

  // Embedded uncorrelated uncertainty

  cb.cp()
      .channel({"mt"})
      .process({"EMB", "jetFakes"})
      .AddSyst(cb, "CMS_scale_emb_t_1prong", "shape", SystMap<>::init(0.71));

  cb.cp()
      .channel({"mt"})
      .process({"EMB", "jetFakes"})
      .AddSyst(cb, "CMS_scale_emb_t_1prong1pizero", "shape", SystMap<>::init(0.71));

  cb.cp()
      .channel({"mt"})
      .process({"EMB", "jetFakes"})
      .AddSyst(cb, "CMS_scale_emb_t_3prong", "shape", SystMap<>::init(0.71));

  // MC + embedded correlated uncertainty

  cb.cp()
      .channel({"mt"})
      .process(JoinStr({{"ZTT", "TTT", "TTL", "VVT", "VVL", "jetFakes"}, {"EMB"}}))
      .AddSyst(cb, "CMS_scale_t_1prong", "shape", SystMap<>::init(0.71));

  cb.cp()
      .channel({"mt"})
      .process(JoinStr({{"ZTT", "TTT", "TTL", "VVT", "VVL", "jetFakes"}, {"EMB"}}))
      .AddSyst(cb, "CMS_scale_t_1prong1pizero", "shape", SystMap<>::init(0.71));

  cb.cp()
      .channel({"mt"})
      .process(JoinStr({{"ZTT", "TTT", "TTL", "VVT", "VVL", "jetFakes"}, {"EMB"}}))
      .AddSyst(cb, "CMS_scale_t_3prong", "shape", SystMap<>::init(0.71));

  // ##########################################################################
  // Uncertainty: Jet energy scale
  // References:
  // - Talk in CMS Htt meeting by Daniel Winterbottom about regional JES splits:
  //   https://indico.cern.ch/event/740094/contributions/3055870/
  // Notes:
  // ##########################################################################

  if (!regional_jec) {
  cb.cp()
      .channel({"mt"})
      .process(mc_processes)
      .AddSyst(cb, "CMS_scale_j", "shape", SystMap<>::init(1.00));
  }

  // Regional JES
  else {
    cb.cp()
        .channel({"mt"})
        .process(mc_processes)
        .AddSyst(cb, "CMS_scale_j_eta0to3", "shape", SystMap<>::init(1.00));

    cb.cp()
        .channel({"mt"})
        .process(mc_processes)
        .AddSyst(cb, "CMS_scale_j_eta0to5", "shape", SystMap<>::init(1.00));

    cb.cp()
        .channel({"mt"})
        .process(mc_processes)
        .AddSyst(cb, "CMS_scale_j_eta3to5", "shape", SystMap<>::init(1.00));

    cb.cp()
        .channel({"mt"})
        .process(mc_processes)
        .AddSyst(cb, "CMS_scale_j_RelativeBal", "shape", SystMap<>::init(1.00));
        
    if (era == 2017) {
      cb.cp()
	  .channel({"mt"})
          .process(mc_processes)
          .AddSyst(cb, "CMS_scale_j_RelativeSample_$ERA", "shape", SystMap<>::init(1.00));
    }
  }

  // ##########################################################################
  // Uncertainty: MET energy scale and Recoil
  // References:
  // Notes:
  // - FIXME: Clustered vs unclustered MET? Inclusion of JES splitting?
  // - FIXME: References?
  // ##########################################################################

  cb.cp()
      .channel({"mt"})
      .process({"ZTT", "TT", "TTT", "TTL", "TTJ", "W", "ZJ", "ZL", "VV", "VVT", "VVL", "VVJ", "ST"})  //Z and W processes are only included due to the EWK fraction. Make sure that there is no contribution to the shift from the DY or Wjets samples.
      .AddSyst(cb, "CMS_scale_met_unclustered", "shape", SystMap<>::init(1.00));
  cb.cp()
      .channel({"mt"})
      .process({"ZTT", "ZL", "ZJ", "W"})
      .AddSyst(cb, "CMS_htt_boson_scale_met_$ERA", "shape", SystMap<>::init(1.00));
  cb.cp()
      .channel({"mt"})
     .process({"ZTT", "ZL", "ZJ", "W"})
      .AddSyst(cb, "CMS_htt_boson_reso_met_$ERA", "shape", SystMap<>::init(1.00));

  // ##########################################################################
  // Uncertainty: Background normalizations
  // References:
  // Notes:
  // - FIXME: Remeasure QCD extrapolation factors for SS and ABCD methods?
  //          Current values are measured by KIT.
  // - FIXME: Adapt for fake factor and embedding
  // - FIXME: W uncertainties: Do we need lnN uncertainties based on the Ersatz
  //          study in Run1 (found in HIG-16043 uncertainty model)
  // - FIXME: References?
  // ##########################################################################

  // VV
  cb.cp()
      .channel({"mt", "mm"})
      .process({"VVT", "VVJ", "VVL", "VV"})
      .AddSyst(cb, "CMS_htt_vvXsec", "lnN", SystMap<>::init(1.06));

  // TT
  cb.cp()
      .channel({"mt", "mm"})
      .process({"TTT", "TTL", "TTJ", "TT"})
      .AddSyst(cb, "CMS_htt_tjXsec", "lnN", SystMap<>::init(1.055));

  cb.cp()
      .channel({"mt", "mm"})
      .process({"ST"})
      .AddSyst(cb, "CMS_htt_stXsec", "lnN", SystMap<>::init(1.055));

  // W
  cb.cp()
      .channel({"mt"})
      .process({"W"})
      .AddSyst(cb, "CMS_htt_wjXsec", "lnN", SystMap<>::init(1.05));
  // TODO: Check if this really should be uncorrelated
  cb.cp()
      .channel({"mm"})
      .process({"W"})
      .AddSyst(cb, "CMS_htt_wjXsec_mm", "lnN", SystMap<>::init(1.05));

  // Z
  cb.cp()
      .channel({"mt"})
      .process({"ZTT", "ZL", "ZJ"})
      .AddSyst(cb, "CMS_htt_zjXsec", "rateParam", SystMap<>::init(0.976460));
  cb.cp()
      .channel({"mm"})
      .process({"ZLL"})
      .AddSyst(cb, "CMS_htt_zjXsec", "rateParam", SystMap<>::init(0.976460));

  // QCD
  cb.cp()
      .channel({"mt"})
      .process({"QCD"})
      .AddSyst(cb, "CMS_ExtrapSSOS_$CHANNEL_$ERA", "lnN", SystMap<>::init(1.30));

  // ##########################################################################
  // Uncertainty: Drell-Yan LO->NLO reweighting
  // References:
  // Notes:
  // - FIXME: References?
  // ##########################################################################

  cb.cp()
      .channel({"mt"})
      .process({"ZTT", "ZL", "ZJ"})
      .AddSyst(cb, "CMS_htt_dyShape_$ERA", "shape", SystMap<>::init(0.10));

  // ##########################################################################
  // Uncertainty: TT shape reweighting
  // References:
  // Notes:
  // - FIXME: References?
  // ##########################################################################

  cb.cp()
      .channel({"mt"})
      .process({"TTT", "TTL", "TTJ", "TT"})
      .AddSyst(cb, "CMS_htt_ttbarShape", "shape", SystMap<>::init(1.00));

  // ##########################################################################
  // Uncertainty: Electron/muon to tau fakes and ZL energy scale
  // References:
  // Notes:
  // - FIXME: References?
  // ##########################################################################

  // ZL energy scale splitted by decay mode
  cb.cp()
      .channel({"mt"})
      .process({"ZL"})
      .AddSyst(cb, "CMS_ZLShape_$CHANNEL_1prong_$ERA", "shape",
               SystMap<>::init(1.00));

  cb.cp()
      .channel({"mt"})
      .process({"ZL"})
      .AddSyst(cb, "CMS_ZLShape_$CHANNEL_1prong1pizero_$ERA", "shape",
               SystMap<>::init(1.00));

  // Muon fakes
  if (era == 2016 || era == 2017) {
      cb.cp()
          .channel({"mt"})
          .process({"ZL"})
          .AddSyst(cb, "CMS_mFakeTau", "lnN", SystMap<>::init(1.500));
  } else if (era == 2018) {
      cb.cp()
          .channel({"mt"})
          .process({"ZL"})
          .AddSyst(cb, "CMS_mFakeTau", "lnN", SystMap<>::init(2.000));
  }

  // ##########################################################################
  // Uncertainty: Jet to tau fakes
  // References:
  // Notes:
  // - FIXME: Adapt for fake factor and embedding
  // - FIXME: References?
  // ##########################################################################

  cb.cp()
      .channel({"mt"})
      .process({"W", "TTJ", "ZJ", "VVJ"})
      .AddSyst(cb, "CMS_htt_jetToTauFake_$ERA", "shape", SystMap<>::init(1.00));

  // ##########################################################################
  // Uncertainty: Embedded events
  // References:
  // - https://twiki.cern.ch/twiki/bin/viewauth/CMS/TauTauEmbeddingSamples2016
  // Notes:
  // ##########################################################################

  // Embedded Normalization: No Lumi, Zjxsec information used, instead derived from data using dimuon selection efficiency
  cb.cp()
      .channel({"mt"})
      .process({"EMB"})
      .AddSyst(cb, "CMS_htt_doublemutrg_$ERA", "lnN", SystMap<>::init(1.04));

  // TTbar contamination in embedded events: 10% shape uncertainty of assumed ttbar->tautau event shape
  cb.cp()
    .channel({"mt"})
    .process({"EMB"})
    .AddSyst(cb, "CMS_htt_emb_ttbar_$ERA", "shape", SystMap<>::init(1.00));

  // Uncertainty of hadronic tau track efficiency correction
  // uncorrelated between eras
  cb.cp()
    .channel({"mt"})
    .process({"EMB"})
    .AddSyst(cb, "CMS_3ProngEff", "shape", SystMap<>::init(1.00));

  cb.cp()
    .channel({"mt"})
    .process({"EMB"})
    .AddSyst(cb, "CMS_1ProngPi0Eff", "shape", SystMap<>::init(1.00));

  // ##########################################################################
  // Uncertainty: Jet fakes
  // References:
  // - https://twiki.cern.ch/twiki/bin/viewauth/CMS/HiggsToTauTauJet2TauFakes
  // Notes:
  // - FIXME: add 2017 norm uncertainties, and properly correlate across years
  // ##########################################################################

  // QCD shape stat.
  cb.cp()
      .channel({"mt"})
      .process({"jetFakes"})
      .AddSyst(cb, "CMS_ff_qcd_njet0_$CHANNEL_stat_$ERA", "shape", SystMap<>::init(1.00));
  cb.cp()
      .channel({"mt"})
      .process({"jetFakes"})
      .AddSyst(cb, "CMS_ff_qcd_njet1_$CHANNEL_stat_$ERA", "shape", SystMap<>::init(1.00));

  // W shape stat.
  cb.cp()
      .channel({"mt"})
      .process({"jetFakes"})
      .AddSyst(cb, "CMS_ff_w_njet0_$CHANNEL_stat_$ERA", "shape", SystMap<>::init(1.00));
  cb.cp()
      .channel({"mt"})
      .process({"jetFakes"})
      .AddSyst(cb, "CMS_ff_w_njet1_$CHANNEL_stat_$ERA", "shape", SystMap<>::init(1.00));

  // TT shape stat.
  cb.cp()
      .channel({"mt"})
      .process({"jetFakes"})
      .AddSyst(cb, "CMS_ff_tt_njet1_stat_$ERA", "shape", SystMap<>::init(1.00));

  // Shape syst. of different contributions (QCD/W/tt)
  cb.cp()
      .channel({"mt"})
      .process({"jetFakes"})
      .AddSyst(cb, "CMS_ff_qcd_$CHANNEL_syst", "shape", SystMap<>::init(1.00));
  cb.cp()
      .channel({"mt"})
      .process({"jetFakes"})
      .AddSyst(cb, "CMS_ff_w_syst", "shape", SystMap<>::init(1.00));
  cb.cp()
      .channel({"mt"})
      .process({"jetFakes"})
      .AddSyst(cb, "CMS_ff_tt_syst", "shape", SystMap<>::init(1.00));

  //below: jetFakes norm uncertainties. Current values are for 2016, which are probably a good approx. for 2017. To be updated.

  // Stat. norm (uncorrelated across years)
 //  cb.cp()
 //      .channel({"et", "mt", "tt"})
 //      .process({"jetFakes"})
 //      .AddSyst(cb, "CMS_ff_norm_stat_$CHANNEL_$BIN_$ERA", "lnN", SystMap<channel, bin_id>::init
 //               ({"mt"}, {1},   1.035) //ggh
 //               ({"mt"}, {2},   1.048) //qqh
 //               ({"mt"}, {11},  1.028) //w
 //               ({"mt"}, {12},  1.037) //ztt
 //               ({"mt"}, {13},  1.036) //tt
 //               ({"mt"}, {14},  1.033) //ss
 //               ({"mt"}, {15},  1.028) //zll
 //               ({"mt"}, {16},  1.042) //misc
 //               ({"mt"}, {100}, 1.026) //incl
 //               ({"et"}, {1},   1.052) //ggh
 //               ({"et"}, {2},   1.079) //qqh
 //               ({"et"}, {11},  1.047) //w
 //               ({"et"}, {12},  1.067) //ztt
 //               ({"et"}, {13},  1.059) //tt
 //               ({"et"}, {14},  1.038) //ss
 //               ({"et"}, {15},  1.067) //zll
 //               ({"et"}, {16},  1.076) //misc
 //               ({"et"}, {100}, 1.046) //incl
 //               ({"tt"}, {1},   1.029) //ggh
 //               ({"tt"}, {2},   1.037) //qqh
 //               ({"tt"}, {12},  1.035) //ztt
 //               ({"tt"}, {16},  1.020) //misc
 //               ({"tt"}, {17},  1.029) //noniso
 //               ({"tt"}, {100}, 1.029) //incl
 //               );

 //  // Syst. norm: Bin-correlated
 //  // uncorrelated between eras
 //  cb.cp()
 //      .channel({"et", "mt", "tt"})
 //      .process({"jetFakes"})
 //      .AddSyst(cb, "CMS_ff_norm_syst_$CHANNEL_$ERA", "lnN", SystMap<channel, bin_id>::init
 //               ({"mt"}, {1},   1.049) //ggh
 //               ({"mt"}, {2},   1.041) //qqh
 //               ({"mt"}, {11},  1.038) //w
 //               ({"mt"}, {12},  1.069) //ztt
 //               ({"mt"}, {13},  1.037) //tt
 //               ({"mt"}, {14},  1.064) //ss
 //               ({"mt"}, {15},  1.048) //zll
 //               ({"mt"}, {16},  1.064) //misc
 //               ({"mt"}, {100}, 1.042) //incl
 //               ({"et"}, {1},   1.042) //ggh
 //               ({"et"}, {2},   1.040) //qqh
 //               ({"et"}, {11},  1.037) //w
 //               ({"et"}, {12},  1.062) //ztt
 //               ({"et"}, {13},  1.040) //tt
 //               ({"et"}, {14},  1.045) //ss
 //               ({"et"}, {15},  1.051) //zll
 //               ({"et"}, {16},  1.041) //misc
 //               ({"et"}, {100}, 1.042) //incl
 //               ({"tt"}, {1},   1.068) //ggh
 //               ({"tt"}, {2},   1.067) //qqh
 //               ({"tt"}, {12},  1.067) //ztt
 //               ({"tt"}, {16},  1.078) //misc
 //               ({"tt"}, {17},  1.070) //noniso
 //               ({"tt"}, {100}, 1.067) //incl
 //               );
 //  // correlated between eras
 //  cb.cp()
 //      .channel({"et", "mt", "tt"})
 //      .process({"jetFakes"})
 //      .AddSyst(cb, "CMS_ff_norm_syst_$CHANNEL", "lnN", SystMap<channel, bin_id>::init
 //               ({"mt"}, {1},   1.049) //ggh
 //               ({"mt"}, {2},   1.041) //qqh
 //               ({"mt"}, {11},  1.038) //w
 //               ({"mt"}, {12},  1.069) //ztt
 //               ({"mt"}, {13},  1.037) //tt
 //               ({"mt"}, {14},  1.064) //ss
 //               ({"mt"}, {15},  1.048) //zll
 //               ({"mt"}, {16},  1.064) //misc
 //               ({"mt"}, {100}, 1.042) //incl
 //               ({"et"}, {1},   1.042) //ggh
 //               ({"et"}, {2},   1.040) //qqh
 //               ({"et"}, {11},  1.037) //w
 //               ({"et"}, {12},  1.062) //ztt
 //               ({"et"}, {13},  1.040) //tt
 //               ({"et"}, {14},  1.045) //ss
 //               ({"et"}, {15},  1.051) //zll
 //               ({"et"}, {16},  1.041) //misc
 //               ({"et"}, {100}, 1.042) //incl
 //               ({"tt"}, {1},   1.068) //ggh
 //               ({"tt"}, {2},   1.067) //qqh
 //               ({"tt"}, {12},  1.067) //ztt
 //               ({"tt"}, {16},  1.078) //misc
 //               ({"tt"}, {17},  1.070) //noniso
 //               ({"tt"}, {100}, 1.067) //incl
 //               );

 //  // Syst. norm: Bin-dependent, correlated across years
 //  // uncorrelated between eras
 //  cb.cp()
 //      .channel({"et", "mt", "tt"})
 //      .process({"jetFakes"})
 //      .AddSyst(cb, "CMS_ff_sub_syst_$CHANNEL_$BIN_$ERA", "lnN", SystMap<channel, bin_id>::init
 //               ({"mt"}, {1},   1.028) //ggh
 //               ({"mt"}, {2},   1.028) //qqh
 //               ({"mt"}, {11},  1.018) //w
 //               ({"mt"}, {12},  1.032) //ztt
 //               ({"mt"}, {13},  1.021) //tt
 //               ({"mt"}, {14},  1.014) //ss
 //               ({"mt"}, {15},  1.028) //zll
 //               ({"mt"}, {16},  1.025) //misc
 //               ({"mt"}, {100}, 1.025) //incl
 //               ({"et"}, {1},   1.028) //ggh
 //               ({"et"}, {2},   1.025) //qqh
 //               ({"et"}, {11},  1.014) //w
 //               ({"et"}, {12},  1.028) //ztt
 //               ({"et"}, {13},  1.021) //tt
 //               ({"et"}, {14},  1.014) //ss
 //               ({"et"}, {15},  1.028) //zll
 //               ({"et"}, {16},  1.025) //misc
 //               ({"et"}, {100}, 1.025) //incl
 //               ({"tt"}, {1},   1.021) //ggh
 //               ({"tt"}, {2},   1.021) //qqh
 //               ({"tt"}, {12},  1.025) //ztt
 //               ({"tt"}, {16},  1.021) //misc
 //               ({"tt"}, {17},  1.014) //noniso
 //               ({"tt"}, {100}, 1.021) //incl
 //               );
 //  // correlated between eras
 //  cb.cp()
 //      .channel({"et", "mt", "tt"})
 //      .process({"jetFakes"})
 //      .AddSyst(cb, "CMS_ff_sub_syst_$CHANNEL_$BIN", "lnN", SystMap<channel, bin_id>::init
 //               ({"mt"}, {1},   1.028) //ggh
 //               ({"mt"}, {2},   1.028) //qqh
 //               ({"mt"}, {11},  1.018) //w
 //               ({"mt"}, {12},  1.032) //ztt
 //               ({"mt"}, {13},  1.021) //tt
 //               ({"mt"}, {14},  1.014) //ss
 //               ({"mt"}, {15},  1.028) //zll
 //               ({"mt"}, {16},  1.025) //misc
 //               ({"mt"}, {100}, 1.025) //incl
 //               ({"et"}, {1},   1.028) //ggh
 //               ({"et"}, {2},   1.025) //qqh
 //               ({"et"}, {11},  1.014) //w
 //               ({"et"}, {12},  1.028) //ztt
 //               ({"et"}, {13},  1.021) //tt
 //               ({"et"}, {14},  1.014) //ss
 //               ({"et"}, {15},  1.028) //zll
 //               ({"et"}, {16},  1.025) //misc
 //               ({"et"}, {100}, 1.025) //incl
 //               ({"tt"}, {1},   1.021) //ggh
 //               ({"tt"}, {2},   1.021) //qqh
 //               ({"tt"}, {12},  1.025) //ztt
 //               ({"tt"}, {16},  1.021) //misc
 //               ({"tt"}, {17},  1.014) //noniso
 //               ({"tt"}, {100}, 1.021) //incl
 //               );
}
} // namespace ch
