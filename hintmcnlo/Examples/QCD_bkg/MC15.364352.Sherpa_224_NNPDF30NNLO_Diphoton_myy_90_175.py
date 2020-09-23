include("MC15JobOptions/Sherpa_2.2.4_NNPDF30NNLO_Common.py")

evgenConfig.description = "Sherpa 2.2.4 gamma gamma + 0,1j@NLO + 2,3j@LO with 90<m_yy<175."
evgenConfig.keywords = ["SM", "diphoton", "NLO" ]
evgenConfig.contact  = [ "atlas-generators-sherpa@cern.ch", "frank.siegert@cern.ch" ]
evgenConfig.minevents = 500
evgenConfig.inputconfcheck = "Diphoton_myy_90_175"

genSeq.Sherpa_i.RunCard="""
(run){
  % scales, tags for scale variations
  FSF:=1.; RSF:=1.; QSF:=1.;
  SCALES STRICT_METS{FSF*MU_F2}{RSF*MU_R2}{QSF*MU_Q2};

  ALPHAQED_DEFAULT_SCALE=0.0

  % tags for process setup
  NJET:=3; LJET:=2,3; QCUT:=10;

  % me generator settings
  ME_SIGNAL_GENERATOR Amegic Comix LOOPGEN Internal;
  LOOPGEN:=OpenLoops;
}(run)

(processes){
  Process 21 21 -> 22 22
  ME_Generator Internal;
  Loop_Generator gg_yy
  Scales VAR{FSF*Abs2(p[2]+p[3])}{RSF*Abs2(p[2]+p[3])}{QSF*Abs2(p[2]+p[3])};
  End process;

  Process 93 93 -> 22 22 93{NJET};
  Order (*,2);
  CKKW sqr(QCUT/E_CMS)/(1.0+sqr(QCUT/0.6)/Abs2(p[2]+p[3]));
  NLO_QCD_Mode MC@NLO {LJET};
  ME_Generator Amegic {LJET};
  RS_ME_Generator Comix {LJET};
  Loop_Generator LOOPGEN {LJET};
  PSI_ItMin 50000 {3}
  Integration_Error 0.99 {3}
  PSI_ItMin 100000 {4,5}
  Integration_Error 0.02 {4,5}
  End process;
}(processes)

(selector){
  "PT" 22 20,E_CMS:18,E_CMS [PT_UP]
  RapidityNLO  22  -2.7  2.7
  IsolationCut  22  0.1  2  0.10;
  Mass 22 22  90  175
  DeltaRNLO 22 22 0.2 1000.0;
}(selector)
"""

# switch to CutTools reduction within OpenLoops
genSeq.Sherpa_i.Parameters += [ "OL_PARAMETERS=redlib1=5=redlib2=5=write_parameters=1" ]

# add PDF4LHC15
genSeq.Sherpa_i.Parameters += [
  "PDF_VARIATIONS=NNPDF30NNLO[all] NNPDF30_nnlo_as_0117 NNPDF30_nnlo_as_0119 MMHT2014nnlo68cl CT14nnlo PDF4LHC15_nnlo_30[all]"
]

genSeq.Sherpa_i.OpenLoopsLibs = [ "ppaa", "ppaaj" ]
genSeq.Sherpa_i.NCores = 240

