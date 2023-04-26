#include"Plotter.cc"
class PDFStudy:public Plotter{
public:
  PDFStudy(TString mode_="");
};
PDFStudy::PDFStudy(TString mode_){
  ScanFiles("hists");
  samples["nnpdf31nlo"]=Sample("NNPDF31 nlo luxqed","SAMPLE nominal nnpdf31",Style(kBlack,-1,3266,"e2"),Style(kBlack,-1,3266,"e2"))+"pdfstudy";
  samples["nnpdf31nnlo"]=Sample("NNPDF31 nnlo luxqed","SAMPLE nnpdf31",Style(kRed,-1,3277,"e2"),Style(kRed,-1,3277,"e2"))%"_nnpdf31nnlo"+"pdfstudy";
  samples["mmht2015"]=Sample("MMHT2015 qed nlo","SAMPLE mmht2015",Style(kGreen,-1,3344,"e2"),Style(kGreen,-1,3344,"e2"))%"_mmht2015"+"pdfstudy";
  samples["pdf4lhc15"]=Sample("LUXqed17 + PDF4LHC15","SAMPLE pdf4lhc15",Style(kBlue,-1,3333,"e2"),Style(kBlue,-1,3333,"e2"))%"_mmht2015"+"pdfstudy";
  samples["nnpdf31nlo_that"]=Sample("Q^{2}=#hat{t}","SAMPLE nominal nnpdf31",Style(kOrange,-1,3277,"e2"),Style(kOrange,-1,3277,"e2"))%"_that"+"pdfstudy";
  samples["nnpdf31nlo_shat"]=Sample("Q^{2}=#hat{s}","SAMPLE nominal nnpdf31",Style(kMagenta,-1,3233,"e2"),Style(kMagenta,-1,3233,"e2"))%"_shat"+"pdfstudy";
  samples["mmht2015qed1"]=Sample("MMHT2015 qed1 nlo","SAMPLE mmht2015",Style(kRed,22,-1,"e"),Style(kRed,-1,-1,"e"))+"MMHT2015QED_0";
  samples["mmht2015qed2"]=Sample("MMHT2015 qed2 nlo","SAMPLE mmht2015",Style(kBlue,23,-1,"e"),Style(kBlue,-1,-1,"e"))+"MMHT2015QED_1"+"MMHT2015QED_2"+"MMHT2015QED_3"+"MMHT2015QED_4";
  samples["mmht2015qed3"]=Sample("MMHT2015 qed2 nlo","SAMPLE mmht2015",Style(kGreen,24,-1,"e"),Style(kGreen,-1,-1,"e"))+"MMHT2015QED_1"+"MMHT2015QED_5"+"MMHT2015QED_6"+"MMHT2015QED_7";

  AddSystematic("muf","muf",Systematic::Type::ENVELOPE,"_muf_up _muf_down","nominal");
  AddSystematic("mur","mur",Systematic::Type::ENVELOPE,"_mur_up _mur_down","nominal");

  AddSystematic("pdf_nnpdf31","pdf_nnpdf31",Systematic::Type::GAUSSIAN,FormRange("_pdf%d",Range(100)),"nnpdf31");

  AddSystematic("pdf_ct14","pdf_ct14",Systematic::Type::ENVELOPE,"_pdf5->_pdf11 _pdf5->_pdf0","ct14");

  AddSystematic("mmht2015_qcd","mmht2015_qcd",Systematic::Type::HESSIAN,FormRange("_pdf%d",Range(1,51,2)),"mmht2015");
  for(int i=0;i<6;i++){
    AddSystematic(Form("mmht2015_qed%d",i),Form("mmht2015_qed%d",i),Systematic::Type::ENVELOPE,Form("_pdf%d _pdf%d",51+2*i,52+2*i),"mmht2015");
  }
  AddSystematic("mmht2015_qed","mmht2015_qed",Systematic::Type::MULTI,FormRange("mmht2015_qed%d",Range(6)),"mmht2015");
  AddSystematic("pdf_mmht2015","pdf_mmht2015",Systematic::Type::MULTI,"mmht2015_qcd mmht2015_qed","mmht2015");

  AddSystematic("pdf4lhc15_qcd","pdf4lhc15_qcd",Systematic::Type::HESSIAN,FormRange("_pdf%d",Range(1,31)),"pdf4lhc15");
  AddSystematic("pdf4lhc15_qed","pdf4lhc15_qed",Systematic::Type::HESSIAN,FormRange("_pdf%d",Range(31,38)),"pdf4lhc15");
  AddSystematic("pdf_pdf4lhc15","pdf_pdf4lhc15",Systematic::Type::MULTI,"pdf4lhc15_qcd pdf4lhc15_qed","pdf4lhc15");

  //AddSystematic("pdf","PDF",Systematic::Type::MULTI,"pdf_nnpdf31 pdf_ct14 pdf_mmht2015 pdf_pdf4lhc15");
  AddSystematic("pdf","PDF",Systematic::Type::MULTI,"pdf_nnpdf31 pdf_ct14 mmht2015_qcd pdf4lhc15_qcd");

  SetupEntries(mode_);
}
