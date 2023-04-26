#include "PDFStudy.cc"
void plot(){
  PDFStudy aa;
  TString common="widthweight nostat type:2 base:0";
  aa.SetupEntries("nnpdf31nlo mmht2015");
  aa.DrawPlot("mass_ac","xmin:50 xmax:3000 logx xtitle:m(ll) sysname:pdf "+common);
  aa.SetupEntries("nnpdf31nlo pdf4lhc15");
  aa.DrawPlot("mass_ac","xmin:50 xmax:3000 logx xtitle:m(ll) sysname:pdf "+common);
  aa.SetupEntries("nnpdf31nlo nnpdf31nnlo");
  aa.DrawPlot("mass_ac","xmin:50 xmax:3000 logx xtitle:m(ll) sysname:pdf "+common);
  aa.SetupEntries("nnpdf31nlo");
  aa.DrawPlot("mass_ac","xmin:50 xmax:3000 logx xtitle:m(ll) sysname:pdf widey "+common);
  aa.DrawPlot("mass_ac","xmin:50 xmax:3000 logx xtitle:m(ll) sysname:muf widey "+common);
  aa.DrawPlot("mass_ac","xmin:50 xmax:3000 logx xtitle:m(ll) sysname:mur widey "+common);
  aa.DrawPlot("rapidity_ac","xtitle:y(ll) xmin:-3 xmax:3 sysname:muf widey "+common);
  aa.DrawPlot("rapidity_ac","xtitle:y(ll) xmin:-3 xmax:3 sysname:pdf widey "+common);
  aa.DrawPlot("pt_ac","xtitle:p_{T}(ll) logx sysname:muf widey "+common);
  aa.DrawPlot("pt_ac","xtitle:p_{T}(ll) logx sysname:pdf widey "+common);
  aa.DrawPlot("cost_ac","xtitle:cos#theta_{CS} sysname:muf widey "+common);
  aa.DrawPlot("cost_ac","xtitle:cos#theta_{CS} sysname:pdf widey "+common);
  aa.SetupEntries("nnpdf31nlo nnpdf31nlo_that nnpdf31nlo_shat");
  aa.DrawPlot("mass_ac","xmin:50 xmax:3000 logx xtitle:m(ll) sysname:muf widey "+common);

}
void plot_old(){
  Plotter aa;
  aa.ScanFiles("hists");  
  aa.SetupEntries("NNPDF31Luxqed NNPDF31Luxqed_GlobalRecoil");
  //aa.SetupEntries("NNPDF31Luxqed NNPDF23QED_DD");
  //aa.SetupEntries("NNPDF23QED NNPDF23QED_DD NNPDF23QED_SD NNPDF23QED_ND NNPDF31Luxqed_GlobalRecoil NNPDF31Luxqed ");
  aa.PrintEntries();
  TString common="Xmin:50 Xmax:3000 logx widthweight 2:widewidey";
  //TString common="Xmin:50 Xmax:3000 xmin:0 xmax:100 widthweight 2:widewidey";
  aa.DrawPlot("mpt","project:x xmin:50 xmax:3000 1:logy "+common);
  aa.DrawPlot("mpt","project:y 1:logy "+common);
  aa.DrawPlot("mpt_ac","project:x xmin:50 xmax:3000 1:logy "+common);
  aa.DrawPlot("mpt_ac","project:y 1:logy "+common);
}
  
  
