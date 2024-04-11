import argparse
import awkward as ak
import numpy as np
import ROOT
import uproot

def createCMSLabel():
	cmsLatex = ROOT.TLatex()
	cmsLatex.SetTextSize(0.04)
	cmsLatex.SetNDC(True)
	cmsLatex.SetTextAlign(11)

	return cmsLatex

def main(file_path, out_dir):
    
    inFile = ROOT.TFile.Open(file_path, "READ")
   
    c = ROOT.TCanvas( 'c', 'c', 100, 10, 600, 600 )
    h_SingleJet180_den = inFile.Get("demo/h_SingleJet180_den")
    h_SingleJet180_num = inFile.Get("demo/h_SingleJet180_num")
    h_SingleJet180_den.Sumw2()
    h_SingleJet180_num.Sumw2()
    eff_SingleJet180 = h_SingleJet180_num.Clone("eff_SingleJet180")
    eff_SingleJet180.Divide(h_SingleJet180_den)
    eff_SingleJet180.SetLineColor(4)
    eff_SingleJet180.SetMarkerColor(4)
    eff_SingleJet180.SetStats(0)
    eff_SingleJet180.SetTitle("SingleJet180 Efficiency; Offline Jet pT [GeV]; Efficiency")
    eff_SingleJet180.Draw("PE")


    c.Draw()
    c.SaveAs(f"{out_dir}/eff_SingleJet180.png")
    c.Close() 
    '''
    for l1thresh in ["30", "60", "120", "180"]:
        
        c1 = ROOT.TCanvas( 'c1', 'c1', 100, 10, 600, 600 )
        h_SingleJet180_num = inFile.Get(f"h_SingleJet180_L1{l1thresh}_num")
        h_SingleJet180_num.Sumw2()
        eff_SingleJet180 = h_SingleJet180_num.Clone(f"eff_SingleJet180_L1{l1thresh}")
        eff_SingleJet180.Divide(h_SingleJet180_den)
        eff_SingleJet180.SetLineColor(4)
        eff_SingleJet180.SetMarkerColor(4)
        eff_SingleJet180.SetStats(0)
        eff_SingleJet180.SetTitle(f"SingleJet180 Efficiency (L1 Jet pT > {l1thresh} GeV); Offline Jet pT [GeV]; Efficiency")
        eff_SingleJet180.Draw("PE")

        c1.Draw()
        c1.SaveAs(f"{out_dir}/eff_SingleJet180_L1{l1thresh}.png")
        c1.Close()
    '''
    c2 = ROOT.TCanvas( 'c2', 'c2', 100, 10, 600, 600 )
    h_SingleJet180_central_den = inFile.Get("demo/h_SingleJet180_central_den")
    h_SingleJet180_central_num = inFile.Get("demo/h_SingleJet180_central_num")
    h_SingleJet180_central_den.Sumw2()
    h_SingleJet180_central_num.Sumw2()
    eff_SingleJet180_central = h_SingleJet180_central_num.Clone("eff_SingleJet180_central")
    eff_SingleJet180_central.Divide(h_SingleJet180_den)
    eff_SingleJet180_central.SetLineColor(4)
    eff_SingleJet180_central.SetMarkerColor(4)
    eff_SingleJet180_central.SetStats(0)
    eff_SingleJet180_central.SetTitle("SingleJet180 Efficiency (|eta| <= 3.0); Offline Jet pT [GeV]; Efficiency")
    eff_SingleJet180_central.Draw("PE")
    c2.Draw()
    c2.SaveAs(f"{out_dir}/eff_SingleJet180_central")
    c2.Close()

    c3 = ROOT.TCanvas( 'c3', 'c3', 100, 10, 600, 600 )
    h_SingleJet180_forward_den = inFile.Get("demo/h_SingleJet180_forward_den")
    h_SingleJet180_forward_num = inFile.Get("demo/h_SingleJet180_forward_num")
    h_SingleJet180_forward_den.Sumw2()
    h_SingleJet180_forward_num.Sumw2()
    eff_SingleJet180_forward = h_SingleJet180_forward_num.Clone("eff_SingleJet180_forward")
    eff_SingleJet180_forward.Divide(h_SingleJet180_den)
    eff_SingleJet180_forward.SetLineColor(4)
    eff_SingleJet180_forward.SetMarkerColor(4)
    eff_SingleJet180_forward.SetStats(0)
    eff_SingleJet180_forward.SetTitle("SingleJet180 Efficiency (|eta| > 3.0); Offline Jet pT [GeV]; Efficiency")
    eff_SingleJet180_forward.Draw("PE")
    c3.Draw()
    c3.SaveAs(f"{out_dir}/eff_SingleJet180_forward")
    c3.Close()

    c4 = ROOT.TCanvas( 'c4', 'c4', 100, 10, 600, 600 )
    h_HT280_num = inFile.Get("demo/h_HT280_num")
    h_HT280_den = inFile.Get("demo/h_HT280_den")
    h_HT280_num.Sumw2()
    h_HT280_den.Sumw2()
    eff_HT280 = h_HT280_num.Clone("eff_HT280")
    eff_HT280.Divide(h_HT280_den)
    eff_HT280.SetLineColor(4)
    eff_HT280.SetMarkerColor(4)
    eff_HT280.SetStats(0)
    eff_HT280.SetTitle("HT280 Efficiency; $H_T$ [GeV]; Efficiency")
    eff_HT280.Draw("PE")
    c4.Draw()
    c4.SaveAs(f"{out_dir}/eff_HT280")
    c4.Close()

    c5 = ROOT.TCanvas( 'c5', 'c5', 100, 10, 600, 600 )
    h_ETMHF90_num = inFile.Get("demo/h_ETMHF90_num")
    h_ETMHF90_den = inFile.Get("demo/h_ETMHF90_den")
    h_ETMHF90_num.Sumw2()
    h_ETMHF90_den.Sumw2()
    eff_ETMHF90 = h_ETMHF90_num.Clone("eff_ETMHF90")
    eff_ETMHF90.Divide(h_ETMHF90_den)
    eff_ETMHF90.SetLineColor(4)
    eff_ETMHF90.SetMarkerColor(4)
    eff_ETMHF90.SetStats(0)
    eff_ETMHF90.SetTitle("ETMHF90 Efficiency; Missing $p_T$ [GeV]; Efficiency")
    eff_ETMHF90.Draw("PE")
    c5.Draw()
    c5.SaveAs(f"{out_dir}/eff_ETMHF90")
    c5.Close()
    
if __name__ == "__main__":
    
    # parse arguments
    parser = argparse.ArgumentParser(description="This program creates plots to study L1 jets.")
    parser.add_argument("-p", "--file_path", help="path to input ROOT file containing histograms")
    parser.add_argument("-o", "--out_dir", help="path to directory to save plots", default=".")
    
    args = parser.parse_args()
    
    main(args.file_path, args.out_dir)
