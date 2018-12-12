import ROOT
ROOT.gStyle.SetOptStat(0)
from math import sqrt
import random
import DataSetInfo as info
import optparse
import os

class WriteNJetPlots:
    def __init__(self):    
        self.colors = [ROOT.kCyan+1, ROOT.kMagenta+1, ROOT.kYellow+1,
                       ROOT.kRed+1, ROOT.kGreen+1, ROOT.kBlue+1,
                       ROOT.kBlack]
        self.signalcolors = [ROOT.kViolet-6, ROOT.kTeal-6, ROOT.kPink-6, ROOT.kAzure-6]
        
    def calcUnc(self,a,unc_a,b,unc_b):
        return  a/b * sqrt( (unc_a/a)**2+(unc_b/b)**2)
    
    def writeHistos(self, data, basenameIn, basenameOut, bin, sys):
        histos = []
        for key, dsi in data.iteritems():
            h = dsi.getHisto(basenameIn+"_"+bin+sys)
            h.SetName(bin+"_"+dsi.label+"_"+basenameOut+sys)
            h.SetTitle(bin+"_"+dsi.label+"_"+basenameOut+sys)
            h.Write()
            histos.append(h)
        return histos

    def writeStatHistos(self, data, basenameIn, basenameOut, bin, sys):
        for key, dsi in data.iteritems():
            h = dsi.getHisto(basenameIn+"_"+bin+sys)
            for i in range(1, h.GetNbinsX()+1):
                name = bin+"_"+dsi.label+"_mcStatBin"+str(i)
                hStatUp = h.Clone(name+"Up")
                hStatUp.SetTitle(name+"Up")
                hStatDown = h.Clone(name+"Down")
                hStatDown.SetTitle(name+"Down")
                content = h.GetBinContent(i)
                error = h.GetBinError(i)
                up = content+error if content+error > 0.0 else 0.0
                down = content-error if content-error > 0.0 else 0.0
                hStatUp.SetBinContent(i, up)
                hStatDown.SetBinContent(i, down)
                hStatUp.Write()
                hStatDown.Write()

    def makePseudoData(self, histos, signalhistos, sgData, basename, bin, sys):
        #make some pseudo_data
        name = bin+"_pseudodata_"+basename+sys
        mynewh = ROOT.TH1D(name, name, histos[0].GetNbinsX(), histos[0].GetBinLowEdge(1), histos[0].GetBinLowEdge(1)+histos[0].GetNbinsX())
        
        #make some pseudo_data with signal
        pseudodataS_histos = []
        for key, dsi in sgData.iteritems():
            sig = dsi.label
            name = bin+"_pseudodataS_"+sig+"_"+basename+sys
            pseudodataS_histos.append( ROOT.TH1D(name, name, histos[0].GetNbinsX(), histos[0].GetBinLowEdge(1), histos[0].GetBinLowEdge(1)+histos[0].GetNbinsX()) )

        for bin in range(histos[0].GetNbinsX()):
            sumdata = 0
            #sumunc = 0
            for h in histos:
                sumdata += h.GetBinContent(bin+1)
                #sumunc += h.GetBinError(bin+1)
            #sumdata = sumdata + sumunc*2*(random.random()-0.5)
            mynewh.SetBinContent(bin+1, int(round(sumdata)))
        
            for i, sigh in enumerate(signalhistos):
                sumsig = sigh.GetBinContent(bin+1)
                #sumsigunc = sigh.GetBinError(bin+1)
                sumdatasigh = sumdata + sumsig #+ sumsigunc*2*(random.random()-0.5)
                pseudodataS_histos[i].SetBinContent(bin+1, int(round(sumdatasigh)))
        
        mynewh.Write()
        for h in pseudodataS_histos:
            h.Write()

class MakeDataCard:
    def __init__(self, outFile, bgData, otData, sgData, basename, mvaBin):
        self.nMVABins = len(mvaBin)
        self.outFile = outFile
        self.bgData = bgData
        self.otData = otData
        self.sgData = sgData
        self.basename = basename
        self.mvaBin = mvaBin
        bgMVAHistos = self.getMVAHists(bgData, mvaBin, basename)
        otMVAHistos = self.getMVAHists(otData, mvaBin, basename)
        sgMVAHistos = self.getMVAHists(sgData, mvaBin, basename)        
        self.dataSets = [(sgMVAHistos, sgData), (bgMVAHistos, bgData), (otMVAHistos, otData)]

    def getMVAHists(self, data, mvaBin, basename):
        histos = []
        for cut in self.mvaBin:
            h = data.getHisto(self.basename+"_"+cut)
            histos.append(h)
        return histos

    def writeDataCard(self, name):
        f = open(name, 'w')

        f.write( "Signal Region 1 Datacard -- signal category" )
        f.write( "\n" )
        f.write( "imax * number of bins\n" )
        f.write( "jmax * number of processes minus 1\n" )
        f.write( "kmax * number of nuisance parameters\n" )
        f.write( "\n" )
        f.write( "-------------------------------------------------------------------------------------------------------------------------------------------\n" )
        f.write( "\n" )
        for i in range(self.nMVABins):
            f.write( "shapes data_obs    sigD{0}    {1} wspace:data_obs_D{0}\n".format(i+1, self.outFile) )
            f.write( "shapes bkg_tt      sigD{0}    {1} wspace:background_tt_D{0}\n".format(i+1, self.outFile) )
            f.write( "shapes bkg_other   sigD{0}    {1} wspace:background_other_D{0}\n".format(i+1, self.outFile) )
            f.write( "shapes signal      sigD{0}    {1} wspace:signal_D{0}\n".format(i+1, self.outFile) )
            f.write( "\n" )
        f.write( "\n" )
        f.write( "-------------------------------------------------------------------------------------------------------------------------------------------\n" )
        bins        = "bin         "
        observation = "observation "
        for i in range(self.nMVABins):
            bins        += " sigD{0} ".format(i+1)
            observation += "  -1   "
        f.write( bins+"\n" )
        f.write( observation+"\n" )
        f.write( "-------------------------------------------------------------------------------------------------------------------------------------------\n" )
        f.write( "# background rate must be taken from _norm param x 1\n" )
        bin      = "bin              "
        process1 = "process          "
        process2 = "process          "
        rate     = "rate             "
        for i in range(self.nMVABins):
            for d in self.dataSets:
                r = 1
                if(d[1].rate):
                    r = d[0][i].Integral()
                bin      += "{0: <16} ".format("sigD"+str(i+1))
                process1 += "{0: <16} ".format(d[1].processName)
                process2 += "{0: <16} ".format(d[1].process)
                rate     += "{0: <16} ".format(r)
        f.write( bin+"\n" )
        f.write( process1+"\n" )
        f.write( process2+"\n" )
        f.write( rate+"\n" )
        f.write( "-------------------------------------------------------------------------------------------------------------------------------------------\n" )
        f.write( "# Normal uncertainties in the signal region\n" )
        lumiSys = "lumi_13TeV      lnN  "
        for i in range(self.nMVABins):
            for d in self.dataSets:
                lumiSys += "{0: <5} ".format(d[1].lumiSys)
        f.write( lumiSys+"\n" )
        f.write( "-------------------------------------------------------------------------------------------------------------------------------------------\n" )
        for i in range(self.nMVABins):
            f.write( "np_tt_D{0} param 0.0 1.0\n".format(i+1) )

if __name__ == "__main__":
    parser = optparse.OptionParser("usage: %prog [options]\n")
    parser.add_option('-d',         dest='basedir',  type='string', default='condor/JECPlots',     help="Path to input rootfiles")
    parser.add_option('--rootFile', dest='rootFile', type='string', default='njets_for_Aron.root', help="Output root file")
    parser.add_option('--dataCard', dest='dataCard', type='string', default='dataCard.txt',        help="Output data card file")
    parser.add_option('-H',         dest='outDir',   type='string', default='',                    help="Can pass in the output directory name")
    options, args = parser.parse_args()

    # Where the root files are stored
    basedir = options.basedir + "/"
    outDir = options.outDir
    if os.path.exists(outDir):
        print "Failed: Output directory %s already exits" % ('"'+outDir+'"')
        exit(0)        
    else:
        os.makedirs(outDir)    
    jettypes = ["pt30_1l"]#, "pt45_0l"]
    mvaBin = [ "D1", "D2","D3","D4"]
    systypes = ["", "_JECUp", "_JECDown", "_JERUp", "_JERDown", "_btgUp", "_btgDown", "_lepUp", "_lepDown"]
    outputfile = ROOT.TFile.Open(outDir + "/" + options.rootFile,"RECREATE")
    outputDataCard = options.dataCard

    # I hadd my ttbar files into TT.root, and I hadd all other backgrounds into BG_noTT.root
    bgData = {
        "TT"    : info.DataSetInfo(basedir=basedir, fileName="TT.root",      label="TT",    processName="bkg_tt",    process="1", rate=False, lumiSys="-"),
        "OTHER" : info.DataSetInfo(basedir=basedir, fileName="BG_noTT.root", label="OTHER", processName="bkg_other", process="2", rate=True,  lumiSys="-"),
    }

    sgData = {
        "SYY_350" : info.DataSetInfo(basedir=basedir, fileName="stealth_stop_350_SYY.root", label="SYY_350", processName="signal", process="0", rate=True, lumiSys="1.05"),
        "SYY_450" : info.DataSetInfo(basedir=basedir, fileName="stealth_stop_450_SYY.root", label="SYY_450", processName="signal", process="0", rate=True, lumiSys="1.05"),
        "SYY_550" : info.DataSetInfo(basedir=basedir, fileName="stealth_stop_550_SYY.root", label="SYY_550", processName="signal", process="0", rate=True, lumiSys="1.05"),
        "SYY_650" : info.DataSetInfo(basedir=basedir, fileName="stealth_stop_650_SYY.root", label="SYY_650", processName="signal", process="0", rate=True, lumiSys="1.05"),
        "SYY_750" : info.DataSetInfo(basedir=basedir, fileName="stealth_stop_750_SYY.root", label="SYY_750", processName="signal", process="0", rate=True, lumiSys="1.05"),
        "SYY_850" : info.DataSetInfo(basedir=basedir, fileName="stealth_stop_850_SYY.root", label="SYY_850", processName="signal", process="0", rate=True, lumiSys="1.05"),
        "RPV_350" : info.DataSetInfo(basedir=basedir, fileName="rpv_stop_350.root",         label="RPV_350", processName="signal", process="0", rate=True, lumiSys="1.05"),
        "RPV_450" : info.DataSetInfo(basedir=basedir, fileName="rpv_stop_450.root",         label="RPV_450", processName="signal", process="0", rate=True, lumiSys="1.05"),
        "RPV_550" : info.DataSetInfo(basedir=basedir, fileName="rpv_stop_550.root",         label="RPV_550", processName="signal", process="0", rate=True, lumiSys="1.05"),
        "RPV_650" : info.DataSetInfo(basedir=basedir, fileName="rpv_stop_650.root",         label="RPV_650", processName="signal", process="0", rate=True, lumiSys="1.05"),
        "RPV_750" : info.DataSetInfo(basedir=basedir, fileName="rpv_stop_750.root",         label="RPV_750", processName="signal", process="0", rate=True, lumiSys="1.05"),
        "RPV_850" : info.DataSetInfo(basedir=basedir, fileName="rpv_stop_850.root",         label="RPV_850", processName="signal", process="0", rate=True, lumiSys="1.05"),
    }

    Data = {
        "data" : info.DataSetInfo(basedir=basedir, fileName="Data.root", label="data", processName="data", process="1", rate=False, lumiSys="-")
    }

    #Write all histos to outputfile
    outputfile.cd()
    wp = WriteNJetPlots()
    for jettype in jettypes:
        for bin in mvaBin:
            for sys in systypes:
                basenameIn  = "h_njetsShifted_" + jettype
                basenameOut = "h_njets_" + jettype                
                histos = wp.writeHistos(bgData, basenameIn, basenameOut, bin, sys)
                signalhistos = wp.writeHistos(sgData, basenameIn, basenameOut, bin, sys)
                if sys == "":
                    wp.writeHistos(Data, basenameIn, basenameOut, bin, sys)
                    wp.writeStatHistos({"OTHER" : bgData["OTHER"]}, basenameIn, basenameOut, bin, sys)
                    wp.writeStatHistos(sgData, basenameIn, basenameOut, bin, sys)
                    wp.makePseudoData(histos, signalhistos, sgData, basenameOut, bin, sys)
    

    #Close outfile
    outputfile.Close()
    
    #Make data card
    for key in sgData:
            md = MakeDataCard(outFile="MVA_ws.root", bgData=bgData["TT"], otData=bgData["OTHER"], sgData=sgData[key], basename="h_njets_pt30_1l", mvaBin=mvaBin)
            md.writeDataCard(outDir+"/"+key+"_"+outputDataCard)
            
    import os
    os.system("cat "+outDir+"/RPV_550_"+outputDataCard)
