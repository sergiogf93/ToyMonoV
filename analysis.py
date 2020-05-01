from ROOT import gROOT, TFile, TH1, TH2D, TCanvas

from optparse import OptionParser
from functions import *   # Here the functions like printError are defined

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

#######################################
parser = OptionParser()

parser.add_option('-i', '--input-file',type='string',help='Path to the input root file', dest = "input_file", default="")
parser.add_option('-t', '--tree-name',type='string',help='Name of the tree in the input file', dest = "tree_name", default="ntuple")
parser.add_option("--do-CNN", action="store_true", dest="doCNN", help="Use the CNN in the selection" , default=False)

(options, args) = parser.parse_args()
#######################################
# -- unit conversion units, and image parameters
GeV = 1000
DegToRad = np.pi / 180.
step = 0.05
N = int(1/step)

#######################################

def main(options):
    #Making sure input file is provided
    if options.input_file == "":
        printError("No input file provided. Run with the -i option and point to a root file")
        exit()
    os.system('mkdir -p Plots')
    runANALYSIS(options)


def runANALYSIS(options):
    rootFile = TFile(options.input_file)
    tree = gROOT.FindObject(options.tree_name)

    met_for_histogram = []  # Array for making a histogram of MET before the cuts

    yield_VTagger= {"raw_events":0,     # Number of raw events in the signal selection (Merged 0b HP)
                    "sum_weights":0}    # Sum of weights in the signal selection (Merged 0b HP)

    yield_CNN= {"raw_events":0,
                "sum_weights":0}

    n_events = tree.GetEntries()
    for i_event in range(n_events):
        tree.GetEntry(i_event) # Update the entry with the properties of the event
        if i_event % 5000 == 0 or i_event == n_events - 1: print("Processing event {} ({:.2f} %)".format(i_event, 100*(i_event+1)/n_events))

        # Filling arrays before selection
        met_for_histogram.append(tree.met_tst_et / 1000.)

        # Start selection:
        if tree.met_fire == 0: continue                                     # MET triggers
        if tree.n_el_baseline != 0 or tree.n_mu_baseline != 0: continue     # Lepton veto
        if tree.n_TCCJet == 0: continue                                     # Ask for one Large-R Jet
        if tree.dPhiTCCJetMet <= 120*DegToRad: continue                     # Delta phi between large jet and MET has to be > 120º
        if tree.n_trackTCCSeparatedBjet != 0: continue                      # No b-jets out of the Large jet
        if tree.n_trackTCCAssociatedBjet != 0: continue                     # Doing the 0 b-jet selection
        if tree.met_tst_et <= 250*GeV: continue                             # MET > 250 GeV

        # V-Tagger
        D2_isWJet = tree.TCCJet_passD2_W50[0] == 1
        D2_isZJet = tree.TCCJet_passD2_Z50[0] == 1
        Mass_isWJet = tree.TCCJet_passMass_W50[0] == 1
        Mass_isZJet = tree.TCCJet_passMass_Z50[0] == 1
        if D2_isWJet and D2_isZJet and Mass_isWJet and Mass_isZJet:         # Pass V-Tagger
            # If you reach here it means that the event passes the selection! (merged 0b HP)
            yield_VTagger["raw_events"] += 1
            yield_VTagger["sum_weights"] += tree.weight

        # CNN
        # Implement applying the CNN to the image per event
        if not options.doCNN: continue

        large_jet_pt = tree.TCCJet_pt[0] / GeV # pT of the leading large R jet in GeV
        th2_image = GetTH2DImage(tree.TCCJet_eta[0], tree.TCCJet_phi[0], tree.TCCJet_lead_clusters_pt, tree.TCCJet_lead_clusters_eta, tree.TCCJet_lead_clusters_phi)            
        image_as_array = ConvertToArray(th2_image)

    #out of looping tree
    makeHistogram(met_for_histogram, "MET", "MET [GeV]", 100, 0, 1200)
    print("-"*50)
    print("Yield = {:.2f}".format(36100 * yield_VTagger["sum_weights"])) # Extrapolate to 36.1 fb-1

def GetTH2DImage(TCCJet_eta, TCCJet_phi, clusters_pt, clusters_eta, clusters_phi):
    hist2D = TH2D("", "", 2*N, -1, 1, 2*N, -1, 1)
    for i_cluster in range(0,len(clusters_pt)):
        d_eta = clusters_eta[i_cluster] - TCCJet_eta
        d_phi = clusters_phi[i_cluster] - TCCJet_phi
        if d_phi > np.pi: d_phi = -2*np.pi + d_phi
        if d_phi < -1*np.pi: d_phi = 2*np.pi - d_phi

        bin = Get2DBin(hist2D, d_eta, d_phi)
        hist2D.SetBinContent(bin, hist2D.GetBinContent(bin) + clusters_pt[i_cluster]/GeV)
    
    # hist2D.Draw("colz")  # This Draws per call the image
    return hist2D

def ConvertToArray(hist):
    array = []
    for c_phi in np.arange(1-step/2., -1, -1*step):
        row = []
        for c_eta in np.arange(-1 + step/2., 1, step):
            bin = Get2DBin(hist, c_eta, c_phi)
            row.append(hist.GetBinContent(bin))
        array.append(row)
    return np.array(array)


def makeHistogram(array,name,xLabel,nbins,xmin,xmax):
    plt.hist(array,bins=nbins,range=[xmin,xmax])
    plt.yscale('log')
    plt.xlabel(xLabel)
    plt.savefig("Plots/" + name + ".png")
    print("Saving image: Plots/" + name + ".png")
    plt.clf()


if __name__ == '__main__':
    main(options)