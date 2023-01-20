#!/usr/bin/env python

## Offline btag weighted samples compared to btag selections
## use with bEnriched and BGenFilter samples

from __future__ import print_function

import sys,os
import array
from argparse import ArgumentParser
from argparse import HelpFormatter

from Analysis.Tools.utils import Process, AnalysisHistograms

from ROOT import TFile,TCanvas, TMultiGraph,TLine, TRatioPlot, gStyle, TH1, Double, gROOT
from ROOT import kRed, kBlue, kBlack, kMagenta, kGreen, kCyan
from rootpy.interactive import wait

import Analysis.Tools.CMS_lumi as CMS_lumi
import Analysis.Tools.tdrstyle as tdrstyle

TH1.SetDefaultSumw2()

#set the tdr style (why some setting don't work? using gStyle??? - not all settings work w/ gStyle???)
tdrstyle.setTDRStyle()

# More styles

#change the CMS_lumi variables (see CMS_lumi.py)
CMS_lumi.lumi_7TeV = "4.8 fb^{-1}"
CMS_lumi.lumi_8TeV = "18.3 fb^{-1}"
CMS_lumi.writeExtraText = 1
CMS_lumi.extraText = "Simulation"
CMS_lumi.lumi_sqrtS = "13 TeV" # used with iPeriod = 0, e.g. for simulation-only plots (default is an empty string)

iPos = 0
if( iPos==0 ): CMS_lumi.relPosX = 0.12

# 
# Simple example of macro: plot with CMS name and lumi text
#  (this script does not pretend to work in all configurations)
# iPeriod = 1*(0/1 7 TeV) + 2*(0/1 8 TeV)  + 4*(0/1 13 TeV) 
# For instance: 
#               iPeriod = 3 means: 7 TeV + 8 TeV
#               iPeriod = 7 means: 7 TeV + 8 TeV + 13 TeV 
#               iPeriod = 0 means: free form (uses lumi_sqrtS)
# Initiated by: Gautier Hamel de Monchenault (Saclay)
# Translated in Python by: Joshua Hardenbrook (Princeton)
# Updated by:   Dinko Ferencek (Rutgers)
#

iPeriod = 0

xbins = []
for i in xrange(0,120,10):
   xbins.append(i)
for i in xrange(120,500,20):
   xbins.append(i)
for i in xrange(500,800,30):
   xbins.append(i)
for i in xrange(800,1200,50):
   xbins.append(i)
for i in xrange(1200,2000,100):
   xbins.append(i)               
for i in xrange(2000,2600,150):
   xbins.append(i)
for i in xrange(2600,3001,400):
   xbins.append(i)

xbins = []   
for i in xrange(0,3000,50):
   xbins.append(i)
   
axbins = array.array('d', xbins)

# -------------------------------------------------------      

def make_plots(histos1,histos2=None,legend1=None,legend2=None,legend_title=None,combined=False,ratio=None,ratio_range=None,flavour=False):

   global output_name
   global x_range, y_range
   global iPos
   global H_ref,W_ref,H,W,T,B,L,R
   global textSize
   
   flavours = ['b','c','udsg','bb','cc']
   canvas = {}
   ratioplots  = {}
   legends = {}
   histos1_rb = {}
   histos2_rb = {}
   
   xmin,xmax = get_range(x_range)
   ymin = -1
   ymax = -1
   if y_range:
      ymin,ymax = get_range(y_range)

   if ratio:
      H_ref = 840
      textSize = 0.036
   H = H_ref
   W = W_ref
   # references for T, B, L, R
   T = 0.080*H_ref
   B = 0.150*H_ref 
   L = 0.135*W_ref
   R = 0.065*W_ref

   if not histos1 and not histos2:
      print("No histograms to plot")
      quit()
   
   if not histos2:
      ratio = False
      
   for var, value1 in histos1.iteritems():
      canvas[var]  = {}
      ratioplots[var]   = {}
      legends[var] = {}
      histos1_rb[var] = {}
      histos2_rb[var] = {}
      for obj, value2 in value1.iteritems():
         if obj == 13 or obj == 23:
            continue
         canvas[var][obj]  = {}
         ratioplots[var][obj]   = {}
         legends[var][obj] = {}
         histos1_rb[var][obj] = {}
         histos2_rb[var][obj] = {}
         for flv,th1 in value2.iteritems():
            if not flavour and flv != 'all':
               continue
            c_title = 'c_{}{}_{}'.format(var,obj,flv)
            canvas[var][obj][flv]  = TCanvas( c_title, c_title, 2000, 10, W, H )
            c = canvas[var][obj][flv]
            ratioplots[var][obj][flv] = None
            hmax_sf = 1.5
            if var == 'pt_jet' or var == 'm_jet':
               canvas[var][obj][flv].SetLogy()
               hmax_sf = 2.
            if var == 'm_jet':
               histos1_rb[var][obj][flv] = histos1[var][obj][flv].Rebin(len(axbins)-1,'{}{}'.format(var,obj),axbins)
               h1 = histos1_rb[var][obj][flv]
#               h1 = histos1[var][obj][flv]
#               h1.Rebin(25)
            else:
               h1 = histos1[var][obj][flv]
               h1.Rebin(10)
            h1.SetName('{}_h1'.format(h1.GetName()))
            h1max = h1.GetMaximum()
            h2 = None
            prep_histogram(h1,title=legend1)
            if histos2:
               if var == 'm_jet':
                  histos2_rb[var][obj][flv] = histos2[var][obj][flv].Rebin(len(axbins)-1,'{}{}'.format(var,obj),axbins)
                  h2 = histos2_rb[var][obj][flv]
#                  h2 = histos2[var][obj][flv]
#                  h2.Rebin(25)
               else:
                  h2 = histos2[var][obj][flv]
                  h2.Rebin(10)
               h2.SetName('{}_h2'.format(h2.GetName()))
               h2max = h2.GetMaximum()
               if y_range:
                  h1.SetMaximum(ymax)
                  h1.SetMinimum(ymin)
               else:
                  if h1max < h2max:
                     h1.SetMaximum(h2max*hmax_sf)
               prep_histogram(h2,color=kRed,title=legend2)
            # Simple plotting of histograms
            if not ratio:
               h1.GetXaxis().SetRangeUser(xmin,xmax)
               h1.Draw()
               if h2:
                  h2.Draw('same')
               legends[var][obj][flv] = prep_canvas(c,legend_title=legend_title)
               #draw the lumi text on the canvas
               CMS_lumi.CMS_lumi(c, iPeriod, iPos)
               c.cd()
               c.Update()
               c.RedrawAxis()
            else:
               ratioplots[var][obj][flv] = TRatioPlot(h1,h2,"divsym")
               rp = ratioplots[var][obj][flv]
               prep_ratioplots(rp)
               rp.Draw()
               legends[var][obj][flv] = mod_ratioplots(rp,c,title=ratio,legend_title=legend_title)
               ylmin,ylmax = get_range(ratio_range)
               rp.GetLowerRefYaxis().SetRangeUser(ylmin,ylmax)
               rp.GetUpperRefXaxis().SetRangeUser(xmin,xmax)
               rp.GetLowerRefXaxis().SetRangeUser(xmin,xmax)
#               c.Update()
#               c.Draw()
               
               #draw the lumi text on the canvas
               CMS_lumi.CMS_lumi(rp.GetUpperPad(), iPeriod, iPos)
               c.cd()
               c.Update()
               c.RedrawAxis()
            
            
            c.SaveAs('{}_{}_{}_{}.png'.format(output_name,var,obj,flv))
               
            if 'm_jet12_h1' in h1.GetName(): 
               output = TFile.Open('{}_m_jet12.root'.format(output_name),'recreate')
               if ratio:
                  m12_graph = rp.GetLowerRefGraph()
                  m12_graph.Write()
               h1.Write()
               h2.Write()
               output.Close()
            
            
#   wait()
            
# ------------------------------------------------------- 
 
def get_range(r):
   r=r.replace(' ','')
   rr = r.split(',')
   rmin = 0.
   rmax = 2.
   if len(rr) > 1:
      rmin = float(rr[0])
      rmax = float(rr[1])
   if rmin > rmax:
      rmin = float(rr[1])
      rmax = float(rr[0])
   return rmin,rmax
         
# ------------------------------------------------------- 
 
def prep_ratioplots(rp,title=None):
   rp.SetH1DrawOpt("e")
   rp.SetH2DrawOpt("e")
   rp.GetLowYaxis().SetNdivisions(6,5,0)

# ------------------------------------------------------- 
 
def mod_ratioplots(rp, canvas,title=None,legend_title=None):
   canvas.SetFillColor(0)
   canvas.SetBorderMode(0)
   canvas.SetFrameFillStyle(0)
   canvas.SetFrameBorderMode(0)
   canvas.SetLeftMargin( L/W )
   canvas.SetRightMargin( R/W )
   canvas.SetTopMargin( T/H )
   canvas.SetBottomMargin( B/H )
#   canvas.SetTickx(0)
#   canvas.SetTicky(0)
   
   rp.GetLowerRefYaxis().SetTitle(title)
   rp.GetLowerRefYaxis().SetTitleOffset(1.95)
   rp.GetUpperRefYaxis().SetTitleOffset(1.95)
   rp.GetLowerRefXaxis().SetTitleOffset(1.2)
   rp.GetLowerRefGraph().SetLineColor(kBlack)
   rp.GetLowerRefGraph().SetLineWidth(2)
   rp.GetLowerRefGraph().SetMarkerStyle(20)
   rp.GetLowerRefGraph().SetMarkerColor(kBlack)
   rp.RangeAxisChanged()
   
   rp.SetLeftMargin(L/W)
   rp.SetRightMargin(R/W)
   rp.SetUpTopMargin(T/H)
#    rp.SetLowBottomMargin(0.35)
# #   rp.SetUpBottomMargin(0.4)
#    rp.SetLowTopMargin(0.01)

   
   rp.GetLowerRefXaxis().SetLabelFont(42)
   rp.GetLowerRefXaxis().SetTitleFont(42)
   rp.GetLowerRefYaxis().SetLabelFont(42)
   rp.GetLowerRefYaxis().SetTitleFont(42)
   rp.GetUpperRefYaxis().SetLabelFont(42)
   rp.GetUpperRefYaxis().SetTitleFont(42)
   rp.GetLowerRefXaxis().SetLabelSize(textSize)
   rp.GetLowerRefXaxis().SetTitleSize(textSize)
   rp.GetLowerRefYaxis().SetLabelSize(textSize)
   rp.GetLowerRefYaxis().SetTitleSize(textSize)
   rp.GetUpperRefYaxis().SetLabelSize(textSize)
   rp.GetUpperRefYaxis().SetTitleSize(textSize)
   
#   rp.SetSplitFraction(0.5)
#   rp.SetSeparationMargin(0.)
   # SetSplitFraction does not work!? Brute force...
   x1= rp.GetLowerPad().GetXlowNDC()
   x2= x1+rp.GetLowerPad().GetWNDC()
   y1= rp.GetLowerPad().GetYlowNDC()
   y2= y1+rp.GetLowerPad().GetHNDC()
   rp.GetLowerPad().SetPad(x1,y1,x2,y2*1.2)
   
   ## LEGEND
   if legend_title:
      legend = rp.GetUpperPad().BuildLegend(0.50,0.70,0.935,0.895,legend_title)
   else:
      legend = rp.GetUpperPad().BuildLegend(0.50,0.76,0.935,0.895)
   legend.SetTextFont(42)
   legend.SetTextSize(0.04)
   
   return legend
   

   
# ------------------------------------------------------- 
 
def prep_canvas(canvas,legend_title=None):
   canvas.SetFillColor(0)
   canvas.SetBorderMode(0)
   canvas.SetFrameFillStyle(0)
   canvas.SetFrameBorderMode(0)
   canvas.SetLeftMargin( L/W )
   canvas.SetRightMargin( R/W )
   canvas.SetTopMargin( T/H )
   canvas.SetBottomMargin( B/H )
#   canvas.SetTickx(0)
#   canvas.SetTicky(0)
   
   if legend_title:
      legend = canvas.BuildLegend(0.50,0.70,0.935,0.895,legend_title)
   else:
      legend = canvas.BuildLegend(0.50,0.76,0.935,0.895)
   legend.SetTextFont(42)
   legend.SetTextSize(0.035)
   
   return legend


# ------------------------------------------------------- 
 
    
def prep_histogram(h,title=None,color=kBlack):
   xmin,xmax = get_range(x_range)
   h.GetXaxis().SetRangeUser(xmin,xmax)

   h.SetLineColor(color)
   h.SetMarkerColor(color)
   h.SetMarkerStyle(20)
   h.SetLineWidth(2)
   
   xAxis = h.GetXaxis()
   xmin,xmax = get_range(x_range)
   xAxis.SetRangeUser(xmin,xmax)
   xAxis.SetNdivisions(6,5,0)
   xAxis.SetTitleOffset(1.2)
   xAxis.SetLabelFont(42)
   xAxis.SetTitleFont(42)
   xAxis.SetLabelSize(textSize)
   xAxis.SetTitleSize(textSize)

   yAxis = h.GetYaxis()
   yAxis.SetNdivisions(6,5,0)
   yAxis.SetTitle('Events')
   yAxis.SetTitleOffset(1.4)
   yAxis.SetLabelFont(42)
   yAxis.SetTitleFont(42)
   yAxis.SetLabelSize(textSize)
   yAxis.SetTitleSize(textSize)
   
   if title:
      h.SetTitle(title)

# ------------------------------------------------------- 
     
def histograms(path,process,hdir,observables,combined=False):
   # combined means objects of the observable are combined, eg dijets m_jet12 is the mass of jet1 and jet2 combined
   process_list = [process]
   if process == 'bEnriched':
      process_list = ['bEnriched','BGenFilter']
   proc_histos = []
   nobjs = []
   for i,p in enumerate(process_list):
      proc = Process(p,path)
      ah = AnalysisHistograms(process=proc,hdir=hdir,observables=observables)
      ah.readTH1(lumi_scale=True,combined=combined)
      proc_histos.append(ah.histograms())
      nobjs.append(len(ah.objects()))
      
   
   histos = proc_histos[0]
   if process == 'bEnriched':
      if nobjs[0] != nobjs[1]:
         print('-e-: bEnriched and BGenFilter samples do not have the same number of jets')
         quit()
         
      for i in range(1,len(proc_histos)):
         for var, value1 in histos.iteritems():
            for obj, value2 in value1.iteritems():
               for flv,th1 in value2.iteritems():
                  histos[var][obj][flv].Add(proc_histos[i][var][obj][flv])
                  
   return histos
      
# -------------------------------------------------------      
      

def main():
   global output_name
   global x_range, y_range
   global textSize
   
   global tdrStyle
   
   global H_ref,W_ref,H,W,T,B,L,R
   
   
   output_name = 'qcd_plots_output'
   H_ref = 600
   W_ref = 800
   
   textSize = 0.05
      
   obs1 = ['pt_jet','eta_jet','phi_jet']
   obs2 = ['m_jet']
   flavours = ['b','c','udsg','bb','cc']
   # From the command line
   # (maybe in the future use ConfigParser from a config file)
   parser = ArgumentParser(prog='qcd_plots.py', formatter_class=lambda prog: HelpFormatter(prog,indent_increment=6,max_help_position=80,width=280), description='QCD plots',add_help=True)
   parser.add_argument("--processes"       , dest="processes"                        , help="processes to be used (NB: bEnriched includes BGenFilter) - MANDATORY")
   parser.add_argument("--paths"           , dest="paths"                            , help="paths of the files w/ histograms (one for each process) - MANDATORY")
   parser.add_argument("--hdirs"           , dest="hdirs"                            , help="histograms directories (one for each process) - MANDATORY")
   parser.add_argument("--single_vars"     , dest="svars"      , default=obs1        , help="observables single object")
   parser.add_argument("--combined_vars"   , dest="cvars"      , default=obs2        , help="observables combined objects")
   parser.add_argument("--legends"         , dest="legs"                             , help="legends")
   parser.add_argument("--legend_title"    , dest="leg_title"  , default="QCD plots" , help="legend title")
   parser.add_argument("--ratio"           , dest="ratio"                            , help="ratio plots title")   
   parser.add_argument("--ratio_range"     , dest="ratio_range", default="0,2"       , help="ratio plots range")   
   parser.add_argument("--x_range"         , dest="x_range"    , default="260,2000"  , help="x-axis range")   
   parser.add_argument("--y_range"         , dest="y_range"                          , help="x-axis range")   
   parser.add_argument("--output"          , dest="output"                           , help="name for output files")   
   args = parser.parse_args()

   # Check and manipulate arguments   
   if not args.processes or not args.paths or not args.hdirs:
      parser.print_help()
      quit()
   if not args.svars and not args.cvars:
      print('No variables to be plotted')
      quit()
   processes = args.processes.split(',')
   files_paths = ['.']
   if args.paths:
      files_paths = args.paths.split(',')
   hdirs = ['/']
   if args.hdirs:
      hdirs = args.hdirs.split(',')
   
   # GLOBALS
   if args.output:
      output_name = args.output
   x_range = args.x_range
   y_range = args.y_range
      
   legs = []   
   histos = []
   for i in range(len(files_paths)):
      legs.append(None)
      histos.append({})
      if args.svars:
         histos[i]['single']   = histograms(path=files_paths[i],process=processes[i],hdir=hdirs[i],observables=args.svars,combined=False)
      if args.cvars:
         histos[i]['combined'] = histograms(path=files_paths[i],process=processes[i],hdir=hdirs[i],observables=args.cvars,combined=True)
      
   if args.legs:
      legs = args.legs.split(',')
      
   # SINGLE HISTOGRAM
   if len(files_paths) == 1:
      if args.svars:
         make_plots(histos1=histos[0]['single']  , legend1=legs[0], legend_title=args.leg_title, combined=False, ratio=args.ratio, flavour=False)
      if args.cvars:
         make_plots(histos1=histos[0]['combined'], legend1=legs[0], legend_title=args.leg_title, combined=True , ratio=args.ratio, flavour=False)
         
   # COMPARISON OF 2 HISTOGRAMS
   if len(files_paths) == 2:
      if args.svars:
         make_plots(histos1=histos[0]['single']  , histos2=histos[1]['single']  ,legend1=legs[0], legend2=legs[1], legend_title=args.leg_title, combined=False, ratio=args.ratio, ratio_range=args.ratio_range, flavour=False)
      if args.cvars:
         make_plots(histos1=histos[0]['combined'], histos2=histos[1]['combined'],legend1=legs[0], legend2=legs[1], legend_title=args.leg_title, combined=True , ratio=args.ratio, ratio_range=args.ratio_range, flavour=False)
   

   

##################

main()


# qcd_plots.py \
# --processes=bEnriched,bEnriched \
# --paths=results/10_btagsel/mssmhbb_fh_2018_cr,results/10_btagweight/mssmhbb_fh_2018_cr \
# --hdirs=final_selection,final_selection \
# --legends="btag selection","btag weight" \
# --legend_title="2018 FH - CR" \
# --single_vars="" \
# --ratio="no veto/veto" \
# --ratio_range="0.51,1.99" \
# --x_range="200,2000" \
# --output=mssmhbb_fh_2018_cr_bsel_x_bw