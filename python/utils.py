import os
import glob
from ROOT import TFile, TH1F, TH2F, TGraphAsymmErrors, TH1



class Process:
   def __init__(self,alias,directory):
      self.m_alias = alias
      self.m_dir = directory
      xs = {}
      if self.m_alias=='QCD_HT':
         xs['100to200']   = 23700000.0
         xs['200to300']   = 1547000.0
         xs['300to500']   = 322600.0
         xs['500to700']   = 29980.0
         xs['700to1000']  = 6334.0
         xs['1000to1500'] = 1088.0
         xs['1500to2000'] = 99.11
         xs['2000toInf']  = 20.23
      if self.m_alias=='bEnriched':
         xs['100to200']   = 1117000.0
         xs['200to300']   = 80430.0
         xs['300to500']   = 16620.0
         xs['500to700']   = 1487.0
         xs['700to1000']  = 296.5
         xs['1000to1500'] = 46.61
         xs['1500to2000'] = 3.72
         xs['2000toInf']  = 0.6462
      if self.m_alias=='BGenFilter':
         xs['100to200']   = 1282000.0
         xs['200to300']   = 111800.0
         xs['300to500']   = 28070.0
         xs['500to700']   = 3082.0
         xs['700to1000']  = 724.2
         xs['1000to1500'] = 138.2
         xs['1500to2000'] = 13.61
         xs['2000toInf']  = 2.909
      if self.m_alias=='MuEnriched':
         xs['50To80']    = 376600.0
         xs['80To120']   = 88930.0
         xs['120To170']  = 21230.0
         xs['170To300']  = 7055.0
         xs['300To470']  = 619.3
         xs['470To600']  = 59.24
         xs['600To800']  = 18.21
         xs['800to1000'] = 3.275
         xs['1000toInf'] = 1.078
      self.m_xsections = xs
      self.m_bins = xs.keys()
      #root files
      rf = {}
      ngen = {}
      lumi_scale = {}
      if len(xs) > 0:
         keyword = self.m_alias
         ld = glob.glob(directory+'/*'+keyword+'*.root')
         for b in self.m_bins:
            f = [s for s in ld if b in s]
            if keyword == 'QCD_HT':
               f = [s for s in ld if b in s and (not 'BGenFilter' in s)]
            if len(f) == 0:
               print(alias+': No file for bin '+b+'. Skipping!')
               continue
            if len(f) > 1:
               print('There is more than 1 file with bin '+b+'. The first one, '+f[0]+' will be considered')
            rf[b] = f[0]
            # get number of generator+pileup weighted generated events
            hfile = TFile( f[0], 'OLD' )
            h_wf = hfile.Get('workflow')
            ngen[b] = h_wf.GetBinContent(3)
            lumi_scale[b] = xs[b]/ngen[b]
            hfile.Close()
      self.m_rootfiles = rf
      self.m_ngen = ngen
      self.m_lumi_scale = lumi_scale
      # scale to same lumi
      
      
      
   def alias(self):
      return self.m_alias
      
   def crossSections(self):
      return self.m_xsections
      
   def bins(self):
      return self.m_bins
      
   def neventsGenerated(self):
      return self.m_ngen
         
   def rootFiles(self):
      return self.m_rootfiles
      
   def luminosityScale(self,lumi=1000.):
      ls = {}
      for b in self.m_bins:
         ls[b] = self.m_lumi_scale[b]*lumi
      return ls


###################################################

class AnalysisHistograms:
   def __init__(self,process,hdir='/',observables=['pt_jet'],flavours=[]):
      self.m_proc = process
      self.m_histo_dir = hdir
      self.m_flvs = flavours
      self.m_vars = observables
      self.m_objs = []
      self.m_dijets = []
      self.m_histos = {}
      self.m_histos_dijet = {}
      self.m_histos_proc = {}
      self.m_histos_dijet_proc = {}
      
# -------------------------------------------------------      

   # read histograms from the files   
   def readTH1(self,lumi_scale=True,combined=False):
      h1 = {}
      pbins = []
      nobjs = 0
      objs = []
      h_types = self.m_vars
      flv_all = self.flavours()
      flv_all.append('all')
      comb = 1
      if combined:
         comb = 2
         # TO DO: flavour combinations
         if len(flv_all) > 1:
            flv_all = ['all']
      for pbin, pfile in self.m_proc.rootFiles().items():  # need to think better these loops to avoid more loops below
         ls = self.m_proc.luminosityScale()[pbin]
         h1[pbin] = {}
         pbins.append(pbin)
         tfile = TFile( pfile, 'OLD' )
         h_dir = tfile.GetDirectory(self.m_histo_dir)
         h_dir_keys = h_dir.GetListOfKeys()
         histo_names = {}
         for typ in h_types:
            histo_names = [x.GetName() for x in h_dir_keys if typ in x.GetName() and x.GetName().count('_')==1 and len(x.GetName())==(len(typ)+comb)]
            h1[pbin][typ] = {}
            for h_name in histo_names:
               obj = int(h_name[-1*comb:])
               h1[pbin][typ][obj] = {}
               for flv in flv_all:
                  hname_flv = h_name
                  if flv != 'all':
                     hname_flv = '{}_{}'.format(h_name,flv)
                  h1[pbin][typ][obj][flv] = h_dir.Get(hname_flv)
                  h1[pbin][typ][obj][flv].SetName('{}_{}'.format(hname_flv,pbin))
                  h1[pbin][typ][obj][flv].SetDirectory(0) # Important!
                  if lumi_scale:
                     h1[pbin][typ][obj][flv].Scale(ls)
         tfile.Close()
                 
      self.m_histos_proc = h1
         
#      self.m_jets = objs
      
      # SUM OF ALL (HT) BINS OF THE PROCESS
      h1_add_proc = {}
      if lumi_scale:
         for typ in h_types:
            h1_add_proc[typ] = {}
            for obj in h1[pbin][typ].keys():
               if nobjs == 0:
                  objs.append(obj)
               h1_add_proc[typ][obj] = {}
               for flv in flv_all:
                  h1_add_proc[typ][obj][flv] = {}
                  for b,pbin in enumerate(pbins):
                     h_name = h1[pbin][typ][obj][flv].GetName().replace('_'+pbin,'')
                     if b == 0:
                        h1_add_proc[typ][obj][flv] = h1[pbin][typ][obj][flv].Clone(h_name)
                     else:
                        h1_add_proc[typ][obj][flv].Add(h1[pbin][typ][obj][flv])
                     h1_add_proc[typ][obj][flv].SetDirectory(0) # Important!
            if nobjs == 0:
               nobjs = len(objs)
   
      self.m_histos = h1_add_proc
      self.m_objs = objs
      
# -------------------------------------------------------      
            
   def histogram_directory(self):
      return self.m_histo_dir
      
   def flavours(self):
      return self.m_flvs
      
   def histograms(self,perProcessBin=False):
      if perProcessBin:
         return self.m_histos_proc
      return self.m_histos
      
   def objects(self):
      return self.m_objs
      
      
############################################################################      
