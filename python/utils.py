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
         xs['100to200']   = 23700000.0
         xs['200to300']   = 1547000.0
         xs['300to500']   = 322600.0
         xs['500to700']   = 29980.0
         xs['700to1000']  = 6334.0
         xs['1000to1500'] = 1088.0
         xs['1500to2000'] = 99.11
         xs['2000toInf']  = 20.23
      if self.m_alias=='BGenFilter':
         xs['100to200']   = 23700000.0
         xs['200to300']   = 1547000.0
         xs['300to500']   = 322600.0
         xs['500to700']   = 29980.0
         xs['700to1000']  = 6334.0
         xs['1000to1500'] = 1088.0
         xs['1500to2000'] = 99.11
         xs['2000toInf']  = 20.23
      if self.m_alias=='QCD_Mu':
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
         if self.m_alias == 'QCD_Mu':
            keyword = 'MuEnrichedPt5'
         ld = glob.glob(directory+'/*'+keyword+'*.root')
         for b in self.m_bins:
            f = [s for s in ld if b in s]
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

