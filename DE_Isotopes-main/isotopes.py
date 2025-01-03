#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 13 00:16:58 2024

@author: eunicekoo
"""


############### Imports ###############

from brainpy import isotopic_variants
import re
from pyteomics import mass
import pandas as pd
import numpy as np
from pyteomics import mzml
import os
import matplotlib.pyplot as plt
from tqdm import tqdm
plt.rcParams['figure.dpi'] = 500


############### Initializations ###############

# all unimod modifications are stored here
unimods = mass.Unimod()

class Spectrum:

    def __init__(self,scan=None):
        self.id = None
        self.level=None
        self.RT=None
        self.mz=None
        self.intens=None
        self.collision_energy = None
        self.TIC=None

        if scan:
            self.get_vals(scan)

    def get_vals(self,scan):
        # extract values from mzml spectrum
        self.id = scan["id"]
        self.scan_num = int(re.search("scan=(\d+)",self.id)[1])
        self.level=scan["ms level"]
        self.RT = scan['scanList']['scan'][0]["scan start time"]
        self.mz = scan["m/z array"]
        self.intens = scan["intensity array"]
        if self.level==2:
            self.collision_energy = scan["precursorList"]["precursor"][0]["activation"]
            self.scanwindow = [scan["scanList"]["scan"][0]["scanWindowList"]["scanWindow"][0][i] for i in ["scan window lower limit","scan window upper limit"]]
            isolationWindow = scan["precursorList"]["precursor"][0]["isolationWindow"]
            self.prec_mz = isolationWindow["isolation window target m/z"]
        self.TIC = scan["total ion current"]

    def peak_list(self):
        return(np.array([self.mz,self.intens]))
    
class SpectrumFile:
    

    def __init__(self,mzml_file=None):

        self.filename = None
        self.ms1scans = {}
        self.ms2scans = {}
        
        if mzml_file:
            self.load_spectra(mzml_file)

    def load_spectra(self, mzml_file):
        self.filename = mzml_file
        with mzml.MzML(mzml_file) as reader:
            for scan in reader:
                spectrum = Spectrum(scan)
                if spectrum.level == 1:
                    self.ms1scans[spectrum.scan_num] = spectrum
                elif spectrum.level == 2:
                    self.ms2scans[spectrum.scan_num] = spectrum

def loadSpectra(mzml_file):
    print("Loading Spectra",end=" ")
    python_spec_file = mzml_file+"_pythonspec"
    if not os.path.exists(python_spec_file):
        print("... from file")
        spectra = SpectrumFile(mzml_file)
            
    print(f"Loaded {len(spectra.ms1scans)} MS1 spectra")
    print(f"Loaded {len(spectra.ms2scans)} MS2 spectra")
    print("finished")
    
    return spectra


############### Functions ###############

# split up the fragment name (b/y)(-loss)(frag index)_charge
def split_frag_name(ion_type):
    
    frag_name, frag_z = ion_type.split("_")
    loss_check = frag_name.split("-")
    loss = ""
    
    if len(loss_check)>1:
        frag_name,loss = loss_check
        
    frag_type = frag_name[0]
    frag_idx = int(frag_name[1:])
    
    return frag_type, frag_idx, loss, frag_z


# get the compostion of the fragment
def get_seq_comp(split_seq, ion_type):
    
    stripped_seq = "".join([i[0] for i in split_seq]) # assumes AA comes first before mods
    
    mods = [int(j) for i in split_seq for j in re.findall("\([A-z]+\:(\d+)\)",i) if len(i)>1]
    seq_comp = mass.Composition(sequence=stripped_seq, ion_type=ion_type)
    
    for unimod_idx in mods:
        seq_comp += unimods.by_id(unimod_idx)["composition"]
        
    return seq_comp


# first get the AA sequence and modifications of the fragment
def fragment_seq(peptide, ion_type):
    
    peptide = "".join(peptide)
    split_peptide = re.findall("([A-Z](?:\(.*?\))?)", peptide)
    
    mods = re.finditer("\((.*?)\)", peptide)

    stripped_peptide = re.sub("\(.*?\)","",peptide)
    
    frag_type,frag_idx,loss,frag_z = split_frag_name(ion_type)
    
    assert int(frag_idx)<len(stripped_peptide)
    if frag_type in 'abc':
        seq = split_peptide[:int(frag_idx)]
    elif frag_type in 'xyz':
        seq = split_peptide[-int(frag_idx):]
    else:
        raise(ValueError("Invalid ion type"))
        
    return seq, [frag_type,frag_idx,loss,frag_z]


def frag_isotope(frag, seq):
    
    split_frag_seq, frag_info = fragment_seq(seq,frag)
    loss = "-" + frag_info[2] if frag_info[2] else frag_info[2]
    ion_type = frag_info[0] + loss
    frag_comp = get_seq_comp(split_frag_seq, ion_type)
    
    isotopes = isotopic_variants(frag_comp, npeaks=10, charge=int(frag_info[3]))
    mono_iso_peak = isotopes[0]
    
    return isotopes


def ppm(ppm, mass):
    
    return (mass * ppm) / 1000000


def get_intensities(theoretical_mz_list, observed_mz_list, intensity_list, ppm_val):
    
    result_dict = {}
    
    for theoretical_mz in theoretical_mz_list:
        ppm_calculated_dict = {}
        for index, observed_mz in enumerate(observed_mz_list):
            if abs(theoretical_mz - observed_mz) <= ppm(ppm_val, observed_mz):
                ppm_calculated_dict[observed_mz] = [theoretical_mz, intensity_list[index], (((theoretical_mz - observed_mz) / observed_mz) * 1000000)]
                result_dict[theoretical_mz] = intensity_list[index]
                
            
        if len(ppm_calculated_dict) > 1:
            print('does this ever get reached')
            
            items = list(ppm_calculated_dict.items())
            min_key, min_value = min(items, key=lambda x: abs(x[1][1]))
            result_dict[min_value[0]] = min_value[1]
            
    return result_dict


############### Analysis ###############

d8_fixed_msms = pd.read_csv("/Volumes/Lab/SS/Explorii/2024-10-31_de-aging/d8/fixed-dpMQ/combined/txt/msms.txt", sep='\t', low_memory=False) 
d8_240822_msms = d8_fixed_msms[d8_fixed_msms['Raw file'] == '2024-08-22_DE-d8_50ng_DDA'] 
d8_241030_msms = d8_fixed_msms[d8_fixed_msms['Raw file'] == '2024-10-30_DE-d8_50ng_DDA'] 

d8_240822_msms['precursor'] = d8_240822_msms['Sequence'] + d8_240822_msms['Charge'].astype(str)
d8_241030_msms['precursor'] = d8_241030_msms['Sequence'] + d8_241030_msms['Charge'].astype(str)

d8_240822_msms_grouped = d8_240822_msms.loc[d8_240822_msms.groupby('precursor')['Score'].idxmax()]
d8_241030_msms_grouped = d8_241030_msms.loc[d8_241030_msms.groupby('precursor')['Score'].idxmax()]

common_precursors = set(d8_240822_msms_grouped['precursor']) & set(d8_241030_msms_grouped['precursor'])

d8_240822_msms_filtered = d8_240822_msms_grouped[d8_240822_msms_grouped['precursor'].isin(common_precursors)]
d8_241030_msms_filtered = d8_241030_msms_grouped[d8_241030_msms_grouped['precursor'].isin(common_precursors)]

d8_240822_mzml = loadSpectra("/Volumes/Lab/SS/Explorii/2024-10-31_de-aging/d8/fixed-dpMQ/2024-08-22_DE-d8_50ng_DDA.mzML")
d8_240822_ms1 = d8_240822_mzml.ms1scans
d8_240822_ms2 = d8_240822_mzml.ms2scans

d8_241030_mzml = loadSpectra("/Volumes/Lab/SS/Explorii/2024-10-31_de-aging/d8/fixed-dpMQ/2024-10-30_DE-d8_50ng_DDA.mzML")
d8_241030_ms1 = d8_241030_mzml.ms1scans
d8_241030_ms2 = d8_241030_mzml.ms2scans

d8_241030_ms1_spec_idxs = np.array(list(d8_241030_ms1.keys()))
d8_241030_ms2_spec_idxs = np.array(list(d8_241030_ms2.keys()))

d8_240822_ms1_spec_idxs = np.array(list(d8_240822_ms1.keys()))
d8_240822_ms2_spec_idxs = np.array(list(d8_240822_ms2.keys()))


############### Iterating through MSMS File ###############

for (index1, row1), (index2, row2) in zip(d8_240822_msms_filtered.iloc[:100].iterrows(), d8_241030_msms_filtered.iloc[:100].iterrows()):
    
    ############### Precursor Information + Theoretical Calculation ###############
    
    aug_prec = row1['precursor']
    aug_ms2_scan_num = row1['Scan number']
    aug_prec_mz = row1['m/z']
    
    oct_prec = row2['precursor']
    oct_ms2_scan_num = row2['Scan number']
    oct_prec_mz = row2['m/z']
    
    precursor_comp = get_seq_comp(aug_prec[:-1], "M")
    precursor_comp['C'] += 8 # adding for diethyl
    precursor_comp['H[2]'] = 16 # adding for diethyl
    theor_isotope_patterns = isotopic_variants(precursor_comp, npeaks=5, charge=int(aug_prec[-1]))

    theor_mz_values = np.array([peak.mz for peak in theor_isotope_patterns])
    theor_intens = np.array([peak.intensity for peak in theor_isotope_patterns])
    
    monoisotopic_mz = theor_isotope_patterns[0].mz
    ppm_difference = abs(aug_prec_mz - monoisotopic_mz) / monoisotopic_mz * 1e6
    
    if ppm_difference <= 5:
        print(aug_prec_mz)
        print(theor_isotope_patterns)
        
        ############### August ###############
        
        aug_ms1_idx = np.searchsorted(d8_240822_ms1_spec_idxs, aug_ms2_scan_num)
        aug_ms1_scan_num_at_idx = d8_240822_ms1_spec_idxs[aug_ms1_idx]
        aug_ms1_scan_num_at_idx_minus_one = d8_240822_ms1_spec_idxs[aug_ms1_idx - 1]
     
        aug_ms1_scan_idx = aug_ms1_scan_num_at_idx if abs(aug_ms1_scan_num_at_idx - aug_ms2_scan_num) < abs(aug_ms1_scan_num_at_idx_minus_one - aug_ms2_scan_num) else aug_ms1_scan_num_at_idx_minus_one
        aug_ms1_scan = d8_240822_ms1[aug_ms1_scan_idx]
        
        aug_observed_mz = aug_ms1_scan.mz
        aug_observed_intens = aug_ms1_scan.intens
        
        aug_mz_index = np.abs(aug_observed_mz - aug_prec_mz).argmin()
        
        aug_last_10_mz = aug_observed_mz[max(0, aug_mz_index - 10):aug_mz_index]
        aug_next_10_mz = aug_observed_mz[aug_mz_index + 1:aug_mz_index + 11]
        aug_last_10_intens = aug_observed_intens[max(0, aug_mz_index - 10):aug_mz_index]
        aug_next_10_intens = aug_observed_intens[aug_mz_index + 1:aug_mz_index + 11]
        
        aug_combined_mz = np.concatenate((aug_last_10_mz, [aug_observed_mz[aug_mz_index]], aug_next_10_mz))
        aug_combined_intens = np.concatenate((aug_last_10_intens, [aug_observed_intens[aug_mz_index]], aug_next_10_intens))
        
        ############### October ###############
        
        oct_ms1_idx = np.searchsorted(d8_241030_ms1_spec_idxs, oct_ms2_scan_num)
        oct_ms1_scan_num_at_idx = d8_241030_ms1_spec_idxs[oct_ms1_idx]
        oct_ms1_scan_num_at_idx_minus_one = d8_241030_ms1_spec_idxs[oct_ms1_idx - 1]
    
        oct_ms1_scan_idx = oct_ms1_scan_num_at_idx if abs(oct_ms1_scan_num_at_idx - oct_ms2_scan_num) < abs(oct_ms1_scan_num_at_idx_minus_one - oct_ms2_scan_num) else oct_ms1_scan_num_at_idx_minus_one
        oct_ms1_scan = d8_241030_ms1[oct_ms1_scan_idx]
        
        oct_observed_mz = oct_ms1_scan.mz
        oct_observed_intens = oct_ms1_scan.intens
        
        oct_mz_index = np.abs(oct_observed_mz - oct_prec_mz).argmin()
        
        oct_last_10_mz = oct_observed_mz[max(0, oct_mz_index - 10):oct_mz_index]
        oct_next_10_mz = oct_observed_mz[oct_mz_index + 1:oct_mz_index + 11]
        oct_last_10_intens = oct_observed_intens[max(0, oct_mz_index - 10):oct_mz_index]
        oct_next_10_intens = oct_observed_intens[oct_mz_index + 1:oct_mz_index + 11]
        
        oct_combined_mz = np.concatenate((oct_last_10_mz, [oct_observed_mz[oct_mz_index]], oct_next_10_mz))
        oct_combined_intens = np.concatenate((oct_last_10_intens, [oct_observed_intens[oct_mz_index]], oct_next_10_intens))
        
        ############### Graph Theoretical, August, October ###############
        
        plt.figure(figsize=(12, 6))
        
        max_combined_intensity = max(np.max(aug_combined_intens), np.max(oct_combined_intens))
        scaled_theor_intens = theor_intens * (max_combined_intensity / np.max(theor_intens))
        
        bar_width = 0.075
        
        plt.bar(aug_combined_mz, aug_combined_intens, width=bar_width, color='blue', alpha=0.5, label='August')
        plt.bar(oct_combined_mz, oct_combined_intens, width=bar_width, color='red', alpha=0.5, label='October')
        plt.bar(theor_mz_values, scaled_theor_intens, width=bar_width, align='center', color='green', alpha=0.3, label='Theoretical (scaled)')
        
        plt.xlabel('m/z')
        plt.ylabel('Intensity')
        plt.title('Spectra Comparison: August, October, and Theoretical\nSequence: {aug_prec}')
        plt.legend()
        
        plt.show()
