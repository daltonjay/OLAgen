# Copyright 2024, Dalton J. Nelson

# Primer and ligation probe assessment tools

import primer3
from Bio.Seq import Seq

# E484_cp = 'aaggttttaattgttactttccttta'
# E484K_vp = 'ggtagcacaccttgtaatggtgtta'
# E484WT_vp = 'ggtagcacaccttgtaatggtgttg'
# E484_VP_Primer_Region = 'GTTAAGGGAGTGAAGACGATCAGA' # Equal to the forward primer
# E484_CP_Primer_Region = 'TTATGAGAAATCAAAGTCTTTGGGTT'
# E484_Hydr_Region = 'AGTCATCTTTCGAGGTGACTTTTAGATTGCT'

def probe_spurious_ligation(variable_probe, common_probe):
    
    vp_spurious = primer3.calc_heterodimer(variable_probe, variable_probe + common_probe)
    cp_spurious = primer3.calc_heterodimer(common_probe, variable_probe + common_probe)
    
    success_indicator = 1
    
    if vp_spurious.dg < -9000:
        success_indicator = 0 
    elif cp_spurious.dg < -9000:
        success_indicator = 0
    else:
        pass
    
    return success_indicator
        
# probe_spurious_ligation(E484K_vp, E484_cp)

def primer_spurious_ligation(vp_primer_region, hydr_probe_region, cp_primer_region, variable_probe, common_probe):
    
    vp_primer_spurious = primer3.calc_heterodimer(vp_primer_region, variable_probe + common_probe)
    hydr_probe_spurious = primer3.calc_heterodimer(hydr_probe_region, variable_probe + common_probe)
    cp_primer_spurious = primer3.calc_heterodimer(cp_primer_region, variable_probe + common_probe)
    
    success_indicator = 1
    
    if vp_primer_spurious.dg < -9000:
        success_indicator = 0
    elif hydr_probe_spurious.dg < -9000:
        success_indicator = 0
    elif cp_primer_spurious.dg < -9000:
        success_indicator = 0
    else:
        pass
    
    return success_indicator

# primer_spurious_ligation(E484_VP_Primer_Region, E484_Hydr_Region, E484_CP_Primer_Region, E484K_vp, E484_cp)

def spurious_priming(vp_primer_region, hydr_probe_region, cp_primer_region, variable_probe, common_probe):
    hydrolysis_probe = str(Seq(hydr_probe_region).reverse_complement())
    rev_primer = str(Seq(cp_primer_region).reverse_complement())
    
    nucleic_acids = [variable_probe, common_probe, hydr_probe_region, cp_primer_region, vp_primer_region, hydrolysis_probe, rev_primer]
    
    success_indicator = 1
    
    for na in nucleic_acids:
        primed = primer3.calc_heterodimer(variable_probe, na)
        
        if primed.dg < -9000:
            success_indicator = 0
            break
    
    return success_indicator        
    
    
# val = spurious_priming(E484_VP_Primer_Region, E484_Hydr_Region, E484_CP_Primer_Region, E484K_vp, E484_cp)
# print(val)