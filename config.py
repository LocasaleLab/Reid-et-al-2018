"""
Biomass coefficient dict from T. Shlomi, 2011 and Recon2 model.
For further analysis, coefficients are divided to amino acid and nucleotide part.
"""

aa_biomass_coefficient_dict_shlomi = {
    "glycine": 0.62815197,
    "serine": 0.167507192,
    "aspartate": 3.128196812,
    "asparagine": 0.054439837,
    "arginine": 0.054439837
}
nucleotide_biomass_coefficient_dict_shlomi = {
    "AMP": 0.023164655,
    "CMP": 0.038607758,
    "GMP": 0.043755459,
    "UMP": 0.023164655,
    "dAMP": 0.0094984,
    "dCMP": 0.006332267,
    "dGMP": 0.006332267,
    "dTMP": 0.0094984
}

aa_biomass_coefficient_dict_recon2 = {
    "glycine": 0.62815197,
    "serine": 0.167507192,
    "aspartate": 3.128196812,
    "asparagine": 0.054439837,
    "arginine": 0.054439837
}

nucleotide_biomass_coefficient_dict_recon2 = {
    "AMP": 0.023164655,
    "CMP": 0.038607758,
    "GMP": 0.043755459,
    "UMP": 0.023164655,
    "dAMP": 0.0094984,
    "dCMP": 0.006332267,
    "dGMP": 0.006332267,
    "dTMP": 0.0094984
}

biomass_constant_dict_shlomi = {
    "aminoacids": aa_biomass_coefficient_dict_shlomi,
    "nucleotide": nucleotide_biomass_coefficient_dict_shlomi
}
