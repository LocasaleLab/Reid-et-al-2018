#!usr/bin/python3


import itertools as it

import scipy.misc
import scipy.linalg
import numpy as np


class Mid:
    """
    A class that contains MID data of EMU. It can be set to natural abundance or
    uniformly labeled. It also supports convolution between two EMUs.
    """
    def __init__(self, name, carbon_number, initial_mid=list()):
        self.name = name
        self.carbon_number = carbon_number
        if len(initial_mid) == 0:
            self.mid = [0] * (carbon_number + 1)
        else:
            self.mid = initial_mid

    def update_mid_list(self, mid_list):
        if len(mid_list) == len(self.mid):
            self.mid = list(mid_list)
        else:
            raise Exception("Not standard mid!")

    def update(self, other_mid):
        self.mid = list(other_mid.mid)
        self.name = other_mid.name
        self.carbon_number = other_mid.carbon_number

    def set_natural_distribution(self):
        self.mid = natural_carbon_mid(self.carbon_number)

    def set_uniformly_labeled_distribution(self):
        self.mid = uniformly_labeled_carbon_mid(self.carbon_number)

    def __mul__(self, other):
        new_mid = Mid(self.name, self.carbon_number, list(self.mid))
        for i in range(len(self.mid)):
            new_mid.mid[i] = self.mid[i] * other
        return new_mid

    def __add__(self, other_mid):
        new_mid = Mid(self.name, self.carbon_number, list(self.mid))
        for i in range(len(self.mid)):
            new_mid.mid[i] = self.mid[i] + other_mid.mid[i]
        return new_mid

    def __sub__(self, other_mid):
        new_mid = Mid(self.name, self.carbon_number, list(self.mid))
        for i in range(len(self.mid)):
            new_mid.mid[i] = self.mid[i] - other_mid.mid[i]
        return new_mid

    def __str__(self):
        return "MID:{}({})".format(self.name, self.mid)


def natural_carbon_mid(carbon_number):
    """
    Generate MID of EMU based on natural abundance of 13-carbon (i.e. unlabeled).
    :param carbon_number: number of carbon atoms in the EMU
    :return: Mid class that contain the natural distribution of target EMU.
    """
    c13_ratio = 0.01109
    c12_ratio = 1 - c13_ratio
    new_mid = [0] * (carbon_number + 1)
    for i in range(carbon_number + 1):
        new_mid[i] = scipy.misc.comb(carbon_number, i) * (c13_ratio ** i) * (c12_ratio ** (carbon_number - i))
    return new_mid


def uniformly_labeled_carbon_mid(carbon_number):
    """
    Generate Mid class of EMU which is labeled by 13-carbon uniformly.
    :param carbon_number: Carbon number
    :return: Mid class that contain the uniform 13-carbon distribution of target EMU.
    """
    new_mid = [0] * (carbon_number + 1)
    new_mid[-1] = 1.0
    return new_mid


class LeastSquareSolution:
    """
    Class that store the least squares result.
    """
    def __init__(self, number, residual, original_b, r_square):
        self.number = number
        self.residual = residual
        self.original_b = original_b
        self.r_square = r_square


"""
Initialize basic EMUs.
"""
mid_3pg = Mid('3PG', 3)
mid_3pg_110 = Mid('3PG_110', 2)
mid_3pg_001 = Mid('3PG_001', 1)
mid_g6p = Mid('G6P', 6)
mid_g6p_011111 = Mid('G6P_011111', 5)
mid_r5p = Mid('R5P', 5)
mid_r5p0 = Mid('R5P_0', 5)
mid_ser = Mid('Ser', 3)
mid_ser0 = Mid('Ser_0', 3)
mid_ser0_110 = Mid('Ser_0_110', 2)
mid_ser0_001 = Mid('Ser_0_001', 1)
mid_gly = Mid('Gly', 2)
mid_gly0 = Mid('Gly_0', 2)
mid_imp = Mid('IMP', 10)
mid_imp0 = Mid('IMP_0', 10)
mid_co2 = Mid('CO2', 1)
mid_pyr = Mid('PYR', 3)
mid_pyr_001 = Mid('PYR_100', 1)
mid_pyr_010 = Mid('PYR_010', 1)
mid_pyr_011 = Mid('PYR_110', 2)
mid_oaa_asp = Mid('OAA_Asp', 4)
mid_asp0 = Mid('Asp_0', 4)
mid_asp0_0111 = Mid('ASP_0_0111', 3)
mid_ump = Mid('UMP', 9)
mid_ump0 = Mid('UMP', 9)


def initialize_mid(current_data_dict):
    """
    Initialize MID of EMUs based on data from current_data_dict.
    :param current_data_dict: Data dict of all metabolites from experiment.
    :return: None
    """
    def update_list_of_mid(mid_list, updated_mid_list):
        for mid_name_, updated_mid in zip(mid_list, updated_mid_list):
            mid_name_.update(updated_mid)

    metabolite_name_dict = {
        "3PG/2PG": "mid_3pg",
        "G6P or F6P": "mid_g6p",
        "D-ribose-5-phosphate": "mid_r5p",
        "serine": "mid_ser",
        "glycine": "mid_gly",
        "IMP": "mid_imp",
        "UMP": "mid_ump",
        "D-aspartate(1-)": "mid_oaa_asp",
        "pyruvate": "mid_pyr"
    }

    for metabolite_name, mid_name in metabolite_name_dict.items():
        exec("{}.update_mid_list({})".format(mid_name, list(current_data_dict[metabolite_name])))

    mid_ser0.set_natural_distribution()
    update_list_of_mid([mid_ser0_001, mid_ser0_110], calculate_emu_mid_from_original_mid(mid_ser0, ["001", "110"]))

    update_list_of_mid([mid_3pg_110, mid_3pg_001], calculate_emu_mid_from_original_mid(mid_3pg, ["110", "001"]))
    update_list_of_mid([mid_g6p_011111], calculate_emu_mid_from_original_mid(mid_g6p, ["011111"]))
    update_list_of_mid(
        [mid_pyr_001, mid_pyr_001, mid_pyr_011], calculate_emu_mid_from_original_mid(mid_pyr, ["001", "001", "011"]))

    mid_r5p0.set_natural_distribution()
    mid_gly0.set_natural_distribution()
    mid_imp0.set_natural_distribution()
    mid_co2.set_natural_distribution()
    mid_asp0.set_natural_distribution()
    mid_asp0_0111.set_natural_distribution()
    mid_ump0.set_natural_distribution()


def calculate_model(biomass_coefficient_dict):
    """
    Calculate model based on initialized MID and biomass coefficient dict
    :param biomass_coefficient_dict: Data dict that records the biomass coefficients
    :return: Flux value vector
    """

    aa_constant_dict = biomass_coefficient_dict["aminoacids"]
    nu_constant_dict = biomass_coefficient_dict["nucleotide"]
    c_gly_out = aa_constant_dict["glycine"]
    c_ser_out = aa_constant_dict["serine"]
    c_asp_out = aa_constant_dict["aspartate"]
    c_amp = nu_constant_dict["AMP"]
    c_cmp = nu_constant_dict["CMP"]
    c_gmp = nu_constant_dict["GMP"]
    c_ump = nu_constant_dict["UMP"]
    c_damp = nu_constant_dict["dAMP"]
    c_dcmp = nu_constant_dict["dCMP"]
    c_dgmp = nu_constant_dict["dGMP"]
    c_dtmp = nu_constant_dict["dTMP"]
    c_imp_out = c_amp + c_gmp + c_damp + c_dgmp
    c_ump_out = c_ump + c_cmp + c_dtmp + c_dcmp

    biomass_ratio = 1
    gly_out_flux = c_gly_out * biomass_ratio
    ser_out_flux = c_ser_out * biomass_ratio
    asp_out_flux = c_asp_out * biomass_ratio
    imp_out_flux = c_imp_out * biomass_ratio
    ump_out_flux = c_ump_out * biomass_ratio
    ser_in_result = mid_mix_least_square_solver((mid_ser0, mid_3pg), mid_ser)[0]
    ser_in = ser_in_result.number
    mid_ser_110 = mid_3pg_110 * (1 - ser_in) + mid_ser0_110 * ser_in
    mid_ser_001 = mid_3pg_001 * (1 - ser_in) + mid_ser0_001 * ser_in
    gly_in_result = mid_mix_least_square_solver((mid_gly0, mid_ser_110), mid_gly)[0]
    gly_in = gly_in_result.number
    mid_thf = mid_ser_001
    mid_syn_imp = mid_conv((mid_gly, mid_thf, mid_thf, mid_co2, mid_r5p), 'Mid_IMP_Syn')
    imp_in_result = mid_mix_least_square_solver((mid_imp0, mid_syn_imp), mid_imp)[0]
    imp_in = imp_in_result.number
    r5p_in_result = mid_mix_least_square_solver((mid_r5p0, mid_g6p_011111), mid_r5p)[0]
    r5p_in = r5p_in_result.number
    if r5p_in < 0:
        r5p_in = 0
    mid_oaa_pyr_tca_2 = mid_conv((mid_pyr, mid_co2), 'Mid_OAA_TCA_2')
    mid_oaa_pyr_tca_1 = mid_conv((mid_pyr_001, mid_pyr_001, mid_pyr_011), 'Mid_OAA_TCA_1')
    asp_in_result, oaa_pyr_tca_2_ratio_result = mid_mix_least_square_solver(
        (mid_asp0, mid_oaa_pyr_tca_2, mid_oaa_pyr_tca_1), mid_oaa_asp)[0:2]
    asp_in = asp_in_result.number
    oaa_pyr_tca_2_ratio = oaa_pyr_tca_2_ratio_result.number
    pyr_tca_1_ratio = 1 - asp_in - oaa_pyr_tca_2_ratio
    if pyr_tca_1_ratio < 1e-10:
        pyr_tca_1_ratio = 0
    mid_oaa_asp_0111 = mid_conv((mid_pyr_011, mid_co2)) * (1 - asp_in - oaa_pyr_tca_2_ratio) + mid_conv((
        mid_pyr_001, mid_pyr_001, mid_pyr_001)) * oaa_pyr_tca_2_ratio + mid_asp0_0111 * asp_in

    mid_cap = mid_co2
    mid_syn_ump = mid_conv((mid_oaa_asp_0111, mid_cap, mid_r5p), 'Mid_UMP_Syn')
    ump_in_result = mid_mix_least_square_solver((mid_ump0, mid_syn_ump), mid_ump)[0]
    ump_in = ump_in_result.number

    imp_syn = (1 - imp_in) * imp_out_flux
    imp_input = imp_in * imp_out_flux
    ser_gly = (1 - gly_in) * (imp_syn + gly_out_flux)
    gly_input = gly_in * (imp_syn + gly_out_flux)
    pg_ser = (1 - ser_in) * (ser_gly + ser_out_flux)
    ser_input = ser_in * (ser_gly + ser_out_flux)
    ump_syn = (1 - ump_in) * ump_out_flux
    ump_input = ump_in * ump_out_flux
    g6p_r5p = (1 - r5p_in) * (imp_syn + ump_syn)
    pyr_tca_1 = pyr_tca_1_ratio * (ump_syn + asp_out_flux)
    pyr_tca_2 = oaa_pyr_tca_2_ratio * (ump_syn + asp_out_flux)
    asp_input = asp_in * (ump_syn + asp_out_flux)

    final_flux_vector = [
        pg_ser, ser_input, ser_gly, gly_input, g6p_r5p, pyr_tca_1, pyr_tca_2, asp_input, imp_syn,
        imp_input, ump_syn, ump_input]
    return final_flux_vector


def mid_conv(mid_tuple, name=""):
    """
    Calculate convolution of multiple Mid objects. Return the new Mid object.
    :param mid_tuple: A tuple that contains all Mid object of substrate EMUs.
    :param name: Name of new Mid object.
    :return: Mid object. Generated convolution result.
    """

    def mid_double_conv(mid1, mid2, name_=""):
        carbon_number = mid1.carbon_number + mid2.carbon_number
        return Mid(name_, carbon_number, convolution(mid1.mid, mid2.mid))

    def convolution(list1, list2):
        new_list = [0] * (len(list1) + len(list2) - 1)
        for index_2 in range(len(list2)):
            for index_1 in range(len(list1)):
                new_list[index_1 + index_2] += list1[index_1] * list2[index_2]
        return new_list

    if len(mid_tuple) == 1:
        mid_tuple[0].name = name
        return mid_tuple[0]
    elif len(mid_tuple) == 2:
        return mid_double_conv(mid_tuple[0], mid_tuple[1], name)
    else:
        return mid_conv((mid_double_conv(mid_tuple[0], mid_tuple[1]),) + mid_tuple[2:], name)


def mid_mix_least_square_solver(mid_substrate_tuple, mid_mixed):
    """
    Solve the mix ratio of multiple substrate by least squares method.

    :param mid_substrate_tuple: A tuple of Mid object that contains all original source metabolite MID.
    :param mid_mixed: A Mid object that represents mixed MID.
    :return: List of mix ratio for each original source metabolite. Sum of those ratios equals to 1.
    """

    substrate_num = len(mid_substrate_tuple)
    mid_length = len(mid_mixed.mid)
    coefficient_matrix = np.zeros([mid_length, 0])
    for current_substrate_mid in mid_substrate_tuple:
        current_mid_vector = np.reshape(np.array(current_substrate_mid.mid), [mid_length, 1])
        coefficient_matrix = np.hstack([coefficient_matrix, current_mid_vector])
    converted_coefficient_matrix = np.zeros([mid_length, substrate_num - 1])
    for column_index in range(np.size(coefficient_matrix, 1) - 1):
        converted_coefficient_matrix[:, column_index] = coefficient_matrix[:, column_index] - coefficient_matrix[:, -1]
    result_vector = np.array(mid_mixed.mid)
    converted_result_vector = result_vector - np.array(mid_substrate_tuple[-1].mid)
    solution, residual = np.linalg.lstsq(converted_coefficient_matrix, converted_result_vector)[0:2]
    valid_solution_modification = False
    solution = np.hstack([solution, 1 - np.sum(solution)])
    for solution_index in range(len(solution)):
        if solution[solution_index] < 0:
            solution[solution_index] = 0
            valid_solution_modification = True
    if valid_solution_modification:
        solution = solution / np.sum(solution)
        residual = [np.sum((result_vector - np.dot(coefficient_matrix, solution)) ** 2)]
    r_square = 1 - residual[0] / (len(result_vector) * result_vector.var())
    solution_object_list = []
    for current_solution in solution:
        solution_object_list.append(LeastSquareSolution(current_solution, residual, result_vector, r_square))
    return solution_object_list


def calculate_emu_mid_from_original_mid(original_mid, required_emu_string_list):
    """
    Calculate MID of EMUs of input substrate like PYR, 3PG and G6P.
    :param original_mid: MID of original input substrate from experiment data.
    :param required_emu_string_list: Required EMU list represented by string, like "00011"
    :return: List of Mid object that contain MID of specific EMU.
    """
    def flat_list(current_list, depth=-1):
        final_list = []
        for item in current_list:
            if type(item) == list:
                if depth == 0:
                    final_list.append(item)
                else:
                    final_list += flat_list(item, depth - 1)
            else:
                final_list.append(item)
        return final_list

    carbon_number = original_mid.carbon_number
    required_emu_carbon_number_list = []
    required_emu_mid_list = []
    required_emu_list = []
    required_mid_class_list = []
    for required_emu_string in required_emu_string_list:
        if len(required_emu_string) != carbon_number:
            raise RuntimeError("Unmatched EMU string length and original MID when calculating EMU MID")
        required_emu_carbon_number_list.append(required_emu_string.count("1"))
        required_emu_mid_list.append([0] * (required_emu_carbon_number_list[-1] + 1))
        required_emu_list.append([bool(int(i)) for i in list(required_emu_string)])

    original_percentage_per_carbon_number = [
        original_mid.mid[i] / (scipy.misc.comb(carbon_number, i, exact=True)) for i in range(carbon_number + 1)]
    emu_combination_for_each_carbon_number = [0] * (carbon_number + 1)
    original_percentage_for_each_emu = [0] * (carbon_number + 1)
    for carbon_number_index in range(carbon_number + 1):
        final_result = []
        comb_itertools = it.combinations(range(carbon_number), carbon_number_index)
        for set_true_locations in comb_itertools:
            basic_valid_location_list = [False] * carbon_number
            for current_loc in set_true_locations:
                basic_valid_location_list[current_loc] = True
            final_result.append(basic_valid_location_list)
        emu_combination_for_each_carbon_number[carbon_number_index] = final_result
        original_percentage_for_each_emu[carbon_number_index] = (
            [original_percentage_per_carbon_number[carbon_number_index]] * len(final_result))
    emu_combination = flat_list(emu_combination_for_each_carbon_number, 1)
    original_percentage = flat_list(original_percentage_for_each_emu, 1)

    for required_emu_index in range(len(required_emu_list)):
        for current_emu, each_original_percent in zip(emu_combination, original_percentage):
            overlap_emu = [i and j for i, j in zip(required_emu_list[required_emu_index], current_emu)]
            overlap_emu_carbon_number = overlap_emu.count(True)
            required_emu_mid_list[required_emu_index][overlap_emu_carbon_number] += each_original_percent
        required_mid_class_list.append(Mid(
            "{}_{}".format(original_mid.name, required_emu_string_list[required_emu_index]),
            required_emu_carbon_number_list[required_emu_index], required_emu_mid_list[required_emu_index]))

    return required_mid_class_list
