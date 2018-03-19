#!usr/bin/python3


import itertools as it

import numpy as np
import pandas as pd

import config as configure
import model


def main():
    """
    This script is used to analyze a simple flux model. The original Excel data file with 13C tracing data
    is used as input. Each experiment should be in one sheet. The results of relative fluxes in each 
    experiment are exported to one .csv file. In each experiment, the average and standard deviation of 
    each condition are calculated separately.
    """
    xlsx_file_name = "13C-Glucose_tracing_Mike.xlsx"
    xlsx_file_path = "./{}".format(xlsx_file_name)
    # experiment_sheet_list = [
    #     "HCT116_WQ2101_PHGDH", "HCT116_Raze755_PHGDH", "BT20_WQ2101_PHGDH"]
    experiment_sheet_list = [
        "HCT116_WQ2101_PHGDH", "BT20_WQ2101_PHGDH"]
    flux_analysis(xlsx_file_path, experiment_sheet_list)


def flux_analysis(xlsx_data_file_path, experiment_sheet_list):
    """
    For each experiment (one sheet in 13C-Glucose_tracing_Mike.xlsx), the MS data are extracted and 
    separated by condition prefix. Replicates from the same condition are combined together to calculate
    average and standard deviation. The results for different experiments are exported to separated .csv 
    files.

    :param xlsx_data_file_path: Excel data file path
    :param experiment_sheet_list: List of required sheet name.
    :return: None
    """
    flux_name_list = [
        "3PG-Ser", "Ser-input", "Ser-Gly", "Gly-input", "G6P-R5P", "PYR-TCA-1", "PYR-TCA-2", "Asp-input", "IMP-Syn",
        "IMP-input", "UMP-Syn", "UMP-input"]
    for xlsx_sheet_name in experiment_sheet_list:
        group_data_dict_list = load_and_group_excel_data(xlsx_data_file_path, xlsx_sheet_name)
        # final_result_list_with_stderr = []
        final_result_list_of_replicates = []
        for condition_prefix, replicates_data_dict_list in group_data_dict_list:
            # mean_std_collection_for_current_condition = calculate_replicates_in_one_condition(
            #     replicates_data_dict_list, configure.biomass_constant_dict_shlomi)
            # final_result_list_with_stderr.append([condition_prefix] + mean_std_collection_for_current_condition)
            flux_replicates_for_current_condition = calculate_replicates_in_one_condition(
                replicates_data_dict_list, configure.biomass_constant_dict_shlomi)
            final_result_list_of_replicates.append([condition_prefix] + flux_replicates_for_current_condition)

        output_file_name = "{}/{}.csv".format(
            xlsx_data_file_path[:xlsx_data_file_path.rindex('/')], xlsx_sheet_name)

        with open(output_file_name, 'w') as output_file_object:
            # output_csv_data(final_result_list_with_stderr, flux_name_list, output_file_object)
            output_csv_data(final_result_list_of_replicates, flux_name_list, output_file_object)


def load_and_group_excel_data(xlsx_file_path, xlsx_sheet_name):
    """
    Load all MS data from one experiment. Classify data based on experiment condition. Merge
    isotopomer of each metabolite to generate dict mapping metabolite name to isotopomer data.

    :param xlsx_file_path: Excel data file.
    :param xlsx_sheet_name: The sheet name that contain specific experiment.
    :return: Data dict grouped by condition. All replicates under same condition are grouped together.
    """
    data_frame = pd.read_excel(xlsx_file_path, sheet_name=xlsx_sheet_name, index_col="Name")
    condition_name_list = data_frame.columns
    group_condition_name_list = []
    for condition_prefix, condition_list in it.groupby(condition_name_list, lambda x: x[:x.rindex('_')]):
        group_condition_name_list.append((condition_prefix, list(condition_list)))
    group_metabolite_row_list = []
    current_metabolite_name = ""
    current_metabolite_group_list = []
    for row_index, metabolite_name in enumerate(data_frame.index):
        metabolite_name = metabolite_name.strip()
        if (
                metabolite_name == current_metabolite_name or "[13C]" in metabolite_name or
                ('-' in metabolite_name and metabolite_name[:metabolite_name.rindex('-')]) == current_metabolite_name):
            current_metabolite_group_list.append(row_index)
        else:
            group_metabolite_row_list.append((current_metabolite_name, current_metabolite_group_list))
            if '-' in metabolite_name:
                last_minus_loc = metabolite_name.rindex('-')
                if metabolite_name[last_minus_loc + 1:].isdigit():
                    current_metabolite_name = metabolite_name[:last_minus_loc]
                else:
                    current_metabolite_name = metabolite_name
            else:
                current_metabolite_name = metabolite_name
            current_metabolite_group_list = [row_index]
    group_metabolite_row_list.append((current_metabolite_name, current_metabolite_group_list))
    group_data_dict_list = []
    for condition_prefix, condition_list in group_condition_name_list:
        current_repeat_data_list = []
        for repeat_condition in condition_list:
            current_data_dict = {}
            for metabolite_name, metabolite_row_list in group_metabolite_row_list:
                data_array = np.array(data_frame[repeat_condition].iloc[metabolite_row_list])
                array_sum = data_array.sum()
                if array_sum > 100:
                    current_data_dict[metabolite_name] = data_array / array_sum
            current_repeat_data_list.append(current_data_dict)
        group_data_dict_list.append((condition_prefix, current_repeat_data_list))
    return group_data_dict_list


def calculate_replicates_in_one_condition(current_replicates_data_dict_list, biomass_constant_dict):
    """
    Use data list of replicates in one condition as input. Calculate fluxes for each replicate separately
    and then combine them to get the average and standard deviation

    :param current_replicates_data_dict_list: data dict list of replicates in same condition
    :param biomass_constant_dict: biomass constant dict
    :return: average and standard deviation for each flux.
    """
    flux_list_by_repeat = []
    for data_dict in current_replicates_data_dict_list:
        model.initialize_mid(data_dict)
        flux_list_by_repeat.append(model.calculate_model(biomass_constant_dict))
    return flux_list_by_repeat
    # final_mean_std_collection = [[], []]
    # flux_list_by_flux = zip(*flux_list_by_repeat)
    # for flux_collection in flux_list_by_flux:
    #     flux_collection_nparray = np.array(flux_collection)
    #     final_mean_std_collection[0].append(flux_collection_nparray.mean())
    #     final_mean_std_collection[1].append(flux_collection_nparray.std())
    # return final_mean_std_collection


def output_csv_data(final_result_list_of_replicates, header_row, file_object):
    """
    Export the average flux values and standard deviations in one sheet to a .csv file.

    :param final_result_list_of_replicates: Average and standard deviation of flux value result.
    :param header_row: Flux name list.
    :param file_object: Output file object.
    :return:
    """
    f_out = file_object
    f_out.write("{}\n".format(",".join([""] + header_row)))
    for final_result_in_each_condition in final_result_list_of_replicates:
        # mean_val = [str(i) for i in final_result_in_each_condition[1]]
        # stdev_val = [str(i) for i in final_result_in_each_condition[2]]
        # output_string = "{}\n{}\n{}\n".format(
        #     final_result_in_each_condition[0], ",".join(["Average"] + mean_val),
        #     ",".join(["Standard deviation"] + stdev_val))
        # f_out.write(output_string)
        output_string_list = [final_result_in_each_condition[0]] + [
            ",".join(["Replicate {}".format(index + 1)] + [str(flux) for flux in flux_list])
            for index, flux_list in enumerate(final_result_in_each_condition[1:])]
        f_out.write("\n".join(output_string_list + [""]))


if __name__ == '__main__':
    main()
