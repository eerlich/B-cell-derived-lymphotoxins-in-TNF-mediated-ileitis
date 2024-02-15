#!/usr/bin/env python

import sys
import os.path
import csv
import datetime
from itertools import combinations

class Barcode:
    """Stores the data associated with a barcode retaining sample info."""

    def __init__(self, barcode, raw_row_content, meta):
        self.barcode = barcode
        self.all_row_list = [raw_row_content]
        self.sample = self.barcode.split("_")[1]
        self.metadata = meta
        self.valid_chains = True
        self.clonal = True
        self.prelim_clonotype = None
        self.final_clonotype = None

    def set_chains(self):
        pass

class TcrBarcode(Barcode):
    """Stores the data associated with TCR data."""

    def __init__(self, barcode, raw_row_content, meta):
        super(TcrBarcode, self).__init__(barcode, raw_row_content, meta)
        self.tcra_list = []
        self.tcra_output = []
        self.tcrb = ""

    def set_chains(self):
        temp_tcrb_list = []
        #valid tcr barcodes can have either 2 or 3 rows associated with it
        if len(self.all_row_list) > 1 and len(self.all_row_list) < 4:
            for row in self.all_row_list:
                if row[5] == "TRA":
                    #tcrav_tcraj_cdr1nt_cdr2nt_cdr3nt
                    tcra = f"{row[6]}_{row[8]}_{row[15]}_{row[19]}_{row[23]}"
                    self.tcra_list.append(tcra)
                    self.tcra_output.append(f"{row[23]}")
                elif row[5] == "TRB":
                    #tcrbv_tcrbj_tcrbj_cdr1nt_cdr2nt_cdr3nt
                    tcrb = f"{row[6]}_{row[7]}_{row[8]}_{row[15]}_{row[19]}_{row[23]}"
                    temp_tcrb_list.append(tcrb)
        #there can be only one TCRB associated with a barcode
        if len(temp_tcrb_list) != 1:
            self.valid_chains = False
        #there can be 1 or 2 TCRAs associated with a barcode
        elif len(self.tcra_list) < 1 or len(self.tcra_list) > 2:
            self.valid_chains = False
        else:
            self.tcrb = temp_tcrb_list.pop()
            self.tcra_output.sort()
            if len(self.tcra_output) < 2:
                self.tcra_output.append("")


class BcrBarcode(Barcode):
    """Stores the data associated with BCR data."""

    def __init__(self, barcode, raw_row_content, meta):
        super(BcrBarcode, self).__init__(barcode, raw_row_content, meta)
        self.igh_genes = ""
        self.igh_cdr3nn = ""
        self.igl_genes = ""
        self.igl_cdr3nn = ""
        self.heavy_constant = ""
        self.vdj_compare = ""
        self.cdr3nns = ""

    def set_chains(self):
        #valid bcr barcode must have exactly one igh and one igl chain
        #so must have only two rows
        if len(self.all_row_list) == 2:
            for row in self.all_row_list:
                if row[5] in ["IGK", "IGL"]:
                    #iglvgene_igljgene
                    self.igl_genes = f"{row[6]}_{row[8]}"
                    self.igl_cdr3nn = row[23]
                elif row[5] == "IGH":
                    #ighvgene_ighvgene
                    self.igh_genes = f"{row[6]}_{row[8]}"
                    self.igh_cdr3nn = row[23]
                    self.heavy_constant = row[9]
        else:
            self.valid_chains = False
        #now check that there's a igh_constant and a igl_constant
        if not self.igh_genes and self.igl_genes:
            self.valid_chains = False

        #make the string with all the vdj_compare info
        self.vdj_compare = (f"{self.igh_genes}_{len(self.igh_cdr3nn)}_"
                           f"{self.igl_genes}_{len(self.igl_cdr3nn)}")
        #make the string with the cdr3nn info
        self.cdr3nns = f"{self.igh_cdr3nn}{self.igl_cdr3nn}"


def ask_path():
    """Ask the user for a configuration file with the contig.csv file(s)."""

    path_ask_text = ("What is the path to a config file with the "
                    "filtered_contig_annotations.csv file(s)?")
    print(path_ask_text)
    while True:
        path = input("> ")
        path = path.strip('"')
        file_type = os.path.splitext(path)[1]
        if not os.path.isfile(path):
            print("This isn't a file.")
            sys.exit()
        if not file_type == ".csv":
            print("The config file should be a .csv file.")
            sys.exit()
        break
    return(path)

def get_data_from_all_contig_files(path):
    """Read the config file and get the file contents from the contig files.

    filtered_contig_annotations.csv files are generated from Cell Ranger from
    10X. Tested with version cellranger-7.1+. Configuration files for this 
    script will have the path to the filtered_contig_annotations.csv file(s) in 
    the first column and the number to append to the barcodes in the second 
    column as a sample id.

    When integrating, Seurat adds "_\d" to the end of barcodes in
    order of how they were in the object.list in the SelectIntegrationFeatures/
    FindIntegrationAnchors calls. To aid adding the output to the Seurat object 
    metadata, this script allows for the number to be appended to be specified.
    The number must be unique for each type of data.

    Returns a dictonary where the key is user assigned number and the value is
    the filtered_contig_annotations file contents (minus the header).
    """
    contig_file_content = []
    with open(path, "r", encoding = "utf-8-sig") as f:
      temp_content = csv.reader(f)
      for x in temp_content:
        contig_file_content.append(x)
    header = contig_file_content.pop(0)
    if header != ['path', '#', 'type'] and header != ['path', '#', 'type', 'metadata']:
        print("The config file isn't formatted correctly. The header row is "
              "wrong.\nIt should be path (first column), # (second column), "
              "type (third column, bcr or tcr), (optional) metadata (fourth)")
        sys.exit()
    
    contig_tcr_dict = {}
    contig_bcr_dict = {}
    
    for line in contig_file_content:
        meta = ""
        try:
            meta = line[3]
        except IndexError:
            pass
        content, immuno = read_contig_file(line[0])
        if immuno == "bcr" and line[2] == "bcr":
            contig_bcr_dict[line[1]] = [content, meta]
        elif immuno == "tcr" and line[2] == "tcr":
            contig_tcr_dict[line[1]] = [content, meta]
        else:
            print(f"The data in file {line[0]} appears to be {immuno}, but "
                  f"it is {line[2]} in the configuration file.")
            sys.exit()

    contig_dict = {"tcr": contig_tcr_dict, "bcr": contig_bcr_dict}

    return contig_dict

def read_contig_file(file_info):
    """Read contig file and test if truly bcr or tcr."""

    file_content = []
    try:
        with open(file_info, "r", encoding = "utf-8-sig") as f:
            temp_content = csv.reader(f)
            for line in temp_content:
                file_content.append(line)
    except IOError:
        print(f"{file_info} can't be opened.")
        sys.exit()
    
    #test the header has all the needed information in the right order
    #also remove the header because it's not necessary
    try:
        row = file_content.pop(0)
    except IndexError:
        print(f"{file_info} appears to be an empty file.")
        sys.exit()
    try:
        test = f"{row[6]}_{row[7]}_{row[8]}_{row[9]}_{row[15]}_{row[19]}_{row[23]}"
    except IndexError:
        print(f"{file_info} does not have the expected header for a "
               "filtered_contig_annotation.csv file from CellRanger.")
        sys.exit()
    if not test == "v_gene_d_gene_j_gene_c_gene_cdr1_nt_cdr2_nt_cdr3_nt":
        print(f"{file_info} does not have the expected header for a "
               "filtered_contig_annotation.csv file from CellRanger.")
        sys.exit()

    
    #test if the file has tcr or bcr data
    if file_content[0][5] in ["IGK", "IGL", "IGH"]:
        immuno_type = "bcr"
    elif file_content[0][5] in ["TRA", "TRB", "TRG", "TRD"]:
        immuno_type = "tcr"
    else:
        print(f"It's unclear if '{file_info}' contains tcr or bcr data. The "
              f"first data has '{file_content[0][5]}' in the first chain "
                "column.")
        sys.exit()
    return file_content, immuno_type

def filter_barcodes(content_dict, is_tcr):
    """Process the input data to retain only barcodes with useful reads.

    Filters barcodes to retain only barcodes with a biologically appropriate
    number of chains. Also appends the user provided id number to the barcodes.

    Output is a filtered dictionary {barcode_sample: vdj barcode objects}
    """
    all_bar_dict = {}
    if is_tcr:
        content_dict = content_dict["tcr"]
    else:
        content_dict = content_dict["bcr"]
    for id, data in content_dict.items():
        for row in data[0]:
            bar_id = f"{row[0]}_{id}"
            try:
                all_bar_dict[bar_id].all_row_list.append(row)
            except KeyError:
                if is_tcr:
                    all_bar_dict[bar_id] = TcrBarcode(bar_id, row, data[1])
                else:
                    all_bar_dict[bar_id] = BcrBarcode(bar_id, row, data[1])
    filtered_bar_dict = {}
    for vdj_obj in all_bar_dict.values():
        vdj_obj.set_chains()
        if vdj_obj.valid_chains:
            filtered_bar_dict[vdj_obj.barcode] = vdj_obj
    return filtered_bar_dict

def check_barcode_overlap(bcr_dict, tcr_dict):
    """Delete barcodes with both tcr and bcr data."""
    overlap = bcr_dict.keys() & tcr_dict.keys()
    for invalid_bar in overlap:
        del bcr_dict[invalid_bar]
        del tcr_dict[invalid_bar]
    return bcr_dict, tcr_dict

def find_tcr_clones(tcr_list):
    """Group all of the TCR chains that could be clonal.

    First, make groups of TCRs where all of the TCRB data is completely
    identical. Then iterate through the tcrb_dict and determine if they are
    clonal by comparing the TCRA chain(s). Since T cells can express two TCRA
    chains and since scRNAseq can have dropouts where only one, instead of 
    both, TCRA chains is associated with a barcode, it's possible that a valid
    TCR clonal group with have some barcodes with TCRA-1, some with TCRA-2, and
    some with both TCRA-1 and TCRA-2.

    Assigns arbitrary preliminary clonotype #s to the objects.
    """
    #key is tcrb, value is a list of tcr object(s) associated with it
    inital_tcr_dict = {}
    #group tcr objects by identical tcrb chains
    for vdj_obj in tcr_list:
        try:
            inital_tcr_dict[vdj_obj.tcrb].append(vdj_obj)
        except KeyError:
            inital_tcr_dict[vdj_obj.tcrb] = [vdj_obj]
    for count, tcr_obj_list in enumerate(inital_tcr_dict.values()):
        if len(tcr_obj_list) == 1:
            tcr_obj_list[0].prelim_clonotype = f"{count}"
        else:
            #key is are tcra(s) as a frozenset
            #value are the associated tcr objects
            temp_tcra_dict = {}
            for tcr in tcr_obj_list:
                tcra_set = frozenset(tcr.tcra_list)
                try:
                    temp_tcra_dict[tcra_set].append(tcr)
                except KeyError:
                    temp_tcra_dict[tcra_set] = [tcr]
            #if there's only one TCRA, then everything is completely identical
            if len(temp_tcra_dict) == 1:
                for tcr_obj in tcr_obj_list:
                    tcr_obj.prelim_clonotype = f"{count}"
            else:
                #make a dict of the tcra_sets with two tcras as a key,
                #value is a list of the one_tcra sets that are subseet of key
                two_tcras = {}
                #make a dictonary of the tcra_sets with one tcra as key,
                #value is a list of two_tcras sets that are a superset of key
                one_tcras = {}
                for tcra_set in temp_tcra_dict.keys():
                    if len(tcra_set) > 1:
                        two_tcras[tcra_set] = []
                    else:
                        one_tcras[tcra_set] = []
                if not two_tcras:
                    for n, tcr_obj in enumerate(tcr_obj_list):
                        tcr_obj.prelim_clonotype = f"{count}_{n}"
                else:
                    for one_tcra in one_tcras.keys():
                        for two_tcra in two_tcras.keys():
                            if one_tcra <= two_tcra:
                                one_tcras[one_tcra].append(two_tcra)
                                two_tcras[two_tcra].append(one_tcra)
                    final_grouping = set([])
                    for single, double in one_tcras.items():
                        #the single tcra is associated with more than one set
                        #of two tcra sets, no ground truth resolution
                        if len(double) > 1:
                            final_grouping.add(frozenset([single]))
                            for duo in double:
                                final_grouping.add(frozenset([duo]))
                                del two_tcras[duo]
                        #the single tcra isn't associated with any other tcras
                        elif len(double) < 1:
                            final_grouping.add(frozenset([single]))
                    for two, single_list in two_tcras.items():
                        final_grouping.add(frozenset([two, *single_list]))
                    for n, group in enumerate(final_grouping):
                        for subset in group:
                            for tcr_obj in temp_tcra_dict[subset]:
                                tcr_obj.prelim_clonotype = f"{count}_{n}"

def assign_mixed_tcra(a_list, long):
    """Assign clonotypes when there's a shared TCRB and mixed TCRAs.

    Specifically, this is for the rare cases where multiple barcodes share
    identical TCRBs, at least one of the barcodes are associated with more than
    one TCRA, but not all of the barcodes should be the same clonotype. For
    example, if there's a group of A, B, C, and AB, A, B, and AB should be one
    clonotype and C should be another.
    """
    subset_of_long = []
    not_in_long = []
    for x in a_list:
        if x <= long:
            subset_of_long.append(x)
        else:
            not_in_long.append(x)
    return subset_of_long, not_in_long

def find_bcr_clones(bcr_list):
    """Group all of the BCR clones that could be clonal.

    Start by grouping all barcodes that have the same igh_v, igh_j,
    igh_cdr3 length, igl_v, igl_j, and igl_cdr3_nn length. To be considered 
    clonal, all 6 must be identical. Then, the amino acid sequences of the 
    candidate cdr3s are compared. If the ratio of amino acid changes/cdr3 
    length is less than 0.1 when compared to at least one other member of the 
    group, then they are defined as part of the same clonal group.

    Assigns arbitrary preliminary clonotype #s to the objects.
    """
    #key is vdj compare, value is a list of bcr object(s) associated with it
    inital_bcr_dict = {}
    #group bcr objects by identical vdj chains and cdr3 lengths
    for vdj_obj in bcr_list:
        try:
            inital_bcr_dict[vdj_obj.vdj_compare].append(vdj_obj)
        except KeyError:
            inital_bcr_dict[vdj_obj.vdj_compare] = [vdj_obj]
    
    for count, bcr_obj_list in enumerate(inital_bcr_dict.values()):
        #if there's only one bcr object associated with the info, then it 
        #cannot be clonal with anything else
        if len(bcr_obj_list) == 1:
            bcr_obj_list[0].prelim_clonotype = f"{count}"
        else:
            #key is cdr3nns, value is a list of bcr object(s)
            cdr3_dict = {}
            for bcr_obj in bcr_obj_list:
                try:
                    cdr3_dict[bcr_obj.cdr3nns].append(bcr_obj)
                except KeyError:
                    cdr3_dict[bcr_obj.cdr3nns] = [bcr_obj]
            #the cdr3 AA sequences of all the associated bcrs are identical
            if len(cdr3_dict) == 1:
                for bcr_obj in bcr_obj_list:
                    bcr_obj.prelim_clonotype = f"{count}"
            else:
                clonal_list = []
                for combo in combinations(cdr3_dict.keys(), 2):
                    differences = 0
                    clonal = True
                    #compare amnio acid by amnio acid
                    for cdr3 in zip(*combo):
                        if (differences/len(combo[0])) > 0.1:
                            clonal = False
                            break
                        elif cdr3[0] != cdr3[1]:
                            differences += 1
                    #get the lists of associated bcr objects
                    bcr_list_a = cdr3_dict[combo[0]]
                    bcr_list_b = cdr3_dict[combo[1]]
                    if clonal:
                        bcr_list_a.extend(bcr_list_b)
                        clonal_list.append(set(bcr_list_a))
                    else:
                        clonal_list.append(set(bcr_list_a))
                        clonal_list.append(set(bcr_list_b))
                
                no_more_matches = []
                while len(clonal_list) > 1:
                    end = clonal_list.pop()
                    temp_list = []
                    no_match = 0
                    for test in clonal_list:
                        if test.isdisjoint(end):
                            temp_list.append(test)
                            no_match += 1
                        else:
                            # | is for bcrs in test, end, or both
                            temp_list.append(test | end)
                    if no_match == len(clonal_list):
                        no_more_matches.append(end)
                    clonal_list = temp_list
                clonal_list.extend(no_more_matches)
                for num, clonal_group in enumerate(clonal_list):
                    for bcr_obj in clonal_group:
                        bcr_obj.prelim_clonotype = f"{count}_{num}"

def sort_clonal_groups(vdj_list):
    """Renumber the clonotypes by frequency and return the sorted list."""
    prelim_clono_dict = {}
    for vdj_obj in vdj_list:
        try:
            prelim_clono_dict[vdj_obj.prelim_clonotype].append(vdj_obj)
        except KeyError:
            prelim_clono_dict[vdj_obj.prelim_clonotype] = [vdj_obj]
    sorted_prelim = {k: v for k, v in sorted(prelim_clono_dict.items(),
                     key=lambda item: len(item[1]), reverse = True)}
    sorted_vdj_list = []
    for count, vdj_obj_list in enumerate(sorted_prelim.values(), start = 1):
        for vdj_obj in sorted(vdj_obj_list, key=lambda bar: bar.sample):
            vdj_obj.final_clonotype = count
            sorted_vdj_list.append(vdj_obj)
            #if it's not clonal, modify object accordingly
            if len(vdj_obj_list) == 1:
                vdj_obj.clonal = False
    return sorted_vdj_list

def write_output_file(a_path, is_tcr, output_list):
    """Write the output .csv file."""
    header = ["barcode", "clonotype", "is_clonal"]
    if is_tcr:
        output_file_name = "tcr_output.csv"
        header.extend(["tcrb_cdr3_nn", "tcra_cdr3_nn_1", "tcra_cdr3_nn_2"])
    else:
        output_file_name = "bcr_output.csv"
        header.extend(["igh_cdr3_nn", "igl_cdr3_nn", "igh_constant"])
    if output_list[0].metadata:
        meta = True
        header.append("metadata")
    else:
        meta = False
    output_path = os.path.join(os.path.dirname(a_path), output_file_name)
    data = []
    for vdj_obj in output_list:
        row_output = []
        row_output.extend([vdj_obj.barcode, 
                          vdj_obj.final_clonotype, 
                          vdj_obj.clonal])
        if is_tcr:
            row_output.append(vdj_obj.tcrb.split("_")[-1])
            row_output.extend(vdj_obj.tcra_output)
        else:
            row_output.extend([vdj_obj.igh_cdr3nn, 
                              vdj_obj.igl_cdr3nn, 
                              vdj_obj.heavy_constant])
        if meta:
            row_output.append(vdj_obj.metadata)
        data.append(row_output)
    with open(output_path, "w") as f:
        file_content = csv.writer(f, lineterminator="\n")
        file_content.writerow(header)
        file_content.writerows(data)

def main():
    config_path = ask_path()
    computation_start = datetime.datetime.now()
    raw_data_dict = get_data_from_all_contig_files(config_path)
    
    filtered_bcr_dict = False
    filtered_tcr_dict = False
    if raw_data_dict["bcr"]:
        filtered_bcr_dict = filter_barcodes(raw_data_dict, False)
    if raw_data_dict["tcr"]:
        filtered_tcr_dict = filter_barcodes(raw_data_dict, True)
    if filtered_bcr_dict and filtered_tcr_dict:
        filtered_bcr_dict, filtered_tcr_dict = check_barcode_overlap(filtered_bcr_dict, filtered_tcr_dict)
    
    if filtered_tcr_dict:
        filtered_list = list(filtered_tcr_dict.values())
        find_tcr_clones(filtered_list)
        sorted_list = sort_clonal_groups(filtered_list)
        write_output_file(config_path, True, sorted_list)
    if filtered_bcr_dict:
        filtered_list = list(filtered_bcr_dict.values())
        find_bcr_clones(filtered_list)
        sorted_list = sort_clonal_groups(filtered_list)
        write_output_file(config_path, False, sorted_list)
    print("Done. Script execution time: "
         f"{datetime.datetime.now() - computation_start}")

if __name__ == '__main__':
    main()
