#!/usr/bin/env python3

import os
from typing import Dict
import Bio.Phylo
import Bio.Phylo.Consensus
from Bio import SeqIO
import requests
import ete3
import re
import glob
from Bio.Align.Applications import MuscleCommandline
import statistics
import pandas as pd

import subprocess


class PhyloAnalyser:
    def __init__(self, input_file_path, analysis_dir_path="./"):
        """
        :param input_file_path: path to file with taxids, 1 taxid per line
        :param analysis_dir_path: path to directory where taxid_to_analyse.txt needs to be stored, where analysis will be performed
        """
        # Reading input file
        self.in_taxid_file_path = input_file_path
        self.taxid_dict = self.read_taxid_input()

        # Setting up necessary directories
        self.analysis_loc = analysis_dir_path
        self.path_archive = f"{self.analysis_loc}/archive"
        self.path_ft_all_inferred = f"{self.analysis_loc}/tree/ft/all/inferred"
        self.path_ft_all_bootstrapped = f"{self.analysis_loc}/tree/ft/all/bootstrapped"

        self.path_ft_all_boot_done_supertree = f"{self.analysis_loc}/tree/ft/all/boot_done/supertree"
        self.path_ft_all_boot_nondone_supertree = f"{self.analysis_loc}/tree/ft/all/boot_nondone/supertree"

        self.path_ft_filtered_inferred = f"{self.analysis_loc}/tree/ft/filtered/inferred"
        self.path_ft_filtered_bootstrapped = f"{self.analysis_loc}/tree/ft/filtered/bootstrapped"

        self.path_ft_filtered_boot_done_consensus = f"{self.analysis_loc}/tree/ft/filtered/boot_done/consensus"
        self.path_ft_filtered_boot_nondone_consensus = f"{self.analysis_loc}/tree/ft/filtered/boot_nondone/consensus"
        self.path_ft_filtered_boot_done_supertree = f"{self.analysis_loc}/tree/ft/filtered/boot_done/supertree"
        self.path_ft_filtered_boot_nondone_supertree = f"{self.analysis_loc}/tree/ft/filtered/boot_nondone/supertree"

        self.path_tree_calculated = f"{self.analysis_loc}/tree/calculated"

        self.path_tree_concat = f"{self.analysis_loc}/tree/concatenated"

        self.path_tree_figures = f"{self.analysis_loc}/tree_figures"

        self.list_of_consensus_dirs = [self.path_ft_filtered_boot_done_consensus,
                                       self.path_ft_filtered_boot_nondone_consensus]

        self.list_of_supertree_dirs = [self.path_ft_all_boot_done_supertree,
                                       self.path_ft_all_boot_nondone_supertree,
                                       self.path_ft_filtered_boot_done_supertree,
                                       self.path_ft_filtered_boot_nondone_supertree]

        self.list_of_tree_dirs = [self.path_ft_all_inferred,
                                  self.path_ft_all_bootstrapped,
                                  self.path_ft_all_boot_done_supertree,
                                  self.path_ft_all_boot_nondone_supertree,
                                  self.path_ft_filtered_inferred,
                                  self.path_ft_filtered_bootstrapped,
                                  self.path_ft_filtered_boot_done_consensus,
                                  self.path_ft_filtered_boot_nondone_consensus,
                                  self.path_ft_filtered_boot_done_supertree,
                                  self.path_ft_filtered_boot_nondone_supertree,
                                  self.path_tree_concat,
                                  self.path_tree_calculated]

        self.setup_dirs()  # needs to be run after self.analysis_loc

    def setup_dirs(self) -> None:
        """Method that sets up necessary directories for whole project to work"""
        if self.analysis_loc[-1] == "/":
            self.analysis_loc = self.analysis_loc[:-1]
        if not os.path.exists(self.path_archive):
            os.makedirs(self.path_archive)
        # Creating "tree" directories
        if not os.path.exists(f"{self.analysis_loc}/tree/ft/all"):
            os.makedirs(f"{self.analysis_loc}/tree/ft/all")
        if not os.path.exists(f"{self.analysis_loc}/tree/ft/filtered"):
            os.makedirs(f"{self.analysis_loc}/tree/ft/filtered")
        if not os.path.exists(f"{self.analysis_loc}/tree/concatenated"):
            os.makedirs(f"{self.analysis_loc}/tree/concatenated")
        for dir_path in self.list_of_tree_dirs:
            if not os.path.exists(dir_path):
                os.makedirs(dir_path)
        if not os.path.exists(self.path_tree_figures):
            os.makedirs(self.path_tree_figures)
        # Creating "seq" directories
        if not os.path.exists(f"{self.analysis_loc}/seq/organisms/stock_ids"):
            os.makedirs(f"{self.analysis_loc}/seq/organisms/stock_ids")
        if not os.path.exists(f"{self.analysis_loc}/seq/organisms/renamed"):
            os.makedirs(f"{self.analysis_loc}/seq/organisms/renamed")
        if not os.path.exists(f"{self.analysis_loc}/seq/organisms/concatenated"):
            os.makedirs(f"{self.analysis_loc}/seq/organisms/concatenated")
        # Creating "cluster" directories:
        if not os.path.exists(f"{self.analysis_loc}/cluster/all"):  # used for supertree
            os.makedirs(f"{self.analysis_loc}/cluster/all")
        if not os.path.exists(f"{self.analysis_loc}/cluster/filtered"):  # after filtering
            os.makedirs(f"{self.analysis_loc}/cluster/filtered")
        # Creating "alignment" directories
        if not os.path.exists(f"{self.analysis_loc}/alignment/all"):  # used for supertree
            os.makedirs(f"{self.analysis_loc}/alignment/all")
        if not os.path.exists(f"{self.analysis_loc}/alignment/filtered"):  # after filtering
            os.makedirs(f"{self.analysis_loc}/alignment/filtered")

    def read_taxid_input(self) -> Dict:
        """
        Method which is used to translate taxids into species names

        :return: Returns taxid2name dictionary object
        """
        ncbi = ete3.NCBITaxa()
        with open(self.in_taxid_file_path) as taxid_input:
            taxids_list = [int(taxid.strip("\n")) for taxid in taxid_input.readlines() if
                           taxid is not "\n"]
        taxid2name = ncbi.get_taxid_translator(taxids_list)
        return taxid2name

    def shorten_taxid_name(self, taxid_name: str) -> str:
        """Method which shortens taxon name to 6 uppercase letters, taking 2 to words as label

        :param taxid_name: Taxid name to shorten
        :return: 6 letter uppercase string
        """
        name_to_list = taxid_name.split(" ")
        first_two_names = [name[:3].upper() for name in
                           name_to_list[:2]]  # taking just 2 first words of name
        return "".join(first_two_names)

    def download_proteome(self, taxon_id: int, min_seq_len: int = 30, max_seq_len: int = 1200):
        """
        Method uses uniprot url with taxids to retrieve proteomes for further analysis

        :param taxon_id: Taxid id retrieved from input file
        :param min_seq_len: minimum length of accepted sequences
        :param max_seq_len: maximum length of accepted sequences
        :return: Nothing, but saves fitlered proteomes in specified directories
        """
        taxon_name = self.taxid_dict[taxon_id]
        short_taxon_name = self.shorten_taxid_name(taxon_name)
        filename = f"proteome_{short_taxon_name}.fasta"
        print(f"Starting to download {taxon_name}.")
        url = f"https://www.uniprot.org/uniprot/?query=organism:{str(taxon_id)}&format=fasta"
        requests_object = requests.get(url)
        # Writing unmodified file
        if os.path.exists(f"{self.analysis_loc}/seq/organisms/stock_ids/{filename}"):
            print(f"FASTA file for {taxon_name} already exists.")
            return
        with open(f"{self.analysis_loc}/seq/organisms/stock_ids/{filename}", "wb") as proteome:
            proteome.write(requests_object.content)
            print(f"{filename} saved in /seq/organisms/stock_ids/")
        # Rename seq with incrementing number
        print("Renaming seq in file")
        with open(f"{self.analysis_loc}/seq/organisms/stock_ids/{filename}") as to_correct, open(
                f"{self.analysis_loc}/seq/organisms/renamed/{filename}", "w") as corrected:
            list_of_records = []
            all_seq_number = 0
            accepted_seq_number = 0
            for record in SeqIO.parse(to_correct, "fasta"):
                all_seq_number += 1
                if min_seq_len < len(record.seq) < max_seq_len:
                    accepted_seq_number += 1  # Implementing counter this due to generator
                    record.id = f"{short_taxon_name}_{accepted_seq_number}"
                    record.description = f"{short_taxon_name}_{accepted_seq_number}"
                    list_of_records.append(record)
            SeqIO.write(list_of_records, corrected, 'fasta')
        print(
            f"Discarded {all_seq_number - accepted_seq_number} out of {all_seq_number} sequences. {accepted_seq_number} sequences left.")
        print("Renaming completed. \nFile saved in seq/organisms/renamed/")

    def concat_files(self, extension: str, path_to_files: str, out_path: str):
        """
        Method used to combine multiple files, enabling to choose extension of file

        :param extension: extension of file
        :param path_to_files: Path to directory with input files
        :param out_path: path to output file inside working directory, where file will be saved
        :return: None, saved concatenated file in chosen path
        """
        infile_paths = f"{self.analysis_loc}/{path_to_files}"
        outfile_path = f"{self.analysis_loc}/{out_path}.{extension}"
        process = ['for f in ', infile_paths, '/*', '; do cat $f >> ', outfile_path, "; done"]
        # subprocess.run(process, shell=True, check=True)
        os.system("".join(process))

    def concat_trees(self, path_to_files: str, out_path: str):
        """
        Method used to combine multiple tree files, created specifically for  combining .nwk files

        :param path_to_files: Path to directory with input files
        :param out_path: path to output file inside working directory, where file will be saved
        :return: None, saved concatenated file in chosen path
        """
        process = ['for f in ', path_to_files, '/*', '; do (cat "${f}"; echo) >> ', out_path,
                   "; done"]
        # subprocess.run(process, shell=True, check=True)
        os.system("".join(process))

    def perform_mmseq_clustering(self, path_to_concat_fasta: str):
        """
        Method performing clustering using mmseq2 easy-cluster method with fixed values of parameters

        :param path_to_concat_fasta: Path to concatenated fastafile with all sequences
        :return: None, saved 3 cluster files which are used in future analysis
        """
        process = ["mmseqs easy-cluster ", self.analysis_loc, "/", path_to_concat_fasta,
                   " clusterRes tmp --min-seq-id 0.5 -c 0.8 --cov-mode 1"]
        #  mmseqs easy-cluster proper.fasta clusterRes tmp --min-seq-id 0.5 -c 0.8 --cov-mode 0
        # subprocess.run(process, shell=True, check=True)
        os.system("".join(process))

    def split_clusters(self):
        """
        Method applies small filter - min and max number of proteins in cluster, discards clasters,
        that consist only of one species.

        :return:
        """
        # Sidenote for author - Had to do precleaning step due to 108360 clusters, after precleaning 7322 clusters
        print("Splitting cluster file.")
        with open(f"clusterRes_all_seqs.fasta", "r") as cluster_file:
            clusters_to_save = []
            clusters = re.split(">.+\n>", cluster_file.read())
            cluster_counter = 0
            for cluster in clusters:
                temp_cluster = f">{cluster}"
                if len(set(re.findall(">[A-Z]{6}", temp_cluster))) > 1:
                    with open(
                            f"{self.analysis_loc}/cluster/all/cluster_{str(cluster_counter).zfill(5)}",
                            "w") as file_to_write:
                        file_to_write.write(temp_cluster)
                    cluster_counter += 1

    def fold_rmdup_cluster(self, path_to_cluster, min_size_cluster: int = 1,
                           max_size_cluster: int = 100):
        """
        Method used to fold single line fasta files and discard clusters of sizes which are out of
        chosen ranges

        :param path_to_cluster: Path to file with clusters
        :param min_size_cluster: minimum number of proteins in cluster
        :param max_size_cluster: maximum number of proteins in cluster
        :return: Splitted clusters in cluster/all/ dir
        """
        if path_to_cluster.split(".")[-1] != "fasta":
            process_one = f"seqkit seq -w 60 {path_to_cluster} > {path_to_cluster}_folded"
            os.system(process_one)
            os.system(f"rm {path_to_cluster}")
            process_two = f"seqkit rmdup --quiet -s < {path_to_cluster}_folded > {path_to_cluster}_nodup.fasta"
            os.system(process_two)
            os.system(f"rm {path_to_cluster}_folded")
            with open(f"{path_to_cluster}_nodup.fasta", "r") as to_read, open(
                    f"{path_to_cluster}.fasta", "w") as to_write:
                temp_cluster = to_read.read()
                if min_size_cluster < temp_cluster.count(">") < max_size_cluster:
                    to_write.write(temp_cluster)
            os.system(f"rm {path_to_cluster}_nodup.fasta")



    def remove_paralogs(self, path_to_cluster: str, out_dir_path: str):
        """
        Method used to remove paralogs in order to obtain one-to-one clusters
        :param path_to_cluster: Path to single cluster
        :param out_dir_path: output directory path
        :return:
        """
        list_of_sequences = []
        list_of_taxons = []
        cluster_name = path_to_cluster.split("/")[-1]
        for number, record in enumerate(SeqIO.parse(path_to_cluster, "fasta")):
            taxon_acronym = record.id[:6]
            if taxon_acronym not in list_of_taxons:
                list_of_sequences.append(record)
                list_of_taxons.append(taxon_acronym)
        if len(list_of_taxons) == len(self.taxid_dict):
            SeqIO.write(list_of_sequences, f"{out_dir_path}/{cluster_name}", "fasta")

    def make_alignment(self, path_to_cluster: str, out_dir_path: str):
        """
        Method using Muscle MAFFT plugin to perform aligning of clusters with sequences

        :param path_to_cluster: Path to chosen cluster
        :param out_dir_path: Output dir path, where alignment will be saved
        :return: None, fasta alignment
        """
        # Align sequences using MAFFT
        cluster_name = path_to_cluster.split("/")[-1]
        align = MuscleCommandline(input=path_to_cluster, out=f"{out_dir_path}/{cluster_name}",
                                  quiet=True)
        subprocess.run(str(align).split())

    def infer_ft_tree(self, path_to_aligned_file: str, out_dir_path: str):
        """
        Method used to infer trees using FastTree tool

        :param path_to_aligned_file: Path to alignment
        :param out_dir_path:  Path to output directory where tree will be stored
        :return: None, saving file in chosen dir
        """
        cluster_name = path_to_aligned_file.split("/")[-1]
        cluster_name_without_extension = cluster_name.split(".")[0]
        cluster_number = cluster_name_without_extension.split("_")[-1]
        os.system(
            f"FastTree -quiet {path_to_aligned_file} > {out_dir_path}/ft_tree_{cluster_number}.nwk")

    def filter_bootstrapped(self, path_to_tree: str, out_dir: str, threshold: float = 0.8):
        """
        Method used to filter trees with bootstrap value lower than chosen.

        :param path_to_tree: Path to .nwk tree file
        :param out_dir: Output directory where file will be saved
        :param threshold: Bootstrap treshold
        :return: None, trees saved in directory if they had greater value than bootstrap threshold
        """
        with open(path_to_tree) as tree_file:
            filename = path_to_tree.split("/")[-1]
            t = ete3.Tree(tree_file.read())
            support_list = []
            for node in t.traverse():
                support_list.append(node.support)
            if statistics.mean(support_list) > threshold:
                t.write(outfile=f"{out_dir}/{filename}")

    def reformat_nwk(self, path_to_tree: str, out_dir: str, nwk_format: int = 0):
        """
        Method used to reformat newick file (chaning species names, resolving polytomy) in inner
        processing of trees

        :param path_to_tree: Path tro tree file
        :param out_dir: Path to output directory
        :param nwk_format: for supertree use 9, for consensus use 0
        :return: None, reformated tree file
        """
        with open(path_to_tree) as tree_file:
            filename = path_to_tree.split("/")[-1]
            t = ete3.Tree(tree_file.read())
            for node in t.traverse():
                node.name = node.name[:6]
            t.resolve_polytomy(recursive=True)
            try:
                t.unroot()
            except:
                pass
            t.write(format=nwk_format, outfile=f"{out_dir}/{filename}")

    def rm_semicolons_from_nwk(self, path_to_file: str):
        """
        Method used to remove semicolons from .nwk file. Some of tools demand .nwk trees without
        semicolon at the end to work properly.

        :param path_to_file: Path to .nwk tree file
        :return: Modified input .nwk file without semicolons at the end of line
        """
        os.system(f'sed -i "s/;//" {path_to_file}')

    def perform_consensus_calc(self, path_input_trees: str, path_output_dir: str,
                               strictnes: float = 0.3):
        """
        Method used to calculate consensus trees, by default strictness is set to 0.3

        :param path_input_trees: Path to input file with concatenated .nwk trees, each in newline
                                 with semicolon at the end
        :param path_output_dir: Path to output directory
        :param strictnes: 0.0 - 1.0 scale strictness, the lower calue, the more binary divisions will be seen
        :return: Computed consensus tree saved in file.
        """
        tree_name_with_ext = path_input_trees.split("/")[-1]
        tree_name_no_ext = tree_name_with_ext.split(".")[0]
        tree_name_out = tree_name_no_ext + ".nwk"
        trees = list(Bio.Phylo.parse(path_input_trees, "newick"))
        majority_tree = Bio.Phylo.Consensus.majority_consensus(trees, 0.3)
        Bio.Phylo.write(majority_tree, f"{path_output_dir}/{tree_name_out}", "newick")

    def fu_to_nwk(self, out_dir: str):
        """
        Method used to extract from fasturec output file tree with best score and convert fasturec
        this tree to regular newick format (with semicolons at the end)

        :param out_dir:
        :return: None, saved in file supertree with best score
        """
        fu_file_path_list = [fu_file_path for fu_file_path in
                             glob.glob(f"./fu.*")]
        for number, file in enumerate(fu_file_path_list):
            with open(file, "r") as fu_file:
                best_tree = fu_file.readline()

            out_filename = file.split("/")[-1]
            out_filename = out_filename.split(".")[1]

            best_tree = best_tree.split(" ")[-1].strip("\n") + ";"
            with open(f"{out_dir}/{out_filename}.nwk", "w") as to_write:
                to_write.write(best_tree)
            os.system(f"mv fu.* {self.path_archive}/")

    def perform_supertree_calc(self, path_to_concat_trees, out_dir):
        """
        Method used to perform supertree calculations using fasturec2 software

        :param path_to_concat_trees: Path to input file with concatenated .nwk trees, each in newline without semicolon at the end
        :param out_dir: Output directory
        :return: None, calculated supertree
        """
        os.system(f"fasturec2/fasturec2 -Y -G {path_to_concat_trees}")
        self.fu_to_nwk(out_dir)

    def calculate_rf_all_trees(self, trees_dir_path: str, known_topology_tree: str):
        """
        Method used to calculate rfs of all trees and comparing it with 'known_topology.nwk' file
        which needs to be in the directory where input file and script is run.

        :param trees_dir_path: Path to directory with computed trees.
        :param known_topology_tree: Path to known topology newick file
        :return:
        """
        trees_to_compare_paths_list = [tree_path for tree_path in
                                       glob.glob(f"{trees_dir_path}/*.nwk")]
        trees_to_compare_paths_list.append(known_topology_tree)
        trees_names_list = [tree_name.split("/")[-1] for tree_name in trees_to_compare_paths_list]
        trees_path_to_name_dict = dict(zip(trees_to_compare_paths_list, trees_names_list))
        df = pd.DataFrame(columns=trees_names_list, index=trees_names_list)
        for row_tree_path in trees_to_compare_paths_list:
            with open(row_tree_path) as tree_from_row:
                row_tree = ete3.Tree(tree_from_row.read())
            row_name = trees_path_to_name_dict[row_tree_path]
            for column_tree_path in trees_to_compare_paths_list:
                with open(column_tree_path) as tree_from_column:
                    column_tree = ete3.Tree(tree_from_column.read())
                column_name = trees_path_to_name_dict[column_tree_path]
                rf, rf_max, x1, x2, x3, x4, x5 = row_tree.robinson_foulds(column_tree,
                                                                          unrooted_trees=True)
                df.at[row_name, column_name] = rf
        # Removing already existing files
        if os.path.exists("all_trees_rf.html"):
            os.remove("all_trees_rf.html")
        if os.path.exists("all_trees_rf"):
            os.remove("all_trees_rf")
        if os.path.exists("all_trees_rf.png"):
            os.remove("all_trees_rf.png")

        df.to_html('all_trees_rf.html')
        subprocess.call(
            'wkhtmltoimage -f png --width 0 all_trees_rf.html all_trees_rf.png', shell=True)
        print("Dataframe saved a picture 'all_trees_rf.png'")
        return df

    def draw_trees_to_png(self, trees_dir_path: str, out_dir: str):
        """
        Method uses ete3 Tree object to render png pictures of trees. In order to function properly
        in directory where script runs it's obligatory to put newick file named "known_topology.nwk".

        :param trees_dir_path: Path to directory with computed trees
        :param out_dir: Path to output directory
        :return:
        """
        trees_to_draw = [tree_path for tree_path in glob.glob(f"{trees_dir_path}/*.nwk")]
        try:
            trees_to_draw.append("known_topology.nwk")
        except:
            pass
        for tree in trees_to_draw:
            out_filename = tree.split("/")[-1]
            out_filename = out_filename.split(".")[0]
            with open(tree) as tree_file:
                t = ete3.Tree(tree_file.read())
            tree_name = tree.split("/")[-1]
            tree_name = tree_name.split(".")[-2]
            t.render(f"{out_dir}/{out_filename}.png", w=500, units="mm")

    def draw_trees_to_ascii(self, trees_dir_path: str, out_dir: str):
        """
        Method uses Bio.Phylo Tree object to create ascii representation of trees. In order to function properly
        in directory where script runs it's obligatory to put newick file named "known_topology.nwk".

        :param trees_dir_path:
        :return:
        """
        trees_to_draw = [tree_path for tree_path in glob.glob(f"{trees_dir_path}/*.nwk")]
        try:
            trees_to_draw.append("known_topology.nwk")
        except:
            pass
        for tree in trees_to_draw:
            out_filename = tree.split("/")[-1]
            out_filename = out_filename.split(".")[0]
            known_topology_tree = Bio.Phylo.read(tree, "newick")
            with open(f"{out_dir}/{out_filename}.txt", "w") as k_t:
                Bio.Phylo.draw_ascii(known_topology_tree, k_t)
