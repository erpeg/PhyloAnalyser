#!/usr/bin/env python3

import os
from typing import Dict

from Bio import SeqIO
import requests
import ete3
import re
import glob
from Bio.Align.Applications import MuscleCommandline
import statistics

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
        self.path_ft_all_inferred = f"{self.analysis_loc}/tree/ft/all/inferred"
        self.path_ft_all_bootstrapped = f"{self.analysis_loc}/tree/ft/all/bootstrapped"

        self.path_ft_all_boot_done_consensus = f"{self.analysis_loc}/tree/ft/all/boot_done/consensus"
        self.path_ft_all_boot_nondone_consensus = f"{self.analysis_loc}/tree/ft/all/boot_nondone/consensus"
        self.path_ft_all_boot_done_supertree = f"{self.analysis_loc}/tree/ft/all/boot_done/supertree"
        self.path_ft_all_boot_nondone_supertree = f"{self.analysis_loc}/tree/ft/all/boot_nondone/supertree"

        self.path_ft_filtered_inferred = f"{self.analysis_loc}/tree/ft/filtered/inferred"
        self.path_ft_filtered_bootstrapped = f"{self.analysis_loc}/tree/ft/filtered/bootstrapped"

        self.path_ft_filtered_boot_done_consensus = f"{self.analysis_loc}/tree/ft/filtered/boot_done/consensus"
        self.path_ft_filtered_boot_nondone_consensus = f"{self.analysis_loc}/tree/ft/filtered/boot_nondone/consensus"
        self.path_ft_filtered_boot_done_supertree = f"{self.analysis_loc}/tree/ft/filtered/boot_done/supertree"
        self.path_ft_filtered_boot_nondone_supertree = f"{self.analysis_loc}/tree/ft/filtered/boot_nondone/supertree"

        self.path_ft_tree_concat = f"{self.analysis_loc}/tree/concatenated"

        self.list_of_consensus_dirs = [self.path_ft_all_boot_done_consensus,
                                       self.path_ft_all_boot_nondone_consensus,
                                       self.path_ft_filtered_boot_done_consensus,
                                       self.path_ft_filtered_boot_nondone_consensus]

        self.list_of_supertree_dirs = [self.path_ft_all_boot_done_supertree,
                                       self.path_ft_all_boot_nondone_supertree,
                                       self.path_ft_filtered_boot_done_supertree,
                                       self.path_ft_filtered_boot_nondone_supertree]

        self.list_of_tree_dirs = [self.path_ft_all_inferred,
                                  self.path_ft_all_bootstrapped,
                                  self.path_ft_all_boot_done_consensus,
                                  self.path_ft_all_boot_nondone_consensus,
                                  self.path_ft_all_boot_done_supertree,
                                  self.path_ft_all_boot_nondone_supertree,
                                  self.path_ft_filtered_inferred,
                                  self.path_ft_filtered_bootstrapped,
                                  self.path_ft_filtered_boot_done_consensus,
                                  self.path_ft_filtered_boot_nondone_consensus,
                                  self.path_ft_filtered_boot_done_supertree,
                                  self.path_ft_filtered_boot_nondone_supertree,
                                  self.path_ft_tree_concat]

        self.setup_dirs()  # needs to be run after self.analysis_loc

    def setup_dirs(self) -> None:
        """Method that sets up necessary directories for whole project to work"""
        if self.analysis_loc[-1] == "/":
            self.analysis_loc = self.analysis_loc[:-1]
        # Creating "downloads" directory
        if not os.path.exists(f"{self.analysis_loc}/downloads"):
            os.makedirs(f"{self.analysis_loc}/downloads")
        # Creating "tree" directories
        if not os.path.exists(f"{self.analysis_loc}/tree/ft/cleaned"):
            os.makedirs(f"{self.analysis_loc}/tree/ft/cleaned")
        if not os.path.exists(f"{self.analysis_loc}/tree/ft/all"):
            os.makedirs(f"{self.analysis_loc}/tree/ft/all")
        if not os.path.exists(f"{self.analysis_loc}/tree/ft/filtered"):
            os.makedirs(f"{self.analysis_loc}/tree/ft/filtered")
        if not os.path.exists(f"{self.analysis_loc}/tree/concatenated"):
            os.makedirs(f"{self.analysis_loc}/tree/concatenated")
        for dir_path in self.list_of_tree_dirs:
            if not os.path.exists(dir_path):
                os.makedirs(dir_path)
        # Creating "seq" directories
        if not os.path.exists(f"{self.analysis_loc}/seq/organisms/stock_ids"):
            os.makedirs(f"{self.analysis_loc}/seq/organisms/stock_ids")
        if not os.path.exists(f"{self.analysis_loc}/seq/organisms/renamed"):
            os.makedirs(f"{self.analysis_loc}/seq/organisms/renamed")
        if not os.path.exists(f"{self.analysis_loc}/seq/organisms/concatenated"):
            os.makedirs(f"{self.analysis_loc}/seq/organisms/concatenated")
        if not os.path.exists(f"{self.analysis_loc}/seq/aligned"):
            os.makedirs(f"{self.analysis_loc}/seq/aligned")
        # Creating "cluster" directories:
        if not os.path.exists(f"{self.analysis_loc}/cluster/all"):  # used for supertree
            os.makedirs(f"{self.analysis_loc}/cluster/all")
        if not os.path.exists(f"{self.analysis_loc}/cluster/cleaned"):  # used for supertree
            os.makedirs(f"{self.analysis_loc}/cluster/cleaned")
        if not os.path.exists(f"{self.analysis_loc}/cluster/filtered"):  # after filtering
            os.makedirs(f"{self.analysis_loc}/cluster/filtered")
        # Creating "alignment" directories
        if not os.path.exists(f"{self.analysis_loc}/alignment/cleaned"):  # used for supertree
            os.makedirs(f"{self.analysis_loc}/alignment/cleaned")
        if not os.path.exists(f"{self.analysis_loc}/alignment/all"):  # used for supertree
            os.makedirs(f"{self.analysis_loc}/alignment/all")
        if not os.path.exists(f"{self.analysis_loc}/alignment/filtered"):  # after filtering
            os.makedirs(f"{self.analysis_loc}/alignment/filtered")

    def read_taxid_input(self) -> Dict:
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

    def concat_files(self, extension, path_to_files, out_path):
        infile_paths = f"{self.analysis_loc}/{path_to_files}"
        outfile_path = f"{self.analysis_loc}/{out_path}.{extension}"
        process = ['for f in ', infile_paths, '/*', '; do cat $f >> ', outfile_path, "; done"]
        # subprocess.run(process, shell=True, check=True)
        os.system("".join(process))

    def concat_trees(self, path_to_files, out_path):
        process = ['for f in ', path_to_files, '/*', '; do (cat "${f}"; echo) >> ', out_path, "; done"]
        # subprocess.run(process, shell=True, check=True)
        os.system("".join(process))

    def perform_mmseq_clustering(self, path_to_concat_fasta):
        process = ["mmseqs easy-cluster ", self.analysis_loc, "/", path_to_concat_fasta,
                   " clusterRes tmp --min-seq-id 0.5 -c 0.8 --cov-mode 1"]
        #  mmseqs easy-cluster proper.fasta clusterRes tmp --min-seq-id 0.5 -c 0.8 --cov-mode 0
        # subprocess.run(process, shell=True, check=True)
        os.system("".join(process))

    def split_clusters(self):
        """
        Method applies small filter - min and max number of proteins in cluster, discards clasters,
        that consist only of one species.


        :param min_size_cluster:
        :param max_size_cluster:
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

    def remove_same_seqs(self, path_to_cluster):
        name_of_cluster = path_to_cluster.split("/")[-1]
        with open(path_to_cluster) as to_correct, open(
                f"{self.analysis_loc}/cluster/cleaned/{name_of_cluster}", "w") as corrected:
            list_of_sequences = []
            list_of_taxons = []
            for record in SeqIO.parse(to_correct, "fasta"):
                if record.id not in list_of_taxons:
                    list_of_sequences.append(record)
                    list_of_taxons.append(record.id)
                SeqIO.write(list_of_sequences, corrected, 'fasta')

    def remove_paralogs(self, path_to_cluster, out_dir_path):
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

    def make_alignment(self, path_to_cluster, out_dir_path):
        # Align sequences using MAFFT
        cluster_name = path_to_cluster.split("/")[-1]
        align = MuscleCommandline(input=path_to_cluster, out=f"{out_dir_path}/{cluster_name}",
                                  quiet=True)
        subprocess.run(str(align).split())

    def infer_ft_tree(self, path_to_aligned_file, out_dir_path):
        cluster_name = path_to_aligned_file.split("/")[-1]
        cluster_name_without_extension = cluster_name.split(".")[0]
        cluster_number = cluster_name_without_extension.split("_")[-1]
        os.system(
            f"FastTree -quiet {path_to_aligned_file} > {out_dir_path}/ft_tree_{cluster_number}.nwk")

    def filter_bootstrapped(self, path_to_tree: str, out_dir: str, threshold: float = 0.8):
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

        :param path_to_tree:
        :param out_dir:
        :param nwk_format: for supertree use 9, for consensus use 0
        :return:
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

    def rm_semicolons_from_nwk(self, path_to_file):
        os.system(f'sed -i "s/;//" {path_to_file}')
