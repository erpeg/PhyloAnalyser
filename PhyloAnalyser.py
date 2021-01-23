#!/usr/bin/env python3

import os
from typing import Dict

from Bio import SeqIO
import requests
import ete3
import re
import glob
from Bio.Align.Applications import MuscleCommandline


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
        self.setup_dirs()  # needs to be run after self.analysis_loc

    def setup_dirs(self) -> None:
        """Method that sets up necessary directories for whole project to work"""
        if self.analysis_loc[-1] == "/":
            self.analysis_loc = self.analysis_loc[:-1]
        # Creating "downloads" directory
        if not os.path.exists(f"{self.analysis_loc}/downloads"):
            os.makedirs(f"{self.analysis_loc}/downloads")
        # Creating "tree" directories
        if not os.path.exists(f"{self.analysis_loc}/tree/ft/all"):
            os.makedirs(f"{self.analysis_loc}/tree/ft/all")
        if not os.path.exists(f"{self.analysis_loc}/tree/ft/filtered"):
            os.makedirs(f"{self.analysis_loc}/tree/ft/filtered")
        if not os.path.exists(f"{self.analysis_loc}/tree/concatenated"):
            os.makedirs(f"{self.analysis_loc}/tree/concatenated")
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
        if not os.path.exists(f"{self.analysis_loc}/cluster/filtered"):  # after filtering
            os.makedirs(f"{self.analysis_loc}/cluster/filtered")
        # Creating "alignment" directories
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
        outfile_path = f"{self.analysis_loc}/{out_path}.{extension}"
        process = ['for f in ', path_to_files, '/*.', extension, '; do cat $f >> ', outfile_path, "; done"]
        # subprocess.run(process, shell=True, check=True)
        os.system("".join(process))

    def perform_mmseq_clustering(self, path_to_concat_fasta):
        process = ["mmseqs easy-cluster ", path_to_concat_fasta, " clusterRes tmp --min-seq-id 0.5 -c 0.8 --cov-mode 0"]
        #  mmseqs easy-cluster proper.fasta clusterRes tmp --min-seq-id 0.5 -c 0.8 --cov-mode 0
        # subprocess.run(process, shell=True, check=True)
        os.system("".join(process))

    def split_clusters(self, min_size_cluster: int = 3, max_size_cluster: int = 60):
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
                if min_size_cluster < cluster.count(">") < max_size_cluster:
                    if len(set(re.findall(">[A-Z]{6}", temp_cluster))) > 1:
                        with open(f"{self.analysis_loc}/cluster/all/cluster_{str(cluster_counter).zfill(5)}", "w") as file_to_write:
                            file_to_write.write(temp_cluster)
                        cluster_counter += 1

    def fold_cluster(self, path_to_cluster):
        if path_to_cluster.split(".")[-1] != "fasta":
            process = f"seqkit seq -w 60 {path_to_cluster} > {path_to_cluster}.fasta"
            os.system(process)
            os.system(f"rm {path_to_cluster}")


    def remove_paralogs(self, path_to_cluster, out_dir_path):
        list_of_sequences = []
        list_of_taxons = []
        cluster_name = path_to_cluster.split("/")[-1]
        for number, record in enumerate(SeqIO.parse(path_to_cluster, "fasta")):
            taxon_acronym = record.id[:6]
            if taxon_acronym not in list_of_taxons:
                list_of_sequences.append(record)
                list_of_taxons.append(taxon_acronym)
        SeqIO.write(list_of_sequences, f"{out_dir_path}/{cluster_name}", "fasta")

    def make_alignment(self, path_to_cluster, out_dir_path):
        # Align sequences using MAFFT
        cluster_name = path_to_cluster.split("/")[-1]
        align = MuscleCommandline(input=path_to_cluster, out=f"{out_dir_path}/{cluster_name}", quiet=True)
        subprocess.run(str(align).split())

    def infer_ft_tree(self, path_to_aligned_file, out_dir_path):
        cluster_name = path_to_aligned_file.split("/")[-1]
        cluster_name_without_extension = cluster_name.split(".")[0]
        cluster_number = cluster_name_without_extension.split("_")[-1]
        os.system(
            f"FastTree -quiet {path_to_aligned_file} > {out_dir_path}/ft_tree_{cluster_number}.nwk")




