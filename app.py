#!/usr/bin/env python3

from PhyloAnalyser import PhyloAnalyser
import glob
import joblib


def main():
    analysis_tool = PhyloAnalyser("taxid_to_analyse.txt", "test_dir")
    taxon_dictionary = analysis_tool.read_taxid_input()

    """
    for taxid in taxon_dictionary.keys():
        analysis_tool.download_proteome(taxid)

    # Concatenate renamed .fasta files into one for clustering
    analysis_tool.concat_files(extension="fasta", path_to_files="test_dir/seq/organisms/renamed",
                               out_path="test_dir/seq/organisms/concatenated/concat_renamed")

    Performing clustering of concatenated .fasta file
    analysis_tool.perform_mmseq_clustering(concatenated_fasta_path="test_dir/seq/organisms/concatenated/concat_renamed.fasta")

    analysis_tool.split_clusters()

    path_to_all_dir = f"{analysis_tool.analysis_loc}/cluster/all"
    clusters_paths_in_all = [cluster_path for cluster_path in
                             glob.glob(f"{path_to_all_dir}/cluster*")]

    # print("Starting folding clusters")
    joblib.Parallel(n_jobs=7)(
        joblib.delayed(analysis_tool.fold_cluster)(cluster_path) for cluster_path in
        clusters_paths_in_all)
    print(f"Clusters saved in in /cluster/all directory.")

    path_to_filtered_dir = f"{analysis_tool.analysis_loc}/cluster/filtered"
    # actually, we could use list_of_clusters_paths, since we are taking it as input in below
    # parallel job, saved it as new variable because using list_of_clusters_paths could be misleading

    # Removing paralogs
    print("Starting removing paralogs from clusters.")
    joblib.Parallel(n_jobs=7)(
        joblib.delayed(analysis_tool.remove_paralogs)(cluster_path, path_to_filtered_dir) for
        cluster_path in clusters_paths_in_all)
    print(f"Clusters without paralogs saved in {path_to_filtered_dir}")

    list_of_filtered_clusters_paths = [cluster_path for cluster_path in
                                       glob.glob(f"{path_to_filtered_dir}/cluster*")]
    """
    path_to_all_alignments_dir = f"{analysis_tool.analysis_loc}/alignment/all"
    """
    # Aligning all sequences
    print("Starting aligning clusters from all.")
    joblib.Parallel(n_jobs=7)(
        joblib.delayed(analysis_tool.make_alignment)(cluster_path, path_to_all_alignments_dir) for
        cluster_path in clusters_paths_in_all)
    print(f"Clusters without paraloogs saved in {path_to_all_alignments_dir}")

    clusters_paths_in_filtered = [cluster_path for cluster_path in
                                  glob.glob(f"{path_to_filtered_dir}/cluster*")]
    """
    path_to_filtered_alignments_dir = f"{analysis_tool.analysis_loc}/alignment/filtered"
    """
    # Aligning filtered sequences
    print("Starting aligning clusters from filtered.")
    joblib.Parallel(n_jobs=7)(
        joblib.delayed(analysis_tool.make_alignment)(cluster_path, path_to_filtered_alignments_dir) for
        cluster_path in list_of_filtered_clusters_paths)
    print(f"Clusters without paralogs saved in {path_to_filtered_alignments_dir}")
    
    """
    path_to_all_alignments_dir_updated = f"{analysis_tool.analysis_loc}/alignment/all"
    list_of_all_alignments_paths = [alignment_path for alignment_path in
                                    glob.glob(f"{path_to_all_alignments_dir_updated}/cluster*")]
    list_of_filtered_alignments_paths = [alignment_path for alignment_path in
                                         glob.glob(f"{path_to_filtered_alignments_dir}/cluster*")]

    path_to_all_tree_dir = f"{analysis_tool.analysis_loc}/tree/ft/all"
    path_to_filtered_tree_dir = f"{analysis_tool.analysis_loc}/tree/ft/filtered"

    # Infer all trees
    print("Starting inferring trees from all.")
    joblib.Parallel(n_jobs=7)(
        joblib.delayed(analysis_tool.infer_ft_tree)(aligned_file_path, path_to_all_tree_dir) for
        aligned_file_path in list_of_all_alignments_paths)
    print(f"Trees saved in trees/all directory")

    # Infer filtered trees
    print("Starting inferring trees from filtered.")
    joblib.Parallel(n_jobs=7)(
        joblib.delayed(analysis_tool.infer_ft_tree)(aligned_file_path, path_to_filtered_tree_dir)
        for
        aligned_file_path in list_of_filtered_alignments_paths)
    print(f"Trees saved in trees/all directory")


    # To do: napsiac funkcje, ktora usuwa te same sekwencje z all sekwencji i od nowa zapuscic calosc

if __name__ == '__main__':
    main()