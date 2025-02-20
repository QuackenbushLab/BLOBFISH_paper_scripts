# Replicating BLOBFISH Paper Results

1.  Download the PPI file, the motif file, the expression data, and the LIONESS networks for males ages 20-29 in each tissue from GRAND: https://grand.networkmedicine.org/tissues/.
2.  Download the curated gene set GMT file from Human MSigDB: https://www.gsea-msigdb.org/gsea/msigdb/collections.jsp.
3.  Run GenerateNullPANDA.R to generate the null models, specifying the following variables:

    -   **ppiFilePath:** The path to the file containing the protein-protein interaction network.
    -   **motifFilePath:** The path to the file containing the transcription factor binding motif.
    -   **nullFilePath:** The path to the directory where you wish to store the null model.

2.  Run analysisPipeline.R to run BLOBFISH on the gene subsets and to generate plots and perform pathway enrichment analysis results, specifying the following variables:

    -   **dir_subcutaneous_adipose:** The directory containing the subcutaneous adipose sample networks.
    -   **dir_skeletal_muscle:** The directory containing the skeletal muscle sample networks.
    -   **dir_skin:** The directory containing the skin sample networks.
    -   **dir_lung:** The directory containing the lung sample networks.
    -   **dir_aorta:** The directory containing the aorta sample networks.
    -   **dir_expression:** The directory containing the gene expression data.
    -   **outdir:** The directory where you wish to save your BLOBFISH networks, plots, and pathway enrichment analysis results.
    -   **gmt_pathway_file:** The path to the file where you have stored the GMT pathway file to use with FGSEA (e.g., c2.cp.v2023.2.Hs.symbols.gmt).
    -   **null_file:** The path to the file where the null PANDA distribution is stored.
    -   **null_output_file**: The path to the file where you wish to store the sampled PANDA distribution to use for BLOBFISH.

1.  Run generate_rand_sets.R, specifying the following variables:

    -   **dir_input:** The directory containing the compiled sample-specific networks in RDS format, generated using the previous step.
    -   **output_file:** The directory where you wish to save the random subsets of genes.

1.  Run runBLOBFISH_parallel.sh, specifying the following variables in single_blobfish_run.R:

    -   **outdir:** The directory where you wish to save your BLOBFISH results.
    -   **null_file:** The path to the file where the BLOBFISH null distribution is stored.
    -   **tissue_dir:** The directory containing the compiled sample-specific networks in RDS format.
    -   **randset_file:** The directory containing the random subsets of genes.
    -   **randoutputs:** The directory where you wish to save your BLOBFISH results (should be identical to **outdir**) in single_blobfish_run.R.
    -   **pvalsdir:** The directory where you wish to store the p-values from the BLOBFISH runs.

1.  Run tfcounts.R, specifying the following variables:

    -   **dir_randoutputs:** The directory where the BLOBFISH results are stored.
    -   **output_file:** The directory where you wish to save the random subsets of genes.
    -   **randset_file:** The directory containing the random subsets of genes.



