# Replicating BLOBFISH Paper Results

## GTEx Tissue Analysis

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

1. Run ensemblTFsAndDescriptions.R to add ENSEMBL IDs and descriptive gene names to BLOBFISH results, specifying the following variables:

    -   **blobfishDir:** The directory where the BLOBFISH results are stored.
    -   **motifFile:** The directory containing the motif prior file from GRAND.

## Comparison to Genes2Networks

1. Run prepareOutputsForGenes2Networks.R.
1. Clone the Genes2Networks repository using the command `git clone git@github.com:MaayanLab/Genes2Networks.git`.
2. [Install Java](https://www.java.com/en/download/).
3. [Install Gradle](https://docs.gradle.org/current/userguide/installation.html).
4. Open the `build.gradle` file in the repository where you cloned Genes2Networks and replace the following lines:
     - `jcenter()` with `mavenCentral()`
     - `apply plugin: 'maven'` with `apply plugin: 'maven-publish'`
     - `sourceCompatibility = 1.7` with `sourceCompatibility = JavaVersion.VERSION_1_8` 
     - `compile 'com.github.MaayanLab.common:common-core:master-SNAPSHOT'` with `implementation 'com.github.MaayanLab.common:common-core:master-SNAPSHOT'`
     - `compile 'com.github.MaayanLab.common:common-swing:master-SNAPSHOT'` with `implementation 'com.github.MaayanLab.common:common-swing:master-SNAPSHOT'`
     - `compile 'com.github.MaayanLab.common:common-graph:master-SNAPSHOT'` with `implementation 'com.github.MaayanLab.common:common-graph:master-SNAPSHOT'`
     - `compile group: 'org.tinyjee.jgraphx', name: 'jgraphx', version:'1.10.1.3'` with `implementation group: 'org.tinyjee.jgraphx', name: 'jgraphx', version:'1.10.1.3'`
     - `testCompile group: 'junit', name: 'junit', version:'3.8.1'` with `testImplementation 'junit:junit:4.13.2'`
     - This block:
       ```java
          task sourcesJar(type: Jar, dependsOn: classes) {
            classifier = 'sources'
            from sourceSets.main.allSource
          }
       ```
       with this block:
       ```java
          tasks.register('sourcesJar', Jar) {
            dependsOn classes
            archiveClassifier.set('sources')
            from sourceSets.main.allSource
          }
       ```
     - `classpath 'de.undercouch:gradle-download-task:3.4.3'` with `classpath 'de.undercouch:gradle-download-task:5.6.0'`
5. Add the following lines to the bottom of `build.gradle`:
  ```java
     compileJava.dependsOn downloadAndUnzipResources
     processTestResources.dependsOn downloadAndUnzipResources
  ```
5. Add the following lines after `apply plugin: 'de.undercouch.download'`:
  ```java
     apply plugin: 'application'

     application {
       mainClass = 'edu.mssm.pharm.maayanlab.Genes2Networks.Genes2Networks'
     }
  ```
5. Run `gradle run --args="<input_dir>/inputGenesGenes2Network.txt <output_dir>/skinOut.sig <input_dir>/skin.sig"`, and the analogous runs for `subcutaneous_adipose.sig`, `skeletal_muscle.sig`, `lung.sig`, and `aorta.sig`, where <input_dir> is the directory containing the files generated by prepareOutputsForGenes2Networks.R, and <output_dir> is an existing directory where you wish to store the output subnetworks.
6. Run analysisGenes2Networks.R, modifying the following:

    -   **indir:** The directory where your Genes2Networks result files are saved.
    -   **outdir:** The directory where you wish to save your results from the analysis of Genes2Networks.
    -   **blobfishDir:** The directory where the BLOBFISH results are stored.

## Comparison to SP
1. Run prepareOutputsForGenes2Networks.R.
2. Run analysisPipelineSP.R.


## eQTL Analysis

## Sparse Network Analysis
       


