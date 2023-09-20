# IRF4-in-Treg-and-Th17
Repository containing the scripts for data processing of our Paper "Unveiling IRF4-steered regulation of context-dependent effector programs in Th17 and Treg cells"


## ABSTRACT

The transcription factor interferon regulatory factor 4 (IRF4) is crucial for the differentiation and fate determination of pro-inflammatory T helper (Th)17 and the functionally opposing group of immunomodulatory regulatory T (Treg) cells. Despite its central role in Th lineage determination, molecular mechanisms of how IRF4 steers these diverse transcriptional programs in Th17 and Treg cells are far from being definitive. In the present study, we integrated data derived from affinity-purification and full mass spectrometry-based proteome analysis with chromatin immune precipitation sequencing (ChIP-Seq) to unveil IRF4-driven lineage determination in Treg and Th17 cells. This allowed the characterization of subtype-specific molecular programs and identification of novel, yet uncharted IRF4 interactors in the Th17/Treg context including master transcription factors FOXP3 and RORG as well as other players such as AHR, TBX21, IRF8, BACH2, SATB1, GTF2I, TCF7, FLI1 and RFX1. Moreover, our data reveal that most of these transcription factors are recruited to IRF composite elements for the regulation of cell type-specific transcriptional programs. To our knowledge, this is the first study investigating the molecular mechanisms of IRF4-steered gene expression and IRF4 interactions in an unbiased approach, directly comparing fully differentiated Th17 and Treg cells.


## Dependencies and preliminary setup
To Run the ChIP-script you'll need to download the mm10 Enhancer database here: 
http://www.enhanceratlas.org/
und Use UCSCs liftOver: https://genome.ucsc.edu/cgi-bin/hgLiftOver
The resulting file should be put here: databases/mm10_lifted.bed"
Similarly, you need to get the "Silencers_Mus_Musculus_blood" data from here:
http://health.tsinghua.edu.cn/silencerdb/
And lift it to mm10 and put it here: databases/silencers_Mus_musculus_Blood_lifted_mm10.bed

The Proteomic Input files are provided with the paper, you'll need to change the file names/locations possibly the sheet number.