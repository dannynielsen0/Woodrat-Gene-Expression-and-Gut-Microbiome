
### Reading the htseq count files into R was having an issue related to columns in the counts file and not matching, to fix this, just duplicate the first column (gene name). This code should work recursively within a directory over each tissue type. Substitute lepSecond_counts_files with bryFirst_counts_files to modify those

```
find lepSecond_counts_files -type f -name '*.tsv' -exec zsh -c '
    for file do
        awk -F"\t" '\''BEGIN {OFS = FS} $2 == "" { $2 = $1 } { print $0 }'\'' "$file" > "$file.tmp" && mv "$file.tmp" "$file"
    done
' zsh {} +
```

### The annotations (gff files) for both N. lepida and N. bryanti were generated from using Liftoff (Shumate and Salzberg 2021) with the following commands

 ```
    #For bryanti
        liftoff -g Neotoma_macrotis_annotations.gff bry_tig.fasta ../Mac_fus_genomes/Nmacrotis_1.1.fasta -u "unmapped_mac_bry.txt" -f other_features.txt -copies -sc 0.90 -flank 0.1 -o mac_to_bry.gff

    #For lepida, use the -db flag to call the liftoff database already made above
        liftoff -db ../mac_bry/Neotoma_macrotis_annotations.gff_db lep.fasta ../Mac_fus_genomes/Nmacrotis_1.1.fasta -u "unmapped_mac_lep.txt" -f ../mac_bry/other_features.txt -copies -sc 0.90 -flank 0.1 -o          mac_to_lep.gff
 ```

### The RNAseq data were aligned with hisat using the following SLURM script

```
        #!/bin/bash -l
    #SBATCH --nodes=1  --ntasks-per-node=1 --cpus-per-task=8 --mem-per-cpu=3500M
    #SBATCH --time=3-00:00:00
    #SBATCH --job-name rna_align_$SLURM_ARRAY_TASK_ID
    #SBATCH --output=/data/gpfs/assoc/matocqlab/Danny/bryFirst_lepLast_final_alignments/message_logs/%A.%a.out              # The output file name: <job_name>.<job_id>.out
    #SBATCH --error=/data/gpfs/assoc/matocqlab/Danny/bryFirst_lepLast_final_alignments/message_logs/%A.%a.err               # The error file name: <job_name>.<job_id>.err
    #SBATCH --account=cpu-s2-matocqlab-0            # The account to charge
    #SBATCH --partition=cpu-s2-core-0               # The parition
     
    
    
    module load unr-rc
    
    cd /data/gpfs/assoc/matocqlab/Danny/bryFirst_lepLast_final_alignments
    
    echo "All jobs in this array have:"
    echo "- SLURM_ARRAY_JOB_ID=${SLURM_ARRAY_JOB_ID}"
    echo "- SLURM_ARRAY_TASK_COUNT=${SLURM_ARRAY_TASK_COUNT}"
    echo "- SLURM_ARRAY_TASK_MIN=${SLURM_ARRAY_TASK_MIN}"
    echo "- SLURM_ARRAY_TASK_MAX=${SLURM_ARRAY_TASK_MAX}"
     
    echo "This job in the array has:"
    echo "- SLURM_JOB_ID=${SLURM_JOB_ID}"
    echo "- SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}"
    
    # grab our filename from a directory listing
    NAMES=($(cat /data/gpfs/assoc/matocqlab/Neotoma_transcriptomics/sample_names.txt))
    FILENAME=${NAMES[$SLURM_ARRAY_TASK_ID]}
    echo "My sample is ${FILENAME}"
    
    # make new directory, change into it, and run
    mkdir ${FILENAME}_align
    cd ${FILENAME}_align
    hisat2 -q --phred33 --no-temp-splicesite --no-mixed --no-discordant --max-intronlen 150000 \
    --rna-strandness RF --no-unal -p 8 \
    --un-conc /data/gpfs/assoc/matocqlab/Danny/bryFirst_lepLast_final_alignments/${FILENAME}_align/${FILENAME} \
    -x /data/gpfs/assoc/matocqlab/bry_DB/bry_DB \
    -1 /data/gpfs/assoc/matocqlab/Danny/FG_Cecum_RNA/Trimmed/${FILENAME}_R1_001_val_1.fq.gz \
    -2 /data/gpfs/assoc/matocqlab/Danny/FG_Cecum_RNA/Trimmed/${FILENAME}_R2_001_val_2.fq.gz \
    -S ${FILENAME}_align.sam
    
    samtools view -bh -@ 8 -f 3 -F 256 -q 40 ${FILENAME}_align.sam | \
    samtools sort -@ 8 > ${FILENAME}_align_sorted.bam
    samtools index -b ${FILENAME}_align_sorted.bam
    rm ${FILENAME}_align.sam
```

### The alignments were then used to count reads using HTSeq

```
#!/bin/bash
#SBATCH --nodes=1  --ntasks-per-node=1 --cpus-per-task=32 --constraint=intelv4 --mem-per-cpu=6500M
#SBATCH --time=1-00:00:00
#SBATCH --job-name htseq_$SLURM_ARRAY_TASK_ID
#SBATCH --output=bryFirst_htseq_counts/message_logs/%A.%a.out              # The output file name: <job_name>.<job_id>.out
#SBATCH --error=bryFirst_htseq_counts/message_logs/%A.%a.err               # The error file name: <job_name>.<job_id>.err
#SBATCH --account=cpu-s2-matocqlab-0            # The account to charge
#SBATCH --partition=cpu-s2-core-0               # The parition
 


module load unr-rc

cd  /data/gpfs/assoc/matocqlab/Danny

echo "All jobs in this array have:"
echo "- SLURM_ARRAY_JOB_ID=${SLURM_ARRAY_JOB_ID}"
echo "- SLURM_ARRAY_TASK_COUNT=${SLURM_ARRAY_TASK_COUNT}"
echo "- SLURM_ARRAY_TASK_MIN=${SLURM_ARRAY_TASK_MIN}"
echo "- SLURM_ARRAY_TASK_MAX=${SLURM_ARRAY_TASK_MAX}"
 
echo "This job in the array has:"
echo "- SLURM_JOB_ID=${SLURM_JOB_ID}"
echo "- SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}"

# grab our filename from a directory listing
NAMES=($(cat /data/gpfs/assoc/matocqlab/Danny/redo.txt))
FILENAME=${NAMES[$SLURM_ARRAY_TASK_ID]}
echo "My sample is ${FILENAME}"

# make new directory, change into it, and run
mkdir /data/gpfs/assoc/matocqlab/Danny/bryFirst_htseq_counts/${FILENAME}_htseq
cd /data/gpfs/assoc/matocqlab/Danny/bryFirst_htseq_counts/${FILENAME}_htseq

htseq-count --format bam \
--stranded=reverse \
--type=exon \
--idattr=gene_id \
--add-chromosome-info \
--additional-attr=Name \
--order=pos \
--mode=union \
--nonunique=random \
--secondary-alignments=ignore \
--samout-format=bam \
--samout=/data/gpfs/assoc/matocqlab/Danny/bryFirst_htseq_counts/${FILENAME}_htseq/${FILENAME}.bam \
--counts_output=/data/gpfs/assoc/matocqlab/Danny/bryFirst_htseq_counts/${FILENAME}_htseq/${FILENAME}.counts.tsv \
-n 32 \
/data/gpfs/assoc/matocqlab/Danny/bryFirst_lepLast_final_alignments/${FILENAME}_align/${FILENAME}_align_sorted.bam \
/data/gpfs/assoc/matocqlab/Danny/final_detox_gffs/bry_final.gff
```

### Counts were then imported into R

    ####We created a DESeq object using the DESeqDataSetFromHTSeqCount function

```
        #This is the command, and we adjust the inputs for each tissue (i.e., foregut, caecum, etc...)
        
            bryFirst_liver_dds <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                  directory = directory,
                                  design= ~ 0 + group)

        #Then, we performed variance stabalization and scaled and centered the data

            vsd <- vst(bry_first_liver_dds, blind=FALSE)
            CEN = scale(assay(vsd), center = T, scale = T)
            pca <- prcomp(t(CEN))
            loadings <- as.data.frame(pca$rotation)

```

### With this, plot PCA
```
            biplot <- factoextra::fviz_pca_ind(pca, repel = TRUE, axes = c(1,2),
                          #select.var = list(contrib = 25), #draw top 25 arrows
                          #select.var = list(name = c("Sult2b1-1", "APOA4_0001", "FADS3_0001","Sult2a3-7")),  #alternative to draw specific substitution loadings
                          labelsize=20,
                          addEllipses = TRUE,
                          habillage = bry_first_caecum_dds$group,
                          col.ind = bry_first_caecum_dds$group,
                          ellipse.level=0.95,
                          palette = c("maroon","pink", "lightgreen", "darkgreen"),
                          geom=c("point"), pointsize = 15,   #change to geom=c("point","text") for sample ID
                          ind.fill = bry_first_caecum_dds$group,
                          invisible = c( "quali"), #remove enlarged symbol for group mean
                          title = "Caecum") +
                          scale_shape_manual(values=c(18,15,17,16)) + 
                          guides(colour = guide_legend(override.aes = list(size = 10))) +
                          theme(plot.title=element_text(hjust=0.5),
                                text = element_text(size = 50),
                                axis.title = element_text(size = 40),
                                axis.text = element_text(size = 30),
                                legend.text = element_text(size = 30),  # Increase legend text size
                                legend.title = element_blank(),
                                legend.position = "none",  # Move legend to bottom
                                legend.direction = "vertical")

```

![PCA of Gene Expression in the Caecum](https://github.com/dannynielsen0/Lab-diet-trial-16S-analysis/blob/main/figures/caecum_gene_expression_pca.jpg?raw=true)



