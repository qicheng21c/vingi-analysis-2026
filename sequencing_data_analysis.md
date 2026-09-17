# Data analysis of junction-enrichment sequencing

## Annotate genomic regions

Genomic regions were summarized from raw data: `ncbiRefSeq.txt.gz` obtained 
from the UCSC Genome Browser using script `build_gene_annotation_features.py`.

## Process sequencing data

Before running the pipeline, create the required working directories, place 
the Python scripts in the scripts directory, and prepare the sequencing data 
and a file_list containing the input file names.

```bash
# Set the working directory.
export work_dir="WORK_DIR"
export scripts="${work_dir}/scripts"
export barcode_dir="${work_dir}/barcode"
export bwa_search_index="${work_dir}/bwa_index"
export demux_data="${work_dir}/demux_data"
export umi_deduped_data="${work_dir}/umi_deduped_data"
export bwa_result="${work_dir}/bwa_result"
export genome_mapping_results="${work_dir}/genome_mapping_results"
export genome_insertion_pos="${work_dir}/genome_insertion_pos"
export rebuild_utr3="${work_dir}/rebuild_utr3"
export utr3_remap="${work_dir}/utr3_remap"
export insertion_report="${work_dir}/insertion_report"

mkdir -p $work_dir $scripts $barcode_dir $bwa_search_index \
    $demux_data $umi_deduped_data $bwa_result $genome_mapping_results \
    $genome_insertion_pos $rebuild_utr3 $utr3_remap $insertion_report

# Specify the input data and BWA-MEM index used for the analysis.
export file_list="DATA_FILE_LIST"
export index_id="pe150_utr3"
export utr_realign_index_id="utr3_tail70"
```

1. **Deduplicate based on UMIs, retaining the reads with the highest sequencing quality.**

```bash
awk -F '\t' 'NR>1{print $1}' "${file_list}" | while read data_series; do
    python "${scripts}/umi_deduplication_of_demultiplexed_data.py" \
        "${demux_data}/${data_series}_1.fq.gz" \
        "${demux_data}/${data_series}_2.fq.gz" \
        "${umi_deduped_data}/${data_series}_1.fq.gz" \
        "${umi_deduped_data}/${data_series}_2.fq.gz"
done
```

2. **BWA-MEM Search against index of the fixed sequence.**

```bash
awk -F '\t' 'NR>1{print $1}' "${file_list}" | while read data_series; do
    bwa mem -t 12 "${bwa_search_index}/${index_id}.fa.index" \
        "${umi_deduped_data}/${data_series}_1.fq.gz" \
        "${umi_deduped_data}/${data_series}_2.fq.gz" \
        | tee "${bwa_result}/${data_series}.bam" \
        | samtools sort -@ 12 | samtools view -@ 12 -f 0x9 -F 0xF04 -h > "${bwa_result}/${data_series}.sam"
done
```

3. **Filter qualified reads (39M, S≥50).**

```bash
awk -F '\t' 'NR>1{print $1}' "${file_list}" | while read data_series; do
    python "${scripts}/filter_utr_bwa_result.py" \
        -f "${bwa_search_index}/${index_id}.fa" \
        -i "${bwa_result}/${data_series}.sam" \
        -o "${bwa_result}/${data_series}.filtered.tsv"
done
```

4. **Lentivirus sequencing data decontamination.**

```bash
export contamination_seq_id="lenti_contamination"
export contamination_removal="${batch_dir}/contamination_removal"

bwa index -a is lenti_contamination.fa -p lenti_contamination.fa.index

mkdir -p $contamination_removal ${contamination_removal}/cleaned_data
for file in ${bwa_result}/*.fa; do
    base_name=$(basename $file .fa)
    bwa mem -t 12 "${bwa_search_index}/${contamination_seq_id}.fa.index" \
    $file \
    | tee ${contamination_removal}/${base_name}.bam \
    | samtools sort -@ 12 | samtools view -@ 12 -F 4 -h > ${contamination_removal}/${base_name}.sam
done

for file in ${contamination_removal}/*.sam; do
    base_name=$(basename $file .sam)
    contaminated_list=${contamination_removal}/${base_name}.id
    raw_tbl=${bwa_result}/${base_name}.filtered.tsv
    filtered_tbl=${contamination_removal}/cleaned_data/${base_name}.tsv

    awk '$0!~/^@/{print $1}' $file > $contaminated_list
    grep -vf $contaminated_list $raw_tbl > $filtered_tbl
done

for file in ${contamination_removal}/cleaned_data/*.tsv; do
    base=$(basename $file .tsv)
    outfa="${contamination_removal}/cleaned_data/${base}.fa"
    awk -F '\t' '{printf ">%s %s\n%s\n", $1, $2, $4}' $file > $outfa
done
```

5. **Characterize the precise insertion boundary of RNA donor and junction scars.**

```bash
# 5.1. Genome alignment of soft-clipped sequences.
for file in ${bwa_result}/*.filtered.tsv; do
    base_name=$(basename $file .filtered.tsv)
    awk -F '\t' '{printf ">%s\n%s\n", $1, $4}' $file > ${bwa_result}/${base_name}.fa
done

# For Vingi data
parallel -j 10 --bar \
    "bwa mem -t 12 hg38.fa.index {} > ${genome_mapping_results}/{/.}.bam" \
    ::: ${bwa_result}/*.fa

# For Lentivirus data
parallel -j 10 --bar \
    "bwa mem -t 12 hg38.fa.index {} > ${genome_mapping_results}/{/.}.bam" \
    ::: ${contamination_removal}/cleaned_data/*.fa

# 5.2. Process genome alignment results to identify junction sequences.
for file in ${genome_mapping_results}/*.bam; do
    base=$(basename $file .bam)
    outsam=${genome_mapping_results}/${base}.sam
    samtools sort -@ 12 $file | samtools view -@ 12 -F 0x904 -h > $outsam
done

for file in "${genome_mapping_results}"/*.sam; do
    python "${scripts}/process_genome_mapping_result.py" "$file" "${file%.sam}_junc3.tsv"
done

for junc_data in ${genome_mapping_results}/*_junc3.tsv; do
    base_name=$(basename $junc_data "_junc3.tsv")
    utr_bwa_data=${bwa_result}/${base_name}.filtered.tsv
    outfile=${rebuild_utr3}/${base_name}.fa
    awk -F '\t' '
    NR==FNR{
        seq1[$1] = $3
        next
    }
    {
        if ($1 in seq1)
            seq2[$1] = seq1[$1] tolower($5)
    }
    END{
        for (i in seq2)
            printf ">%s\n%s\n", i, seq2[i] > outf
    }
    ' outf="$outfile" "$utr_bwa_data" "$junc_data"
done

# 5.3. Realign junctions to the donor sequence and calculate the per-position mapping depth.
for file in ${rebuild_utr3}/*.fa; do
    base_name=$(basename $file ".fa")
    bwa mem -t 12 ${bwa_search_index}/${utr_realign_index_id}.fa.index $file \
    | samtools sort | samtools view -h > ${utr3_remap}/${base_name}.sam
    python ${scripts}/process_utr_remap_result.py ${utr3_remap}/${base_name}.sam ${utr3_remap}/${base_name}.tsv
    python ${scripts}/utr3_map_depth_count.py ${utr3_remap}/${base_name}.tsv ${utr3_remap}/depth_${base_name}.tsv
done
```

6. **Get genomic insertion sites.**

```bash
for file in $genome_mapping_results/*_junc3.tsv; do
    base=$(basename $file _junc3.tsv)
    outf=${genome_insertion_pos}/${base}_insertion_pos.tsv
    python ${scripts}/get_insertion_pos.py $file $outf
done
```

7. **Map insertion sites onto annotated genomic regions.**

```bash
for file in ${genome_insertion_pos}/*_insertion_pos.tsv; do
    base=$(basename $file _insertion_pos.tsv)
    python ~/vingi/insertion_mapping/report_insertion_site.py $file > ${insertion_report}/${base}_insertion_report.tsv
done
```