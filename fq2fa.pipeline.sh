
#!/bin/bash
set -euo pipefail

##### ##### ##### Creating corrected Nanopore fastq reads ##### ##### #####

/u/ykaya/miniconda3/envs/path/bin/pod5 convert fast5 -o sampleX.pod5 -t 40 -r /ptmp/ykaya/fast5/EA/69-11/20210505_1603_X2_FAQ08451_969f9923/fast5

/u/ykaya/dorado-0.8.1-linux-x64/bin/dorado basecaller hac,5mCG_5hmCG sampleX.pod5 --models-directory /u/ykaya/dorado-0.8.1-linux-x64/bin/ --device cuda:0 --emit-moves --reference /u/ykaya/TAIR10.fa  > sampleX.bam

/u/ykaya/miniconda3/envs/path/bin/samtools fastq sampleX.bam > sampleX.fastq

/u/ykaya/dorado-0.8.1-linux-x64/bin/dorado download --model herro-v1
/u/ykaya/dorado-0.8.1-linux-x64/bin/dorado correct -m herro-v1 sampleX.fastq > corrected_sampleX.fasta

##### ##### ########## ##### 1st Draft Assembly w/ Flye  ##### ##### ########## #####

#!/usr/bin/env bash

/u/ykaya/Flye/bin/flye --nano-corr corrected_sampleX.fasta --threads 128 -g 135m --out-dir sampleX_Asm

########	########## 2nd Assembly w/ MABS-flye (10 iteration) ##############	#########


#!/usr/bin/env bash

# ------------- Configuration -------------
MABS_FLYE="/netscratch/irg/grp_hancock/Mabs-2.28/mabs-flye.py"
# ^-- Adjust this to the actual path or command you'd like to run if it is not mabs-flye.py

HIFI_DIR="/netscratch/irg/grp_hancock/Raw_HiFi_Data/mpgc_Arabidopsis_EA/fastq/raw"
ONT_DIR="/netscratch/irg/grp_hancock/Mehmet/Arabidopsis_thaliana_GlobalSet/rawdata_nanopore/Other_Samples"
OUTPUT_DIR="/netscratch/irg/grp_hancock/Raw_HiFi_Data/mpgc_Arabidopsis_EA/fastq/Mabs_assembly_Flye"

mkdir -p "$OUTPUT_DIR"

SAMPLES=(
    "ET108.1"
    "ET131.2"
    "ET148.10"
    "ET3.1"
    "ET33.1"
    "ET49.2"
    "ET53.7"
)

for SAMPLE in "${SAMPLES[@]}"; do
    fastq_file="${HIFI_DIR}/${SAMPLE}.fastq.gz"
    ultra_long_input="${ONT_DIR}/${SAMPLE}.fastq"

    # Optional checks
    if [[ ! -f "$fastq_file" ]]; then
        echo "WARNING: HiFi FASTQ not found: $fastq_file"
        continue
    fi
    if [[ ! -f "$ultra_long_input" ]]; then
        echo "WARNING: ONT FASTQ not found: $ultra_long_input"
        continue
    fi

    # Create a dedicated directory for each sample
    sample_dir="${OUTPUT_DIR}/${SAMPLE}"
    mkdir -p "$sample_dir"

    # The output prefix within the sample directory
    out_prefix="${sample_dir}/${SAMPLE}.asm"

    # Create the run script inside the sample directory
    run_script="${sample_dir}/run_mabs_flye.sh"
    cat <<EOF > "$run_script"
#!/usr/bin/env bash

echo "Running mabs-flye for sample: $SAMPLE"

"${MABS_FLYE}" \\
    --nanopore_reads "${ultra_long_input}" \\
    --download_busco_dataset eudicots_odb10.2020-09-10.tar.gz \\
    --genome_size 135m \\
    --output_folder "${out_prefix}" \\
    --threads 40

echo "Done with ${SAMPLE}"
EOF

    # Make the run script executable
    chmod +x "$run_script"

    # (Optional) Immediately run the script from within the sample directory.
    echo "Running script for sample: ${SAMPLE}"
    (cd "$sample_dir" && ./run_mabs_flye.sh)
done

echo "All samples processed."


############# Merge Assemblies: Assembly with the longest contiguity (N50) selected as query, the assembly with second longest N50 was used as reference to join query assembly

merge_wrapper.py hybrid_assembly.fasta self_assembly.fasta


########################## Polishing, Purging, and Scaffolding ######################

#!/bin/bash
set -euo pipefail

# -------------------------
# Define directories & files
# -------------------------
# Base directory for assemblies (each sample has its own subdirectory)
ASSEMBLIES_DIR="/path/sampleX_Asm"

# Directory with paired-end NGS reads (must have *_R1.fastq.gz and *_R2.fastq.gz files)
NGS_DATA_DIR="/netscratch/irg/grp_hancock/African_Genomes_Project/path/"

# Col-CEN reference genome for scaffolding
COLCEN_REF="/netscratch/irg/grp_hancock/Raw_HiFi_Data/mpgc_Arabidopsis_EA/path/ColCEN.fasta"

# -------------------------
# Define tool paths
# -------------------------
MINIMAP2="/opt/share/software/bin/minimap2"
RACON="/opt/share/software/bin/racon"
SAMTOOLS="/opt/share/software/packages/miniconda3-4.12.0/bin/samtools"
# Pilon – adjust if you need to run via java; here we assume it's directly executable:
PILON="/opt/share/software/bin/pilon"
# Tools for purge_dups
PBCSTAT="/opt/share/software/bin/pbcstat"
CALCUTS="/opt/share/software/bin/calcuts"
SPLIT_FA="/opt/share/software/bin/split_fa"
GET_SEQS="/opt/share/software/bin/get_seqs"
PURGE_DUPS="/opt/share/software/bin/purge_dups"
# Scaffolding tools
RAGTAG="/opt/share/software/bin/ragtag.py"
NUCMER="/opt/share/software/bin/nucmer"

# -------------------------
# Process each sample
# -------------------------
for sampleDir in "$ASSEMBLIES_DIR"/*; do
    if [ -d "$sampleDir" ]; then
        sampleId=$(basename "$sampleDir")
        echo "Processing sample: $sampleId"

        # Define the directory holding the final assembly and output files
        ASM_DIR="${sampleDir}/${sampleId}.asm/The_best_assembly"
        rawAssembly="${ASM_DIR}/assembly.fasta"

        # Corrected reads (if available) assumed to be in sample folder
        correctedReads="${sampleDir}/corrected_${sampleId}.fasta"

        # Create subdirectories inside ASM_DIR for intermediate outputs
        RACON_DIR="${ASM_DIR}/racon"
        PILON_DIR="${ASM_DIR}/pilon"
        PURGE_DIR="${ASM_DIR}/purge_dups"
        mkdir -p "$RACON_DIR" "$PILON_DIR" "$PURGE_DIR"

        # -------------------------
        # Racon polishing (using long reads)
        # -------------------------
        if [ ! -f "$correctedReads" ]; then
            echo "Corrected reads not found for sample $sampleId. Skipping Racon polishing..."
            raconPolished="$rawAssembly"
        else
            echo "Running Racon polishing for $sampleId..."
            # First round
            $MINIMAP2 -ax map-ont -t 20 "$rawAssembly" "$correctedReads" > "${RACON_DIR}/${sampleId}.1.sam"
            $RACON -m 8 -x -6 -g -8 -w 500 -t 20 "$correctedReads" "${RACON_DIR}/${sampleId}.1.sam" "$rawAssembly" > "${RACON_DIR}/${sampleId}.assembly_polished1.fasta"
            # Second round
            $MINIMAP2 -ax map-ont -t 20 "${RACON_DIR}/${sampleId}.assembly_polished1.fasta" "$correctedReads" > "${RACON_DIR}/${sampleId}.2.sam"
            $RACON -m 8 -x -6 -g -8 -w 500 -t 20 "$correctedReads" "${RACON_DIR}/${sampleId}.2.sam" "${RACON_DIR}/${sampleId}.assembly_polished1.fasta" > "${RACON_DIR}/${sampleId}.assembly_polished2.fasta"
            echo "Racon polishing completed for $sampleId."
            raconPolished="${RACON_DIR}/${sampleId}.assembly_polished2.fasta"
        fi

        # -------------------------
        # Pilon polishing (using NGS data)
        # -------------------------
        if [ -f "${NGS_DATA_DIR}/${sampleId}_R1.fastq.gz" ]; then
            ngsR1="${NGS_DATA_DIR}/${sampleId}_R1.fastq.gz"
            ngsR2="${NGS_DATA_DIR}/${sampleId}_R2.fastq.gz"
        elif [ -f "${NGS_DATA_DIR}/${sampleId}_R1.fastq" ]; then
            ngsR1="${NGS_DATA_DIR}/${sampleId}_R1.fastq"
            ngsR2="${NGS_DATA_DIR}/${sampleId}_R2.fastq"
        else
            echo "NGS data not found for sample $sampleId. Skipping Pilon polishing..."
            pilonPolished="$raconPolished"
        fi

        if [ -n "${ngsR1:-}" ]; then
            echo "Running Pilon polishing for $sampleId..."
            inputAssembly="$raconPolished"
            # Run three rounds of Pilon polishing
            for round in {1..3}; do
                samFile="${PILON_DIR}/${sampleId}.NGS.${round}.sam"
                bamFile="${PILON_DIR}/${sampleId}.NGS.${round}.sorted.bam"
                echo "Pilon round $round for sample $sampleId..."
                $MINIMAP2 -ax sr -t 50 "$inputAssembly" "$ngsR1" "$ngsR2" > "$samFile"
                $SAMTOOLS view --threads 50 -Sb "$samFile" | $SAMTOOLS sort --threads 50 -o "$bamFile"
                $SAMTOOLS index -@ 50 "$bamFile"
                # Run Pilon polishing. Adjust memory/parameters if needed.
                $PILON --genome "$inputAssembly" --frags "$bamFile" --outdir "$PILON_DIR" --output "${sampleId}.pilon${round}"
                inputAssembly="${PILON_DIR}/${sampleId}.pilon${round}.fasta"
                echo "Pilon round $round completed for $sampleId."
            done
            pilonPolished="$inputAssembly"
            echo "Pilon polishing completed for $sampleId."
        fi

        # -------------------------
        # Purge duplicates
        # -------------------------
        echo "Mapping reads to assembly for purging for sample $sampleId..."
        $MINIMAP2 -x map-ont -t 20 "$pilonPolished" "$correctedReads" -o "${PURGE_DIR}/${sampleId}.nano.paf"

        echo "Calculating depth statistics for $sampleId..."
        $PBCSTAT -O "$PURGE_DIR" "${PURGE_DIR}/${sampleId}.nano.paf"

        echo "Calculating cutoffs for $sampleId..."
        $CALCUTS "${PURGE_DIR}/PB.stat" > "${PURGE_DIR}/cutoffs.auto"

        echo "Splitting assembly for $sampleId..."
        $SPLIT_FA "$pilonPolished" > "${PURGE_DIR}/${sampleId}.split"

        echo "Self-mapping split assembly for $sampleId..."
        $MINIMAP2 -xasm5 -DP -t 20 "${PURGE_DIR}/${sampleId}.split" "${PURGE_DIR}/${sampleId}.split" -o "${PURGE_DIR}/${sampleId}.split.self.paf"

        echo "Identifying duplicates for $sampleId..."
        $PURGE_DUPS -2 -T "${PURGE_DIR}/cutoffs.auto" -c "${PURGE_DIR}/PB.base.cov" "${PURGE_DIR}/${sampleId}.split.self.paf" > "${PURGE_DIR}/dups.bed"

        echo "Extracting non-duplicate sequences for $sampleId..."
        $GET_SEQS -e -p "${ASM_DIR}/${sampleId}.purged" "${PURGE_DIR}/dups.bed" "$pilonPolished"
        purgedAssembly="${ASM_DIR}/${sampleId}.purged.fa"
        echo "Purging completed for $sampleId. Final purged file: $purgedAssembly"

        # -------------------------
        # Scaffolding with Ragtag
        # -------------------------
        echo "Running Ragtag scaffolding for sample $sampleId..."
        $RAGTAG scaffold "$COLCEN_REF" "$purgedAssembly" --aligner $NUCMER -o "${ASM_DIR}/${sampleId}_ragtag"
        echo "Scaffolding completed for $sampleId."

        echo "All steps completed for sample $sampleId."
        echo "---------------------------------------------------------"
    fi
done

echo "All samples processed."
~