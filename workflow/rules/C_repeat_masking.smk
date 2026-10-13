# ===============================================================================
# GABi - Repeat Masking Rules
# by Diego De Panis, 2026
# note: AI tools may have been used to improve, clean and/or comment this version of the code
# ===============================================================================

# -------------------------------------------------------------------------------
#  HELPER FUNCTIONS
# -------------------------------------------------------------------------------

#


# -------------------------------------------------------------------------------
#  REPEAT MODELER + MASKER
# ===============================================================================

rule C01_repeat_masking:
    input:
        assembly=lambda wildcards: get_assembly_path(wildcards.species, wildcards.assembly),
        compleasm=os.path.join(config["OUT_FOLDER"], "GABi_results", "{species}", "{assembly}",
                              "assembly_stats", "compleasm", "summary.txt")
    output:
        # Main masking outputs (direct to masking/)
        masked_soft=os.path.join(config["OUT_FOLDER"], "GABi_results", "{species}", "{assembly}",
                                "masking", "{assembly}.masked.fa"),
        repeat_out=os.path.join(config["OUT_FOLDER"], "GABi_results", "{species}", "{assembly}",
                               "masking", "{assembly}.fa.out"),
        repeat_gff=os.path.join(config["OUT_FOLDER"], "GABi_results", "{species}", "{assembly}",
                               "masking", "{assembly}.fa.out.gff"),
        repeat_align=os.path.join(config["OUT_FOLDER"], "GABi_results", "{species}", "{assembly}",
                                 "masking", "{assembly}.fa.align"),
        # Remove .cat from required outputs since it's not always generated
        repeat_tbl=os.path.join(config["OUT_FOLDER"], "GABi_results", "{species}", "{assembly}",
                               "masking", "{assembly}.fa.tbl"),
        # Repeat landscape analysis outputs
        align_divsum=os.path.join(config["OUT_FOLDER"], "GABi_results", "{species}", "{assembly}",
                                 "masking", "{assembly}.fa.align.divsum"),
        align_with_div=os.path.join(config["OUT_FOLDER"], "GABi_results", "{species}", "{assembly}",
                                   "masking", "{assembly}.fa.align_with_div"),
        landscape_html=os.path.join(config["OUT_FOLDER"], "GABi_results", "{species}", "{assembly}",
                                   "masking", "{assembly}.fa.align.divsum.html"),
        # RepeatModeler outputs
        families_fa=os.path.join(config["OUT_FOLDER"], "GABi_results", "{species}", "{assembly}",
                                "masking", "modeler", "{assembly}-families.fa"),
        families_stk=os.path.join(config["OUT_FOLDER"], "GABi_results", "{species}", "{assembly}",
                                 "masking", "modeler", "{assembly}-families.stk"),
        # Library outputs
        combined_lib=os.path.join(config["OUT_FOLDER"], "GABi_results", "{species}", "{assembly}",
                                 "masking", "libraries", "{assembly}_combined.fa"),
        species_lib=os.path.join(config["OUT_FOLDER"], "GABi_results", "{species}", "{assembly}",
                                "masking", "libraries", "{assembly}-rm.fa")
    container: CONTAINERS["tetools"]
    params:
        output_dir=os.path.join(config["OUT_FOLDER"], "GABi_results", "{species}", "{assembly}", "masking"),
        # RepeatMasker 4.2.x with the RMBlast engine treats -pa N as "N parallel jobs,
        # 4 cores each", so -pa N really consumes N*4 cores. Passing -pa {threads} (the
        # old behaviour) oversubscribed the allocation 4x. Keep at least 1.
        pa=lambda wildcards, threads: max(1, threads // 4)
    benchmark:
        os.path.join(config["OUT_FOLDER"], "GABi_results", "{species}", "{assembly}", "logs",
                    "C01_repeat_masking_benchmark.txt")
    log:
        os.path.join(config["OUT_FOLDER"], "GABi_results", "{species}", "{assembly}", "logs",
                    "C01_repeat_masking.log")
    threads: cpu_func("repeat_masking")
    resources:
        mem_mb=mem_func("repeat_masking"),
        runtime=time_func("repeat_masking")
    shell:
        """
        set -euo pipefail
        export LC_ALL=C
        mkdir -p $(dirname {log})

        exec > {log} 2>&1

        echo "[GABi] Starting repeat masking for {wildcards.species} {wildcards.assembly}"


        # FAST TEMP DIRECTORY SETUP
        # -------------------------------------------------------------------
        # try GABi_FAST_TMP first, then fall back to regular temp

        REQUIRED_GB=150
        PENALTY_PER_JOB=50

        get_fast_workdir() {{
            local REQUIRED_GB="${{1:-150}}"

            echo "[GABi] Fast temp directory selection for RepeatMasker:" >&2
            echo "[GABi]   Required: ${{REQUIRED_GB}} GB" >&2

            # First try GABi_FAST_TMP (RAM-based or fast local SSD)
            if [ -n "${{GABi_FAST_TMP:-}}" ]; then
                echo "[GABi]   Checking FAST_TMP: $GABi_FAST_TMP" >&2

                if mkdir -p "$GABi_FAST_TMP" 2>/dev/null; then
                    AVAIL_KB=$(df -k "$GABi_FAST_TMP" 2>/dev/null | awk 'NR==2 {{print $4}}' || echo "0")
                    AVAIL_GB=$((AVAIL_KB / 1024 / 1024))

                    # Count existing GABi masking jobs
                    EXISTING_JOBS=$(find "$GABi_FAST_TMP" -maxdepth 1 -type d -name "GABi_mask_*" 2>/dev/null | wc -l || echo "0")
                    PENALTY=$((EXISTING_JOBS * PENALTY_PER_JOB))
                    ADJUSTED_GB=$((AVAIL_GB - PENALTY))

                    echo "[GABi]     Available: ${{AVAIL_GB}} GB (adjusted: ${{ADJUSTED_GB}} GB after $EXISTING_JOBS concurrent jobs)" >&2

                    if [ "${{ADJUSTED_GB:-0}}" -ge "$REQUIRED_GB" ]; then
                        echo "[GABi] ⚡ Using RAM-mounted FAST_TMP: $GABi_FAST_TMP" >&2
                        echo "$GABi_FAST_TMP"
                        return 0
                    fi
                    echo "[GABi]     ❌ Insufficient space in FAST_TMP" >&2
                fi
            fi

            # Fall back to regular gabi_get_workdir
            echo "[GABi]   Falling back to regular temp directory..." >&2
            gabi_get_workdir "$REQUIRED_GB"
        }}

        WORK_DIR="$(get_fast_workdir $REQUIRED_GB)"
        TEMP_WORK_DIR="$(mktemp -d "$WORK_DIR/GABi_mask_{wildcards.species}_{wildcards.assembly}_XXXXXX")"
        trap 'echo "[GABi] Cleaning up temp directory..."; rm -rf "$TEMP_WORK_DIR"' EXIT

        mkdir -p "$TEMP_WORK_DIR"/modeler "$TEMP_WORK_DIR"/libraries "$TEMP_WORK_DIR"/masker

        echo "[GABi] Working in temporary directory: $TEMP_WORK_DIR"

        # RepeatMasker/RepeatModeler 4.2.x keep a per-user cache (.RepeatMaskerCache)
        # and scratch dotfiles in $HOME. GABi may have multiple C01 jobs in flight at once,
        # all sharing one $HOME, so give each job a private one inside its temp dir.
        # The EXIT trap above removes it, so nothing is left behind between runs.
        export HOME="$TEMP_WORK_DIR"


        # EXTRACT LINEAGE FROM COMPLEASM
        # -------------------------------------------------------------------
        # Both greps below are allowed to match nothing, hence the `|| true`:
        # under `set -euo pipefail` a grep that finds no line exits 1 and takes
        # the whole rule down. With LINEAGE: eukaryota the summary holds exactly
        # one lineage, so `grep -v eukaryota` legitimately matches nothing and
        # used to kill the job here - silently, and before the two fallbacks
        # below (which exist for precisely that case) could ever run.
        LINEAGES=$(grep "^## lineage:" {input.compleasm} | sed 's/## lineage: //g' | sed 's/_odb.*//g' || true)
        if [ -z "$LINEAGES" ]; then
            echo "[GABi] ERROR: no '## lineage:' line in {input.compleasm}" >&2
            echo "[GABi]        The compleasm summary is empty or truncated, so the" >&2
            echo "[GABi]        repeat library cannot be chosen. Re-run A02 for this" >&2
            echo "[GABi]        assembly." >&2
            exit 1
        fi

        # Prefer non-eukaryota lineage if available, otherwise use first lineage
        LINEAGE=$(echo "$LINEAGES" | grep -v "^eukaryota$" | head -n1 || true)
        if [ -z "$LINEAGE" ]; then
            LINEAGE=$(echo "$LINEAGES" | head -n1)
        fi
        if [ -z "$LINEAGE" ]; then
            LINEAGE="eukaryota"  # fallback
        fi

        echo "[GABi] Using lineage: $LINEAGE"

        # Create output directories
        mkdir -p {params.output_dir}/modeler
        mkdir -p {params.output_dir}/libraries


        # 1. PREPARE ASSEMBLY FILE
        # -------------------------------------------------------------------
        cd $TEMP_WORK_DIR/modeler

        if [[ "{input.assembly}" == *.gz ]]; then
            echo "[GABi] Decompressing assembly..."
            gunzip -c {input.assembly} > {wildcards.assembly}.fa
        else
            echo "[GABi] Linking assembly..."
            ln -s {input.assembly} {wildcards.assembly}.fa
        fi


        # 2. REPEATMODELER
        # -------------------------------------------------------------------
        echo "[GABi] Building RepeatModeler database..."
        BuildDatabase -name {wildcards.assembly} {wildcards.assembly}.fa

        echo "[GABi] Running RepeatModeler..."
        # -LTRStruct stays on: the TETools 2.00 image still ships the whole LTR pipeline
        # RepeatModeler needs for it (LTR_retriever, LtrDetector, GenomeTools/LtrHarvest,
        # cd-hit and Ninja), and RepModelConfig.pm in the image already points at those
        # bundled copies, so no -LTRStruct path has to be configured here. If a future
        # image drops one of them, RepeatModeler fails at startup with the missing
        # tool's name -- check RepModelConfig.pm in the image before touching this line.
        # NB -threads (RepeatModeler 2.0.2+) is a real core count, unlike RepeatMasker's
        # -pa below.
        RepeatModeler -database {wildcards.assembly} -threads {threads} -LTRStruct


        # 3. PREPARE REPEAT LIBRARIES
        # -------------------------------------------------------------------
        cd $TEMP_WORK_DIR/libraries

        echo "[GABi] Extracting species-specific repeat library for: $LINEAGE"

        # Locate famdb.py and the FamDB partition directory. Dfam TETools moved both
        # between image releases (1.8x keeps them under /opt/RepeatMasker, 2.00 ships
        # FamDB separately), so probe the known layouts instead of hardcoding one path.
        FAMDB_BIN="$(command -v famdb.py 2>/dev/null || true)"
        if [ -z "$FAMDB_BIN" ]; then
            for f in /opt/RepeatMasker/famdb.py /opt/FamDB/famdb.py /opt/FamDB-Dfam-*/famdb.py; do
                if [ -f "$f" ]; then FAMDB_BIN="$f"; break; fi
            done
        fi
        if [ -z "$FAMDB_BIN" ]; then
            echo "[GABi] ERROR: famdb.py not found in the TETools container." >&2
            exit 1
        fi

        FAMDB_DIR=""
        for d in /opt/RepeatMasker/Libraries/famdb /opt/FamDB/Libraries/famdb \
                 /opt/FamDB/famdb /opt/FamDB-Dfam-*/Libraries/famdb; do
            if [ -d "$d" ] && ls "$d"/*.h5 >/dev/null 2>&1; then FAMDB_DIR="$d"; break; fi
        done
        if [ -z "$FAMDB_DIR" ]; then
            echo "[GABi] ERROR: no FamDB partition (*.h5) found in the TETools container." >&2
            exit 1
        fi

        echo "[GABi] famdb.py:  $FAMDB_BIN"
        echo "[GABi] FamDB dir: $FAMDB_DIR"
        echo "[GABi] FamDB partitions: $(cd "$FAMDB_DIR" && echo *.h5)"

        if [ -x "$FAMDB_BIN" ]; then
            FAMDB_CMD="$FAMDB_BIN"
        else
            FAMDB_CMD="python3 $FAMDB_BIN"
        fi

        $FAMDB_CMD -i "$FAMDB_DIR" families \
            --format fasta_name --include-class-in-name --ancestors --descendants \
            "$LINEAGE" > {wildcards.assembly}-rm.fa

        N_FAM=$(grep -c '^>' {wildcards.assembly}-rm.fa || true)
        echo "[GABi] FamDB families extracted for '$LINEAGE': ${{N_FAM:-0}}"
        if [ "${{N_FAM:-0}}" -eq 0 ]; then
            echo "[GABi] ERROR: FamDB returned 0 families for lineage '$LINEAGE'." >&2
            echo "[GABi]        Masking would silently fall back to the RepeatModeler" >&2
            echo "[GABi]        library alone, so this is treated as a failure." >&2
            echo "[GABi]        Check that the FamDB partition covering '$LINEAGE' is" >&2
            echo "[GABi]        present in $FAMDB_DIR (Dfam is shipped partitioned; the" >&2
            echo "[GABi]        root partition alone covers only a few clades)." >&2
            exit 1
        fi

        echo "[GABi] Combining repeat libraries..."
        cat {wildcards.assembly}-rm.fa $TEMP_WORK_DIR/modeler/{wildcards.assembly}-families.fa > {wildcards.assembly}_combined.fa


        # 4. REPEATMASKER
        # -------------------------------------------------------------------
        cd $TEMP_WORK_DIR/masker
        ln -s $TEMP_WORK_DIR/modeler/{wildcards.assembly}.fa {wildcards.assembly}.fa

        echo "[GABi] Running RepeatMasker..."
        RepeatMasker -pa {params.pa} -a -s -gccalc -xsmall -gff \
            -lib $TEMP_WORK_DIR/libraries/{wildcards.assembly}_combined.fa \
            {wildcards.assembly}.fa

        # Rename masked file
        mv {wildcards.assembly}.fa.masked {wildcards.assembly}.masked.fa


        # 5. GENERATE REPEAT LANDSCAPE ANALYSIS
        # -------------------------------------------------------------------
        echo "[GABi] Generating repeat landscape analysis..."

        # Calculate divergence from alignment
        calcDivergenceFromAlign.pl -s {wildcards.assembly}.fa.align.divsum \
            -a {wildcards.assembly}.fa.align_with_div {wildcards.assembly}.fa.align

        # Extract genome size from table file.
        # Same pipefail trap as the lineage greps above: a truncated .tbl made
        # this assignment fail and killed the job here with no message, after
        # RepeatModeler and RepeatMasker had already done all the work.
        SIZE=$(grep "total length:" {wildcards.assembly}.fa.tbl | awk '{{print $3}}' || true)
        if [ -z "$SIZE" ]; then
            echo "[GABi] ERROR: no 'total length:' line in {wildcards.assembly}.fa.tbl" >&2
            echo "[GABi]        RepeatMasker's table is incomplete, so the repeat" >&2
            echo "[GABi]        landscape cannot be scaled to the genome size." >&2
            exit 1
        fi
        echo "[GABi] Genome size detected: $SIZE bp"

        # Create repeat landscape plot
        createRepeatLandscape.pl -div {wildcards.assembly}.fa.align.divsum \
            -g $SIZE > {wildcards.assembly}.fa.align.divsum.html


        # 6. COPY AND ORGANIZE ALL OUTPUTS
        # -------------------------------------------------------------------
        echo "[GABi] Organizing output files..."

        # 6a. Copy RepeatMasker outputs

        # Copy main RepeatMasker outputs to masking/
        cp $TEMP_WORK_DIR/masker/{wildcards.assembly}.masked.fa {output.masked_soft}
        cp $TEMP_WORK_DIR/masker/{wildcards.assembly}.fa.out {output.repeat_out}
        cp $TEMP_WORK_DIR/masker/{wildcards.assembly}.fa.out.gff {output.repeat_gff}
        cp $TEMP_WORK_DIR/masker/{wildcards.assembly}.fa.align {output.repeat_align}
        cp $TEMP_WORK_DIR/masker/{wildcards.assembly}.fa.tbl {output.repeat_tbl}

        # Copy repeat landscape analysis outputs
        cp $TEMP_WORK_DIR/masker/{wildcards.assembly}.fa.align.divsum {output.align_divsum}
        cp $TEMP_WORK_DIR/masker/{wildcards.assembly}.fa.align_with_div {output.align_with_div}
        cp $TEMP_WORK_DIR/masker/{wildcards.assembly}.fa.align.divsum.html {output.landscape_html}

        # Note: We're not copying the .cat file as it's not used in the pipeline

        # 6b. Organize RepeatModeler outputs
        cd $TEMP_WORK_DIR/modeler

        # Tar the RM_* directory and remove it
        if ls -d RM_* >/dev/null 2>&1; then
            echo "[GABi] Archiving RepeatModeler work directory..."
            RM_DIR=$(ls -d RM_* | head -n1)
            tar -cf ${{RM_DIR}}.tar ${{RM_DIR}}
            rm -rf ${{RM_DIR}}  # Remove only the directory, not the tar
        fi

        # Remove the assembly file (we don't need it in the output)
        rm -f {wildcards.assembly}.fa

        # Now copy everything that's left in the modeler directory
        cp -r * {params.output_dir}/modeler/

        # 6c. Copy library files
        cp $TEMP_WORK_DIR/libraries/{wildcards.assembly}_combined.fa {output.combined_lib}
        cp $TEMP_WORK_DIR/libraries/{wildcards.assembly}-rm.fa {output.species_lib}

        echo "[GABi] RepeatModeler/RepeatMasker completed successfully!"
        """


# Rule 2: Post-processing (compression and final organization)
rule C02_masking_post:
    input:
        # Main files from repeat_masking
        masked_soft=os.path.join(config["OUT_FOLDER"], "GABi_results", "{species}", "{assembly}",
                                "masking", "{assembly}.masked.fa"),
        repeat_out=os.path.join(config["OUT_FOLDER"], "GABi_results", "{species}", "{assembly}",
                               "masking", "{assembly}.fa.out"),
        repeat_gff=os.path.join(config["OUT_FOLDER"], "GABi_results", "{species}", "{assembly}",
                               "masking", "{assembly}.fa.out.gff"),
        repeat_align=os.path.join(config["OUT_FOLDER"], "GABi_results", "{species}", "{assembly}",
                                 "masking", "{assembly}.fa.align"),
        # Large landscape file that should be compressed in post-processing
        align_with_div=os.path.join(config["OUT_FOLDER"], "GABi_results", "{species}", "{assembly}",
                                   "masking", "{assembly}.fa.align_with_div"),
    output:
        # Compressed final outputs
        masked_soft_gz=os.path.join(config["OUT_FOLDER"], "GABi_results", "{species}", "{assembly}",
                                   "masking", "{assembly}.masked.fa.gz"),
        repeat_out_gz=os.path.join(config["OUT_FOLDER"], "GABi_results", "{species}", "{assembly}",
                                  "masking", "{assembly}.fa.out.gz"),
        repeat_gff_gz=os.path.join(config["OUT_FOLDER"], "GABi_results", "{species}", "{assembly}",
                                  "masking", "{assembly}.fa.out.gff.gz"),
        bed_file_gz=os.path.join(config["OUT_FOLDER"], "GABi_results", "{species}", "{assembly}",
                                "masking", "{assembly}.mask_from_rm.bed.gz"),
        bed_file_tbi=os.path.join(config["OUT_FOLDER"], "GABi_results", "{species}", "{assembly}",
                                "masking", "{assembly}.mask_from_rm.bed.gz.tbi"),
        # Compressed landscape file (large file)
        align_with_div_gz=os.path.join(config["OUT_FOLDER"], "GABi_results", "{species}", "{assembly}",
                                      "masking", "{assembly}.fa.align_with_div.gz"),
        # Compressed library files
        families_fa_gz=os.path.join(config["OUT_FOLDER"], "GABi_results", "{species}", "{assembly}",
                                   "masking", "modeler", "{assembly}-families.fa.gz"),
        combined_lib_gz=os.path.join(config["OUT_FOLDER"], "GABi_results", "{species}", "{assembly}",
                                    "masking", "libraries", "{assembly}_combined.fa.gz"),
        species_lib_gz=os.path.join(config["OUT_FOLDER"], "GABi_results", "{species}", "{assembly}",
                                   "masking", "libraries", "{assembly}-rm.fa.gz")
    container: CONTAINERS["gabi_base"]
    params:
        masking_dir=os.path.join(config["OUT_FOLDER"], "GABi_results", "{species}", "{assembly}", "masking")
    benchmark:
        os.path.join(config["OUT_FOLDER"], "GABi_results", "{species}", "{assembly}", "logs",
                    "C02_masking_post_benchmark.txt")
    log:
        os.path.join(config["OUT_FOLDER"], "GABi_results", "{species}", "{assembly}", "logs",
                    "C02_masking_post.log")
    threads: cpu_func("masking_post")
    resources:
        mem_mb=mem_func("masking_post"),
        runtime=time_func("masking_post")
    shell:
        r'''
        set -euo pipefail
        export LC_ALL=C
        mkdir -p "$(dirname "{log}")"

        exec > "{log}" 2>&1

        echo "[GABi] Masking post-processing for {wildcards.species} {wildcards.assembly}"


        # TEMP DIRECTORY SETUP
        # -------------------------------------------------------------------
        WORK_DIR="$(gabi_get_workdir 50)"
        TEMP_DIR="$(mktemp -d "$WORK_DIR/GABi_maskpost_{wildcards.species}_{wildcards.assembly}_XXXXXX")"
        trap 'rm -rf "$TEMP_DIR"' EXIT
        cd "$TEMP_DIR"

        echo "[GABi] Working in: $TEMP_DIR"

        # SAFE COMPRESSION FUNCTION
        # Compresses a file safely without losing data on failure
        compress_safely() {{
            local infile="$1"
            local outfile="$2"
            local threads="{threads}"

            echo "[GABi]   Compressing: $(basename "$infile")"

            # Check input file exists and is readable
            if [[ ! -f "$infile" ]] || [[ ! -r "$infile" ]]; then
                echo "[GABi]   ERROR: Input file not found or not readable: $infile"
                return 1
            fi

            # Get original size for verification
            local orig_size=$(stat -c%s "$infile")
            echo "[GABi]     Original size: $orig_size bytes"

            # Compress to temporary file (atomic operation)
            local temp_gz="${{outfile}}.tmp.$$"
            if ! pigz -p "$threads" -c "$infile" > "$temp_gz" 2>&1; then
                echo "[GABi]   ERROR: Compression failed for $infile"
                rm -f "$temp_gz"
                return 1
            fi

            # Verify compressed file integrity
            echo "[GABi]     Verifying compressed file..."
            if ! pigz -t "$temp_gz" 2>&1; then
                echo "[GABi]   ERROR: Compressed file failed validation: $infile"
                rm -f "$temp_gz"
                return 1
            fi

            # Verify decompressed size matches (optional but recommended for critical data)
            local decomp_size=$(pigz -dc "$temp_gz" 2>&1 | wc -c)
            if [[ "$orig_size" -ne "$decomp_size" ]]; then
                echo "[GABi]   ERROR: Size mismatch! Original: $orig_size, Decompressed: $decomp_size"
                rm -f "$temp_gz"
                return 1
            fi

            # Move to final location (atomic on same filesystem)
            if ! mv "$temp_gz" "$outfile" 2>&1; then
                echo "[GABi]   ERROR: Failed to move compressed file to final location"
                rm -f "$temp_gz"
                return 1
            fi

            # Only NOW remove the original file (after all checks passed)
            if ! rm -f "$infile" 2>&1; then
                echo "[GABi]   WARNING: Failed to remove original file: $infile"
                echo "[GABi]   WARNING: Compressed file is safe at: $outfile"
                # Don't return error - compressed file is created successfully
            fi

            local comp_size=$(stat -c%s "$outfile")
            local ratio=$(awk "BEGIN {{printf \"%.1f\", 100.0 * $comp_size / $orig_size}}")
            echo "[GABi]   ✓ Success: $(basename "$infile") (${{ratio}}% of original)"

            return 0
        }}

        # 1. BUILD BED FILE
        echo "[GABi] Converting RepeatMasker .out to 3-col BED..."

        if ! command -v rmsk2bed >/dev/null 2>&1; then
          echo "[GABi] ERROR: rmsk2bed not found in PATH"
          exit 127
        fi

        set -o pipefail
        rmsk2bed < "{input.repeat_out}" \
          | cut -f1-3 \
          | sort -k1,1 -k2,2n \
          | bedtools merge -i - \
          | tee {wildcards.assembly}.3cols.bed \
          | bgzip -@ {threads} -c > "{output.bed_file_gz}"

        # Verify BED file was created successfully before indexing
        if [[ ! -f "{output.bed_file_gz}" ]] || [[ ! -s "{output.bed_file_gz}" ]]; then
            echo "[GABi] ERROR: BED file creation failed"
            exit 1
        fi

        tabix -p bed "{output.bed_file_gz}"

        if [[ ! -f "{output.bed_file_tbi}" ]]; then
            echo "[GABi] ERROR: BED indexing failed"
            exit 1
        fi

        # 2. COMPRESS ALL FILES SAFELY
        echo "[GABi] Compressing files safely (with verification)..."

        # Track failures
        FAILED=0

        # Compress each file - stop on first failure
        compress_safely "{input.masked_soft}" "{output.masked_soft_gz}" || FAILED=1
        [[ $FAILED -eq 0 ]] && compress_safely "{input.repeat_out}" "{output.repeat_out_gz}" || FAILED=1
        [[ $FAILED -eq 0 ]] && compress_safely "{input.repeat_gff}" "{output.repeat_gff_gz}" || FAILED=1
        [[ $FAILED -eq 0 ]] && compress_safely "{input.align_with_div}" "{output.align_with_div_gz}" || FAILED=1
        [[ $FAILED -eq 0 ]] && compress_safely "{params.masking_dir}/modeler/{wildcards.assembly}-families.fa" \
                                               "{output.families_fa_gz}" || FAILED=1
        [[ $FAILED -eq 0 ]] && compress_safely "{params.masking_dir}/libraries/{wildcards.assembly}_combined.fa" \
                                               "{output.combined_lib_gz}" || FAILED=1
        [[ $FAILED -eq 0 ]] && compress_safely "{params.masking_dir}/libraries/{wildcards.assembly}-rm.fa" \
                                               "{output.species_lib_gz}" || FAILED=1

        if [[ $FAILED -ne 0 ]]; then
            echo "[GABi] ERROR: One or more compression operations failed!"
            echo "[GABi] Original files have NOT been deleted for failed compressions."
            exit 1
        fi

        echo "[GABi] Done - all files compressed successfully."

        '''


# Rule 3: Plot repeat landscape as SVG
rule C03_plot_repeat_landscape:
    input:
        html=os.path.join(config["OUT_FOLDER"], "GABi_results", "{species}", "{assembly}",
                         "masking", "{assembly}.fa.align.divsum.html")
    output:
        svg=os.path.join(config["OUT_FOLDER"], "GABi_results", "{species}", "{assembly}",
                        "masking", "{assembly}.repeat_landscape.svg"),
        tsv=os.path.join(config["OUT_FOLDER"], "GABi_results", "{species}", "{assembly}",
                        "masking", "{assembly}.repeat_landscape.tsv")
    container: CONTAINERS["gabi_base"]  # Assumes R and ggplot2 are available
    params:
        html2tsv_script=str(scripts_dir / "RMhtml2TSV.py"),
        plot_script=str(scripts_dir / "plotRM.R")
    benchmark:
        os.path.join(config["OUT_FOLDER"], "GABi_results", "{species}", "{assembly}", "logs",
                    "C03_plot_repeat_landscape_benchmark.txt")
    log:
        os.path.join(config["OUT_FOLDER"], "GABi_results", "{species}", "{assembly}", "logs",
                    "C03_plot_repeat_landscape.log")
    threads: cpu_func("split_intervals")
    resources:
        mem_mb=mem_func("split_intervals"),
        runtime=time_func("split_intervals")
    shell:
        """
        set -euo pipefail
        mkdir -p $(dirname {log})

        exec > {log} 2>&1

        echo "[GABi] Converting HTML to TSV for {wildcards.species} {wildcards.assembly}"

        # Convert HTML to TSV
        python {params.html2tsv_script} -html {input.html} -out {output.tsv}

        # Create config file for R script (format: sample_id<TAB>tsv_path)
        TEMP_CFG=$(mktemp)
        trap 'rm -f "$TEMP_CFG"' EXIT

        echo -e "{wildcards.assembly}\t$(realpath {output.tsv})" > "$TEMP_CFG"

        echo "[GABi] Creating SVG plot"

        # Run R plotting script
        Rscript {params.plot_script} "$TEMP_CFG" {output.svg}

        echo "[GABi] Done: SVG created at {output.svg}"

        """

# Rule 4: for external bed file
rule C04_process_external_bed:
    """Process and validate a sibling BED file found next to the assembly"""
    input:
        bed=lambda w: get_sibling_bed_path(w.species, w.assembly) or []
    output:
        bed_gz=os.path.join(
            config["OUT_FOLDER"], "GABi_results", "{species}", "{assembly}",
            "masking", "{assembly}.mask_from_file.bed.gz"
        ),
        bed_tbi=os.path.join(
            config["OUT_FOLDER"], "GABi_results", "{species}", "{assembly}",
            "masking", "{assembly}.mask_from_file.bed.gz.tbi"
        )
    threads: cpu_func("masking_post")
    resources:
        mem_mb=mem_func("masking_post"),
        runtime=time_func("masking_post")
    container: CONTAINERS["gabi_base"]
    log:
        os.path.join(
            config["OUT_FOLDER"], "GABi_results", "{species}", "{assembly}",
            "logs", "C04_process_external_bed.{assembly}.log"
        )
    shell:
        r'''
        set -euo pipefail
        mkdir -p $(dirname {log})

        exec > "{log}" 2>&1

        echo "[GABi] Processing external BED file: {input.bed}"


        # TEMP DIRECTORY SETUP
        # -------------------------------------------------------------------
        WORK_DIR="$(gabi_get_workdir 10)"
        TEMP_DIR="$(mktemp -d "$WORK_DIR/GABi_extbed_{wildcards.species}_{wildcards.assembly}_XXXXXX")"
        trap 'rm -rf "$TEMP_DIR"' EXIT
        cd "$TEMP_DIR"

        # Process BED file - sort, merge, compress
        if [[ "{input.bed}" == *.gz ]]; then
            zcat "{input.bed}" | \
                awk 'BEGIN{{OFS="\t"}} NF>=3 {{print $1,$2,$3}}' | \
                sort -k1,1 -k2,2n | \
                bedtools merge -i - | \
                tee temp.bed | \
                bgzip -@ {threads} -c > "{output.bed_gz}"
        else
            cat "{input.bed}" | \
                awk 'BEGIN{{OFS="\t"}} NF>=3 {{print $1,$2,$3}}' | \
                sort -k1,1 -k2,2n | \
                bedtools merge -i - | \
                tee temp.bed | \
                bgzip -@ {threads} -c > "{output.bed_gz}"
        fi

        tabix -p bed "{output.bed_gz}"

        echo "[GABi] External BED processing complete: $(wc -l < temp.bed) regions"

        '''


rule C04_extract_softmask_bed:
    """Extract masked regions from softmasked assembly (fails if not softmasked)"""
    input:
        assembly=lambda w: get_assembly_path(w.species, w.assembly)
    output:
        bed=os.path.join(
            config["OUT_FOLDER"], "GABi_results", "{species}", "{assembly}",
            "masking", "{assembly}.mask_from_ref.bed.gz"
        ),
        tbi=os.path.join(
            config["OUT_FOLDER"], "GABi_results", "{species}", "{assembly}",
            "masking", "{assembly}.mask_from_ref.bed.gz.tbi"
        )
    params:
        script = str(scripts_dir / "softmasked_to_bed.py")
    threads: cpu_func("split_intervals")
    resources:
        mem_mb=mem_func("split_intervals"),
        runtime=time_func("split_intervals")
    container: CONTAINERS["gabi_base"]
    log:
        os.path.join(
            config["OUT_FOLDER"], "GABi_results", "{species}", "{assembly}",
            "logs", "C04_extract_softmask.{assembly}.log"
        )
    shell:
        r'''
        set -euo pipefail
        mkdir -p $(dirname {log})

        exec > {log} 2>&1

        echo "[GABi] Extracting masked regions from softmasked assembly"


        # TEMP DIRECTORY SETUP
        # -------------------------------------------------------------------
        WORK_DIR="$(gabi_get_workdir 25)"
        TEMP_DIR="$(mktemp -d "$WORK_DIR/GABi_softmask_{wildcards.species}_{wildcards.assembly}_XXXXXX")"
        trap 'rm -rf "$TEMP_DIR"' EXIT
        cd "$TEMP_DIR"

        # Copy input to local temp to avoid Lustre filesystem issues
        echo "[GABi] Copying assembly to temp directory..."
        if [[ "{input.assembly}" == *.gz ]]; then
            gunzip -c "{input.assembly}" > local_assembly.fa
        else
            cp "{input.assembly}" local_assembly.fa
        fi

        # Extract softmasked regions to BED
        echo "[GABi] Running softmasked_to_bed.py..."
        python {params.script} local_assembly.fa temp.bed

        # Check if any regions were found
        if [[ ! -s temp.bed ]]; then
            echo "[GABi] ------------------------------------------------------------"
            echo "[GABi] ❌ ERROR: No softmasked regions found in assembly!"
            echo "[GABi] ============================================================"
            echo "[GABi]"
            echo "[GABi] MASKING is set to 'auto' and no sibling BED file was found"
            echo "[GABi] next to the assembly, so this rule attempted to extract"
            echo "[GABi] masked regions directly from the FASTA — but the assembly"
            echo "[GABi] does not appear to be softmasked."
            echo "[GABi]"
            echo "[GABi] Please either:"
            echo "[GABi]   1. Set MASKING='on' in the control panel to run RepeatMasker"
            echo "[GABi]   2. Place a sibling BED file (same base name as the assembly,"
            echo "[GABi]      with extension .bed or .bed.gz) next to the assembly file"
            echo "[GABi]   3. Use a properly softmasked assembly with MASKING='auto'"
            echo "[GABi]   4. Set MASKING='off' to skip masking entirely"
            echo "[GABi]"
            exit 1
        fi

        # Sort, merge, and compress
        echo "[GABi] Sorting and merging regions..."
        sort -k1,1 -k2,2n temp.bed | \
        bedtools merge -i - | \
        tee merged.bed | \
        bgzip -@ {threads} -c > "{output.bed}"

        tabix -p bed "{output.bed}"

        echo "[GABi] Extracted $(wc -l < merged.bed) masked regions"

        '''
