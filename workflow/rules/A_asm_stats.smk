# ===============================================================================
# GABi - Assembly Stats Rules
# by Diego De Panis, 2026
# note: AI tools may have been used to improve, clean and/or comment this version of the code
# ===============================================================================

# -------------------------------------------------------------------------------
#  ORTHODB VERSION (compleasm, not a user setting)
# -------------------------------------------------------------------------------

# OrthoDB version used for gene completeness (same for compleasm and BUSCO)
# Also written to the A00 and A02 logs.
ODB_VERSION = "odb12.2"
if IS_MAIN and _as_bool(config.get("RUN_COMPL", True)):
    print(f"[GABi] OrthoDB release for compleasm: {ODB_VERSION}")


# -------------------------------------------------------------------------------
#  HELPER FUNCTIONS
# -------------------------------------------------------------------------------

# Report input collectors (GEP2-style: each returns [] if disabled) ------------

def get_report_gfastats_input(wildcards):
    """Gfastats output path for this assembly (always present)."""
    return os.path.join(
        config["OUT_FOLDER"], "GABi_results", wildcards.species,
        wildcards.assembly, "assembly_stats", "gfastats", "stats.txt"
    )


def get_report_compleasm_input(wildcards):
    """Compleasm tar.gz path (only if RUN_COMPL is enabled).
    The tar.gz contains {lineage}_odb{version}/full_table.tsv which is
    extracted at runtime for --compleasm-full."""
    if not _as_bool(config.get("RUN_COMPL", True)):
        return []
    return [os.path.join(
        config["OUT_FOLDER"], "GABi_results", wildcards.species,
        wildcards.assembly, "assembly_stats", "compleasm",
        "compleasm_results.tar.gz"
    )]


def get_all_report_inputs(wildcards):
    """Collect all inputs for the mini-report rule."""
    inputs = [get_report_gfastats_input(wildcards)]
    inputs.extend(get_report_compleasm_input(wildcards))
    return inputs


# -------------------------------------------------------------------------------
#  COMPLEASM EUKA DB AND PLACEMENT FILES (will be run only once!)
# ===============================================================================

rule A00_download_compleasm_db:
    output:
        flag=f"busco_lineages/eukaryota_{ODB_VERSION}.done",
        placement_flag=f"busco_lineages/placement_files_{ODB_VERSION}.done"
    params:
        db_dir="busco_lineages",
        odb=ODB_VERSION
    container: CONTAINERS["compleasm"]
    benchmark:
        "logs/A00_download_compleasm_db_benchmark.txt"
    log:
        "logs/A00_download_db.log"
    threads: cpu_func("download_db")
    resources:
        mem_mb=mem_func("download_db"),
        runtime=time_func("download_db")
    shell:
        """
        mkdir -p $(dirname {log})
        exec > {log} 2>&1

        echo "[GABi] OrthoDB release: {params.odb}"
        mkdir -p {params.db_dir}
        compleasm download eukaryota --odb {params.odb} --library_path {params.db_dir}
        touch {output.flag}
        touch {output.placement_flag}
        """


# -------------------------------------------------------------------------------
#  ASM STATS
# ===============================================================================

rule A01_run_gfastats:
    input:
        assembly=lambda wildcards: get_assembly_path(wildcards.species, wildcards.assembly)
    output:
        stats=os.path.join(config["OUT_FOLDER"], "GABi_results", "{species}", "{assembly}", "assembly_stats", "gfastats", "stats.txt")
    container: CONTAINERS["gabi_base"]
    benchmark:
        os.path.join(config["OUT_FOLDER"], "GABi_results", "{species}", "{assembly}", "logs", "A01_gfastats_benchmark.txt")
    log:
        os.path.join(config["OUT_FOLDER"], "GABi_results", "{species}", "{assembly}", "logs", "A01_gfastats.log")
    threads: cpu_func("gfastats")
    resources:
        mem_mb=mem_func("gfastats"),
        runtime=time_func("gfastats")
    shell:
        """
        mkdir -p $(dirname {output.stats})
        mkdir -p $(dirname {log})
        gfastats -f {input.assembly} -j {threads} --nstar-report > {output.stats} 2> {log}
        """


rule A02_run_compleasm:
    input:
        assembly=lambda w: get_assembly_path(w.species, w.assembly),
        db_flag=f"busco_lineages/eukaryota_{ODB_VERSION}.done",
        placement_flag=f"busco_lineages/placement_files_{ODB_VERSION}.done"
    output:
        summary=os.path.join(
            config["OUT_FOLDER"], "GABi_results", "{species}", "{assembly}",
            "assembly_stats", "compleasm", "summary.txt"
        ),
        targz=os.path.join(
            config["OUT_FOLDER"], "GABi_results", "{species}", "{assembly}",
            "assembly_stats", "compleasm", "compleasm_results.tar.gz"
        )
    container: CONTAINERS["compleasm"]
    params:
        outdir=os.path.join(
            config["OUT_FOLDER"], "GABi_results", "{species}", "{assembly}",
            "assembly_stats", "compleasm"
        ),
        shared_db="busco_lineages",
        lineage=config.get("LINEAGE", "eukaryota"),
        odb=ODB_VERSION
    benchmark:
        os.path.join(
            config["OUT_FOLDER"], "GABi_results", "{species}", "{assembly}",
            "logs", "A02_compleasm_benchmark.txt"
        )
    log:
        os.path.join(
            config["OUT_FOLDER"], "GABi_results", "{species}", "{assembly}",
            "logs", "A02_compleasm.log"
        )
    threads: cpu_func("compleasm")
    resources:
        mem_mb=mem_func("compleasm"),
        runtime=time_func("compleasm")
    shell:
        r'''
        set -euo pipefail
        mkdir -p "{params.outdir}" "$(dirname {log})"

        exec > "{log}" 2>&1

        LINEAGE="{params.lineage}"
        echo "[GABi] Configured LINEAGE: $LINEAGE"
        echo "[GABi] OrthoDB release:   {params.odb}"

        # TEMP DIRECTORY SETUP
        # -------------------------------------------------------------------
        WORK_DIR="$(gabi_get_workdir 50)"
        TEMP_DIR="$(mktemp -d "$WORK_DIR/GABi_compleasm_{wildcards.species}_{wildcards.assembly}_XXXXXX")"
        trap 'rm -rf "$TEMP_DIR"' EXIT

        TMPDB="$TEMP_DIR/compleasm_db"
        mkdir -p "$TMPDB"

        echo "[GABi] Preparing isolated lineage DB in $TMPDB"
        # Using cp -rL instead of tar piping to avoid "short read" errors on some filesystems
        # This brings in eukaryota_{params.odb} + placement files (for --autolineage)
        cp -rL "{params.shared_db}"/* "$TMPDB/"

        # FETCH ADDITIONAL LINEAGE IF NEEDED (lives only in $TMPDB, dies with trap)
        # -------------------------------------------------------------------
        if [ "$LINEAGE" != "eukaryota" ] && [ "$LINEAGE" != "auto" ]; then
            echo "[GABi] Downloading additional lineage '$LINEAGE' into temp DB"
            compleasm download "$LINEAGE" --odb {params.odb} --library_path "$TMPDB"
        fi

        # BUILD LINEAGE FLAG
        # -------------------------------------------------------------------
        if [ "$LINEAGE" = "auto" ]; then
            LINEAGE_FLAG="--autolineage"
            echo "[GABi] Using --autolineage (placement files required)"
        else
            LINEAGE_FLAG="--lineage $LINEAGE"
            echo "[GABi] Using --lineage $LINEAGE"
        fi

        # RUN COMPLEASM
        # -------------------------------------------------------------------
        TMPOUT="$TEMP_DIR/output"
        mkdir -p "$TMPOUT"

        compleasm run \
          -a "{input.assembly}" \
          -o "$TMPOUT" \
          -t {threads} \
          --library_path "$TMPDB" \
          --odb {params.odb} \
          $LINEAGE_FLAG

        # PACKAGE OUTPUTS
        # -------------------------------------------------------------------
        echo "[GABi] Packaging outputs"

        # Archive must live OUTSIDE the tree being archived
        ARCHIVE_TMP="$TEMP_DIR/compleasm_results.tar.gz"

        cd "$TMPOUT"
        tar -czf "$ARCHIVE_TMP" --exclude='summary.txt' .

        gzip -t "$ARCHIVE_TMP"
        tar -tzf "$ARCHIVE_TMP" >/dev/null
        echo "[GABi] Archive verified"

        # Locate summary.txt
        if [ -f "summary.txt" ]; then
            SUMMARY_SRC="summary.txt"
            echo "[GABi] Found summary.txt at top level"
        else
            echo "[GABi] Warning: summary.txt not at top level, searching..."
            SUMMARY_SRC=$(find . -name "summary.txt" -type f | head -1 || true)
            if [ -z "$SUMMARY_SRC" ]; then
                echo "[GABi] ERROR: No summary.txt found anywhere in output"
                ls -la
                exit 1
            fi
            echo "[GABi] Found summary at: $SUMMARY_SRC"
        fi

        # Copy only the two final artifacts to the output directory
        mkdir -p "{params.outdir}"
        cp "$SUMMARY_SRC" "{output.summary}"
        cp "$ARCHIVE_TMP" "{output.targz}"

        echo "[GABi] Done: archive at {output.targz}"
        '''


rule A03_create_mini_report:
    input:
        deps = get_all_report_inputs
    output:
        report = os.path.join(
            config["OUT_FOLDER"], "GABi_results", "{species}", "{assembly}",
            "assembly_stats", "{assembly}_mini_report.md"
        )
    params:
        species   = lambda w: w.species,
        assembly  = lambda w: w.assembly,
        gfastats  = lambda w: get_report_gfastats_input(w),
        compleasm = lambda w: get_report_compleasm_input(w),
        lineage   = config.get("LINEAGE", "eukaryota"),
        script_path = str(scripts_dir / "make_gep2_report.py")
    container: CONTAINERS["gabi_base"]
    benchmark:
        os.path.join(config["OUT_FOLDER"], "GABi_results", "{species}", "{assembly}", "logs", "A03_mini_report_benchmark.txt")
    log:
        os.path.join(config["OUT_FOLDER"], "GABi_results", "{species}", "{assembly}", "logs", "A03_mini_report.log")
    threads: cpu_func("light_task")
    resources:
        mem_mb=mem_func("light_task"),
        runtime=time_func("light_task")
    shell:
        """
        mkdir -p $(dirname {log})
        exec > {log} 2>&1

        cmd="python {params.script_path} -s {params.species} -a {params.assembly} -g {params.gfastats}"

        # Compleasm: extract full_table.tsv from tar.gz files (if any)
        COMPLEASM_FULLS=""
        COMPLEASM_CLEANUP=""
        LINEAGE="{params.lineage}"

        for targz in {params.compleasm}; do
            if [ -n "$targz" ]; then
                extract_dir=$(dirname "$targz")
                tar -xzf "$targz" -C "$extract_dir"

                # Pick the right lineage's full_table.tsv:
                #   - specific lineage (e.g., aves): find that one
                #   - auto: prefer non-eukaryota (the auto-detected result)
                #   - eukaryota: find eukaryota
                if [ "$LINEAGE" = "auto" ]; then
                    # Prefer the non-eukaryota lineage (the interesting auto-detected one)
                    full_table=$(find "$extract_dir" -name "full_table.tsv" -path "*_odb*" \
                        ! -path "*eukaryota_odb*" | head -1)
                    # Fallback to eukaryota if no other lineage found
                    if [ -z "$full_table" ]; then
                        full_table=$(find "$extract_dir" -name "full_table.tsv" -path "*eukaryota_odb*" | head -1)
                    fi
                else
                    # Specific lineage: find exactly that one
                    full_table=$(find "$extract_dir" -name "full_table.tsv" -path "*${{LINEAGE}}_odb*" | head -1)
                fi

                if [ -n "$full_table" ]; then
                    COMPLEASM_FULLS="$COMPLEASM_FULLS $full_table"
                    echo "[GABi] Using compleasm full_table.tsv from: $(dirname "$full_table")"
                else
                    echo "[GABi] WARNING: No full_table.tsv found for lineage '$LINEAGE'"
                fi

                # Track ALL extracted odb directories for cleanup
                for odb_dir in $(find "$extract_dir" -maxdepth 1 -type d -name "*_odb*"); do
                    COMPLEASM_CLEANUP="$COMPLEASM_CLEANUP $odb_dir"
                done
            fi
        done

        if [ -n "$COMPLEASM_FULLS" ]; then
            cmd="$cmd --compleasm-full $COMPLEASM_FULLS"
        fi

        cmd="$cmd -o {output.report}"

        echo "[GABi] Command: $cmd"
        $cmd

        # Footer: point to GEP2 for a full QC of the reference
        (
            echo ""
            echo "<br>"
            echo ""
            echo "*This is a quick check of the reference assembly. For a thorough quality control (k-mer QV, contamination screening and more), run a full [GEP2](https://github.com/diegomics/GEP2) analysis.*"
        ) >> {output.report}

        # Cleanup ALL extracted compleasm directories (eukaryota + target lineage)
        for cleanup_dir in $COMPLEASM_CLEANUP; do
            if [ -d "$cleanup_dir" ]; then
                rm -rf "$cleanup_dir"
                echo "[GABi] Cleaned up extracted compleasm: $cleanup_dir"
            fi
        done
        """
