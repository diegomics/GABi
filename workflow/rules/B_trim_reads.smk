# ===============================================================================
# GABi - Read Trimming and QC Rules
# by Diego De Panis, 2026
# note: AI tools may have been used to improve, clean and/or comment this version of the code
# ===============================================================================
#
# Read processing is SAMPLE-CENTRIC (not assembly-specific):
#   - Raw reads are symlinked to: GABi_results/data/{sample_id}/{read_type}/
#   - Processed reads go to:      GABi_results/data/{sample_id}/{read_type}/processed/
#   - QC reports go to:           GABi_results/data/{sample_id}/{read_type}/reports/
#   - Logs go to:                 GABi_results/data/{sample_id}/{read_type}/logs/
#
# This allows the same processed reads to be used with multiple assemblies
# without redundant processing.
# ===============================================================================


# Disambiguate producers that yield the same reads/* targets
ruleorder: B00_link_long_reads > B01_compress_long_reads
ruleorder: B00_link_pe_reads > B01_compress_pe_reads


# -------------------------------------------------------------------------------
#  PATH TEMPLATES
# -------------------------------------------------------------------------------

# Base path for sample data (sample-centric, no species/assembly)
_DATA_BASE = os.path.join(config["OUT_FOLDER"], "GABi_results", "data")


# -------------------------------------------------------------------------------
#  HELPER FUNCTIONS (Sample-Centric)
# -------------------------------------------------------------------------------

def _find_sample_reads(sample_id, read_type_key):
    """
    Find read files for a sample, across every species/assembly in the config.
    Returns: (normalized_read_type, reads_dict) or (None, {}) if not found.

    This used to return the FIRST (species, assembly) node containing the
    sample. _collect_sample_reads() gathers them all instead and raises when
    two assemblies declare different files for the same sample and technology,
    so the answer no longer depends on dict iteration order.
    """
    target_norm = normalize_read_type(read_type_key)
    found = _collect_sample_reads(sample_id, read_type_key)
    if target_norm in found:
        return target_norm, found[target_norm]
    return None, {}


def _pick_group_sample(w):
    """Pick the read group based on sample_id and read_type wildcards (no species/assembly needed)"""
    rt, reads = _find_sample_reads(w.sample_id, w.read_type)

    if rt is None:
        raise ValueError(f"No valid read data for sample {w.sample_id} with read_type={w.read_type}")

    expected_rt = normalize_read_type(w.read_type)
    if rt != expected_rt:
        raise ValueError(f"read_type mismatch: expected {expected_rt}, got {rt}")

    for grp in _enumerate_groups(rt, reads):
        if str(grp["idx"]) == str(w.idx):
            expected_base = grp["base"]
            if expected_base != w.base:
                raise ValueError(
                    f"Base mismatch for {w.sample_id} "
                    f"Path{w.idx}: expected {expected_base}, got {w.base}"
                )
            return rt, grp

    raise ValueError(f"Could not find group Path{w.idx} for {w.sample_id}/{w.read_type}")


def _mk_dirs(log):
    os.makedirs(os.path.dirname(str(log)), exist_ok=True)


def _fastp_adapter_args():
    """fastp options for the configured ADAPTERS_PE (validated at load)."""
    mode, fasta = adapters_pe_mode()
    if mode == "auto":
        return "--detect_adapter_for_pe"
    if mode == "overlap":
        return ""
    return f'--adapter_fasta "{fasta}"'


# A guess fastp took from its built-in list is printed with its name line
# (">Illumina TruSeq Adapter Read 1") before the sequence; a guess it made up
# from over-represented sequence is printed bare. That bare line is the signal.
_FASTP_DETECT_RE = re.compile(r"^Detecting adapter sequence for (read[12])\.\.\.\s*$")


def _tandem_repeat_unit(seq, max_period=12, max_mismatch=0.1):
    """Repeat unit if seq is a short tandem repeat (telomere, satellite, ...), else None.
    """
    seq = seq.upper()
    for k in range(1, max_period + 1):
        n = len(seq) - k
        if n >= 2 * k and sum(seq[i] != seq[i + k] for i in range(n)) <= max_mismatch * n:
            return seq[:k]
    return None


def _fastp_unknown_adapter_guesses(text):
    """[(read, sequence)] for each adapter fastp guessed without matching its known list."""
    lines = text.splitlines()
    out = []
    for i, line in enumerate(lines):
        m = _FASTP_DETECT_RE.match(line.strip())
        if not m:
            continue
        nxt = lines[i + 1].strip() if i + 1 < len(lines) else ""
        if not nxt or nxt.startswith(">") or nxt.startswith("No adapter detected"):
            continue
        out.append((m.group(1), nxt))
    return out


_MISSING_IN = os.path.join(
    config["OUT_FOLDER"], "GABi_results", ".missing_input_placeholder_DO_NOT_CREATE"
)


# -------------------------------------------------------------------------------
#  INPUT FUNCTIONS (Sample-Centric)
# -------------------------------------------------------------------------------

def _linkable_long_src(w):
    """Check if long read source is already gzipped (can be symlinked)"""
    rt, g = _pick_group_sample(w)
    return g["r"] if str(g["r"]).endswith(".gz") else _MISSING_IN


def _compressible_long_src(w):
    """Check if long read source needs compression"""
    rt, g = _pick_group_sample(w)
    return g["r"] if not str(g["r"]).endswith(".gz") else _MISSING_IN


def _linkable_pe_r1(w):
    """Check if PE R1 can be symlinked (not trimming and already gzipped)"""
    if _as_bool(config.get("TRIM_PE", True)):
        return _MISSING_IN
    rt, g = _pick_group_sample(w)
    return g["r1"] if str(g["r1"]).endswith(".gz") else _MISSING_IN


def _linkable_pe_r2(w):
    """Check if PE R2 can be symlinked (not trimming and already gzipped)"""
    if _as_bool(config.get("TRIM_PE", True)):
        return _MISSING_IN
    rt, g = _pick_group_sample(w)
    return g["r2"] if str(g["r2"]).endswith(".gz") else _MISSING_IN


def _compressible_pe_r1(w):
    """Check if PE R1 needs compression"""
    if _as_bool(config.get("TRIM_PE", True)):
        return _MISSING_IN
    rt, g = _pick_group_sample(w)
    needs = (not str(g["r1"]).endswith(".gz")) or (not str(g["r2"]).endswith(".gz"))
    return g["r1"] if needs else _MISSING_IN


def _compressible_pe_r2(w):
    """Check if PE R2 needs compression"""
    if _as_bool(config.get("TRIM_PE", True)):
        return _MISSING_IN
    rt, g = _pick_group_sample(w)
    needs = (not str(g["r1"]).endswith(".gz")) or (not str(g["r2"]).endswith(".gz"))
    return g["r2"] if needs else _MISSING_IN


def _get_ont_reads_for_correction(w):
    """Get ONT reads that need correction"""
    if not _as_bool(config.get("CORRECT_ONT", False)):
        return _MISSING_IN

    rt, reads = _find_sample_reads(w.sample_id, "ont")
    if rt != "ont":
        return _MISSING_IN

    # Return the intermediate linked/compressed ONT file
    for grp in _enumerate_groups(rt, reads):
        if str(grp["idx"]) == str(w.idx) and grp["base"] == w.base:
            return os.path.join(
                _DATA_BASE, w.sample_id, "ont",
                f"Path{w.idx}_{w.base}.fq.gz"
            )

    return _MISSING_IN


def _get_hifi_input_for_filtering(w):
    """Get linked/compressed HiFi reads for filtering"""
    return os.path.join(
        _DATA_BASE, w.sample_id, "hifi",
        f"Path{w.idx}_{w.base}.fq.gz"
    )


def _get_fastqc_input(wildcards):
    """Get the correct input file for FastQC based on processing settings"""
    filename = wildcards.filename
    sample_id = wildcards.sample_id
    read_type = wildcards.read_type

    # Check if it's a PE read (has _R1 or _R2)
    if "_R1" in filename or "_R2" in filename:
        if _as_bool(config.get("TRIM_PE", True)):
            # FastQC on original PE reads (before trimming)
            match = re.match(r"Path(\d+)_(.+)_(R[12])$", filename)
            if match:
                idx, base, read_num = match.groups()
                rt, reads = _find_sample_reads(sample_id, read_type)
                for grp in _enumerate_groups(rt, reads):
                    if str(grp["idx"]) == idx and grp["base"] == base:
                        if read_num == "R1":
                            return grp["r1"]
                        else:
                            return grp["r2"]

    # Default: use the file in the reads directory
    return os.path.join(_DATA_BASE, sample_id, read_type, f"{filename}.fq.gz")


def _get_nanoplot_input(w):
    if w.suffix in ("_corrected", "_filtered"):
        return os.path.join(
            _DATA_BASE, w.sample_id, w.read_type, "processed",
            f"Path{w.idx}_{w.base}{w.suffix}.fq.gz"
        )
    return os.path.join(
        _DATA_BASE, w.sample_id, w.read_type,
        f"Path{w.idx}_{w.base}{w.suffix}.fq.gz"
    )

# -------------------------------------------------------------------------------
#  WILDCARD CONSTRAINTS
# -------------------------------------------------------------------------------

wildcard_constraints:
    sample_id = r"[^/]+",
    read_type = r"(hifi|ont|pacbio|nanopore|illumina|10x)",
    idx = r"\d+",
    base = r"[^/]+"


# -------------------------------------------------------------------------------
#  LONG READS - LINK/COMPRESS
# ===============================================================================

rule B00_link_long_reads:
    """Symlink compressed long reads to sample data directory"""
    input:
        src = _linkable_long_src
    output:
        linked = os.path.join(
            _DATA_BASE, "{sample_id}", "{read_type}",
            "Path{idx}_{base}.fq.gz"
        )
    wildcard_constraints:
        read_type = r"(hifi|ont|pacbio|nanopore)",
        base = r"[^/]+(?<!_corrected)(?<!_filtered)"
    threads: 1
    resources:
        mem_mb = mem_func("light_task"),
        runtime = time_func("light_task")
    log:
        os.path.join(_DATA_BASE, "{sample_id}", "{read_type}", "logs",
                     "B00_link_long_reads.Path{idx}.{base}.log")
    run:
        rt, grp = _pick_group_sample(wildcards)
        if grp["kind"] != "long":
            raise ValueError("This target belongs to a paired-end job, not long reads.")
        if input.src == _MISSING_IN:
            raise ValueError("link_long_reads should only run when the source is .gz")
        _mk_dirs(log)

        src = input.src
        dest = output.linked
        os.makedirs(os.path.dirname(dest), exist_ok=True)

        cmd = f'ln -sf "$(realpath "{src}")" "{dest}" && echo "Linked {src} to {dest}" > {log} 2>&1'
        shell(cmd)


rule B01_compress_long_reads:
    """Compress uncompressed long reads to sample data directory"""
    input:
        src = _compressible_long_src
    output:
        compressed = os.path.join(
            _DATA_BASE, "{sample_id}", "{read_type}",
            "Path{idx}_{base}.fq.gz"
        )
    wildcard_constraints:
        read_type = r"(hifi|ont|pacbio|nanopore)",
        base = r"[^/]+(?<!_corrected)(?<!_filtered)"
    threads: cpu_func("compress_reads")
    resources:
        mem_mb = mem_func("compress_reads"),
        runtime = time_func("compress_reads")
    container: CONTAINERS["gabi_base"]
    log:
        os.path.join(_DATA_BASE, "{sample_id}", "{read_type}", "logs",
                     "B01_compress_long_reads.Path{idx}.{base}.log")
    run:
        rt, grp = _pick_group_sample(wildcards)
        if grp["kind"] != "long":
            raise ValueError("This target belongs to a paired-end job, not long reads.")
        if input.src == _MISSING_IN:
            raise ValueError("compress_long_reads should only run when source is NOT .gz")
        _mk_dirs(log)

        os.makedirs(os.path.dirname(output.compressed), exist_ok=True)
        cmd = f'pigz -p {threads} -c "{input.src}" > "{output.compressed}" 2> {log}'
        shell(cmd)


# -------------------------------------------------------------------------------
#  ONT CORRECTION (if CORRECT_ONT=True)
# ===============================================================================

rule B03_correct_ont_reads:
    """Error-correct ONT reads using hifiasm"""
    input:
        reads = _get_ont_reads_for_correction
    output:
        corrected = os.path.join(
            _DATA_BASE, "{sample_id}", "ont", "processed",
            "Path{idx}_{base}_corrected.fq.gz"
        )
    wildcard_constraints:
        base = r"[^/]+(?<!_corrected)"
    threads: cpu_func("ont_correction")
    resources:
        mem_mb = mem_func("ont_correction"),
        runtime = time_func("ont_correction")
    container: CONTAINERS["gabi_base"]
    log:
        os.path.join(_DATA_BASE, "{sample_id}", "ont", "logs",
                     "B03_correct_ont.Path{idx}.{base}.log")
    shell:
        r'''
        set -euo pipefail
        mkdir -p "$(dirname {output.corrected})" "$(dirname {log})"

        exec > {log} 2>&1

        if [[ "{input.reads}" == *"_placeholder_"* ]]; then
            echo "ONT correction skipped (CORRECT_ONT=False or not ONT reads)"
            exit 1
        fi

        echo "[GABi] Starting ONT error correction with hifiasm"
        echo "[GABi] Input: {input.reads}"


        # TEMP DIRECTORY SETUP
        # -------------------------------------------------------------------
        WORK_DIR="$(gabi_get_workdir 20)"
        TEMP_DIR="$(mktemp -d "$WORK_DIR/GABi_ont_correct_{wildcards.sample_id}_{wildcards.idx}_XXXXXX")"
        trap 'rm -rf "$TEMP_DIR"' EXIT
        cd "$TEMP_DIR"

        echo "[GABi] Working directory: $TEMP_DIR"

        hifiasm -o tmp --ont -t {threads} -l0 --write-ec {input.reads} >> hifi.log 2>&1 &
        HIFI_PID=$!

        LAST_SIZE=0
        STAGNANT_COUNT=0
        LOOP_COUNT=0

        while kill -0 $HIFI_PID 2>/dev/null; do
            LOOP_COUNT=$((LOOP_COUNT + 1))

            if grep -q "Reads has been written" hifi.log 2>/dev/null; then
                echo "Log indicates EC reads complete" >> {log}
                sleep 30
                kill -TERM $HIFI_PID 2>/dev/null || true
                break
            fi

            if [[ -f "tmp.ec.fq" ]]; then
                CURRENT_SIZE=$(stat -c%s "tmp.ec.fq" 2>/dev/null || echo 0)
                echo "Loop $LOOP_COUNT: EC file size = $CURRENT_SIZE bytes" >> {log}

                if [[ $CURRENT_SIZE -eq $LAST_SIZE ]] && [[ $CURRENT_SIZE -gt 1000000 ]]; then
                    STAGNANT_COUNT=$((STAGNANT_COUNT + 1))
                    echo "File size unchanged for $((STAGNANT_COUNT * 30))s" >> {log}

                    if [[ $STAGNANT_COUNT -ge 4 ]]; then
                        echo "File stagnant for 120s, terminating hifiasm" >> {log}
                        kill -TERM $HIFI_PID 2>/dev/null || true
                        break
                    fi
                else
                    STAGNANT_COUNT=0
                    LAST_SIZE=$CURRENT_SIZE
                fi
            fi

            if [[ $LOOP_COUNT -gt 2880 ]]; then
                echo "Timeout reached (24h), terminating" >> {log}
                kill -TERM $HIFI_PID 2>/dev/null || true
                break
            fi

            sleep 30
        done

        wait $HIFI_PID || EXIT_CODE=$?
        echo "Hifiasm exit code: ${{EXIT_CODE:-0}}" >> {log}

        if [[ -f "tmp.ec.fq" ]] && [[ $(stat -c%s "tmp.ec.fq") -gt 1000000 ]]; then
            echo "Compressing corrected reads..." >> {log}
            pigz -p {threads} -c tmp.ec.fq > {output.corrected}
            echo "Successfully created corrected reads" >> {log}
        else
            echo "ERROR: No valid EC file produced" >> {log}
            ls -la >> {log}
            exit 1
        fi

        echo "[GABi] ONT correction complete"
        '''


# -------------------------------------------------------------------------------
#  HIFI PROCESSING
# ===============================================================================

rule B04_check_uli_primers:
    """Check for ULI primers in HiFi reads"""
    input:
        reads = os.path.join(
            _DATA_BASE, "{sample_id}", "hifi",
            "Path{idx}_{base}.fq.gz"
        )
    output:
        yaml = os.path.join(
            _DATA_BASE, "{sample_id}", "hifi", "reports",
            "Path{idx}_{base}_uli_mqc.yaml"
        )
    wildcard_constraints:
        # As in B00/B01/B03/B05: without the guard, base can bind a name that
        # ends in a processed-file suffix.
        base = r"[^/]+(?<!_filtered)"
    threads: cpu_func("reads_qc")
    resources:
        mem_mb = mem_func("reads_qc"),
        runtime = time_func("reads_qc")
    container: CONTAINERS["gabi_base"]
    log:
        os.path.join(_DATA_BASE, "{sample_id}", "hifi", "logs",
                     "B04_check_uli.Path{idx}.{base}.log")
    shell:
        r'''
        set -euo pipefail
        mkdir -p "$(dirname {output.yaml})" "$(dirname {log})"
        exec > {log} 2>&1

        echo "Checking for ULI primers in HiFi reads..."

        # Work in temp directory
        WORK_DIR="$(gabi_get_workdir 5)"
        TEMP_DIR="$(mktemp -d "$WORK_DIR/GABi_uli_check_{wildcards.sample_id}_XXXXXX")"
        trap 'rm -rf "$TEMP_DIR"' EXIT
        cd "$TEMP_DIR"

        # Subsample reads
        seqtk sample -s789 {input.reads} 100000 > subset.fq

        # Count ULI primer occurrences
        hits=$(awk 'NR%4==2' subset.fq | grep -E -c 'AAGCAGTGGTATCAACGCAGAGTACT|AGTACTCTGCGTTGATACCACTGCTT' || echo 0)
        reads=$(awk 'NR%4==2' subset.fq | wc -l)
        ulipct=$(awk -v h="$hits" -v r="$reads" 'BEGIN{{printf("%.2f", r>0 ? 100*h/r : 0)}}')

        echo "ULI primer hits: $hits / $reads reads ($ulipct%)"

        # Generate MultiQC custom content YAML
        cat > {output.yaml} << EOF
id: 'uli_primer_check'
section_name: 'ULI Primer Check'
description: 'Detection of ULI adapter sequences in HiFi reads'
plot_type: 'generalstats'
pconfig:
    - uli_pct:
        title: 'ULI %'
        description: 'Percentage of reads containing ULI primer sequences'
        min: 0
        max: 100
        suffix: '%'
        format: '{{:,.2f}}'
data:
    {wildcards.sample_id}_Path{wildcards.idx}:
        uli_pct: $ulipct
EOF

        echo "[GABi] ULI check complete"
        '''


rule B05_filter_hifi_adapters:
    """Filter HiFi reads containing SMRTbell adapter sequences using bbduk"""
    input:
        reads = _get_hifi_input_for_filtering,
        adapters = workflow.source_path("../scripts/hifi_blocklist.fa")
    output:
        filtered = os.path.join(
            _DATA_BASE, "{sample_id}", "hifi", "processed",
            "Path{idx}_{base}_filtered.fq.gz"
        ),
        stats = os.path.join(
            _DATA_BASE, "{sample_id}", "hifi", "reports",
            "Path{idx}_{base}_bbduk.stats"
        )
    wildcard_constraints:
        base = r"[^/]+(?<!_filtered)"
    threads: cpu_func("filter_hifi")
    resources:
        mem_mb = mem_func("filter_hifi"),
        runtime = time_func("filter_hifi")
    container: CONTAINERS["gabi_base"]
    log:
        os.path.join(_DATA_BASE, "{sample_id}", "hifi", "logs",
                     "B05_filter_hifi.Path{idx}.{base}.log")
    shell:
        r'''
        set -euo pipefail
        mkdir -p "$(dirname {output.filtered})" "$(dirname {output.stats})" "$(dirname {log})"
        exec > {log} 2>&1

        echo "[GABi] Starting HiFi adapter filtering with bbduk"
        echo "[GABi] Input: {input.reads}"
        echo "[GABi] Ref:   {input.adapters}"

        bbduk.sh -Xmx$(({resources.mem_mb} * 85 / 100))m \
            in={input.reads} \
            out={output.filtered} \
            ref={input.adapters} \
            stats={output.stats} \
            k=27 hdist=1 rcomp=t \
            threads={threads} \
            overwrite=t

        echo "[GABi] HiFi filtering complete"
        '''


# -------------------------------------------------------------------------------
#  PE READS - LINK/COMPRESS
# ===============================================================================

rule B00_link_pe_reads:
    """Symlink compressed PE reads (only if TRIM_PE=False)"""
    input:
        r1 = _linkable_pe_r1,
        r2 = _linkable_pe_r2
    output:
        r1 = os.path.join(
            _DATA_BASE, "{sample_id}", "{read_type}",
            "Path{idx}_{base}_R1.fq.gz"
        ),
        r2 = os.path.join(
            _DATA_BASE, "{sample_id}", "{read_type}",
            "Path{idx}_{base}_R2.fq.gz"
        )
    wildcard_constraints:
        read_type = r"(illumina|10x)",
        base = r"[^/]+(?<!_trimmed)"
    threads: 1
    resources:
        mem_mb = mem_func("light_task"),
        runtime = time_func("light_task")
    log:
        os.path.join(_DATA_BASE, "{sample_id}", "{read_type}", "logs",
                     "B00_link_pe_reads.Path{idx}.{base}.log")
    run:
        rt, grp = _pick_group_sample(wildcards)
        if grp["kind"] != "pe":
            raise ValueError("This target belongs to a long-read job, not PE.")
        if input.r1 == _MISSING_IN or input.r2 == _MISSING_IN:
            raise ValueError("link_pe_reads should only run when TRIM_PE=False and source is .gz")
        _mk_dirs(log)

        os.makedirs(os.path.dirname(output.r1), exist_ok=True)
        cmd = (
            f'ln -sf "$(realpath "{input.r1}")" "{output.r1}" && '
            f'ln -sf "$(realpath "{input.r2}")" "{output.r2}" && '
            f'echo "Linked PE reads" > {log} 2>&1'
        )
        shell(cmd)


rule B01_compress_pe_reads:
    """Compress PE reads (only if TRIM_PE=False and source not gzipped)"""
    input:
        r1 = _compressible_pe_r1,
        r2 = _compressible_pe_r2
    output:
        r1 = os.path.join(
            _DATA_BASE, "{sample_id}", "{read_type}",
            "Path{idx}_{base}_R1.fq.gz"
        ),
        r2 = os.path.join(
            _DATA_BASE, "{sample_id}", "{read_type}",
            "Path{idx}_{base}_R2.fq.gz"
        )
    wildcard_constraints:
        read_type = r"(illumina|10x)",
        base = r"[^/]+(?<!_trimmed)"
    threads: cpu_func("compress_reads")
    resources:
        mem_mb = mem_func("compress_reads"),
        runtime = time_func("compress_reads")
    container: CONTAINERS["gabi_base"]
    log:
        os.path.join(_DATA_BASE, "{sample_id}", "{read_type}", "logs",
                     "B01_compress_pe_reads.Path{idx}.{base}.log")
    run:
        rt, grp = _pick_group_sample(wildcards)
        if grp["kind"] != "pe":
            raise ValueError("This target belongs to a long-read job, not PE.")
        if input.r1 == _MISSING_IN or input.r2 == _MISSING_IN:
            raise ValueError("compress_pe_reads should only run when TRIM_PE=False and source is NOT .gz")
        _mk_dirs(log)

        os.makedirs(os.path.dirname(output.r1), exist_ok=True)
        cmd = (
            f'if [[ "{input.r1}" == *.gz ]]; then ln -sf "$(realpath "{input.r1}")" "{output.r1}"; '
            f'else pigz -p {threads} -c "{input.r1}" > "{output.r1}"; fi && '
            f'if [[ "{input.r2}" == *.gz ]]; then ln -sf "$(realpath "{input.r2}")" "{output.r2}"; '
            f'else pigz -p {threads} -c "{input.r2}" > "{output.r2}"; fi 2> {log}'
        )
        shell(cmd)


# -------------------------------------------------------------------------------
#  PE READS - TRIMMING (if TRIM_PE=True)
# ===============================================================================

rule B02_trim_pe_reads:
    """Trim PE reads with fastp (only if TRIM_PE=True)"""
    input:
        r1 = lambda w: _pick_group_sample(w)[1]["r1"],
        r2 = lambda w: _pick_group_sample(w)[1]["r2"]
    output:
        r1 = os.path.join(
            _DATA_BASE, "{sample_id}", "{read_type}", "processed",
            "Path{idx}_{base}_R1_trimmed.fq.gz"
        ),
        r2 = os.path.join(
            _DATA_BASE, "{sample_id}", "{read_type}", "processed",
            "Path{idx}_{base}_R2_trimmed.fq.gz"
        ),
        json = os.path.join(
            _DATA_BASE, "{sample_id}", "{read_type}", "reports",
            "Path{idx}_{base}_fastp.json"
        ),
        html = os.path.join(
            _DATA_BASE, "{sample_id}", "{read_type}", "reports",
            "Path{idx}_{base}_fastp.html"
        )
    wildcard_constraints:
        read_type = r"(illumina|10x)",
        base = r"[^/]+"
    threads: cpu_func("trim_pe")
    resources:
        mem_mb = mem_func("trim_pe"),
        runtime = time_func("trim_pe")
    container: CONTAINERS["gabi_base"]
    log:
        os.path.join(_DATA_BASE, "{sample_id}", "{read_type}", "logs",
                     "B02_trim_pe_reads.Path{idx}.{base}.log")
    params:
        adapter_mode = lambda w: config.get("ADAPTERS_PE", "illumina"),
        adapter_args = lambda w: _fastp_adapter_args()
    run:
        rt, grp = _pick_group_sample(wildcards)
        if grp["kind"] != "pe":
            raise ValueError("This target belongs to a long-read job, not PE.")
        _mk_dirs(log)

        os.makedirs(os.path.dirname(output.r1), exist_ok=True)
        os.makedirs(os.path.dirname(output.json), exist_ok=True)

        tenx = _is_10x(rt)
        extra = "--trim_front1 23 --trim_front2 0" if tenx else ""

        log_path = str(log)
        with open(log_path, "a") as fh:
            fh.write(f"[GABi] ADAPTERS_PE: {params.adapter_mode}\n")
        # The log is appended to, so it may hold earlier attempts. Remember where
        # this one starts, so the check below reads only this run's output.
        start = os.path.getsize(log_path)

        cmd = (
            f'fastp -i "{input.r1}" -I "{input.r2}" -o "{output.r1}" -O "{output.r2}" -w {threads} '
            f'{params.adapter_args} -l 50 --overrepresentation_analysis {extra} '
            f'-j "{output.json}" -h "{output.html}" >> "{log_path}" 2>&1'
        )
        shell(cmd)

        if params.adapter_args == "--detect_adapter_for_pe":
            with open(log_path) as fh:
                fh.seek(start)
                guesses = _fastp_unknown_adapter_guesses(fh.read())
            for read, seq in guesses:
                unit = _tandem_repeat_unit(seq)
                hint = (f" It is a short tandem repeat (unit {unit}), e.g. a telomere or "
                        f"satellite repeat." if unit else "")
                msg = (f"[GABi] WARNING: {wildcards.sample_id} Path{wildcards.idx}: fastp guessed "
                       f"the {read} adapter as {seq[:30]}..., which is not a known adapter. "
                       f"It may be genome sequence, now trimmed from every read that "
                       f"contains it.{hint} Check the fastp report, or set "
                       f"ADAPTERS_PE: illumina and re-trim this sample.")
                with open(log_path, "a") as fh:
                    fh.write(msg + "\n")
                print(msg)


# -------------------------------------------------------------------------------
#  QC RULES
# ===============================================================================

rule B06_fastqc_reads:
    """Run FastQC on subsampled reads"""
    input:
        reads = _get_fastqc_input
    output:
        html = os.path.join(
            _DATA_BASE, "{sample_id}", "{read_type}", "reports",
            "{filename}_fastqc.html"
        ),
        zip = os.path.join(
            _DATA_BASE, "{sample_id}", "{read_type}", "reports",
            "{filename}_fastqc.zip"
        )
    wildcard_constraints:
        filename = r"[^/]+"
    threads: cpu_func("reads_qc")
    resources:
        mem_mb = mem_func("reads_qc"),
        runtime = time_func("reads_qc")
    container: CONTAINERS["gabi_base"]
    log:
        os.path.join(_DATA_BASE, "{sample_id}", "{read_type}", "logs",
                     "B06_fastqc.{filename}.log")
    shell:
        r'''
        set -euo pipefail
        mkdir -p "$(dirname {output.html})" "$(dirname {log})"


        # TEMP DIRECTORY SETUP
        # -------------------------------------------------------------------
        WORK_DIR="$(gabi_get_workdir 5)"
        TEMP_DIR="$(mktemp -d "$WORK_DIR/GABi_fastqc_{wildcards.sample_id}_XXXXXX")"
        trap 'rm -rf "$TEMP_DIR"' EXIT
        cd "$TEMP_DIR"

        # fastqc is a java wrapper - set memory
        export _JAVA_OPTIONS="-Xmx$(({resources.mem_mb} * 90 / 100))m"

        # Subsample for FastQC (250k reads)
        echo "Subsampling for FastQC..." > {log}
        TEMP_SUB="{wildcards.filename}.fq"
        seqtk sample -s123 {input.reads} 250000 > $TEMP_SUB 2>> {log}

        # Run FastQC on subsample
        fastqc -t {threads} -o . $TEMP_SUB >> {log} 2>&1

        # Move outputs to correct location
        mv "{wildcards.filename}_fastqc.html" "{output.html}"
        mv "{wildcards.filename}_fastqc.zip" "{output.zip}"

        echo "[GABi] FastQC complete"
        '''


rule B07_nanoplot_long_reads:
    """Run NanoPlot on long reads"""
    input:
        reads = _get_nanoplot_input
    output:
        dir = directory(os.path.join(
            _DATA_BASE, "{sample_id}", "{read_type}", "reports",
            "Path{idx}_{base}{suffix}_nanoplot"
        ))
    wildcard_constraints:
        read_type = r"(hifi|ont|pacbio|nanopore)",
        suffix = r"(|_corrected|_filtered)",
        base = r"[^/]+?(?<!_corrected)(?<!_filtered)"
    threads: cpu_func("reads_qc")
    resources:
        mem_mb = mem_func("reads_qc"),
        runtime = time_func("reads_qc")
    container: CONTAINERS["gabi_base"]
    log:
        os.path.join(_DATA_BASE, "{sample_id}", "{read_type}", "logs",
                     "B07_nanoplot.Path{idx}.{base}{suffix}.log")
    shell:
        r'''
        set -euo pipefail
        mkdir -p "$(dirname {log})"
        exec > {log} 2>&1

        echo "[GABi] Running NanoPlot for {wildcards.sample_id} {wildcards.read_type}"

        # Create output directory
        mkdir -p {output.dir}

        # Run NanoPlot with subsampling for speed
        NanoPlot \
            --fastq {input.reads} \
            --outdir {output.dir} \
            --prefix {wildcards.sample_id}_{wildcards.read_type}_ \
            --threads {threads} \
            --downsample 100000 \
            --plots hex dot \
            --N50 \
            --loglength

        echo "[GABi] NanoPlot complete"
        '''


# -------------------------------------------------------------------------------
#  AGGREGATE QC
# ===============================================================================

def _get_qc_reports_for_multiqc(w):
    """
    Paths are dedup because the same report is now reached once per assembly.
    """
    reports = []
    seen = set()

    def want(path):
        if path not in seen:
            seen.add(path)
            reports.append(path)

    for rt, reads in _collect_sample_reads(w.sample_id).items():
        reports_dir = os.path.join(_DATA_BASE, w.sample_id, rt, "reports")

        for grp in _enumerate_groups(rt, reads):
            idx = grp["idx"]
            base = grp["base"]

            if grp["kind"] == "long":
                # FastQC for raw reads
                want(os.path.join(reports_dir, f"Path{idx}_{base}_fastqc.zip"))
                # NanoPlot
                want(os.path.join(reports_dir, f"Path{idx}_{base}_nanoplot"))

                # ULI check for HiFi
                if rt == "hifi":
                    want(os.path.join(reports_dir, f"Path{idx}_{base}_uli_mqc.yaml"))
                    # Filtering stats if enabled
                    if _as_bool(config.get("FILTER_HIFI", True)):
                        want(os.path.join(reports_dir, f"Path{idx}_{base}_bbduk.stats"))

                # Corrected ONT NanoPlot
                if rt == "ont" and _as_bool(config.get("CORRECT_ONT", False)):
                    want(os.path.join(reports_dir, f"Path{idx}_{base}_corrected_nanoplot"))

            elif grp["kind"] == "pe":
                # Per-mate FastQC in both cases. With TRIM_PE on, this used to
                # ask for the fastp JSON alone, so the default run had no
                # per-mate FastQC at all - while _get_fastqc_input() carried a
                # whole branch for reading the raw R1/R2 that nothing ever
                # requested. B06 subsamples to 250k reads, so the two extra
                # jobs are cheap, and the report then shows pre-trim per-mate
                # quality beside the trimming summary.
                want(os.path.join(reports_dir, f"Path{idx}_{base}_R1_fastqc.zip"))
                want(os.path.join(reports_dir, f"Path{idx}_{base}_R2_fastqc.zip"))
                if _as_bool(config.get("TRIM_PE", True)):
                    # Fastp reports
                    want(os.path.join(reports_dir, f"Path{idx}_{base}_fastp.json"))

    return reports


rule B08_aggregate_read_qc:
    """Aggregate all read QC reports with MultiQC"""
    input:
        reports = _get_qc_reports_for_multiqc
    output:
        html = os.path.join(_DATA_BASE, "{sample_id}", "multiqc_report.html"),
        # MultiQC names its data directory after the report: multiqc_report.html
        # -> multiqc_report_data. The rule used to declare 'multiqc_data', pass
        # --no-data-dir, and then mkdir the empty directory to satisfy the
        # declaration.
        data = directory(os.path.join(_DATA_BASE, "{sample_id}", "multiqc_report_data"))
    threads: 1
    resources:
        mem_mb = mem_func("light_task"),
        runtime = time_func("light_task")
    container: CONTAINERS["gabi_base"]
    log:
        os.path.join(_DATA_BASE, "{sample_id}", "logs", "B08_multiqc.log")
    shell:
        r'''
        set -euo pipefail
        mkdir -p "$(dirname {output.html})" "$(dirname {log})"
        exec > {log} 2>&1

        echo "[GABi] Running MultiQC for sample {wildcards.sample_id}"

        # Create file list for MultiQC
        FILELIST=$(mktemp)
        trap 'rm -f "$FILELIST"' EXIT

        for f in {input.reports}; do
            if [[ -e "$f" ]]; then
                echo "$f" >> "$FILELIST"
            fi
        done

        if [[ ! -s "$FILELIST" ]]; then
            echo "[GABi] ERROR: no QC reports found for {wildcards.sample_id}" >&2
            exit 1
        fi

        # Run MultiQC. No '|| true': a MultiQC that fails outright used to
        # still mark the job complete.
        multiqc \
            --file-list "$FILELIST" \
            --outdir "$(dirname {output.html})" \
            --filename "$(basename {output.html})" \
            --force \
            --data-dir 2>&1

        echo "[GABi] MultiQC complete"
        '''
