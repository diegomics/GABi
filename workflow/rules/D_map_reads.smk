# ===============================================================================
# GABi - Mapping Rules
# by Diego De Panis, 2026
# note: AI tools may have been used to improve, clean and/or comment this version of the code
# ===============================================================================

# consistent idx/base separator for BAM filenames
_SEP = "_"

# -------------------------------------------------------------------------------
#  HELPER FUNCTIONS
# -------------------------------------------------------------------------------

def _norm_rt(rt):
    return normalize_read_type(rt)

def _maybe_temp(path, condition):
    """Wrap path in Snankemake's directive temp() if condition is True to automatically delete intermediate files"""
    return temp(path) if condition else path

def _idx_dir(sp, asm):
    return os.path.join(config["OUT_FOLDER"], "GABi_results", sp, asm, "idx")

def _bams_dir(sp, asm, sid):
    return os.path.join(config["OUT_FOLDER"], "GABi_results", sp, asm, "samples", sid, "BAMs")

# def _ref_fa(sp, asm):
#    return samples_config["sp_name"][sp]["asm_id"][asm]["asm_file"]

def _platform_from_rt(rt):
    rt = _norm_rt(rt)
    if rt in {"illumina","pe","sr","10x","tenx","10xgenomics","hic","hi-c","linked","linkedreads"}:
        return "ILLUMINA"
    if rt in {"hifi","pacbio","pacbiohifi","pb_hifi","pbhifi"}:
        return "PACBIO"
    if rt in {"ont","nanopore","oxfordnanopore"}:
        return "ONT"
    return "ILLUMINA"

def _rg_string(w, idx, base):
    # read group: ID + SM + PL (+ LB, PU for disambiguation)
    sm = w.sample_id
    rt = _norm_rt(w.read_type)
    pl = _platform_from_rt(rt)
    lb = f"{sm}.{rt}"
    pu = f"{sm}.P{idx}"
    rgid = f"{sm}.{rt}.P{idx}.{base}"
    return f'@RG\\tID:{rgid}\\tSM:{sm}\\tPL:{pl}\\tLB:{lb}\\tPU:{pu}'

def _minimap_preset(rt):
    rt = _norm_rt(rt)
    if rt in {"hifi","pacbio","pacbiohifi","pb_hifi","pbhifi"}:
        return "map-hifi"
    if rt in {"ont","nanopore","oxfordnanopore"}:
        return "map-ont"
    return "map-hifi"

def _reads_for_group(w):
    """
    Return (kind, paths) for mapping inputs of the requested group,
    using the standardized files in reads/.
    """
    sp, asm, sid = w.species, w.assembly, w.sample_id
    # Pass the read_type from wildcards to get the correct data
    rt, trim, read_dict = _get_sample_node(sp, asm, sid, w.read_type)
    if rt is None:
        raise ValueError(f"No valid read data for {sp}/{asm}/{sid}/{w.read_type}")

    # Verify the read type matches
    if rt != normalize_read_type(w.read_type):
        raise ValueError(f"read_type mismatch: got {rt}, expected {normalize_read_type(w.read_type)}")

    for grp in _enumerate_groups(rt, read_dict):
        if str(grp["idx"]) == str(w.idx) and grp["base"] == w.base:
            outs = _outputs_for_group(sp, asm, sid, grp, trim)
            if grp["kind"] == "long":
                return "long", (outs[0],)
            else:
                return "pe", (outs[0], outs[1])
    raise ValueError(f"Could not find group Path{w.idx}_{w.base} for {sp}/{asm}/{sid}/{w.read_type}")

def _minibwa_prefix(sp, asm):
    return os.path.join(_idx_dir(sp, asm), asm)  # prefix

def _minimap_index_path(sp, asm):
    return os.path.join(_idx_dir(sp, asm), f"{asm}.mmi")

def _final_group_bam(sp, asm, sid, grp):
    rt = grp["rt"]; idx = grp["idx"]; base = grp["base"]
    bdir = _bams_dir(sp, asm, sid)
    stem = f"{rt}_Path{idx}{_SEP}{base}"
    if grp["kind"] == "pe":
        return os.path.join(
            bdir,
            f"{stem}.markdup.bam" if _as_bool(config.get("DEDUP_PE", True)) else f"{stem}.sorted.bam",
        )
    else:
        return os.path.join(
            bdir,
            f"{stem}.markdup.bam" if _as_bool(config.get("DEDUP_LONG", True)) else f"{stem}.sorted.bam",
        )

def _qc_dir(sp, asm, sid):
    return os.path.join(_bams_dir(sp, asm, sid), "qc")

def _reads_for_pe_group(w):
    kind, paths = _reads_for_group(w)
    if kind != "pe":
        raise ValueError(f"Requested PE group but found {kind} for {w.species}/{w.assembly}/{w.sample_id}")
    return paths  # (r1, r2)

def _reads_for_long_group(w):
    kind, paths = _reads_for_group(w)
    if kind != "long":
        raise ValueError(f"Requested long group but found {kind} for {w.species}/{w.assembly}/{w.sample_id}")
    return paths  # (r1,)

def _qc_summaries_for_assembly(w):
    sp, asm = w.species, w.assembly
    ins = []
    for _sp, _asm, sid in all_samples_iter():
        if _sp != sp or _asm != asm:
            continue

        # Get all read types for this sample
        try:
            sample_data = samples_config["sp_name"][_sp]["asm_id"][_asm]["sample_id"][sid]
            if not sample_data.get("read_type"):
                continue

            # Process each read type
            for rt_key in sample_data["read_type"]:
                if rt_key == "None" or rt_key is None:
                    continue

                rt_norm = normalize_read_type(rt_key)
                _, _, reads = _get_sample_node(_sp, _asm, sid, rt_key)

                if reads:  # Only if there are actual read files
                    qcdir = os.path.join(config["OUT_FOLDER"], "GABi_results", sp, asm,
                                       "samples", sid, "BAMs", "qc")
                    ins.append(os.path.join(qcdir, f"{sid}.{rt_norm}.summary.txt"))
        except (KeyError, TypeError, AttributeError):
            continue

    return ins


# -------------------------------------------------------------------------------
#  INDEXING
# ===============================================================================
# one canonical reference FASTA in idx/: symlink if uncompressed, else decompress once
rule D00_ref_fa_canonical:
    input:
        ref=lambda w: get_assembly_path(w.species, w.assembly)
    output:
        fa=os.path.join(config["OUT_FOLDER"], "GABi_results", "{species}", "{assembly}", "idx", "{assembly}.fa")
    threads: cpu_func("split_intervals")
    resources:
        mem_mb=mem_func("split_intervals"),
        runtime=time_func("split_intervals")
    container: CONTAINERS["gabi_base"]
    shell:
        r'''
        set -euo pipefail
        mkdir -p "$(dirname {output.fa})"
        if [[ "{input.ref}" == *.gz ]]; then
            gunzip -c "{input.ref}" > "{output.fa}"
        else
            ln -sf "$(realpath "{input.ref}")" "{output.fa}"
        fi
        '''

rule D01_minimap2_index:
    """Build minimap2 index once per assembly."""
    input:
        fa=rules.D00_ref_fa_canonical.output.fa
    output:
        mmi=os.path.join(config["OUT_FOLDER"], "GABi_results", "{species}", "{assembly}", "idx", "{assembly}.mmi")
    threads: cpu_func("minimap2")
    resources:
        mem_mb=mem_func("minimap2"),
        runtime=time_func("minimap2")
    container: CONTAINERS["gabi_base"]
    log:
        os.path.join(config["OUT_FOLDER"], "GABi_results", "{species}", "{assembly}", "logs", "D01_minimap2_index.log")
    shell:
        r'''
        set -euo pipefail
        mkdir -p $(dirname {log})

        exec > {log} 2>&1

        echo "[GABi] Creating minimap2 index for {wildcards.assembly}"


        # TEMP DIRECTORY SETUP
        # -------------------------------------------------------------------
        WORK_DIR="$(gabi_get_workdir 75)"
        TEMP_DIR="$(mktemp -d "$WORK_DIR/GABi_mm2index_{wildcards.species}_{wildcards.assembly}_XXXXXX")"
        trap 'rm -rf "$TEMP_DIR"' EXIT
        cd "$TEMP_DIR"

        echo "[GABi] Copying assembly to temp directory..."
        cp "{input.fa}" local.fa

        echo "[GABi] Running minimap2 -d..."
        minimap2 -d local.mmi local.fa

        echo "[GABi] Copying index back to output..."
        mkdir -p "$(dirname {output.mmi})"
        cp local.mmi "{output.mmi}"

        echo "[GABi] Done"

        '''

rule D01_minibwa_index:
    """Create the minibwa index once per assembly (prefix-based).
    """
    input:
        fa=rules.D00_ref_fa_canonical.output.fa
    output:
        flag=os.path.join(config["OUT_FOLDER"], "GABi_results", "{species}", "{assembly}", "idx", "minibwa.done")
    params:
        prefix=lambda w: os.path.join(_idx_dir(w.species, w.assembly), w.assembly),
        idx_dir=lambda w: _idx_dir(w.species, w.assembly)
    threads: cpu_func("minibwa_index")
    resources:
        mem_mb=mem_func("minibwa_index"),
        runtime=time_func("minibwa_index")
    container: CONTAINERS["gabi_base"]
    log:
        os.path.join(config["OUT_FOLDER"], "GABi_results", "{species}", "{assembly}", "logs", "D01_minibwa_index.log")
    shell:
        r'''
        set -euo pipefail
        mkdir -p $(dirname {log})

        exec > {log} 2>&1

        echo "[GABi] Creating minibwa index for {wildcards.assembly}"

        mkdir -p "{params.idx_dir}"

        # Built directly at the final prefix: minibwa takes an explicit
        # out.prefix, so there is nothing to rename afterwards.
        minibwa index -t {threads} "{input.fa}" "{params.prefix}"

        echo "[GABi] Index files:"
        ls -la "{params.prefix}".* || true

        for ext in l2b mbw; do
            if [[ ! -s "{params.prefix}.$ext" ]]; then
                echo "[GABi] ERROR: expected index file {params.prefix}.$ext is missing or empty." >&2
                echo "[GABi]        minibwa produced:" >&2
                ls -la "$(dirname {params.prefix})" >&2
                exit 1
            fi
        done

        touch "{output.flag}"

        echo "[GABi] Done"
        '''


# -------------------------------------------------------------------------------
#  LONG READS MAPPING
# ===============================================================================

rule D02_map_long_reads:
    """minimap2 mapping of one long-read group -> sorted BAM + index + flagstat."""
    input:
        r1=lambda w: _reads_for_long_group(w)[0],
        idx=rules.D01_minimap2_index.output.mmi
    output:
        bam=_maybe_temp(
            os.path.join(
                config["OUT_FOLDER"], "GABi_results", "{species}", "{assembly}",
                "samples", "{sample_id}", "BAMs",
                "{read_type}_Path{idx}" + _SEP + "{base}.sorted.bam"
            ),
            _as_bool(config.get("DEDUP_LONG", True))
        ),
        bai=_maybe_temp(
            os.path.join(
                config["OUT_FOLDER"], "GABi_results", "{species}", "{assembly}",
                "samples", "{sample_id}", "BAMs",
                "{read_type}_Path{idx}" + _SEP + "{base}.sorted.bam.bai"
            ),
            _as_bool(config.get("DEDUP_LONG", True))
        )
    wildcard_constraints:
        read_type=r"(hifi|ont)",
        idx=r"\d+"
    threads: cpu_func("minimap2")
    resources:
        mem_mb=mem_func("minimap2"),
        runtime=time_func("minimap2")
    container: CONTAINERS["gabi_base"]
    log:
        os.path.join(
            config["OUT_FOLDER"], "GABi_results", "{species}", "{assembly}", "samples", "{sample_id}",
            "logs", "D02_map_long_reads.{read_type}.Path{idx}" + _SEP + "{base}.log"
        )
    run:
        r1 = input.r1
        preset = _minimap_preset(wildcards.read_type)
        rg = _rg_string(wildcards, wildcards.idx, wildcards.base)
        ref_mmi = input.idx
        os.makedirs(os.path.dirname(output.bam), exist_ok=True)
        os.makedirs(os.path.dirname(str(log)), exist_ok=True)
        # Streaming pipe - no temp dir needed
        shell(
            'set -euo pipefail\n'
            'minimap2 -t {t} -a -x {preset} -R "{rg}" "{mmi}" "{r1}" | '
            'samtools sort -@ {t} -o "{bam}"\n'
            'samtools index -@ {t} "{bam}"\n'
            'samtools flagstat -@ {t} "{bam}" >> "{log}"\n'
            .format(t=threads, preset=preset, rg=rg, mmi=ref_mmi, r1=r1, bam=output.bam, log=str(log))
        )

rule D03_markdup_long_reads:
    """Optional: mark duplicates for long reads group (enabled via DEDUP_LONG=true)."""
    input:
        bam=rules.D02_map_long_reads.output.bam
    output:
        bam=os.path.join(
            config["OUT_FOLDER"], "GABi_results", "{species}", "{assembly}", "samples", "{sample_id}", "BAMs",
            "{read_type}_Path{idx}" + _SEP + "{base}.markdup.bam",
        ),
        bai=os.path.join(
            config["OUT_FOLDER"], "GABi_results", "{species}", "{assembly}", "samples", "{sample_id}", "BAMs",
            "{read_type}_Path{idx}" + _SEP + "{base}.markdup.bam.bai",
        ),
        metrics=os.path.join(
            config["OUT_FOLDER"], "GABi_results", "{species}", "{assembly}", "samples", "{sample_id}", "BAMs",
            "{read_type}_Path{idx}" + _SEP + "{base}.markdup_metrics.txt",
        )
    wildcard_constraints:
        read_type=r"(hifi|ont)",
        idx=r"\d+"
    threads: cpu_func("markdup")
    resources:
        mem_mb=mem_func("markdup"),
        runtime=time_func("markdup")
    container: CONTAINERS["gabi_base"]
    log:
        os.path.join(
            config["OUT_FOLDER"], "GABi_results", "{species}", "{assembly}", "samples", "{sample_id}",
            "logs", "D03_markdup_long_reads.{read_type}.Path{idx}" + _SEP + "{base}.log"
        )
    shell:
        r'''
        set -euo pipefail
        mkdir -p "$(dirname {output.bam})" "$(dirname {log})"

        exec > "{log}" 2>&1

        echo "[GABi] Starting duplicate marking for {input.bam}"


        # TEMP DIRECTORY SETUP
        # -------------------------------------------------------------------
        WORK_DIR="$(gabi_get_workdir 100)"
        TEMP_DIR="$(mktemp -d "$WORK_DIR/GABi_markduplong_{wildcards.species}_{wildcards.sample_id}_XXXXXX")"
        trap 'rm -rf "$TEMP_DIR"' EXIT

        echo "[GABi] Using temp directory: $TEMP_DIR"

        echo "[GABi] Running sambamba markdup"
        sambamba markdup -p -t {threads} --tmpdir "$TEMP_DIR" \
            "{input.bam}" "{output.bam}" 2>&1 | tee "{output.metrics}"

        echo "[GABi] Indexing BAM"
        samtools index -@ {threads} "{output.bam}"

        echo "[GABi] Done"

        '''


# -------------------------------------------------------------------------------
#  SHORT READS MAPPING
# ===============================================================================

rule D02_map_pe_reads:
    """minibwa mapping of one PE group -> sorted BAM + index + flagstat.
    """
    input:
        r1=lambda w: _reads_for_pe_group(w)[0],
        r2=lambda w: _reads_for_pe_group(w)[1],
        idx=rules.D01_minibwa_index.output.flag
    output:
        bam=_maybe_temp(
            os.path.join(
                config["OUT_FOLDER"], "GABi_results", "{species}", "{assembly}",
                "samples", "{sample_id}", "BAMs",
                "{read_type}_Path{idx}" + _SEP + "{base}.sorted.bam"
            ),
            _as_bool(config.get("DEDUP_PE", True))
        ),
        bai=_maybe_temp(
            os.path.join(
                config["OUT_FOLDER"], "GABi_results", "{species}", "{assembly}",
                "samples", "{sample_id}", "BAMs",
                "{read_type}_Path{idx}" + _SEP + "{base}.sorted.bam.bai"
            ),
            _as_bool(config.get("DEDUP_PE", True))
        )
    wildcard_constraints:
        read_type=r"(illumina|10x)",
        idx=r"\d+"
    threads: cpu_func("minibwa_map")
    resources:
        mem_mb=mem_func("minibwa_map"),
        runtime=time_func("minibwa_map")
    container: CONTAINERS["gabi_base"]
    log:
        os.path.join(
            config["OUT_FOLDER"], "GABi_results", "{species}", "{assembly}", "samples", "{sample_id}",
            "logs", "D02_map_pe_reads.{read_type}.Path{idx}" + _SEP + "{base}.log"
        )
    run:
        r1, r2 = input.r1, input.r2
        ref_prefix = _minibwa_prefix(wildcards.species, wildcards.assembly)
        rg = _rg_string(wildcards, wildcards.idx, wildcards.base)
        os.makedirs(os.path.dirname(output.bam), exist_ok=True)
        os.makedirs(os.path.dirname(str(log)), exist_ok=True)
        # Streaming pipe - no temp dir needed
        shell(
            'set -euo pipefail\n'
            '(\n'
            '  minibwa map -x sr -t {t} -R "{rg}" "{ref}" "{r1}" "{r2}" | \n'
            '  samtools sort -@ {t} -o "{bam}"\n'
            ') 2>> "{log}"\n'
            'samtools index -@ {t} "{bam}" 2>> "{log}"\n'
            'samtools flagstat -@ {t} "{bam}" >> "{log}"\n'
            .format(t=threads, rg=rg, ref=ref_prefix, r1=r1, r2=r2,
                    bam=output.bam, log=str(log))
        )

rule D03_markdup_pe_reads:
    """Mark duplicates for PE group (enabled by DEDUP_PE)."""
    input:
        bam=rules.D02_map_pe_reads.output.bam
    output:
        bam=os.path.join(
            config["OUT_FOLDER"], "GABi_results", "{species}", "{assembly}", "samples", "{sample_id}", "BAMs",
            "{read_type}_Path{idx}" + _SEP + "{base}.markdup.bam",
        ),
        bai=os.path.join(
            config["OUT_FOLDER"], "GABi_results", "{species}", "{assembly}", "samples", "{sample_id}", "BAMs",
            "{read_type}_Path{idx}" + _SEP + "{base}.markdup.bam.bai",
        ),
        metrics=os.path.join(
            config["OUT_FOLDER"], "GABi_results", "{species}", "{assembly}", "samples", "{sample_id}", "BAMs",
            "{read_type}_Path{idx}" + _SEP + "{base}.markdup_metrics.txt",
        )
    wildcard_constraints:
        read_type=r"(illumina|10x)",
        idx=r"\d+"
    threads: cpu_func("markdup")
    resources:
        mem_mb=mem_func("markdup"),
        runtime=time_func("markdup")
    container: CONTAINERS["gabi_base"]
    log:
        os.path.join(
            config["OUT_FOLDER"], "GABi_results", "{species}", "{assembly}", "samples", "{sample_id}",
            "logs", "D03_markdup_pe_reads.{read_type}.Path{idx}" + _SEP + "{base}.log"
        )
    shell:
        r'''
        set -euo pipefail
        mkdir -p "$(dirname {output.bam})" "$(dirname {log})"

        exec > "{log}" 2>&1

        echo "[GABi] Starting duplicate marking for {input.bam}"


        # TEMP DIRECTORY SETUP
        # -------------------------------------------------------------------
        WORK_DIR="$(gabi_get_workdir 100)"
        TEMP_DIR="$(mktemp -d "$WORK_DIR/GABi_markdup_{wildcards.species}_{wildcards.sample_id}_XXXXXX")"
        trap 'rm -rf "$TEMP_DIR"' EXIT

        echo "[GABi] Using temp directory: $TEMP_DIR"

        echo "[GABi] Running sambamba markdup"
        sambamba markdup -p -t {threads} --tmpdir "$TEMP_DIR" \
            "{input.bam}" "{output.bam}" 2>&1 | tee "{output.metrics}"

        echo "[GABi] Indexing BAM"
        samtools index -@ {threads} "{output.bam}"

        echo "[GABi] Done"

        '''



# -------------------------------------------------------------------------------
#  PRE-ALIGNED INPUT
# ===============================================================================
#  Alternative entry point: the data table gives a .BAM in read_files instead of
#  FASTQ. read_type still declares the technology, exactly as for reads - the
#  header is read only to CHECK that declaration, never to supply it.
#
#  The BAM is registered as {sid}.{rt}.imported.bam and handed to D04, whose run
#  block already symlinks rather than merges when given a single file. From
#  {sid}.{rt}.merged.bam onwards nothing knows or cares how the BAM was made:
#  D05 QC, D06 sex inference, D07 selection, D08 linking and all of E and F are
#  unchanged.
#
#  The symlink is absolute because the source usually lives outside GABi_results.
#  As with reads, the path must be on a partition bound into the container.

def _input_bam_for(w):
    """The user's BAM for this sample/read_type, as a rule input."""
    path = bam_entry_path(w.species, w.assembly, w.sample_id, w.read_type)
    return [path] if path else []


rule D02_import_bam:
    """Validate a pre-aligned BAM and register it for this sample/read_type."""
    input:
        bam=_input_bam_for,
        fai=os.path.join(
            config["OUT_FOLDER"], "GABi_results", "{species}", "{assembly}",
            "idx", "{assembly}.fa.fai"
        ),
    output:
        bam=os.path.join(
            config["OUT_FOLDER"], "GABi_results", "{species}", "{assembly}",
            "samples", "{sample_id}", "BAMs", "{sample_id}.{read_type}.imported.bam"
        ),
        report=os.path.join(
            config["OUT_FOLDER"], "GABi_results", "{species}", "{assembly}",
            "samples", "{sample_id}", "BAMs", "qc",
            "{sample_id}.{read_type}.imported_bam.md"
        ),
    params:
        script=str(scripts_dir / "validate_input_bam.py"),
    threads: cpu_func("light_task")
    resources:
        mem_mb=mem_func("light_task"),
        runtime=time_func("light_task")
    container: CONTAINERS["gabi_base"]
    log:
        os.path.join(
            config["OUT_FOLDER"], "GABi_results", "{species}", "{assembly}",
            "samples", "{sample_id}", "logs",
            "D02_import_bam.{sample_id}.{read_type}.log"
        )
    shell:
        r'''
        set -euo pipefail
        mkdir -p "$(dirname {output.bam})" "$(dirname {output.report})" "$(dirname {log})"
        exec > "{log}" 2>&1

        echo "[GABi] Pre-aligned input for {wildcards.sample_id} ({wildcards.read_type})"
        echo "[GABi] Source: {input.bam}"

        # Header only - no pysam needed, and no pass over the alignments.
        #
        # mktemp, not a bare filename: this rule does not cd anywhere, so a
        # relative path lands in Snakemake's working directory and EVERY sample's
        # import job would write to the same file concurrently
        HEADER="$(mktemp "${{TMPDIR:-/tmp}}/GABi_bamhdr_{wildcards.sample_id}_{wildcards.read_type}_XXXXXX")"
        trap 'rm -f "$HEADER"' EXIT
        samtools view -H "{input.bam}" > "$HEADER"

        python {params.script:q} \
            --header "$HEADER" \
            --fai "{input.fai}" \
            --bam "{input.bam}" \
            --sample-id "{wildcards.sample_id}" \
            --read-type "{wildcards.read_type}" \
            --report "{output.report}"

        # Absolute: the source normally sits outside GABi_results, so a relative
        # link would depend on how the two paths are bound into the container.
        ln -sfn "$(realpath "{input.bam}")" "{output.bam}"
        [ -e "{output.bam}" ] || {{ echo "[GABi] ERROR: symlink does not resolve" >&2; exit 1; }}

        echo "[GABi] Registered as $(basename "{output.bam}")"
        '''


# ----------------------------- merge per sample & read type -----------------------------

def _imported_bam(sp, asm, sid, rt):
    """Path of the validated symlink D02_import_bam produces."""
    return os.path.join(_bams_dir(sp, asm, sid), f"{sid}.{rt}.imported.bam")


def _sample_group_bams_by_type(w):
    """What D04 merges for this sample&read_type.

    For read input that is the per-group BAMs from D02/D03. For pre-aligned
    input it is the single validated symlink from D02_import_bam - D04's run
    block already symlinks rather than merges when handed one file, so the
    import path needs no special case there.
    """
    sp, asm, sid = w.species, w.assembly, w.sample_id
    rt_req = normalize_read_type(w.read_type)

    # Get data for the specific read type requested
    rt, trim, reads = _get_sample_node(sp, asm, sid, w.read_type)
    if rt is None or rt != rt_req:
        return []

    if is_bam_entry(sp, asm, sid, w.read_type):
        return [_imported_bam(sp, asm, sid, rt)]

    outs = []
    for grp in _enumerate_groups(rt, reads):
        outs.append(_final_group_bam(sp, asm, sid, grp))
    return outs

rule D04_merge_sample_readtype_bams:
    input:
        _sample_group_bams_by_type
    output:
        bam=os.path.join(
            config["OUT_FOLDER"], "GABi_results", "{species}", "{assembly}", "samples", "{sample_id}", "BAMs",
            "{sample_id}.{read_type}.merged.bam",
        ),
        bai=os.path.join(
            config["OUT_FOLDER"], "GABi_results", "{species}", "{assembly}", "samples", "{sample_id}", "BAMs",
            "{sample_id}.{read_type}.merged.bam.bai",
        ),
        stat=os.path.join(
            config["OUT_FOLDER"], "GABi_results", "{species}", "{assembly}", "samples", "{sample_id}", "BAMs",
            "{sample_id}.{read_type}.merged.flagstat.txt",
        )
    threads: cpu_func("minibwa_map")
    resources:
        mem_mb=mem_func("minibwa_map"),
        runtime=time_func("minibwa_map")
    container: CONTAINERS["gabi_base"]
    log:
        os.path.join(
            config["OUT_FOLDER"], "GABi_results", "{species}", "{assembly}", "samples", "{sample_id}",
            "logs", "D04_merge_bams.{sample_id}.{read_type}.log"
        )
    run:
        os.makedirs(os.path.dirname(output.bam), exist_ok=True)
        if len(input) == 0:
            raise ValueError(f"No BAMs to merge for {wildcards.sample_id}/{wildcards.read_type}")
        if len(input) == 1:
            import os as _os
            rel = _os.path.relpath(input[0], _os.path.dirname(output.bam))
            shell('ln -sf "{src}" "{dst}"'.format(src=rel, dst=output.bam))
        else:
            shell('samtools merge -@ {t} -o "{out}" {ins}'.format(
                t=threads, out=output.bam, ins=" ".join(f'"{x}"' for x in input)
            ))
        shell(
            'set -euo pipefail\n'
            'samtools index -@ {t} "{bam}"\n'
            'samtools flagstat -@ {t} "{bam}" > "{stat}"\n'
            .format(t=threads, bam=output.bam, stat=output.stat)
        )

# ----------------------------- QC of merged BAM -----------------------------
rule D05_qc_merged_bam:
    """QC of merged BAM: flagstat, idxstats, coverage, concise summary (incl. dispersion & outliers)."""
    input:
        bam=rules.D04_merge_sample_readtype_bams.output.bam,
        bai=rules.D04_merge_sample_readtype_bams.output.bai
    output:
        flagstat=os.path.join(
            config["OUT_FOLDER"], "GABi_results", "{species}", "{assembly}",
            "samples", "{sample_id}", "BAMs", "qc",
            "{sample_id}.{read_type}.flagstat.txt",
        ),
        idxstats=os.path.join(
            config["OUT_FOLDER"], "GABi_results", "{species}", "{assembly}",
            "samples", "{sample_id}", "BAMs", "qc",
            "{sample_id}.{read_type}.idxstats.tsv",
        ),
        coverage=os.path.join(
            config["OUT_FOLDER"], "GABi_results", "{species}", "{assembly}",
            "samples", "{sample_id}", "BAMs", "qc",
            "{sample_id}.{read_type}.coverage.tsv",
        ),
        summary=os.path.join(
            config["OUT_FOLDER"], "GABi_results", "{species}", "{assembly}",
            "samples", "{sample_id}", "BAMs", "qc",
            "{sample_id}.{read_type}.summary.txt",
        ),
        summary_md=os.path.join(
            config["OUT_FOLDER"], "GABi_results", "{species}", "{assembly}",
            "samples", "{sample_id}", "BAMs", "qc",
            "{sample_id}.{read_type}.summary.md",
        )
    params:
        qc_bam_summary = str(scripts_dir / "qc_bam_summary.py"),
        rt=lambda w: _norm_rt(w.read_type),   # canonical read type for summary
        disp_pct=0.95,                        # for Depth dispersion (95%)
        hi_mult=5.0                           # High-depth outliers threshold: depth >= hi_mult * WMD
    threads: cpu_func("bam_qc")
    resources:
        mem_mb=mem_func("bam_qc"),
        runtime=time_func("bam_qc")
    container: CONTAINERS["gabi_base"]
    log:
        os.path.join(
            config["OUT_FOLDER"], "GABi_results", "{species}", "{assembly}",
            "samples", "{sample_id}", "logs",
            "D05_qc_merged_bam.{sample_id}.{read_type}.log"
        )
    shell:
        r'''
        set -euo pipefail
        mkdir -p "$(dirname {output.flagstat})" "$(dirname {log})"

        exec > "{log}" 2>&1

        echo "[GABi] QC for sample={wildcards.sample_id} read_type={params.rt}"

        # samtools outputs
        samtools flagstat -@ {threads} "{input.bam}" > "{output.flagstat}"
        samtools idxstats "{input.bam}" > "{output.idxstats}"
        samtools coverage "{input.bam}" > "{output.coverage}"

        # Primary mapped (exclude secondary+supplementary+QCfail+unmapped) = -F 2308
        PRIMARY=$(samtools view -@ {threads} -c -F 2308 "{input.bam}")
        # Primary mapped, excluding duplicates too = -F 3332
        PRIMARY_NODUP=$(samtools view -@ {threads} -c -F 3332 "{input.bam}")

        # Compute summary & dispersion metrics
        python {params.qc_bam_summary:q} \
            --sample "{wildcards.sample_id}" \
            --read-type "{params.rt}" \
            --coverage "{output.coverage}" \
            --flagstat "{output.flagstat}" \
            --primary "$PRIMARY" \
            --primary-nodup "$PRIMARY_NODUP" \
            --disp-pct {params.disp_pct} \
            --hi-mult {params.hi_mult} \
            --out "{output.summary}" \
            --md-out "{output.summary_md}"

        echo "[GABi] QC complete"
        '''


# ===============================================================================
#  D06_INFER_SEX_PLOIDY — coverage-based ploidy inference for sex chromosomes
# ===============================================================================
#  Runs per sample. Reads the priority tech's coverage.tsv (from D05) and
#  classifies each declared sex contig as haploid / diploid / ambiguous.
#  Output drives DeepVariant's --haploid_contigs flag in E10 (wired from
#  E's DV input function).
#
#  Inference happens regardless of sample_sex — the rule always produces
#  a TSV documenting what coverage says about each declared sex contig.
#  The decision of how to *use* that TSV (haploid flagging, ignoring, etc.)
#  lives in E's DV input function.


# --- Helpers: sex inference reads the SELECTED technology's coverage ---
# Which technology that is only becomes known once D07 has run, so the rule
# depends on every callable technology's coverage.tsv and resolves the right one
# from tech_selection.json inside the shell. The files are a few hundred bytes
# each and D05 produces them for every technology regardless, so depending on
# all of them costs nothing and buys a correct dependency edge.
#
# This has to agree with the technology variant calling actually runs on: if sex
# inference read one technology's coverage while DeepVariant called on another,
# the haploid-contig flags would describe the wrong BAM!

_TECH_SELECTION_TMPL = os.path.join(
    config["OUT_FOLDER"], "GABi_results", "{species}", "{assembly}",
    "samples", "{sample_id}", "BAMs", "tech_selection.json"
)


def _d_declared_haploid(w):
    """Haploid contigs implied by a DECLARED sample_sex_karyo, else [].

    Used only for the cross-check below - the calling rules resolve ploidy
    through haploid_contigs_for_sample() in E_variant_call.smk.
    """
    contigs, source = resolve_haploid_contigs(
        samples_config, w.species, w.assembly, w.sample_id
    )
    return contigs if source == "declared" else []


def _d_all_tech_coverages(w):
    """coverage.tsv for every callable technology of this sample."""
    qc = _qc_dir(w.species, w.assembly, w.sample_id)
    return [
        os.path.join(qc, f"{w.sample_id}.{t}.coverage.tsv")
        for t in callable_techs(w.species, w.assembly, w.sample_id)
    ]


rule D06_infer_sex_ploidy:
    """
    Infer per-contig ploidy of sex chromosomes from coverage depth.

    Reads the SELECTED technology's coverage.tsv (per tech_selection.json) and
    classifies each declared asm_sex_chr contig as haploid / diploid /
    ambiguous. Always runs when its output is requested by a downstream rule.
    If sex_chr is empty for the assembly, produces a header-only TSV — keeps
    the DAG consistent.
    """
    input:
        selection=_TECH_SELECTION_TMPL,
        coverages=_d_all_tech_coverages,
    output:
        tsv=os.path.join(
            config["OUT_FOLDER"], "GABi_results", "{species}", "{assembly}",
            "samples", "{sample_id}", "BAMs", "qc",
            "{sample_id}.sex_inference.tsv",
        ),
    params:
        # sex_chr passed as comma-separated string; may be empty.
        sex_chr=lambda w: ",".join(
            get_sex_chr(samples_config, w.species, w.assembly)
        ),
        tiny_contig_bp=int(config.get("TINY_CONTIG_BP", 1_000_000)),
        script=str(scripts_dir / "infer_sex_ploidy.py"),
        qc_dir=lambda w: _qc_dir(w.species, w.assembly, w.sample_id),
        # For the declared-vs-inferred cross-check below.
        sample_sex=lambda w: get_sample_sex(
            samples_config, w.species, w.assembly, w.sample_id
        ),
        declared_haploid=lambda w: ",".join(_d_declared_haploid(w)),
    threads: 1
    resources:
        mem_mb=mem_func("split_intervals"),
        runtime=time_func("split_intervals"),
    container: CONTAINERS["gabi_base"]
    log:
        os.path.join(
            config["OUT_FOLDER"], "GABi_results", "{species}", "{assembly}",
            "samples", "{sample_id}", "logs",
            "D06_infer_sex_ploidy.{sample_id}.log",
        )
    shell:
        r'''
        set -euo pipefail
        mkdir -p "$(dirname {output.tsv})" "$(dirname {log})"
        exec > "{log}" 2>&1

        echo "[GABi] Sex-chromosome ploidy inference for {wildcards.sample_id}"

        # Resolve the selected technology at run time. Every candidate
        # coverage.tsv is a declared input, so whichever one this picks is
        # already present and up to date.
        TECH=$(python -c "import json,sys; print(json.load(open(sys.argv[1]))['chosen'])" "{input.selection}")
        COVERAGE="{params.qc_dir}/{wildcards.sample_id}.${{TECH}}.coverage.tsv"

        if [ ! -s "$COVERAGE" ]; then
            echo "[GABi] ERROR: selected technology '$TECH' but $COVERAGE is missing or empty" >&2
            exit 1
        fi

        echo "[GABi] Selected technology: $TECH"
        echo "[GABi] Coverage source: $COVERAGE"

        python {params.script:q} \
            --coverage "$COVERAGE" \
            --sex-chr "{params.sex_chr}" \
            --tiny-contig-bp {params.tiny_contig_bp} \
            --output "{output.tsv}" \
            --sample-id "{wildcards.sample_id}"

        # DECLARED vs INFERRED CROSS-CHECK
        # -------------------------------------------------------------------
        # When sample_sex_karyo is declared, the karyotype decides ploidy and
        # this TSV is not consulted by the callers. It is still worth knowing
        # when coverage disagrees with the declaration - a swapped sample, a
        # mis-entered karyotype or a mis-declared sex contig all look like this.
        # The declaration still wins; this only reports.
        if [ "{params.sample_sex}" != "unknown" ]; then
            INFERRED=$(awk -F'\t' '
                NR==1 {{ for(i=1;i<=NF;i++){{ if($i=="contig") c=i; if($i=="call") k=i }}; next }}
                (k>0 && $k=="haploid") {{ print $c }}
            ' "{output.tsv}" | LC_ALL=C sort | paste -sd, - || true)
            DECLARED=$(printf '%s' "{params.declared_haploid}" | tr ',' '\n' \
                       | grep -v '^$' | LC_ALL=C sort | paste -sd, - || true)

            echo "[GABi] Declared sample_sex_karyo: {params.sample_sex}"
            echo "[GABi]   haploid by declaration: ${{DECLARED:-(none)}}"
            echo "[GABi]   haploid by coverage   : ${{INFERRED:-(none)}}"

            if [ "$INFERRED" != "$DECLARED" ]; then
                echo "[GABi] ⚠️  SEX MISMATCH for {wildcards.sample_id}: coverage disagrees with the declared karyotype."
                echo "[GABi]     Using the DECLARED karyotype for variant calling ({params.sample_sex})."
                echo "[GABi]     Check the sample's sex, its asm_sex_chr entry, and {output.tsv}"
            else
                echo "[GABi] Declaration and coverage agree"
            fi
        fi

        echo "[GABi] Done"
        '''




rule D06_aggregate_bam_qc:
    """Build a Markdown table combining per-sample BAM QC summaries."""
    input:
        _qc_summaries_for_assembly
    output:
        md=os.path.join(
            config["OUT_FOLDER"], "GABi_results", "{species}", "{assembly}",
            "samples", "bam_qc.md"
        )
    threads: cpu_func("light_task")
    resources:
        mem_mb=mem_func("light_task"),
        runtime=time_func("light_task")
    container: CONTAINERS["gabi_base"]
    log:
        os.path.join(
            config["OUT_FOLDER"], "GABi_results", "{species}", "{assembly}",
            "logs", "D06_aggregate_bam_qc.log"
        )
    params:
        aggregate_bam_qc = str(scripts_dir / "aggregate_bam_qc.py")
    shell:
        r'''
        set -euo pipefail
        mkdir -p "$(dirname {output.md})" "$(dirname {log})"
        python {params.aggregate_bam_qc:q} --out "{output.md}" {input} > "{log}" 2>&1
        '''


# ===============================================================================
#  TECHNOLOGY SELECTION
# ===============================================================================
#  Exactly one sequencing technology is chosen per sample. Technologies are
#  never merged.
#
#  The choice depends on per-technology coverage, which only exists after
#  mapping and D05 QC - so it cannot be made while the DAG is being built.
#  D07_select_tech writes it to tech_selection.json, and D08_link_selected_bam
#  points {sid}.selected.bam at the winner.
#
#  That symlink is the mechanism that keeps the DAG static: downstream rules
#  declare ONE fixed input path. If they instead declared
#  {sid}.{tech}.merged.bam, the technology would be part of a filename that is
#  unknown at DAG-build time, and an input function cannot resolve it - input
#  functions are evaluated exactly once, before any job has run. (params
#  functions ARE re-evaluated once inputs exist, which is why the caller rules
#  can read the chosen technology out of the JSON via params.)

def _d_callable_techs(w):
    """Technologies for this sample that can actually be mapped and called."""
    return callable_techs(w.species, w.assembly, w.sample_id)


def _d_tech_summaries(w):
    """D05 QC summaries for every callable technology of this sample."""
    qc = _qc_dir(w.species, w.assembly, w.sample_id)
    return [
        os.path.join(qc, f"{w.sample_id}.{t}.summary.txt")
        for t in _d_callable_techs(w)
    ]


def _d_tech_bams(w):
    """Per-technology merged BAMs (plus indexes) for every callable technology."""
    bdir = _bams_dir(w.species, w.assembly, w.sample_id)
    out = []
    for t in _d_callable_techs(w):
        out.append(os.path.join(bdir, f"{w.sample_id}.{t}.merged.bam"))
        out.append(os.path.join(bdir, f"{w.sample_id}.{t}.merged.bam.bai"))
    return out


rule D07_select_tech:
    """Pick one technology per sample, by DATA_PRIORITY and coverage.

    Reads MEAN_DEPTH from each technology's D05 QC summary. Technologies the
    configured CALLER cannot handle are excluded before ranking - that list is
    resolved at DAG-build time and passed in, so the script itself knows nothing
    about callers.
    """
    input:
        summaries=_d_tech_summaries,
    output:
        selection=os.path.join(
            config["OUT_FOLDER"], "GABi_results", "{species}", "{assembly}",
            "samples", "{sample_id}", "BAMs", "tech_selection.json"
        ),
    params:
        script=str(scripts_dir / "tech_selection.py"),
        # TECH=PATH pairs: no filename parsing, no reliance on argument order
        pairs=lambda w: " ".join(
            f"{t}={os.path.join(_qc_dir(w.species, w.assembly, w.sample_id), f'{w.sample_id}.{t}.summary.txt')}"
            for t in _d_callable_techs(w)
        ),
        eligible=",".join(caller_eligible_techs()),
        priority=config.get("DATA_PRIORITY", "hifi>illumina>ont"),
        good_cov=config.get("GOOD_TECH_COV", 15),
        min_cov=config.get("MIN_TECH_COV", 1),
    container: CONTAINERS["gabi_base"]
    threads: cpu_func("light_task")
    resources:
        mem_mb=mem_func("light_task"),
        runtime=time_func("light_task")
    log:
        os.path.join(
            config["OUT_FOLDER"], "GABi_results", "{species}", "{assembly}",
            "samples", "{sample_id}", "logs", "D07_select_tech.{sample_id}.log"
        )
    shell:
        r'''
        set -euo pipefail
        mkdir -p "$(dirname {output.selection})" "$(dirname {log})"

        python {params.script:q} \
            --sample-id "{wildcards.sample_id}" \
            --summaries {params.pairs} \
            --eligible "{params.eligible}" \
            --priority "{params.priority}" \
            --good-cov {params.good_cov} \
            --min-cov {params.min_cov} \
            --output "{output.selection}" \
            --log "{log}"
        '''


rule D08_link_selected_bam:
    """Point {sid}.selected.bam at the technology chosen by D07.

    A symlink, not a copy - these BAMs are large and the per-technology file is
    kept anyway. Links are relative (basename) so they survive the different
    mount points inside and outside the container.
    """
    input:
        selection=rules.D07_select_tech.output.selection,
        bams=_d_tech_bams,
    output:
        bam=os.path.join(
            config["OUT_FOLDER"], "GABi_results", "{species}", "{assembly}",
            "samples", "{sample_id}", "BAMs", "{sample_id}.selected.bam"
        ),
        bai=os.path.join(
            config["OUT_FOLDER"], "GABi_results", "{species}", "{assembly}",
            "samples", "{sample_id}", "BAMs", "{sample_id}.selected.bam.bai"
        ),
    container: CONTAINERS["gabi_base"]
    threads: cpu_func("light_task")
    resources:
        mem_mb=mem_func("light_task"),
        runtime=time_func("light_task")
    log:
        os.path.join(
            config["OUT_FOLDER"], "GABi_results", "{species}", "{assembly}",
            "samples", "{sample_id}", "logs", "D08_link_selected_bam.{sample_id}.log"
        )
    shell:
        r'''
        set -euo pipefail
        mkdir -p "$(dirname {output.bam})" "$(dirname {log})"
        exec > "{log}" 2>&1

        BAM_DIR="$(dirname "{output.bam}")"

        TECH=$(python -c "import json,sys; print(json.load(open(sys.argv[1]))['chosen'])" "{input.selection}")
        if [ -z "$TECH" ]; then
            echo "[GABi] ERROR: no 'chosen' technology in {input.selection}" >&2
            exit 1
        fi

        SRC="{wildcards.sample_id}.${{TECH}}.merged.bam"
        if [ ! -e "$BAM_DIR/$SRC" ]; then
            echo "[GABi] ERROR: selected technology '$TECH' but $BAM_DIR/$SRC does not exist" >&2
            exit 1
        fi

        echo "[GABi] {wildcards.sample_id}: selected technology '$TECH'"
        echo "[GABi] Linking $(basename "{output.bam}") -> $SRC"

        ln -sfn "$SRC"       "{output.bam}"
        ln -sfn "$SRC.bai"   "{output.bai}"

        # -e follows symlinks, so this confirms the targets actually resolve
        if [ ! -e "{output.bam}" ] || [ ! -e "{output.bai}" ]; then
            echo "[GABi] ERROR: symlink created but does not resolve" >&2
            exit 1
        fi

        echo "[GABi] Done"
        '''
