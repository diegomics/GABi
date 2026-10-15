# ===============================================================================
# GABi - Variant Calling Rules
# by Diego De Panis, 2026
# note: AI tools may have been used to improve, clean and/or comment this version of the code
# ===============================================================================

# -------------------------------------------------------------------------------
#  OUTPUTH PATH TEMPLATES (literal patterns only)
# -------------------------------------------------------------------------------

IDX_DIR_TMPL        = os.path.join(config["OUT_FOLDER"], "GABi_results", "{species}", "{assembly}", "idx")
REF_FA_TMPL         = os.path.join(IDX_DIR_TMPL, "{assembly}.fa")
REF_FAI_TMPL        = os.path.join(IDX_DIR_TMPL, "{assembly}.fa.fai")
REF_DICT_TMPL       = os.path.join(IDX_DIR_TMPL, "{assembly}.dict")
INTERVALS_DIR_TMPL  = os.path.join(IDX_DIR_TMPL, "intervals")

SAMPLE_BAMS_TMPL    = os.path.join(config["OUT_FOLDER"], "GABi_results", "{species}", "{assembly}", "samples", "{sample_id}", "BAMs")
SAMPLE_VCFS_TMPL    = os.path.join(config["OUT_FOLDER"], "GABi_results", "{species}", "{assembly}", "samples", "{sample_id}", "VCFs")

INTERVAL_LIST_TMPL  = os.path.join(INTERVALS_DIR_TMPL, "interval_{shard}.list")
GVCF_PART_TMPL      = os.path.join(SAMPLE_VCFS_TMPL, "tmp", "{sample_id}.bp.{shard}.g.vcf.gz")
GVCF_PART_TBI_TMPL  = GVCF_PART_TMPL + ".tbi"
GENO_PART_TMPL      = os.path.join(SAMPLE_VCFS_TMPL, "tmp", "{sample_id}.bp.{shard}.genotyped.vcf.gz")
GENO_PART_TBI_TMPL  = GENO_PART_TMPL + ".tbi"

# FINAL_* templates are defined after _BP_MODE is set (see "Resolution-aware final paths" below)

# Caller-specific intermediates (routed to FINAL_* by E20/E21)
GATK_GVCF_TMPL     = os.path.join(SAMPLE_VCFS_TMPL, "{sample_id}.gatk.g.vcf.gz")
GATK_GVCF_TBI_TMPL = GATK_GVCF_TMPL + ".tbi"
GATK_BCF_TMPL      = os.path.join(SAMPLE_VCFS_TMPL, "{sample_id}.gatk.raw.bcf")
GATK_BCF_CSI_TMPL  = GATK_BCF_TMPL + ".csi"

DV_GVCF_TMPL       = os.path.join(SAMPLE_VCFS_TMPL, "{sample_id}.dv.g.vcf.gz")
DV_GVCF_TBI_TMPL   = DV_GVCF_TMPL + ".tbi"
DV_BCF_TMPL        = os.path.join(SAMPLE_VCFS_TMPL, "{sample_id}.dv.raw.bcf")
DV_BCF_CSI_TMPL    = DV_BCF_TMPL + ".csi"

BCFT_GVCF_TMPL     = os.path.join(SAMPLE_VCFS_TMPL, "{sample_id}.bcft.g.vcf.gz")
BCFT_GVCF_TBI_TMPL = BCFT_GVCF_TMPL + ".tbi"
BCFT_BCF_TMPL      = os.path.join(SAMPLE_VCFS_TMPL, "{sample_id}.bcft.raw.bcf")
BCFT_BCF_CSI_TMPL  = BCFT_BCF_TMPL + ".csi"
BCFT_PART_TMPL     = os.path.join(SAMPLE_VCFS_TMPL, "tmp", "{sample_id}.bcft.{shard}.bcf")
BCFT_PART_CSI_TMPL = BCFT_PART_TMPL + ".csi"

# Sex-aware ploidy for bcftools call, written per sample by E29.
BCFT_PLOIDY_TMPL   = os.path.join(SAMPLE_VCFS_TMPL, "{sample_id}.ploidy.txt")
BCFT_SAMPLES_TMPL  = os.path.join(SAMPLE_VCFS_TMPL, "{sample_id}.samples.txt")

# Joint genotyping output (per species/assembly, not per sample)
JOINT_VCFS_DIR_TMPL = os.path.join(config["OUT_FOLDER"], "GABi_results", "{species}", "{assembly}", "jointVCFs")
JOINT_BCF_TMPL      = os.path.join(JOINT_VCFS_DIR_TMPL, "{assembly}.joint.bcf")
JOINT_BCF_CSI_TMPL  = JOINT_BCF_TMPL + ".csi"

# Cohort-level files for JOINT_GENO: multi (bcftools multi-sample calling)
COHORT_PLOIDY_TMPL  = os.path.join(JOINT_VCFS_DIR_TMPL, "{assembly}.cohort.ploidy.txt")
COHORT_SAMPLES_TMPL = os.path.join(JOINT_VCFS_DIR_TMPL, "{assembly}.cohort.samples.txt")
COHORT_PART_TMPL    = os.path.join(JOINT_VCFS_DIR_TMPL, "tmp", "{assembly}.multi.{shard}.bcf")
COHORT_PART_CSI_TMPL = COHORT_PART_TMPL + ".csi"


# -------------------------------------------------------------------------------
#  HELPER FUNCTIONS (prefixed to avoid collisions)
# -------------------------------------------------------------------------------

def vc_intervals_dir(sp, asm):
    return os.path.join(config["OUT_FOLDER"], "GABi_results", sp, asm, "idx", "intervals")

# NOTE: priority_list(), sample_techs() and pick_priority_tech() are defined in
# the Snakefile HELPERS section and shared with D_map_reads.smk.

# Path templates for the technology-selection artefacts produced by D07/D08.
# Every caller reads the same two paths, whatever technology was selected: the
# selected BAM is a symlink, so the technology never appears in a filename that
# the DAG has to know in advance.
TECH_SELECTION_TMPL = os.path.join(
    SAMPLE_BAMS_TMPL, "tech_selection.json"
)
SELECTED_BAM_TMPL = os.path.join(
    SAMPLE_BAMS_TMPL, "{sample_id}.selected.bam"
)


def vc_selection_json(sp, asm, sid):
    """Path to this sample's tech_selection.json (written by D07_select_tech)."""
    return TECH_SELECTION_TMPL.format(species=sp, assembly=asm, sample_id=sid)


def vc_get_calling_bam(sp, asm, sid):
    """BAM to call variants on: always the selected-technology symlink.
    """
    return SELECTED_BAM_TMPL.format(species=sp, assembly=asm, sample_id=sid)


def vc_selected_tech(sp, asm, sid, default=None):
    """Read the selected technology out of tech_selection.json.

    Only safe from a `params` function or a shell block, never from an `input`
    function: Snakemake evaluates input functions exactly once while building
    the DAG, before any job has run, whereas params are re-evaluated once a
    job's inputs are present (Scheduler.reset_params_and_resources). The
    `default` fallback covers that first, pre-run evaluation and `--dry-run`,
    where params are deliberately not refreshed.
    """
    try:
        with open(vc_selection_json(sp, asm, sid)) as fh:
            return json.load(fh).get("chosen") or default
    except (OSError, ValueError, KeyError):
        return default


def vc_interval_path(sp, asm, shard):
    return os.path.join(vc_intervals_dir(sp, asm), f"interval_{shard}.list")

# Assemblies whose interval layout has already been checked this session.
_INTERVALS_CHECKED = set()


def _vc_assert_intervals_current(sp, asm):
    """Fail if intervals/ was built for a different asm_sex_chr than is declared.

    The shard layout depends on asm_sex_chr (E00 gives each sex contig its own
    shard so E01 can set ploidy per sample). That dependency is declared as a
    param, which Snakemake uses as a rerun trigger - but only when `params` is
    among the configured rerun-triggers. GABi's execution profiles set
    `rerun-triggers: [mtime]`, under which an existing intervals/ directory is
    reused no matter what asm_sex_chr says.

    Left undetected that is a silent correctness bug: sex contigs stay bundled
    with autosomes, E01 finds a mixed shard, falls back to diploid, and emits
    perfectly ordinary-looking 0/0 0/1 1/1 genotypes on a haploid chromosome.
    E01 refuses to run in that state, but only after the job has started; this
    catches it while the DAG is still being built.
    """
    if (sp, asm) in _INTERVALS_CHECKED:
        return
    manifest = os.path.join(vc_intervals_dir(sp, asm), "MANIFEST.json")
    declared = sorted(get_sex_chr(samples_config, sp, asm))
    try:
        with open(manifest) as fh:
            recorded = sorted(json.load(fh).get("sex_contigs", []))
    except (OSError, ValueError):
        # No manifest yet (first build) - nothing to compare against.
        return
    if recorded != declared:
        raise ValueError(
            f"[GABi] The interval layout for {sp}/{asm} is out of date.\n"
            f"       built for sex contigs : {', '.join(recorded) or '(none)'}\n"
            f"       currently declared    : {', '.join(declared) or '(none)'}\n"
            f"       Sex contigs need a shard of their own so GATK can call them\n"
            f"       at the right ploidy per sample. Rebuild the intervals:\n"
            f"         rm -rf {vc_intervals_dir(sp, asm)}\n"
            f"       or:\n"
            f"         snakemake --forcerun E00_make_intervals ...\n"
            f"       (GABi's profiles use rerun-triggers: [mtime], so Snakemake\n"
            f"        will not notice this on its own.)"
        )
    _INTERVALS_CHECKED.add((sp, asm))


def vc_discover_shards(sp, asm):
    # force interval checkpoint to run
    _ = checkpoints.E00_make_intervals.get(species=sp, assembly=asm)
    _vc_assert_intervals_current(sp, asm)
    files = sorted(glob.glob(os.path.join(vc_intervals_dir(sp, asm), "interval_*.list")))
    return [os.path.splitext(os.path.basename(x))[0].split("_", 1)[1] for x in files]

# config knobs
# GABi assumes diploid samples, so ploidy is fixed rather than configurable.
# Sex chromosomes are handled per caller, not through this value (DeepVariant
# uses --haploid_contigs, fed by D06's coverage-based sex inference).
# This is also GATK HaplotypeCaller's own default for --sample-ploidy, so it is
# passed explicitly only to keep the assumption visible in the command line.
_SAMPLE_PLOIDY   = 2
_KEEP_BP_GVCF    = bool(_as_bool(config.get("KEEP_gVCF", True)))

# Resolution-aware final paths
# basepair -> {sid}.bp.g.vcf.gz / {sid}.bp.raw.bcf  (name reflects per-site content)
# block    -> {sid}.g.vcf.gz    / {sid}.raw.bcf      (no prefix needed)
##_BP_TAG = ".bp" if _BP_MODE else ""

FINAL_GVCF_TMPL     = os.path.join(SAMPLE_VCFS_TMPL, "{sample_id}" + _BP_TAG + ".g.vcf.gz")
FINAL_GVCF_TBI_TMPL = FINAL_GVCF_TMPL + ".tbi"
FINAL_BCF_TMPL      = os.path.join(SAMPLE_VCFS_TMPL, "{sample_id}" + _BP_TAG + ".raw.bcf")
FINAL_BCF_CSI_TMPL  = FINAL_BCF_TMPL + ".csi"

# Variant caller: "deep" (DeepVariant, default), "bcft" (bcftools) or "gatk".
# normalize_caller() lives in the Snakefile HELPERS section and accepts the
# common spellings/capitalisations, raising on anything unrecognised.
_CALLER = normalize_caller(config.get("CALLER", "deep"))

_CALLER_LABEL = {
    "deep": "DeepVariant",
    "bcft": "bcftools",
    "gatk": "GATK HaplotypeCaller",
}[_CALLER]
print(f"[GABi] Variant caller: {_CALLER_LABEL}")

# Joint genotyping via GLnexus (requires CALLER: deep)
# YAML parses on/off as True/False — handle both booleans and strings
# JOINT_GENO parsing lives in the Snakefile alongside normalize_caller(), so
# the output collector and this module cannot disagree about what the setting
# means
_JOINT_GENO_MODE = joint_geno_mode()
_JOINT_GENO = _JOINT_GENO_MODE != "off"

# How each caller's per-sample gVCFs are combined into a cohort call set:
#   deep -> GLnexus --config DeepVariant
#   gatk -> GLnexus --config gatk        (GLnexus ships presets for both)
#   bcft -> bcftools merge --gvcf        (GLnexus has no bcftools preset)
#
# All three MERGE the per-sample gVCFs rather than re-calling from BAMs, so the
# joint call set is always derived from the same evidence as the per-sample
# outputs. It also means a cohort can mix technologies: each sample was already
# called with the profile appropriate to its own data.

# JOINT_GENO: multi asks for true multi-sample calling instead, which only
# bcftools offers - DeepVariant and GATK have no equivalent mode.
if _JOINT_GENO_MODE == "multi" and _CALLER != "bcft":
    raise ValueError(
        f"[GABi] JOINT_GENO: multi requires CALLER: bcft (got '{_CALLER}').\n"
        f"       True multi-sample calling means `bcftools mpileup -b <bams> | "
        f"bcftools call -m`;\n"
        f"       DeepVariant and GATK have no multi-sample calling mode.\n"
        f"       Use JOINT_GENO: merge to combine their per-sample gVCFs instead."
    )

if _JOINT_GENO_MODE == "multi":
    _JOINT_METHOD = "bcftools-multisample"
elif _CALLER == "bcft":
    _JOINT_METHOD = "bcftools-merge"
else:
    _JOINT_METHOD = "glnexus"

# GLnexus reads a per-sample reference depth from every non-variant record, and
# which FORMAT field holds it depends on how the gVCF was made:
#
#   deep, any resolution   FORMAT/MIN_DP  - E10 always emits DeepVariant's own
#                                           banded gVCF; only E11's BCF expands,
#                                           so RESOLUTION never changes it
#   gatk + block           FORMAT/MIN_DP  - -ERC GVCF writes reference bands
#   gatk + basepair        FORMAT/DP      - -ERC BP_RESOLUTION writes one record
#                                           per position, and a single-position
#                                           record has no band minimum
#
# The stock "gatk" preset hardcodes ref_dp_format: MIN_DP, so basepair gVCFs
# need a preset of their own. It is NOT used for block mode: GATK bands carry
# both DP and MIN_DP, so a DP-based config would run there too but would take
# each band's average instead of its minimum - less conservative, and a silent
# change to results.
if _CALLER == "gatk" and _BP_MODE:
    _GLNEXUS_CONFIG = str(scripts_dir / "glnexus_gatk_bp.yml")
else:
    _GLNEXUS_CONFIG = {"deep": "DeepVariant", "gatk": "gatk"}.get(_CALLER, "DeepVariant")

# Precomputed rather than written inline as a conditional expression in the
# rule's `container:` directive. Snakemake's parser reads a directive value
# token by token and treats the line break as the end of the value, so a
# parenthesised expression continued onto a line beginning with `else` makes it
# read `else` as a bare rule keyword: "Unexpected keyword else in rule
# definition". Keep directive values to a single line, or hoist them here.
_JOINT_CONTAINER = (
    CONTAINERS["glnexus"] if _JOINT_METHOD == "glnexus" else CONTAINERS["gabi_base"]
)

# How the cohort call set was produced, for the joint reports in F.
_JOINT_CALLER_LABEL = {
    ("deep", "glnexus"): "DeepVariant + GLnexus",
    ("gatk", "glnexus"): "GATK HaplotypeCaller + GLnexus",
    ("bcft", "bcftools-merge"): "bcftools + bcftools merge --gvcf",
    ("bcft", "bcftools-multisample"): "bcftools multi-sample calling",
}.get((_CALLER, _JOINT_METHOD), f"{_CALLER_LABEL} + {_JOINT_METHOD}")

if _JOINT_GENO:
    # Auto-force KEEP_gVCF on (the per-sample gVCFs are the input either way)
    _keep_raw = _as_bool(config.get("KEEP_gVCF", True))
    _keep_on = _keep_raw is True or str(_keep_raw).strip().lower() in ("true", "yes", "on", "1")
    if not _keep_on:
        print("[GABi] JOINT_GENO=on -> setting KEEP_gVCF=on (required for joint genotyping)")
        config["KEEP_gVCF"] = True
    _KEEP_BP_GVCF = True
    _JOINT_LABEL = {
        "glnexus": f"GLnexus, --config {_GLNEXUS_CONFIG}",
        "bcftools-merge": "bcftools merge --gvcf",
        "bcftools-multisample": "bcftools multi-sample calling (mpileup -b | call -mv)",
    }[_JOINT_METHOD]
    print(f"[GABi] Joint genotyping: enabled ({_JOINT_LABEL})")
else:
    print("[GABi] Joint genotyping: disabled")


# -------------------------------------------------------------------------------
#  CALLING CONFIG VALIDATION (fail fast, at DAG-build time)
# -------------------------------------------------------------------------------
# Rejects caller / technology / resolution combinations that cannot produce
# correct results, so a run stops here rather than after hours of compute

def _validate_calling_config():
    """Check CALLER against each sample's available technologies."""

    # --- note on caller x resolution (nothing is rejected here) -----------------
    # Every caller is valid at either resolution, but they reach basepair by
    # different routes:
    #   bcft - omit `call --gvcf`, so mpileup/call emit one record per position
    #   gatk - -ERC BP_RESOLUTION, one record per position from HaplotypeCaller
    #   deep - DeepVariant only emits reference blocks, so E11 expands them with
    #          `bcftools convert --gvcf2vcf`. That loop (gvcf_to_vcf() in
    #          vcfconvert.c) rewrites a single record once per position, altering
    #          only POS and the REF base: FORMAT/GT, GQ, MIN_DP, PL and VAF are
    #          replicated unchanged across the whole block. Per-site depth is
    #          therefore the block MINIMUM rather than the depth at that exact
    #          position - conservative for F01's depth filters, and intended.
    #          It is also the only route to the MASKED INFO tag, which F01
    #          applies in basepair mode only.

    # Nothing below matters if we are not calling at all.
    if not _as_bool(config.get("RUN_CALLING", True)):
        return

    # --- caller x per-sample technologies --------------------------------------
    unusable, narrowed, gatk_hifi = [], [], []

    for sp, asm, sid in all_samples_iter():
        techs = [t for t in sample_techs(sp, asm, sid) if t not in NON_CALLABLE_TECHS]
        if not techs:
            continue
        eligible = [t for t in techs if t in CALLER_TECHS[_CALLER]]
        if not eligible:
            unusable.append((sp, asm, sid, techs))
            continue
        excluded = [t for t in techs if t not in eligible]
        if excluded:
            narrowed.append((sid, excluded, eligible))
        if _CALLER == "gatk" and "hifi" in eligible:
            gatk_hifi.append(sid)

    if unusable:
        detail = "\n".join(
            f"         {sid} ({sp}/{asm}): has {', '.join(techs)}"
            for sp, asm, sid, techs in unusable
        )
        raise ValueError(
            f"[GABi] CALLER: {_CALLER} cannot be used with any technology available\n"
            f"       for these samples:\n{detail}\n"
            f"       {_CALLER} supports: {', '.join(sorted(CALLER_TECHS[_CALLER]))}.\n"
            f"       Use CALLER: deep (recommended for every technology) or\n"
            f"       CALLER: bcft, or drop these samples from the data table."
        )

    # --- multi-sample joint calling x cohort technology uniformity ------------
    # `bcftools mpileup` applies ONE -X profile to every BAM it is given, so
    # true multi-sample calling is only sound when every sample in a cohort was
    # called on the same technology. Uniformity only has to hold WITHIN a
    # (species, assembly) cohort - different species may use different techs.
    #
    # Selection happens at run time, so this cannot always be settled here. What
    # can: if every sample in a cohort has exactly ONE candidate technology,
    # the outcome is already determined. Those cases are rejected now; the rest
    # are re-checked against the actual tech_selection.json files in EJ21.
    if _JOINT_METHOD == "bcftools-multisample":
        cohorts = {}
        for sp, asm, sid in all_samples_iter():
            techs = callable_techs(sp, asm, sid, _CALLER)
            if techs:
                cohorts.setdefault((sp, asm), []).append((sid, techs))

        for (sp, asm), members in cohorts.items():
            if len(members) < 2:
                continue          # joint genotyping needs >=2 samples anyway
            if not all(len(t) == 1 for _, t in members):
                continue          # outcome depends on coverage - checked at run time
            chosen = {t[0] for _, t in members}
            if len(chosen) > 1:
                detail = "\n".join(
                    f"         {sid}: {t[0]}" for sid, t in sorted(members)
                )
                raise ValueError(
                    f"[GABi] JOINT_GENO: multi needs every sample in a cohort to use the\n"
                    f"       same technology, but {sp}/{asm} mixes "
                    f"{', '.join(sorted(chosen))}:\n{detail}\n"
                    f"       bcftools mpileup applies one -X profile to all BAMs, so the\n"
                    f"       wrong error model would be used for some samples.\n"
                    f"       Use JOINT_GENO: merge, which calls each sample with its own\n"
                    f"       profile and then combines the gVCFs."
                )

    if IS_MAIN:
        for sid, excluded, eligible in narrowed:
            print(f"[GABi] Sample {sid}: {', '.join(excluded)} cannot be used with "
                  f"CALLER: {_CALLER} -> calling restricted to {' / '.join(eligible)}")
        for sid in gatk_hifi:
            print(f"[GABi] ⚠️  Sample {sid}: GATK HaplotypeCaller on HiFi reads is "
                  f"discouraged. DeepVariant (PACBIO model) or bcftools "
                  f"(-X pacbio-ccs) handle HiFi error profiles better.")


_validate_calling_config()


def _get_sample_read_category(species, assembly, sample_id):
    """'short' or 'long', for the technology actually being called on.
    Only one technology is used per sample now.
    Read from `params`, so tech_selection.json is present by the time the value
    is used; the fallback covers DAG build and dry runs.
    """
    tech = vc_selected_tech(species, assembly, sample_id)
    if tech is None:
        tech = pick_priority_tech(
            species, assembly, sample_id,
            techs=callable_techs(species, assembly, sample_id, _CALLER),
        )
    if tech is None:
        return "short"
    return "long" if _is_long(tech) else "short"


def _vc_pairhmm_impl(read_category):
    """
    Determine PairHMM implementation based on read category and USE_AVX
    (set in the Snakefile, next to ODB_VERSION).

    USE_AVX options:
      auto: FASTEST_AVAILABLE for short reads, LOGLESS_CACHING for long reads
      on:   FASTEST_AVAILABLE for all (may crash on some CPUs or with long reads)
      off:  LOGLESS_CACHING for all (safest, slowest)
    """
    mode = str(USE_AVX).lower()

    if mode in {"off", "false", "no"}:
        # Disabled: always use safe implementation
        return "LOGLESS_CACHING"

    if mode in {"on", "true", "yes"}:
        # Forced on: use AVX for everything (user takes responsibility)
        return "FASTEST_AVAILABLE"

    # Auto mode (default): AVX for short reads only
    if read_category == "long":
        return "LOGLESS_CACHING"
    else:
        return "FASTEST_AVAILABLE"


def vc_hc_params(wildcards):
    """
    Return HaplotypeCaller parameters based on read type.

    Short reads: optimized for high coverage, PCR artifacts present
    Long reads: optimized for longer reads, no PCR, safer pairhmm
    """
    read_category = _get_sample_read_category(
        wildcards.species, wildcards.assembly, wildcards.sample_id
    )

    pairhmm = _vc_pairhmm_impl(read_category)

    if read_category == "long":
        return {
            "read_category": "long",
            "pairhmm": pairhmm,
            "pcr_indel_model": "NONE",            # No PCR artifacts in long reads
            "max_haplotypes": 64,                 # Fewer haplotypes needed
            "max_reads_start": 40,                # Each read covers more
            "max_asm_region": 200,                # Smaller assembly regions
        }
    else:  # short reads
        return {
            "read_category": "short",
            "pairhmm": pairhmm,
            "pcr_indel_model": "CONSERVATIVE",    # Handle PCR artifacts
            "max_haplotypes": 128,                # More haplotypes for high coverage
            "max_reads_start": 50,                # Standard for short reads
            "max_asm_region": 300,                # Larger assembly regions
        }



# -------------------------------------------------------------------------------
#  DeepVariant helpers
# -------------------------------------------------------------------------------

# Samples already announced by _caller_for_sample(). That function is called
# once per routed file - E20 (gVCF + tbi) and E21 (BCF + csi), plus F01's params
# - so without this guard a multi-tech sample logs the same line 4-5 times.
_ANNOUNCED_MULTITECH = set()


def _caller_for_sample(sp, asm, sid):
    """Caller to use for a sample, decided at DAG-build time.
    This is simply the CALLER setting.

    Technology/caller compatibility is enforced up front by
    _validate_calling_config(), so by the time this runs every sample is known
    to have at least one technology the configured caller can handle.
    """
    if IS_MAIN and (sp, asm, sid) not in _ANNOUNCED_MULTITECH:
        _ANNOUNCED_MULTITECH.add((sp, asm, sid))
        techs = [t for t in sample_techs(sp, asm, sid) if t not in NON_CALLABLE_TECHS]
        if len(techs) > 1:
            eligible = [t for t in techs if t in CALLER_TECHS[_CALLER]]
            chosen = pick_priority_tech(sp, asm, sid, techs=eligible)
            print(f"[GABi] Sample {sid}: {len(techs)} technologies "
                  f"({', '.join(techs)}) -> {_CALLER_LABEL} on '{chosen}'")
    return _CALLER


def _tech_to_dv_model(tech):
    """Map a read technology onto a DeepVariant --model_type."""
    if tech == "hifi":
        return "PACBIO"
    if tech == "ont":
        return "ONT_R104"
    return "WGS"          # illumina, 10x


def _dv_model_type(sp, asm, sid):
    """DeepVariant --model_type for the SELECTED technology.

    Used from `params`, so by the time the value reaches the shell,
    tech_selection.json exists and this returns the real answer. The
    DAG-build-time fallback below is only ever seen during dry runs; E10's
    shell re-reads the JSON and refuses to proceed if the two disagree.
    """
    tech = vc_selected_tech(sp, asm, sid)
    if tech is None:
        # Not yet selected (DAG build / dry run) - fall back to the priority
        # technology purely so the dry-run output shows something sensible.
        tech = pick_priority_tech(
            sp, asm, sid, techs=callable_techs(sp, asm, sid, _CALLER)
        )
    return _tech_to_dv_model(tech) if tech else "WGS"


def _dv_get_bam(sp, asm, sid):
    """BAM for DeepVariant: the selected-technology symlink from D08."""
    return SELECTED_BAM_TMPL.format(species=sp, assembly=asm, sample_id=sid)


# -------------------------------------------------------------------------------
#  bcftools helpers
# -------------------------------------------------------------------------------

# mpileup -X platform profiles (bcftools >= 1.20). The versioned names are used
# deliberately: the bare aliases ("illumina", "ont-sup", "pacbio-ccs") float to
# whatever the newest tuning is, which would silently change calls between
# bcftools releases.
#
# NOTE on ont-sup-1.20: it is tuned for SUP (super-accuracy) basecalled R10
# data and DOES call indels via the CNS model. The older plain "ont" profile
# sets -I (skip indels entirely). If ONT reads were basecalled with fast/hac,
# or come from R9, indel calls from this profile should be treated with
# suspicion. GABi's data table does not record basecaller or chemistry, so this
# is an assumption rather than something we can check.
_BCFT_PROFILES = {
    "illumina": "illumina-1.20",
    "10x":      "illumina-1.20",
    "hifi":     "pacbio-ccs-1.20",
    "ont":      "ont-sup-1.20",
}

# DP thresholds delimiting gVCF reference blocks (block resolution only).
# Sites below the first value are emitted as individual records. Deliberately
# not a config knob: it is a file-size/resolution trade-off, not a scientific
# one, and the control panel is already long enough.
_BCFT_GVCF_BINS = "4,8,16,24,48,96"


def _bcft_profile(tech):
    """mpileup -X profile for a read technology."""
    return _BCFT_PROFILES.get(tech, "illumina-1.20")


def _bcft_profile_for_sample(sp, asm, sid):
    """mpileup -X profile for the SELECTED technology.

    Read from `params`, so tech_selection.json exists by the time the value
    reaches the shell; the fallback only shows up in dry runs. E30's shell
    re-reads the JSON and refuses to run if the two disagree.
    """
    tech = vc_selected_tech(sp, asm, sid)
    if tech is None:
        tech = pick_priority_tech(
            sp, asm, sid, techs=callable_techs(sp, asm, sid, _CALLER)
        )
    return _bcft_profile(tech)


def _bcft_haploid_contigs(w):
    """Contigs to call at ploidy 1, via the shared resolver.

    Defined in terms of haploid_contigs_for_sample() rather than repeating its
    logic, so bcftools and DeepVariant cannot disagree about a sample's ploidy.
    """
    return haploid_contigs_for_sample(w.species, w.assembly, w.sample_id)


# -------------------------------------------------------------------------------
#  Caller routing (select caller-specific intermediate -> final output)
# -------------------------------------------------------------------------------

# -------------------------------------------------------------------------------
#  Path template for D06's output (already defined in D, but re-state here
#  for clarity since E references it). Keep in sync with D06_infer_sex_ploidy!!
# -------------------------------------------------------------------------------
SEX_INFER_TSV_TMPL = os.path.join(
    config["OUT_FOLDER"], "GABi_results", "{species}", "{assembly}",
    "samples", "{sample_id}", "BAMs", "qc", "{sample_id}.sex_inference.tsv"
)


# -------------------------------------------------------------------------------
#  Reader: extract haploid contigs from D06's sex_inference.tsv
# -------------------------------------------------------------------------------
def _read_inferred_haploid_contigs(tsv_path):
    """
    Parse sex_inference.tsv and return the list of contigs whose call is
    "haploid". Ambiguous and diploid rows are excluded — they are NOT passed
    to --haploid_contigs (safe default per design).

    Returns [] if the file is empty (header only) or unreadable.
    """
    contigs = []
    try:
        with open(tsv_path) as f:
            header = f.readline().strip().split("\t")
            if "contig" not in header or "call" not in header:
                return []
            i_contig = header.index("contig")
            i_call = header.index("call")
            for line in f:
                parts = line.rstrip("\n").split("\t")
                if len(parts) <= max(i_contig, i_call):
                    continue
                if parts[i_call].strip() == "haploid":
                    contigs.append(parts[i_contig].strip())
    except (FileNotFoundError, IOError):
        # Should not happen at rule-execution time (Snakemake guarantees the
        # input file exists), but degrade safely if it does.
        return []
    return contigs


# -------------------------------------------------------------------------------
#  Input function: conditionally depend on sex_inference.tsv
# -------------------------------------------------------------------------------
def _dv_sex_inference_input(w):
    """
    Return the path to sex_inference.tsv as a list — empty list if this
    sample does not need inference (declared sample_sex, missing sex_chr,
    or "none" case).

    Returning [] tells Snakemake "no extra dependency" — D06 will not be
    invoked for this sample.
    """
    contigs, source = resolve_haploid_contigs(
        samples_config, w.species, w.assembly, w.sample_id
    )
    if source == "inferred":
        return [SEX_INFER_TSV_TMPL.format(
            species=w.species, assembly=w.assembly, sample_id=w.sample_id
        )]
    # "declared", "missing_sex_chr", or "none" — no TSV needed.
    return []


# -------------------------------------------------------------------------------
#  Params function: build the --haploid_contigs flag string for make_examples
# -------------------------------------------------------------------------------
def haploid_contigs_for_sample(sp, asm, sid):
    """Contigs to call at ploidy 1 for one sample. THE single source of truth.

    Every caller must resolve haploid contigs through this function.

    Resolution, in order:
      - sample_sex_karyo declared hetero -> every contig in asm_sex_chr.
        The karyotype wins; sex_inference.tsv is NOT consulted.
      - sample_sex_karyo declared homo   -> none. Karyotype wins.
      - sample_sex_karyo unknown/empty          -> read sex_inference.tsv and
        take only rows whose call is "haploid". "ambiguous" and "diploid" are
        both left diploid, which is the safe direction.
      - no asm_sex_chr on the assembly          -> none.
    """
    contigs, source = resolve_haploid_contigs(samples_config, sp, asm, sid)

    if source == "inferred":
        tsv = SEX_INFER_TSV_TMPL.format(
            species=sp, assembly=asm, sample_id=sid
        )
        contigs = _read_inferred_haploid_contigs(tsv)

    return contigs or []


def _dv_haploid_contigs_flag(w):
    """--haploid_contigs flag for DeepVariant's make_examples, or "" for none."""
    contigs = haploid_contigs_for_sample(w.species, w.assembly, w.sample_id)
    if not contigs:
        return ""
    return f"--haploid_contigs={','.join(contigs)}"



def vc_shard_ploidy(w):
    """--sample-ploidy for one sample on one shard: 1 or _SAMPLE_PLOIDY.

    E00 gives every sex contig a shard of its own, so a shard is either wholly
    sex-linked or wholly autosomal. A shard is called haploid only when EVERY
    contig in it is haploid for this sample; anything else stays diploid, which
    keeps the safe direction for a sample whose inference came back ambiguous.

    Read from `params`, so the interval file and sex_inference.tsv both exist by
    the time the value is used. The fallback below is only reached during DAG
    construction and dry runs, and E01's shell re-derives the value to check.

    GATK cannot express ploidy 0 (a W in a ZZ male, a Y in an XX female).
    haploid_contigs_for_sample() only ever reports contigs to call at ploidy 1,
    so the situation never arises here - those contigs are called diploid and
    their near-zero coverage leaves the calls to be filtered downstream.
    """
    haploid = set(haploid_contigs_for_sample(w.species, w.assembly, w.sample_id))
    if not haploid:
        return _SAMPLE_PLOIDY
    try:
        with open(vc_interval_path(w.species, w.assembly, w.shard)) as fh:
            shard_contigs = {ln.rsplit(":", 1)[0] for ln in fh if ln.strip()}
    except OSError:
        return _SAMPLE_PLOIDY
    if shard_contigs and shard_contigs <= haploid:
        return 1
    return _SAMPLE_PLOIDY


def _routed_caller_file(w, tmpl_map):
    """Route to the correct caller-specific intermediate file.

    tmpl_map is keyed by canonical caller name (deep / bcft / gatk).
    """
    caller = _caller_for_sample(w.species, w.assembly, w.sample_id)
    if caller not in tmpl_map:
        raise ValueError(
            f"[GABi] No output template registered for caller '{caller}'."
        )
    return tmpl_map[caller].format(
        species=w.species, assembly=w.assembly, sample_id=w.sample_id
    )


# -------------------------------------------------------------------------------
#  REF READINESS (depend on map_reads rule)
# ===============================================================================

rule E00_ref_faidx_for_calling:
    input:
        fa=rules.D00_ref_fa_canonical.output.fa
    output:
        fai=REF_FAI_TMPL
    threads: cpu_func("samtools_index")
    resources:
        mem_mb=mem_func("samtools_index"),
        runtime=time_func("samtools_index")
    container: CONTAINERS["gabi_base"]
    log:
        os.path.join(config["OUT_FOLDER"], "GABi_results", "{species}", "{assembly}", "logs", "E00_ref_faidx.log")
    shell:
        r'''
        set -euo pipefail
        mkdir -p $(dirname {log})

        exec > {log} 2>&1

        echo "[GABi] Creating FASTA index for {input.fa}"


        # TEMP DIRECTORY SETUP
        # -------------------------------------------------------------------
        WORK_DIR="$(gabi_get_workdir 25)"
        TEMP_DIR="$(mktemp -d "$WORK_DIR/GABi_faidx_{wildcards.species}_{wildcards.assembly}_XXXXXX")"
        trap 'rm -rf "$TEMP_DIR"' EXIT
        cd "$TEMP_DIR"

        echo "[GABi] Copying assembly to temp directory..."
        cp "{input.fa}" local.fa

        echo "[GABi] Running samtools faidx..."
        samtools faidx local.fa

        echo "[GABi] Copying index back to output..."
        cp local.fa.fai "{output.fai}"

        echo "[GABi] Done"

        '''

rule E00_ref_dict_for_calling:
    input:
        fa=rules.D00_ref_fa_canonical.output.fa
    output:
        dict=REF_DICT_TMPL
    threads: cpu_func("samtools_index")
    resources:
        mem_mb=mem_func("samtools_index"),
        runtime=time_func("samtools_index")
    container: CONTAINERS["gatk"]
    log:
        os.path.join(config["OUT_FOLDER"], "GABi_results", "{species}", "{assembly}", "logs", "E00_ref_dict.log")
    shell:
        r'''
        set -euo pipefail
        mkdir -p $(dirname {log})

        exec > {log} 2>&1

        echo "[GABi] Creating sequence dictionary for {input.fa}"


        # TEMP DIRECTORY SETUP
        # -------------------------------------------------------------------
        WORK_DIR="$(gabi_get_workdir 25)"
        TEMP_DIR="$(mktemp -d "$WORK_DIR/GABi_dict_{wildcards.species}_{wildcards.assembly}_XXXXXX")"
        trap 'rm -rf "$TEMP_DIR"' EXIT
        cd "$TEMP_DIR"

        echo "[GABi] Copying assembly to temp directory..."
        cp "{input.fa}" local.fa

        echo "[GABi] Running gatk CreateSequenceDictionary..."
        gatk CreateSequenceDictionary -R local.fa -O local.dict

        echo "[GABi] Copying dictionary back to output..."
        cp local.dict "{output.dict}"

        echo "[GABi] Done"

        '''


# -------------------------------------------------------------------------------
#  PER ASM INTERVALS (no contig splitting)
# ===============================================================================

checkpoint E00_make_intervals:
    input:
        fa=rules.D00_ref_fa_canonical.output.fa,
        fai=rules.E00_ref_faidx_for_calling.output.fai,
        dict=rules.E00_ref_dict_for_calling.output.dict
    output:
        idir=directory(INTERVALS_DIR_TMPL)
    params:
        target_bp=int(50_000_000),
        # Declared as a param so Snakemake can SEE it. The shard layout depends
        # on asm_sex_chr, but that is read from samples_config inside run:,
        # which is invisible to the rerun triggers - so an existing intervals/
        # directory would be reused after asm_sex_chr changed, silently keeping
        # sex contigs bundled with autosomes. E01 would then find a mixed shard,
        # fall back to diploid, and call sex chromosomes at the wrong ploidy
        # with nothing in the logs to say so. As a param, any change to the
        # declared sex contigs re-shards.
        sex_chr=lambda w: ",".join(
            sorted(get_sex_chr(samples_config, w.species, w.assembly))
        ),
    threads: cpu_func("split_intervals")
    resources:
        mem_mb=mem_func("split_intervals"), runtime=time_func("split_intervals")
    container: CONTAINERS["gabi_base"]
    run:
        from pathlib import Path
        idir = Path(output.idir); idir.mkdir(parents=True, exist_ok=True)
        contigs = []
        with open(input.fai) as fh:
            for ln in fh:
                toks = ln.rstrip("\n").split("\t")
                if len(toks) >= 2: contigs.append((toks[0], int(toks[1])))

        #
        # GATK's --sample-ploidy applies to a whole HaplotypeCaller invocation,
        # but ploidy is a property of the sample (a hetero sample needs
        # ploidy 1 on its sex chromosomes) while intervals are a property of the
        # assembly. Segregating here is what lets E01 pick the ploidy per sample
        # per shard: sharding stays assembly-level, ploidy stays sample-level.
        #
        # One shard per sex contig rather than one shard for all of them,
        # because inference is per contig - a sample can come back haploid on Z
        # and ambiguous (hence diploid) on W, and a shared shard could only take
        # one ploidy for both.
        #
        # With no asm_sex_chr declared this is a no-op and the layout is exactly
        # as before. The other callers are unaffected in behaviour - bcftools
        # gets its ploidy from a file, not from shards - though they do see a
        # slightly different shard count.
        sex_chr = set(get_sex_chr(samples_config, wildcards.species, wildcards.assembly))
        sex_contigs = [(n, L) for n, L in contigs if n in sex_chr]
        contigs     = [(n, L) for n, L in contigs if n not in sex_chr]

        shard, current, acc, tgt = 1, [], 0, params.target_bp

        def _write_shard(entries, n_shard):
            # GATK reads the .list form (CHROM:FROM-TO, 1-based inclusive);
            # bcftools mpileup -R wants a regions file, for which BED is the
            # unambiguous choice (0-based, half-open). Contig names containing
            # ':' or '-' would be ambiguous in the .list form but are safe in
            # BED, which is another reason to emit both rather than convert.
            with open(idir / f"interval_{n_shard:03d}.list", "w") as f:
                for n, L in entries:
                    f.write(f"{n}:1-{L}\n")
            with open(idir / f"interval_{n_shard:03d}.bed", "w") as f:
                for n, L in entries:
                    f.write(f"{n}\t0\t{L}\n")

        def _flush():
            nonlocal shard, current, acc
            if not current: return
            _write_shard(current, shard)
            shard += 1; current.clear(); acc = 0

        for name, L in contigs:
            if L > tgt and not current:
                _write_shard([(name, L)], shard)
                shard += 1
            elif acc + L <= tgt:
                current.append((name, L)); acc += L
            else:
                _flush(); current.append((name, L)); acc = L
        _flush()

        # Sex contigs last, one shard each
        for name, L in sex_contigs:
            _write_shard([(name, L)], shard)
            shard += 1

        with open(idir / "MANIFEST.json", "w") as mf:
            json.dump({
                "intervals": sorted(str(p) for p in idir.glob("interval_*.list")),
                "sex_contigs": sorted(n for n, _ in sex_contigs),
            }, mf, indent=2)


# -------------------------------------------------------------------------------
#  HaplotypeCaller (BP_RESOLUTION), scattered
# ===============================================================================

rule E01_haplotypecaller_scattered:
    input:
        bam=lambda w: vc_get_calling_bam(w.species, w.assembly, w.sample_id),
        bai=lambda w: vc_get_calling_bam(w.species, w.assembly, w.sample_id) + ".bai",
        # Declared unconditionally, for every sample.
        selection=TECH_SELECTION_TMPL,
        ref=rules.D00_ref_fa_canonical.output.fa,
        fai=rules.E00_ref_faidx_for_calling.output.fai,
        dict=rules.E00_ref_dict_for_calling.output.dict,
        interval=lambda w: vc_interval_path(w.species, w.assembly, w.shard),
        # Needed by vc_shard_ploidy() when this sample's ploidy is inferred
        # rather than declared. [] for declared / no-sex-chr samples.
        sex_inference=_dv_sex_inference_input,
    output:
        gvcf=temp(GVCF_PART_TMPL),
        tbi=temp(GVCF_PART_TBI_TMPL)
    wildcard_constraints:
        shard=r"\d{3}"
    params:
        ploidy=vc_shard_ploidy,
        haploid_list=lambda w: ",".join(
            haploid_contigs_for_sample(w.species, w.assembly, w.sample_id)
        ),
        hc_params=vc_hc_params,
        bp_mode=_BP_MODE,
    threads: cpu_func("gatk_haplotypecaller")
    resources:
        mem_mb=mem_func("gatk_haplotypecaller"), runtime=time_func("gatk_haplotypecaller")
    container: CONTAINERS["gatk"]
    log:
        os.path.join(config["OUT_FOLDER"], "GABi_results", "{species}", "{assembly}", "samples", "{sample_id}", "logs",
                     "E01_hc.{sample_id}.{shard}.log")
    shell:
        r'''
        set -euo pipefail
        mkdir -p "$(dirname {output.gvcf})" "$(dirname {log})"

        exec > {log} 2>&1

        # Ensure BAM has a valid index (bai or csi); repair if empty/missing
        if [ -f "{input.bam}.bai" ] && [ ! -s "{input.bam}.bai" ]; then rm -f "{input.bam}.bai"; fi
        if [ ! -s "{input.bam}.bai" ] && [ ! -s "{input.bam}.csi" ]; then
            echo "[GABi] Re-indexing BAM before HC"
            samtools index -@ {threads} "{input.bam}" || true
        fi


        # TEMP DIRECTORY SETUP
        # Write output to local fast storage, then copy to final location
        # -------------------------------------------------------------------
        WORK_DIR="$(gabi_get_workdir 25)"
        TEMP_DIR="$(mktemp -d "$WORK_DIR/GABi_hc_{wildcards.species}_{wildcards.assembly}_{wildcards.sample_id}_{wildcards.shard}_XXXXXX")"
        trap 'rm -rf "$TEMP_DIR"' EXIT
        cd "$TEMP_DIR"
        echo "[GABi] Working directory: $TEMP_DIR"


        # JAVA/GATK MEMORY SETTINGS
        # -------------------------------------------------------------------
        # samjdk.buffer_size (bytes) based on available RAM
        #  < 8 GB  -> 256 KiB
        #  8–16 GB -> 512 KiB
        # 16–64 GB ->   1 MiB
        # ≥64 GB   ->   2 MiB
        MEM_MB={resources.mem_mb}
        XMX=$(( MEM_MB - 2048 )); [ "$XMX" -lt 1024 ] && XMX=1024
        if   [ "$MEM_MB" -ge 65536 ]; then BUF=2097152
        elif [ "$MEM_MB" -ge 16384 ]; then BUF=1048576
        elif [ "$MEM_MB" -ge  8192 ]; then BUF=524288
        else                               BUF=262144
        fi


        # READ-TYPE AWARE PARAMETERS
        # -------------------------------------------------------------------
        READ_CAT="{params.hc_params[read_category]}"
        PAIRHMM="{params.hc_params[pairhmm]}"
        PCR_MODEL="{params.hc_params[pcr_indel_model]}"
        MAX_HAP="{params.hc_params[max_haplotypes]}"
        MAX_READS="{params.hc_params[max_reads_start]}"
        MAX_ASM="{params.hc_params[max_asm_region]}"

        # Use native PairHMM threads only when a native implementation is selected
        NTH={threads}
        if [ "$PAIRHMM" = "LOGLESS_CACHING" ]; then NTH=1; fi

        # PLOIDY - verify rather than trust, as in E10/E30. params were resolved
        # from the interval file and sex_inference.tsv; if that ever stops being
        # re-evaluated, a sex chromosome called at the wrong ploidy is silent.
        SHARD_CONTIGS=$(cut -d: -f1 "{input.interval}" | LC_ALL=C sort -u | paste -sd, -)
        echo "[GABi] Shard contigs: $SHARD_CONTIGS"
        echo "[GABi] Haploid for this sample: {params.haploid_list:q}"
        echo "[GABi] Sample ploidy for this shard: {params.ploidy}"

        # A shard holding a haploid contig must be called at ploidy 1. E00 gives
        # every sex contig a shard of its own, so the only way to see one mixed
        # in with autosomes is a stale intervals/ directory built before
        # asm_sex_chr was declared. Diploid calls on a sex chromosome look
        # perfectly normal in the output, so refuse instead of proceeding.
        if [ -n "{params.haploid_list}" ] && [ "{params.ploidy}" != "1" ]; then
            for C in $(printf '%s' "{params.haploid_list}" | tr ',' ' '); do
                if cut -d: -f1 "{input.interval}" | grep -qxF "$C"; then
                    echo "[GABi] ERROR: shard {wildcards.shard} contains haploid contig '$C'" >&2
                    echo "[GABi]        but is being called at ploidy {params.ploidy}." >&2
                    echo "[GABi]        The interval layout is stale - it was built before" >&2
                    echo "[GABi]        asm_sex_chr was set. Re-run E00_make_intervals:" >&2
                    echo "[GABi]          snakemake --forcerun E00_make_intervals ..." >&2
                    echo "[GABi]        or delete $(dirname "{input.interval}")" >&2
                    exit 1
                fi
            done
        fi

        echo "[GABi] HaplotypeCaller for $READ_CAT reads (resolution={params.bp_mode})"
        echo "[GABi] Memory: XMX=${{XMX}}m, buffer=${{BUF}}, threads={threads}"
        echo "[GABi] Parameters: pairhmm=$PAIRHMM, pcr=$PCR_MODEL, haplotypes=$MAX_HAP, reads_start=$MAX_READS, asm_region=$MAX_ASM"

        # RESOLUTION-DEPENDENT FLAGS
        # -------------------------------------------------------------------
        # basepair: -ERC BP_RESOLUTION + --output-mode EMIT_ALL_CONFIDENT_SITES
        #           (emits a record for every single base pair)
        # block:    -ERC GVCF  (collapses non-variant regions into ref-blocks)
        #           --output-mode EMIT_ALL_CONFIDENT_SITES is incompatible with
        #           GVCF mode and must NOT be used
        if [ "{params.bp_mode}" = "True" ]; then
            ERC_FLAG="-ERC BP_RESOLUTION --output-mode EMIT_ALL_CONFIDENT_SITES"
            echo "[GABi] Resolution: basepair (BP_RESOLUTION)"
        else
            ERC_FLAG="-ERC GVCF"
            echo "[GABi] Resolution: block (GVCF)"
        fi

        # Run HaplotypeCaller - output to local temp
        gatk --java-options "-Xmx${{XMX}}m -XX:ParallelGCThreads=1 -Dsamjdk.compression_level=1 -Dsamjdk.buffer_size=${{BUF}} \
          -Dsamjdk.use_async_io_read_samtools=true -Dsamjdk.use_async_io_write_samtools=true" \
          HaplotypeCaller \
          -R "{input.ref}" \
          -I "{input.bam}" \
          -L "{input.interval}" \
          -O "output.g.vcf.gz" \
          $ERC_FLAG \
          --sample-ploidy {params.ploidy} \
          --native-pair-hmm-threads $NTH \
          --pair-hmm-implementation $PAIRHMM \
          --pcr-indel-model $PCR_MODEL \
          --max-num-haplotypes-in-population $MAX_HAP \
          --max-reads-per-alignment-start $MAX_READS \
          --max-assembly-region-size $MAX_ASM \
          --interval-padding 0 \
          --tmp-dir "$TEMP_DIR"

        gatk IndexFeatureFile -I "output.g.vcf.gz"

        # Copy results to final location
        echo "[GABi] Copying output to final location..."
        cp "output.g.vcf.gz" "{output.gvcf}"
        cp "output.g.vcf.gz.tbi" "{output.tbi}"

        echo "[GABi] Done"

        '''


# -------------------------------------------------------------------------------
#  GATHER BP gVCF per sample
# ===============================================================================

rule E02_gather_bp_gvcf:
    input:
        gvcfs=lambda w: [
            GVCF_PART_TMPL.format(
                species=w.species, assembly=w.assembly,
                sample_id=w.sample_id, shard=s
            )
            for s in vc_discover_shards(w.species, w.assembly)
        ],
        tbis=lambda w: [
            GVCF_PART_TBI_TMPL.format(
                species=w.species, assembly=w.assembly,
                sample_id=w.sample_id, shard=s
            )
            for s in vc_discover_shards(w.species, w.assembly)
        ]
    output:
        gvcf=temp(GATK_GVCF_TMPL),
        tbi=temp(GATK_GVCF_TBI_TMPL)
    threads: cpu_func("gatk_gathervcfs")
    resources:
        mem_mb=mem_func("gatk_gathervcfs"), runtime=time_func("gatk_gathervcfs")
    container: CONTAINERS["gatk"]
    log:
        os.path.join(
            config["OUT_FOLDER"], "GABi_results", "{species}", "{assembly}",
            "samples", "{sample_id}", "logs", "E02_gather_gvcf.{sample_id}.log"
        )
    shell:
        r'''
        set -euo pipefail
        mkdir -p "$(dirname {output.gvcf})" "$(dirname {log})"

        exec > {log} 2>&1


        # TEMP DIRECTORY SETUP
        # Write output to local fast storage, then copy to final location
        # -------------------------------------------------------------------
        WORK_DIR="$(gabi_get_workdir 50)"
        TEMP_DIR="$(mktemp -d "$WORK_DIR/GABi_gather_gvcf_{wildcards.species}_{wildcards.assembly}_{wildcards.sample_id}_XXXXXX")"
        trap 'rm -rf "$TEMP_DIR"' EXIT
        cd "$TEMP_DIR"
        echo "[GABi] Working directory: $TEMP_DIR"

        MEM_MB={resources.mem_mb}
        XMX=$(( MEM_MB - 2048 )); [ "$XMX" -lt 1024 ] && XMX=1024
        if   [ "$MEM_MB" -ge 65536 ]; then BUF=2097152
        elif [ "$MEM_MB" -ge 16384 ]; then BUF=1048576
        elif [ "$MEM_MB" -ge  8192 ]; then BUF=524288
        else                               BUF=262144
        fi

        # Count the actual number of gVCF files (not including .tbi)
        set -- {input.gvcfs}
        N=$#
        echo "[GABi] gather_bp_gvcf: N=$N gVCFs to gather"

        if [ "$N" -eq 1 ]; then
            # Fast path: single shard -> just copy + reindex
            SRC="$1"
            echo "[GABi] Single gVCF, using fast path (copy)"
            cp "$SRC" "output.g.vcf.gz"
            gatk IndexFeatureFile -I "output.g.vcf.gz"
        else
            echo "[GABi] Multiple gVCFs, using GatherVcfs"
            INLIST=""
            for f in {input.gvcfs}; do
                INLIST="$INLIST -I $f"
            done

            gatk --java-options "-Xmx${{XMX}}m -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=${{BUF}} \
                                 -Dsamjdk.compression_level=2 \
                                 -Dsamjdk.use_async_io_read_samtools=true \
                                 -Dsamjdk.use_async_io_write_samtools=true" \
                 GatherVcfs $INLIST -O "output.g.vcf.gz" --TMP_DIR "$TEMP_DIR"

            gatk IndexFeatureFile -I "output.g.vcf.gz"
        fi

        # Copy results to final location
        echo "[GABi] Copying output to final location..."
        cp "output.g.vcf.gz" "{output.gvcf}"
        cp "output.g.vcf.gz.tbi" "{output.tbi}"

        echo "[GABi] Done"

        '''

# -------------------------------------------------------------------------------
#  SCATTERED SINGLE-SAMPLE GENOTYPING + FINAL BCF
# ===============================================================================

rule E03_genotype_scattered:
    input:
        ref=rules.D00_ref_fa_canonical.output.fa,
        fai=rules.E00_ref_faidx_for_calling.output.fai,
        dict=rules.E00_ref_dict_for_calling.output.dict,
        gvcf=GVCF_PART_TMPL,
        tbi=GVCF_PART_TBI_TMPL,
        interval=lambda w: vc_interval_path(w.species, w.assembly, w.shard),
        sex_inference=_dv_sex_inference_input,
    output:
        vcf=temp(GENO_PART_TMPL),
        tbi=temp(GENO_PART_TBI_TMPL)
    wildcard_constraints:
        shard=r"\d{3}"
    threads: cpu_func("gatk_genotypegvcfs")
    resources:
        mem_mb=mem_func("gatk_genotypegvcfs"), runtime=time_func("gatk_genotypegvcfs")
    container: CONTAINERS["gatk"]
    params:
        bp_mode=_BP_MODE,
        ploidy=vc_shard_ploidy,
    log:
        os.path.join(config["OUT_FOLDER"], "GABi_results", "{species}", "{assembly}", "samples", "{sample_id}", "logs",
                     "E03_genotype.{sample_id}.{shard}.log")
    shell:
        r'''
        set -euo pipefail
        mkdir -p "$(dirname {output.vcf})" "$(dirname {log})"

        exec > {log} 2>&1


        # TEMP DIRECTORY SETUP
        # Write output to local fast storage, then copy to final location
        # -------------------------------------------------------------------
        WORK_DIR="$(gabi_get_workdir 20)"
        TEMP_DIR="$(mktemp -d "$WORK_DIR/GABi_genotype_{wildcards.species}_{wildcards.assembly}_{wildcards.sample_id}_{wildcards.shard}_XXXXXX")"
        trap 'rm -rf "$TEMP_DIR"' EXIT
        cd "$TEMP_DIR"
        echo "[GABi] Working directory: $TEMP_DIR"

        MEM_MB={resources.mem_mb}
        XMX=$(( MEM_MB - 2048 )); [ "$XMX" -lt 1024 ] && XMX=1024
        if   [ "$MEM_MB" -ge 65536 ]; then BUF=2097152
        elif [ "$MEM_MB" -ge 16384 ]; then BUF=1048576
        elif [ "$MEM_MB" -ge  8192 ]; then BUF=524288
        else                               BUF=262144
        fi

        # RESOLUTION-DEPENDENT FLAGS
        # -------------------------------------------------------------------
        # basepair: emit every site individually (threshold 0 + non-variant sites)
        # block:    emit variant sites + reference blocks (default QUAL threshold)
        if [ "{params.bp_mode}" = "True" ]; then
            EXTRA_FLAGS="--include-non-variant-sites --standard-min-confidence-threshold-for-calling 0"
            echo "[GABi] Resolution: basepair (all sites emitted individually)"
        else
            EXTRA_FLAGS="--include-non-variant-sites"
            echo "[GABi] Resolution: block (variant sites + reference blocks)"
        fi

        echo "[GABi] Sample ploidy for this shard: {params.ploidy}"

        # Run GenotypeGVCFs - output to local temp
        gatk --java-options "-Xmx${{XMX}}m -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=${{BUF}} \
          -Dsamjdk.compression_level=1 -Dsamjdk.use_async_io_read_samtools=true -Dsamjdk.use_async_io_write_samtools=true" \
          GenotypeGVCFs -R "{input.ref}" -V "{input.gvcf}" -L "{input.interval}" \
          --sample-ploidy {params.ploidy} \
          --interval-padding 0 \
          -O "output.vcf.gz" \
          $EXTRA_FLAGS \
          --tmp-dir "$TEMP_DIR"

        gatk IndexFeatureFile -I "output.vcf.gz"

        # Copy results to final location
        echo "[GABi] Copying output to final location..."
        cp "output.vcf.gz" "{output.vcf}"
        cp "output.vcf.gz.tbi" "{output.tbi}"

        echo "[GABi] Done"

        '''

rule E04_gather_genotyped_to_bcf:
    input:
        vcfs=lambda w: [ GENO_PART_TMPL.format(species=w.species, assembly=w.assembly, sample_id=w.sample_id, shard=s)
                    for s in vc_discover_shards(w.species, w.assembly) ],
        tbis=lambda w: [ GENO_PART_TBI_TMPL.format(species=w.species, assembly=w.assembly, sample_id=w.sample_id, shard=s)
                    for s in vc_discover_shards(w.species, w.assembly) ]
    output:
        bcf=temp(GATK_BCF_TMPL),
        csi=temp(GATK_BCF_CSI_TMPL)
    threads: cpu_func("bcftools_concat")
    resources:
        mem_mb=mem_func("bcftools_concat"), runtime=time_func("bcftools_concat")
    container: CONTAINERS["gabi_base"]
    log:
        os.path.join(config["OUT_FOLDER"], "GABi_results", "{species}", "{assembly}", "samples", "{sample_id}", "logs",
                     "E04_gather_bcf.{sample_id}.log")
    shell:
        r'''
        set -euo pipefail
        mkdir -p "$(dirname {output.bcf})" "$(dirname {log})"

        exec > {log} 2>&1


        # TEMP DIRECTORY SETUP
        # Write output to local fast storage, then copy to final location
        # -------------------------------------------------------------------
        WORK_DIR="$(gabi_get_workdir 50)"
        TEMP_DIR="$(mktemp -d "$WORK_DIR/GABi_gather_bcf_{wildcards.species}_{wildcards.assembly}_{wildcards.sample_id}_XXXXXX")"
        trap 'rm -rf "$TEMP_DIR"' EXIT
        cd "$TEMP_DIR"
        echo "[GABi] Working directory: $TEMP_DIR"

        # How many genotyped VCFs did we get?
        set -- {input.vcfs}
        N=$#
        echo "[GABi] gather_genotyped_to_bcf: inputs=$N, threads={threads}"

        if [ "$N" -eq 1 ]; then
            SRC="$1"
            echo "[GABi] Single input, converting to BCF"
            # single input -> just convert to BCF
            export HTS_COMPRESSION_LEVEL=1
            bcftools view -O b -o "output.bcf" "$SRC" --threads {threads}
        else
            echo "[GABi] Multiple inputs, indexing and concatenating"
            # Index all input files first (in case they're not indexed)
            for vcf in {input.vcfs}; do
                if [ ! -f "${{vcf}}.tbi" ] && [ ! -f "${{vcf}}.csi" ]; then
                    echo "[GABi] Indexing $vcf"
                    bcftools index -t "$vcf"
                fi
            done

            # multiple inputs -> concatenate to BCF
            export HTS_COMPRESSION_LEVEL=1
            bcftools concat -a -D -O b -o "output.bcf" {input.vcfs} --threads {threads}
        fi

        echo "[GABi] Creating CSI index for BCF"
        bcftools index -f --csi --threads {threads} "output.bcf"

        # Copy results to final location
        echo "[GABi] Copying output to final location..."
        cp "output.bcf" "{output.bcf}"
        cp "output.bcf.csi" "{output.csi}"

        echo "[GABi] Done"

        # Clean up empty tmp directory if all temp files are gone
        TMP_DIR="$(dirname {input.vcfs[0]})"
        if [ -d "$TMP_DIR" ] && [ -z "$(ls -A "$TMP_DIR" 2>/dev/null)" ]; then
            rmdir "$TMP_DIR" 2>/dev/null || true
            echo "[GABi] Removed empty tmp directory"
        fi

        '''
# ===============================================================================
#  DeepVariant PATH (single rule replaces the GATK scatter-gather chain)
# ===============================================================================

rule E10_deepvariant:
    """Run DeepVariant as three separate binaries with per-step env.

    Why split from run_deepvariant: the three steps need different thread
    settings and run_deepvariant sets env once for all of them.
      make_examples  N parallel processes, each pinned to 1 thread, to avoid
                     RLIMIT_NPROC / cgroup pids.max blowup on clusters
      call_variants  1 process, capped to `threads` (see the note at that step -
                     it previously ran uncapped, which is only safe under a
                     scheduler that enforces the allocation)
      postprocess    multi-CPU but light, --cpus 4
    """
    input:
        ref=rules.D00_ref_fa_canonical.output.fa,
        fai=rules.E00_ref_faidx_for_calling.output.fai,
        bam=lambda w: _dv_get_bam(w.species, w.assembly, w.sample_id),
        bai=lambda w: _dv_get_bam(w.species, w.assembly, w.sample_id) + ".bai",
        selection=TECH_SELECTION_TMPL,
        sex_inference=_dv_sex_inference_input,
    output:
        gvcf=temp(DV_GVCF_TMPL),
        tbi=temp(DV_GVCF_TBI_TMPL),
    params:
        model_type=lambda w: _dv_model_type(w.species, w.assembly, w.sample_id),
        model_dir=lambda w: {
            "WGS": "/opt/models/wgs",
            "PACBIO": "/opt/models/pacbio",
            "ONT_R104": "/opt/models/ont_r104",
        }[_dv_model_type(w.species, w.assembly, w.sample_id)],
        haploid_contigs=_dv_haploid_contigs_flag,
    threads: cpu_func("deepvariant")
    resources:
        mem_mb=mem_func("deepvariant"),
        runtime=time_func("deepvariant")
    container: CONTAINERS["deepvariant"]
    log:
        os.path.join(
            config["OUT_FOLDER"], "GABi_results", "{species}", "{assembly}",
            "samples", "{sample_id}", "logs", "E10_deepvariant.{sample_id}.log"
        )
    shell:
        r'''
        set -euo pipefail
        mkdir -p "$(dirname {output.gvcf})" "$(dirname {log})"
        exec > {log} 2>&1

        # Raise thread limit (best-effort, applies to all child processes)
        ulimit -u unlimited 2>/dev/null || ulimit -u "$(ulimit -Hu)" 2>/dev/null || true

        echo "[GABi] DeepVariant (split-steps) for {wildcards.sample_id}"

        # The model type comes from a params function that read
        # tech_selection.json. Snakemake re-evaluates params once a job's inputs
        # exist, so it should already reflect the selected technology - but this
        # rule runs for days, so verify rather than assume. If Snakemake ever
        # stops refreshing params, this turns a silently wrong model into a
        # loud failure.
        SELECTED_TECH=$(python3 -c "import json,sys; print(json.load(open(sys.argv[1]))['chosen'])" "{input.selection}")
        case "$SELECTED_TECH" in
            hifi) EXPECTED_MODEL=PACBIO ;;
            ont)  EXPECTED_MODEL=ONT_R104 ;;
            *)    EXPECTED_MODEL=WGS ;;
        esac
        if [ "$EXPECTED_MODEL" != "{params.model_type}" ]; then
            echo "[GABi] ERROR: model type does not match the selected technology." >&2
            echo "[GABi]        tech_selection.json says '$SELECTED_TECH' (-> $EXPECTED_MODEL)" >&2
            echo "[GABi]        but this job was configured with {params.model_type}" >&2
            exit 1
        fi
        echo "[GABi] Selected technology: $SELECTED_TECH"
        echo "[GABi] Model type: {params.model_type} ({params.model_dir})"
        # Validate model type — this rule has been tested only for these.
        case "{params.model_type}" in
            WGS|PACBIO|ONT_R104)
                echo "[GABi] Model type validated: {params.model_type}"
                ;;
            *)
                echo "[GABi] ERROR: Unsupported model type '{params.model_type}'"
                echo "[GABi] This rule supports: WGS, PACBIO, ONT_R104"
                echo "[GABi] WES, HYBRID_PACBIO_ILLUMINA, and others have not been validated."
                exit 1
                ;;
        esac
        echo "[GABi] Threads/shards: {threads}"
        echo "[GABi] ulimit -u (soft/hard): $(ulimit -Su) / $(ulimit -Hu)"

        WORK_DIR="$(gabi_get_workdir 150)"
        TEMP_DIR="$(mktemp -d "$WORK_DIR/GABi_deepvariant_{wildcards.species}_{wildcards.assembly}_{wildcards.sample_id}_XXXXXX")"
        trap 'rm -rf "$TEMP_DIR"' EXIT
        echo "[GABi] Working directory: $TEMP_DIR"

        DV_TMP="$TEMP_DIR/dv_tmp"
        mkdir -p "$DV_TMP"

        EXAMPLES="$DV_TMP/make_examples.tfrecord@{threads}.gz"
        GVCF_RECORDS="$DV_TMP/gvcf.tfrecord@{threads}.gz"
        CV_OUTPUT="$DV_TMP/call_variants_output.tfrecord.gz"
        MODEL_JSON="{params.model_dir}/model.example_info.json"

        # Shared JAX compile cache (helps small-model in make_examples
        # AND reused in call_variants).
        export JAX_COMPILATION_CACHE_DIR="$TEMP_DIR/jax_cache"
        mkdir -p "$JAX_COMPILATION_CACHE_DIR"

        # ==========================================================
        # Mirror run_deepvariant.py decisions that depend on the
        # model.example_info.json and on the model type.
        # ==========================================================

        # Small model: only enable if the JSON declares one
        # (matches run_deepvariant.py:_use_small_model).
        SMALL_MODEL_ME_FLAG=""
        SMALL_MODEL_PP_FLAG=""
        if python3 - "$MODEL_JSON" <<'PYEOF' 2>/dev/null
import json, sys
try:
    with open(sys.argv[1]) as f:
        info = json.load(f)
    flags = info.get('flags_for_calling') or dict()
    sys.exit(0 if flags.get('trained_small_model_path') else 1)
except Exception:
    sys.exit(1)
PYEOF
        then
            SMALL_MODEL_ME_FLAG="--call_small_model_examples"
            SMALL_MODEL_PP_FLAG="--small_model_cvo_records $DV_TMP/make_examples_call_variant_outputs.tfrecord@{threads}.gz"
            echo "[GABi] Small model: ENABLED (declared in model.example_info.json)"
        else
            echo "[GABi] Small model: disabled (not declared in model.example_info.json)"
        fi

        # Phasing: enable for long-read models
        # (matches run_deepvariant.py:_should_phase_vcf).
        PHASE_ME_FLAGS=""
        PHASE_PP_FLAGS=""
        case "{params.model_type}" in
            PACBIO|ONT_R104)
                PHASE_TSV="$DV_TMP/read-phasing_debug@{threads}.tsv"
                PHASE_ME_FLAGS="--output_local_read_phasing $PHASE_TSV --output_phase_info"
                PHASE_PP_FLAGS="--phased_reads_input_path $PHASE_TSV"
                echo "[GABi] Phasing: ENABLED for {params.model_type}"
                ;;
            *)
                echo "[GABi] Phasing: disabled for {params.model_type}"
                ;;
        esac

        export HAPLOID_FLAG="{params.haploid_contigs}"
        if [ -n "$HAPLOID_FLAG" ]; then
            echo "[GABi] Haploid contigs: $HAPLOID_FLAG"
        else
            echo "[GABi] Haploid contigs: (none — all sites called diploid)"
        fi

        # ==========================================================
        # STEP 1: make_examples  --  THREAD-PINNED, PARALLEL
        # ==========================================================
        echo "[GABi] ===== make_examples starting ====="
        (
            export OPENBLAS_NUM_THREADS=1
            export MKL_NUM_THREADS=1
            export OMP_NUM_THREADS=1
            export NUMEXPR_NUM_THREADS=1
            export TF_NUM_INTEROP_THREADS=1
            export TF_NUM_INTRAOP_THREADS=1
            export TF_ENABLE_ONEDNN_OPTS=0
            export TF_CPP_MIN_LOG_LEVEL=2

            time seq 0 $(( {threads} - 1 )) | parallel -q -j {threads} --halt now,fail=1 --line-buffer \
                /opt/deepvariant/bin/make_examples \
                    --mode calling \
                    --ref "{input.ref}" \
                    --reads "{input.bam}" \
                    --examples "$EXAMPLES" \
                    --gvcf "$GVCF_RECORDS" \
                    --checkpoint "{params.model_dir}" \
                    --checkpoint_json "$MODEL_JSON" \
                    $SMALL_MODEL_ME_FLAG \
                    $PHASE_ME_FLAGS \
                    ${{HAPLOID_FLAG:-}} \
                    --sample_name "{wildcards.sample_id}" \
                    --task {{}}
        )
        echo "[GABi] ===== make_examples done ====="

        # ==========================================================
        # STEP 2: call_variants  --  SINGLE PROCESS, THREAD-CAPPED
        # ==========================================================
        # This step used to run with no thread limits, on the reasoning that a
        # single inference process should be free to use the whole machine. That
        # holds under SLURM, where the cgroup clamps it to the allocation. Run
        # locally on a shared server there is nothing to clamp it: measured at
        # 12534% CPU (~125 of 160 cores) for a job Snakemake was told uses 18.

        echo "[GABi] ===== call_variants starting ====="
        (
            export TF_CPP_MIN_LOG_LEVEL=2

            export OMP_NUM_THREADS={threads}
            export OPENBLAS_NUM_THREADS={threads}
            export MKL_NUM_THREADS={threads}
            export NUMEXPR_NUM_THREADS={threads}
            export TF_NUM_INTRAOP_THREADS={threads}
            export TF_NUM_INTEROP_THREADS=2

            echo "[GABi] Thread cap: {threads} (intra-op), 2 (inter-op)"
            time /opt/deepvariant/bin/call_variants \
                --outfile "$CV_OUTPUT" \
                --examples "$EXAMPLES" \
                --checkpoint "{params.model_dir}"
        )
        echo "[GABi] ===== call_variants done ====="
        echo "[GABi] Check the cap took: user/real above should be ~{threads},"
        echo "[GABi] not the machine's core count. If it is still far higher,"
        echo "[GABi] the TF_* variables are being ignored (JAX backend)."

        # ==========================================================
        # STEP 3: postprocess_variants  --  MULTI-CPU, LIGHTWEIGHT
        # ==========================================================
        echo "[GABi] ===== postprocess_variants starting ====="
        (
            export TF_CPP_MIN_LOG_LEVEL=2
            time /opt/deepvariant/bin/postprocess_variants \
                --ref "{input.ref}" \
                --infile "$CV_OUTPUT" \
                --outfile "$TEMP_DIR/output.vcf.gz" \
                --nonvariant_site_tfrecord_path "$GVCF_RECORDS" \
                --gvcf_outfile "$TEMP_DIR/output.g.vcf.gz" \
                --checkpoint_json "$MODEL_JSON" \
                $SMALL_MODEL_PP_FLAG \
                $PHASE_PP_FLAGS \
                ${{HAPLOID_FLAG:-}} \
                --cpus 4 \
                --sample_name "{wildcards.sample_id}"
        )
        echo "[GABi] ===== postprocess_variants done ====="

        echo "[GABi] Copying gVCF to output..."
        cp "$TEMP_DIR/output.g.vcf.gz"     "{output.gvcf}"
        cp "$TEMP_DIR/output.g.vcf.gz.tbi" "{output.tbi}"
        echo "[GABi] Done"
        '''


rule E11_dv_to_bcf:
    """Convert DeepVariant gVCF to BCF.
    block mode:    gVCF -> BCF  (reference blocks preserved)
    basepair mode: gVCF -> expand ref-blocks to per-site records -> BCF
    """
    input:
        gvcf=DV_GVCF_TMPL,
        tbi=DV_GVCF_TBI_TMPL,
        ref=rules.D00_ref_fa_canonical.output.fa,
        fai=rules.E00_ref_faidx_for_calling.output.fai,
    output:
        bcf=temp(DV_BCF_TMPL),
        csi=temp(DV_BCF_CSI_TMPL),
    params:
        bp_mode=_BP_MODE,
    threads: cpu_func("bcftools_concat")
    resources:
        mem_mb=mem_func("bcftools_concat"),
        runtime=time_func("bcftools_concat")
    container: CONTAINERS["gabi_base"]
    log:
        os.path.join(
            config["OUT_FOLDER"], "GABi_results", "{species}", "{assembly}",
            "samples", "{sample_id}", "logs", "E11_dv_to_bcf.{sample_id}.log"
        )
    shell:
        r'''
        set -euo pipefail
        mkdir -p "$(dirname {output.bcf})" "$(dirname {log})"

        exec > {log} 2>&1

        echo "[GABi] Converting DeepVariant gVCF to BCF for {wildcards.sample_id}"


        # TEMP DIRECTORY SETUP
        # -------------------------------------------------------------------
        WORK_DIR="$(gabi_get_workdir 50)"
        TEMP_DIR="$(mktemp -d "$WORK_DIR/GABi_dv2bcf_{wildcards.species}_{wildcards.assembly}_{wildcards.sample_id}_XXXXXX")"
        trap 'rm -rf "$TEMP_DIR"' EXIT
        cd "$TEMP_DIR"
        echo "[GABi] Working directory: $TEMP_DIR"

        export HTS_COMPRESSION_LEVEL=1

        if [ "{params.bp_mode}" = "True" ]; then
            echo "[GABi] Resolution: basepair (expanding reference blocks)"
            # REF case: `convert --gvcf2vcf` copies the reference base verbatim
            # from the FASTA, so on a soft-masked assembly every expanded record
            # inside a repeat comes back lowercase, while variant records (whose
            # REF DeepVariant already uppercased) stay uppercase - two
            # conventions in one file, split by whether a position was inside a
            # reference block. The awk pass normalises field 4 only; VCF is
            # tab-delimited, so no other field can be touched. Same as E31.
            bcftools convert --gvcf2vcf "{input.gvcf}" \
                --fasta-ref "{input.ref}" -Ov \
              | awk 'BEGIN{{OFS="\t"}} /^#/{{print; next}} {{$4=toupper($4); print}}' \
              | bcftools view -A -O b -o "output.bcf" --threads {threads}
        else
            echo "[GABi] Resolution: block (preserving reference blocks)"
            bcftools view -A -O b -o "output.bcf" "{input.gvcf}" --threads {threads}
        fi

        echo "[GABi] Creating CSI index for BCF"
        bcftools index -f --csi --threads {threads} "output.bcf"

        # Copy results to final location
        echo "[GABi] Copying output to final location..."
        cp "output.bcf"     "{output.bcf}"
        cp "output.bcf.csi" "{output.csi}"

        echo "[GABi] Done"

        '''


# ===============================================================================
#  bcftools PATH  (mpileup | call, scattered over the same shards as GATK)
# ===============================================================================

rule E29_bcftools_ploidy:
    """Write the per-sample ploidy and sample-sex files for `bcftools call`.

    bcftools expresses ploidy as (CHROM, FROM, TO, SEX, PLOIDY) rows, with each
    sample assigned a SEX label via -S. GABi uses the SAMPLE ID ITSELF as that
    label rather than a karyotype such as XY/ZW. Two reasons:

      - Sex is usually not declared. sample_sex_karyo is empty for most data,
        and D06 infers ploidy per contig from coverage instead - a per-sample
        answer that no fixed karyotype vocabulary can express.
      - It makes the single-sample and cohort cases the same mechanism. A joint
        ploidy file is just these per-sample rows concatenated.

    Coordinates are written explicitly (1..contig_length from the .fai) rather
    than as '*', which is only special-cased for the default row.
    """
    input:
        fai=rules.E00_ref_faidx_for_calling.output.fai,
        # Present only when this sample's ploidy has to be inferred; returns []
        # for declared / no-sex-contig samples so D06 is not pulled in for them.
        sex_inference=_dv_sex_inference_input,
    output:
        ploidy=BCFT_PLOIDY_TMPL,
        samples=BCFT_SAMPLES_TMPL,
    threads: 1
    resources:
        mem_mb=mem_func("split_intervals"),
        runtime=time_func("split_intervals")
    container: CONTAINERS["gabi_base"]
    log:
        os.path.join(
            config["OUT_FOLDER"], "GABi_results", "{species}", "{assembly}",
            "samples", "{sample_id}", "logs", "E29_bcftools_ploidy.{sample_id}.log"
        )
    run:
        from pathlib import Path

        sid = wildcards.sample_id
        Path(output.ploidy).parent.mkdir(parents=True, exist_ok=True)
        Path(str(log)).parent.mkdir(parents=True, exist_ok=True)

        lengths = {}
        with open(input.fai) as fh:
            for ln in fh:
                toks = ln.rstrip("\n").split("\t")
                if len(toks) >= 2:
                    lengths[toks[0]] = int(toks[1])

        haploid = _bcft_haploid_contigs(wildcards)
        missing = [c for c in haploid if c not in lengths]
        if missing:
            raise ValueError(
                f"[GABi] {sid}: sex contigs not present in the reference index: "
                f"{', '.join(missing)}. Check asm_sex_chr against the assembly."
            )

        with open(output.ploidy, "w") as f:
            # Default row: everything not listed below is diploid.
            f.write(f"*\t*\t*\t{sid}\t2\n")
            for c in haploid:
                f.write(f"{c}\t1\t{lengths[c]}\t{sid}\t1\n")

        with open(output.samples, "w") as f:
            f.write(f"{sid}\t{sid}\n")

        with open(str(log), "w") as f:
            f.write(f"[GABi] Ploidy file for {sid}\n")
            f.write(f"[GABi] Default ploidy: 2\n")
            if haploid:
                f.write(f"[GABi] Haploid contigs: {', '.join(haploid)}\n")
            else:
                f.write("[GABi] Haploid contigs: (none - all contigs diploid)\n")


rule E30_bcftools_call_scattered:
    """bcftools mpileup | call for one shard of one sample."""
    input:
        bam=lambda w: vc_get_calling_bam(w.species, w.assembly, w.sample_id),
        bai=lambda w: vc_get_calling_bam(w.species, w.assembly, w.sample_id) + ".bai",
        selection=TECH_SELECTION_TMPL,
        ref=rules.D00_ref_fa_canonical.output.fa,
        fai=rules.E00_ref_faidx_for_calling.output.fai,
        bed=lambda w: os.path.join(
            vc_intervals_dir(w.species, w.assembly), f"interval_{w.shard}.bed"
        ),
        ploidy=rules.E29_bcftools_ploidy.output.ploidy,
        samples=rules.E29_bcftools_ploidy.output.samples,
    output:
        bcf=temp(BCFT_PART_TMPL),
        csi=temp(BCFT_PART_CSI_TMPL),
    wildcard_constraints:
        shard=r"\d{3}"
    params:
        profile=lambda w: _bcft_profile_for_sample(w.species, w.assembly, w.sample_id),
        bins=_BCFT_GVCF_BINS,
        bp_mode=_BP_MODE,
    threads: cpu_func("gatk_haplotypecaller")   # TODO: add a bcftools_call key to resources.yaml
    resources:
        mem_mb=mem_func("gatk_haplotypecaller"),
        runtime=time_func("gatk_haplotypecaller")
    container: CONTAINERS["gabi_base"]
    log:
        os.path.join(
            config["OUT_FOLDER"], "GABi_results", "{species}", "{assembly}",
            "samples", "{sample_id}", "logs", "E30_bcftools_call.{sample_id}.{shard}.log"
        )
    shell:
        r'''
        set -euo pipefail
        mkdir -p "$(dirname {output.bcf})" "$(dirname {log})"
        exec > {log} 2>&1

        # Same guard as E10: params were resolved from tech_selection.json, and
        # should already reflect the selected technology, but verify rather than
        # trust - an mpileup profile mismatched to the data is silent, not fatal.
        SELECTED_TECH=$(python3 -c "import json,sys; print(json.load(open(sys.argv[1]))['chosen'])" "{input.selection}")
        case "$SELECTED_TECH" in
            hifi)         EXPECTED_PROFILE=pacbio-ccs-1.20 ;;
            ont)          EXPECTED_PROFILE=ont-sup-1.20 ;;
            illumina|10x) EXPECTED_PROFILE=illumina-1.20 ;;
            *)            EXPECTED_PROFILE=illumina-1.20 ;;
        esac
        if [ "$EXPECTED_PROFILE" != "{params.profile}" ]; then
            echo "[GABi] ERROR: mpileup profile does not match the selected technology." >&2
            echo "[GABi]        tech_selection.json says '$SELECTED_TECH' (-> $EXPECTED_PROFILE)" >&2
            echo "[GABi]        but this job was configured with {params.profile}" >&2
            exit 1
        fi

        echo "[GABi] bcftools calling {wildcards.sample_id} shard {wildcards.shard}"
        echo "[GABi] Selected technology: $SELECTED_TECH  |  profile: {params.profile}"
        echo "[GABi] Ploidy file:"
        sed 's/^/[GABi]   /' "{input.ploidy}"

        # RESOLUTION
        # -------------------------------------------------------------------
        # block:    --gvcf collapses runs of hom-ref into blocks binned by DP
        # basepair: omit --gvcf (and never pass -v) so `call` emits one record
        #           per position - genuinely per-site, unlike expanding blocks
        if [ "{params.bp_mode}" = "True" ]; then
            GVCF_FLAG=""
            echo "[GABi] Resolution: basepair (one record per position)"
        else
            GVCF_FLAG="--gvcf {params.bins}"
            echo "[GABi] Resolution: block (gVCF, DP bins {params.bins})"
        fi

        # ANNOTATIONS - both of these are load-bearing
        # -------------------------------------------------------------------
        # -a FORMAT/DP : `call --gvcf` bins reference blocks by per-sample DP,
        #                and DP is NOT in mpileup's default annotation set.
        #                Verified on bcftools 1.23: omitting it aborts with
        #                "--gvcf output mode requires FORMAT/DP tag, which is
        #                not present in the input header". gvcf.c also has a
        #                silent fallback (a missing DP flushes every site as its
        #                own record), so on a build without that up-front header
        #                check the result would be a basepair-sized file wearing
        #                a .g.vcf name instead of an error. Required either way.
        # -a FORMAT/AD : needed by F01 for allele-balance tags. On by default,
        #                requested explicitly so it survives future changes.
        # call -a GQ   : FORMAT/GQ is an optional `call` annotation. Without it
        #                MIN_GQ tagging is inert and F01's NO_GQ filter fires on
        #                every variant record.
        bcftools mpileup \
            -X {params.profile} \
            -f "{input.ref}" \
            -a FORMAT/AD,FORMAT/DP \
            -R "{input.bed}" \
            -Ou "{input.bam}" \
        | bcftools call \
            -m \
            -a GQ \
            $GVCF_FLAG \
            --ploidy-file "{input.ploidy}" \
            -S "{input.samples}" \
            -Ob -o "{output.bcf}" \
            --threads {threads}

        bcftools index -f --csi --threads {threads} "{output.bcf}"

        echo "[GABi] Done"
        '''


rule E31_bcftools_gather:
    """Concatenate per-shard bcftools BCFs, fill uncovered positions, emit the
    per-sample gVCF and raw BCF.

    `bcftools mpileup` emits a record only where at least one read is present,
    so positions with no coverage are absent from the call set entirely.
    DeepVariant and GATK gVCFs instead span the whole reference, representing
    those regions as blocks with DP=0. fill_uncovered.py closes that difference
    so all three callers describe the same genome - which matters for anything
    downstream that walks positions rather than records, and for F01, which can
    only tag NO_DATA on a record that exists.
    """
    input:
        bcfs=lambda w: [
            BCFT_PART_TMPL.format(
                species=w.species, assembly=w.assembly,
                sample_id=w.sample_id, shard=s
            )
            for s in vc_discover_shards(w.species, w.assembly)
        ],
        csis=lambda w: [
            BCFT_PART_CSI_TMPL.format(
                species=w.species, assembly=w.assembly,
                sample_id=w.sample_id, shard=s
            )
            for s in vc_discover_shards(w.species, w.assembly)
        ],
        ref=rules.D00_ref_fa_canonical.output.fa,
        fai=rules.E00_ref_faidx_for_calling.output.fai,
    output:
        gvcf=temp(BCFT_GVCF_TMPL),
        tbi=temp(BCFT_GVCF_TBI_TMPL),
        bcf=temp(BCFT_BCF_TMPL),
        csi=temp(BCFT_BCF_CSI_TMPL),
    params:
        fill_script=str(scripts_dir / "fill_uncovered.py"),
        bp_mode=_BP_MODE,
    threads: cpu_func("bcftools_concat")
    resources:
        mem_mb=mem_func("bcftools_concat"),
        runtime=time_func("bcftools_concat")
    container: CONTAINERS["gabi_base"]
    log:
        os.path.join(
            config["OUT_FOLDER"], "GABi_results", "{species}", "{assembly}",
            "samples", "{sample_id}", "logs", "E31_bcftools_gather.{sample_id}.log"
        )
    shell:
        r'''
        set -euo pipefail
        mkdir -p "$(dirname {output.gvcf})" "$(dirname {log})"
        exec > {log} 2>&1

        WORK_DIR="$(gabi_get_workdir 50)"
        TEMP_DIR="$(mktemp -d "$WORK_DIR/GABi_bcftgather_{wildcards.species}_{wildcards.assembly}_{wildcards.sample_id}_XXXXXX")"
        trap 'rm -rf "$TEMP_DIR"' EXIT
        cd "$TEMP_DIR"

        export HTS_COMPRESSION_LEVEL=1

        set -- {input.bcfs}
        echo "[GABi] Concatenating $# shard(s) for {wildcards.sample_id}"

        # Shards are disjoint whole contigs, so this is a straight concatenation;
        # -a -D matches the GATK gather and is harmless insurance against any
        # future change to how intervals are built.
        bcftools concat -a -D -O b -o "concat.bcf" {input.bcfs} --threads {threads}

        # ONE RECORD PER POSITION (basepair mode only)
        # -------------------------------------------------------------------
        # mpileup emits an indel CANDIDATE alongside the substitution record at
        # the same position. When `call` genotypes that candidate as hom-ref you
        # get two records for one position, both with ALT=".", e.g.
        #     364  T  .  DP=8;...
        #     364  t  .  INDEL;IDV=2;...
        # DeepVariant and GATK emit exactly one record per position, and
        # duplicates would be double-counted by anything walking positions
        # (heterozygosity windows, missingness fractions, PSMC input).
        #
        # The filter drops ONLY indel records with no ALT allele. Real indel
        # variants have N_ALT>0 and are kept, as are all substitution records.
        # Not applied in block mode: `call --gvcf` already resolves this
        # collision itself (see the comment in bcftools' gvcf.c), and those
        # single-position records carry no INDEL flag to match on.
        if [ "{params.bp_mode}" = "True" ]; then
            echo "[GABi] Removing hom-ref indel duplicates (basepair mode)"
            bcftools view -e 'INFO/INDEL=1 && N_ALT=0' "concat.bcf" \
                -O b -o "called.bcf" --threads {threads}
        else
            mv "concat.bcf" "called.bcf"
        fi
        bcftools index -f --csi --threads {threads} "called.bcf"

        # FILL UNCOVERED POSITIONS
        # -------------------------------------------------------------------
        # %END is INFO/END where present (block mode) and POS otherwise
        # (basepair), so one query describes the covered span in both modes.
        echo "[GABi] Locating positions with no coverage"
        bcftools query -f '%CHROM\t%POS\t%END\n' "called.bcf" > covered.tsv
        bcftools view -h "called.bcf" > called_header.txt

        python {params.fill_script:q} \
            --fai "{input.fai}" \
            --fasta "{input.ref}" \
            --covered covered.tsv \
            --header called_header.txt \
            --output gaps.vcf

        bgzip -f gaps.vcf
        bcftools index -f -t gaps.vcf.gz
        FILL_SRC="gaps.vcf.gz"

        if [ "{params.bp_mode}" = "True" ]; then
            # Expand the no-data intervals to one record per position, so
            # basepair really does mean every position in the reference.
            echo "[GABi] Expanding no-data intervals to per-position records"
            bcftools convert --gvcf2vcf -f "{input.ref}" gaps.vcf.gz \
                -Oz -o gaps_expanded.vcf.gz
            bcftools index -f -t gaps_expanded.vcf.gz
            FILL_SRC="gaps_expanded.vcf.gz"
        fi

        # REF NORMALISATION - case, and length on hom-ref records
        # -------------------------------------------------------------------
        # (a) CASE. Three sources of lowercase REF on a soft-masked assembly:
        #       1. `convert --gvcf2vcf` copies the FASTA verbatim when expanding
        #       2. mpileup INDEL records carry the reference sequence verbatim
        #          (substitution records ARE uppercased, so one position can
        #          carry both conventions - e.g. "364 T" and "364 t INDEL")
        #       3. `call --gvcf` single-position records derived from (2)
        #
        # (b) LENGTH. A hom-ref record (ALT=".") must carry a single-base REF.
        #     mpileup emits a hom-ref INDEL CANDIDATE whose REF spans the whole
        #     indel, and `call --gvcf` opens a reference block by copying the
        #     REF of whichever record started it (gvcf.c: kputs(rec->d.allele[0],
        #     &gvcf->als)). A block opened by such a candidate therefore inherits
        #     a multi-base REF - e.g. "22925051 ATTT . END=22925099" - which is
        #     not valid gVCF. `bcftools merge --gvcf` then has to split that
        #     block at a position covered by another sample's indel, rebuilds a
        #     REF from the four-base allele, and aborts with
        #     "The REF prefixes differ".
        #     Truncating to the first base loses nothing: on a hom-ref record
        #     every base after the first is unused. Dropping these records
        #     instead would lose real depth, since some positions are covered by
        #     nothing else.
        echo "[GABi] Merging call set with no-data records, normalising REF"
        bcftools concat -a -Ov "called.bcf" "$FILL_SRC" --threads {threads} \
          | awk 'BEGIN{{OFS="\t"}}
                 /^#/ {{print; next}}
                 {{ $4 = toupper($4)
                    if ($5 == "." && length($4) > 1) $4 = substr($4, 1, 1)
                    print }}' \
          | bgzip -c -@ {threads} > "output.g.vcf.gz"
        bcftools index -f -t --threads {threads} "output.g.vcf.gz"

        # Tagging BCF. -A trims an unseen symbolic allele (<*>/<NON_REF>) if one
        # is present. `bcftools call` already resolves mpileup's <*>, so this is
        # a no-op here - verified on bcftools 1.23, 0/177 records carried <*>
        # before or after. Kept for symmetry with E11 (DeepVariant gVCFs DO
        # carry <*>, where trimming matters because F01's MULTI_ALT filter tests
        # N_ALT>1 and would otherwise count the symbolic allele as a real ALT).
        bcftools view -A -O b -o "output.raw.bcf" "output.g.vcf.gz" --threads {threads}
        bcftools index -f --csi --threads {threads} "output.raw.bcf"

        cp "output.g.vcf.gz"     "{output.gvcf}"
        cp "output.g.vcf.gz.tbi" "{output.tbi}"
        cp "output.raw.bcf"      "{output.bcf}"
        cp "output.raw.bcf.csi"  "{output.csi}"

        echo "[GABi] Done"

        # Drop the shard tmp dir once its contents are gone
        TMP_PARTS="$(dirname {input.bcfs[0]})"
        if [ -d "$TMP_PARTS" ] && [ -z "$(ls -A "$TMP_PARTS" 2>/dev/null)" ]; then
            rmdir "$TMP_PARTS" 2>/dev/null || true
        fi
        '''


# ===============================================================================
#  ROUTER RULES — caller-specific intermediates -> final output paths
# ===============================================================================
#  These rules exist so that downstream modules (F, G, H) always see the same
#  filenames ({sid}[.bp].g.vcf.gz, {sid}[.bp].raw.bcf) regardless of which caller
#  produced them.  The input functions select the GATK or DV intermediate based
#  on _caller_for_sample(), which inspects the per-sample tech count at DAG time.

rule E20_route_gvcf:
    """Route caller-specific gVCF to final output path."""
    input:
        gvcf=lambda w: _routed_caller_file(w, {
            "deep": DV_GVCF_TMPL, "bcft": BCFT_GVCF_TMPL, "gatk": GATK_GVCF_TMPL,
        }),
        tbi=lambda w: _routed_caller_file(w, {
            "deep": DV_GVCF_TBI_TMPL, "bcft": BCFT_GVCF_TBI_TMPL, "gatk": GATK_GVCF_TBI_TMPL,
        }),
    output:
        gvcf=(FINAL_GVCF_TMPL if _KEEP_BP_GVCF else temp(FINAL_GVCF_TMPL)),
        tbi=(FINAL_GVCF_TBI_TMPL if _KEEP_BP_GVCF else temp(FINAL_GVCF_TBI_TMPL))
    threads: 1
    resources:
        mem_mb=mem_func("split_intervals"),
        runtime=time_func("split_intervals")
    container: CONTAINERS["gabi_base"]
    log:
        os.path.join(
            config["OUT_FOLDER"], "GABi_results", "{species}", "{assembly}",
            "samples", "{sample_id}", "logs", "E20_route_gvcf.{sample_id}.log"
        )
    shell:
        r'''
        set -euo pipefail
        mkdir -p "$(dirname {log})"
        echo "[GABi] Routing gVCF for {wildcards.sample_id}" > {log}
        cp "{input.gvcf}" "{output.gvcf}"
        cp "{input.tbi}"  "{output.tbi}"
        echo "[GABi] Done" >> {log}
        '''


rule E21_route_bcf:
    """Route caller-specific BCF to final output path."""
    input:
        bcf=lambda w: _routed_caller_file(w, {
            "deep": DV_BCF_TMPL, "bcft": BCFT_BCF_TMPL, "gatk": GATK_BCF_TMPL,
        }),
        csi=lambda w: _routed_caller_file(w, {
            "deep": DV_BCF_CSI_TMPL, "bcft": BCFT_BCF_CSI_TMPL, "gatk": GATK_BCF_CSI_TMPL,
        }),
    output:
        bcf=FINAL_BCF_TMPL,
        csi=FINAL_BCF_CSI_TMPL
    threads: 1
    resources:
        mem_mb=mem_func("split_intervals"),
        runtime=time_func("split_intervals")
    container: CONTAINERS["gabi_base"]
    log:
        os.path.join(
            config["OUT_FOLDER"], "GABi_results", "{species}", "{assembly}",
            "samples", "{sample_id}", "logs", "E21_route_bcf.{sample_id}.log"
        )
    shell:
        r'''
        set -euo pipefail
        mkdir -p "$(dirname {log})"
        echo "[GABi] Routing BCF for {wildcards.sample_id}" > {log}
        cp "{input.bcf}" "{output.bcf}"
        cp "{input.csi}" "{output.csi}"
        echo "[GABi] Done" >> {log}
        '''


# ===============================================================================
#  JOINT GENOTYPING  (JOINT_GENO: on, any caller, >=2 samples per assembly)
# ===============================================================================
#  One rule serves all three callers because they all take the same input - the
#  per-sample gVCFs - and differ only in the tool that merges them. The tool and
#  its container are chosen from _CALLER, which is a global setting resolved at
#  DAG-build time, so no routing or runtime branching is needed.

def _joint_gvcf_indexes_for_assembly(w):
    """.tbi for each gVCF. `bcftools merge` requires indexed inputs; GLnexus
    does not, but declaring them costs nothing and keeps the two paths alike."""
    return [g + ".tbi" for g in _joint_gvcfs_for_assembly(w)]


def _cohort_sample_ids(sp, asm):
    """Sample IDs of one cohort, in data-table order, restricted to samples
    that actually have reads. Shared by every cohort-level rule so the ploidy
    file, the sample list and the BAM list can never fall out of step."""
    sids = []
    try:
        asm_data = samples_config["sp_name"][sp]["asm_id"][asm]
        for sid in asm_data.get("sample_id", {}):
            if sid in ("None", None):
                continue
            if callable_techs(sp, asm, sid, _CALLER):
                sids.append(sid)
    except (KeyError, TypeError, AttributeError):
        pass
    return sids


def _cohort_selected_bams(w):
    return [SELECTED_BAM_TMPL.format(species=w.species, assembly=w.assembly,
                                     sample_id=s)
            for s in _cohort_sample_ids(w.species, w.assembly)]


def _cohort_selected_bais(w):
    return [b + ".bai" for b in _cohort_selected_bams(w)]


def _cohort_selections(w):
    return [TECH_SELECTION_TMPL.format(species=w.species, assembly=w.assembly,
                                       sample_id=s)
            for s in _cohort_sample_ids(w.species, w.assembly)]


def _cohort_sex_inference(w):
    """sex_inference.tsv for the cohort members whose ploidy is inferred."""
    out = []
    for sid in _cohort_sample_ids(w.species, w.assembly):
        _, source = resolve_haploid_contigs(
            samples_config, w.species, w.assembly, sid
        )
        if source == "inferred":
            out.append(SEX_INFER_TSV_TMPL.format(
                species=w.species, assembly=w.assembly, sample_id=sid
            ))
    return out


def _joint_gvcfs_for_assembly(w):
    """Collect all per-sample gVCFs under a species/assembly for joint genotyping."""
    sp, asm = w.species, w.assembly
    gvcfs = []
    try:
        asm_data = samples_config["sp_name"][sp]["asm_id"][asm]
        if "sample_id" not in asm_data:
            return gvcfs
        for sid in asm_data["sample_id"]:
            if sid in ("None", None):
                continue
            sample_data = asm_data["sample_id"][sid]
            # Check sample has reads
            has_reads = False
            for rt_key in sample_data.get("read_type", {}):
                if rt_key not in ("None", None):
                    _, _, reads = _get_sample_node(sp, asm, sid, rt_key)
                    if reads:
                        has_reads = True
                        break
            if has_reads:
                gvcfs.append(
                    FINAL_GVCF_TMPL.format(
                        species=sp, assembly=asm, sample_id=sid
                    )
                )
    except (KeyError, TypeError, AttributeError):
        pass
    return gvcfs


# ------------------------------------------------------------------------------
#  Only ONE of the two paths below may define a rule producing JOINT_BCF_TMPL,
#  otherwise Snakemake reports an ambiguous rule. _JOINT_METHOD is a global
#  setting resolved at DAG-build time, so guarding the definitions is enough -
#  and it avoids a router rule that would copy the whole cohort BCF.
# ------------------------------------------------------------------------------

if _JOINT_METHOD == "bcftools-multisample":

    rule EJ20_cohort_ploidy:
        """Cohort ploidy and sample-sex files for multi-sample bcftools calling.

        Exactly the per-sample files from E29, concatenated: each sample is its
        own "sex" label, carrying its own default row plus a row per haploid
        contig. Grouping samples by shared karyotype would produce a shorter
        file, but the grouping key is the SET of haploid contigs, and since
        inference classifies each contig independently (haploid / diploid /
        ambiguous) one ambiguous call creates a new group - so the compact form
        is not reliably compact. This form is trivially correct and lets you
        diff a sample's cohort rows against its own {sid}.ploidy.txt.
        """
        input:
            fai=rules.E00_ref_faidx_for_calling.output.fai,
            sex_inference=_cohort_sex_inference,
        output:
            ploidy=COHORT_PLOIDY_TMPL,
            samples=COHORT_SAMPLES_TMPL,
        threads: 1
        resources:
            mem_mb=mem_func("split_intervals"),
            runtime=time_func("split_intervals")
        container: CONTAINERS["gabi_base"]
        log:
            os.path.join(
                config["OUT_FOLDER"], "GABi_results", "{species}", "{assembly}",
                "logs", "EJ20_cohort_ploidy.{assembly}.log"
            )
        run:
            from pathlib import Path

            sp, asm = wildcards.species, wildcards.assembly
            Path(output.ploidy).parent.mkdir(parents=True, exist_ok=True)
            Path(str(log)).parent.mkdir(parents=True, exist_ok=True)

            lengths = {}
            with open(input.fai) as fh:
                for ln in fh:
                    toks = ln.rstrip("\n").split("\t")
                    if len(toks) >= 2:
                        lengths[toks[0]] = int(toks[1])

            sids = _cohort_sample_ids(sp, asm)
            per_sample = {
                sid: haploid_contigs_for_sample(sp, asm, sid) for sid in sids
            }

            missing = sorted({c for cs in per_sample.values() for c in cs
                              if c not in lengths})
            if missing:
                raise ValueError(
                    f"[GABi] {asm}: sex contigs absent from the reference index: "
                    f"{', '.join(missing)}"
                )

            with open(output.ploidy, "w") as f:
                for sid in sids:
                    f.write(f"*\t*\t*\t{sid}\t2\n")
                for sid in sids:
                    for c in per_sample[sid]:
                        f.write(f"{c}\t1\t{lengths[c]}\t{sid}\t1\n")

            with open(output.samples, "w") as f:
                for sid in sids:
                    f.write(f"{sid}\t{sid}\n")

            with open(str(log), "w") as f:
                f.write(f"[GABi] Cohort ploidy for {asm}: {len(sids)} sample(s)\n")
                for sid in sids:
                    hap = ", ".join(per_sample[sid]) or "(none)"
                    f.write(f"[GABi]   {sid}: haploid -> {hap}\n")


    rule EJ21_multisample_call_scattered:
        """bcftools mpileup over all cohort BAMs at once, for one shard."""
        input:
            bams=_cohort_selected_bams,
            bais=_cohort_selected_bais,
            selections=_cohort_selections,
            ref=rules.D00_ref_fa_canonical.output.fa,
            fai=rules.E00_ref_faidx_for_calling.output.fai,
            bed=lambda w: os.path.join(
                vc_intervals_dir(w.species, w.assembly), f"interval_{w.shard}.bed"
            ),
            ploidy=rules.EJ20_cohort_ploidy.output.ploidy,
            samples=rules.EJ20_cohort_ploidy.output.samples,
        output:
            bcf=temp(COHORT_PART_TMPL),
            csi=temp(COHORT_PART_CSI_TMPL),
        wildcard_constraints:
            shard=r"\d{3}"
        threads: cpu_func("gatk_haplotypecaller")   # TODO: add a bcftools_call key
        resources:
            mem_mb=mem_func("gatk_haplotypecaller"),
            runtime=time_func("gatk_haplotypecaller")
        container: CONTAINERS["gabi_base"]
        log:
            os.path.join(
                config["OUT_FOLDER"], "GABi_results", "{species}", "{assembly}",
                "logs", "EJ21_multisample_call.{assembly}.{shard}.log"
            )
        shell:
            r'''
            set -euo pipefail
            mkdir -p "$(dirname {output.bcf})" "$(dirname {log})"
            exec > {log} 2>&1

            # TECHNOLOGY UNIFORMITY - the assertion behind JOINT_GENO: multi
            # ---------------------------------------------------------------
            # mpileup takes a single -X profile for every BAM. The DAG-time
            # check rejected cohorts that were provably mixed, but where a
            # sample had more than one candidate technology the outcome was
            # only decided when D07 ran. This is that check, on real results.
            TECHS=$(for j in {input.selections}; do
                        python3 -c "import json,sys; print(json.load(open(sys.argv[1]))['chosen'])" "$j"
                    done | LC_ALL=C sort -u)
            N_TECHS=$(printf '%s\n' "$TECHS" | grep -c . || true)

            if [ "$N_TECHS" -ne 1 ]; then
                echo "[GABi] ERROR: JOINT_GENO: multi needs one technology across the cohort," >&2
                echo "[GABi]        but {wildcards.assembly} selected: $(printf '%s ' $TECHS)" >&2
                echo "[GABi]        bcftools mpileup applies one -X profile to every BAM." >&2
                echo "[GABi]        Use JOINT_GENO: merge instead." >&2
                exit 1
            fi

            case "$TECHS" in
                hifi)         PROFILE=pacbio-ccs-1.20 ;;
                ont)          PROFILE=ont-sup-1.20 ;;
                illumina|10x) PROFILE=illumina-1.20 ;;
                *)            echo "[GABi] ERROR: unknown technology '$TECHS'" >&2; exit 1 ;;
            esac

            # mktemp, not a bare filename: the rule does not cd anywhere, so a
            # relative path lands in Snakemake's working directory - and all
            # shards of a cohort would race on the same file. $TMPDIR is set to
            # GABi's temp dir by the shell prefix.
            BAM_LIST="$(mktemp "${{TMPDIR:-/tmp}}/GABi_bamlist_{wildcards.assembly}_{wildcards.shard}_XXXXXX")"
            trap 'rm -f "$BAM_LIST"' EXIT
            printf '%s\n' {input.bams} > "$BAM_LIST"
            N_BAMS=$(wc -l < "$BAM_LIST")

            echo "[GABi] Multi-sample calling: {wildcards.assembly} shard {wildcards.shard}"
            echo "[GABi] Samples: $N_BAMS  |  technology: $TECHS  |  profile: $PROFILE"

            # -mv: variant sites only. A cohort call set is a variant matrix -
            # this is what GLnexus emits too - and an all-sites cohort file
            # would be n_samples x genome in size for no downstream gain.
            bcftools mpileup \
                -X "$PROFILE" \
                -f "{input.ref}" \
                -a FORMAT/AD,FORMAT/DP \
                -R "{input.bed}" \
                -b "$BAM_LIST" \
                -Ou \
            | bcftools call \
                -m -v \
                -a GQ \
                --ploidy-file "{input.ploidy}" \
                -S "{input.samples}" \
                -Ob -o "{output.bcf}" \
                --threads {threads}

            bcftools index -f --csi --threads {threads} "{output.bcf}"
            echo "[GABi] Done"
            '''


    rule EJ22_multisample_gather:
        """Concatenate the per-shard cohort BCFs into the joint call set."""
        input:
            bcfs=lambda w: [
                COHORT_PART_TMPL.format(species=w.species, assembly=w.assembly, shard=s)
                for s in vc_discover_shards(w.species, w.assembly)
            ],
            csis=lambda w: [
                COHORT_PART_CSI_TMPL.format(species=w.species, assembly=w.assembly, shard=s)
                for s in vc_discover_shards(w.species, w.assembly)
            ],
        output:
            bcf=JOINT_BCF_TMPL,
            csi=JOINT_BCF_CSI_TMPL,
        threads: cpu_func("bcftools_concat")
        resources:
            mem_mb=mem_func("bcftools_concat"),
            runtime=time_func("bcftools_concat")
        container: CONTAINERS["gabi_base"]
        log:
            os.path.join(
                config["OUT_FOLDER"], "GABi_results", "{species}", "{assembly}",
                "logs", "EJ22_multisample_gather.{assembly}.log"
            )
        shell:
            r'''
            set -euo pipefail
            mkdir -p "$(dirname {output.bcf})" "$(dirname {log})"
            exec > {log} 2>&1

            set -- {input.bcfs}
            echo "[GABi] Concatenating $# cohort shard(s) for {wildcards.assembly}"

            export HTS_COMPRESSION_LEVEL=1
            bcftools concat -a -D -O b -o "{output.bcf}" {input.bcfs} --threads {threads}
            bcftools index -f --csi --threads {threads} "{output.bcf}"

            echo "[GABi] Done"

            TMP_PARTS="$(dirname {input.bcfs[0]})"
            if [ -d "$TMP_PARTS" ] && [ -z "$(ls -A "$TMP_PARTS" 2>/dev/null)" ]; then
                rmdir "$TMP_PARTS" 2>/dev/null || true
            fi
            '''

# Second `if` rather than `else`: a top-level `else` immediately after a rule
# definition asks Snakemake's parser to dedent out of rule context onto a
# keyword, which it does not always handle. Two independent conditions are
# unambiguous and mutually exclusive by construction.
if _JOINT_METHOD != "bcftools-multisample":

    rule EJ10_joint_genotype:
        """Combine the per-sample gVCFs of one assembly into a cohort call set.

        Runs when JOINT_GENO: on and the assembly has >=2 samples, for any caller.
        GLnexus for DeepVariant and GATK (it ships a preset for each); bcftools
        merge --gvcf for bcftools, which GLnexus has no preset for. The container
        follows the tool: GLnexus is not present in gabi_base.
        """
        input:
            gvcfs=_joint_gvcfs_for_assembly,
            tbis=_joint_gvcf_indexes_for_assembly,
            ref=rules.D00_ref_fa_canonical.output.fa,
            fai=rules.E00_ref_faidx_for_calling.output.fai,
        output:
            bcf=JOINT_BCF_TMPL,
            csi=JOINT_BCF_CSI_TMPL,
        params:
            method=_JOINT_METHOD,
            glnexus_config=_GLNEXUS_CONFIG,
        threads: cpu_func("glnexus")
        resources:
            mem_mb=mem_func("glnexus"),
            runtime=time_func("glnexus")
        container: _JOINT_CONTAINER
        log:
            os.path.join(
                config["OUT_FOLDER"], "GABi_results", "{species}", "{assembly}",
                "logs", "EJ10_joint_genotype.{assembly}.log"
            )
        shell:
            r'''
            set -euo pipefail
            mkdir -p "$(dirname {output.bcf})" "$(dirname {log})"

            exec > {log} 2>&1

            N_SAMPLES=$(echo {input.gvcfs} | wc -w)
            echo "[GABi] Joint genotyping for {wildcards.assembly}"
            echo "[GABi] Method: {params.method}"
            echo "[GABi] Samples: $N_SAMPLES"
            echo "[GABi] Threads: {threads}"


            # TEMP DIRECTORY SETUP
            # -------------------------------------------------------------------
            WORK_DIR="$(gabi_get_workdir 100)"
            TEMP_DIR="$(mktemp -d "$WORK_DIR/GABi_joint_{wildcards.species}_{wildcards.assembly}_XXXXXX")"
            trap 'rm -rf "$TEMP_DIR"' EXIT
            echo "[GABi] Working directory: $TEMP_DIR"

            if [ "{params.method}" = "glnexus" ]; then
                # Convert mem_mb to GB for GLnexus (floor, minimum 4)
                MEM_GB=$(( {resources.mem_mb} / 1024 ))
                [ "$MEM_GB" -lt 4 ] && MEM_GB=4
                echo "[GABi] GLnexus config: {params.glnexus_config}"

                # A value containing "/" is a custom preset FILE (gatk +
                # basepair); anything else is a built-in preset name. Check the
                # file up front - glnexus_cli's own error for a missing config
                # is not obvious, and this runs after a long bulk load.
                case "{params.glnexus_config}" in
                    */*)
                        if [ ! -s "{params.glnexus_config}" ]; then
                            echo "[GABi] ERROR: GLnexus config file not found:" >&2
                            echo "[GABi]        {params.glnexus_config}" >&2
                            echo "[GABi]        It should sit in workflow/scripts/" >&2
                            exit 1
                        fi
                        echo "[GABi] Custom preset (BP_RESOLUTION gVCFs carry FORMAT/DP, not MIN_DP)"
                        ;;
                esac

                glnexus_cli \
                    --config {params.glnexus_config} \
                    --dir "$TEMP_DIR/glnexus_db" \
                    --threads {threads} \
                    --mem-gbytes "$MEM_GB" \
                    {input.gvcfs} \
                > "$TEMP_DIR/joint.bcf"
            else
                # bcftools merge --gvcf needs the reference to fill in bases at the
                # points where reference blocks from different samples are split.
                echo "[GABi] Merging gVCF blocks with bcftools merge --gvcf"
                bcftools merge --gvcf "{input.ref}" \
                    -Ob -o "$TEMP_DIR/joint.bcf" \
                    --threads {threads} \
                    {input.gvcfs}
            fi

            echo "[GABi] Creating CSI index"
            bcftools index -f --csi --threads {threads} "$TEMP_DIR/joint.bcf"

            # Copy results to final location
            echo "[GABi] Copying output to final location..."
            cp "$TEMP_DIR/joint.bcf"     "{output.bcf}"
            cp "$TEMP_DIR/joint.bcf.csi" "{output.csi}"

            echo "[GABi] Joint genotyping complete ($N_SAMPLES samples)"

            '''
