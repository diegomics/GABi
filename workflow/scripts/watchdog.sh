#!/bin/bash

# Download helpers
# by Diego De Panis, 2026
# This script is part of the GABi pipeline
# note: AI tools may have been used to improve, clean and/or comment this version of the code

gabi_download_with_timeout() {
    local STALL_TIMEOUT=${1:-3600}
    shift

    echo "[GABi] Running with ${STALL_TIMEOUT}s timeout: $@"
    timeout --signal=TERM --kill-after=30 "$STALL_TIMEOUT" "$@"
    return $?
}


# ---------------------------------------------------------------------
# HTTP fetch (to stdout): curl preferred, wget fallback, python urllib
# last resort. Used for short API responses where we need the body
# inline. Retries built in. Diagnostics to stderr.
# ---------------------------------------------------------------------
gabi_http_fetch_stdout() {
    # gabi_http_fetch_stdout <url>
    local URL=$1

    if command -v curl >/dev/null 2>&1; then
        # --fail: without it an HTTP error page arrives with exit status 0,
        # and the caller parses the error page as data. wget and urllib
        # below already fail on HTTP errors.
        curl -sS --fail --max-time 60 --retry 3 --retry-delay 10 "$URL"
    elif command -v wget >/dev/null 2>&1; then
        wget -qO- --tries=3 --waitretry=10 --timeout=60 "$URL"
    elif command -v python3 >/dev/null 2>&1 || command -v python >/dev/null 2>&1; then
        local PY
        PY=$(command -v python3 || command -v python)
        "$PY" -c '
import sys, urllib.request
req = urllib.request.Request(sys.argv[1], headers={"User-Agent": "GABi/urllib"})
sys.stdout.write(urllib.request.urlopen(req, timeout=60).read().decode())
' "$URL"
    else
        echo "[GABi] No HTTP client available in container" >&2
        return 127
    fi
}

# ---------------------------------------------------------------------
# HTTP fetch (to file): same backend preference as above, designed for
# large downloads. Retries and generous timeouts built in.
# ---------------------------------------------------------------------
gabi_http_fetch() {
    # gabi_http_fetch <url> <output_file>
    # Downloads with retries, returns tool exit status.
    local URL=$1
    local OUT=$2

    if command -v curl >/dev/null 2>&1; then
        curl -sS -L --fail \
            --connect-timeout 30 \
            --max-time 14400 \
            --retry 5 \
            --retry-delay 60 \
            --retry-max-time 10800 \
            -o "$OUT" "$URL"
    elif command -v wget >/dev/null 2>&1; then
        wget -q --tries=5 --waitretry=60 --timeout=14400 -O "$OUT" "$URL"
    elif command -v python3 >/dev/null 2>&1 || command -v python >/dev/null 2>&1; then
        local PY
        PY=$(command -v python3 || command -v python)
        # Stdlib-only: urllib.request with manual retries. No deps.
        "$PY" - "$URL" "$OUT" <<'PYEOF'
import sys, time, urllib.request, urllib.error, shutil
url, out = sys.argv[1], sys.argv[2]
for attempt in range(1, 6):
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "GABi/urllib"})
        with urllib.request.urlopen(req, timeout=14400) as r, open(out, "wb") as f:
            shutil.copyfileobj(r, f, length=1024 * 1024)
        sys.exit(0)
    except (urllib.error.URLError, OSError) as e:
        sys.stderr.write(f"[GABi]   python urllib attempt {attempt} failed: {e}\n")
        if attempt < 5:
            time.sleep(60)
sys.exit(1)
PYEOF
    else
        echo "[GABi] No HTTP client (curl/wget/python) available in container" >&2
        return 127
    fi
}


# ---------------------------------------------------------------------
# True if an ENA submitted_format value names FASTQ ("FASTQ", "fastq",
# "FASTQ;FASTQ"). A case, not `echo | grep -q`: under pipefail, grep -q
# exiting at the first match can SIGPIPE the producer and flip the answer.
# ---------------------------------------------------------------------
_gabi_is_fastq_format() {
    case "$1" in
        *[Ff][Aa][Ss][Tt][Qq]*) return 0 ;;
    esac
    return 1
}


# ---------------------------------------------------------------------
# Could this run's SUBMITTED files be FASTQ? Asked before
# `enaDataGet.py -f submitted`, which otherwise fetches whatever was
# submitted - for PacBio HiFi usually the unaligned BAM (BAM;BAI), often
# larger than the FASTQ it stands in for - only for the rule to discard it.
#
# Fails OPEN: returns 1 ("skip") only when ENA gives a well-formed answer
# naming a non-FASTQ format. A failed query, an unexpected answer or an
# empty field all return 0, keeping the old behaviour, so an ENA hiccup can
# never remove a fallback that might have worked.
# ---------------------------------------------------------------------
gabi_ena_submitted_may_be_fastq() {
    local ACC=$1 META="" HEADER="" ROW="" RUN="" FMT=""
    META=$(gabi_http_fetch_stdout \
        "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${ACC}&result=read_run&fields=submitted_format&format=tsv") || return 0
    { IFS= read -r HEADER; IFS= read -r ROW; } <<< "$META" || true
    HEADER=${HEADER%$'\r'}
    ROW=${ROW%$'\r'}
    [ "$HEADER" = "$(printf 'run_accession\tsubmitted_format')" ] || return 0
    IFS=$'\t' read -r RUN FMT <<< "$ROW" || true
    [ "$RUN" = "$ACC" ] && [ -n "$FMT" ] || return 0
    _gabi_is_fastq_format "$FMT" && return 0
    echo "[GABi] ENA lists the submitted files of ${ACC} as '${FMT}', not FASTQ:" >&2
    echo "[GABi]   skipping their download - the rule can only use FASTQ" >&2
    return 1
}


# ---------------------------------------------------------------------
# Query ENA portal API filereport for a run accession and print the
# HTTPS URLs (one per line) to stdout. Diagnostics go to stderr.
#
# Prefers fastq_ftp; falls back to submitted_ftp when submitted_format
# contains FASTQ.
#
# Exit codes:
#   0 success
#   2 network error reaching ENA portal
#   3 empty response / no data row (accession not found?)
#   4 no FASTQ available (submission may be CRAM/BAM-only)
#   5 URL count doesn't match requested mode
# ---------------------------------------------------------------------
gabi_ena_get_urls() {
    local ACC=$1
    local MODE=${2:-paired}   # paired | single
    local API_URL="https://www.ebi.ac.uk/ena/portal/api/filereport"
    local FIELDS="fastq_ftp,submitted_ftp,submitted_format"
    local META ROW URLS SOURCE COUNT FETCH_EXIT

    echo "[GABi] Querying ENA portal API for ${ACC} (mode=${MODE})..." >&2

    # Not `META=$(...)` + `FETCH_EXIT=$?`: under set -e a failing
    # substitution exits on the assignment line, before $? is read.
    FETCH_EXIT=0
    META=$(gabi_http_fetch_stdout \
        "${API_URL}?accession=${ACC}&result=read_run&fields=${FIELDS}&format=tsv") || FETCH_EXIT=$?
    if [ "$FETCH_EXIT" -ne 0 ]; then
        echo "[GABi] ENA portal API query failed (exit ${FETCH_EXIT})" >&2
        return 2
    fi

    # Judge the answer by its content, not the exit status: a 200 can still
    # carry a maintenance page. ENA always puts run_accession first.
    local HEADER=""
    { IFS= read -r HEADER; } <<< "$META" || true
    HEADER=${HEADER%$'\r'}
    case "$HEADER" in
        run_accession$'\t'*) ;;
        *)
            echo "[GABi] ENA portal API returned an unexpected response for ${ACC}" >&2
            echo "[GABi]   first line: ${HEADER:0:200}" >&2
            return 2
            ;;
    esac

    ROW=$(echo "$META" | awk 'NR==2')
    ROW=${ROW%$'\r'}
    if [ -z "$ROW" ]; then
        echo "[GABi] ENA portal returned no data row for ${ACC}" >&2
        echo "[GABi] Response was:" >&2
        echo "$META" >&2
        return 3
    fi

    local FASTQ_FTP SUBMITTED_FTP SUBMITTED_FORMAT
    # ENA portal always prepends run_accession as column 1, regardless
    # of what fields= requests. Our requested fields start at column 2.
    FASTQ_FTP=$(echo "$ROW"        | awk -F'\t' '{print $2}')
    SUBMITTED_FTP=$(echo "$ROW"    | awk -F'\t' '{print $3}')
    SUBMITTED_FORMAT=$(echo "$ROW" | awk -F'\t' '{print $4}')

    if [ -n "$FASTQ_FTP" ]; then
        URLS="$FASTQ_FTP"
        SOURCE="fastq_ftp"
    elif [ -n "$SUBMITTED_FTP" ] && _gabi_is_fastq_format "$SUBMITTED_FORMAT"; then
        URLS="$SUBMITTED_FTP"
        SOURCE="submitted_ftp (FASTQ)"
    else
        echo "[GABi] No FASTQ available via ENA portal for ${ACC}" >&2
        echo "[GABi]   fastq_ftp='${FASTQ_FTP}'" >&2
        echo "[GABi]   submitted_ftp='${SUBMITTED_FTP}'" >&2
        echo "[GABi]   submitted_format='${SUBMITTED_FORMAT}'" >&2
        return 4
    fi

    COUNT=$(echo "$URLS" | tr ';' '\n' | grep -c .)
    if [ "$MODE" = "paired" ] && [ "$COUNT" -ne 2 ]; then
        echo "[GABi] Expected 2 URLs for paired-end ${ACC}, got ${COUNT}" >&2
        return 5
    fi
    if [ "$MODE" = "single" ] && [ "$COUNT" -lt 1 ]; then
        echo "[GABi] No URLs for ${ACC}" >&2
        return 5
    fi

    echo "[GABi] Found ${COUNT} URL(s) for ${ACC} via ${SOURCE}" >&2
    # Emit URLs one per line, prepending https:// (API returns them
    # as bare hostnames, e.g. ftp.sra.ebi.ac.uk/vol1/...)
    echo "$URLS" | tr ';' '\n' | awk 'NF > 0 {print "https://" $0}'
    return 0
}


# ---------------------------------------------------------------------
# Download paired-end reads via portal API + HTTPS.
# CWD must be the destination directory.
# Produces <acc>_1.fastq.gz and <acc>_2.fastq.gz on success.
# ---------------------------------------------------------------------
gabi_ena_download_paired() {
    local ACC=$1
    local URLS URL1 URL2 API_EXIT

    API_EXIT=0
    URLS=$(gabi_ena_get_urls "$ACC" paired) || API_EXIT=$?
    if [ "$API_EXIT" -ne 0 ]; then
        return "$API_EXIT"
    fi

    URL1=$(echo "$URLS" | sed -n '1p')
    URL2=$(echo "$URLS" | sed -n '2p')

    gabi_http_fetch "$URL1" "${ACC}_1.fastq.gz" || return 10
    gabi_http_fetch "$URL2" "${ACC}_2.fastq.gz" || return 11
    return 0
}


# ---------------------------------------------------------------------
# Download single-end/long reads via portal API + HTTPS.
# CWD must be the destination directory.
# Produces <acc>.fastq.gz on success.
# ---------------------------------------------------------------------
gabi_ena_download_single() {
    local ACC=$1
    local URLS URL1 API_EXIT

    API_EXIT=0
    URLS=$(gabi_ena_get_urls "$ACC" single) || API_EXIT=$?
    if [ "$API_EXIT" -ne 0 ]; then
        return "$API_EXIT"
    fi

    URL1=$(echo "$URLS" | sed -n '1p')
    gabi_http_fetch "$URL1" "${ACC}.fastq.gz" || return 10
    return 0
}
