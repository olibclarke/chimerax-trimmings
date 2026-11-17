"""Dependency-free sequence features for Trimmings; no ChimeraX/Qt imports.

SEG is a Python port of NCBI's public-domain blast_seg.c (Ilya Dondoshansky),
using Wootton and Federhen's entropy triggering and probability-based trimming.
Reference: https://github.com/ncbi/ncbi-cxx-toolkit-public/blob/master/src/algo/blast/core/blast_seg.c
PROSITE definitions are downloaded unchanged from SIB and retained with their
copyright/license notice. They are data, not bundled or modified definitions.
"""

from collections import Counter
from functools import lru_cache
from hashlib import sha256
import json
import math
from pathlib import Path
import re
from time import monotonic, time
from urllib.parse import urlencode
from urllib.request import Request, urlopen
from uuid import uuid4


VERSION = 1
PROSITE_URL = "https://ftp.expasy.org/databases/prosite/prosite.dat"
SCAN_URL = "https://prosite.expasy.org/cgi-bin/prosite/scanprosite/PSScan.cgi"
AMINO_ACIDS = frozenset("ACDEFGHIKLMNPQRSTVWY")
CACHE_SECONDS = 14 * 86400


def _check(cancel, deadline=None):
    if cancel.is_set():
        raise InterruptedError("Sequence-feature lookup cancelled.")
    if deadline is not None and monotonic() >= deadline:
        raise TimeoutError("Sequence-feature lookup timed out.")


def _atomic_write(path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_name(path.name + "." + uuid4().hex + ".tmp")
    try:
        temp.write_bytes(data)
        temp.replace(path)
    finally:
        temp.unlink(missing_ok=True)


def _request(url, timeout, data=None):
    request = Request(url, data=data, headers={"User-Agent": "ChimeraX-Trimmings/1.0",
                      "Content-Type": "application/x-www-form-urlencoded"})
    with urlopen(request, timeout=timeout) as response:
        content = response.read(40 * 1024 * 1024 + 1)
    if len(content) > 40 * 1024 * 1024:
        raise ValueError("Sequence-feature response exceeded the download limit.")
    return content


@lru_cache(maxsize=16384)
def _lnfact(n):
    # Match the six-decimal NCBI lookup table, including tie-breaking in trimming.
    return round(math.lgamma(n + 1), 6) if n <= 10000 else ((n + .5) * math.log(n) - n + .9189385332)


def _seg_probability(counts, length):
    state = sorted((n for n in counts.values() if n), reverse=True)
    classes = Counter(state + [0] * (20 - len(state)))
    assignments = _lnfact(20) - sum(_lnfact(n) for n in classes.values())
    permutations = _lnfact(length) - sum(_lnfact(n) for n in state)
    return assignments + permutations - length * 2.9957322735539909


def _seg_trim(sequence, cancel):
    size, best, left, right = len(sequence), 1.0, 0, len(sequence) - 1
    for length in range(size, max(1, size - 50), -1):
        _check(cancel)
        counts = Counter(aa for aa in sequence[:length] if aa in AMINO_ACIDS)
        for start in range(size - length + 1):
            probability = _seg_probability(counts, length)
            if probability < best:
                best, left, right = probability, start, start + length - 1
            if start + length < size:
                old, new = sequence[start], sequence[start + length]
                if old in AMINO_ACIDS:
                    counts[old] -= 1
                if new in AMINO_ACIDS:
                    counts[new] += 1
    return left, right


def seg(sequence, cancel):
    """Return inclusive one-based SEG intervals (12, 2.2, 2.5; maxtrim=50).

    Ambiguous residues are excluded from the 20-letter composition; windows
    containing >2 ambiguous residues cannot trigger or extend a segment.
    An explicit stack replaces the reference implementation's left-trim recursion.
    """
    pending, regions = [(sequence, 0)], []
    while pending:
        seq, offset = pending.pop()
        if len(seq) < 12:
            continue
        _check(cancel)
        entropy = [-1.0] * len(seq)
        for start in range(len(seq) - 11):
            counts = Counter(aa for aa in seq[start:start + 12] if aa in AMINO_ACIDS)
            total = sum(counts.values())
            if total >= 10:
                entropy[start + 5] = abs(sum(n * math.log(n / total) / math.log(2)
                                                           for n in sorted(counts.values(), reverse=True)) / total)
        i, last, lowlim = 5, len(seq) - 7, 5
        while i <= last:
            if 0 <= entropy[i] <= 2.2:
                lo, hi = i, i
                while lo >= lowlim and 0 <= entropy[lo] <= 2.5:
                    lo -= 1
                lo += 1
                while hi <= last and 0 <= entropy[hi] <= 2.5:
                    hi += 1
                hi -= 1
                left, right = lo - 5, hi + 6
                trim_left, trim_right = _seg_trim(seq[left:right + 1], cancel)
                left, right = left + trim_left, left + trim_right
                if i + 6 < left:
                    pending.append((seq[lo - 5:left], offset + lo - 5))
                regions.append((left + offset + 1, right + offset + 1))
                i = min(hi, right + 5)
                lowlim = i + 1
            i += 1
    return sorted(set(regions))


def low_complexity(sequence, cancel):
    hits = []
    for start, end in seg(sequence, cancel):
        composition = Counter(sequence[start - 1:end])
        summary = ", ".join(f"{aa} {n / (end - start + 1):.0%}" for aa, n in composition.most_common(3))
        hits.append({"origin": "sequence", "feature_kind": "low_complexity", "type": "Low complexity",
                     "name": "Low-complexity region", "accession": "SEG", "source": "SEG",
                     "start": start, "end": end, "evalue": None, "bitscore": None, "outcompeted": False,
                     "entry_url": "https://www.ncbi.nlm.nih.gov/books/NBK279684/",
                     "structures_url": "", "description": "Composition: " + summary + ".",
                     "evidence": "SEG prediction; window 12, trigger 2.2, extension 2.5, maxtrim 50. "
                                 "Low complexity does not establish intrinsic disorder."})
    return hits


def compile_pattern(pattern):
    """Translate PROSITE's constrained pattern language; reject arbitrary regex.

    Character sets containing < or > allow that terminus instead of a residue.
    Positive lookahead preserves overlapping occurrences. Conserved/excluded
    positions never accept an ambiguous residue unless explicitly listed.
    """
    pattern = pattern.strip().removesuffix(".")
    tokens = re.findall(r"\[[A-Z<>]+\]|\{[A-Z]+\}|[A-Zx<>]|\(\d+(?:,\d+)?\)|-", pattern)
    if "".join(tokens) != pattern:
        raise ValueError("Unsupported PROSITE pattern: " + pattern)
    pieces = []
    repeatable = False
    for token in tokens:
        if token == "-":
            repeatable = False
        elif token.startswith("("):
            if not repeatable:
                raise ValueError("Invalid PROSITE repetition.")
            bounds = [int(n) for n in token[1:-1].split(",")]
            if not 0 <= bounds[0] <= bounds[-1] <= 1000:
                raise ValueError("Invalid PROSITE repetition bounds.")
            pieces.append("{" + token[1:-1] + "}")
            repeatable = False
        else:
            repeatable = token not in ("<", ">")
            if token.startswith("["):
                letters = token[1:-1].replace("<", "").replace(">", "")
                alternatives = (["[" + letters + "]"] if letters else [])
                alternatives += (["\\A"] if "<" in token else []) + (["\\Z"] if ">" in token else [])
                pieces.append("(?:" + "|".join(alternatives) + ")")
            elif token.startswith("{"):
                allowed = "".join(sorted(AMINO_ACIDS - set(token[1:-1])))
                pieces.append("[" + allowed + "]")
            else:
                pieces.append({"x": "[A-Z]", "<": "\\A", ">": "\\Z"}.get(token, token))
    regex = re.compile("(?=(" + "".join(pieces) + "))")
    if regex.match(""):
        raise ValueError("A PROSITE pattern must consume at least one residue.")
    return regex


def parse_prosite(text):
    release = re.search(r"^CC\s+Release (.+?)\.$", text, re.M)
    entries = {}
    for block in text.split("//"):
        identifier = re.search(r"^ID\s+(\w+); (PATTERN|MATRIX)\.", block, re.M)
        if not identifier:
            continue
        accession = re.search(r"^AC\s+(PS\d{5});", block, re.M)
        if not accession:
            raise ValueError("Invalid PROSITE accession.")
        acc = accession[1]
        if acc in entries:
            raise ValueError("Duplicate PROSITE accession.")
        pattern = "".join(re.findall(r"^PA\s+(.+)$", block, re.M))
        if identifier[2] == "PATTERN" and not pattern:
            raise ValueError("Missing PROSITE pattern.")
        entries[acc] = {"accession": acc, "name": identifier[1], "kind": identifier[2],
                        "description": " ".join(re.findall(r"^DE\s+(.+)$", block, re.M)),
                        "common": bool(re.search(r"/SKIP-FLAG=TRUE;", block)), "pattern": pattern,
                        "prosite_release": release[1] if release else "local file",
                        "regex": compile_pattern(pattern) if pattern else None}
    if not entries:
        raise ValueError("No PROSITE entries found; the response may be an error page.")
    return entries


@lru_cache(maxsize=2)
def _local_prosite(path, modified, size):
    """Keep compiled definitions between commands; file edits invalidate them."""
    data = Path(path).read_bytes()
    return parse_prosite(data.decode("utf-8")), sha256(data).hexdigest()


def load_prosite(cache_dir, cancel, deadline, database=None, refresh=False, report=lambda text: None):
    path = Path(database).expanduser() if database else Path(cache_dir) / "prosite.dat"
    _check(cancel, deadline)
    if database or path.is_file() and not refresh:
        stamp = path.stat()
        entries, digest = _local_prosite(str(path.resolve()), stamp.st_mtime_ns, stamp.st_size)
    else:
        report("downloading PROSITE definitions (once; sequences stay local)")
        _check(cancel, deadline)
        data = _request(PROSITE_URL, min(60, max(.1, deadline - monotonic())))
        entries, digest = parse_prosite(data.decode("utf-8")), sha256(data).hexdigest()
        _check(cancel, deadline)
        _atomic_write(path, data)
    _check(cancel, deadline)
    return entries, digest


def pattern_hits(sequence, entries, common, cancel):
    results = []
    for entry in entries.values():
        _check(cancel)
        if entry["kind"] != "PATTERN" or entry["common"] and not common:
            continue
        for match in entry["regex"].finditer(sequence):
            start, end = match.start(1) + 1, match.end(1)
            results.append(_prosite_hit(entry, sequence, start, end, "pattern"))
    return results


def _prosite_hit(entry, sequence, start, end, kind):
    return {
        **{key: value for key, value in entry.items() if key not in ("regex", "kind")},
        "origin": "prosite", "feature_kind": "motif" if kind == "pattern" else "profile",
        "database": "prosite" if kind == "pattern" else "profile", "source": "PROSITE",
        "type": "Pattern match" if kind == "pattern" else "Profile-supported region",
        "start": start, "end": end, "matched_sequence": sequence[start - 1:end],
        "evalue": None, "bitscore": None, "outcompeted": False,
        "entry_url": "https://prosite.expasy.org/" + entry["accession"],
        "structures_url": "", "evidence": "Sequence pattern match; function is not established by this match alone."
        if kind == "pattern" else "ScanProsite profile match at the standard level-0 cutoff. "
                                  "A profile region alone does not establish a particular residue's function.",
    }


def parse_scanprosite(payload, sequence, identifier, entries, common=False):
    try:
        if (not isinstance(payload, dict) or type(payload["n_match"]) is not int
                or type(payload["n_seq"]) is not int or payload["n_seq"] not in (0, 1)
                or not isinstance(payload["matchset"], list)
                or payload["n_match"] != len(payload["matchset"])):
            raise ValueError()
        results = []
        for match in payload["matchset"]:
            if match.get("sequence_id") != identifier and match.get("sequence_ac") != identifier:
                raise ValueError()
            acc, start, end = match["signature_ac"], match["start"], match["stop"]
            if (acc not in entries or type(start) is not int or type(end) is not int
                    or not 1 <= start <= end <= len(sequence)):
                raise ValueError()
            entry = entries[acc]
            if entry["common"] and not common:
                continue
            if entry["kind"] == "MATRIX":
                score, level = float(match["score"]), int(match["level"])
                if not math.isfinite(score) or level < 0:
                    raise ValueError()
                hit = _prosite_hit(entry, sequence, start, end, "profile")
                hit.update(score=score, profile_level=level, source="ScanProsite")
                results.append(hit)
        return results
    except (KeyError, TypeError, ValueError, OverflowError):
        raise ValueError("Invalid ScanProsite results or sequence coordinates; local features were retained.") from None


def search(sequences, cache_dir, timeout, cancel, report, *, database=None,
           common=False, profiles=False, refresh=False):
    """Local SEG + PROSITE, optionally enriched online; successful zeroes cache too.

    Local features survive a failed definition download or optional online search.
    Failures are warnings, never cached as valid empty remote responses.
    """
    deadline = monotonic() + timeout
    cache_dir = Path(cache_dir)
    results, warnings = {}, []
    try:
        entries, digest = load_prosite(cache_dir, cancel, deadline, database, refresh, report)
    except Exception as error:
        _check(cancel, deadline)
        if database:
            raise
        entries, digest = {}, "unavailable"
        warnings.append("PROSITE patterns unavailable; low-complexity features retained: " + str(error))
    for index, sequence in enumerate(dict.fromkeys(sequences), 1):
        _check(cancel, deadline)
        report(f"local sequence {index}/{len(sequences)}")
        fingerprint = json.dumps([VERSION, sequence, digest, common])
        path = cache_dir / "results" / (sha256(fingerprint.encode()).hexdigest() + ".json")
        local = _read_cache(path, len(sequence)) if entries and not refresh else None
        if local is None:
            local = low_complexity(sequence, cancel) + pattern_hits(sequence, entries, common, cancel)
            if entries:
                _atomic_write(path, json.dumps({"fetched": time(), "hits": local}).encode())
        results[sequence] = list(local)
    # Complete every local scan before waiting on the optional network service.
    # One slow/failed remote request must not prevent other chains' local results.
    if profiles and entries:
        for index, sequence in enumerate(results, 1):
            _check(cancel)
            if monotonic() >= deadline:
                warnings.append("ScanProsite profile time budget exhausted; local features retained.")
                break
            fingerprint = json.dumps([VERSION, sequence, digest, common])
            path = cache_dir / "results" / (sha256(fingerprint.encode()).hexdigest() + ".json")
            online_path = path.with_name(path.stem + "-profiles.json")
            online = _read_cache(online_path, len(sequence)) if not refresh else None
            if online is None:
                report(f"ScanProsite profiles {index}/{len(sequences)} (sending sequence)")
                identifier = "trimmings_" + sha256(sequence.encode()).hexdigest()
                data = urlencode({"seq": ">" + identifier + "\n" + sequence, "output": "json",
                                  "skip": "0" if common else "1", "lowscore": "0"}).encode()
                try:
                    payload = json.loads(_request(SCAN_URL, min(60, max(.1, deadline - monotonic())), data))
                    online = parse_scanprosite(payload, sequence, identifier, entries, common)
                    _check(cancel, deadline)
                    _atomic_write(online_path, json.dumps({"fetched": time(), "hits": online}).encode())
                except Exception as error:
                    _check(cancel)
                    warnings.append("ScanProsite profiles unavailable; local features retained: " + str(error))
                    online = []
            results[sequence].extend(online)
    return results, warnings


def _read_cache(path, sequence_length):
    try:
        payload = json.loads(path.read_text())
        if (0 <= time() - payload["fetched"] < CACHE_SECONDS and isinstance(payload["hits"], list)
                and all(isinstance(hit, dict) and {"start", "end", "accession", "origin", "feature_kind",
                                                  "name", "type", "source", "outcompeted", "evalue", "bitscore"} <= hit.keys()
                        and type(hit["start"]) is int and type(hit["end"]) is int
                        and 1 <= hit["start"] <= hit["end"] <= sequence_length for hit in payload["hits"])):
            return payload["hits"]
    except (OSError, ValueError, TypeError, KeyError):
        pass
    return None
