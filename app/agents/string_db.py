"""STRING-db association scores between candidate genes and the seed proteins.

Keyless REST API (caller_identity only). An optional STRING_API_KEY env var is
sent only as a courtesy/fallback. For each candidate we fetch BOTH the functional
and physical networks restricted to [candidate + seeds], keep the edges touching
the candidate, and aggregate across seeds (max / mean / best-seed / count >= cutoff).
Pairs STRING has no evidence for are treated as score 0.

All read-only, cached, graceful-empty on failure. Never writes user_data/.
"""
from __future__ import annotations

import os
import urllib.parse
import urllib.request
from functools import lru_cache

STRING_BASE = "https://string-db.org/api"
_CALLER = "scrap-ai"


def _post(path: str, data: dict, timeout: float = 40.0) -> bytes:
    key = os.environ.get("STRING_API_KEY")
    if key:
        data = {**data, "api_key": key}
    body = urllib.parse.urlencode(data).encode()
    req = urllib.request.Request(f"{STRING_BASE}/{path}", data=body,
                                 headers={"User-Agent": "SCRAP-AI/1.0"})
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        return resp.read()


def _network_edges(candidate: str, seeds: tuple[str, ...], physical: bool) -> dict[str, float]:
    """Return {seed: combined_score} for edges STRING reports touching candidate."""
    ids = "\r".join([candidate, *seeds])
    data = {"identifiers": ids, "species": 9606, "required_score": 0,
            "caller_identity": _CALLER}
    if physical:
        data["network_type"] = "physical"
    try:
        txt = _post("tsv/network", data).decode()
    except Exception:
        return {}
    out: dict[str, float] = {}
    lines = txt.strip().splitlines()
    seedset = set(seeds)
    for line in lines[1:]:
        c = line.split("\t")
        if len(c) < 6:
            continue
        a, b, score = c[2], c[3], c[5]
        pair = {a, b}
        if candidate not in pair:
            continue
        other = (pair - {candidate}).pop() if len(pair) == 2 else None
        if other in seedset:
            try:
                out[other] = max(out.get(other, 0.0), float(score))
            except ValueError:
                pass
    return out


def _aggregate(scores: dict[str, float], seeds: tuple[str, ...], cutoff: float) -> dict:
    vals = [scores.get(s, 0.0) for s in seeds]
    best_seed, best = "", 0.0
    for s in seeds:
        if scores.get(s, 0.0) > best:
            best, best_seed = scores[s], s
    n = len(seeds) or 1
    return {
        "max": round(max(vals) if vals else 0.0, 3),
        # max counted only if it clears the confidence cutoff (used for scoring)
        "max_thresh": round(best if best >= cutoff else 0.0, 3),
        "mean": round(sum(vals) / n, 3),
        "best_seed": best_seed if best >= cutoff else "",
        "best": round(best, 3),
        "n_seeds_ge_cutoff": sum(1 for v in vals if v >= cutoff),
        "linked_seeds": {k: round(v, 3) for k, v in sorted(scores.items(),
                                                           key=lambda kv: -kv[1])},
    }


@lru_cache(maxsize=512)
def string_scores(candidate: str, seeds: tuple[str, ...], cutoff: float = 0.4) -> dict:
    """Functional + physical STRING summary for a candidate vs the seed set."""
    func = _network_edges(candidate, seeds, physical=False)
    phys = _network_edges(candidate, seeds, physical=True)
    return {
        "gene": candidate,
        "functional": _aggregate(func, seeds, cutoff),
        "physical": _aggregate(phys, seeds, cutoff),
        "func_links": {k: round(v, 3) for k, v in func.items()},
        "phys_links": {k: round(v, 3) for k, v in phys.items()},
        "string_url": ("https://string-db.org/cgi/network?identifiers="
                       + urllib.parse.quote("\r".join([candidate, *seeds]))
                       + "&species=9606"),
        "cutoff": cutoff,
    }


@lru_cache(maxsize=128)
def network_image(candidate: str, seeds: tuple[str, ...], physical: bool = False) -> bytes | None:
    """STRING-rendered PNG of candidate + seeds (bytes), or None on failure."""
    ids = "\r".join([candidate, *seeds])
    data = {"identifiers": ids, "species": 9606, "caller_identity": _CALLER,
            "network_flavor": "confidence"}
    if physical:
        data["network_type"] = "physical"
    try:
        img = _post("image/network", data)
        return img if img[:4] == b"\x89PNG" else None
    except Exception:
        return None
