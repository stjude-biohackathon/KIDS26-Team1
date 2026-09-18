#!/usr/bin/env python3
"""
Shared library for the protein-complex-classifier skill.
Sequence fetching, PLM embedding generation, and biological feature extraction.
"""

import os, json, re, gc, time, math, random, warnings
import numpy as np
import pandas as pd
import requests
warnings.filterwarnings("ignore")

MAX_LEN = 1022
AMINO_ACIDS = list("ACDEFGHIKLMNPQRSTVWY")


# ══════════════════════════════════════════════════════════════════════════════
# GENE LIST + SEQUENCE FETCHING
# ══════════════════════════════════════════════════════════════════════════════
def read_genes(path):
    """Read one gene symbol per line, skipping blanks/comments. Accepts any
    non-whitespace token of letters/digits/hyphen/underscore/dot -- gene
    symbols are NOT all-uppercase-only (e.g. HGNC 'Cxxorf###' open-reading-
    frame names contain lowercase 'orf'; some composite loci use
    underscores, e.g. APOBEC3A_B). An overly strict all-caps regex here
    would silently drop legitimate genes before they ever reach sequence
    lookup."""
    with open(path) as f:
        return [l.strip() for l in f
                if l.strip() and not l.startswith("#")
                and re.match(r'^[A-Za-z0-9_.\-]+$', l.strip())]


def fetch_uniprot_sequence(gene_symbol, retries=3):
    """Map a human gene symbol -> (UniProt accession, sequence) via REST API."""
    url = "https://rest.uniprot.org/uniprotkb/search"
    params = {
        "query": f"gene_exact:{gene_symbol} AND organism_id:9606 AND reviewed:true",
        "format": "json", "size": 1, "fields": "accession,sequence",
    }
    for attempt in range(retries):
        try:
            r = requests.get(url, params=params, timeout=30)
            if r.status_code == 200:
                data = r.json()
                if data.get("results"):
                    entry = data["results"][0]
                    return (entry["primaryAccession"],
                            entry.get("sequence", {}).get("value", ""))
            elif r.status_code == 429:
                time.sleep(2 * (attempt + 1)); continue
        except Exception:
            pass
        time.sleep(1 * (attempt + 1))
    return None, None


def fetch_sequences(genes, cache_path=None):
    """Fetch (and cache) sequences for a list of genes. Returns dict."""
    if cache_path and os.path.exists(cache_path):
        with open(cache_path) as f:
            return json.load(f)
    seq_data = {}
    for i, gene in enumerate(genes, 1):
        print(f"  [{i:>3}/{len(genes)}] {gene}...", end=" ", flush=True)
        uid, seq = fetch_uniprot_sequence(gene)
        if uid and seq:
            seq_data[gene] = {"uniprot_id": uid, "sequence": seq}
            print(f"{uid} ({len(seq)} aa)")
        else:
            print("FAILED")
        if i % 15 == 0:
            time.sleep(0.5)
    if cache_path:
        with open(cache_path, "w") as f:
            json.dump(seq_data, f)
    return seq_data


# ══════════════════════════════════════════════════════════════════════════════
# PLM EMBEDDINGS
# ══════════════════════════════════════════════════════════════════════════════
def embed_esmc(seqs, gene_names, verbose=True):
    """ESM-C 600M mean-pooled embeddings (per EvolutionaryScale esm SDK)."""
    import torch
    from esm.models.esmc import ESMC
    from esm.sdk.api import ESMProtein, LogitsConfig

    device = torch.device("cpu")
    model = ESMC.from_pretrained("esmc_600m", device=device)
    model.eval()
    test = model.logits(model.encode(ESMProtein(sequence="ACDEFGHIKLMNPQRSTVWY")),
                        LogitsConfig(sequence=True, return_embeddings=True))
    dim = test.embeddings.shape[-1]

    emb = np.zeros((len(seqs), dim), dtype=np.float32)
    for i, (seq, gene) in enumerate(zip(seqs, gene_names)):
        out = model.logits(model.encode(ESMProtein(sequence=seq[:MAX_LEN])),
                           LogitsConfig(sequence=True, return_embeddings=True))
        emb[i] = out.embeddings.mean(dim=1).squeeze(0).detach().float().cpu().numpy()
        if verbose and ((i + 1) % 25 == 0 or i == 0):
            print(f"    [{i+1:>3}/{len(seqs)}] {gene}")
    del model; gc.collect()
    return emb, dim


def embed_prott5(seqs, gene_names, verbose=True):
    """ProtT5-XL-UniRef50 mean-pooled embeddings."""
    import torch
    from transformers import T5Tokenizer, T5EncoderModel
    tok = T5Tokenizer.from_pretrained("Rostlab/prot_t5_xl_uniref50", do_lower_case=False)
    model = T5EncoderModel.from_pretrained("Rostlab/prot_t5_xl_uniref50")
    model.eval()
    dim = model.config.d_model

    emb = np.zeros((len(seqs), dim), dtype=np.float32)
    for i, (seq, gene) in enumerate(zip(seqs, gene_names)):
        inp = tok(" ".join(list(seq[:MAX_LEN])), return_tensors="pt",
                  add_special_tokens=True, max_length=MAX_LEN + 2, truncation=True)
        with torch.no_grad():
            out = model(**inp)
        emb[i] = out.last_hidden_state[0, :-1, :].mean(dim=0).numpy()
        if verbose and ((i + 1) % 25 == 0 or i == 0):
            print(f"    [{i+1:>3}/{len(seqs)}] {gene}")
    del model, tok; gc.collect()
    return emb, dim


def embed_proteomelm(esmc_embeddings, verbose=True):
    """ProteomeLM-S contextualized embeddings (manual forward from ESM-C input)."""
    import torch
    import torch.nn.functional as F
    from huggingface_hub import hf_hub_download
    import safetensors.torch

    config_path = hf_hub_download("Bitbol-Lab/ProteomeLM-S", "config.json")
    weights_path = hf_hub_download("Bitbol-Lab/ProteomeLM-S", "model.safetensors")
    with open(config_path) as f:
        cfg = json.load(f)
    dim, n_heads, n_layers = cfg["dim"], cfg["n_heads"], cfg["n_layers"]
    head_size = dim // n_heads
    sd = {k: v.float() for k, v in safetensors.torch.load_file(weights_path).items()}

    def linear(x, w, b): return F.linear(x, w, b)
    def layer_norm(x, w, b, eps=1e-12): return F.layer_norm(x, (x.shape[-1],), w, b, eps)

    def block(x, idx):
        p = f"transformer.transformer.layer.{idx}"
        B, S, D = x.shape
        n = layer_norm(x, sd[f"{p}.sa_layer_norm.weight"], sd[f"{p}.sa_layer_norm.bias"])
        Q = linear(n, sd[f"{p}.attention.q_lin.weight"], sd[f"{p}.attention.q_lin.bias"]).view(B,S,n_heads,head_size).transpose(1,2)
        K = linear(n, sd[f"{p}.attention.k_lin.weight"], sd[f"{p}.attention.k_lin.bias"]).view(B,S,n_heads,head_size).transpose(1,2)
        V = linear(n, sd[f"{p}.attention.v_lin.weight"], sd[f"{p}.attention.v_lin.bias"]).view(B,S,n_heads,head_size).transpose(1,2)
        attn = F.softmax(torch.matmul(Q, K.transpose(-2,-1)) / (head_size**0.5), dim=-1)
        sa = torch.matmul(attn, V).transpose(1,2).contiguous().view(B,S,D)
        sa = linear(sa, sd[f"{p}.attention.out_lin.weight"], sd[f"{p}.attention.out_lin.bias"])
        x = x + sa
        n2 = layer_norm(x, sd[f"{p}.output_layer_norm.weight"], sd[f"{p}.output_layer_norm.bias"])
        ff = F.gelu(linear(n2, sd[f"{p}.ffn.lin1.weight"], sd[f"{p}.ffn.lin1.bias"]))
        return x + linear(ff, sd[f"{p}.ffn.lin2.weight"], sd[f"{p}.ffn.lin2.bias"])

    x = torch.tensor(esmc_embeddings, dtype=torch.float32).unsqueeze(0)
    with torch.no_grad():
        h = (linear(x, sd["embedding_main.weight"], sd["embedding_main.bias"])
           + linear(x, sd["embedding_encoder.weight"], sd["embedding_encoder.bias"]))
        for i in range(n_layers):
            h = block(h, i)
    emb = h[0].numpy()
    del sd, h, x; gc.collect()
    return emb, emb.shape[1]


def embed_proteinbert(seqs, gene_names, verbose=True):
    """ProteinBERT (Rostlab/prot_bert) mean-pooled embeddings."""
    import torch
    from transformers import BertTokenizer, BertModel
    tok = BertTokenizer.from_pretrained("Rostlab/prot_bert", do_lower_case=False)
    model = BertModel.from_pretrained("Rostlab/prot_bert")
    model.eval()
    dim = model.config.hidden_size

    emb = np.zeros((len(seqs), dim), dtype=np.float32)
    for i, (seq, gene) in enumerate(zip(seqs, gene_names)):
        inp = tok(" ".join(list(seq[:MAX_LEN])), return_tensors="pt",
                  add_special_tokens=True, max_length=MAX_LEN + 2, truncation=True)
        with torch.no_grad():
            out = model(**inp)
        emb[i] = out.last_hidden_state[0, 1:-1, :].mean(dim=0).numpy()
        if verbose and ((i + 1) % 25 == 0 or i == 0):
            print(f"    [{i+1:>3}/{len(seqs)}] {gene}")
    del model, tok; gc.collect()
    return emb, dim


# ══════════════════════════════════════════════════════════════════════════════
# BIOLOGICAL FEATURES
# ══════════════════════════════════════════════════════════════════════════════
def physicochemical_features(seq):
    from Bio.SeqUtils.ProtParam import ProteinAnalysis
    clean = "".join(c for c in seq if c in AMINO_ACIDS)
    if len(clean) < 10:
        return {}
    pa = ProteinAnalysis(clean)
    f = {"length": len(clean), "log_length": math.log10(len(clean)),
         "molecular_weight": pa.molecular_weight(), "isoelectric_point": pa.isoelectric_point(),
         "gravy": pa.gravy(), "aromaticity": pa.aromaticity(),
         "instability_index": pa.instability_index()}
    try: f["charge_at_pH7"] = pa.charge_at_pH(7.0)
    except: f["charge_at_pH7"] = 0
    ext = pa.molar_extinction_coefficient()
    f["extinction_reduced"], f["extinction_oxidized"] = ext[0], ext[1]
    ss = pa.secondary_structure_fraction()
    f["helix_fraction"], f["turn_fraction"], f["sheet_fraction"] = ss[0], ss[1], ss[2]
    aa_pct = pa.amino_acids_percent
    for aa in AMINO_ACIDS:
        f[f"aa_pct_{aa}"] = aa_pct.get(aa, 0.0)
    return f


def sequence_composition_features(seq):
    clean = "".join(c for c in seq if c in AMINO_ACIDS)
    L = len(clean); f = {}
    groups = [("polar", set("STNQ")), ("nonpolar", set("GAVILMFPW")),
              ("positive_charged", set("KRH")), ("negative_charged", set("DE")),
              ("aromatic", set("FWY")), ("tiny", set("AGCS")),
              ("small", set("AGCSDNPSTV")), ("aliphatic", set("AILV"))]
    for name, g in groups:
        f[f"frac_{name}"] = sum(1 for c in clean if c in g) / max(L, 1)
    f["net_charge_frac"] = f["frac_positive_charged"] - f["frac_negative_charged"]
    lc, i = 0, 0
    while i < L:
        j = i
        while j < L and clean[j] == clean[i]: j += 1
        if j - i >= 5: lc += j - i
        i = j
    f["low_complexity_frac"] = lc / max(L, 1)
    dis, ordr = set("AEGRQSPDK"), set("WFYILMVNCH")
    f["disorder_promoting_frac"] = sum(1 for c in clean if c in dis) / max(L, 1)
    f["order_promoting_frac"] = sum(1 for c in clean if c in ordr) / max(L, 1)
    f["disorder_order_ratio"] = f["disorder_promoting_frac"] / max(f["order_promoting_frac"], 0.01)
    key_dp = ["LL","LE","EL","EE","AL","LA","KL","LK","AA","KE","EK","VL","LV","SS","PP","GG",
              "RR","KK","DD","QQ","NN","HH","CC","WW","PG","GP","RS","SR","DE","ED","KR","RK",
              "SG","GS","AE","EA","DK","KD","LI","IL"]
    dps = {}
    for k in range(L - 1):
        dp = clean[k:k+2]; dps[dp] = dps.get(dp, 0) + 1
    for dp in key_dp:
        f[f"dp_{dp}"] = dps.get(dp, 0) / max(L - 1, 1)
    return f


CHROMATIN_DOMAINS = {"SNF2_N":"PF00176","Helicase_C":"PF00271","YEATS":"PF03366",
    "Bromodomain":"PF00439","HSA":"PF07529","Actin":"PF00022","AAA_ATPase":"PF00004",
    "SANT":"PF00249","Chromo":"PF00385","SET":"PF00856","PHD":"PF00628","RING":"PF00097",
    "ZnF_C2H2":"PF00096","WD40":"PF00400","ARID":"PF01388","HMG_box":"PF00505",
    "Histone":"PF00125","RuvB":"PF17774"}


def fetch_interpro(uid, retries=3):
    url = f"https://www.ebi.ac.uk/interpro/api/entry/all/protein/uniprot/{uid}/"
    for attempt in range(retries):
        try:
            r = requests.get(url, timeout=30, headers={"Accept": "application/json"})
            if r.status_code == 200: return r.json().get("results", [])
            elif r.status_code == 204: return []
            elif r.status_code == 429: time.sleep(2 * (attempt + 1))
        except: time.sleep(1 * (attempt + 1))
    return []


def domain_features(domains):
    f = {}
    pfam = [d["metadata"]["accession"] for d in domains
            if d.get("metadata", {}).get("source_database") == "pfam"]
    f["n_interpro_entries"] = len(domains)
    f["n_pfam_domains"] = len(pfam)
    f["n_unique_pfam"] = len(set(pfam))
    pfam_set = set(pfam)
    name_lower = " ".join((d.get("metadata", {}).get("name") or "").lower() for d in domains)
    for dname, pid in CHROMATIN_DOMAINS.items():
        f[f"has_{dname}"] = 1 if pid in pfam_set else 0
    f["has_any_chromatin_domain"] = 1 if any(f[f"has_{d}"] for d in CHROMATIN_DOMAINS) else 0
    f["has_dna_binding_kw"] = 1 if any(k in name_lower for k in ["dna-binding","dna binding","nucleic acid"]) else 0
    f["has_histone_kw"] = 1 if any(k in name_lower for k in ["histone","nucleosome","chromatin"]) else 0
    f["has_atpase_kw"] = 1 if any(k in name_lower for k in ["atpase","atp-binding","helicase"]) else 0
    return f


KEY_GO = {"chromatin_remodeling":"GO:0006338","nucleosome_binding":"GO:0031491",
    "histone_binding":"GO:0042393","ATP_dependent_chromatin":"GO:0043044",
    "helicase_activity":"GO:0004386","ATPase_activity":"GO:0016887","DNA_binding":"GO:0003677",
    "transcription_coactivator":"GO:0003713","chromatin_binding":"GO:0003682",
    "protein_complex_binding":"GO:0032403","nucleus":"GO:0005634","nucleoplasm":"GO:0005654",
    "chromatin":"GO:0000785","chromosome":"GO:0005694","histone_modification":"GO:0016570",
    "DNA_repair":"GO:0006281","transcription_regulation":"GO:0006355"}


def fetch_go_text(gene, uniprot_service, retries=3):
    """Query UniProt for a gene's GO annotations, with retry + backoff on
    transient failures. FIX (2026-09-18, ported from the load-tested
    user_data/pipeline_lib_gpu.py): the original version was a bare
    try/except with NO retry -- a single SSL error or timeout silently
    returned "" and zeroed that gene's 28 GO features with no warning.
    Empirically this can happen even in modest, non-concurrent single-run
    usage; see references/methodology.md for the case study. A broken SSL
    connection can stay broken on the same client instance, so a fresh
    client is created before each retry."""
    client = uniprot_service
    for attempt in range(retries):
        try:
            result = client.search(
                f"gene_exact:{gene} AND organism_id:9606 AND reviewed:true",
                frmt="tsv", columns="accession,go_id,go_p,go_f,go_c")
            if result and len(result.strip().split("\n")) > 1:
                parts = result.strip().split("\n")[1].split("\t")
                return " ".join(parts[1:5]) if len(parts) > 1 else ""
            return ""
        except Exception:
            if attempt < retries - 1:
                time.sleep(min(20, 2 ** attempt) + random.uniform(0, 1))
                try:
                    client = type(uniprot_service)(verbose=False)
                except Exception:
                    pass
    return ""


def go_features(go_text):
    f = {}; t = go_text.lower()
    for name, gid in KEY_GO.items():
        f[f"go_{name}"] = 1 if gid.lower() in t else 0
    kws = [("go_kw_chromatin",["chromatin","nucleosome"]),("go_kw_helicase",["helicase","atpase","atp-dependent"]),
        ("go_kw_histone",["histone"]),("go_kw_dna_repair",["dna repair","dna damage"]),
        ("go_kw_transcription",["transcription"]),("go_kw_kinase",["kinase","phosphorylation"]),
        ("go_kw_replication",["replication"]),("go_kw_splicing",["splicing","spliceosome"]),
        ("go_kw_nuclear_pore",["nuclear pore","nucleoporin"])]
    for kw_name, kwlist in kws:
        f[kw_name] = 1 if any(kw in t for kw in kwlist) else 0
    f["n_go_annotations"] = len([x for x in go_text.split(";") if x.strip()])
    f["chromatin_go_score"] = sum([f.get("go_chromatin_remodeling",0),f.get("go_nucleosome_binding",0),
        f.get("go_histone_binding",0),f.get("go_ATP_dependent_chromatin",0),
        f.get("go_chromatin_binding",0),f.get("go_histone_modification",0),
        f.get("go_kw_chromatin",0),f.get("go_kw_histone",0)])
    return f


def fetch_string(gene, retries=5):
    """Query STRING for interaction partners, with retry + backoff on both
    exceptions AND HTTP 429/5xx responses. FIX (2026-09-18, ported from the
    load-tested user_data/pipeline_lib_gpu.py): the original version only
    backed off on exceptions -- a 429 (rate-limited, a valid HTTP response,
    not an exception) fell straight through with no retry and silently
    returned [], which under sustained load caused the majority of genes to
    incorrectly show 0 STRING partners (confirmed empirically: re-fetching
    the same genes individually afterward found real partners). See
    references/methodology.md for the case study."""
    for attempt in range(retries):
        try:
            r = requests.get("https://string-db.org/api/json/interaction_partners",
                             params={"identifiers": gene, "species": 9606,
                                     "required_score": 700, "limit": 200}, timeout=30)
            if r.status_code == 200:
                return [{"partner": e.get("preferredName_B",""), "score": e.get("score",0)}
                        for e in r.json()]
            elif r.status_code == 429 or r.status_code >= 500:
                time.sleep(min(30, 2 ** attempt) + random.uniform(0, 1))
            else:
                return []
        except Exception:
            time.sleep(min(30, 2 ** attempt) + random.uniform(0, 1))
    return []


def ppi_features(gene, partners, known_positives, string_cache):
    partner_names = set(p["partner"] for p in partners)
    scores = [p["score"] for p in partners]
    f = {}
    f["n_string_partners"] = len(partners)
    f["mean_string_score"] = float(np.mean(scores)) if scores else 0
    f["max_string_score"] = max(scores) if scores else 0
    shared = partner_names & known_positives
    f["n_shared_with_positives"] = len(shared)
    f["frac_shared_with_positives"] = len(shared) / max(len(partner_names), 1)
    n_links = sum(1 for pg in known_positives
                  if gene in set(p["partner"] for p in string_cache.get(pg, [])))
    f["n_positive_genes_linking_to_this"] = n_links
    f["frac_positives_linking"] = n_links / max(len(known_positives), 1)
    f["n_high_conf_partners"] = sum(1 for s in scores if s >= 0.9)
    chromatin = {"SMARCA4","SMARCA2","SMARCC1","SMARCC2","ARID1A","SMARCB1","BRG1","CHD1","CHD4","EP400","TRRAP"}
    f["n_chromatin_remodeler_links"] = len(partner_names & chromatin)
    return f
