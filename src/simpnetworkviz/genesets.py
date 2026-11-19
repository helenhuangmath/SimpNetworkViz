"""Loading gene-set libraries (GMT files, Enrichr libraries, bundled demo)."""
from __future__ import annotations

import os
from importlib import resources
from pathlib import Path
from typing import Dict, Iterable, List, Optional

import requests

GeneSets = Dict[str, List[str]]

ENRICHR_URL = "https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName={name}"

#: A curated selection of widely used Enrichr libraries. Any Enrichr library
#: name works; see https://maayanlab.cloud/Enrichr/#libraries for the full list.
POPULAR_LIBRARIES = {
    "MSigDB_Hallmark_2020": "MSigDB Hallmark gene sets (50 broad processes)",
    "KEGG_2021_Human": "KEGG pathways (human)",
    "Reactome_2022": "Reactome pathways (human)",
    "WikiPathways_2023_Human": "WikiPathways (human)",
    "GO_Biological_Process_2023": "Gene Ontology biological process",
    "GO_Molecular_Function_2023": "Gene Ontology molecular function",
    "GO_Cellular_Component_2023": "Gene Ontology cellular component",
    "KEGG_2019_Mouse": "KEGG pathways (mouse)",
    "WikiPathways_2019_Mouse": "WikiPathways (mouse)",
    "BioPlanet_2019": "NCATS BioPlanet pathways",
}

DEMO = "demo"

# ---------------------------------------------------------------------------
# Species and MSigDB collections
# ---------------------------------------------------------------------------

#: Supported species: MSigDB code, NCBI taxon (for STRING), Enrichr naming hint.
SPECIES = {
    "human": {"msigdb": "Hs", "taxon": 9606, "aliases": ("human", "hs", "homo sapiens", "9606")},
    "mouse": {"msigdb": "Mm", "taxon": 10090, "aliases": ("mouse", "mm", "mus musculus", "10090")},
}

MSIGDB_RELEASE = "2025.1"
MSIGDB_URL = "https://data.broadinstitute.org/gsea-msigdb/msigdb/release/{release}.{code}/{file}.v{release}.{code}.symbols.gmt"

#: MSigDB collections by name: (human file stem, mouse file stem); ``None`` = not available for that species.
MSIGDB_COLLECTIONS = {
    "Hallmark": ("h.all", "mh.all"),
    "GO": ("c5.go", "m5.go"),
    "GO_BP": ("c5.go.bp", "m5.go.bp"),
    "GO_CC": ("c5.go.cc", "m5.go.cc"),
    "GO_MF": ("c5.go.mf", "m5.go.mf"),
    "Reactome": ("c2.cp.reactome", "m2.cp.reactome"),
    "WikiPathways": ("c2.cp.wikipathways", "m2.cp.wikipathways"),
    "BioCarta": ("c2.cp.biocarta", "m2.cp.biocarta"),
    "KEGG": ("c2.cp.kegg_medicus", None),
    "KEGG_legacy": ("c2.cp.kegg_legacy", None),
    "PID": ("c2.cp.pid", None),
    "CP": ("c2.cp", "m2.cp"),
    "CGP": ("c2.cgp", "m2.cgp"),
    "TFT": ("c3.tft", "m3.gtrd"),
    "miRNA": ("c3.mir", "m3.mirdb"),
    "HPO": ("c5.hpo", None),
    "TumorPhenotype": (None, "m5.mpt"),
    "Oncogenic": ("c6.all", None),
    "ImmuneSigDB": ("c7.immunesigdb", None),
    "Vaccine": ("c7.vax", None),
    "CellType": ("c8.all", "m8.all"),
    "Positional": ("c1.all", "m1.all"),
    "All": ("msigdb", "msigdb"),
}


def normalize_species(species) -> str:
    """``"human"`` or ``"mouse"`` from a name, MSigDB code or NCBI taxon id."""
    key = str(species).strip().lower()
    for name, info in SPECIES.items():
        if key == name or key in info["aliases"]:
            return name
    raise ValueError(f"Unknown species {species!r}: use 'human' or 'mouse' (or the taxon id 9606 / 10090)")


def species_taxon(species) -> int:
    """NCBI taxon id of a species name (for STRING)."""
    try:
        return int(species)  # already a taxon id
    except (TypeError, ValueError):
        return SPECIES[normalize_species(species)]["taxon"]


def msigdb_file(collection: str, species: str = "human", release: str = MSIGDB_RELEASE) -> str:
    """File stem of an MSigDB collection for a species, e.g. ``GO_BP`` + mouse -> ``m5.go.bp.v2025.1.Mm.symbols``."""
    sp = normalize_species(species)
    names = {k.lower(): k for k in MSIGDB_COLLECTIONS}
    key = names.get(collection.strip().lower())
    if key is None:
        raise ValueError(f"Unknown MSigDB collection {collection!r}. Choose from: {', '.join(MSIGDB_COLLECTIONS)}")
    stem = MSIGDB_COLLECTIONS[key][0 if sp == "human" else 1]
    if stem is None:
        avail = [k for k, v in MSIGDB_COLLECTIONS.items() if v[0 if sp == "human" else 1]]
        raise ValueError(f"MSigDB has no {key} collection for {sp}. Available for {sp}: {', '.join(avail)}")
    code = SPECIES[sp]["msigdb"]
    return f"{stem}.v{release}.{code}.symbols"


def fetch_msigdb(collection: str, species: str = "human", release: str = MSIGDB_RELEASE, force: bool = False,
                 timeout: float = 120) -> GeneSets:
    """Download an MSigDB collection (gene symbols) for ``species`` from the Broad release server.

    ``collection`` is a key of :data:`MSIGDB_COLLECTIONS` (``Hallmark``, ``GO_BP``,
    ``Reactome``, ``KEGG`` ...); mouse collections are MSigDB's native mouse
    sets (``mh``, ``m2``, ``m5`` ...), not orthology-mapped human ones. Files are
    cached under ``~/.cache/simpnetworkviz/msigdb``.
    """
    stem = msigdb_file(collection, species, release)
    cache = _cache_dir() / "msigdb" / f"{stem}.gmt"
    if cache.exists() and not force:
        return read_gmt(cache)
    sp = normalize_species(species)
    url = MSIGDB_URL.format(release=release, code=SPECIES[sp]["msigdb"], file=stem.split(".v")[0])
    resp = requests.get(url, timeout=timeout)
    if resp.status_code != 200:
        raise ValueError(f"MSigDB download failed ({resp.status_code}) for {url}. Check the release ({release}) "
                         "or download the GMT from https://www.gsea-msigdb.org/gsea/msigdb and pass it with --gmt.")
    gene_sets = parse_gmt_lines(resp.text.splitlines())
    if not gene_sets:
        raise ValueError(f"No gene sets in {url}")
    cache.parent.mkdir(parents=True, exist_ok=True)
    write_gmt(gene_sets, cache, description=f"MSigDB {collection} {sp} {release}")
    return gene_sets


def _cache_dir() -> Path:
    base = os.environ.get("SIMPNETWORKVIZ_CACHE") or os.path.join(
        os.environ.get("XDG_CACHE_HOME", os.path.expanduser("~/.cache")), "simpnetworkviz"
    )
    path = Path(base)
    path.mkdir(parents=True, exist_ok=True)
    return path


def normalize_gene(gene: str) -> str:
    """Canonical form used for matching: stripped, upper case, weights removed."""
    return gene.strip().split(",")[0].strip().upper()


def parse_gmt_lines(lines: Iterable[str]) -> GeneSets:
    """Parse GMT-formatted lines: ``name<TAB>description<TAB>gene1<TAB>gene2...``.

    Also accepts Enrichr's text format where genes may carry a ``,weight`` suffix.
    """
    gene_sets: GeneSets = {}
    for line in lines:
        line = line.rstrip("\r\n")
        if not line.strip() or line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) < 3:
            continue
        name = fields[0].strip()
        genes = []
        seen = set()
        for g in fields[2:]:
            g = normalize_gene(g)
            if g and g not in seen:
                seen.add(g)
                genes.append(g)
        if genes:
            gene_sets[name] = genes
    return gene_sets


def read_gmt(path) -> GeneSets:
    with open(path, encoding="utf-8") as fh:
        return parse_gmt_lines(fh)


def write_gmt(gene_sets: GeneSets, path, description: str = "") -> None:
    with open(path, "w", encoding="utf-8") as fh:
        for name, genes in gene_sets.items():
            fh.write("\t".join([name, description] + list(genes)) + "\n")


def fetch_enrichr_library(name: str, force: bool = False, timeout: float = 60) -> GeneSets:
    """Download an Enrichr gene-set library (cached under ``~/.cache/simpnetworkviz``)."""
    cache = _cache_dir() / f"{name}.gmt"
    if cache.exists() and not force:
        return read_gmt(cache)
    resp = requests.get(ENRICHR_URL.format(name=name), timeout=timeout)
    resp.raise_for_status()
    gene_sets = parse_gmt_lines(resp.text.splitlines())
    if not gene_sets:
        raise ValueError(
            f"Enrichr returned no gene sets for library {name!r}. "
            "Check the name at https://maayanlab.cloud/Enrichr/#libraries"
        )
    write_gmt(gene_sets, cache, description=name)
    return gene_sets


def demo_path(filename: str) -> str:
    return str(resources.files("simpnetworkviz").joinpath("data", filename))


def load_gene_sets(source, min_size: int = 1, max_size: Optional[int] = None, species: str = "human") -> GeneSets:
    """Load gene sets from a GMT path, an Enrichr library name, an MSigDB collection, ``"demo"`` or a dict.

    ``"msigdb:GO_BP"`` (or ``"msigdb:Hallmark,msigdb:Reactome"``) downloads the
    MSigDB collection for ``species`` (``"human"`` or ``"mouse"``; see
    :data:`MSIGDB_COLLECTIONS` and :func:`fetch_msigdb`). Several sources may be
    combined by passing a list or a comma separated string, e.g.
    ``"KEGG_2021_Human,Reactome_2022"``. Pathway names are prefixed with the
    library name only when names collide.
    """
    if isinstance(source, dict):
        gene_sets = {k: [normalize_gene(g) for g in v] for k, v in source.items()}
    else:
        if isinstance(source, str):
            sources = [s.strip() for s in source.split(",") if s.strip()]
        else:
            sources = [str(s) for s in source]
        gene_sets = {}
        for src in sources:
            if src == DEMO:
                loaded = read_gmt(demo_path("demo_genesets.gmt"))
            elif src.lower().startswith("msigdb:"):
                loaded = fetch_msigdb(src.split(":", 1)[1], species=species)
            elif os.path.exists(src):
                loaded = read_gmt(src)
            else:
                loaded = fetch_enrichr_library(src)
            for k, v in loaded.items():
                key = k if k not in gene_sets else f"{src}: {k}"
                gene_sets[key] = v
    return {
        k: v
        for k, v in gene_sets.items()
        if len(v) >= min_size and (max_size is None or len(v) <= max_size)
    }


def read_gene_list(path) -> Dict[str, Optional[float]]:
    """Read genes from a text/CSV/TSV file.

    Accepts one gene per line, or a table whose first column is the gene and an
    optional numeric column (e.g. log2 fold change). Returns ``{gene: value}``.
    """
    import pandas as pd

    with open(path, encoding="utf-8") as fh:
        text = fh.read()
    lines = [l for l in text.splitlines() if l.strip() and not l.startswith("#")]
    if not lines:
        return {}
    sep = "\t" if "\t" in lines[0] else ("," if "," in lines[0] else None)
    if sep is None:
        genes: Dict[str, Optional[float]] = {}
        for l in lines:
            for tok in l.replace(";", " ").split():
                genes.setdefault(normalize_gene(tok), None)
        return genes
    df = pd.read_csv(path, sep=sep, header=None, comment="#", dtype=str)
    first = [str(x).strip() for x in df.iloc[0]]
    has_header = first[0].lower() in _HEADER_WORDS or (
        df.shape[1] > 1 and not _is_number(first[1]) and len(df) > 1 and _is_number(df.iloc[1, 1])
    )
    if has_header:
        df.columns = first
        df = df.iloc[1:].reset_index(drop=True)
    for c in df.columns[1:]:
        converted = pd.to_numeric(df[c], errors="coerce")
        if converted.notna().any():
            df[c] = converted
    return gene_values_from_frame(df)


_HEADER_WORDS = {"gene", "genes", "symbol", "gene_symbol", "gene_name", "genesymbol", "id", "name", "gene_id"}


def _is_number(x) -> bool:
    try:
        float(x)
        return True
    except (TypeError, ValueError):
        return False


def gene_values_from_frame(df) -> Dict[str, Optional[float]]:
    import pandas as pd

    gene_col = df.columns[0]
    value_col = None
    preferred = ["log2fc", "logfc", "log2foldchange", "log2_fold_change", "avg_log2fc", "fc", "score", "value", "t", "stat"]
    lower = {str(c).lower(): c for c in df.columns[1:]}
    for p in preferred:
        if p in lower:
            value_col = lower[p]
            break
    if value_col is None:
        for c in df.columns[1:]:
            if pd.api.types.is_numeric_dtype(df[c]):
                value_col = c
                break
    out: Dict[str, Optional[float]] = {}
    for _, row in df.iterrows():
        g = row[gene_col]
        if pd.isna(g):
            continue
        v = None if value_col is None or pd.isna(row[value_col]) else float(row[value_col])
        out.setdefault(normalize_gene(str(g)), v)
    return out
