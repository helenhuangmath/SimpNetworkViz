# Install

Python 3.9 or newer is required.

```bash
git clone https://github.com/helenhuangmath/SimpNetworkViz
cd SimpNetworkViz
pip install -e .             # command line + Python API
pip install -e ".[gui]"      # optional: the point-and-click app ([The point-and-click app](commands/gui.md))
pip install -e ".[map]"      # optional: UMAP/t-SNE for the pathway map ([The whole-library view and the pathway map](guide/library_map.md))
pip install -e ".[test]"     # optional: run the tests with  pytest -q
snv --version
```

Dependencies are numpy, pandas, scipy, networkx, matplotlib and requests; they are installed automatically. Internet is needed only to download a library the first time and for STRING interactions. The demo runs offline.

## Check the install

```bash
snv --version
snv demo          # runs offline, writes figures and tables to demo_output/
```

If `snv` is not found after installing, use `python -m simpnetworkviz.cli` wherever this documentation says `snv`.

## Optional extras

| Extra | Install with | What it adds |
|---|---|---|
| `gui` | `pip install -e ".[gui]"` | The point-and-click app, [`snv gui`](commands/gui.md) |
| `map` | `pip install -e ".[map]"` | UMAP and t-SNE for the [pathway map](guide/library_map.md) |
| `test` | `pip install -e ".[test]"` | pytest and the reference implementations used by the validation tests |

## Downloaded files

Enrichr libraries are cached in `~/.cache/simpnetworkviz`, and MSigDB collections in `~/.cache/simpnetworkviz/msigdb`. Delete these folders to free the space or to force a fresh download.
