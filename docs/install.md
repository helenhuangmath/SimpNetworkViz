# Install SimpNetworkViz

## Requirements

* Python 3.9 or newer, on Linux, macOS or Windows.
* These packages, which pip installs for you: numpy, pandas, scipy, networkx, matplotlib and requests.
* For the web app only: streamlit.

Cytoscape is not required. An internet connection is needed only to download gene-set libraries by name and to query STRING; see [Network access](qa.md#which-steps-need-the-internet).

## Install from GitHub

A fresh virtual environment is recommended:

```bash
python -m venv snv-env
source snv-env/bin/activate        # Windows: snv-env\Scripts\activate
```

or with conda:

```bash
conda create -n snv python=3.12
conda activate snv
```

Then install the package:

```bash
pip install git+https://github.com/helenhuangmath/SimpNetworkViz
```

## Install from a clone

Clone the repository if you also want the examples and tests, or if you plan to change the code:

```bash
git clone https://github.com/helenhuangmath/SimpNetworkViz
cd SimpNetworkViz
pip install -e .
```

`-e` makes an editable install: changes to the source take effect without reinstalling.

## Optional extras

| Extra | Install with | What it adds |
|---|---|---|
| `gui` | `pip install -e ".[gui]"` | The point-and-click web app, [`snv gui`](commands/gui.md) |
| `test` | `pip install -e ".[test]"` | pytest, to run the test suite |

## Check the install

```bash
snv --version
snv demo
```

`snv demo` runs offline on the bundled demo data and writes its figures and tables to `demo_output/`. It should end like this:

```text
Pathway network: 17 pathways, 16 edges, 5 clusters
...
Gene network: 88 genes, 893 edges

Results written to .../demo_output/
```

`simpnetworkviz` is a longer alias for the `snv` command. If the shell cannot find `snv`, the environment you installed into is not active.

## Update or remove

```bash
pip install --upgrade git+https://github.com/helenhuangmath/SimpNetworkViz   # update
pip uninstall simpnetworkviz                                                 # remove
```

Downloaded gene-set libraries are cached in `~/.cache/simpnetworkviz`. Delete that folder to free the space or to force a fresh download.
