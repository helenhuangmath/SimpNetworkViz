# Contributing

Questions, bug reports and pull requests are welcome.

## Ask a question or report a bug

Open an [issue](https://github.com/helenhuangmath/SimpNetworkViz/issues). For a bug, please include:

* the command or code you ran
* the full error message
* the output of `snv --version` and `python --version`
* a small input file that reproduces the problem, if you can share one

## Set up for development

```bash
git clone https://github.com/helenhuangmath/SimpNetworkViz
cd SimpNetworkViz
pip install -e ".[test,gui]"
pytest -q
```

Tests run on GitHub Actions with Python 3.9 and 3.12 on every push and pull request.

## Make a change

1. Create a branch.
2. Make the change, and add or update a test in `tests/` for any change in behavior.
3. Run `pytest -q`.
4. If the change affects what users see, update the README and the matching page in `docs/`.
5. Open a pull request that says what changed and why.

## Where things are

| Part | File |
|---|---|
| Gene sets and gene lists | `src/simpnetworkviz/genesets.py` |
| Enrichment test, reading result tables | `src/simpnetworkviz/enrichment.py` |
| Similarity, clustering, labels, gene networks | `src/simpnetworkviz/network.py` |
| Layout | `src/simpnetworkviz/layout.py` |
| Colors | `src/simpnetworkviz/style.py` |
| Figures | `src/simpnetworkviz/plotting.py` |
| HTML, Cytoscape, GraphML and CSV export | `src/simpnetworkviz/export.py` |
| One-call workflow | `src/simpnetworkviz/pipeline.py` |
| Command line | `src/simpnetworkviz/cli.py` |
| Web app | `src/simpnetworkviz/app.py` |

## Build this documentation

The documentation is written in Markdown and built with Sphinx.

```bash
pip install -e . -r docs/requirements.txt
sphinx-build -b html docs docs/_build/html
```

Open `docs/_build/html/index.html` in a browser to review the generated documentation.

## Versioning

The package follows [semantic versioning](https://semver.org/).

## License

By contributing, you agree that your contribution is released under the project's [MIT License](https://github.com/helenhuangmath/SimpNetworkViz/blob/main/LICENSE).
