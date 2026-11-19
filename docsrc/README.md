# Documentation sources

This folder holds the sources of the documentation site, <https://helenhuangmath.github.io/SimpNetworkViz/>. The pages are Markdown, built with Sphinx and the PyData theme. The front page is a hand-written HTML page.

GitHub Pages serves the `docs/` folder of the `gh-pages` branch, so the built site is committed there.

## Rebuild the site

```bash
pip install -e . -r docsrc/requirements.txt
bash docsrc/build.sh
```

`build.sh` builds into a temporary folder and copies the result into `docs/`. It leaves `docs/images/` alone. Open `docs/index.html` to check the result, then commit `docsrc/` and `docs/` to the `gh-pages` branch and push.

## Where things are

| Path | Contents |
|---|---|
| `_extra/index.html`, `_extra/styles.css` | The front page of the site, hand-written HTML and CSS |
| `documentation.md` | Home page of the documentation |
| `install.md`, `workflow.md`, `fileformats.md`, `gallery.md`, `theory.md`, `qa.md` | Top-level pages |
| `commands/` | One page per `snv` subcommand. The option lists are generated from the parser in `src/simpnetworkviz/cli.py`. |
| `guide/` | The whole-library view and pathway map, grouping, figures, Python API |
| `tutorial/` | The step-by-step tutorials; `tutorial/index.md` is the "Get started" page |
| `api.md` | Reference of the Python API, generated from the docstrings |
| `images/tutorial/` | PNG previews of the Cl13 vs Arm tutorial's PDF figures |
| `images/vignettes/` | Figures of the other tutorials, made by running their code on the demo data |
| `_extra/examples/` | Interactive example maps, copied to the site as they are |

Most pages started as sections of the top-level `README.md`, and `tutorial/cl13_vs_arm.md` as `examples/Cl13_vs_Arm/TUTORIAL.md`. When you change one of those, change the matching page here too.

## The PathwayAgent tab

The navigation bar has a tab that links to the companion package PathwayAgent. Its address is set in two places: `PATHWAYAGENT_URL` in `conf.py` for the documentation pages, and the `PathwayAgent` link in `_extra/index.html` for the front page. Change both when PathwayAgent gets its own site.
