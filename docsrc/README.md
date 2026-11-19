# Documentation sources

This folder holds the sources of the documentation site, <https://helenhuangmath.github.io/SimpNetworkViz/>. The pages are Markdown, built with Sphinx and the Read the Docs theme.

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
| `index.md` | Landing page |
| `install.md`, `workflow.md`, `fileformats.md`, `gallery.md`, `theory.md`, `qa.md` | Top-level pages |
| `commands/` | One page per `snv` subcommand. The option lists are generated from the parser in `src/simpnetworkviz/cli.py`. |
| `guide/` | The whole-library view and pathway map, grouping, figures, Python API |
| `tutorial/` | The demo walkthrough and the Cl13 vs Arm tutorial |
| `api.md` | API reference, generated from the docstrings |
| `images/tutorial/` | PNG previews of the tutorial's PDF figures |
| `_extra/` | Files copied to the site as they are: the interactive example maps |

Most pages started as sections of the top-level `README.md`, and `tutorial/cl13_vs_arm.md` as `examples/Cl13_vs_Arm/TUTORIAL.md`. When you change one of those, change the matching page here too.
