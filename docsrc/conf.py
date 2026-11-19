"""Sphinx configuration for the SimpNetworkViz documentation.

The sources are in this folder (docsrc/). The built site is committed to docs/ on the
gh-pages branch, which GitHub Pages serves. Rebuild with:

    pip install -e . -r docsrc/requirements.txt
    bash docsrc/build.sh
"""
import os
import sys
from importlib.metadata import PackageNotFoundError, version as _version

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))  # for _cli_for_docs.py

project = "SimpNetworkViz"
author = "Hua Huang"
copyright = "2026, Hua Huang"

try:
    release = _version("simpnetworkviz")
except PackageNotFoundError:  # package not installed: the API and option tables will be empty
    release = "1.0.0"
version = release

extensions = [
    "myst_parser",
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "sphinx_copybutton",
    "sphinxarg.ext",
    "sphinx_design",
]

source_suffix = {".md": "markdown"}
# The site root (index.html) is the hand-written landing page in _extra/.
# The documentation starts at documentation.html.
root_doc = "documentation"
exclude_patterns = ["_build", "_extra", "README.md", "Thumbs.db", ".DS_Store"]

# MyST (Markdown) options
myst_enable_extensions = ["colon_fence", "deflist"]
myst_heading_anchors = 3

# Keep "--option" as typed instead of turning the two hyphens into a dash
smartquotes = False

# API reference
autodoc_member_order = "bysource"
autodoc_typehints = "none"
napoleon_google_docstring = False
napoleon_numpy_docstring = True

# Copy button: strip shell and Python prompts
copybutton_prompt_text = r">>> |\.\.\. |\$ "
copybutton_prompt_is_regexp = True

# HTML output
SITE = "https://helenhuangmath.github.io/SimpNetworkViz/"
html_theme = "pydata_sphinx_theme"
html_title = "SimpNetworkViz"
html_static_path = ["_static"]
html_css_files = ["custom.css"]
# Files copied to the site root as they are: the landing page (index.html, styles.css),
# the interactive example maps and .nojekyll
html_extra_path = ["_extra"]
html_show_sourcelink = False
html_theme_options = {
    "logo": {"text": "SimpNetworkViz", "link": SITE},
    "header_links_before_dropdown": 5,
    "navbar_align": "left",
    "navbar_end": ["theme-switcher", "navbar-icon-links"],
    "icon_links": [
        {"name": "GitHub", "url": "https://github.com/helenhuangmath/SimpNetworkViz", "icon": "fa-brands fa-github"},
    ],
    "secondary_sidebar_items": ["page-toc"],
    "show_toc_level": 2,
    "navigation_with_keys": False,
    "footer_start": ["copyright"],
    "footer_end": ["sphinx-version"],
}
html_context = {"default_mode": "light"}


def _tidy_docstring(app, what, name, obj, options, lines):
    """Make two docstrings valid reStructuredText without touching the package source."""
    for i, line in enumerate(lines):
        lines[i] = line.replace("|log2FC|", "abs(log2FC)")
    if name.endswith("build_gene_network"):
        out = []
        for line in lines:
            text = line.strip()
            if text == "mode":
                out += ["``mode`` is one of:", ""]
            elif line.startswith("    ``"):          # an option of the list: ``"shared"``  description
                out.append("* " + text)
            elif line.startswith("          "):      # continuation of the description
                out.append("  " + text)
            else:
                out.append(line)
        lines[:] = out


def setup(app):
    app.connect("autodoc-process-docstring", _tidy_docstring)
