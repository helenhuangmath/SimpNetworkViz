"""Sphinx configuration for the SimpNetworkViz documentation.

Build locally with:

    pip install -e . -r docs/requirements.txt
    sphinx-build -b html docs docs/_build/html
"""
from simpnetworkviz import __version__

project = "SimpNetworkViz"
author = "Hua Huang"
copyright = "2026, Hua Huang"

release = __version__
version = release

extensions = [
    "myst_parser",
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "sphinx_copybutton",
]

source_suffix = {".md": "markdown"}
root_doc = "index"
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

# MyST (Markdown) options
myst_enable_extensions = ["colon_fence", "deflist"]
myst_heading_anchors = 3

# API reference
autodoc_member_order = "bysource"
autodoc_typehints = "none"
napoleon_google_docstring = False
napoleon_numpy_docstring = True

# Copy button: strip shell and Python prompts
copybutton_prompt_text = r">>> |\.\.\. |\$ "
copybutton_prompt_is_regexp = True

# HTML output
html_theme = "sphinx_rtd_theme"
html_title = f"SimpNetworkViz {release} documentation"
html_static_path = ["_static"]
html_css_files = ["custom.css"]
html_theme_options = {
    "navigation_depth": 3,
    "collapse_navigation": False,
    "style_external_links": True,
}
html_context = {
    "display_github": True,
    "github_user": "helenhuangmath",
    "github_repo": "SimpNetworkViz",
    "github_version": "main",
    "conf_py_path": "/docs/",
}
