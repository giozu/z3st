# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
# Z3ST: An open-source FEniCSx framework for thermo-mechanical analysis
# Author: Giovanni Zullo
# Version: 0.1.0 (2025)
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---

# ============================================================================
# IMPORTS AND PATHS
# ============================================================================

import os
import sys

# Add project paths so autodoc finds your Python modules
sys.path.insert(0, os.path.abspath("../.."))
sys.path.insert(0, os.path.abspath("../../core"))
sys.path.insert(0, os.path.abspath("../../materials"))
sys.path.insert(0, os.path.abspath("../../models"))
sys.path.insert(0, os.path.abspath("../../utils"))

# ============================================================================
# PROJECT INFORMATION
# ============================================================================

project = "Z3ST"
copyright = "2025, Giovanni Zullo"
author = "Giovanni Zullo"
release = "0.1.0"

# ============================================================================
# GENERAL CONFIGURATION
# ============================================================================

extensions = [
    # Core Sphinx extensions
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.napoleon",
    "sphinx.ext.mathjax",
    "sphinx.ext.viewcode",
    "sphinx.ext.todo",
    "sphinx.ext.ifconfig",
    "sphinx.ext.imgmath",

    # Type hints and doc improvements
    "sphinx_autodoc_typehints",
    "myst_parser",  # enables Markdown support
]

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

language = "en"

# ============================================================================
# AUTODOC & NAPOLEON SETTINGS
# ============================================================================

autodoc_member_order = "bysource"  # preserves order of code
autodoc_typehints = "description"  # show type hints in function docs
autodoc_default_options = {
    "members": True,
    "undoc-members": True,
    "show-inheritance": True,
    "private-members": False,
}

autosummary_generate = True
napoleon_google_docstring = True
napoleon_numpy_docstring = True
napoleon_include_init_with_doc = True
napoleon_include_special_with_doc = True
add_module_names = False

# ============================================================================
# HTML OUTPUT
# ============================================================================

# You can use either sphinx_book_theme (modern) or sphinx_rtd_theme (classic)
# If you prefer the modern style, uncomment the book theme lines below.

# html_theme = "sphinx_book_theme"
# html_theme_options = {
#     "repository_url": "https://github.com/giozu/z3st",
#     "use_repository_button": True,
#     "use_download_button": False,
#     "use_fullscreen_button": True,
#     "path_to_docs": "docs/",
#     "use_sidenotes": True,
# }

# ============================================================================
# HTML OUTPUT (modern theme)
# ============================================================================

html_theme = "sphinx_book_theme"

html_theme_options = {
    "repository_url": "https://github.com/giozu/z3st",
    "use_repository_button": True,
    "use_issues_button": True,
    "use_edit_page_button": False,
    "use_fullscreen_button": True,
    "use_download_button": True,
    "path_to_docs": "docs/",
    "home_page_in_toc": True,
    "show_toc_level": 2,
    "use_sidenotes": True,
    "show_navbar_depth": 2,
    "launch_buttons": {
        "binderhub_url": "",
    },
}

html_title = "z3st documentation"
# html_static_path = ["_static"]

# ============================================================================
# MATH AND TODO OPTIONS
# ============================================================================

todo_include_todos = True
math_number_all = True
math_eqref_format = "Eq.{number}"

# ============================================================================
# MARKDOWN (MyST) CONFIGURATION
# ============================================================================

myst_enable_extensions = [
    "colon_fence",        # for ::: fenced blocks
    "deflist",            # for definition lists
    "html_admonition",    # for notes, warnings, etc.
    "html_image",         # for Markdown images
    "substitution",       # for variable substitution
    "replacements",
]
myst_heading_anchors = 3  # automatic anchors up to H3
