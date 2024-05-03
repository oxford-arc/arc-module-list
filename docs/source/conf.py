# Configuration file for the Sphinx documentation builder.

# -- Project information

project = 'arc-module-list'

author = 'The ARC Team'

release = '0.1'
version = '0.1.0'

# -- General configuration

extensions = [
    'sphinx.ext.duration',
    'sphinx.ext.doctest',
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.intersphinx',
#    "sphinx_favicon",
]

intersphinx_mapping = {
    'python': ('https://docs.python.org/3/', None),
    'sphinx': ('https://www.sphinx-doc.org/en/master/', None),
}
intersphinx_disabled_domains = ['std']

templates_path = ['_templates']
pygments_style = 'sphinx'

# -- Options for HTML output
#favicons = [
#    "favicon16.png",
#    "favicon32.png",
#    "favicon96.png",
#    "favicon160.png",
#    "arc_icon.svg",
#]

html_theme = 'sphinx_rtd_theme'
html_favicon = 'favicon.ico'

# -- Options for EPUB output
epub_show_urls = 'footnote'
