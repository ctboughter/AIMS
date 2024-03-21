# Configuration file for the Sphinx documentation builder.

# -- Project information

project = 'AIMS'
copyright = '2024, Boughter'
author = 'Boughter'

release = '0.9'
version = 'AIMS_0.9'

# -- General configuration

# At some point, write a CSS file to change colors
#def setup(app):
#    app.add_css_file('custom.css')

extensions = [
    'sphinx.ext.duration',
    'sphinx.ext.doctest',
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.intersphinx',
]

intersphinx_mapping = {
    'python': ('https://docs.python.org/3/', None),
    'sphinx': ('https://www.sphinx-doc.org/en/master/', None),
}
intersphinx_disabled_domains = ['std']

templates_path = ['_templates']

html_css_files = [
    "style.css",
]

# -- Options for HTML output

html_theme = 'sphinx_rtd_theme'

# -- Options for EPUB output
epub_show_urls = 'footnote'
