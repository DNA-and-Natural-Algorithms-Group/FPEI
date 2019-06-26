# -*- coding: utf-8 -*-
#
# Multistrand documentation build configuration file, created by
# sphinx-quickstart on Tue Sep 28 22:48:53 2010.
#
# This file is execfile()d with the current directory set to its containing dir.
#
# Note that not all possible configuration values are present in this
# autogenerated file.
#
# All configuration values have a default; values that are commented out
# serve to show the default.

import sys, os

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#sys.path.insert(0, os.path.abspath('.'))

# by including the following function, we can add some event handlers to modify docstrings.
def setup(app):

    app.connect('autodoc-process-docstring',proc_docstring)

def proc_docstring(app, what, name, obj, options, lines):
    """
    Processes docstrings to make them readable both by pydoc and by sphinx.

    Examples below, for this function. Blocks are marked with exactly
    2 spaces extra indent in this example
    
    Keyword Arguments:
    app             -- the Sphinx application object
    what [type=str] -- the type of the object which the docstring belongs to (one of 'module', 'class', 'exception', 'function', 'method', 'attribute')
    name [type=str] -- fully qualified name of object
    obj             -- the object itself
    options         -- the options given to the directive: an object with
                       attributes inherited_members, undoc_members,
                       show_inheritance and noindex that are true if the flag
                       option of same name was given to the auto directive
    lines [type=list] -- the lines of the docstring

    Return Value: None

    Results [omitting opening section which is unchanged]:
      :param app: the Sphinx application object
      :param what: the type of the object which the docstring belongs to (one of 'module', 'class', 'exception', 'function', 'method', 'attribute')
      :type what: str
      :param name: fully qualified name of object
      :type name: str
      :param obj: the object itself
      :param options: the options given to the directive: an object with attributes inherited_members, undoc_members, show_inheritance and noindex that are true if the flag option of same name was given to the auto directive
      :param lines: the lines of the docstring
      :type lines: list
      
      :rtype: None
    """
    if what != 'function' and what != 'method' and what != 'class' and what != 'attribute':
        return
    
    import re
    args_loc = -1
    retval_loc = -1
    for l in lines:
        k = l.lower()
        if k.startswith("arguments:") or k.startswith("keyword arguments:"):
            args_loc = lines.index(l)
        if k.startswith("return value:"):
            retval_loc = lines.index(l)
            retval_start = len(k) - 13

    replacedata = []
    # list of index, replacevalue pairs, where replacevalue is a
    # dictionary with the following items:
    # "argname":str
    # 'argtext':str
    # 'argkeywords': list of two items, the first is the 'default'
    #                value, second is the type.  Note that these are
    #                either strings or None
    def add_to_val( idx, currentval, line ):
        #helper to add more values to a specific line, following our parse rules.
        # Returns True if it thinks we need to stop, False otherwise.
        
        parsed_data = re.match(r'^\s*(.*?)\s*((?:\[(.*)\]\s*)?--\s*(.*?)\s*)?$',line)
        # note: if we don't match the above, parsed_data will be
        # None. But this shouldn't happen, as the above can match
        # /any/ string including a blank one, due to the conditionals
        # and *. Note that a blank string will have group 0 == '', as
        # that non-greedy group is the weakest possible match.
        groups = parsed_data.groups()

        if len(groups[0]) == 0 and groups[1] == None:
            if len(currentval) > 0:
                replacedata.append(currentval.copy())
                currentval.clear()
            return True

        if groups[1] == None:
            currentval['argtext'] = currentval['argtext'].strip() + ' ' + groups[0]
            currentval['range'] = [currentval['range'][0]] + [idx]
            return False

        #If we were already processing a line and get a new starter, add the old one.
        if len(currentval) > 0:
            replacedata.append(currentval.copy())

        #process the two parts. Left is the argument name and data
        #about it, right side is all the text.
        currentval.clear()
        currentval['argname'] = groups[0]
        currentval['argtext'] = groups[3]
        currentval['argkeywords'] = [None,None]
        currentval['range'] = [idx]
        if groups[2] != None:
            keypairs = [i.split('=') for i in groups[2].split(',')]
            for item in keypairs:
                stripped_name = item[0].strip().lower()
                if stripped_name == 'default':
                    currentval['argkeywords'][0] = item[1].strip()
                elif stripped_name == 'type':
                    currentval['argkeywords'][1] = item[1].strip()
        # done with processing this line
        return False

    currentval = {}
    if args_loc >= 0:
        for i in range(args_loc+1,len(lines)):
            flag = add_to_val( i, currentval, lines[i] )
            if flag: break
    if len(currentval) > 0:
        replacedata.append( currentval.copy() )
    #now we grab the return value info. 
    if retval_loc >= 0:
        currentval.clear()
        if retval_start > 0:
            add_to_val( retval_loc, currentval, '-- ' + lines[retval_loc][13:])
        for i in range(retval_loc+1,len(lines)):
            flag = add_to_val( i, currentval, lines[i] )
            if flag: break

    if len(currentval) > 0:
        replacedata.append( currentval.copy() )

    def get_line(replace_info):
        if replace_info['argname'] == '':
            return [':rtype: {0}{1}'.format('' if replace_info['argkeywords'] == [None,None] else '[type=' + replace_info['argkeywords'][1] + '] ', replace_info['argtext'])]
        return [':param {0}: {1}{2}'.format( replace_info['argname'],  \
                                             '' if replace_info['argkeywords'][0] == None else '[default=' + replace_info['argkeywords'][0] + '] ', \
                                             replace_info['argtext'])] + \
                ([] if replace_info['argkeywords'][1] == None else \
                [':type {0}: {1}'.format( replace_info['argname'], \
                                                    replace_info['argkeywords'][1])])
    
    
    #Finally, we process the data and spit out appropriate lines.
    replacedata.sort( key=lambda x:x['range'][-1], reverse=True ) # biggest final item is last, to avoid clobbering when we do slice replacement.
    for i in replacedata:
        if len(i['range']) == 1:
            lines[i['range'][0]:i['range'][0]+1] = get_line(i)[:]
        else:
            lines[i['range'][0]:i['range'][1]+1] = get_line(i)[:]

    for l in lines:
        k = l.lower()
        if k.startswith("arguments:") or k.startswith("keyword arguments:") or k.startswith("return value:"):
            lines.remove(l)
    return
# -- General configuration -----------------------------------------------------

# If your documentation needs a minimal Sphinx version, state it here.
#needs_sphinx = '1.0'

# Add any Sphinx extension module names here, as strings. They can be extensions
# coming with Sphinx (named 'sphinx.ext.*') or your custom ones.
extensions = ['sphinx.ext.autodoc','rst2pdf.pdfbuilder', 'sphinx.ext.doctest', 'sphinx.ext.viewcode']

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# The suffix of source filenames.
source_suffix = '.rst'

# The encoding of source files.
#source_encoding = 'utf-8-sig'

# The master toctree document.
master_doc = 'index'

# General information about the project.
project = u'Multistrand'
copyright = u'Caltech 2010-2017 <help@multistran.org>'

# The version info for the project you're documenting, acts as replacement for
# |version| and |release|, also used in various other places throughout the
# built documents.
#
# The short X.Y version.
version = '2.1'
# The full version, including alpha/beta/rc tags.
release = '2.1'

# The language for content autogenerated by Sphinx. Refer to documentation
# for a list of supported languages.
#language = None

# There are two options for replacing |today|: either, you set today to some
# non-false value, then it is used:
#today = ''
# Else, today_fmt is used as the format for a strftime call.
#today_fmt = '%B %d, %Y'

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
exclude_patterns = []

# The reST default role (used for this markup: `text`) to use for all documents.
#default_role = None

# If true, '()' will be appended to :func: etc. cross-reference text.
#add_function_parentheses = True

# If true, the current module name will be prepended to all description
# unit titles (such as .. function::).
#add_module_names = True

# If true, sectionauthor and moduleauthor directives will be shown in the
# output. They are ignored by default.
#show_authors = False

# The name of the Pygments (syntax highlighting) style to use.
pygments_style = 'sphinx'

# A list of ignored prefixes for module index sorting.
#modindex_common_prefix = []

# Include __init__ docstrings.
autoclass_content = 'both'

autodoc_default_flags = ['members']
# add this for other defaults: ['undoc-members', 'inherited-members', 'show-inheritance']
autodoc_member_order = 'groupwise'

# -- Options for HTML output ---------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
html_theme = 'default'

# Theme options are theme-specific and customize the look and feel of a theme
# further.  For a list of options available for each theme, see the
# documentation.
#html_theme_options = {}

# Add any paths that contain custom themes here, relative to this directory.
#html_theme_path = []

# The name for this set of Sphinx documents.  If None, it defaults to
# "<project> v<release> documentation".
#html_title = None

# A shorter title for the navigation bar.  Default is the same as html_title.
#html_short_title = None

# The name of an image file (relative to this directory) to place at the top
# of the sidebar.
#html_logo = None

# The name of an image file (within the static path) to use as favicon of the
# docs.  This file should be a Windows icon file (.ico) being 16x16 or 32x32
# pixels large.
#html_favicon = None

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']

# If not '', a 'Last updated on:' timestamp is inserted at every page bottom,
# using the given strftime format.
#html_last_updated_fmt = '%b %d, %Y'

# If true, SmartyPants will be used to convert quotes and dashes to
# typographically correct entities.
#html_use_smartypants = True

# Custom sidebar templates, maps document names to template names.
#html_sidebars = {}

# Additional templates that should be rendered to pages, maps page names to
# template names.
#html_additional_pages = {}

# If false, no module index is generated.
#html_domain_indices = True

# If false, no index is generated.
#html_use_index = True

# If true, the index is split into individual pages for each letter.
#html_split_index = False

# If true, links to the reST sources are added to the pages.
#html_show_sourcelink = True

# If true, "Created using Sphinx" is shown in the HTML footer. Default is True.
#html_show_sphinx = True

# If true, "(C) Copyright ..." is shown in the HTML footer. Default is True.
#html_show_copyright = True

# If true, an OpenSearch description file will be output, and all pages will
# contain a <link> tag referring to it.  The value of this option must be the
# base URL from which the finished HTML is served.
#html_use_opensearch = ''

# This is the file name suffix for HTML files (e.g. ".xhtml").
#html_file_suffix = None

# Output file base name for HTML help builder.
htmlhelp_basename = 'Multistranddoc'


# -- Options for LaTeX output --------------------------------------------------

# The paper size ('letter' or 'a4').
#latex_paper_size = 'letter'

# The font size ('10pt', '11pt' or '12pt').
#latex_font_size = '10pt'

# Grouping the document tree into LaTeX files. List of tuples
# (source start file, target name, title, author, documentclass [howto/manual]).
latex_documents = [
  ('index', 'Multistrand.tex', u'Multistrand Documentation',
   u'Caltech \\textless{}help@multistrand.org\\textgreater{}', 'manual'),
]

# The name of an image file (relative to this directory) to place at the top of
# the title page.
#latex_logo = None

# For "manual" documents, if this is true, then toplevel headings are parts,
# not chapters.
#latex_use_parts = False

# If true, show page references after internal links.
#latex_show_pagerefs = False

# If true, show URL addresses after external links.
#latex_show_urls = False

# Additional stuff for the LaTeX preamble.
#latex_preamble = ''

# Documents to append as an appendix to all manuals.
#latex_appendices = []

# If false, no module index is generated.
#latex_domain_indices = True


# -- Options for manual page output --------------------------------------------

# One entry per manual page. List of tuples
# (source start file, name, description, authors, manual section).
man_pages = [
    ('index', 'multistrand', u'Multistrand Documentation',
     [u'Caltech <help@multistrand.org>'], 1),
    ('objects','objects', u'Objects',
     [u'Caltech <help@multistrand.org>'], 1)
]

# -- Options for PDF output --------------------------------------------------
# Grouping the document tree into PDF files. List of tuples
# (source start file, target name, title, author, options).
#
pdf_documents = [    ('index', 'multistrand', u'Multistrand Documentation',
                      u'Caltech'),
                     ('objects','objects', u'Objects',u'Caltech')]

pdf_language = "en_US"

# Section level that forces a break page.
# For example: 1 means top-level sections start in a new page
# 0 means disabled
pdf_break_level = 0

# When a section starts in a new page, force it to be 'even', 'odd',
# or just use 'any'
pdf_breakside = 'any'

# Insert footnotes where they are defined instead of
# at the end.
#pdf_inline_footnotes = True

# verbosity level. 0 1 or 2
#pdf_verbosity = 0

# If false, no index is generated.
pdf_use_index = False

# If false, no modindex is generated.
pdf_use_modindex = False

# If false, no coverpage is generated.
#pdf_use_coverpage = True

# Name of the cover page template to use
#pdf_cover_template = 'sphinxcover.tmpl'

# Documents to append as an appendix to all manuals.
#pdf_appendices = []

# Enable experimental feature to split table cells. Use it
# if you get "DelayedTable too big" errors
#pdf_splittables = False

# Set the default DPI for images
#pdf_default_dpi = 72

# Enable rst2pdf extension modules (default is empty list)
# you need vectorpdf for better sphinx's graphviz support
#pdf_extensions = ['vectorpdf']

# Page template name for "regular" pages
#pdf_page_template = 'cutePage'
