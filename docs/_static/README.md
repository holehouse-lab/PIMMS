# Static Doc Directory

Custom static files for the docs (style sheets, JavaScript, images) go here. Sphinx copies them into the built site after the theme's own static files, so a file named `default.css` here would overwrite the builtin `default.css`.

The path to this folder is set in the Sphinx `conf.py` file in the line:
```python
html_static_path = ['_static']
```

At present the directory holds `custom.css`, which `conf.py` adds to every page with:
```python
html_css_files = ['custom.css']
```

## Examples of files to add to this directory
* Custom Cascading Style Sheets
* Custom JavaScript code
* Static logo images
