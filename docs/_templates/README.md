# Templates Doc Directory

Custom Jinja templates for the HTML builder go here. Sphinx looks in this directory before the theme's own templates, so a file named `page.html` here would take the place of the builtin `page.html`. The directory currently holds no templates.

The path to this folder is set in the Sphinx `conf.py` file in the line:
```python
templates_path = ['_templates']
```

## Examples of files to add to this directory
* HTML extensions of stock pages like `page.html` or `layout.html`
