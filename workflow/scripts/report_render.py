"""report_render.py - render the HTML run report from a Jinja template

input:   a Jinja template, and a params YAML of template variables
output:  one self-contained HTML file at --output
usage:   report_render.py --params params.yaml --output report.html TEMPLATE
env:     click and jinja2, plus pyyaml for the params file
notes:   three filters are registered on the environment. a plot is inlined
         as base64 by embed. the pdf beside a png is named by pdf. a value is
         rendered as yaml text by to_yaml
"""
import base64
import mimetypes
import os

import click
import jinja2
import yaml


def embed(path):
    """<img> tag with the file inlined, or a placeholder when it is absent.

    A missing plot appears in the page as text. Every plot that reached the
    template is rendered.
    """
    if not path or not os.path.exists(path):
        return f'<p class="missing">missing: {path}</p>'
    mime = mimetypes.guess_type(path)[0] or "application/octet-stream"
    with open(path, "rb") as fh:
        data = base64.b64encode(fh.read()).decode("ascii")
    return (f'<img src="data:{mime};base64,{data}" '
            f'alt="{os.path.basename(path)}">')


def pdf(path):
    """Absolute path of the pdf written next to a rendered png.

    The path is absolute because the report is read from anywhere. The pdf
    itself stays in the run directory.
    """
    if not path:
        return ""
    return os.path.abspath(os.path.splitext(path)[0] + ".pdf")


@click.command()
@click.option("-o", "--output", type=str, required=True)
@click.option("-p", "--params", type=click.Path(exists=True), required=True)
@click.argument("template", type=click.Path(exists=True))
def render(template, output, params):
    """Render TEMPLATE with PARAMS into a self-contained HTML report."""
    with open(params, encoding="utf-8") as fh:
        d = yaml.safe_load(fh) or {}

    env = jinja2.Environment(
        loader=jinja2.FileSystemLoader(os.path.dirname(os.path.abspath(template))),
        autoescape=False,
    )
    env.filters["embed"] = embed
    env.filters["pdf"] = pdf
    env.filters["to_yaml"] = yaml.dump

    html = env.get_template(os.path.basename(template)).render(**d)
    with open(output, "w", encoding="utf-8") as fh:
        fh.write(html)


if __name__ == "__main__":
    render()
