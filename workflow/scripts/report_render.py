"""Render the HTML report from a Jinja template.

Plots are inlined as base64 so the result is one file that can be attached to
an email or dropped in a shared folder.

usage: python report_render.py --params params.yaml --output report.html template.jinja
"""
import base64
import mimetypes
import os

import click
import jinja2
import yaml


def embed(path):
    """<img> tag with the file inlined, or a placeholder when it is absent.

    A missing plot is reported in the page rather than raised, so one failed
    plot does not cost the whole report.
    """
    if not path or not os.path.exists(path):
        return f'<p class="missing">missing: {path}</p>'
    mime = mimetypes.guess_type(path)[0] or "application/octet-stream"
    with open(path, "rb") as fh:
        data = base64.b64encode(fh.read()).decode("ascii")
    return (f'<img src="data:{mime};base64,{data}" '
            f'alt="{os.path.basename(path)}">')


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
    env.filters["to_yaml"] = yaml.dump

    html = env.get_template(os.path.basename(template)).render(**d)
    with open(output, "w", encoding="utf-8") as fh:
        fh.write(html)


if __name__ == "__main__":
    render()
