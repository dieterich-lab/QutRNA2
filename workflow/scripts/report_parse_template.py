import click
import jinja2
import yaml


@click.command()
@click.option("-o", "--output", type=str, required=True)
@click.option("-p", "--params", type=click.Path(exists=True))
@click.argument("template", type=click.Path(exists=True))
def parse(template, output, params):
  """Parse jinja template"""

  d = {}
  if params:
    with open(params, "r") as f:
      d = yaml.safe_load(f)

  with open(template, "r") as in_file:
    s = in_file.read()
    jinja2.filters.FILTERS.update(
      {
        "read_yaml": yaml.safe_load,
        "to_yaml": yaml.dump,
      })
    t = jinja2.Template(s)
    r = t.render(**d)
    with open(output, "w") as out_file:
      out_file.write(r)


if __name__ == "__main__":
  parse()