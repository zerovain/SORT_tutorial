"""Format percent-cell sources into notebooks and sync documentation copies."""
from pathlib import Path
import nbformat

root = Path(__file__).resolve().parent
for source in sorted((root / "notebooks").glob("[0-9]*.py")):
    notebook_path = source.with_suffix(".ipynb")
    previous = nbformat.read(notebook_path, as_version=4) if notebook_path.exists() else None
    cells, kind, lines = [], None, []

    def append_cell():
        if kind is None:
            return
        text = "\n".join(lines).strip()
        if kind == "markdown":
            text = "\n".join(line[2:] if line.startswith("# ") else line[1:] if line.startswith("#") else line for line in text.splitlines())
            cells.append(nbformat.v4.new_markdown_cell(text))
        else:
            cells.append(nbformat.v4.new_code_cell(text))

    for line in source.read_text().splitlines():
        if line.startswith("# %%"):
            append_cell()
            kind, lines = ("markdown" if "[markdown]" in line else "code"), []
        else:
            lines.append(line)
    append_cell()
    if previous is not None:
        for index, cell in enumerate(cells):
            if index < len(previous.cells):
                old = previous.cells[index]
                if cell.cell_type == old.cell_type == "code" and cell.source == old.source:
                    cell.outputs = old.outputs
                    cell.execution_count = old.execution_count
    book = nbformat.v4.new_notebook(cells=cells, metadata={"kernelspec": {
        "display_name": "Python 3 (ipykernel)", "language": "python", "name": "python3"},
        "language_info": {"name": "python"}})
    nbformat.validate(book)
    nbformat.write(book, notebook_path)
for notebook in (root / "notebooks").glob("*.ipynb"):
    (root / "docs/tutorials" / notebook.name).write_bytes(notebook.read_bytes())
print("Formatted source notebooks; documentation copies synchronized.")
