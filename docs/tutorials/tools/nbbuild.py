"""Execute a list of notebook cells in one namespace and write an .ipynb with their outputs (nbformat 4, by hand).

stdout / stderr of each cell, matplotlib figures and pyvista plots (Plotter.show() is captured as a PNG, like the
"static" Jupyter backend) become the outputs. Run with the pytec3d Python; the notebook paths are relative to the
directory of the notebook (the runner changes to it first).
"""
import base64
import getpass
import io
import json
import os
import socket
import sys
import traceback

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402


class Notebook:
    def __init__(self):
        self.cells = []

    def md(self, text):
        self.cells.append(("markdown", text.strip("\n")))

    def code(self, text):
        self.cells.append(("code", text.strip("\n")))


def _lines(text):
    lines = text.split("\n")
    return [line + "\n" for line in lines[:-1]] + ([lines[-1]] if lines[-1] else [])


def run(nb, path):
    path = os.path.abspath(path)
    os.chdir(os.path.dirname(path))
    # Neutral provenance in the stored outputs (no host, user or script path of the build machine).
    socket.gethostname = lambda: "workstation"
    getpass.getuser = lambda: "user"
    sys.argv = [os.path.basename(path)]
    # the open figures are captured after each cell: plt.show() is a no-op (with Agg it only warns)
    plt.show = lambda *args, **kwargs: None

    images = []
    try:
        import pyvista as pv

        pv.OFF_SCREEN = True

        def show(self, *args, **kwargs):
            png = io.BytesIO()
            img = self.screenshot(return_img=True)
            plt.imsave(png, img, format="png")
            images.append(png.getvalue())
            self.close()

        pv.Plotter.show = show
        pv.set_jupyter_backend = lambda *args, **kwargs: None   # needs IPython; the show() above stands for 'static'
        # pyvista sends VTK messages to Python logging; give them back to VTK (they then go to the terminal, as in a
        # notebook where pyvista is imported later).
        from vtkmodules.vtkCommonCore import vtkOutputWindow
        vtkOutputWindow.SetInstance(vtkOutputWindow())
    except ImportError:
        pass

    class Capture(io.TextIOBase):
        def __init__(self):
            self.buffer_ = io.StringIO()

        def write(self, s):
            self.buffer_.write(s)
            return len(s)

    namespace = {"__name__": "__main__"}
    nb_cells, count = [], 0
    out, err = Capture(), Capture()
    real = (sys.stdout, sys.stderr)
    sys.stdout, sys.stderr = out, err
    failed = False
    try:
        for kind, source in nb.cells:
            if kind == "markdown":
                nb_cells.append({"cell_type": "markdown", "metadata": {}, "source": _lines(source)})
                continue
            count += 1
            out.buffer_, err.buffer_ = io.StringIO(), io.StringIO()
            images.clear()
            outputs = []
            try:
                exec(compile(source, f"<cell {count}>", "exec"), namespace)
            except Exception:
                failed = True
                err.write(traceback.format_exc())
                real[1].write(f"cell {count} failed:\n{traceback.format_exc()}\n")
            for name, stream in (("stdout", out), ("stderr", err)):
                text = stream.buffer_.getvalue()
                if text:
                    outputs.append({"output_type": "stream", "name": name, "text": _lines(text)})
            pngs = []
            for num in plt.get_fignums():
                png = io.BytesIO()
                plt.figure(num).savefig(png, format="png", dpi=90, bbox_inches="tight")
                pngs.append(png.getvalue())
            plt.close("all")
            for data in pngs + images:
                outputs.append({"output_type": "display_data", "metadata": {},
                                "data": {"image/png": base64.b64encode(data).decode(), "text/plain": ["<Figure>"]}})
            nb_cells.append({"cell_type": "code", "execution_count": count, "metadata": {}, "outputs": outputs,
                             "source": _lines(source)})
    finally:
        sys.stdout, sys.stderr = real
    for i, c in enumerate(nb_cells):
        c["id"] = f"cell-{i}"
    notebook = {"cells": nb_cells, "metadata": {
        "kernelspec": {"display_name": "Python 3 (pytec3d)", "language": "python", "name": "python3"},
        "language_info": {"name": "python", "version": "3.12"}}, "nbformat": 4, "nbformat_minor": 5}
    with open(path, "w", encoding="utf-8") as handle:
        json.dump(notebook, handle, indent=1)
        handle.write("\n")
    print(("FAILED " if failed else "ok ") + path)
