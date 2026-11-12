"""Small pyvista helpers of the tutorial notebooks: draw the regions of a model and the nodes, edges or vectors a setup
step selects or writes, so that each step can be checked on the mesh before the model is loaded.

    from setup_views import view
    p = view(grid, regions={1: "air", 2: "conductor", 3: ("my label", "#ff8800")})   # cells coloured by TYPE
    p.nodes(ids, "red", "_A = 0")                                  # a set of points
    p.edges(start, end, "blue", "tree edges")                     # edges (point ids of their ends)
    p.show()
"""
from __future__ import annotations

import numpy as np
import pyvista as pv

COLORS = {"air": "#e8eef4", "iron": "#9aa3ad", "coil": "#e39a3b", "conductor": "#c0703a", "layer": "#8fc1e3"}
MARKERS = ["#c0392b", "#27ae60", "#1f5fa8", "#8e44ad", "#d35400", "#16a085"]


class View:
    """A pyvista plotter on a pytec3d grid, with a legend of what is drawn."""

    def __init__(self, grid, regions: dict | None = None, opacity: float | None = None, edges: bool = True,
                 size=(1250, 680), clip: dict | None = None):
        self.grid = pv.wrap(grid).copy()
        self.xyz = np.asarray(self.grid.points)
        self.flat = np.ptp(self.xyz[:, 2]) < 1e-12           # a 2D mesh
        self.p = pv.Plotter(window_size=size)
        self.p.set_background("white")
        self.entries: list[list[str]] = []
        self.n = 0
        g = self.grid.clip(**clip) if clip else self.grid
        for t, name in (regions or {}).items():
            part = g.extract_cells(np.flatnonzero(g.cell_data["TYPE"] == t))
            if part.n_cells == 0:
                continue
            name, color = name if isinstance(name, tuple) else (name, COLORS.get(name, name))
            alpha = opacity if opacity is not None else (1.0 if self.flat else (0.12 if name == "air" else 0.45))
            self.p.add_mesh(part, color=color, opacity=alpha, lighting=not self.flat)
            self.entries.append([f"TYPE {t}: {name}", "#7a8794" if name == "air" else color])
        if edges:
            solid = [t for t, nm in (regions or {}).items() if nm != "air"]
            shown = g if self.flat else g.extract_cells(np.flatnonzero(np.isin(g.cell_data["TYPE"], solid)))
            if shown.n_cells:
                self.p.add_mesh(shown.linear_copy().extract_all_edges(), color="#9aa5b1", line_width=0.5)

    def _color(self, color):
        if color is None:
            color = MARKERS[self.n % len(MARKERS)]
            self.n += 1
        return color

    def nodes(self, ids, color=None, label=None, size=8):
        """The points ``ids`` (point ids, or a boolean mask)."""
        ids = np.flatnonzero(ids) if np.asarray(ids).dtype == bool else np.asarray(ids, dtype=int)
        color = self._color(color)
        if len(ids):
            self.p.add_mesh(pv.PolyData(self.xyz[ids]), color=color, point_size=size, render_points_as_spheres=True)
        if label:
            self.entries.append([f"{label} ({len(ids)})", color])
        return self

    def edges(self, start, end, color=None, label=None, width=3, values=None, cmap="viridis", arrows=False):
        """Edges by the point ids of their ends; ``values`` colours them, ``arrows`` shows their orientation."""
        start, end = np.asarray(start, dtype=int), np.asarray(end, dtype=int)
        n = len(start)
        if n == 0:
            return self
        color = self._color(color) if values is None else None
        if arrows:
            a, b = self.xyz[start], self.xyz[end]
            cloud = pv.PolyData(0.5 * (a + b) - 0.3 * (b - a))
            cloud["v"] = b - a
            length = float(np.linalg.norm(b - a, axis=1).mean())
            self.p.add_mesh(cloud.glyph(orient="v", scale=False, factor=0.6 * length), color=color)
        else:
            lines = np.column_stack((np.full(n, 2), np.arange(n), np.arange(n) + n)).ravel()
            poly = pv.PolyData(np.vstack((self.xyz[start], self.xyz[end])), lines=lines)
            if values is not None:
                poly.cell_data[label or "value"] = values
                self.p.add_mesh(poly, scalars=label or "value", cmap=cmap, line_width=width,
                                scalar_bar_args={"color": "black"})
                return self
            self.p.add_mesh(poly, color=color, line_width=width)
        if label:
            self.entries.append([f"{label} ({n})", color])
        return self

    def vectors(self, positions, vectors, color=None, label=None, factor=None):
        """Arrows of ``vectors`` at ``positions`` (e.g. Js at the cell centroids)."""
        cloud = pv.PolyData(np.asarray(positions, dtype=float))
        cloud["v"] = np.asarray(vectors, dtype=float)
        size = factor or 0.06 * float(np.linalg.norm(np.ptp(self.xyz, axis=0)))
        color = self._color(color)
        self.p.add_mesh(cloud.glyph(orient="v", scale=False, factor=size), color=color)
        if label:
            self.entries.append([label, color])
        return self

    def show(self, title=None, view=None, zoom=1.0, bounds=None):
        """Render: ``view`` "xy", "xz", "yz" or a camera position; ``bounds`` to frame a part of the model."""
        if self.entries:
            width = min(0.6, 0.0065 * max(len(label) for label, _ in self.entries) + 0.03)
            self.p.add_legend(self.entries, bcolor="white", face="circle", loc="lower left",
                              size=(width, 0.032 * len(self.entries) + 0.01))
        if title:
            self.p.add_text(title, font_size=10, color="black")
        if self.flat or view == "xy":
            self.p.view_xy()
        elif view == "xz":
            self.p.view_xz()
        elif view == "yz":
            self.p.view_yz()
        elif view is not None:
            self.p.camera_position = view
        else:
            self.p.view_isometric()
        if bounds is not None:
            self.p.reset_camera(bounds=bounds)
        self.p.camera.zoom(zoom)
        self.p.show()


def view(grid, **kwargs) -> View:
    return View(grid, **kwargs)
