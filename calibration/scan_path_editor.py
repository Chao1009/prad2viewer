#!/usr/bin/env python3
"""
HyCal Path Editor
=================
GUI for manually building scan paths by clicking modules on the HyCal map.
Paths are saved to calibration/paths.json with named profiles.

When two clicked modules are far apart, intermediate modules whose centres
lie on the connecting line (within tolerance) are auto-inserted.

Usage
-----
    python calibration/hycal_path_editor.py
"""

from __future__ import annotations

import json
import math
import os
import tkinter as tk
from tkinter import ttk, messagebox, simpledialog
from typing import Dict, List, Tuple

from scan_utils import C, Module, load_modules, filter_scan_modules, DEFAULT_DB_PATH


# ============================================================================
#  PATH PROFILES
# ============================================================================

PATHS_FILE = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                          "paths.json")
MAX_LG_LAYERS = 6


def load_paths(filepath: str = PATHS_FILE) -> Dict[str, List[str]]:
    if os.path.exists(filepath):
        with open(filepath) as f:
            return json.load(f)
    return {}


def save_paths(profiles: Dict[str, List[str]], filepath: str = PATHS_FILE):
    # Compact arrays: each profile on ~one line
    with open(filepath, "w") as f:
        f.write("{\n")
        items = list(profiles.items())
        for idx, (name, path) in enumerate(items):
            arr = json.dumps(path, separators=(", ", ": "))
            comma = "," if idx < len(items) - 1 else ""
            f.write(f"  {json.dumps(name)}: {arr}{comma}\n")
        f.write("}\n")


# ============================================================================
#  AUTO-INSERT: find modules on the line between two points
# ============================================================================

def _point_to_segment_dist(px, py, ax, ay, bx, by) -> float:
    """Distance from point (px,py) to segment (ax,ay)-(bx,by)."""
    dx, dy = bx - ax, by - ay
    len2 = dx * dx + dy * dy
    if len2 < 1e-9:
        return math.sqrt((px - ax) ** 2 + (py - ay) ** 2)
    t = max(0, min(1, ((px - ax) * dx + (py - ay) * dy) / len2))
    proj_x = ax + t * dx
    proj_y = ay + t * dy
    return math.sqrt((px - proj_x) ** 2 + (py - proj_y) ** 2)


def find_intermediate_modules(
    m_from: Module, m_to: Module,
    candidates: List[Module], tol_frac: float = 0.4,
) -> List[Module]:
    """Return modules whose centre lies on the line from m_from to m_to.

    Tolerance = tol_frac * min(module_sx, module_sy) for each candidate.
    Results are ordered by distance from m_from.
    """
    ax, ay = m_from.x, m_from.y
    bx, by = m_to.x, m_to.y
    seg_len = math.sqrt((bx - ax) ** 2 + (by - ay) ** 2)
    if seg_len < 1e-6:
        return []

    hits: List[Tuple[float, Module]] = []
    for m in candidates:
        if m.name == m_from.name or m.name == m_to.name:
            continue
        tol = tol_frac * min(m.sx, m.sy)
        d = _point_to_segment_dist(m.x, m.y, ax, ay, bx, by)
        if d <= tol:
            dist_from_a = math.sqrt((m.x - ax) ** 2 + (m.y - ay) ** 2)
            hits.append((dist_from_a, m))
    hits.sort(key=lambda h: h[0])
    return [m for _, m in hits]


# ============================================================================
#  GUI
# ============================================================================

class PathEditorGUI:

    CANVAS_SIZE = 680
    CANVAS_PAD  = 8
    MOD_SHRINK  = 0.90

    def __init__(self, root: tk.Tk, all_modules: List[Module],
                 paths_file: str = PATHS_FILE):
        self.root = root
        self.all_modules = all_modules
        self._mod_by_name: Dict[str, Module] = {m.name: m for m in all_modules}
        self._paths_file = paths_file
        self._lg_layers = 0

        # Precompute PbWO4 bbox for layer filter
        pwo4 = [m for m in all_modules if m.mod_type == "PbWO4"]
        self._pwo4_min_x = min(m.x for m in pwo4)
        self._pwo4_max_x = max(m.x for m in pwo4)
        self._pwo4_min_y = min(m.y for m in pwo4)
        self._pwo4_max_y = max(m.y for m in pwo4)
        glass = [m for m in all_modules if m.mod_type == "PbGlass"]
        self._lg_sx = glass[0].sx if glass else 38.15
        self._lg_sy = glass[0].sy if glass else 38.15

        self._selectable: set = set()   # names of selectable modules
        self._update_selectable(0)

        # Current path & profile
        self._path: List[str] = []      # ordered module names
        self._profile_name: str = ""
        self._profiles = load_paths(self._paths_file)

        # Canvas state
        self._cell_ids: Dict[str, int] = {}
        self._scale = 1.0
        self._ox = self._oy = 0.0
        self._x_min = self._y_max = 0.0

        self._build_ui()

    # -- selectable filter ---------------------------------------------------

    def _update_selectable(self, lg_layers: int):
        self._lg_layers = lg_layers
        self._selectable = {m.name for m in
                            filter_scan_modules(self.all_modules, lg_layers,
                                                self._lg_sx, self._lg_sy)}

    # -- UI ------------------------------------------------------------------

    def _build_ui(self):
        self.root.title("HyCal Path Editor")
        self.root.configure(bg=C.BG)
        self.root.resizable(True, True)

        style = ttk.Style()
        style.theme_use("clam")
        style.configure(".", background=C.BG, foreground=C.TEXT,
                         fieldbackground=C.PANEL, bordercolor=C.BORDER)
        style.configure("TLabel", background=C.BG, foreground=C.TEXT)
        style.configure("TLabelframe", background=C.BG, foreground=C.ACCENT)
        style.configure("TLabelframe.Label", background=C.BG,
                         foreground=C.ACCENT, font=("Consolas", 9, "bold"))
        style.configure("TButton", background=C.PANEL, foreground=C.TEXT,
                         padding=4)
        style.map("TButton",
                  background=[("active", C.BORDER)],
                  foreground=[("disabled", "#484f58")])
        style.configure("Accent.TButton", background="#1f6feb",
                         foreground="white")
        style.configure("Danger.TButton", background="#da3633",
                         foreground="white")
        style.configure("Green.TButton", background="#238636",
                         foreground="white")
        style.configure("TCombobox", fieldbackground=C.PANEL,
                         background=C.BORDER, foreground=C.TEXT,
                         selectbackground=C.BORDER,
                         selectforeground=C.TEXT)
        style.map("TCombobox",
                  fieldbackground=[("readonly", C.PANEL),
                                   ("disabled", C.BG)],
                  foreground=[("disabled", "#484f58")])

        # -- top bar ---------------------------------------------------------
        top = tk.Frame(self.root, bg="#0d1520", height=32)
        top.pack(fill="x")
        tk.Label(top, text="  HYCAL PATH EDITOR  ", bg="#0d1520", fg=C.GREEN,
                 font=("Consolas", 13, "bold")).pack(side="left", padx=8)

        # -- main area -------------------------------------------------------
        main = tk.Frame(self.root, bg=C.BG)
        main.pack(fill="both", expand=True, padx=6, pady=4)

        left = tk.Frame(main, bg=C.BG)
        left.grid(row=0, column=0, sticky="nsew", padx=(0, 4))
        right = tk.Frame(main, bg=C.BG)
        right.grid(row=0, column=1, sticky="nsew")
        main.columnconfigure(0, weight=0)
        main.columnconfigure(1, weight=1)
        main.rowconfigure(0, weight=1)

        self._build_canvas(left)
        self._build_controls(right)

    # -- canvas --------------------------------------------------------------

    def _build_canvas(self, parent):
        self._canvas_frame = ttk.LabelFrame(parent, text="")
        self._update_canvas_label()
        self._canvas_frame.pack(fill="both", expand=True)

        sz = self.CANVAS_SIZE
        self._canvas = tk.Canvas(self._canvas_frame, width=sz, height=sz,
                                  bg="#0a0e14", highlightthickness=0)
        self._canvas.pack(padx=4, pady=4)
        self._canvas.bind("<Button-1>", self._on_canvas_click)

        self._compute_canvas_mapping()
        self._draw_modules()

        # legend
        leg = tk.Frame(self._canvas_frame, bg=C.BG)
        leg.pack(fill="x", padx=4, pady=(0, 4))
        for label, colour in [("Selectable", C.MOD_TODO),
                               ("In Path", C.MOD_INPATH),
                               ("Excluded", C.MOD_EXCLUDED),
                               ("PbGlass", C.MOD_GLASS)]:
            tk.Canvas(leg, width=10, height=10, bg=colour,
                      highlightthickness=0).pack(side="left", padx=(6, 1))
            tk.Label(leg, text=label, bg=C.BG, fg=C.DIM,
                     font=("Consolas", 8)).pack(side="left")

    def _update_canvas_label(self):
        n = len(self._path)
        name = self._profile_name or "(unsaved)"
        self._canvas_frame.configure(
            text=f" Module Map | Profile: {name} | Path: {n} modules ")

    def _compute_canvas_mapping(self):
        drawn = [m for m in self.all_modules if m.mod_type != "LMS"]
        if not drawn:
            return
        x_min = min(m.x - m.sx / 2 for m in drawn)
        x_max = max(m.x + m.sx / 2 for m in drawn)
        y_min = min(m.y - m.sy / 2 for m in drawn)
        y_max = max(m.y + m.sy / 2 for m in drawn)
        usable = self.CANVAS_SIZE - 2 * self.CANVAS_PAD
        self._scale = min(usable / (x_max - x_min),
                          usable / (y_max - y_min))
        draw_w = (x_max - x_min) * self._scale
        draw_h = (y_max - y_min) * self._scale
        self._ox = self.CANVAS_PAD + (usable - draw_w) / 2
        self._oy = self.CANVAS_PAD + (usable - draw_h) / 2
        self._x_min = x_min
        self._y_max = y_max

    def _mod_to_canvas(self, m: Module) -> Tuple[float, float, float, float]:
        cx = self._ox + (m.x - self._x_min) * self._scale
        cy = self._oy + (self._y_max - m.y) * self._scale
        hw = m.sx * self._scale * self.MOD_SHRINK / 2
        hh = m.sy * self._scale * self.MOD_SHRINK / 2
        return (cx - hw, cy - hh, cx + hw, cy + hh)

    def _mod_center(self, m: Module) -> Tuple[float, float]:
        cx = self._ox + (m.x - self._x_min) * self._scale
        cy = self._oy + (self._y_max - m.y) * self._scale
        return cx, cy

    def _draw_modules(self):
        self._canvas.delete("all")
        self._cell_ids.clear()
        path_set = set(self._path)

        for m in self.all_modules:
            if m.mod_type == "LMS":
                continue
            x0, y0, x1, y1 = self._mod_to_canvas(m)
            if m.name in path_set:
                color = C.MOD_INPATH
            elif m.name in self._selectable:
                color = C.MOD_TODO
            elif m.mod_type == "PbGlass":
                color = C.MOD_GLASS
            else:
                color = C.MOD_EXCLUDED
            rid = self._canvas.create_rectangle(
                x0, y0, x1, y1, fill=color, outline="", width=0,
                tags=(f"mod_{m.name}",))
            self._cell_ids[m.name] = rid

        self._draw_path_line()

    def _draw_path_line(self):
        self._canvas.delete("path_line")
        if len(self._path) < 2:
            return
        coords = []
        for name in self._path:
            m = self._mod_by_name.get(name)
            if m:
                coords.extend(self._mod_center(m))
        if len(coords) >= 4:
            self._canvas.create_line(
                *coords, fill=C.PATH_LINE, width=1.5,
                smooth=False, tags=("path_line",))
            self._canvas.tag_raise("path_line")

    def _refresh_canvas(self):
        """Lightweight refresh: recolor modules and redraw path line."""
        path_set = set(self._path)
        for m in self.all_modules:
            if m.mod_type == "LMS":
                continue
            rid = self._cell_ids.get(m.name)
            if rid is None:
                continue
            if m.name in path_set:
                color = C.MOD_INPATH
            elif m.name in self._selectable:
                color = C.MOD_TODO
            elif m.mod_type == "PbGlass":
                color = C.MOD_GLASS
            else:
                color = C.MOD_EXCLUDED
            self._canvas.itemconfigure(rid, fill=color)
        self._draw_path_line()
        self._update_canvas_label()
        self._update_path_list()

    # -- click handler -------------------------------------------------------

    def _on_canvas_click(self, event):
        items = self._canvas.find_closest(event.x, event.y)
        if not items:
            return
        tags = self._canvas.gettags(items[0])
        name = None
        for tag in tags:
            if tag.startswith("mod_"):
                name = tag[4:]
                break
        if name is None or name not in self._selectable:
            return

        # If already last in path, ignore (prevent double-click duplicates)
        if self._path and self._path[-1] == name:
            return

        # Auto-insert intermediate modules if line crosses others
        if self._path:
            last_mod = self._mod_by_name[self._path[-1]]
            new_mod = self._mod_by_name[name]
            candidates = [self._mod_by_name[n] for n in self._selectable
                          if n not in set(self._path) and n != name]
            intermediates = find_intermediate_modules(last_mod, new_mod,
                                                      candidates)
            for im in intermediates:
                if im.name not in set(self._path):
                    self._path.append(im.name)

        self._path.append(name)
        self._refresh_canvas()

    # -- controls panel ------------------------------------------------------

    def _build_controls(self, parent):
        # === Profile ===
        pf = ttk.LabelFrame(parent, text=" Profile ")
        pf.pack(fill="x", pady=(0, 4))

        r_pf = tk.Frame(pf, bg=C.BG)
        r_pf.pack(fill="x", padx=6, pady=4)
        tk.Label(r_pf, text="Profile:", bg=C.BG, fg=C.TEXT,
                 font=("Consolas", 9)).pack(side="left")
        self._profile_var = tk.StringVar()
        self._profile_combo = ttk.Combobox(
            r_pf, textvariable=self._profile_var,
            values=sorted(self._profiles.keys()), width=16,
            font=("Consolas", 9), state="readonly")
        self._profile_combo.pack(side="left", padx=(4, 4))
        self._profile_combo.bind("<<ComboboxSelected>>", self._on_load_profile)
        self.root.option_add("*TCombobox*Listbox.background", C.PANEL)
        self.root.option_add("*TCombobox*Listbox.foreground", C.TEXT)
        self.root.option_add("*TCombobox*Listbox.selectBackground", C.ACCENT)
        self.root.option_add("*TCombobox*Listbox.selectForeground", "white")

        ttk.Button(r_pf, text="New", style="Green.TButton",
                   command=self._cmd_new_profile).pack(side="left", padx=2)
        ttk.Button(r_pf, text="Save", style="Accent.TButton",
                   command=self._cmd_save).pack(side="left", padx=2)
        ttk.Button(r_pf, text="Delete", style="Danger.TButton",
                   command=self._cmd_delete).pack(side="left", padx=2)

        # === Settings ===
        sf = ttk.LabelFrame(parent, text=" Settings ")
        sf.pack(fill="x", pady=(0, 4))

        r_lg = tk.Frame(sf, bg=C.BG)
        r_lg.pack(fill="x", padx=6, pady=4)
        tk.Label(r_lg, text="LG layers (0-6):", bg=C.BG, fg=C.TEXT,
                 font=("Consolas", 9)).pack(side="left")
        self._lg_var = tk.IntVar(value=0)
        tk.Spinbox(r_lg, from_=0, to=MAX_LG_LAYERS,
                   textvariable=self._lg_var, width=4,
                   bg=C.PANEL, fg=C.TEXT, font=("Consolas", 9),
                   buttonbackground=C.BORDER, insertbackground=C.TEXT,
                   command=self._on_lg_changed).pack(side="right")

        # === Path editing ===
        ef = ttk.LabelFrame(parent, text=" Path ")
        ef.pack(fill="both", expand=True, pady=(0, 4))

        bf = tk.Frame(ef, bg=C.BG)
        bf.pack(fill="x", padx=6, pady=4)
        ttk.Button(bf, text="Undo Last",
                   command=self._cmd_undo).pack(side="left", padx=2)
        ttk.Button(bf, text="Clear All", style="Danger.TButton",
                   command=self._cmd_clear).pack(side="left", padx=2)

        self._path_list = tk.Text(ef, height=20, bg=C.PANEL, fg=C.TEXT,
                                   font=("Consolas", 9), wrap="word",
                                   state="disabled", borderwidth=0,
                                   highlightthickness=0)
        self._path_list.pack(fill="both", expand=True, padx=6, pady=(0, 4))
        self._update_path_list()

    def _update_path_list(self):
        self._path_list.configure(state="normal")
        self._path_list.delete("1.0", "end")
        if self._path:
            self._path_list.insert("end",
                                    ", ".join(self._path) + f"\n\n({len(self._path)} modules)")
        else:
            self._path_list.insert("end", "(empty)")
        self._path_list.configure(state="disabled")

    # -- commands ------------------------------------------------------------

    def _on_lg_changed(self):
        self._update_selectable(self._lg_var.get())
        self._draw_modules()
        self._refresh_canvas()

    def _on_load_profile(self, _event=None):
        name = self._profile_var.get()
        if name in self._profiles:
            self._profile_name = name
            self._path = list(self._profiles[name])
            self._refresh_canvas()

    def _cmd_new_profile(self):
        name = simpledialog.askstring("New Profile", "Profile name:",
                                       parent=self.root)
        if not name or not name.strip():
            return
        name = name.strip()
        self._profile_name = name
        self._path = []
        self._profiles[name] = []
        self._profile_combo["values"] = sorted(self._profiles.keys())
        self._profile_var.set(name)
        self._refresh_canvas()

    def _cmd_save(self):
        if not self._profile_name:
            self._cmd_new_profile()
            if not self._profile_name:
                return
        self._profiles[self._profile_name] = list(self._path)
        save_paths(self._profiles, self._paths_file)
        self._profile_combo["values"] = sorted(self._profiles.keys())
        messagebox.showinfo("Saved",
                            f"Profile '{self._profile_name}' saved "
                            f"({len(self._path)} modules)")

    def _cmd_delete(self):
        name = self._profile_name or self._profile_var.get()
        if not name or name not in self._profiles:
            return
        if not messagebox.askyesno("Delete", f"Delete profile '{name}'?"):
            return
        del self._profiles[name]
        save_paths(self._profiles, self._paths_file)
        self._profile_combo["values"] = sorted(self._profiles.keys())
        self._profile_var.set("")
        self._profile_name = ""
        self._path = []
        self._refresh_canvas()

    def _cmd_undo(self):
        if self._path:
            self._path.pop()
            self._refresh_canvas()

    def _cmd_clear(self):
        if self._path and messagebox.askyesno("Clear", "Clear entire path?"):
            self._path.clear()
            self._refresh_canvas()


# ============================================================================
#  MAIN
# ============================================================================

def main():
    import argparse
    parser = argparse.ArgumentParser(
        description="HyCal Scan Path Editor")
    parser.add_argument("--database", default=DEFAULT_DB_PATH,
                        help="Path to hycal_map.json")
    parser.add_argument("--paths", default=PATHS_FILE,
                        help="Path to paths.json")
    args = parser.parse_args()

    all_modules = load_modules(args.database)
    print(f"Loaded {len(all_modules)} modules")

    root = tk.Tk()
    PathEditorGUI(root, all_modules, paths_file=args.paths)
    root.mainloop()


if __name__ == "__main__":
    main()
