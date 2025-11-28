# -*- coding: utf-8 -*-
"""
Generate a local HTML viewer that follows the AMK visualization style for
all XYZ structures contained in ./min_xyz/*.xyz.
"""
import glob
import html
import json
import os
import re
from pathlib import Path

from ase.io import read

try:
    from amk_result_viewer import (  # type: ignore
        atoms_to_xyz as amk_atoms_to_xyz,
        xyz_to_jsmol_data_script as amk_xyz_to_jsmol_data_script,
        HTML_TEMPLATE as AMK_HTML_TEMPLATE,
    )
except Exception:  # pragma: no cover - fallback when optional module is absent
    amk_atoms_to_xyz = None
    amk_xyz_to_jsmol_data_script = None
    AMK_HTML_TEMPLATE = None


def _fallback_atoms_to_xyz(atoms, comment=""):
    lines = [str(len(atoms)), comment]
    for atom in atoms:
        x, y, z = atom.position
        lines.append(f"{atom.symbol} {x:.6f} {y:.6f} {z:.6f}")
    return "\n".join(lines)


def _fallback_xyz_to_jsmol_data_script(xyz_text):
    xyz_oneline = xyz_text.replace("\n", "|")
    return "data 'model X'|" + xyz_oneline + "|end 'model X';"


atoms_to_xyz = amk_atoms_to_xyz or _fallback_atoms_to_xyz
xyz_to_jsmol_data_script = amk_xyz_to_jsmol_data_script or _fallback_xyz_to_jsmol_data_script

DEFAULT_HTML_TEMPLATE = r"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8"/>
<meta name="viewport" content="width=device-width,initial-scale=1"/>
<title>XYZ Collection Visualization</title>
<style>
  body{margin:0;padding:24px;background:#f7f7f9;color:#111;
       font-family:-apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,Helvetica,Arial,"PingFang SC","Microsoft YaHei",sans-serif;}
  .grid{display:grid;grid-template-columns:repeat(3, minmax(300px, 1fr));gap:16px}
  .card{background:#fff;border-radius:14px;box-shadow:0 2px 14px rgba(0,0,0,.06);
        padding:12px 12px 16px 12px;display:flex;flex-direction:column}
  .meta{padding:2px 4px 10px 4px}
  .name{font-weight:600;font-size:15px}
  .summary{color:#555;font-size:12px;margin-top:4px}
  .controls{padding:8px 4px;display:flex;gap:8px;align-items:center;flex-wrap:wrap}
  /* 关键：不固定外层高度，仅由 JSmol 控制，避免出现“第二个空框” */
  .viewer{width:100%;border-radius:12px;overflow:hidden;border:1px solid #eee}
  select{padding:4px 8px;border:1px solid #ddd;border-radius:6px;font-size:12px}
  button{padding:4px 12px;border:1px solid #ddd;border-radius:6px;background:#f8f9fa;cursor:pointer;font-size:12px}
  button:hover{background:#e9ecef}
  .freq-info{color:#666;font-size:11px;margin-left:8px}
  .imaginary{color:#dc3545;font-weight:600}
  .warning{color:#ffc107;font-weight:600}
</style>
<script src="https://chemapps.stolaf.edu/jmol/jsmol/JSmol.min.js"></script>
</head>
<body>
<div class="grid">
"""


def _derive_html_template():
    base = AMK_HTML_TEMPLATE or DEFAULT_HTML_TEMPLATE
    if "Gaussian16 Transition State Visualization" in base:
        base = base.replace(
            "Gaussian16 Transition State Visualization",
            "XYZ Collection Visualization",
        )
    return base


HTML_TEMPLATE = _derive_html_template()


def parse_xyz_file(path):
    """Parse a single XYZ file into an ASE Atoms object."""
    try:
        return read(path, format="xyz")
    except Exception as exc:
        print(f"[ERROR] Failed to parse {path}: {exc}")
        return None


def make_structure_card(idx, file_path, atoms):
    safe_title = html.escape(os.path.basename(file_path))
    rel_path = html.escape(os.path.relpath(file_path))
    natoms = len(atoms)
    formula = atoms.get_chemical_formula()
    xyz_text = atoms_to_xyz(atoms, comment=rel_path)
    base_script = xyz_to_jsmol_data_script(xyz_text)

    javascript_code = f"""
      <script>
        (function(){{
          if(!window.Applets) window.Applets = {{}};
          var Info = {{
            width: "100%",
            height: 420,
            debug: false,
            color: "0xFFFFFF",
            use: "HTML5",
            j2sPath: "https://chemapps.stolaf.edu/jmol/jsmol/j2s",
            script: {json.dumps(base_script)},
            disableJ2SLoadMonitor: true,
            disableInitialConsole: true,
            allowJavaScript: true,
            serverURL: "https://chemapps.stolaf.edu/jmol/jsmol/php/jsmol.php",
            addSelectionOptions: false,
            console: "none"
          }};
          var host = document.getElementById("app{idx}");
          host.innerHTML = Jmol.getAppletHtml("app{idx}", Info);
          var applet = window.Applets["app{idx}"] = Jmol._applets["app{idx}"];
          Jmol.script(applet, "background white; set antialiasDisplay on; center all; zoom 100;");

          var currentStyle = "ballstick";
          var representationCommands = {{
            ballstick: "select *; spacefill 23%; wireframe 0.15;",
            stick: "select *; spacefill 0; wireframe 0.2;",
            spacefill: "select *; spacefill 100%; wireframe 0;"
          }};

          function applyStyle(style) {{
            var cmd = representationCommands[style] || representationCommands.ballstick;
            currentStyle = style in representationCommands ? style : "ballstick";
            Jmol.script(applet, cmd);
          }}

          applyStyle(currentStyle);

          window.showBallStick{idx} = function() {{ applyStyle("ballstick"); }};
          window.showStick{idx} = function() {{ applyStyle("stick"); }};
          window.showSpacefill{idx} = function() {{ applyStyle("spacefill"); }};
          window.resetView{idx} = function() {{
            Jmol.script(applet, "reset; center all; zoom 100; spin off;");
            applyStyle(currentStyle);
          }};
        }})();
      </script>
    """

    return f"""
    <div class="card">
      <div class="meta">
        <div class="name">{safe_title}</div>
        <div class="summary">Atoms: {natoms} | Formula: {formula}</div>
        <div class="freq-info">Path: {rel_path}</div>
      </div>
      <div class="controls">
        <button onclick="showBallStick{idx}()">Ball &amp; Stick</button>
        <button onclick="showStick{idx}()">Stick</button>
        <button onclick="showSpacefill{idx}()">Spacefill</button>
        <button onclick="resetView{idx}()">Reset view</button>
      </div>
      <div id="app{idx}" class="viewer"></div>
      {javascript_code}
    </div>
    """


def build_page(cards):
    return HTML_TEMPLATE + "\n".join(cards) + "\n</div>\n</body>\n</html>"


def _natural_key(path):
    """Return a list suitable for natural sorting of filenames."""
    basename = os.path.basename(path)
    return [
        int(chunk) if chunk.isdigit() else chunk.lower()
        for chunk in re.split(r"(\d+)", basename)
    ]


def main():
    pattern = os.path.join("min_xyz", "*.xyz")
    xyz_files = sorted(glob.glob(pattern), key=_natural_key)
    if not xyz_files:
        raise SystemExit("No XYZ files found under ./min_xyz/")

    print("=" * 60)
    print(f"Discovered {len(xyz_files)} XYZ files in ./min_xyz/")
    print("=" * 60)

    cards = []
    total_atoms = 0
    skipped = 0

    for idx, file_path in enumerate(xyz_files):
        rel = os.path.relpath(file_path)
        print(f"[{idx+1:03d}/{len(xyz_files):03d}] {rel}")
        atoms = parse_xyz_file(file_path)
        if atoms is None:
            skipped += 1
            continue
        total_atoms += len(atoms)
        cards.append(make_structure_card(len(cards), file_path, atoms))

    if not cards:
        raise SystemExit("Unable to parse any XYZ files.")

    html_doc = build_page(cards)
    output_file = Path("xyz_visualization.html")
    output_file.write_text(html_doc, encoding="utf-8")

    print("\n" + "=" * 60)
    print("XYZ VISUALIZATION SUMMARY")
    print("=" * 60)
    print(f"Total files discovered : {len(xyz_files)}")
    print(f"Successfully rendered  : {len(cards)}")
    print(f"Skipped (parse errors) : {skipped}")
    print(f"Total atoms displayed  : {total_atoms}")
    print(f"Output file            : {output_file.resolve()}")
    print("Open the HTML file in your browser to view all structures.")


if __name__ == "__main__":
    main()
