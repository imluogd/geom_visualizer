# -*- coding: utf-8 -*-
import os
import glob
import io
import html
from ase.io import read, write
#for path verification

HTML_TEMPLATE = r"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8"/>
<meta name="viewport" content="width=device-width,initial-scale=1"/>
<title>Gaussian16 Visualization</title>
<style>
  body{margin:0;padding:24px;background:#f7f7f9;color:#111;
       font-family:-apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,Helvetica,Arial,"PingFang SC","Microsoft YaHei",sans-serif;}
  h1{margin:0 0 16px 0;font-size:22px}
  .grid{display:grid;grid-template-columns:repeat(auto-fill,minmax(360px,1fr));gap:16px}
  .card{background:#fff;border-radius:14px;box-shadow:0 2px 14px rgba(0,0,0,.06);
        padding:12px 12px 16px 12px;display:flex;flex-direction:column}
  .meta{padding:2px 4px 10px 4px}
  .name{font-weight:600;font-size:15px}
  .summary{color:#555;font-size:12px;margin-top:4px}
  .viewer{width:100%;height:400px;border-radius:12px;overflow:hidden;border:1px solid #eee}
</style>
<!-- JSmol from CDN -->
<script src="https://chemapps.stolaf.edu/jmol/jsmol/JSmol.min.js"></script>
</head>
<body>
<h1>Gaussian16 visualization</h1>
<div class="grid">
"""


def atoms_to_xyz(atoms, comment=""):
    """Return XYZ string (with natoms & comment line)."""
    buf = io.StringIO()
    write(buf, atoms, format="xyz", comment=comment)
    return buf.getvalue()


def xyz_to_jsmol_data_script(xyz_text):
    """
    将多行 XYZ 转成 JSmol 可内联的 'load data' 脚本。
    """
    xyz_oneline = xyz_text.replace("\n", "|")
    return "data 'model X'|" + xyz_oneline + "|end 'model X';show data;"


def make_card(idx, title, natoms, formula, jsmol_script):
    # 每个卡片一个 <div id="appN"> 容器 + 初始化脚本
    safe_title = html.escape(title)
    escaped_jscript = jsmol_script.replace('"', '&quot;').replace("'", "\\'")
    
    return f"""
    <div class="card">
      <div class="meta">
        <div class="name">{safe_title}</div>
        <div class="summary">Atoms: {natoms} | Formula: {formula}</div>
      </div>
      <div id="app{idx}" class="viewer"></div>
      <script>
        (function(){{
          if(!window.Applets) window.Applets = {{}};
          var Info = {{
            width: "100%",
            height: "100%",
            debug: false,
            color: "0xFFFFFF",
            use: "HTML5",
            j2sPath: "https://chemapps.stolaf.edu/jmol/jsmol/j2s",
            script: "{escaped_jscript}",
            disableJ2SLoadMonitor: true,
            disableInitialConsole: true,
            allowJavaScript: true,
            serverURL: "https://chemapps.stolaf.edu/jmol/jsmol/php/jsmol.php",
            addSelectionOptions: false,
            console: "none"
          }};
          window.Applets["app{idx}"] = Jmol.getApplet("app{idx}", Info);
          document.getElementById("app{idx}").innerHTML = Jmol.getAppletHtml(window.Applets["app{idx}"]);
          Jmol.script(window.Applets["app{idx}"], "background white; set antialiasDisplay on;");
          Jmol.script(window.Applets["app{idx}"], "select *; spacefill 23%; wireframe 0.15;");
          Jmol.script(window.Applets["app{idx}"], "zoom 100; center all;");
        }})();
      </script>
    </div>
    """


def main():
    logs = sorted(glob.glob("*.log"))
    if not logs:
        raise SystemExit("No .log files found in current folder.")

    cards_html = []
    n_ok = 0
    for i, f in enumerate(logs):
        try:
            atoms = read(f, format="gaussian-out", index=-1)  # 取最后一帧（最终结构）
        except Exception as e:
            print(f"[WARN] Skip {f}: {e}")
            continue

        xyz = atoms_to_xyz(atoms, comment=os.path.basename(f))
        jsmol_script = xyz_to_jsmol_data_script(xyz)
        formula = atoms.get_chemical_formula()
        card = make_card(i, f, len(atoms), formula, jsmol_script)
        cards_html.append(card)
        n_ok += 1

    if not cards_html:
        raise SystemExit("Parsed 0 final structures from .log files.")

    # 直接闭合HTML文档，不需要额外的tail
    html_doc = HTML_TEMPLATE + "\n".join(cards_html) + "\n</div>\n</body>\n</html>"
    out = "gallery_jsmol.html"
    with open(out, "w", encoding="utf-8") as fh:
        fh.write(html_doc)
    print(
        f"Done. Wrote interactive gallery to: {os.path.abspath(out)}  (parsed {n_ok} structures)")


if __name__ == "__main__":
    main()