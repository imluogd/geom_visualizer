# -*- coding: utf-8 -*-
import os
import glob
import io
import html
import numpy as np
from ase.io import read
import cclib

HTML_TEMPLATE = r"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8"/>
<meta name="viewport" content="width=device-width,initial-scale=1"/>
<title>Gaussian16 Transition State Visualization</title>
<style>
  body{margin:0;padding:24px;background:#f7f7f9;color:#111;
       font-family:-apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,Helvetica,Arial,"PingFang SC","Microsoft YaHei",sans-serif;}
  .grid{display:grid;grid-template-columns:repeat(auto-fill,minmax(400px,1fr));gap:16px}
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

def atoms_to_xyz(atoms, comment=""):
    """Convert ASE atoms to XYZ format string."""
    lines = []
    lines.append(str(len(atoms)))
    lines.append(comment)
    for atom in atoms:
        symbol = atom.symbol
        x, y, z = atom.position
        lines.append(f"{symbol} {x:.6f} {y:.6f} {z:.6f}")
    return "\n".join(lines)

def create_jsmol_vibration_script(atoms, frequencies, normal_modes):
    """
    Create JSmol script for vibrational visualization.
    """
    if not frequencies or not normal_modes:
        return None
    coords = []
    for atom in atoms:
        symbol = atom.symbol
        x, y, z = atom.position
        coords.append(f"{symbol} {x:.6f} {y:.6f} {z:.6f}")
    natoms = len(atoms)
    script_lines = ["load data 'model VIBRATIONS'"]
    for i, (freq, mode) in enumerate(zip(frequencies, normal_modes)):
        script_lines.append(f"{natoms}")
        script_lines.append(f"freq = {freq:.2f} cm-1")
        for j, coord_line in enumerate(coords):
            dx, dy, dz = mode[j]
            script_lines.append(f"{coord_line} 0 {dx:.6f} {dy:.6f} {dz:.6f}")
    script_lines.append("end 'model VIBRATIONS'")
    return "\n".join(script_lines).replace("\n", "|")

def xyz_to_jsmol_data_script(xyz_text):
    """Convert multi-line XYZ to JSmol inline 'load data' script."""
    xyz_oneline = xyz_text.replace("\n", "|")
    # 不使用 'show data;'，避免生成额外信息区/空框
    return "data 'model X'|" + xyz_oneline + "|end 'model X';"

def parse_gaussian_output(log_file):
    """Parse Gaussian16 log file to get final structure and imaginary freqs."""
    try:
        parser = cclib.io.ccopen(log_file)
        data = parser.parse()
        if hasattr(data, 'atomcoords') and len(data.atomcoords) > 0:
            final_coords = data.atomcoords[-1]
            atomic_numbers = data.atomnos
            from ase import Atoms
            atoms = Atoms(numbers=atomic_numbers, positions=final_coords)
        else:
            atoms = read(log_file, format="gaussian-out", index=-1)
            data = None
        imaginary_frequencies, imaginary_modes = [], []
        if data and hasattr(data, 'vibfreqs') and hasattr(data, 'vibdisps'):
            for i, freq in enumerate(data.vibfreqs):
                if freq < 0:
                    imaginary_frequencies.append(freq)
                    if i < len(data.vibdisps):
                        imaginary_modes.append(data.vibdisps[i])
        return atoms, imaginary_frequencies, imaginary_modes
    except Exception as e:
        print(f"Error parsing {log_file}: {e}")
        try:
            atoms = read(log_file, format="gaussian-out", index=-1)
            return atoms, [], []
        except:
            return None, [], []

def _tail_contains(path, needle, n_lines=400):
    """Read last n_lines of a text file and test whether 'needle' appears."""
    try:
        with open(path, 'rb') as f:
            f.seek(0, os.SEEK_END)
            size = f.tell()
            # 估算读取字节数（每行~120字节为粗略估计）
            approx = min(size, n_lines * 200)
            f.seek(max(0, size - approx), os.SEEK_SET)
            tail = f.read().decode(errors='ignore')
        return needle in tail
    except Exception:
        return False

def gaussian_normally_terminated(path):
    """
    判定 Gaussian 日志是否“正常结束”：
    - 包含 'Normal termination of Gaussian'
    - 且不包含 'Error termination'
    """
    has_normal = _tail_contains(path, "Normal termination of Gaussian", n_lines=600)
    has_error  = _tail_contains(path, "Error termination", n_lines=600)
    return has_normal and not has_error

def make_vibration_card(idx, title, atoms, frequencies, normal_modes):
    """Create HTML card with vibrational visualization for imaginary modes."""
    safe_title = html.escape(title)
    natoms = len(atoms)
    formula = atoms.get_chemical_formula()
    n_imaginary = len(frequencies)
    warning_text = ""
    if n_imaginary == 0:
        warning_text = '<div class="warning">⚠ No imaginary frequencies (minimum structure)</div>'
    elif n_imaginary > 1:
        warning_text = f'<div class="warning">⚠ {n_imaginary} imaginary frequencies (higher-order saddle point)</div>'

    # Base structure
    xyz = atoms_to_xyz(atoms, comment=os.path.basename(title))
    base_script = xyz_to_jsmol_data_script(xyz)

    # Vibrations (if any)
    vibr_script = ""
    if n_imaginary > 0:
        vibr_script = create_jsmol_vibration_script(atoms, frequencies, normal_modes) or ""

    # Escape for JS
    esc_base = base_script.replace('"', '&quot;').replace("'", "\\'")
    esc_vibr = vibr_script.replace('"', '&quot;').replace("'", "\\'") if vibr_script else ""

    # Controls
    if n_imaginary == 1:
        freq = frequencies[0]
        vibration_controls = f"""
        <div class="controls">
          <button onclick="showStructure{idx}()">Structure</button>
          <button onclick="showVibration{idx}()">Show vibration</button>
          <span class="freq-info imaginary">{freq:.1f} cm⁻¹</span>
        </div>"""
    elif n_imaginary > 1:
        options = "\n".join(
            f'<option value="{i+1}">{f:.1f} cm⁻¹</option>' for i, f in enumerate(frequencies)
        )
        vibration_controls = f"""
        <div class="controls">
          <button onclick="showStructure{idx}()">Structure</button>
          <button onclick="showVibrations{idx}()">Show vibrations</button>
          <select id="freqSelect{idx}" onchange="changeVibration{idx}()" disabled>
            <option value="0">Select imaginary mode...</option>
            {options}
          </select>
        </div>"""
    else:
        vibration_controls = f"""
        <div class="controls">
          <button onclick="showStructure{idx}()">Structure</button>
          <span style="color:#888;">No imaginary frequencies</span>
        </div>"""

    # 只生成一次 applet：getAppletHtml("id", Info) 注入
    javascript_code = f"""
      <script>
        (function(){{
          if(!window.Applets) window.Applets = {{}};
          var baseScript = "{esc_base}";
          var vibrScript = "{esc_vibr}";
          // 默认振动幅度更明显一些
          var vibScale = 0.8;

          var Info = {{
            width: "100%",
            height: 400,
            debug: false,
            color: "0xFFFFFF",
            use: "HTML5",
            j2sPath: "https://chemapps.stolaf.edu/jmol/jsmol/j2s",
            script: baseScript,
            disableJ2SLoadMonitor: true,
            disableInitialConsole: true,
            allowJavaScript: true,
            serverURL: "https://chemapps.stolaf.edu/jmol/jsmol/php/jsmol.php",
            addSelectionOptions: false,
            console: "none"
          }};
          var html = Jmol.getAppletHtml("app{idx}", Info);
          var host = document.getElementById("app{idx}");
          host.innerHTML = html;
          window.Applets["app{idx}"] = Jmol._applets["app{idx}"];
          // 初始样式
          Jmol.script(window.Applets["app{idx}"], "background white; set antialiasDisplay on;");
          Jmol.script(window.Applets["app{idx}"], "select *; spacefill 23%; wireframe 0.15;");
          Jmol.script(window.Applets["app{idx}"], "zoom 100; center all;");

          // 切换函数
          window.showStructure{idx} = function() {{
            Jmol.script(window.Applets["app{idx}"], baseScript);
            Jmol.script(window.Applets["app{idx}"], "vibration off; select *; spacefill 23%; wireframe 0.15;");
            var s = document.getElementById("freqSelect{idx}");
            if (s) s.disabled = true;
          }};
          window.showVibration{idx} = function() {{
            if (!vibrScript) return;
            Jmol.script(window.Applets["app{idx}"], vibrScript);
            Jmol.script(window.Applets["app{idx}"], "vibration on; vibration scale " + vibScale + "; model 1;");
          }};
          window.showVibrations{idx} = function() {{
            if (!vibrScript) return;
            Jmol.script(window.Applets["app{idx}"], vibrScript);
            Jmol.script(window.Applets["app{idx}"], "vibration on; vibration scale " + vibScale + "; model 1;");
            var s = document.getElementById("freqSelect{idx}");
            if (s) s.disabled = false;
          }};
          window.changeVibration{idx} = function() {{
            var s = document.getElementById("freqSelect{idx}");
            var m = s ? s.value : 0;
            if (m > 0) {{
              Jmol.script(window.Applets["app{idx}"], "model " + m + "; vibration on; vibration scale " + vibScale + ";");
            }}
          }};
        }})();
      </script>"""

    return f"""
    <div class="card">
      <div class="meta">
        <div class="name">{safe_title}</div>
        <div class="summary">Atoms: {natoms} | Formula: {formula} | Imaginary modes: {n_imaginary}</div>
        {warning_text}
      </div>
      {vibration_controls}
      <div id="app{idx}" class="viewer"></div>
      {javascript_code}
    </div>
    """

def main():
    logs = sorted(glob.glob("*.log"))
    if not logs:
        raise SystemExit("No .log files found in current folder.")
    print("=" * 60)

    cards_html = []
    n_valid_ts = n_minimum = n_higher_order = n_errors = n_skipped = 0

    for i, log_file in enumerate(logs):
        print(f"Processing {log_file}...")

        # 新增：未正常结束的日志直接跳过并打印
        if not gaussian_normally_terminated(log_file):
            if _tail_contains(log_file, "Error termination", n_lines=600):
                print(f"[SKIP] {log_file}: error-terminated (Gaussian reported an error).")
            else:
                print(f"[SKIP] {log_file}: not normally terminated (no 'Normal termination of Gaussian').")
            n_skipped += 1
            continue

        atoms, frequencies, normal_modes = parse_gaussian_output(log_file)
        if atoms is None:
            print(f"[ERROR] Failed to parse {log_file}")
            n_errors += 1
            continue

        n_imaginary = len(frequencies)
        if n_imaginary == 1:
            n_valid_ts += 1
            print(f"  -> Valid transition state: {frequencies[0]:.2f} cm⁻¹")
        elif n_imaginary == 0:
            n_minimum += 1
            print(f"  -> WARNING: No imaginary frequencies found (appears to be minimum)")
        else:
            n_higher_order += 1
            print(f"  -> WARNING: {n_imaginary} imaginary frequencies: {[f'{f:.2f}' for f in frequencies]} cm⁻¹")

        card = make_vibration_card(i - n_skipped, log_file, atoms, frequencies, normal_modes)
        cards_html.append(card)

    if not cards_html:
        raise SystemExit("No valid structures found in log files.")

    html_doc = HTML_TEMPLATE + "\n".join(cards_html) + "\n</div>\n</body>\n</html>"

    output_file = "gaussian_ts_analysis.html"
    with open(output_file, "w", encoding="utf-8") as fh:
        fh.write(html_doc)

    print("\n" + "=" * 60)
    print("TRANSITION STATE ANALYSIS SUMMARY")
    print("=" * 60)
    print(f"Valid transition states:     {n_valid_ts}")
    print(f"Minimum structures:          {n_minimum}")
    print(f"Higher-order saddle points:  {n_higher_order}")
    print(f"Parse errors:                {n_errors}")
    print(f"Skipped (not finished):      {n_skipped}")
    print(f"Total processed:             {len(cards_html)}")
    print(f"\nOutput: {os.path.abspath(output_file)}")

if __name__ == "__main__":
    main()
