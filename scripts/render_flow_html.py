#!/usr/bin/env python3
"""One HTML: PO / MMD / MSE flow, OnlineRFPerm check, ADDIS/SAFFRON, boards."""
from __future__ import annotations

import html
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "Python" / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from layer_freeze_cv import attach_layer_dicts  # noqa: E402

OUT = ROOT / "results" / "layer_freeze_online_cv"
DATASETS = (
    "dgp_concept",
    "dgp_covariate",
    "bankmarketing",
    "eeg",
    "electricity",
    "covertype",
)
ACTION_ZH = {
    "keep_training": "接着 train",
    "watch": "再观察",
    "x_shift": "X shift",
    "tricky": "tricky",
    "freeze": "冻层",
}
ACTION_COLOR = {
    "keep_training": "#2a7d4f",
    "watch": "#d48b16",
    "x_shift": "#6b4ea0",
    "tricky": "#555",
    "freeze": "#b33",
}


def _fmt(v, n=3):
    if v is None:
        return ""
    try:
        x = float(v)
    except (TypeError, ValueError):
        return html.escape(str(v))
    if not np_ok(x):
        return ""
    return f"{x:.{n}g}"


def np_ok(x) -> bool:
    try:
        return x == x and abs(x) != float("inf")
    except Exception:
        return False


def load_result(name: str) -> dict | None:
    path = OUT / name / "summary.json"
    if not path.exists():
        return None
    result = json.loads(path.read_text(encoding="utf-8"))
    result.setdefault("dataset", name)
    return attach_layer_dicts(result)


def row_cells(r: dict) -> str:
    act = r.get("action") or "keep_training"
    color = ACTION_COLOR.get(act, "#333")
    hop = "●" if r.get("rfperm_hop") else ""
    tested = "MSE" if r.get("fdr_tested") else "—"
    addis = "reject" if r.get("fdr_reject_addis") or r.get("fdr_reject") else ""
    saff = "reject" if r.get("fdr_reject_saffron") else ""
    if r.get("fdr_discarded") and not r.get("fdr_tested"):
        addis = "discard"
    elif r.get("fdr_discarded"):
        addis = addis or "discard"
    return (
        f"<tr>"
        f"<td>{r.get('t')}</td>"
        f"<td>{_fmt(r.get('po_stream'))}</td><td>{_fmt(r.get('po_ma'))}</td>"
        f"<td>{_fmt(r.get('mse_stream'))}</td><td>{_fmt(r.get('mse_ma'))}</td>"
        f"<td>{_fmt(r.get('mmd_vs_ref') if r.get('mmd_vs_ref') is not None else r.get('mmd_stream'))}</td>"
        f"<td>{_fmt(r.get('rfperm_T'))}</td><td>{hop}</td>"
        f"<td>{_fmt(r.get('rfperm_p'))}</td>"
        f"<td>{tested}</td>"
        f"<td>{_fmt(r.get('fdr_alpha'))}</td>"
        f"<td>{addis}</td><td>{saff}</td>"
        f"<td style='color:{color};font-weight:600'>{ACTION_ZH.get(act, act)}</td>"
        f"</tr>"
    )


def dataset_section(name: str, result: dict) -> str:
    title = html.escape(str(result.get("title") or name))
    board = result.get("board_action") or "keep_training"
    rec = ACTION_ZH.get(board, board)
    color = ACTION_COLOR.get(board, "#1f4e79")
    onset = result.get("onset_hat")
    onset_fdr = result.get("onset_fdr")
    n_test = (result.get("fdr_addis") or {}).get("n_fdr_tested", 0)
    n_rej = (result.get("fdr_addis") or {}).get("n_fdr_reject", 0)
    n_saff = (result.get("fdr_saffron") or {}).get("n_fdr_reject", 0)
    trend = OUT / name / "po_mse_trend.png"
    img = f'<p class="shot"><img src="{name}/po_mse_trend.png" alt="{name} trend"></p>' if trend.exists() else ""
    board_link = f'<a href="{name}/board.html">单独看板</a>' if (OUT / name / "board.html").exists() else ""
    rows = result.get("rows") or []
    body = "".join(row_cells(r) for r in rows)
    has_rf = any(r.get("rfperm_T") is not None for r in rows)
    note = ""
    if not has_rf:
        note = "<p class='muted'>这张表还没有 OnlineRFPerm 序列（后来才接到闭环上）。读 PO × MSE × MMD 即可。</p>"
    return f"""
<section id="{name}" class="dataset">
  <h2>{title}</h2>
  <p class="pill" style="background:{color}">{html.escape(str(rec))}</p>
  <p class="meta">n_ref={result.get('n_ref')} · n_new={result.get('n_new')} · batches={result.get('n_batches')}
  · onset_hat={onset} · ADDIS onset={onset_fdr} · MSE-tested={n_test} · ADDIS reject={n_rej} · SAFFRON reject={n_saff}
  · {board_link}</p>
  {note}
  {img}
  <div class="scroll">
  <table>
    <thead><tr>
      <th>t</th><th>PO</th><th>PO MA</th><th>MSE</th><th>MSE MA</th>
      <th>MMD vs ref</th><th>RFPerm T</th><th>hop</th><th>p</th>
      <th>FDR 输入</th><th>α_t</th><th>ADDIS</th><th>SAFFRON</th><th>看板</th>
    </tr></thead>
    <tbody>{body}</tbody>
  </table>
  </div>
</section>
"""


def overview_row(name: str, result: dict) -> str:
    board = result.get("board_action") or ""
    add = result.get("fdr_addis") or {}
    return (
        f"<tr><td><a href='#{name}'>{html.escape(name)}</a></td>"
        f"<td>{html.escape(str(result.get('n_new', '')))}</td>"
        f"<td>{result.get('onset_true')}</td>"
        f"<td>{result.get('onset_hat')}</td>"
        f"<td>{add.get('onset_fdr')}</td>"
        f"<td>{add.get('n_fdr_tested', 0)}</td>"
        f"<td>{add.get('n_fdr_reject', 0)}</td>"
        f"<td>{(result.get('fdr_saffron') or {}).get('n_fdr_reject', 0)}</td>"
        f"<td>{html.escape(ACTION_ZH.get(board, board))}</td></tr>"
    )


def render(out: Path = OUT) -> Path:
    sections = []
    overview = []
    loaded = {}
    for name in DATASETS:
        result = load_result(name)
        if result is None:
            continue
        loaded[name] = result
        sections.append(dataset_section(name, result))
        overview.append(overview_row(name, result))
    contrast = ""
    if (out / "dgp_contrast.png").exists():
        contrast = '<p class="shot"><img src="dgp_contrast.png" alt="DGP contrast"></p>'
    page = f"""<!DOCTYPE html>
<html lang="zh">
<head>
<meta charset="utf-8"/>
<meta name="viewport" content="width=device-width, initial-scale=1"/>
<title>PO-risk / MMD / MSE → OnlineRFPerm → ADDIS</title>
<style>
:root {{ --ink:#122; --paper:#f7f5f0; --card:#fff; --line:#d9d3c7; --navy:#1f4e79; --muted:#5c5a55; }}
* {{ box-sizing:border-box; }}
body {{ margin:0; font-family:"IBM Plex Sans","Noto Sans SC",sans-serif; color:var(--ink); background:var(--paper); }}
header {{ background:var(--navy); color:#fff; padding:28px 32px 22px; }}
header h1 {{ margin:0 0 8px; font-size:1.55rem; font-weight:600; }}
header p {{ margin:0; max-width:880px; line-height:1.5; opacity:.92; }}
nav {{ position:sticky; top:0; background:#eef2f6; border-bottom:1px solid var(--line); padding:8px 32px; z-index:2; }}
nav a {{ color:var(--navy); margin-right:14px; text-decoration:none; font-size:.92rem; }}
main {{ padding:24px 32px 64px; max-width:1180px; }}
h2 {{ font-size:1.25rem; margin:28px 0 10px; }}
h3 {{ font-size:1.05rem; margin:18px 0 8px; }}
.flow {{ display:grid; grid-template-columns:repeat(5,1fr); gap:10px; margin:16px 0 8px; }}
.step {{ background:var(--card); border:1px solid var(--line); border-radius:10px; padding:12px 12px 14px; min-height:132px; }}
.step b {{ display:block; color:var(--navy); margin-bottom:6px; }}
.step.alert {{ border-color:#c45c26; background:#fff7f0; }}
.step code {{ font-size:.8rem; }}
.muted {{ color:var(--muted); font-size:.92rem; line-height:1.45; }}
.pill {{ display:inline-block; color:#fff; padding:4px 10px; border-radius:999px; font-size:.85rem; }}
table {{ border-collapse:collapse; background:var(--card); width:100%; font-variant-numeric:tabular-nums; }}
th, td {{ border:1px solid var(--line); padding:5px 8px; font-size:.86rem; }}
th {{ background:#eef2f6; text-align:left; }}
.scroll {{ overflow-x:auto; }}
.shot img {{ max-width:100%; background:#fff; border:1px solid var(--line); }}
.dataset {{ margin-top:36px; padding-top:8px; border-top:1px solid var(--line); }}
.meta {{ color:var(--muted); font-size:.9rem; }}
.grid2 {{ display:grid; grid-template-columns:1fr 1fr; gap:16px; }}
@media (max-width: 980px) {{ .flow, .grid2 {{ grid-template-columns:1fr; }} }}
code {{ background:#eee; padding:1px 4px; border-radius:3px; }}
</style>
</head>
<body>
<header>
  <h1>看板流程：PO-risk · MMD · MSE → OnlineRFPerm → online FDR</h1>
  <p>PO-risk 和 MMD 都是对 reference batch 的一句 readout。serving MSE 崩了，才用冻住的
  <code>RandomForestRegressor().predict(X_new)</code> 做 OnlineRFPerm 检验。
  这些 p 值走 <b>ADDIS</b>（主）和 <b>SAFFRON</b>（对照）做 online FDR。不做 online-bootstrap。</p>
</header>
<nav>
  <a href="#flow">流程</a>
  <a href="#what">WHAT</a>
  <a href="#fdr">FDR</a>
  <a href="#overview">对照</a>
  <a href="#dgp_concept">concept DGP</a>
  <a href="#dgp_covariate">covariate DGP</a>
  <a href="#bankmarketing">bank</a>
  <a href="#eeg">EEG</a>
  <a href="#electricity">electricity</a>
  <a href="#covertype">covertype</a>
</nav>
<main>
<section id="flow">
  <h2>1. 每一段 batch 怎么读</h2>
  <div class="flow">
    <div class="step"><b>1 · 新 batch T=1</b>n_ref 钉死。新来的是 T=1。何时 update 仍是业务逻辑。</div>
    <div class="step"><b>2 · PO-risk</b>独立 RF：μ(Y|X)、e(T|X)。<code>φ=(Y−μ)(T−e)</code>。问 P(Y|X) 有没有 hop。RF 不该先突然崩。</div>
    <div class="step"><b>3 · MMD²(X_new, X_ref)</b>只跟 reference 比。不是 pairwise 历史，不是上一层 representation。</div>
    <div class="step"><b>4 · serving MSE</b>现模型在新 batch 上还好不好用。评估在 update 之前。更可能先崩的是这一条。</div>
    <div class="step alert"><b>5 · MSE 崩了 → RFPerm</b><code>pred = rf.predict(X_new)</code><br>T = MSE − E_ref<br>p 进 ADDIS / SAFFRON</div>
  </div>
  <p class="muted">安静 batch 的 p 记成 1，ADDIS 直接 discard（p &gt; τ=0.5）。只在 MSE 崩的那些点花 alpha。</p>
</section>

<section id="what">
  <h2>2. WHAT：PO × MSE × MMD 对照</h2>
  <div class="scroll">
  <table>
    <thead><tr><th>PO-risk</th><th>serving MSE</th><th>MMD vs ref</th><th>读法</th><th>动作</th></tr></thead>
    <tbody>
      <tr><td>安静</td><td>安静</td><td>—</td><td>没有 hop，模型还在拟合</td><td>接着 train</td></tr>
      <tr><td>崩了</td><td>没崩</td><td>—</td><td>机制可能动了，误差还在线内</td><td>再观察</td></tr>
      <tr><td>崩了</td><td>崩了</td><td>—</td><td>hop 可见而且现模型也崩了</td><td>冻那一层的 training</td></tr>
      <tr><td>没崩</td><td>崩了</td><td>崩了</td><td>不是 P(Y|X)，是 P(X)</td><td>X shift；OnlineRFPerm + FDR 确认 WHEN</td></tr>
      <tr><td>没崩</td><td>崩了</td><td>没崩</td><td>不是 concept 也不是 X</td><td>tricky；RFPerm 仍可检验 MSE 跳变</td></tr>
    </tbody>
  </table>
  </div>
</section>

<section id="fdr">
  <h2>3. WHEN 的检验：OnlineRFPerm + online FDR</h2>
  <div class="grid2">
    <div>
      <h3>OnlineRFPerm</h3>
      <p class="muted">参考窗上 fit 一次浅 RF，之后不再更新。</p>
<pre>rf = RandomForestRegressor().fit(X_ref, Y_ref)
pred = rf.predict(np.asarray(X_new))
T = mean((Y_new - pred)**2) - E_ref
p = rank of T against the T pool</pre>
      <p class="muted">FDR 用的是相对 <b>D_ref null pool</b> 的 rank p（不是短流上 1/(t+1) 那种粗 p）。last-two hop 仍标跳变。安静 batch 的 p 记成 1，ADDIS discard。</p>
    </div>
    <div>
      <h3>ADDIS（主）/ SAFFRON（对照）</h3>
      <p class="muted">α = 0.05。ADDIS：p &gt; τ=0.5 discard，p ≤ λ=0.25 才是 candidate，因子 (τ−λ)/τ。SAFFRON：λ=0.5，不 discard。两者都是 infinite-horizon γ_t = 1/(t(t+1))，不看未来有多长。</p>
      <p class="muted">短流上 rank p 很粗（0.2 / 1/3 / 0.5），hop 可能亮、FDR 仍不 reject —— 这是控制，不是漏检 bug。</p>
    </div>
  </div>
</section>

<section id="overview">
  <h2>4. 这几张看板</h2>
  {contrast}
  <div class="scroll">
  <table>
    <thead><tr>
      <th>dataset</th><th>n_new</th><th>labeled</th><th>hop onset</th>
      <th>ADDIS onset</th><th>MSE-tested</th><th>ADDIS rej</th><th>SAFFRON rej</th><th>看板</th>
    </tr></thead>
    <tbody>
      {''.join(overview)}
    </tbody>
  </table>
  </div>
  <p class="muted">concept DGP：MMD 安静、PO/MSE 动 → 再观察。covariate DGP：MMD 走高、PO 安静。
  bank / EEG：MSE 先崩、PO 安静、MMD 过线 → X shift；OnlineRFPerm hop 标 WHEN，FDR 再确认。</p>
</section>

{''.join(sections)}
</main>
</body>
</html>
"""
    dest = out / "flow.html"
    dest.write_text(page, encoding="utf-8")
    return dest


def main() -> int:
    path = render(OUT)
    print("wrote", path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
