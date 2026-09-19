#!/usr/bin/env python3
"""Boss-facing heuristic MVP: PO-risk / MSE / MMD boards in one HTML."""
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
    "electricity",
    "covertype",
    "bankmarketing",
    "eeg",
)
ACTION_ZH = {
    "keep_training": "接着 train",
    "watch": "再观察",
    "x_shift": "X 在动",
    "tricky": "tricky",
    "freeze": "冻这一层",
}
ACTION_COLOR = {
    "keep_training": "#2a7d4f",
    "watch": "#d48b16",
    "x_shift": "#6b4ea0",
    "tricky": "#555",
    "freeze": "#b33",
}
BLURB = {
    "dgp_concept": "机制在慢慢翻（β），X 没动。PO / MSE 抬，MMD 安静 → 再观察，先不冻层。",
    "dgp_covariate": "同一套 P(Y|X)，X 的均值在走。MMD 过线，PO 安静，模型还能用 → 接着 train。",
    "electricity": "电价流。RF 的 PO-risk 没崩，MSE 的移动平均也没过 2× → 每一层接着 backprop。",
    "covertype": "地理顺序。MSE 先崩、PO 安静、MMD vs 参考窗过线 → 不是 concept，是 X 在动。",
    "bankmarketing": "营销活动顺序。t=4 模型误差跳，PO 仍安静，MMD 过线 → X shift。",
    "eeg": "脑电时间序。t=2 误差跳，PO 安静，MMD 过线 → 同样是 X shift。",
}


def _fmt(v, n=3):
    if v is None:
        return ""
    try:
        x = float(v)
    except (TypeError, ValueError):
        return html.escape(str(v))
    if x != x or abs(x) == float("inf"):
        return ""
    return f"{x:.{n}g}"


def load_result(name: str) -> dict | None:
    path = OUT / name / "summary.json"
    if not path.exists():
        return None
    result = json.loads(path.read_text(encoding="utf-8"))
    result.setdefault("dataset", name)
    return attach_layer_dicts(result)


def mmd_of(r: dict):
    v = r.get("mmd_vs_ref")
    if v is None:
        v = r.get("mmd_stream")
    return v


def row_cells(r: dict) -> str:
    act = r.get("action") or "keep_training"
    color = ACTION_COLOR.get(act, "#333")
    hop = "跳" if r.get("rfperm_hop") else ""
    return (
        f"<tr>"
        f"<td>{r.get('t')}</td>"
        f"<td>{_fmt(r.get('po_stream'))}</td>"
        f"<td>{_fmt(r.get('mse_stream'))}</td>"
        f"<td>{_fmt(mmd_of(r))}</td>"
        f"<td>{_fmt(r.get('rfperm_T'))}</td>"
        f"<td>{hop}</td>"
        f"<td style='color:{color};font-weight:600'>{ACTION_ZH.get(act, act)}</td>"
        f"</tr>"
    )


def dataset_section(name: str, result: dict) -> str:
    title = html.escape(str(result.get("title") or name))
    board = result.get("board_action") or "keep_training"
    rec = ACTION_ZH.get(board, board)
    color = ACTION_COLOR.get(board, "#1f4e79")
    blurb = html.escape(BLURB.get(name, ""))
    trend = OUT / name / "po_mse_trend.png"
    img = f'<p class="shot"><img src="{name}/po_mse_trend.png" alt="{name}"></p>' if trend.exists() else ""
    board_link = ""
    if (OUT / name / "board.html").exists():
        board_link = f' · <a href="{name}/board.html">明细看板</a>'
    rows = result.get("rows") or []
    body = "".join(row_cells(r) for r in rows)
    onset = result.get("onset_hat")
    onset_s = "" if onset is None else f" · 跳变点 t={onset}"
    return f"""
<section id="{name}" class="dataset">
  <div class="headrow">
    <h2>{title}</h2>
    <span class="pill" style="background:{color}">{html.escape(str(rec))}</span>
  </div>
  <p class="blurb">{blurb}</p>
  <p class="meta">n_ref={result.get('n_ref')} · 每段 n={result.get('n_new')} · {result.get('n_batches')} 段{onset_s}{board_link}</p>
  {img}
  <div class="scroll">
  <table>
    <thead><tr>
      <th>段</th><th>PO-risk</th><th>模型 MSE</th><th>MMD vs 参考窗</th>
      <th>冻住 RF 的 T</th><th>跳变</th><th>动作</th>
    </tr></thead>
    <tbody>{body}</tbody>
  </table>
  </div>
</section>
"""


def overview_card(name: str, result: dict) -> str:
    board = result.get("board_action") or "keep_training"
    color = ACTION_COLOR.get(board, "#1f4e79")
    rec = ACTION_ZH.get(board, board)
    title = html.escape(str(result.get("title") or name))
    blurb = html.escape(BLURB.get(name, ""))
    return f"""
<a class="card" href="#{name}">
  <span class="pill" style="background:{color}">{html.escape(str(rec))}</span>
  <strong>{title}</strong>
  <span>{blurb}</span>
</a>
"""


def render(out: Path = OUT) -> Path:
    sections = []
    cards = []
    for name in DATASETS:
        result = load_result(name)
        if result is None:
            continue
        sections.append(dataset_section(name, result))
        cards.append(overview_card(name, result))
    contrast = ""
    if (out / "dgp_contrast.png").exists():
        contrast = """
<section id="contrast">
  <h2>对照实验：机制在变 vs X 在走</h2>
  <p class="muted">左：concept（β 慢慢翻）— MMD 安静，PO / MSE 动。右：covariate（均值走）— MMD 过线，PO 安静。灰虚线是标好的 onset，橙线是冻住 RF 标到的跳变。</p>
  <p class="shot"><img src="dgp_contrast.png" alt="concept vs covariate"></p>
</section>
"""
    page = f"""<!DOCTYPE html>
<html lang="zh">
<head>
<meta charset="utf-8"/>
<meta name="viewport" content="width=device-width, initial-scale=1"/>
<title>在线层冻结 · 启发式看板 MVP</title>
<style>
:root {{
  --ink:#15202b; --paper:#f4f1ea; --card:#fff; --line:#ddd6c8;
  --navy:#1f4e79; --muted:#5f5b53; --green:#2a7d4f; --gold:#d48b16; --purple:#6b4ea0;
}}
* {{ box-sizing:border-box; }}
body {{ margin:0; font-family:"IBM Plex Sans","Noto Sans SC","Source Han Sans SC",sans-serif; color:var(--ink); background:var(--paper); }}
header {{ background:linear-gradient(160deg,#163a5c 0%,#1f4e79 60%,#2c6a99 100%); color:#fff; padding:36px 32px 28px; }}
header .kicker {{ letter-spacing:.12em; font-size:.72rem; opacity:.8; text-transform:uppercase; margin:0 0 10px; }}
header h1 {{ margin:0 0 10px; font-size:1.7rem; font-weight:650; line-height:1.25; }}
header p {{ margin:0; max-width:760px; line-height:1.55; opacity:.94; font-size:1.02rem; }}
nav {{ position:sticky; top:0; background:#eef2f6; border-bottom:1px solid var(--line); padding:10px 32px; z-index:3; }}
nav a {{ color:var(--navy); margin-right:16px; text-decoration:none; font-size:.9rem; }}
nav a:hover {{ text-decoration:underline; }}
main {{ padding:8px 32px 72px; max-width:1120px; margin:0 auto; }}
h2 {{ font-size:1.22rem; margin:32px 0 10px; }}
.muted {{ color:var(--muted); font-size:.95rem; line-height:1.5; }}
.flow {{ display:grid; grid-template-columns:1fr 28px 1fr 28px 1fr 28px 1fr; gap:0; align-items:stretch; margin:18px 0 8px; }}
.step {{ background:var(--card); border:1px solid var(--line); border-radius:12px; padding:14px 14px 16px; }}
.step b {{ display:block; color:var(--navy); margin-bottom:6px; font-size:.95rem; }}
.step p {{ margin:0; font-size:.88rem; line-height:1.45; color:#333; }}
.arrow {{ display:flex; align-items:center; justify-content:center; color:#9aa; font-size:1.4rem; }}
.rules {{ display:grid; grid-template-columns:repeat(5,1fr); gap:10px; margin:16px 0; }}
.rule {{ background:var(--card); border-radius:12px; border-top:5px solid #ccc; padding:12px 12px 14px; min-height:128px; }}
.rule b {{ display:block; margin-bottom:6px; }}
.rule span {{ font-size:.84rem; color:#444; line-height:1.4; }}
.cards {{ display:grid; grid-template-columns:repeat(3,1fr); gap:12px; margin:16px 0 8px; }}
a.card {{ display:block; background:var(--card); border:1px solid var(--line); border-radius:12px; padding:14px; text-decoration:none; color:inherit; min-height:148px; }}
a.card:hover {{ border-color:var(--navy); }}
a.card strong {{ display:block; margin:8px 0 6px; font-size:.98rem; }}
a.card span:last-child {{ font-size:.86rem; color:var(--muted); line-height:1.4; }}
.pill {{ display:inline-block; color:#fff; padding:3px 10px; border-radius:999px; font-size:.8rem; }}
table {{ border-collapse:collapse; background:var(--card); width:100%; font-variant-numeric:tabular-nums; }}
th, td {{ border:1px solid var(--line); padding:6px 9px; font-size:.86rem; }}
th {{ background:#eef2f6; text-align:left; }}
.scroll {{ overflow-x:auto; }}
.shot img {{ max-width:100%; background:#fff; border:1px solid var(--line); border-radius:8px; }}
.dataset {{ margin-top:40px; padding-top:4px; border-top:1px solid var(--line); }}
.headrow {{ display:flex; align-items:center; gap:12px; flex-wrap:wrap; }}
.headrow h2 {{ margin:18px 0 0; }}
.blurb {{ margin:8px 0; line-height:1.5; }}
.meta {{ color:var(--muted); font-size:.88rem; }}
.note {{ background:#fff; border-left:4px solid var(--navy); padding:12px 16px; margin:18px 0; }}
code {{ background:#eee; padding:1px 5px; border-radius:3px; font-size:.86em; }}
footer {{ color:var(--muted); font-size:.82rem; margin-top:48px; }}
@media (max-width: 980px) {{
  .flow, .rules, .cards {{ grid-template-columns:1fr; }}
  .arrow {{ display:none; }}
  header, nav, main {{ padding-left:16px; padding-right:16px; }}
}}
</style>
</head>
<body>
<header>
  <p class="kicker">Layer-freeze online CV · MVP</p>
  <h1>新数据来了，这张看板告诉你：<br>机制变了没有、模型还能不能用、该不该冻层。</h1>
  <p>启发式，不是自动停训器。三句话：<b>PO-risk</b> 看 P(Y|X) 有没有 hop，
  <b>模型 MSE</b> 看现模型还付不付得起房租，
  <b>MMD</b> 看是不是 X 自己在走。何时真正去 update，仍归业务。</p>
</header>
<nav>
  <a href="#flow">怎么读</a>
  <a href="#rules">五种动作</a>
  <a href="#boards">六张看板</a>
  <a href="#contrast">对照实验</a>
  <a href="#dgp_concept">concept</a>
  <a href="#dgp_covariate">covariate</a>
  <a href="#electricity">电价</a>
  <a href="#covertype">covertype</a>
  <a href="#bankmarketing">银行</a>
  <a href="#eeg">脑电</a>
</nav>
<main>
<section id="flow">
  <h2>每一段新 batch 怎么读</h2>
  <div class="flow">
    <div class="step"><b>1 · 钉死参考窗</b><p>前面 1 万行是 T=0。新来的这一段是 T=1。参考窗不再改。</p></div>
    <div class="arrow">→</div>
    <div class="step"><b>2 · PO-risk</b><p>独立随机森林估 μ(Y|X)、e(T|X)。问的是机制，不是现模型好不好。RF 不该突然崩。</p></div>
    <div class="arrow">→</div>
    <div class="step"><b>3 · 模型 MSE + MMD</b><p>MSE：现网络在新段上的误差。<br>MMD：只跟参考窗比 X，<code>MMD²(X_new, X_ref)</code>。</p></div>
    <div class="arrow">→</div>
    <div class="step"><b>4 · 动作</b><p>两都崩才冻层。PO 崩、模型没崩 → 再看。MSE 崩、PO 没崩 → 先看是不是 X 在动。</p></div>
  </div>
  <p class="note">口径钉死：MMD 只对 <b>参考窗</b>，不对历史所有 batch 的 pairwise 均值，也不对上一层 representation。
  冻住的 <code>RandomForestRegressor().predict(X_new)</code> 只用来标「何时跳了一下」。不做 online-bootstrap。</p>
</section>

<section id="rules">
  <h2>五种启发式动作</h2>
  <div class="rules">
    <div class="rule" style="border-top-color:var(--green)"><b>接着 train</b><span>PO 安静、MSE 安静。没有 hop，模型还在拟合。一层都不冻。</span></div>
    <div class="rule" style="border-top-color:var(--gold)"><b>再观察</b><span>PO 崩了，MSE 没崩。机制可能动了，误差还在线内。先不冻。</span></div>
    <div class="rule" style="border-top-color:#b33"><b>冻这一层</b><span>PO 和 MSE 都崩。现策略不太行。开 freeze-depth，只训上面几层。</span></div>
    <div class="rule" style="border-top-color:var(--purple)"><b>X 在动</b><span>MSE 崩了，PO 没崩，MMD 过线。不是 concept drift，是 P(X) 在变。不按机制去冻层。</span></div>
    <div class="rule" style="border-top-color:#555"><b>tricky</b><span>MSE 崩了，PO 和 MMD 都安静。不是 P(Y|X) 也不是 X。先不冻。</span></div>
  </div>
  <p class="muted">崩 = 因果移动平均 ≥ 2× 各自在参考窗上的 baseline。重点：RF 的 PO-risk 不该先突然塌；更可能先崩的是模型 MSE。</p>
</section>

<section id="boards">
  <h2>六张看板，一句话结论</h2>
  <div class="cards">
    {''.join(cards)}
  </div>
</section>

{contrast}

{''.join(sections)}

<footer>
  这是启发式 MVP，用来一起看数，不是上线决策器。何时 update 仍是业务逻辑。
  明细图和表都在本页；需要原始 JSON 看各数据集目录下的 summary.json。
</footer>
</main>
</body>
</html>
"""
    dest = out / "mvp.html"
    dest.write_text(page, encoding="utf-8")
    (out / "flow.html").write_text(page, encoding="utf-8")
    return dest


def main() -> int:
    path = render(OUT)
    index = OUT / "index.html"
    index.write_text(
        """<!DOCTYPE html><meta charset="utf-8">
<meta http-equiv="refresh" content="0; url=mvp.html">
<title>启发式看板 MVP</title>
<p><a href="mvp.html">打开给老板看的看板 MVP</a></p>
""",
        encoding="utf-8",
    )
    print("wrote", path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
