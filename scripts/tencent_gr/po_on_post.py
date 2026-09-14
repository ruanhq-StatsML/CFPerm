#!/usr/bin/env python3
"""转化粒 FSDS：满窗 y_post_clk_1d，X=User/Ctx/交叉（SKU 闸关）。

不是预报榜。RF-domain = P(X) 谁来了；PO/LOGO/LOCO = 哪包/哪列扛早/晚对 Y 的差。

  PYTHONPATH=. python3 scripts/tencent_gr/po_on_post.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import kendalltau

sys.path.insert(0, str(Path(__file__).resolve().parent))
from cross_feats import add_cross  # noqa: E402
from fit_post_cnv_model import BRANCHES, daily_mix, prep  # noqa: E402
from po_fs_logo import po_risk_fit, rf_domain  # noqa: E402


def mass_tower(vimp, names, families):
    out = {g: 0.0 for g in families}
    for v, n in zip(vimp, names):
        out[tower_of(n)] += float(v)
    s = sum(out.values()) + 1e-12
    return {k: out[k] / s for k in families}

ROOT = Path(__file__).resolve().parents[2]
POST = ROOT / "results/tencent_gr_fs150/tables/post.parquet"
USER = ROOT / "results/tencent_gr_fs150/user_feats_selected.parquet"
OUT = ROOT / "results/tencent_gr_fs150"
SEED = 0

USER_COLS = {
    "n_clk_before_7d",
    "n_clk_before_1d",
    "n_clk_before_1h",
    "lag_post_clk_1d_rate",
    "lag_empty_any",
    "n_prior_cnv",
}
CTX_COLS = {
    "empty_any",
    "dt_any_min",
    "dt_any_min_miss",
    "sess_clk_before",
    "sess_pos",
    "log1p_price",
}


def tower_of(name: str) -> str:
    if name.startswith("uc_"):
        return "cross"
    if name.startswith("us_"):
        return "sku"
    if name in USER_COLS:
        return "user"
    if name in CTX_COLS:
        return "ctx"
    return "other"


def mix_by_w(post: pd.DataFrame, w: np.ndarray) -> dict:
    out = {}
    for lab, m in (("early", w == 0), ("late", w == 1)):
        g = post.loc[m]
        out[lab] = {
            "n": int(len(g)),
            "pos": float(g["_y"].mean()),
            "empty_any": float(g["empty_any"].mean()),
            "sess0": float((g["sess_clk_before"].fillna(0) == 0).mean()),
            "sku_clk": float((~g["wo_prior_clk"].fillna(True).astype(bool)).mean()),
        }
    return out


def board(X, y, w, names, *, seed: int, do_loco: bool) -> dict:
    families = sorted({tower_of(n) for n in names})
    groups = {g: [i for i, n in enumerate(names) if tower_of(n) == g] for g in families}
    auc, v_rf = rf_domain(X, w, seed=seed)
    full = po_risk_fit(X, y, w, seed=seed)
    logo = {}
    for g, ix in groups.items():
        keep = [j for j in range(X.shape[1]) if j not in set(ix)]
        if not keep:
            logo[g] = {"R_minus": float("nan"), "delta": float("nan")}
            continue
        rm = po_risk_fit(X[:, keep], y, w, seed=seed + 11 + families.index(g))["risk"]
        logo[g] = {"R_minus": rm, "delta": float(full["risk"] - rm)}
        print(f"  LOGO {g:6s} Δ={logo[g]['delta']:+.6f}", flush=True)
    pos = {g: max(logo[g]["delta"] or 0.0, 0.0) for g in families}
    s = sum(pos.values()) + 1e-12
    share = {g: pos[g] / s for g in families}
    po_v = full["vimp"]
    rf_v = v_rf
    tau = float(kendalltau(po_v, rf_v).correlation)
    order = np.argsort(-po_v)
    top_po = [
        {"rank": r + 1, "name": names[j], "tower": tower_of(names[j]), "po_vimp": float(po_v[j]), "rf_vimp": float(rf_v[j])}
        for r, j in enumerate(order)
    ]
    loco = []
    if do_loco:
        for j, n in enumerate(names):
            keep = [k for k in range(X.shape[1]) if k != j]
            rm = po_risk_fit(X[:, keep], y, w, seed=seed + 80 + j)["risk"]
            d = float(full["risk"] - rm)
            loco.append({"name": n, "tower": tower_of(n), "impurity": float(po_v[j]), "loco_dR": d})
            print(f"  LOCO {n:42s} ΔR={d:+.6e}", flush=True)
        loco.sort(key=lambda r: -r["loco_dR"])
    return {
        "n": int(len(y)),
        "p": int(len(names)),
        "pos_rate": float(y.mean()),
        "w1_rate": float(w.mean()),
        "rf_domain_auc": auc,
        "po_risk": full["risk"],
        "family_n": {g: len(groups[g]) for g in families},
        "rf_mass": mass_tower(rf_v, names, families),
        "po_mass": mass_tower(po_v, names, families),
        "logo": logo,
        "logo_share": share,
        "kendall_po_vs_rf": tau,
        "ranked_po": top_po,
        "loco": loco,
        "groups": {g: [names[i] for i in groups[g]] for g in families},
    }


def write_md(path: Path, mix: dict, clocks: dict) -> None:
    lines = [
        "# 转化粒 FSDS（User / Ctx / 交叉）",
        "",
        "Y = 满窗 `y_post_clk_1d`。X 停在 `cnv_ts`。SKU 闸关，同品列不进。",
        "不是预报榜。RF-domain = P(X)；PO/LOGO/LOCO = 哪包/哪列扛早/晚对 Y 的差。",
        "",
        f"daily mix：n_cnv={mix['n_cnv']}  sku_clk={mix['sku_clk_share']:.4f}  empty={mix['empty_any']:.3f}  sess0={mix['sess_clk0']:.3f}  GATE={mix['sku_gate']}",
        "",
    ]
    for cname, rec in clocks.items():
        lines += [
            f"## 时钟 `{cname}`",
            "",
            f"n={rec['n']} pos={rec['pos_rate']:.3f} W1={rec['w1_rate']:.3f}  "
            f"RF-domain AUC **{rec['rf_domain_auc']:.3f}**  PO-risk R **{rec['po_risk']:.6f}**  "
            f"τ(PO, RF)={rec['kendall_po_vs_rf']:.3f}",
            "",
            "构成（early vs late）：",
            "",
            "| | n | pos | empty_any | sess0 | sku_clk |",
            "|---|---:|---:|---:|---:|---:|",
        ]
        for lab in ("early", "late"):
            m = rec["mix"][lab]
            lines.append(
                f"| {lab} | {m['n']} | {m['pos']:.3f} | {m['empty_any']:.3f} | {m['sess0']:.3f} | {m['sku_clk']:.4f} |"
            )
        lines += [
            "",
            "| tower | n | RF-domain | PO-VIMP | LOGO Δ | LOGO share |",
            "|---|---:|---:|---:|---:|---:|",
        ]
        fam = rec["logo_share"]
        for g in sorted(fam, key=lambda x: -fam[x]):
            lg = rec["logo"][g]
            dlt = lg["delta"]
            dlt_s = f"{dlt:+.5f}" if dlt == dlt else "nan"
            lines.append(
                f"| {g} | {rec['family_n'][g]} | {rec['rf_mass'][g]:.3f} | "
                f"{rec['po_mass'][g]:.3f} | {dlt_s} | {fam[g]:.3f} |"
            )
        lines += ["", "PO-VIMP 序：", "", "| rank | feat | tower | PO-VIMP | RF-VIMP |", "|---|---|---|---:|---:|"]
        for r in rec["ranked_po"]:
            lines.append(
                f"| {r['rank']} | `{r['name']}` | {r['tower']} | {r['po_vimp']:.4f} | {r['rf_vimp']:.4f} |"
            )
        if rec.get("loco"):
            lines += [
                "",
                "LOCO ΔR（>0 = 这列在扛 PO-risk；impurity 头名可以是负的）：",
                "",
                "| feat | tower | impurity | LOCO ΔR |",
                "|---|---|---:|---:|",
            ]
            for r in rec["loco"]:
                lines.append(
                    f"| `{r['name']}` | {r['tower']} | {r['impurity']:.4f} | {r['loco_dR']:+.2e} |"
                )
        lines.append("")
    lines += [
        "用户时钟 `_seq_t_end` 和 150 板同一把刀（后来的人）。`cnv_ts` 是单的早晚，跟满窗 follow 缠在一起，只作对照。",
        "SKU 漏斗构成两批都接近 0，不进 X。",
        "",
    ]
    path.write_text("\n".join(lines), encoding="utf-8")
    print("wrote", path)


def main() -> None:
    post = pd.read_parquet(POST)
    mix = daily_mix(post)
    print("daily mix", json.dumps(mix, indent=2), flush=True)

    ycol = "y_post_clk_1d"
    sub = post.loc[post[ycol].notna()].copy()
    sub["_y"] = sub[ycol].astype(int)
    cols = [c for b in ("heat", "empty", "any_path", "sess", "lag") for c in BRANCHES[b]]
    x = add_cross(prep(sub, cols), sku_on=False)
    names = list(x.columns)
    X = x.fillna(0.0).to_numpy(np.float64)
    y = sub["_y"].to_numpy(np.int64)
    print(f"满窗 n={len(sub)} p={len(names)} pos={float(y.mean()):.3f}", flush=True)
    print("towers", {g: sum(1 for n in names if tower_of(n) == g) for g in ("user", "ctx", "cross", "other")}, flush=True)

    u = pd.read_parquet(USER, columns=["user_id", "_seq_t_end"])
    u["user_id"] = u["user_id"].astype(np.int64)
    sub = sub.merge(u, on="user_id", how="left")
    clocks_w = {
        "seq_t_end_median": (sub["_seq_t_end"].to_numpy(np.float64) > np.nanmedian(sub["_seq_t_end"].to_numpy(np.float64))).astype(int),
        "cnv_ts_median": (sub["cnv_ts"].to_numpy(np.float64) > np.median(sub["cnv_ts"].to_numpy(np.float64))).astype(int),
    }

    clocks = {}
    for i, (cname, w) in enumerate(clocks_w.items()):
        print(f"\n== {cname} W1={w.mean():.3f} ==", flush=True)
        rec = board(X, y, w, names, seed=SEED + 17 * i, do_loco=(cname == "seq_t_end_median"))
        rec["mix"] = mix_by_w(sub, w)
        rec["clock"] = cname
        clocks[cname] = rec
        print(
            json.dumps(
                {"clock": cname, "rf_auc": rec["rf_domain_auc"], "po_risk": rec["po_risk"], "logo_share": rec["logo_share"]},
                indent=2,
            ),
            flush=True,
        )

    dump = {k: {kk: vv for kk, vv in rec.items() if kk != "groups"} for k, rec in clocks.items()}
    (OUT / "POST_FSDS.json").write_text(json.dumps({"mix": mix, "clocks": dump}, indent=2), encoding="utf-8")
    write_md(OUT / "POST_FSDS.md", mix, clocks)


if __name__ == "__main__":
    main()
