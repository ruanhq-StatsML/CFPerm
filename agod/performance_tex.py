"""Paste-ready LaTeX performance tables for gated online RFPerm."""
from __future__ import annotations

from pathlib import Path


def _shown(v, fmt):
    return float(fmt % v)


def _bold_min(vals, v, fmt, tol=1e-12):
    shown = [_shown(x, fmt) for x in vals]
    cell = fmt % v
    if abs(_shown(v, fmt) - min(shown)) < tol:
        return r"$\mathbf{%s}$" % cell
    return r"$%s$" % cell


def _bold_max(vals, v, fmt, tol=1e-12):
    shown = [_shown(x, fmt) for x in vals]
    cell = fmt % v
    if abs(_shown(v, fmt) - max(shown)) < tol:
        return r"$\mathbf{%s}$" % cell
    return r"$%s$" % cell


def _rmse_fmt(rmse):
    if rmse >= 100:
        return "%.1f"
    if rmse >= 10:
        return "%.2f"
    return "%.3f"


def synth_tables(summary, gate, learners):
    lines = [
        r"% Online RFPerm + PO-risk. Next-batch MSE. Uniform default.",
        r"% rfperm = last-two, w=sqrt(po_risk0) on T=1 iff consecutive OOS jump.",
        r"% resid = drop old batch (hard lever). oracle = knows the concept cut.",
        "",
    ]
    for learner in learners:
        lines += [
            r"\begin{table}[ht]\centering",
            r"\caption{Next-batch MSE (\texttt{%s}). Last-two uniform is the" % learner,
            r"default. \texttt{rfperm}: online RFPerm + PO-risk, fire iff "
            r"$e_{\mathrm{now}}/e_{\mathrm{prev}}\ge\gamma=%.2f$." % gate,
            r"\texttt{resid}: drop the old batch. \texttt{oracle}: knows the cut.}",
            r"\label{tab:po-refit-mse-%s}" % learner,
            r"\small",
            r"\setlength{\tabcolsep}{4pt}",
            r"\begin{tabular}{@{}lccccc@{}}\toprule",
            r"Scene & uniform & rfperm & drop-old & oracle & fire \\",
            r"\midrule",
        ]
        for scene, cell in summary[learner].items():
            m = cell["methods"]
            vals = [
                m["uniform_pair"]["mse_mean"],
                m["rfperm"]["mse_mean"],
                m["resid"]["mse_mean"],
            ]
            lines.append(
                r"%s & %s & %s & %s & $%.3f$ & $%.2f$ \\"
                % (
                    scene,
                    _bold_min(vals, m["uniform_pair"]["mse_mean"], "%.3f"),
                    _bold_min(vals, m["rfperm"]["mse_mean"], "%.3f"),
                    _bold_min(vals, m["resid"]["mse_mean"], "%.3f"),
                    m["oracle"]["mse_mean"],
                    m["rfperm"]["fire_mean"],
                )
            )
        lines.extend([r"\bottomrule", r"\end{tabular}", r"\end{table}", ""])
        if "concept" not in summary[learner]:
            continue
        m = summary[learner]["concept"]["methods"]
        hops = sorted(int(t) for t in m["uniform_pair"].get("by_t", {}))
        if not hops:
            continue
        lines += [
            r"\begin{table}[ht]\centering",
            r"\caption{Concept hop path (\texttt{%s}), cut at $B_4$." % learner,
            r"Next-batch MSE. The cut hop is unforecastable; annealing pays on the hop after.}",
            r"\label{tab:po-refit-concept-%s}" % learner,
            r"\small",
            r"\setlength{\tabcolsep}{4pt}",
            r"\begin{tabular}{@{}clcccc@{}}\toprule",
            r"$t$ & test & uniform & rfperm & drop-old & oracle \\",
            r"\midrule",
        ]
        for t in hops:
            mark = r" $\leftarrow$ cut" if t + 1 == 4 else ""
            row_vals = [
                m["uniform_pair"]["by_t"][str(t)],
                m["rfperm"]["by_t"][str(t)],
                m["resid"]["by_t"][str(t)],
            ]
            lines.append(
                r"%d & $B_{%d}$%s & %s & %s & %s & $%.3f$ \\"
                % (
                    t,
                    t + 1,
                    mark,
                    _bold_min(row_vals, row_vals[0], "%.3f"),
                    _bold_min(row_vals, row_vals[1], "%.3f"),
                    _bold_min(row_vals, row_vals[2], "%.3f"),
                    m["oracle"]["by_t"][str(t)],
                )
            )
        lines.extend([r"\bottomrule", r"\end{tabular}", r"\end{table}", ""])
    return lines


def _real_index(rows):
    learners = []
    datasets = []
    for r in rows:
        if r["learner"] not in learners:
            learners.append(r["learner"])
        if r["dataset"] not in datasets:
            datasets.append(r["dataset"])
    return learners, datasets


def _real_score(row):
    v = float(row["online_score"])
    if row["task"] == "mse":
        return v ** 0.5
    return v


def real_tables(rows, gate, n_per):
    learners, datasets = _real_index(rows)
    lines = [
        r"%% Real consecutive-batch clocks. Batch size %d, $\gamma=%.2f$."
        % (n_per, gate),
        r"% RMSE $\downarrow$ on continuous tasks; Acc $\uparrow$ on discrete.",
        r"% fire = rfperm consecutive-OOS fire rate (first hop never fires).",
        "",
    ]
    for learner in learners:
        lines += [
            r"\begin{table}[ht]\centering",
            r"\caption{Real consecutive clocks (\texttt{%s}), batch size %d,"
            % (learner, n_per),
            r"24 hops. Last-two uniform is the default. \texttt{rfperm} fires",
            r"iff $e_{\mathrm{now}}/e_{\mathrm{prev}}\ge\gamma=%.2f$." % gate,
            r"Continuous: RMSE ($\downarrow$). Discrete: Acc ($\uparrow$).}",
            r"\label{tab:po-refit-real-%s}" % learner,
            r"\small",
            r"\setlength{\tabcolsep}{3.5pt}",
            r"\begin{tabular}{@{}llcrrrr@{}}\toprule",
            r"Dataset & clock & metric & uniform & rfperm & drop-old & fire \\",
            r"\midrule",
        ]
        for ds in datasets:
            sub = [r for r in rows if r["learner"] == learner and r["dataset"] == ds]
            if not sub:
                continue
            by_m = {r["method"]: r for r in sub}
            task = sub[0]["task"]
            clock = sub[0].get("clock") or ""
            metric = "RMSE" if task == "mse" else "Acc"
            scores = {
                m: _real_score(by_m[m])
                for m in ("uniform_pair", "rfperm", "resid")
            }
            vals = [scores["uniform_pair"], scores["rfperm"], scores["resid"]]
            fmt = _rmse_fmt(scores["uniform_pair"]) if task == "mse" else "%.3f"
            bold = _bold_min if task == "mse" else _bold_max
            lines.append(
                r"%s & %s & %s & %s & %s & %s & $%.2f$ \\"
                % (
                    ds.replace("_", r"\_"),
                    clock,
                    metric,
                    bold(vals, scores["uniform_pair"], fmt),
                    bold(vals, scores["rfperm"], fmt),
                    bold(vals, scores["resid"], fmt),
                    by_m["rfperm"]["fire_rate"],
                )
            )
        lines.extend([r"\bottomrule", r"\end{tabular}", r"\end{table}", ""])
    extra = _periodic_hop_tables(rows)
    if extra:
        lines += extra
    return lines


def _periodic_hop_tables(rows, names=("occupancy", "beijing_pm25")):
    """Fire hops on periodic clocks — the leftover failure mode."""
    lines = []
    for ds in names:
        for learner in ("rf", "xgb"):
            sub = [
                r
                for r in rows
                if r["dataset"] == ds
                and r["learner"] == learner
                and r["method"] == "rfperm"
            ]
            if not sub:
                continue
            rec = sub[0]
            hist = rec.get("history") or []
            fires = [h for h in hist if h.get("fired")]
            if not fires:
                continue
            metric = "Acc" if rec["task"] == "acc" else "RMSE"
            by_m = {
                r["method"]: {int(h["t"]): h for h in (r.get("history") or [])}
                for r in rows
                if r["dataset"] == ds and r["learner"] == learner
            }
            ts = sorted({int(h["t"]) for h in fires})
            neighbors = []
            for t in ts:
                for u in (t, t + 1):
                    if u not in neighbors and any(
                        u in by_m.get(m, {})
                        for m in ("uniform_pair", "rfperm", "resid")
                    ):
                        neighbors.append(u)
            lines += [
                r"\begin{table}[ht]\centering",
                r"\caption{Periodic reheating on \texttt{%s} (\texttt{%s})."
                % (ds.replace("_", r"\_"), learner),
                r"Fire hops and the hop after (anneal). Metric: %s.}" % metric,
                r"\label{tab:po-refit-periodic-%s-%s}"
                % (ds.replace("_", "-"), learner),
                r"\small",
                r"\begin{tabular}{@{}ccccc@{}}\toprule",
                r"$t$ & fired & uniform & rfperm & drop-old \\",
                r"\midrule",
            ]
            for t in neighbors:
                cells = []
                for method in ("uniform_pair", "rfperm", "resid"):
                    h = by_m.get(method, {}).get(t)
                    if h is None:
                        cells.append(r"---")
                        continue
                    sc = float(h["score"])
                    if rec["task"] == "mse":
                        sc = sc ** 0.5
                        cells.append("$" + (_rmse_fmt(sc) % sc) + "$")
                    else:
                        cells.append(r"$%.3f$" % sc)
                fired = by_m.get("rfperm", {}).get(t, {}).get("fired", False)
                mark = r"yes" if fired else r"no"
                lines.append(
                    r"%d & %s & %s & %s & %s \\"
                    % (t, mark, cells[0], cells[1], cells[2])
                )
            lines.extend([r"\bottomrule", r"\end{tabular}", r"\end{table}", ""])
    return lines


def write_synth_tex(summary, gate, learners, path):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(synth_tables(summary, gate, learners)), encoding="utf-8")
    return path


def write_real_tex(rows, gate, n_per, path):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(real_tables(rows, gate, n_per)), encoding="utf-8")
    return path


def write_all_tex(synth, real, path):
    """One paste-ready bundle: synth MSE + concept hops + real clocks."""
    lines = [
        r"% Paste-ready AGOD performance comparison.",
        r"% Requires booktabs. Last-two uniform default; rfperm = gated",
        r"% last-two $\sqrt{\mathrm{po\_risk0}}$ on $T=1$; resid = drop-old.",
        "",
    ]
    if synth:
        lines += synth_tables(
            synth["summary"],
            synth["gate"],
            synth["learners"],
        )
    if real:
        lines += real_tables(
            real["rows"],
            real["gate"],
            real.get("n_per", 200),
        )
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines), encoding="utf-8")
    return path
