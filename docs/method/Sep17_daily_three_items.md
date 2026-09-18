# Sep17 daily wrap — 三件事状态

> 其他线（CFPerm dual/blend、PO/IPTW、obs-PO、image-OOD、HF efficiency 等）默认 **已到位**；今天主线只盯下面 3 项。

---

## 1. 监控 gradient（Grad-OnlineRFPerm）— **DONE / 可写进报告**

**口径（钉死）**
- 未冻结参数一个标量 \(g_t=\|\nabla_{\theta_U} L\|_2\)（不是 per-layer 多测、不是 norm 平均）
- **一整条** OnlineRFPerm：\(T_t=g_t-e_{\mathrm{ref}}\) → rank/EWMA \(p\) → FDR
- layer share 只做 reject 后 freeze-depth **诊断**

**结果（直观 lead-time，符合预期）**
- MVP 5 datasets × 5 seeds：mean lead \(\approx\mathbf{-3.24}\) batches；P(earlier)=**72%**
- electricity / eeg 最强；bank 持平偏温和
- Extras：grace 推迟假阳；α 方向稳定；electricity 上 freeze 闭环 **更低 MSE + 更少 FLOPs**；stocks_IWM lead −4

**交付物**
- 公式+实验 LaTeX：`docs/method/Grad_OnlineRFPerm.tex`
- 代码：`agod/grad_rfperm.py` · `scripts/run_agod_grad_rfperm_monitor.py` · `scripts/run_grad_rfperm_extras.py`
- 图：`results/agod_grad_rfperm/` · `results/grad_rfperm_extras/`
- PR：https://github.com/ruanhq-StatsML/CFPerm/pull/64

---

## 2. 推荐数据 batch / streaming testing — **部分完成，需定清单**

**已有可推荐 / 已跑过的 stream packs**
| 用途 | packs |
|---|---|
| MVP tabular streams | synthetic（控漂移）、Covertype、bank-marketing、electricity、eeg-eye-state |
| 生产向 hourly / market | metro_interstate、beijing_pm25、stocks_{MSFT,IWM,SPY,…}、waymo_proxy |
| 多模态 / 其他 | MSR-VTT、Affec、Amazon（既有 AGOD 线） |

**建议写进报告的「推荐测试矩阵」**
1. **Batch**：固定 \(f_{\mathrm{ref}}\)，batch=128/256，burn=8，报告 next-MSE / gate duty  
2. **Streaming**：连续 batch，对比 Grad vs MSE vs MMD vs PO 的 \(t_{\mathrm{reject}}\) / Lead  
3. **Null**：stationary synthetic + grace（Type-I / early FPR）  
4. **Freeze 闭环**：reject → freeze_early / freeze_low_share → MSE–FLOPs  

**还差（若今天要收口）**
- 一页「官方推荐 dataset 清单 + 默认超参」表（可从上面直接抽）
- 可选：把推荐矩阵写进 `docs/method/` 单独一节 / README

---

## 3. OnlinePermOOB + LLM / TabPFN software — **待启动（今天未落地代码）**

**目标形态（建议）**
- 同一 OnlinePermOOB / OnlineRFPerm **软件壳**
- 后端可插：现有 RF / MLP，以及 **TabPFN**（tabular）、**LLM**（text / multimodal features）
- API：`fit_ref` → stream `update(batch)` → `p_t, reject` → optional adapt / freeze

**今天状态**
- 核心 OnlineRFPerm 骨架已有：`agod/online_rfperm.py`
- Grad 监控已验证同一 FDR 壳可挂非 MSE 标量
- **尚无** TabPFN / LLM backend 包、安装脚本、统一 CLI

**下一步（最小切片）**
1. `OnlinePermOOB` facade：统一 `Backend` protocol（RF | TabPFN | LLM-embed+head）  
2. TabPFN path：小 tabular stream smoke（electricity / bank）  
3. LLM path：frozen embed → 轻量 head，Grad 或 MSE 标量进同一 OnlineRFPerm  
4. 一页 software README（install / one-liner / backends）

---

## 一句话给老板

今天：**gradient 监控（单流 Lead）做完且符合预期**；**streaming/batch 测试矩阵有数据与脚本可推荐**；**OnlinePermOOB×LLM/TabPFN 软件壳是下一刀**（骨架在，backend 未接）。其余线已就绪。
