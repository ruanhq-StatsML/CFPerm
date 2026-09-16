# 大模型审核：检测审核逻辑的突变

场景不是「模型突然不会了」，也不是「幻觉率变高就不丝滑」。
这里检测的是 **审核逻辑** \(P(Y\mid X)\) 有没有在相邻窗里 **一刀切地换制度**。

\(X\)：待审样本（query、回复、可选上下文）。
\(Y\)：审核决定（放行 / 拒绝 / 转人工，或 RLHF 的 chosen/rejected）。
上一窗拟合的 probe \(\mu_0\) 就是当时的审核逻辑。

**Objective 是映射有没有 hop**（RFPerm / PO 当分布距离、置换检验看预测退化）。  
**不是 ratio，也不是 AUC。**

| 东西 | 角色 | 跟 objective |
|---|---|---|
| last-two fire / 对参考的 RFPerm | 审核映射变没变 | **就是 objective** |
| \(e_{\mathrm{now}}/e_{\mathrm{prev}}\) | 工程 gate，用来扣扳机 | 无关：阈值 \(\gamma\) 可换，不是要优化的量 |
| domain AUC / `po_risk0` AUROC | 旁证（画像、排序好不好用） | **不一致**：那是分类可分性，不是置换距离 |
| 拒绝率 / 幻觉率 / pos-rate | 水平高低 | 无关：高可以仍然丝滑 |

quiet vs fire 才是读法。ratio=1.34 或 AUC=0.999 都不构成这条方法的目标函数。

## 1. 两条线不要混

| 问题 | 统计层 | 审核侧读法 |
|---|---|---|
| 审核链路还顺吗 | 相邻窗 last-two hop | 同一套政策、同一套 judge，连续两批决定还像不像 |
| 审核还信得过吗 | 相对上线参考的退化 | 相对冻结金标 / 上线时的 probe，现在已经漂了多少 |

**逐渐变严、逐渐变松**：拒绝率 10% → 12% → 14%，相邻窗仍像，链路算丝滑。last-two 可以不 fire。相对原参考的退化检验可以在「审核标准已经偏到不能当金标」之前响。

**死光了 / 奇特逻辑**：政策一夜改版、judge 模型换代、系统 prompt 重写、标注团队整体替换、某一类越狱流量突然占满窗口。相邻窗上 probe 对不上新决定，这才是变点。ratio 只是扣扳机的闸。

高拒绝率本身不是变点。稳定地严、稳定地松，都还可以很顺。断的是 **映射突然换了**。

## 2. 线上长什么样

按时间切 batch（小时 / 日 / 发布窗）：

```
probe ← fit(X_{t-1}, Y_{t-1})          # 上一窗的审核逻辑
e_now ← err(probe, X_t, Y_t)
fire  ← (e_now / e_prev ≥ γ)           # 默认 γ=1.25
po_i  ← |Y_i - μ0(X_i)|                # 相对旧审核逻辑的残差
audit ← argsort(-po)[:k]               # 只在 fire 后打开人工对照窗
```

| 状态 | 含义 | 动作 |
|---|---|---|
| quiet | 审核逻辑连续 | 标准路径；DPO/RLHF 数据可按原门禁合并；\(w=1\) |
| fire | 审核逻辑 hop | **不要把当前 judge 当金标**；Top-\(k\) `po_risk0` 人工复审；拒合并或限量合并偏好对；可选 \(\sqrt{\mathrm{PO}}\) 只打在 \(T=1\) |

Fire 打开的是对照窗，不是「这条违规了」的分类器。

## 3. 已有一跑：HH-RLHF，RFPerm-as-Judge

数据：`Anthropic/hh-rlhf` helpful-base。\(n=4000\)，50 窗 × 80 条，cut 在 \(B_4\)。Cite：`bai2022hh-rlhf`。

读法只看 **有没有 fire**（映射 hop 了没有）：

| 窗 | fire | 读法 |
|---|---|---|
| 安静 \(t=2\) | 否 | 审核逻辑还顺 |
| **cut \(t=4\)** | **是**（delay 0） | 偏好映射 hop |
| 之后 | 不持续当金标 | 新审核制度可以重新变顺 |

Gate 上留下的 ratio（1.03 / 1.34）、judge-err 3.67、style AUC 0.999、text AUC 0.494，都是 **日志**。AUC 说的是 \(P(X)\) 好不好分，跟「映射距离 / 置换检验」不是同一个 objective。Style AUC 高只说明语域也动了，不能拿来当这条方法的成绩。

动作仍然：quiet 放行；fire 后 `po_risk0` Top-\(k\) 抽偏好对，不把新 judge 当金标。`po_risk0` 是 hop **之后** 的排序，不是 objective。

对照：HaluEval 幻觉制度同样是 cut 当窗 fire。那边 gate 更宽，不代表 objective 变成了更大的 ratio。

## 4. 审核场景里什么算「突变」

具体会撞上的触发，而不是「模型变蠢」：

1. **政策包切版**：新安全红线、新地区合规、某类内容从 warn 改成 block。
2. **Judge 换代**：外挂 LLM 审核、内部分类器、规则引擎权重一夜间替换。
3. **提示词 / 系统策略重写**：同一批样本，旧 probe 对不上新决定。
4. **人审池更换**：外包团队、指南、抽检比例整体换。
5. **流量成分突变叠在审核上**：越狱、刷评、某产品线突然放量。若只有 mix 变、映射没变，是画像 hop 不是审核逻辑 hop。domain AUC 可以当旁证，**不能替代** fire。

HH 这场：last-two **fire 了** ⇒ 审核映射当变点处理。AUC 再高也只是旁证。

## 5. 逐渐漂 vs 突然死

- 审核员越来越严、每周紧一点：链路仍可丝滑。last-two 不响是对的。相对上线参考的退化可以早响——「这套审核已经不是上线时那套了」，在彻底不能用之前拦下来。
- 某天政策/judge 一刀切：这就是死光了，变点。HH 的 `first_fire_t = cut = 4` 是这个。

两层一起才稳健：不把温水误报成崩盘，也不把真 hop 当成慢慢变差。

## 6. 不声称什么

- 不是道德审核器，不输出「这条该不该过」。
- 不是事实核查。
- 不优化、不汇报 ratio 或 AUC 当主结果。
- `po_risk0` 不是违规分数，是相对旧 probe 的残差强度，只用于 fire 后排序。
- 不能把 hop 归因到某一个词、某一条规则。

## 文件

- 本说明：`docs/reports/LLM_Judge_Audit_Hop.md`
- 可贴表：`docs/reports/LLM_Judge_Audit_Hop_tables_only.tex`
- 数字来源：HF landing `hh_judge_style.json`（branch `cursor/hf-llm-landing-protos-92f8`）
