# HF landing prototypes — run report

Packages: local `agod` OnlineRFPerm + PO-risk; sklearn HashingVectorizer RF.
Datasets: `Anthropic/hh-rlhf` (helpful-base subset), `pminervini/HaluEval` (qa_samples).

## 1) RFPerm-as-Judge + 文风 / preference hop (HH-RLHF)

- judge err pre→post: `0.147` → `0.539` (ratio `3.67`)
- style domain AUC (P(X) register): `0.999`
- text domain AUC: `0.494`
- fire_rate `0.104`, first_fire_t `4`
- hop@cut fired=`True` ratio=`1.342857142857143`

Read: high style AUC ⇒ 文风画像变了；judge ratio / fire ⇒ 偏好映射变了（DPO/RLHF 连更）。
动作：fire 后按 `po_risk0` Top-k 抽检偏好对；安静回均匀。

## 2) Hallucination regime + RAG-hit + router (HaluEval)

- halluc rate before/after: `0.080` / `0.942`
- fire_rate `0.091`, first_fire_t `4`
- instance probe AUROC(post): `1.0`
- hop@cut fired=`True` ranking=`{'n_t1': 100, 'auroc_po_risk0': 1.0, 'precision_at_10': 1.0, 'mean_po_t1': 1.041720734126984, 'halluc_rate': 0.82, 'mean_rag_hit_top10': 1.0}`
- router_task_shift: `{'task0_rate_pre': 0.54, 'task0_rate_post': 0.5665, 'halluc_by_task_post': {'task0': 0.941747572815534, 'task1': 0.9423298731257209}, 'domain_auc_pre_post': 0.811158934252631}`

Read: regime fire = 幻觉标签制度跳变；`po_risk0` 只排序；RAG-hit 低 ⇒ 检索支持缺口；
task_shift ⇒ 下一步刷新哪个 router 专家。

## Stance

Detection / localization / ranking / gated reweight — not fact-checking, not unique causal attribution.
