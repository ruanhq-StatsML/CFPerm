# LOGO trimodal probe

Shares are localization proxies, not a unique decomposition.
Layer 1 = Brier LOGO. Layer 2 = PO-risk / MMD mix. T is the batch label.

n_ref=640, n_new=160, batches=5, onset=2.
Shifted slice is always 3.

## Per-batch two-layer shares and next-batch plan

| kind | t | onset | action | update | Brier | π_loss v/a/t | π_PO v/a/t | π_MMD v/a/t | mix_PO v/a/t | video tower | audio tower | text tower | fusion | top subset | subset π_MMD v/a/t | subset dominant |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| concept_video | 0 |  | keep_training | full_train | 0.197 | 0/1/0 | 0.0302/0.51/0.46 | 0/1/0 | 1/0.338/1 | full_train | full_train | full_train | fusion_full | 1 | 1/0/0 | text |
| concept_video | 1 |  | keep_training | full_train | 0.201 | 0/0/0 | 0/1/0 | 0/0.469/0.531 | 0/0.681/0 | full_train | full_train | full_train | fusion_full | 0 | 0.591/0/0.409 | video |
| concept_video | 2 | yes | keep_training | full_train | 0.177 | 0/0/0 | 0.104/0.0376/0.858 | 0/1/0 | 1/0.0363/1 | full_train | full_train | full_train | fusion_full | 3 | 0/1/0 | audio |
| concept_video | 3 | yes | keep_training | full_train | 0.233 | 1/0/0 | 0.8/0.0367/0.163 | 1/0/0 | 0.444/1/1 | full_train | full_train | full_train | fusion_full | 3 | 0.059/0.941/0 | video |
| concept_video | 4 | yes | keep_training | full_train | 0.227 | 0/0/0 | 0.467/0.231/0.302 | 1/0/0 | 0.318/1/1 | full_train | full_train | full_train | fusion_full | 3 | 0/0.0878/0.912 | video |
| covariate_audio | 0 |  | keep_training | full_train | 0.197 | 0/1/0 | 0.0302/0.51/0.46 | 0/1/0 | 1/0.338/1 | full_train | full_train | full_train | fusion_full | 1 | 1/0/0 | text |
| covariate_audio | 1 |  | keep_training | full_train | 0.201 | 0/0/0 | 0/1/0 | 0/0.469/0.531 | 0/0.681/0 | full_train | full_train | full_train | fusion_full | 0 | 0.591/0/0.409 | video |
| covariate_audio | 2 | yes | keep_training | full_train | 0.168 | 0/0/0 | 0/0.392/0.608 | 0/1/0 | 0/0.282/1 | full_train | full_train | full_train | fusion_full | 2 | 0/1/0 | text |
| covariate_audio | 3 | yes | keep_training | full_train | 0.21 | 0/0/0 | 0.639/0.0247/0.336 | 0/1/0 | 1/0.0241/1 | full_train | full_train | full_train | fusion_full | 2 | 0.888/0.013/0.0993 | video |
| covariate_audio | 4 | yes | keep_training | full_train | 0.169 | 0/0/0 | 0.0184/0/0.982 | 0/1/0 | 1/0/1 | full_train | full_train | full_train | fusion_full | 0 | 0/0/1 | text |
| both | 0 |  | keep_training | full_train | 0.197 | 0/1/0 | 0.0302/0.51/0.46 | 0/1/0 | 1/0.338/1 | full_train | full_train | full_train | fusion_full | 1 | 1/0/0 | text |
| both | 1 |  | keep_training | full_train | 0.201 | 0/0/0 | 0/1/0 | 0/0.469/0.531 | 0/0.681/0 | full_train | full_train | full_train | fusion_full | 0 | 0.591/0/0.409 | video |
| both | 2 | yes | keep_training | full_train | 0.174 | 0/0/0 | 0/0.305/0.695 | 0/1/0 | 0/0.234/1 | full_train | full_train | full_train | fusion_full | 2 | 0/1/0 | text |
| both | 3 | yes | keep_training | full_train | 0.221 | 0/0/0 | 0.831/0/0.169 | 0/1/0 | 1/0/1 | full_train | full_train | full_train | fusion_full | 3 | 0/1/0 | video |
| both | 4 | yes | keep_training | full_train | 0.209 | 0/0/0 | 0.297/0/0.703 | 0/1/0 | 1/0/1 | full_train | full_train | full_train | fusion_full | 3 | 0/1/0 | video |

## online Brier and cumulative regret

| kind | oracle | t | frozen | full update | LOGO policy | regret |
| --- | --- | --- | --- | --- | --- | --- |
| concept_video | full_update | 0 | 0.197 | 0.197 | 0.197 | 0 |
| concept_video | full_update | 1 | 0.201 | 0.201 | 0.201 | 0 |
| concept_video | full_update | 2 | 0.177 | 0.179 | 0.179 | 0 |
| concept_video | full_update | 3 | 0.233 | 0.228 | 0.228 | 0 |
| concept_video | full_update | 4 | 0.227 | 0.214 | 0.214 | 0 |
| covariate_audio | frozen | 0 | 0.197 | 0.197 | 0.197 | 0 |
| covariate_audio | frozen | 1 | 0.201 | 0.201 | 0.201 | 0.000305 |
| covariate_audio | frozen | 2 | 0.168 | 0.169 | 0.169 | 0.000851 |
| covariate_audio | frozen | 3 | 0.21 | 0.208 | 0.208 | -0.00118 |
| covariate_audio | frozen | 4 | 0.169 | 0.163 | 0.163 | -0.00688 |
| both | frozen | 0 | 0.197 | 0.197 | 0.197 | 0 |
| both | frozen | 1 | 0.201 | 0.201 | 0.201 | 0.000305 |
| both | frozen | 2 | 0.174 | 0.175 | 0.175 | 0.00137 |
| both | frozen | 3 | 0.221 | 0.218 | 0.218 | -0.00212 |
| both | frozen | 4 | 0.209 | 0.201 | 0.201 | -0.0102 |
