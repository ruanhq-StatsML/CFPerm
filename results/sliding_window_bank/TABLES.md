# Sliding-window candidate feature bank

Cheap rolling stats vs frozen D_ref. Online PCA is the fast P(X) pointer.
MMD still vs D_ref, not pairwise history.

| kind | t | onset | T | Brier excess | PCA recon | PCA video | PCA audio | PCA text | argmax |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| concept | 0 |  | 0.142 | 0.0954 | 0.0339 | -0.112 | 0.0766 | 0.0299 | audio |
| concept | 1 |  | 0.112 | 0.0728 | 0.0258 | -0.0685 | 0.0462 | -0.00969 | audio |
| concept | 2 |  | 0.091 | 0.0514 | 0.064 | -0.0208 | 0.0792 | 0.0762 | audio |
| concept | 3 | yes | 0.0739 | 0.0321 | 0.1 | 0.0247 | 0.0789 | 0.0956 | text |
| concept | 4 | yes | 0.0873 | 0.0349 | 0.121 | 0.0754 | 0.0843 | 0.124 | text |
| concept | 5 | yes | 0.104 | 0.0474 | 0.107 | 0.075 | 0.13 | 0.0478 | audio |
| concept | 6 | yes | 0.113 | 0.0506 | 0.0723 | 0.0655 | 0.109 | 0.0443 | audio |
| concept | 7 | yes | 0.128 | 0.0687 | 0.0629 | 0.0339 | 0.129 | 0.0313 | audio |
| concept | 8 | yes | 0.152 | 0.0857 | 0.0353 | 0.00121 | 0.0514 | 0.031 | audio |
| concept | 9 | yes | 0.168 | 0.116 | 0.0623 | 0.0563 | 0.0937 | -0.00921 | audio |
| concept | 10 | yes | 0.211 | 0.169 | 0.0524 | 0.0502 | 0.0893 | -0.00584 | audio |
| concept | 11 | yes | 0.24 | 0.22 | 0.0854 | 0.06 | 0.0996 | -0.00155 | audio |
| covariate | 0 |  | 0.114 | 0.0413 | 0.0339 | -0.112 | 0.0766 | 0.0299 | audio |
| covariate | 1 |  | 0.0865 | 0.0162 | 0.0258 | -0.0685 | 0.0462 | -0.00969 | audio |
| covariate | 2 |  | 0.0691 | 0.00416 | 0.064 | -0.0208 | 0.0792 | 0.0762 | audio |
| covariate | 3 | yes | 0.0627 | 0.0083 | 0.103 | 0.0247 | 0.0834 | 0.0956 | text |
| covariate | 4 | yes | 0.0783 | 0.0322 | 0.125 | 0.0754 | 0.0864 | 0.124 | text |
| covariate | 5 | yes | 0.0872 | 0.0418 | 0.114 | 0.075 | 0.144 | 0.0478 | audio |
| covariate | 6 | yes | 0.0934 | 0.0265 | 0.0887 | 0.0655 | 0.138 | 0.0443 | audio |
| covariate | 7 | yes | 0.086 | 0.0149 | 0.102 | 0.0339 | 0.23 | 0.0313 | audio |
| covariate | 8 | yes | 0.0975 | 0.0229 | 0.0778 | 0.00121 | 0.159 | 0.031 | audio |
| covariate | 9 | yes | 0.0783 | 0.0315 | 0.143 | 0.0563 | 0.327 | -0.00921 | audio |
| covariate | 10 | yes | 0.0731 | 0.0305 | 0.18 | 0.0502 | 0.417 | -0.00584 | audio |
| covariate | 11 | yes | 0.0522 | 0.0165 | 0.291 | 0.06 | 0.593 | -0.00155 | audio |
