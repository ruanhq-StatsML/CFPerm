# ChannelAttribution-style + graph localization

- Heuristic: first / last / linear (CA-compatible)
- Markov-1 removal-effect on **transition digraph**
- Localization: `networkx.ego_graph` + Personalized PageRank
- Graphs: transition digraph **and** co-occurrence (windowed)
- Figure: `tencent_channel_graph_attribution.png`

ChannelAttribution pip wheel failed to build (Cython); heuristic+Markov-removal reimplemented to match CA API; localization via networkx.ego_graph + pagerank.
