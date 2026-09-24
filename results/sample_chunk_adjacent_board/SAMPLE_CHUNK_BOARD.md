# Sample-chunk adjacent boards (multi-dataset)

Sort → every N samples as a window → adjacent chunk FS board. Business grain = sample count (1000/2000); joint viz only.

## diffusiondb
- meta: `{'dataset': 'diffusiondb', 'y': 'image_nsfw', 'sort': 'timestamp', 'n': 8000}`
- ![diffusiondb](diffusiondb_board.png)

### every 1000 samples
- pairs=7

| t0→t1 | n0/n1 | ΔȲ | AUC | top_fsds |
|---|---|---:|---:|---|
| 0→1 | 1000/1000 | 0.0054 | 0.602 | mucha, alphonse mucha, alphonse |
| 1→2 | 1000/1000 | -0.0040 | 0.636 | alphonse, alphonse mucha, mucha |
| 2→3 | 1000/1000 | 0.0246 | 0.610 | alphonse, alphonse mucha, mucha |
| 3→4 | 1000/1000 | -0.0088 | 0.604 | artgerm, mucha, by artgerm |
| 4→5 | 1000/1000 | 0.0214 | 0.594 | mucha, girl, female |
| 5→6 | 1000/1000 | 0.0374 | 0.579 | dress, white, woman |
| 6→7 | 1000/1000 | 0.0280 | 0.609 | fashion, body, young |

### every 2000 samples
- pairs=3

| t0→t1 | n0/n1 | ΔȲ | AUC | top_fsds |
|---|---|---:|---:|---|
| 0→1 | 2000/2000 | 0.0110 | 0.611 | alphonse mucha, alphonse, mucha |
| 1→2 | 2000/2000 | 0.0142 | 0.603 | mucha, alphonse, alphonse mucha |
| 2→3 | 2000/2000 | 0.0621 | 0.613 | dress, woman, mucha |

## tencent_gr
- meta: `{'dataset': 'tencent_gr_edges', 'y': 'y_convert', 'sort': 'e_last_ts', 'n': 8000, 'd': 35}`
- ![tencent_gr](tencent_gr_board.png)

### every 1000 samples
- pairs=6

| t0→t1 | n0/n1 | ΔȲ | AUC | top_fsds |
|---|---|---:|---:|---|
| 0→1 | 1000/1000 | 0.0090 | 0.990 | e_log1p_exp, e_n_exp, i_log1p_n_exp |
| 1→2 | 1000/1000 | -0.0270 | 1.000 | e_log1p_exp, e_n_exp, i_log1p_n_exp |
| 2→3 | 1000/1000 | -0.0150 | 0.999 | e_n_exp, e_log1p_exp, u_n_exp |
| 3→4 | 1000/1000 | 0.0040 | 1.000 | e_log1p_exp, e_n_exp, u_n_exp |
| 4→5 | 1000/1000 | 0.0010 | 0.997 | e_log1p_exp, e_n_exp, i_credit_last |
| 5→6 | 1000/1000 | -0.0190 | 0.989 | e_log1p_exp, e_n_exp, i_share_last |

### every 2000 samples
- pairs=3

| t0→t1 | n0/n1 | ΔȲ | AUC | top_fsds |
|---|---|---:|---:|---|
| 0→1 | 2000/2000 | -0.0300 | 1.000 | e_log1p_exp, e_n_exp, i_log1p_n_exp |
| 1→2 | 2000/2000 | -0.0030 | 0.997 | e_n_exp, e_log1p_exp, u_n_exp |
| 2→3 | 2000/2000 | -0.0190 | 0.999 | e_log1p_exp, e_n_exp, i_share_last |

## waymo_proxy
- meta: `{'dataset': 'waymo_proxy', 'y': 'proxy_y', 'sort': 'row_index', 'n': 8000, 'd': 9}`
- ![waymo_proxy](waymo_proxy_board.png)

### every 1000 samples
- pairs=7

| t0→t1 | n0/n1 | ΔȲ | AUC | top_fsds |
|---|---|---:|---:|---|
| 0→1 | 1000/1000 | 0.2131 | 0.943 | x0, x8, x3 |
| 1→2 | 1000/1000 | 0.2036 | 0.928 | x0, x8, x3 |
| 2→3 | 1000/1000 | 0.2359 | 0.937 | x0, x8, x3 |
| 3→4 | 1000/1000 | 0.2322 | 0.941 | x0, x8, x3 |
| 4→5 | 1000/1000 | 0.2373 | 0.938 | x0, x8, x5 |
| 5→6 | 1000/1000 | 0.2465 | 0.944 | x0, x8, x3 |
| 6→7 | 1000/1000 | 0.2465 | 0.947 | x0, x8, x3 |

### every 2000 samples
- pairs=3

| t0→t1 | n0/n1 | ΔȲ | AUC | top_fsds |
|---|---|---:|---:|---|
| 0→1 | 2000/2000 | 0.4281 | 0.963 | x0, x8, x5 |
| 1→2 | 2000/2000 | 0.4688 | 0.972 | x0, x8, x5 |
| 2→3 | 2000/2000 | 0.4884 | 0.973 | x0, x8, x5 |

## metro_interstate
- meta: `{'dataset': 'metro_interstate', 'y': 'traffic_volume', 'sort': 'date_time', 'n': 8000, 'd': 19}`
- ![metro_interstate](metro_interstate_board.png)

### every 1000 samples
- pairs=7

| t0→t1 | n0/n1 | ΔȲ | AUC | top_fsds |
|---|---|---:|---:|---|
| 0→1 | 1000/1000 | -267.1500 | 0.954 | m5, m0, m14 |
| 1→2 | 1000/1000 | -218.6630 | 0.964 | m5, m9, m3 |
| 2→3 | 1000/1000 | 243.3420 | 0.970 | m5, m6, m12 |
| 3→4 | 1000/1000 | -3.7380 | 0.974 | m5, m0, m16 |
| 4→5 | 1000/1000 | 171.8230 | 0.983 | m5, m14, m10 |
| 5→6 | 1000/1000 | 42.4760 | 0.984 | m5, m0, m9 |
| 6→7 | 1000/1000 | -110.6580 | 0.974 | m5, m1, m9 |

### every 2000 samples
- pairs=3

| t0→t1 | n0/n1 | ΔȲ | AUC | top_fsds |
|---|---|---:|---:|---|
| 0→1 | 2000/2000 | -230.5670 | 0.968 | m5, m0, m6 |
| 1→2 | 2000/2000 | 203.8445 | 0.980 | m5, m6, m0 |
| 2→3 | 2000/2000 | 73.0585 | 0.979 | m5, m0, m6 |

## beijing_pm25
- meta: `{'dataset': 'beijing_pm25', 'y': 'pm2.5', 'sort': 'stamp', 'n': 8000, 'd': 13}`
- ![beijing_pm25](beijing_pm25_board.png)

### every 1000 samples
- pairs=7

| t0→t1 | n0/n1 | ΔȲ | AUC | top_fsds |
|---|---|---:|---:|---|
| 0→1 | 1000/1000 | 12.7520 | 0.971 | p8, p0, p2 |
| 1→2 | 1000/1000 | -12.2280 | 0.986 | p8, p2, p0 |
| 2→3 | 1000/1000 | 8.3620 | 0.957 | p8, p0, p10 |
| 3→4 | 1000/1000 | 24.0380 | 0.982 | p8, p0, p7 |
| 4→5 | 1000/1000 | -5.5500 | 0.968 | p8, p0, p10 |
| 5→6 | 1000/1000 | -5.0780 | 0.993 | p8, p0, p10 |
| 6→7 | 1000/1000 | 28.0400 | 0.989 | p8, p0, p2 |

### every 2000 samples
- pairs=3

| t0→t1 | n0/n1 | ΔȲ | AUC | top_fsds |
|---|---|---:|---:|---|
| 0→1 | 2000/2000 | -1.6710 | 0.971 | p8, p0, p2 |
| 1→2 | 2000/2000 | 25.4440 | 0.977 | p8, p0, p10 |
| 2→3 | 2000/2000 | 6.1670 | 0.991 | p8, p0, p10 |

