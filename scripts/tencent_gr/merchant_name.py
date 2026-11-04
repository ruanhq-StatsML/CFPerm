#!/usr/bin/env python3
"""Merchant names: generate vs encode are two different jobs.

generate_merchant_name(seed, i)  — 造一个英文店名（标签，不是特征）
name_matrix(names)               — 真要把店名塞进 X 再用（TF-IDF / hashing）

造名不要 TF-IDF。TF-IDF 是「已经有字符串之后」的稀疏编码。
LabelEncoder 对 400 个互异店名 = 400 个哑变量，店和店之间不共享 token。

  python3 scripts/tencent_gr/merchant_name.py
"""
from __future__ import annotations

from typing import Sequence

ADJ = (
    "Amber Bright Cedar Delta Ember Frost Golden Harbor Ivory Jade "
    "Kinetic Lumen Maple North Oak Pearl Quartz Ridge Silver Tide"
).split()
NOUN = (
    "Basket Bazaar Cart Depot Emporium Forge Goods House Imports Kettle "
    "Mart Outlet Plaza Rack Supply Trader Union Vault Works Yard"
).split()
SFX = "Co Ltd Shop Store Collective Trading Studio Market".split()


def generate_merchant_name(seed: int, i: int = 0, backend: str = "faker") -> str:
    """One English shop name. `i` = merchant index so a seed yields a catalog.

    backends
    --------
    faker     faker.Faker.company()          实务默认，已装
    template  ADJ+NOUN+SFX                   无依赖，可复现
    messy     随机拉丁串                     顶替 vLLM SamplingParams
    """
    b = backend.lower()
    if b == "faker":
        return _faker(seed, i)
    if b == "template":
        return f"{ADJ[i % len(ADJ)]} {NOUN[(i * 7) % len(NOUN)]} {SFX[(i * 13) % len(SFX)]}"
    if b == "messy":
        import numpy as np

        rng = np.random.default_rng(int(seed) * 1_000_003 + int(i))
        k = int(rng.integers(6, 14))
        alpha = np.array(list("abcdefghijklmnopqrstuvwxyz"))
        return "".join(rng.choice(alpha, size=k))
    raise ValueError(f"unknown backend {backend!r}; use faker|template|messy")


def generate_catalog_names(n: int, seed: int, backend: str = "faker") -> list[str]:
    return [generate_merchant_name(seed, i, backend=backend) for i in range(n)]


def _faker(seed: int, i: int) -> str:
    from faker import Faker

    fake = Faker("en_US")
    fake.seed_instance(int(seed) * 1_000_003 + int(i))
    return fake.company()


def name_matrix(names: Sequence[str], *, kind: str = "tfidf", k: int = 8):
    """把店名编成稠密矩阵。只在店名要进 X 时用。

    tfidf   sklearn TfidfVectorizer char_wb 3-5 → TruncatedSVD(k)
    hash    HashingVectorizer n_features=k  （不存词表，可在线）
    """
    import numpy as np
    from sklearn.decomposition import TruncatedSVD
    from sklearn.feature_extraction.text import HashingVectorizer, TfidfVectorizer

    docs = [str(x) for x in names]
    if kind == "hash":
        hv = HashingVectorizer(
            n_features=k, alternate_sign=False, ngram_range=(1, 2), analyzer="word"
        )
        return hv.fit_transform(docs).toarray()
    tf = TfidfVectorizer(analyzer="char_wb", ngram_range=(3, 5), min_df=1)
    x = tf.fit_transform(docs)
    kk = min(int(k), max(1, x.shape[1] - 1), max(1, x.shape[0] - 1))
    if kk < 1:
        return np.zeros((len(docs), k))
    svd = TruncatedSVD(n_components=kk, random_state=0)
    z = svd.fit_transform(x)
    if z.shape[1] < k:
        z = np.hstack([z, np.zeros((z.shape[0], k - z.shape[1]))])
    return z


def _demo() -> None:
    print("generate")
    for b in ("faker", "template", "messy"):
        print(f"  {b:10s}", [generate_merchant_name(0, i, backend=b) for i in range(3)])
    names = generate_catalog_names(6, 0, "faker")
    z = name_matrix(names, kind="tfidf", k=4)
    print("tfidf", z.shape, names[:2])


if __name__ == "__main__":
    _demo()
