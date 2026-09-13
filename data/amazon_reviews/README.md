# Amazon Reviews shards (local cache)

Download from Hugging Face:
[`jingxiang11111/amazon_reviews_for_rec`](https://huggingface.co/datasets/jingxiang11111/amazon_reviews_for_rec)

```bash
mkdir -p data/amazon_reviews/shards
# example: a few train shards for smoke
curl -L -o data/amazon_reviews/shards/data-000000-00cd13e1.tar.gz \
  https://huggingface.co/datasets/jingxiang11111/amazon_reviews_for_rec/resolve/main/train/data-000000-00cd13e1.tar.gz
```

Do **not** commit `.tar.gz` shards to git.
