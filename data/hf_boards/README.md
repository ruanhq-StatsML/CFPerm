Labeled Hugging Face boards used by `scripts/run_hf_block_aware_board.py`.

Parquet is streamed at runtime from Hub `refs/convert/parquet` (not committed):

- `scikit-learn/adult-census-income` — income label, sex as batch; blocks demography / work / hours / capital
- `fancyzhx/yelp_polarity` vs `fancyzhx/amazon_polarity` — sentiment, marketplace as batch; Amazon rows via Dataset Viewer
- `nyu-mll/multi_nli` validation_matched — entailment, fiction vs telephone
