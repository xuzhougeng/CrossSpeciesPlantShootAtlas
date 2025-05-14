

TF_prediction.sh is used for de novo annotation of transcription factors across different species.

## Integration

The notebook in integration follows the [SCGLUE tutorial](https://scglue.readthedocs.io/en/latest/) to perform joint analysis of scATAC-seq and scRNA-seq data, along with regulatory inference.
This process generates `gene2peaks.links` and `gene2peaks.bedpe` files as outputs.


# Resources

gene2symbol.txt : TF geneid and its symbol

ath_cisbp_all.bed: 

```bash
python -m inferelator_prior.pwm_to_meme --motif pwms_all_motifs/* --info TF_Information_all_motifs_plus.txt --out ath.cisbp.meme
python -u run_fimo.py -f Athaliana_447_TAIR10.fa -m ath.cisbp.meme -o ath_cisbp_all.csv -j 7 -a '--thresh 1e-5'
```

Run following code to convert csv to bed

```python
input = "ath_cisbp_all.csv"
output = "ath_cisbp_all.bed"
import pandas as pd
fimo = pd.read_csv(input, usecols=["motif_alt_id", "sequence_name", "start", "stop"])
fimo.columns = ["name", "chrom", "chromStart", "chromEnd"]
fimo = fimo.loc[:, ["chrom", "chromStart", "chromEnd", "name"]]
fimo.sort_values(["chrom", "chromStart"], inplace=True)
fimo.to_csv(output, sep="\t", header=False, index=False)
```

