# Automated ACMG/AMP Variant Pathogenicity Classification for TPP1 (CLN2 Disease): A Reproducible Bayesian Pipeline

Lena Traczuk
Dawid Fleischer



---

## Abstract

Automated, standardized, and reproducible variant classification is critical for clinical genetics. Here, we present **tpp1-cln2-pipeline**, an automated Snakemake workflow implementing the Tavtigian Bayesian point framework for the classification of variants in the *TPP1* gene (associated with late infantile neuronal ceroid lipofuscinosis / CLN2 disease). The pipeline integrates variant normalization, VEP v115, SpliceAI, gnomAD v4.1.0, and ClinVar. Analyzing 337 real ClinVar variants, the pipeline identified 63 clinically actionable (Pathogenic or Likely Pathogenic) variants, demonstrating 75.8% concordance with the commercial Franklin by Genoox platform.

---

## Introduction

Late infantile neuronal ceroid lipofuscinosis (CLN2 disease) is a rare autosomal recessive neurodegenerative disorder caused by mutations in the *TPP1* gene [1]. Accurate classification of *TPP1* variants according to the ACMG/AMP guidelines [2] is critical for early diagnosis and treatment selection (e.g., enzyme replacement therapy).

We present an automated, containerized pipeline that implements the Tavtigian Bayesian point system [3] with ClinGen-recommended gene-specific calibrations for *TPP1* variant interpretation.

---

## Materials and Methods

The pipeline is implemented in Snakemake [4]. Variants are normalized with BCFtools [5], annotated with VEP v115 [6] (incorporating REVEL [7] and BayesDel [8]), and splicing effects are predicted with SpliceAI [9]. Population frequencies are retrieved from gnomAD v4.1.0 [10] and clinical records from ClinVar [11].

ACMG criteria are mapped to points: PVS1 (+8), PS (+4), PM (+2), PP (+1), BP (-1), BM (-2), BS (-4), BA (-8). The final classification thresholds are: Pathogenic (>= 10), Likely Pathogenic (6-9), VUS (0-5), Likely Benign (-1 to -5), Benign (<= -6).

---

## Results and Discussion

We analyzed **337 real ClinVar variants** in the *TPP1* gene. The pipeline classified 37 as Pathogenic, 26 as Likely Pathogenic, 103 as VUS, 57 as Likely Benign, and 114 as Benign. A total of **63 variants (18.7%)** were classified as clinically actionable (Pathogenic or Likely Pathogenic).

Comparison with the commercial Franklin by Genoox platform showed **75.8% concordance**, with discordances primarily driven by conservative splicing predictions (SpliceAI) and in silico predictor (REVEL) thresholds in our pipeline.

---

## References

1. Markham, A. (2017). Cerliponase Alfa: First Global Approval. *Drugs*, 77(11), 1247-1253.
2. Richards, S., et al. (2015). Standards and guidelines for the interpretation of sequence variants. *Genetics in Medicine*, 17(5), 405-424.
3. Tavtigian, S. V., et al. (2020). Modeling the ACMG/AMP variant classification guidelines as a Bayesian classification framework. *Genetics in Medicine*, 22(3), 617-626.
4. Mölder, F., et al. (2021). Sustainable data analysis with Snakemake. *F1000Research*, 10, 490.
5. Danecek, P., et al. (2021). Twelve years of SAMtools and BCFtools. *GigaScience*, 10(2), giab008.
6. McLaren, W., et al. (2016). The Ensembl Variant Effect Predictor. *Genome Biology*, 17(1), 122.
7. Ioannidis, N. M., et al. (2016). REVEL: An Ensemble Method for Predicting the Pathogenicity of Rare Missense Variants. *American Journal of Human Genetics*, 99(4), 877-885.
8. Feng, B. J. (2017). Analysis of Functional Impact of Missense Variants Using BayesDel. *Human Mutation*, 38(10), 1234-1245.
9. Jaganathan, K., et al. (2019). Predicting Splicing from Primary Sequence with Deep Learning. *Cell*, 176(3), 535-548.
10. Karczewski, K. J., et al. (2020). The mutational constraint spectrum quantified from variation in 141,456 humans. *Nature*, 581(7809), 434-443.
11. Landrum, M. J., et al. (2018). ClinVar: improving access to variant interpretations and supporting evidence. *Nucleic Acids Research*, 46(D1), D1062-D1067.
