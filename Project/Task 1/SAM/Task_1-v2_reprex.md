# Task 1 - Differential Expression Analysis

Clean up the data and prepare for DE analysis with edgeR

    #;-) 
    #;-) Attaching package: 'dplyr'
    #;-) The following objects are masked from 'package:stats':
    #;-) 
    #;-)     filter, lag
    #;-) The following objects are masked from 'package:base':
    #;-) 
    #;-)     intersect, setdiff, setequal, union

    #;-) Genes: 60683 | Samples: 6
    #;-) # A tibble: 6 × 6
    #;-)   sample               subject treatment time  run    condition
    #;-)   <chr>                <fct>   <fct>     <chr> <chr>  <fct>    
    #;-) 1 S_2_Rem_5dpi_S70001  S_2     Rem       5dpi  S70001 Other    
    #;-) 2 S_2_SARS_5dpi_S70003 S_2     SARS      5dpi  S70003 COVID    
    #;-) 3 S_2_mock_5dpi_S70002 S_2     mock      5dpi  S70002 Control  
    #;-) 4 S_3_Rem_5dpi_S69995  S_3     Rem       5dpi  S69995 Other    
    #;-) 5 S_3_SARS_5dpi_S69996 S_3     SARS      5dpi  S69996 COVID    
    #;-) 6 S_3_mock_5dpi_S69997 S_3     mock      5dpi  S69997 Control

We start off by importing packages (tidyverse, limma, edgeR etc), and load in the raw data counts from the provided tsv file. Parsing those, we create the counts matrix and sample metadata as separate objects. That makes it easier to work with in later steps.

## Data

After loading the data we’re starting with ~60,000 genes, across the six samples.

The treatments/condition outline the different categories, with subject indicating the donor (experiment) each with replicates. Time and run are not used in the design, as they are confounded with treatment.

Ordinarily, we would filter out the samples belonging to the REM treatment, as our focus is on the COVID vs Control. However, due to the small sample size (n=2 per group), we will retain all samples for the DE analysis to maximize statistical power, and only report on the SARS vs mock contrast.

That way the REM samples can still contribute to controlling the subject variation, and accordingly make the COVID vs Control contrast more robust.

## Filtering Genes

    #;-) Kept 13375 genes after PC + filterByExpr
    #;-) [1] "(Intercept)"   "subjectS_3"    "treatmentSARS" "treatmentRem"

To reduce multiple-testing issues, and increase power, we drop uninformative rows and non-protein-coding genes. This is done by leveraging edgeR’s filterByExpr function, which keeps genes with enough counts given the design/groups.

The effect of this is dropping the gene count from ~60,000 to ~15,000, which is a more manageable number for DE analysis.

## Design Notes

We’ll take a moment here to consider this design, as its a bit unintuitive.

Essentially since we have such a low sample count, we want to ensure that whatever variance we find is due to the factor of interest, and not an artefact of individual variability. Under normal circumstances, we would have many more samples and many more replicates, which would help to average out individual variability.

Since we don’t have that luxury here, we instead try to account for it in the design. This is done by including the subject (donor) as a blocking factor, which helps to control for individual variability.

![](https://i.imgur.com/O69gDrP.png)<!-- -->

    #;-) [1] -0.01648809

Per-donor SARS–mock log2FCs hardly correlate (ρ ≈ −0.02), i.e., donors don’t move in lockstep. With only two donors, that heterogeneity plus BH explains why single-tool FDR finds no hits.

    #;-) # A tibble: 4 × 4
    #;-)   sample               subject treatment weight
    #;-)   <chr>                <fct>   <fct>      <dbl>
    #;-) 1 S_2_SARS_5dpi_S70003 S_2     SARS        3.48
    #;-) 2 S_2_mock_5dpi_S70002 S_2     mock        5.49
    #;-) 3 S_3_SARS_5dpi_S69996 S_3     SARS        5.49
    #;-) 4 S_3_mock_5dpi_S69997 S_3     mock        6.26

![](https://i.imgur.com/fnK1HCl.png)<!-- -->

One SARS sample has the lowest weight by far, indicating that its the noisiest. Voom downweights automatically, but with only 2 donors, it still limits power.

We’ll keep it for now, but we will have to keep it in mind when assessing the robustness of our findings.

## Limma voom

Limma Voom combines the voom transformation (which models the mean–variance relationship of RNA-seq counts) with the limma linear modelling and empirical Bayes framework. This approach converts count data to log-counts per million (logCPM) with precision weights, then fits fast, flexible linear models.

### Benefits

- Handles thousands of genes and potentially complex designs with ease.
- Simple to handle covariates, blocking or paired samples (which is useful in this instance)
- Bayes moderation improves the stability (but really needs the replicates to be equal to or greater than 4, see limitations)
- Contrasts and effect size thresholds are easy to implement.

### Limitations

- Difficulty handling small sample sizes (a significant issue here)
- Difficulty with sparse data sets (not as much an issue)
- Assumes a log normalised data, rather than modelling counts directly. This is surmountable.

### Differential Analysis

With only two samples (effectively) per group, limma will struggle on its own. It will get signal, but it is unlikely that it will be able to reach statistically robust conclusions without help.

### Results

    #;-) [1] 0

Per strict results, no significant genes were found at 5% FDR and 2x fold-change.

Given the requirements, and the donor heterogeneity, this is not entirely unexpected.

### Interpretation

The lack of significant genes under strict criteria suggests that the observed changes may not be robust, potentially due to the limited sample size and high variability.

Limma voom is generally a powerful method, but with only two donors per group and considerable variability, it struggles to identify significant DE genes under stringent thresholds.

One solution is to relax the criteria, for example by lowering the fold-change threshold to 1.5x (log2FC of 0.58). This can help to identify more genes that may be biologically relevant, albeit with a higher risk of false positives.

We can also use the second DE tool, edgeR, to provide an independent check on the findings. If both methods agree on certain genes, we can be more confident in their validity, even if they don’t meet the strictest criteria in either method alone.

## DESeq2

## EdgeR

EdgeR is another widely used tool for differential expression analysis for RNA data. It models the count data directly against a negative binomial distribution, a natural fit for RNA’s dispersed seq counts. It estimates the common, trended and tagwise distributions and only then fits the GLM (Generalised Linear Model) to test for the DE.

In this instance we also use the glmQLFit to implement treast-style testing with a minimum log-fold threshold.

### Benefits

- Negative Binomial works well with small sample sizes (good in this instance)
- Robust dispersion estimates can handle the outliers reliably.
- GLM’s accomodate complex contrasts, covariates and paired designs (useful in this case)
- Also allows for the treat-style tests, letting us check against the minimum meaningul log fold change.

### Limmitations

- It does rely on the NB being a good fit. If there are extreme outliers, or an otherwise distorted distribution, this will challenge it.
- While its possible to handle continous data, it requires a chunk of work. Given counts as our base data, we can somewhat handle this.
- Complexity compounds quickly, once we start venturing past the base into QL vs LRT or treat vs standard. Increases the risk of a poor conclusion.

### Differential Analysis

With our dataset only having the effective 2 samples per group, edgeR is generally more reliable, but will still have issues without more replicates. If the variance isn’t controlled confidently, it will miss true signals.

### Results

    #;-) [1] 4

EdgeR finds four significant genes at 5% FDR and 2x fold-change, which is a modest improvement over limma voom but still limited by the sample size.

### Interpretation

Even with robust NB modelling, the tiny n and high donor variability restrict statistical power.

Again, under normal circumstances with more samples, we would expect more robust findings. Given the current constraints, we can consider relaxing the criteria to identify more potential DE genes.

In the interim, we can combine the results from both methods. If both methods agree on certain genes, we can be more confident in their validity, even if they don’t meet the strictest criteria in either method alone.

## Consensus

Combining results from both methods to find robust DE genes.

    #;-) $strict_n
    #;-) [1] 0
    #;-) 
    #;-) $relaxed_n
    #;-) [1] 25
    #;-) [1] 0
    #;-) [1] 25
    #;-)               gene gene_name logFC_limma      p_limma FDR_limma logFC_edgeR
    #;-) 1  ENSG00000169245    CXCL10    7.674656 0.0019737513  0.834224    8.355664
    #;-) 2  ENSG00000169248    CXCL11    6.646544 0.0030143128  0.834224    6.903596
    #;-) 3  ENSG00000165949     IFI27    4.614557 0.0003199601  0.834224    4.290016
    #;-) 4  ENSG00000111335      OAS2    9.345509 0.0062945131  0.834224   10.969933
    #;-) 5  ENSG00000183486       MX2    4.362118 0.0025178178  0.834224    4.514855
    #;-) 6  ENSG00000185885    IFITM1    2.173705 0.0028657271  0.834224    2.232746
    #;-) 7  ENSG00000137959    IFI44L    4.054755 0.0029700210  0.834224    3.939224
    #;-) 8  ENSG00000130303      BST2    2.458443 0.0056908554  0.834224    2.593005
    #;-) 9  ENSG00000068079     IFI35    1.940231 0.0013875208  0.834224    1.931519
    #;-) 10 ENSG00000184979     USP18    3.664720 0.0040300184  0.834224    3.445659
    #;-)         p_edgeR    FDR_edgeR agree_dir mean_logFC     p_fisher   FDR_fisher
    #;-) 1  1.290094e-08 0.0001725501      TRUE   8.015160 6.466082e-10 2.744205e-06
    #;-) 2  1.273797e-07 0.0008518516      TRUE   6.775070 8.708446e-09 1.847932e-05
    #;-) 3  4.429500e-06 0.0148111403      TRUE   4.452287 3.029335e-08 4.285499e-05
    #;-) 4  2.408210e-06 0.0107366049      TRUE  10.157721 2.880830e-07 3.056561e-04
    #;-) 5  2.759295e-05 0.0738111536      TRUE   4.438487 1.214567e-06 1.030924e-03
    #;-) 6  2.579383e-04 0.4928463138      TRUE   2.203225 1.117473e-05 7.904258e-03
    #;-) 7  3.176810e-04 0.5311229373      TRUE   3.996990 1.403357e-05 7.983990e-03
    #;-) 8  1.787040e-04 0.3983609644      TRUE   2.525724 1.504993e-05 7.983990e-03
    #;-) 9  1.613510e-03 0.9999988525      TRUE   1.935875 3.136436e-05 1.331103e-02
    #;-) 10 5.185079e-04 0.6935043228      TRUE   3.555189 2.941847e-05 1.331103e-02
    #;-)    log2FC_S2 log2FC_S3     CRISPR_action
    #;-) 1  1.5457822  8.608265 CRISPRi (repress)
    #;-) 2  2.7125705  6.873306 CRISPRi (repress)
    #;-) 3  1.9079216  5.096290 CRISPRi (repress)
    #;-) 4  4.2205814 10.535523 CRISPRi (repress)
    #;-) 5  0.4493125  4.763338 CRISPRi (repress)
    #;-) 6  0.9790080  2.994869 CRISPRi (repress)
    #;-) 7  0.6678397  4.692395 CRISPRi (repress)
    #;-) 8  1.1945973  3.384587 CRISPRi (repress)
    #;-) 9  1.0517964  2.375402 CRISPRi (repress)
    #;-) 10 0.6596151  4.384116 CRISPRi (repress)

Under the strict consensus, no genes are found, which is not surprising given the earlier results.

So we take the approach of relaxing the criteria, and use each other as a check.

We require that both methods agree on the direction of change, and use a Fisher combined p-value with a more lenient FDR threshold of 10%. We also lower the fold-change requirement to 1.5x (log2FC of 0.58).

This ‘relaxed’ consensus finds 25 significant genes, which is a more reasonable number to work with, we can output these as potential targets for Task 2, and the CRISPR interventions (CRISPRa for down in COVID, CRISPRi for up). If we can find safe guides to targets these genes specifically, we can look at the impact of increasing and decreasing their expressions respectively.

### Check with Leave One Out (LOO)

To confirm the robustness of our findings, we can perform a Leave One Out analysis, where we repeat the DE analysis multiple times, each time leaving out one sample. This helps to ensure that our results are not overly dependent on any single sample.

It’s potentially overkill for this small dataset, but useful as a safeguard.

    #;-) [1] 25
    #;-)  [1] "ENSG00000169245" "ENSG00000169248" "ENSG00000165949" "ENSG00000111335"
    #;-)  [5] "ENSG00000183486" "ENSG00000185885" "ENSG00000137959" "ENSG00000130303"
    #;-)  [9] "ENSG00000068079" "ENSG00000184979" "ENSG00000197632" "ENSG00000111331"
    #;-) [13] "ENSG00000175445" "ENSG00000126709" "ENSG00000198932" "ENSG00000146433"
    #;-) [17] "ENSG00000187608" "ENSG00000172183" "ENSG00000266412" "ENSG00000134321"
    #;-) [21] "ENSG00000107036" "ENSG00000135913" "ENSG00000178385" "ENSG00000188283"
    #;-) [25] "ENSG00000185900"

LOO confirms that there are 25 robust hits in the relaxed set, which is a good sign. These genes are likely to be truly differentially expressed due to the treatment, rather than being artifacts of individual sample variability.

In essence this gives us some confidence in our findings, and helps to reassure that the genes noted earlier as noisier are not driving the results.

### BH Diagnostics

Finally we confirm the FDR control via Benjamini-Hochberg (BH) diagnostics. This is a useful sanity-check to ensure that our p-values and FDR adjustments are behaving as expected.

If our values are off base, then the robustness checks are moot.

BH asks “is the p-value small enough, given its rank among all p-values, to be called significant at the chosen FDR level q?”

The plots and summaries below confirm that while Limma-treats p values aren’t great, the Fisher-combined p-values are good enough for the first 5-6 ranks. The relaxed hits at 5-10% FDR then are reasonable, and we can trust them.

    #;-) # A tibble: 15 × 7
    #;-)     rank        x p_sorted cutoff_q05 pass_q05 cutoff_q10 pass_q10
    #;-)    <int>    <dbl>    <dbl>      <dbl> <lgl>         <dbl> <lgl>   
    #;-)  1     1 0.000236 0.000320  0.0000118 FALSE     0.0000236 FALSE   
    #;-)  2     2 0.000471 0.000389  0.0000236 FALSE     0.0000471 FALSE   
    #;-)  3     3 0.000707 0.000432  0.0000353 FALSE     0.0000707 FALSE   
    #;-)  4     4 0.000943 0.00105   0.0000471 FALSE     0.0000943 FALSE   
    #;-)  5     5 0.00118  0.00139   0.0000589 FALSE     0.000118  FALSE   
    #;-)  6     6 0.00141  0.00163   0.0000707 FALSE     0.000141  FALSE   
    #;-)  7     7 0.00165  0.00181   0.0000825 FALSE     0.000165  FALSE   
    #;-)  8     8 0.00189  0.00197   0.0000943 FALSE     0.000189  FALSE   
    #;-)  9     9 0.00212  0.00204   0.000106  FALSE     0.000212  FALSE   
    #;-) 10    10 0.00236  0.00227   0.000118  FALSE     0.000236  FALSE   
    #;-) 11    11 0.00259  0.00244   0.000130  FALSE     0.000259  FALSE   
    #;-) 12    12 0.00283  0.00252   0.000141  FALSE     0.000283  FALSE   
    #;-) 13    13 0.00306  0.00287   0.000153  FALSE     0.000306  FALSE   
    #;-) 14    14 0.00330  0.00297   0.000165  FALSE     0.000330  FALSE   
    #;-) 15    15 0.00353  0.00301   0.000177  FALSE     0.000353  FALSE
    #;-) q05 q10 
    #;-)   0   0
    #;-) adj<05 adj<10 
    #;-)      0      0

![](https://i.imgur.com/BT8dG7v.png)<!-- -->

    #;-) # A tibble: 15 × 7
    #;-)     rank        x p_sorted cutoff_q05 pass_q05 cutoff_q10 pass_q10
    #;-)    <int>    <dbl>    <dbl>      <dbl> <lgl>         <dbl> <lgl>   
    #;-)  1     1 0.000236 6.47e-10  0.0000118 TRUE      0.0000236 TRUE    
    #;-)  2     2 0.000471 8.71e- 9  0.0000236 TRUE      0.0000471 TRUE    
    #;-)  3     3 0.000707 3.03e- 8  0.0000353 TRUE      0.0000707 TRUE    
    #;-)  4     4 0.000943 2.88e- 7  0.0000471 TRUE      0.0000943 TRUE    
    #;-)  5     5 0.00118  1.21e- 6  0.0000589 TRUE      0.000118  TRUE    
    #;-)  6     6 0.00141  1.12e- 5  0.0000707 TRUE      0.000141  TRUE    
    #;-)  7     7 0.00165  1.40e- 5  0.0000825 TRUE      0.000165  TRUE    
    #;-)  8     8 0.00189  1.50e- 5  0.0000943 TRUE      0.000189  TRUE    
    #;-)  9     9 0.00212  2.94e- 5  0.000106  TRUE      0.000212  TRUE    
    #;-) 10    10 0.00236  3.14e- 5  0.000118  TRUE      0.000236  TRUE    
    #;-) 11    11 0.00259  3.83e- 5  0.000130  TRUE      0.000259  TRUE    
    #;-) 12    12 0.00283  5.47e- 5  0.000141  TRUE      0.000283  TRUE    
    #;-) 13    13 0.00306  9.74e- 5  0.000153  TRUE      0.000306  TRUE    
    #;-) 14    14 0.00330  1.05e- 4  0.000165  TRUE      0.000330  TRUE    
    #;-) 15    15 0.00353  1.16e- 4  0.000177  TRUE      0.000353  TRUE
    #;-) q05 q10 
    #;-)  15  25
    #;-) adj<05 adj<10 
    #;-)     15     25

![](https://i.imgur.com/qK2aZJf.png)<!-- -->

    #;-) BH critical ranks — limma-treat: q=0.05→0, q=0.10→0; discoveries via p.adjust: adj<05 0, adj<10 0
    #;-) BH critical ranks — Fisher: q=0.05→15, q=0.10→25; discoveries via p.adjust: adj<05 15, adj<10 25

## Export Results (for Re-Use in later tasks)

With the results confirmed, we can export the relevant tables for use in later tasks.
