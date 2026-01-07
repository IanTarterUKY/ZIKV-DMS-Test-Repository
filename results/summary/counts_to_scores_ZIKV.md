# Analyze counts and compute escape scores
This Python Jupyter notebook analyzes the variant counts and looks at mutation coverage and jackpotting.
It then computes an "escape scores" for each variant after grouping by barcode or substitutions as specified in the configuration.

## Set up analysis

This notebook primarily makes use of the Bloom lab's [dms_variants](https://jbloomlab.github.io/dms_variants) package, and uses [plotnine](https://github.com/has2k1/plotnine) for ggplot2-like plotting syntax:


```python
import collections
import math
import os
import warnings

import Bio.SeqIO

import dms_variants.codonvarianttable
from dms_variants.constants import CBPALETTE
import dms_variants.plotnine_themes

from IPython.display import display, HTML

import matplotlib.pyplot as plt

import numpy

import pandas as pd

from plotnine import *

import seaborn

import yaml

%matplotlib inline
```

Set [plotnine](https://github.com/has2k1/plotnine) theme to the gray-grid one defined in [dms_variants](https://jbloomlab.github.io/dms_variants):


```python
theme_set(dms_variants.plotnine_themes.theme_graygrid())
```

Versions of key software:


```python
print(f"Using dms_variants version {dms_variants.__version__}")
```

    Using dms_variants version 1.4.3


Ignore warnings that clutter output:


```python
warnings.simplefilter('ignore')
```

Read the configuration file:


```python
with open('config.yaml') as f:
    config = yaml.safe_load(f)
```

Create output directory:


```python
os.makedirs(config['escape_scores_dir'], exist_ok=True)
```

Read information about the samples:


```python
samples_df = pd.read_csv(config['barcode_runs_ZIKV'])
```

## Initialize codon-variant table
Initialize [CodonVariantTable](https://jbloomlab.github.io/dms_variants/dms_variants.codonvarianttable.html#dms_variants.codonvarianttable.CodonVariantTable) from wildtype gene sequence and the variant counts CSV file.
We will then use the plotting functions of this variant table to analyze the counts per sample:


```python
wt_seqrecord = Bio.SeqIO.read(config['wildtype_sequence_ZIKV'], 'fasta')
geneseq = str(wt_seqrecord.seq)
primary_target = wt_seqrecord.name
print(f"Read sequence of {len(geneseq)} nt for {primary_target} from {config['wildtype_sequence_ZIKV']}")
      
print(f"Initializing CodonVariantTable from gene sequence and {config['variant_counts_ZIKV']}")
      
variants = dms_variants.codonvarianttable.CodonVariantTable.from_variant_count_df(
                geneseq=geneseq,
                variant_count_df_file=config['variant_counts_ZIKV'],
                primary_target=primary_target)
      
print('Done initializing CodonVariantTable.')
```

    Read sequence of 327 nt for ZIKV from data/wildtype_zikv.fasta
    Initializing CodonVariantTable from gene sequence and results/counts/ZIKV/variant_counts.csv.gz
    Done initializing CodonVariantTable.


## Sequencing counts per sample
Average counts per variant for each sample.
Note that these are the **sequencing** counts, in some cases they may outstrip the actual number of sorted cells:


```python
p = variants.plotAvgCountsPerVariant(libraries=variants.libraries,
                                     by_target=False,
                                     orientation='v')
p = p + theme(panel_grid_major_x=element_blank())  # no vertical grid lines
_ = p.draw()
```


    
![png](counts_to_scores_ZIKV_files/counts_to_scores_ZIKV_19_0.png)
    


And the numerical values plotted above:


```python
display(HTML(
 variants.avgCountsPerVariant(libraries=variants.libraries,
                               by_target=False)
 .pivot_table(index='sample',
              columns='library',
              values='avg_counts_per_variant')
 .round(1)
 .to_html()
 ))
```


<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>library</th>
      <th>lib1</th>
    </tr>
    <tr>
      <th>sample</th>
      <th></th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>exp1-none-0-ref</th>
      <td>253.8</td>
    </tr>
    <tr>
      <th>exp1-Msc-1000-ref</th>
      <td>1061.7</td>
    </tr>
    <tr>
      <th>exp1-Z2131-1000-abneg</th>
      <td>279.5</td>
    </tr>
    <tr>
      <th>exp1-Z6161-1000-abneg</th>
      <td>191.3</td>
    </tr>
    <tr>
      <th>exp1-Z6163-1000-abneg</th>
      <td>151.1</td>
    </tr>
    <tr>
      <th>exp1-Z3039-1000-abneg</th>
      <td>308.9</td>
    </tr>
    <tr>
      <th>exp1-Z2029-1000-abneg</th>
      <td>190.8</td>
    </tr>
    <tr>
      <th>exp1-Z2039-1000-abneg</th>
      <td>138.8</td>
    </tr>
    <tr>
      <th>exp1-Z2044-1000-abneg</th>
      <td>200.5</td>
    </tr>
    <tr>
      <th>exp1-Z2103-1000-abneg</th>
      <td>198.8</td>
    </tr>
    <tr>
      <th>exp1-Z3051-1000-abneg</th>
      <td>183.2</td>
    </tr>
    <tr>
      <th>exp1-Z2212-1000-abneg</th>
      <td>49.3</td>
    </tr>
    <tr>
      <th>exp1-ZV2-1000-abneg</th>
      <td>183.7</td>
    </tr>
    <tr>
      <th>exp1-ZV16-1000-abneg</th>
      <td>155.8</td>
    </tr>
    <tr>
      <th>exp1-ZV67-1000-abneg</th>
      <td>124.9</td>
    </tr>
    <tr>
      <th>exp1-Z2102-1000-abneg</th>
      <td>235.6</td>
    </tr>
    <tr>
      <th>exp1-Z2208-1000-abneg</th>
      <td>138.9</td>
    </tr>
    <tr>
      <th>exp1-Z3022-1000-abneg</th>
      <td>143.8</td>
    </tr>
    <tr>
      <th>exp1-Z6147-1000-abneg</th>
      <td>123.5</td>
    </tr>
  </tbody>
</table>


## Mutations per variant
Average number of mutations per gene among all variants of the primary target, separately for each date:


```python
#this plotting is very slow when lots of samples, so for now plots are commented out

for date, date_df in samples_df.groupby('date', sort=False):
   p = variants.plotNumCodonMutsByType(variant_type='all',
                                       orientation='v',
                                       libraries=variants.libraries,
                                       samples=date_df['sample'].unique().tolist(),
                                       widthscale=2)
   p = p + theme(panel_grid_major_x=element_blank())  # no vertical grid lines
   fig = p.draw()
   display(fig)
   plt.close(fig)
```


    
![png](counts_to_scores_ZIKV_files/counts_to_scores_ZIKV_23_0.png)
    



    
![png](counts_to_scores_ZIKV_files/counts_to_scores_ZIKV_23_1.png)
    


Now similar plots but showing mutation frequency across the gene:


```python
# this plotting is very slow when lots of samples, so for now code commented out

for date, date_df in samples_df.groupby('date', sort=False):
   p = variants.plotMutFreqs(variant_type='all',
                             mut_type='codon',
                             orientation='v',
                             libraries=variants.libraries,
                             samples=date_df['sample'].unique().tolist(),
                             widthscale=1.5)
   fig = p.draw()
   display(fig)
   plt.close(fig)
```


    
![png](counts_to_scores_ZIKV_files/counts_to_scores_ZIKV_25_0.png)
    



    
![png](counts_to_scores_ZIKV_files/counts_to_scores_ZIKV_25_1.png)
    


## Jackpotting and mutation coverage in pre-selection libraries
We look at the distribution of counts in the "reference" (pre-selection) libraries to see if they seem jackpotted (a few variants at very high frequency):


```python
pre_samples_df = samples_df.query('selection == "reference"')
```

Distribution of mutations along the gene for the pre-selection samples; big spikes may indicate jackpotting:


```python
# this plotting is very slow when lots of samples, so for now code commented out

p = variants.plotMutFreqs(variant_type='all',
                         mut_type='codon',
                         orientation='v',
                         libraries=variants.libraries,
                         samples=pre_samples_df['sample'].unique().tolist(),
                         widthscale=1.5)
_ = p.draw()
```


    
![png](counts_to_scores_ZIKV_files/counts_to_scores_ZIKV_29_0.png)
    


How many mutations are observed frequently in pre-selection libraries?
Note that the libraries have been pre-selected for ZV-2/16/67 binding, so we expect stop variants to mostly be missing.
Make the plot both for all variants and just single-mutant variants:


```python
# this plotting is very slow when lots of samples, so for now code commented out

for variant_type in ['all', 'single']:
   p = variants.plotCumulMutCoverage(
                         variant_type=variant_type,
                         mut_type='aa',
                         orientation='v',
                         libraries=variants.libraries,
                         samples=pre_samples_df['sample'].unique().tolist(),
                         widthscale=1.8,
                         heightscale=1.2)
   _ = p.draw()
```


    
![png](counts_to_scores_ZIKV_files/counts_to_scores_ZIKV_31_0.png)
    



    
![png](counts_to_scores_ZIKV_files/counts_to_scores_ZIKV_31_1.png)
    


Now make a plot showing what number and fraction of counts are for each variant in each pre-selection sample / library.
If some variants constitute a very high fraction, then that indicates extensive bottlenecking:


```python
for ystat in ['frac_counts', 'count']:
    p = variants.plotCountsPerVariant(ystat=ystat,
                                      libraries=variants.libraries,
                                      samples=pre_samples_df['sample'].unique().tolist(),
                                      orientation='v',
                                      widthscale=1.75,
                                      )
    _ = p.draw()
```


    
![png](counts_to_scores_ZIKV_files/counts_to_scores_ZIKV_33_0.png)
    



    
![png](counts_to_scores_ZIKV_files/counts_to_scores_ZIKV_33_1.png)
    


Now make the same plot breaking down by variant class, which enables determination of which types of variants are at high and low frequencies.
For this plot (unlike one above not classified by category) we only show variants of the primary target (not the homologs), and also group synonymous with wildtype in order to reduce number of plotted categories to make more interpretable:


```python
for ystat in ['frac_counts', 'count']:
    p = variants.plotCountsPerVariant(ystat=ystat,
                                      libraries=variants.libraries,
                                      samples=pre_samples_df['sample'].unique().tolist(),
                                      orientation='v',
                                      widthscale=1.75,
                                      by_variant_class=True,
                                      classifyVariants_kwargs={'syn_as_wt': True},
                                      primary_target_only=True,
                                      )
    _ = p.draw()
```


    
![png](counts_to_scores_ZIKV_files/counts_to_scores_ZIKV_35_0.png)
    



    
![png](counts_to_scores_ZIKV_files/counts_to_scores_ZIKV_35_1.png)
    


We also directly look to see what is the variant in each reference library / sample with the highest fraction counts.
Knowing if the highest frequency variant is shared helps determine **where** in the experiment the jackpotting happened:


```python
frac_counts_per_variant = (
        variants.add_frac_counts(variants.variant_count_df)
        .query(f"sample in {pre_samples_df['sample'].unique().tolist()}")
        )

display(HTML(
    frac_counts_per_variant
    .sort_values('frac_counts', ascending=False)
    .groupby(['library', 'sample'])
    .head(n=1)
    .sort_values(['library', 'sample'])
    .set_index(['library', 'sample'])
    [['frac_counts', 'target', 'barcode', 'aa_substitutions', 'codon_substitutions']]
    .round(4)
    .to_html()
    ))
```


<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th></th>
      <th>frac_counts</th>
      <th>target</th>
      <th>barcode</th>
      <th>aa_substitutions</th>
      <th>codon_substitutions</th>
    </tr>
    <tr>
      <th>library</th>
      <th>sample</th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th rowspan="2" valign="top">lib1</th>
      <th>exp1-none-0-ref</th>
      <td>0.0834</td>
      <td>ZIKV</td>
      <td>GGGCATCGGCAAACCT</td>
      <td></td>
      <td></td>
    </tr>
    <tr>
      <th>exp1-Msc-1000-ref</th>
      <td>0.0898</td>
      <td>ZIKV</td>
      <td>GGGCATCGGCAAACCT</td>
      <td></td>
      <td></td>
    </tr>
  </tbody>
</table>


To further where the jackpotting relative to the generation of the reference samples, we plot the correlation among the fraction of counts for the different reference samples.
If the fractions are highly correlated, that indicates that the jackpotting occurred in some upstream step common to the reference samples:


```python
# this code makes a full matrix of scatter plots, but is REALLY SLOW. So for now,
# it is commented out in favor of code that just makes correlation matrix.
for lib, lib_df in frac_counts_per_variant.groupby('library'):
   wide_lib_df = lib_df.pivot_table(index=['target', 'barcode'],
                                    columns='sample',
                                    values='frac_counts')
   g = seaborn.pairplot(wide_lib_df, corner=True, plot_kws={'alpha': 0.5}, diag_kind='kde')
   _ = g.fig.suptitle(lib, size=18)
   plt.show()
```


    
![png](counts_to_scores_ZIKV_files/counts_to_scores_ZIKV_39_0.png)
    


## Examine counts for wildtype variants
The type of score we use to quantify escape depends on how well represented wildtype is in the selected libraries.
If wildtype is still well represented, we can use a more conventional functional score that gives differential selection relative to wildtype.
If wildtype is not well represented, then we need an alternative score that does not involve normalizing frequencies to wildtype.

First get average fraction of counts per variant for each variant class:


```python
counts_by_class = (
    variants.variant_count_df
    .pipe(variants.add_frac_counts)
    .pipe(variants.classifyVariants,
          primary_target=variants.primary_target,
          non_primary_target_class='homolog',
          class_as_categorical=True)
    .groupby(['library', 'sample', 'variant_class'])
    .aggregate(avg_frac_counts=pd.NamedAgg('frac_counts', 'mean'))
    .reset_index()
    .merge(samples_df[['sample', 'library', 'date', 'antibody', 'concentration', 'selection']],
           on=['sample', 'library'], validate='many_to_one')
    )

display(HTML(counts_by_class.head().to_html(index=False)))
```


<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>library</th>
      <th>sample</th>
      <th>variant_class</th>
      <th>avg_frac_counts</th>
      <th>date</th>
      <th>antibody</th>
      <th>concentration</th>
      <th>selection</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>lib1</td>
      <td>exp1-none-0-ref</td>
      <td>wildtype</td>
      <td>0.000003</td>
      <td>7/29/2025</td>
      <td>none</td>
      <td>0</td>
      <td>reference</td>
    </tr>
    <tr>
      <td>lib1</td>
      <td>exp1-none-0-ref</td>
      <td>synonymous</td>
      <td>0.000001</td>
      <td>7/29/2025</td>
      <td>none</td>
      <td>0</td>
      <td>reference</td>
    </tr>
    <tr>
      <td>lib1</td>
      <td>exp1-none-0-ref</td>
      <td>1 nonsynonymous</td>
      <td>0.000002</td>
      <td>7/29/2025</td>
      <td>none</td>
      <td>0</td>
      <td>reference</td>
    </tr>
    <tr>
      <td>lib1</td>
      <td>exp1-none-0-ref</td>
      <td>&gt;1 nonsynonymous</td>
      <td>0.000005</td>
      <td>7/29/2025</td>
      <td>none</td>
      <td>0</td>
      <td>reference</td>
    </tr>
    <tr>
      <td>lib1</td>
      <td>exp1-none-0-ref</td>
      <td>stop</td>
      <td>0.000001</td>
      <td>7/29/2025</td>
      <td>none</td>
      <td>0</td>
      <td>reference</td>
    </tr>
  </tbody>
</table>


Plot average fraction of all counts per variant for each variant class.
If the values for wildtype are low for the non-reference samples (such as more similar to stop the nonsynonymous), then normalizing by wildtype in calculating scores will probably not work well as wildtype is too depleted:


```python
min_frac = 1e-7  # plot values < this as this

p = (ggplot(counts_by_class
            .assign(avg_frac_counts=lambda x: numpy.clip(x['avg_frac_counts'], min_frac, None))
            ) +
     aes('avg_frac_counts', 'sample', color='selection') +
     geom_point(size=2) +
     scale_color_manual(values=CBPALETTE[1:]) +
     facet_grid('library ~ variant_class') +
     scale_x_log10() +
     theme(axis_text_x=element_text(angle=90),
           figure_size=(2.5 * counts_by_class['variant_class'].nunique(),
                        0.2 * counts_by_class['library'].nunique() * 
                        counts_by_class['sample'].nunique())
           ) +
     geom_vline(xintercept=min_frac, linetype='dotted', color=CBPALETTE[3])
     )

_ = p.draw()
```


    
![png](counts_to_scores_ZIKV_files/counts_to_scores_ZIKV_43_0.png)
    


## Compute escape scores
We use the escape score metric, which does **not** involve normalizing to wildtype and so isn't strongly affected by low wildtype counts.
We compute the scores using the method [dms_variants.codonvarianttable.CodonVariantTable.escape_scores](https://jbloomlab.github.io/dms_variants/dms_variants.codonvarianttable.html?highlight=escape_scores#dms_variants.codonvarianttable.CodonVariantTable.escape_scores).

First, define what samples to compare for each calculation, matching each post-selection (escape) to the pre-selection (reference) sample on the same date:


```python
score_sample_df = (
    samples_df
    .query('selection == "escape"')
    .rename(columns={'sample': 'post_sample',
                     'number_cells': 'pre_cells_sorted'})
    .merge(samples_df
           .query('selection == "reference"')
           [['sample', 'library', 'date', 'number_cells']]
           .rename(columns={'sample': 'pre_sample',
                            'number_cells': 'post_cells_sorted'}),
           how='left', on=['date', 'library'], validate='many_to_one',
           )
    .assign(name=lambda x: x['antibody'] + '_' + x['concentration'].astype(str))
    # add dates to names where needed to make unique
    .assign(n_libs=lambda x: x.groupby(['name', 'date'])['pre_sample'].transform('count'))
    .sort_values(['name', 'date', 'n_libs'], ascending=False)
    .assign(i_name=lambda x: x.groupby(['library', 'name'], sort=False)['pre_sample'].cumcount(),
            name=lambda x: x.apply(lambda r: r['name'] + '_' + str(r['date']) if r['i_name'] > 0 else r['name'],
                                   axis=1),
            )
    .sort_values(['antibody', 'concentration', 'library', 'i_name'])
    # get columns of interest
    [['name', 'library', 'antibody', 'concentration', 'date',
      'pre_sample', 'post_sample', 'frac_escape', 'pre_cells_sorted', 'post_cells_sorted']]
    )

assert len(score_sample_df.groupby(['name', 'library'])) == len(score_sample_df)

print(f"Writing samples used to compute escape scores to {config['escape_score_samples_ZIKV']}\n")
score_sample_df.to_csv(config['escape_score_samples_ZIKV'], index=False)

display(HTML(score_sample_df.to_html(index=False)))
```

    Writing samples used to compute escape scores to results/escape_scores/samples_ZIKV.csv
    



<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>name</th>
      <th>library</th>
      <th>antibody</th>
      <th>concentration</th>
      <th>date</th>
      <th>pre_sample</th>
      <th>post_sample</th>
      <th>frac_escape</th>
      <th>pre_cells_sorted</th>
      <th>post_cells_sorted</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>Z2029_1000</td>
      <td>lib1</td>
      <td>Z2029</td>
      <td>1000</td>
      <td>7/31/2025</td>
      <td>exp1-Msc-1000-ref</td>
      <td>exp1-Z2029-1000-abneg</td>
      <td>0.043530</td>
      <td>400000.0</td>
      <td>NaN</td>
    </tr>
    <tr>
      <td>Z2039_1000</td>
      <td>lib1</td>
      <td>Z2039</td>
      <td>1000</td>
      <td>7/31/2025</td>
      <td>exp1-Msc-1000-ref</td>
      <td>exp1-Z2039-1000-abneg</td>
      <td>0.046397</td>
      <td>394000.0</td>
      <td>NaN</td>
    </tr>
    <tr>
      <td>Z2044_1000</td>
      <td>lib1</td>
      <td>Z2044</td>
      <td>1000</td>
      <td>7/31/2025</td>
      <td>exp1-Msc-1000-ref</td>
      <td>exp1-Z2044-1000-abneg</td>
      <td>0.075850</td>
      <td>550000.0</td>
      <td>NaN</td>
    </tr>
    <tr>
      <td>Z2102_1000</td>
      <td>lib1</td>
      <td>Z2102</td>
      <td>1000</td>
      <td>7/31/2025</td>
      <td>exp1-Msc-1000-ref</td>
      <td>exp1-Z2102-1000-abneg</td>
      <td>0.013522</td>
      <td>99000.0</td>
      <td>NaN</td>
    </tr>
    <tr>
      <td>Z2103_1000</td>
      <td>lib1</td>
      <td>Z2103</td>
      <td>1000</td>
      <td>7/31/2025</td>
      <td>exp1-Msc-1000-ref</td>
      <td>exp1-Z2103-1000-abneg</td>
      <td>0.026778</td>
      <td>209000.0</td>
      <td>NaN</td>
    </tr>
    <tr>
      <td>Z2131_1000</td>
      <td>lib1</td>
      <td>Z2131</td>
      <td>1000</td>
      <td>7/31/2025</td>
      <td>exp1-Msc-1000-ref</td>
      <td>exp1-Z2131-1000-abneg</td>
      <td>0.048761</td>
      <td>475000.0</td>
      <td>NaN</td>
    </tr>
    <tr>
      <td>Z2208_1000</td>
      <td>lib1</td>
      <td>Z2208</td>
      <td>1000</td>
      <td>7/31/2025</td>
      <td>exp1-Msc-1000-ref</td>
      <td>exp1-Z2208-1000-abneg</td>
      <td>0.015506</td>
      <td>122000.0</td>
      <td>NaN</td>
    </tr>
    <tr>
      <td>Z2212_1000</td>
      <td>lib1</td>
      <td>Z2212</td>
      <td>1000</td>
      <td>7/31/2025</td>
      <td>exp1-Msc-1000-ref</td>
      <td>exp1-Z2212-1000-abneg</td>
      <td>0.003806</td>
      <td>27000.0</td>
      <td>NaN</td>
    </tr>
    <tr>
      <td>Z3022_1000</td>
      <td>lib1</td>
      <td>Z3022</td>
      <td>1000</td>
      <td>7/31/2025</td>
      <td>exp1-Msc-1000-ref</td>
      <td>exp1-Z3022-1000-abneg</td>
      <td>0.002809</td>
      <td>21800.0</td>
      <td>NaN</td>
    </tr>
    <tr>
      <td>Z3039_1000</td>
      <td>lib1</td>
      <td>Z3039</td>
      <td>1000</td>
      <td>7/31/2025</td>
      <td>exp1-Msc-1000-ref</td>
      <td>exp1-Z3039-1000-abneg</td>
      <td>0.053299</td>
      <td>475000.0</td>
      <td>NaN</td>
    </tr>
    <tr>
      <td>Z3051_1000</td>
      <td>lib1</td>
      <td>Z3051</td>
      <td>1000</td>
      <td>7/31/2025</td>
      <td>exp1-Msc-1000-ref</td>
      <td>exp1-Z3051-1000-abneg</td>
      <td>0.044180</td>
      <td>385000.0</td>
      <td>NaN</td>
    </tr>
    <tr>
      <td>Z6147_1000</td>
      <td>lib1</td>
      <td>Z6147</td>
      <td>1000</td>
      <td>7/31/2025</td>
      <td>exp1-Msc-1000-ref</td>
      <td>exp1-Z6147-1000-abneg</td>
      <td>0.037429</td>
      <td>297000.0</td>
      <td>NaN</td>
    </tr>
    <tr>
      <td>Z6161_1000</td>
      <td>lib1</td>
      <td>Z6161</td>
      <td>1000</td>
      <td>7/31/2025</td>
      <td>exp1-Msc-1000-ref</td>
      <td>exp1-Z6161-1000-abneg</td>
      <td>0.043332</td>
      <td>413000.0</td>
      <td>NaN</td>
    </tr>
    <tr>
      <td>Z6163_1000</td>
      <td>lib1</td>
      <td>Z6163</td>
      <td>1000</td>
      <td>7/31/2025</td>
      <td>exp1-Msc-1000-ref</td>
      <td>exp1-Z6163-1000-abneg</td>
      <td>0.035633</td>
      <td>360000.0</td>
      <td>NaN</td>
    </tr>
    <tr>
      <td>ZV16_1000</td>
      <td>lib1</td>
      <td>ZV16</td>
      <td>1000</td>
      <td>7/31/2025</td>
      <td>exp1-Msc-1000-ref</td>
      <td>exp1-ZV16-1000-abneg</td>
      <td>0.115498</td>
      <td>754000.0</td>
      <td>NaN</td>
    </tr>
    <tr>
      <td>ZV2_1000</td>
      <td>lib1</td>
      <td>ZV2</td>
      <td>1000</td>
      <td>7/31/2025</td>
      <td>exp1-Msc-1000-ref</td>
      <td>exp1-ZV2-1000-abneg</td>
      <td>0.107901</td>
      <td>869000.0</td>
      <td>NaN</td>
    </tr>
    <tr>
      <td>ZV67_1000</td>
      <td>lib1</td>
      <td>ZV67</td>
      <td>1000</td>
      <td>7/31/2025</td>
      <td>exp1-Msc-1000-ref</td>
      <td>exp1-ZV67-1000-abneg</td>
      <td>0.003442</td>
      <td>28000.0</td>
      <td>NaN</td>
    </tr>
  </tbody>
</table>


Compute the escape scores for variants of the primary target and classify the variants:


```python
print(f"Computing escape scores for {primary_target} variants using {config['escape_score_type']} "
      f"score type with a pseudocount of {config['escape_score_pseudocount']} and "
      f"an escape fraction floor {config['escape_score_floor_E']}, an escape fraction ceiling "
      f"{config['escape_score_ceil_E']}, and grouping variants by {config['escape_score_group_by']}.")

escape_scores = (variants.escape_scores(score_sample_df,
                                        score_type=config['escape_score_type'],
                                        pseudocount=config['escape_score_pseudocount'],
                                        floor_E=config['escape_score_floor_E'],
                                        ceil_E=config['escape_score_ceil_E'],
                                        by=config['escape_score_group_by'],
                                        )
                 .query('target == @primary_target')
                 .pipe(variants.classifyVariants,
                       primary_target=variants.primary_target,
                       syn_as_wt=(config['escape_score_group_by'] == 'aa_substitutions'),
                       )
                 )
print('Here are the first few lines of the resulting escape scores:')
display(HTML(escape_scores.head().to_html(index=False)))
```

    Computing escape scores for ZIKV variants using frac_escape score type with a pseudocount of 0.5 and an escape fraction floor 0, an escape fraction ceiling 1, and grouping variants by barcode.
    Here are the first few lines of the resulting escape scores:



<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>name</th>
      <th>target</th>
      <th>library</th>
      <th>pre_sample</th>
      <th>post_sample</th>
      <th>barcode</th>
      <th>score</th>
      <th>score_var</th>
      <th>pre_count</th>
      <th>post_count</th>
      <th>codon_substitutions</th>
      <th>n_codon_substitutions</th>
      <th>aa_substitutions</th>
      <th>n_aa_substitutions</th>
      <th>variant_class</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>Z2029_1000</td>
      <td>ZIKV</td>
      <td>lib1</td>
      <td>exp1-Msc-1000-ref</td>
      <td>exp1-Z2029-1000-abneg</td>
      <td>GGGCATCGGCAAACCT</td>
      <td>0.038507</td>
      <td>2.261641e-10</td>
      <td>40297687</td>
      <td>6421676</td>
      <td></td>
      <td>0</td>
      <td></td>
      <td>0</td>
      <td>wildtype</td>
    </tr>
    <tr>
      <td>Z2029_1000</td>
      <td>ZIKV</td>
      <td>lib1</td>
      <td>exp1-Msc-1000-ref</td>
      <td>exp1-Z2029-1000-abneg</td>
      <td>GGGGCACTACGCTTCT</td>
      <td>0.048129</td>
      <td>4.385351e-10</td>
      <td>27677288</td>
      <td>5512578</td>
      <td></td>
      <td>0</td>
      <td></td>
      <td>0</td>
      <td>wildtype</td>
    </tr>
    <tr>
      <td>Z2029_1000</td>
      <td>ZIKV</td>
      <td>lib1</td>
      <td>exp1-Msc-1000-ref</td>
      <td>exp1-Z2029-1000-abneg</td>
      <td>ATTCATGAGTTTCAGT</td>
      <td>0.045075</td>
      <td>8.539195e-10</td>
      <td>14166710</td>
      <td>2642613</td>
      <td></td>
      <td>0</td>
      <td></td>
      <td>0</td>
      <td>wildtype</td>
    </tr>
    <tr>
      <td>Z2029_1000</td>
      <td>ZIKV</td>
      <td>lib1</td>
      <td>exp1-Msc-1000-ref</td>
      <td>exp1-Z2029-1000-abneg</td>
      <td>TACTGCACCCGTTGTG</td>
      <td>0.049546</td>
      <td>9.755161e-10</td>
      <td>13793392</td>
      <td>2828173</td>
      <td></td>
      <td>0</td>
      <td></td>
      <td>0</td>
      <td>wildtype</td>
    </tr>
    <tr>
      <td>Z2029_1000</td>
      <td>ZIKV</td>
      <td>lib1</td>
      <td>exp1-Msc-1000-ref</td>
      <td>exp1-Z2029-1000-abneg</td>
      <td>TAGTGCCGGGGCTGAT</td>
      <td>0.052769</td>
      <td>1.145835e-09</td>
      <td>12674369</td>
      <td>2767797</td>
      <td></td>
      <td>0</td>
      <td></td>
      <td>0</td>
      <td>wildtype</td>
    </tr>
  </tbody>
</table>


## Apply pre-selection count filter to variant escape scores
Now determine a pre-selection count filter in order to flag for removal variants with counts that are so low that the estimated score is probably noise.
We know that stop codons should be largely purged pre-selection, and so the counts for them are a good indication of the "noise" threshold.
We therefore set the filter using the number of pre-selection counts for the stop codons.

To do this, we first compute the number of pre-selection counts for stop-codon variants at various quantiles and look at these.
We then take the number of pre-selection counts at the specified quantile as the filter cutoff, and filter scores for all variants with pre-selection counts less than this filter cutoff:

NOTE, since we didn't pre-sort this library, this is not the appropriate method. This has historically beena round 100, so we'll just apply a filter of 100.


```python
filter_precount = config['escape_score_precount_filter']

print(f"\nSetting the pre-count filter cutoff to the absolute value of {config['escape_score_precount_filter']}")
```

    
    Setting the pre-count filter cutoff to the absolute value of 100


Apply the filter to the escape scores, so that scores that fail the pre-selection count filter are now marked with `pass_pre_count_filter` of `False`:


```python
escape_scores = (
    escape_scores
    .assign(pass_pre_count_filter=lambda x: x['pre_count'] >= config['escape_score_precount_filter'])
    )
escape_scores.head()
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>name</th>
      <th>target</th>
      <th>library</th>
      <th>pre_sample</th>
      <th>post_sample</th>
      <th>barcode</th>
      <th>score</th>
      <th>score_var</th>
      <th>pre_count</th>
      <th>post_count</th>
      <th>codon_substitutions</th>
      <th>n_codon_substitutions</th>
      <th>aa_substitutions</th>
      <th>n_aa_substitutions</th>
      <th>variant_class</th>
      <th>pass_pre_count_filter</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>Z2029_1000</td>
      <td>ZIKV</td>
      <td>lib1</td>
      <td>exp1-Msc-1000-ref</td>
      <td>exp1-Z2029-1000-abneg</td>
      <td>GGGCATCGGCAAACCT</td>
      <td>0.038507</td>
      <td>2.261641e-10</td>
      <td>40297687</td>
      <td>6421676</td>
      <td></td>
      <td>0</td>
      <td></td>
      <td>0</td>
      <td>wildtype</td>
      <td>True</td>
    </tr>
    <tr>
      <th>1</th>
      <td>Z2029_1000</td>
      <td>ZIKV</td>
      <td>lib1</td>
      <td>exp1-Msc-1000-ref</td>
      <td>exp1-Z2029-1000-abneg</td>
      <td>GGGGCACTACGCTTCT</td>
      <td>0.048129</td>
      <td>4.385351e-10</td>
      <td>27677288</td>
      <td>5512578</td>
      <td></td>
      <td>0</td>
      <td></td>
      <td>0</td>
      <td>wildtype</td>
      <td>True</td>
    </tr>
    <tr>
      <th>2</th>
      <td>Z2029_1000</td>
      <td>ZIKV</td>
      <td>lib1</td>
      <td>exp1-Msc-1000-ref</td>
      <td>exp1-Z2029-1000-abneg</td>
      <td>ATTCATGAGTTTCAGT</td>
      <td>0.045075</td>
      <td>8.539195e-10</td>
      <td>14166710</td>
      <td>2642613</td>
      <td></td>
      <td>0</td>
      <td></td>
      <td>0</td>
      <td>wildtype</td>
      <td>True</td>
    </tr>
    <tr>
      <th>3</th>
      <td>Z2029_1000</td>
      <td>ZIKV</td>
      <td>lib1</td>
      <td>exp1-Msc-1000-ref</td>
      <td>exp1-Z2029-1000-abneg</td>
      <td>TACTGCACCCGTTGTG</td>
      <td>0.049546</td>
      <td>9.755161e-10</td>
      <td>13793392</td>
      <td>2828173</td>
      <td></td>
      <td>0</td>
      <td></td>
      <td>0</td>
      <td>wildtype</td>
      <td>True</td>
    </tr>
    <tr>
      <th>4</th>
      <td>Z2029_1000</td>
      <td>ZIKV</td>
      <td>lib1</td>
      <td>exp1-Msc-1000-ref</td>
      <td>exp1-Z2029-1000-abneg</td>
      <td>TAGTGCCGGGGCTGAT</td>
      <td>0.052769</td>
      <td>1.145835e-09</td>
      <td>12674369</td>
      <td>2767797</td>
      <td></td>
      <td>0</td>
      <td></td>
      <td>0</td>
      <td>wildtype</td>
      <td>True</td>
    </tr>
  </tbody>
</table>
</div>



## Apply ACE2-binding / expression filter to variant mutations
We also used deep mutational scanning to estimate how each mutation affected ACE2 binding and expression in the B.1.351 background.
Here we flag for removal any variants of the primary target that have (or have mutations) that were measured to decrease ACE2-binding or expression beyond a minimal threshold, in order to avoid these variants muddying the signal as spurious escape mutants.

To do this, we first determine all mutations that do / do-not having binding that exceeds the thresholds.

Note that because we are working on this serum-mapping project at the same time as we are working on the ACE2-binding / RBD-expression project, the scores will be preliminary until all final analyses have been done on the DMS project end. So, we will allow either preliminary or "final" measurements to be used. 


```python
#mut_bind_expr_file = config['mut_bind_expr']
    
#print(f"Reading ACE2-binding and expression for mutations from {mut_bind_expr_file}, "
#      f"and filtering for variants that have single mutations that "
#      f"only have mutations with binding >={config['escape_score_min_bind_mut_Omicron_BA2']} and "
#      f"expression >={config['escape_score_min_expr_mut_Omicron_BA2']}.")

#mut_bind_expr = (pd.read_csv(mut_bind_expr_file)
#                 .query("target == 'Omicron_BA2'")
#                 # need to add back the offset numbering for some silly, circuitous reason 
#                 .assign(RBD_site=lambda x: x['position']-config['site_number_offset'] ,
#                         RBD_mutation=lambda x: x['wildtype']+x['RBD_site'].astype(str)+x['mutant']
#                        )
#                )

#print('Here is what that dataframe looks like:')

#display(HTML(mut_bind_expr.query('delta_bind > -2.5').head().to_html(index=False)))
```


```python
#assert mut_bind_expr['RBD_mutation'].nunique() == len(mut_bind_expr)
#for prop in ['bind', 'expr']:
#    muts_adequate = set(mut_bind_expr
#                        .query(f"delta_{prop} >= {config[f'escape_score_min_{prop}_mut_Omicron_BA2']}")
#                        ['RBD_mutation']
#                        )
#    print(f"{len(muts_adequate)} of {len(mut_bind_expr)} mutations have adequate {prop}.")
#    escape_scores[f"muts_pass_{prop}_filter"] = (
#        escape_scores
#        ['aa_substitutions']
#        .map(lambda s: set(s.split()).issubset(muts_adequate))
#        )

# annotate as passing overall filter if passes all mutation and binding filters:
#escape_scores['pass_ACE2bind_expr_filter'] = (
#        escape_scores['muts_pass_bind_filter'] &
#        escape_scores['muts_pass_expr_filter'] 
#        )

#display(HTML(escape_scores.query('not pass_ACE2bind_expr_filter & variant_class != "stop"').head().to_html(index=False)))
```

Print the number of mutations that pass RBD bind, RBD expression, and are not to sites that are disulfide bonds (if specified in config) 


```python
# if we are excluding all cysteines to remove spurious mutations that break disulfide bonds:
#if config['exclude_cysteines']:
#    print("Here are the number of mutations that pass the bind, express, and disulfide filters:")
#    print(len(mut_bind_expr
#              .assign(pass_cysteine_filter=lambda x: x['mutation'].str[0] != "C")
#              .query(f"delta_bind >= {config[f'escape_score_min_bind_mut_Wuhan_Hu_1']} & \
#                       delta_expr >= {config[f'escape_score_min_expr_mut_Wuhan_Hu_1']} & \
#                       pass_cysteine_filter")
#             ))
#   print("There are these many possible mutations (excluding wildtype and disulfides!):")
#    print(mut_bind_expr.query('wildtype!="C"')['position'].nunique()*19
#         )

#else:
#    print("Here are the number of mutations that pass the bind and express filters:")
#    print(len(mut_bind_expr
#              .assign(pass_cysteine_filter=lambda x: x['mutation'].str[0] != "C")
#              .query(f"delta_bind >= {config[f'escape_score_min_bind_mut_Wuhan_Hu_1']} & \
#                       delta_expr >= {config[f'escape_score_min_expr_mut_Wuhan_Hu_1']}")
#             ))
#    print("There are these many possible mutations (excluding wildtype!):")
#    print(mut_bind_expr['position'].nunique()*19
#         )
```


```python
#print('These are the sites that are involved in disulfide bonds:')
#print(mut_bind_expr.query('wildtype=="C"')['position'].unique())
```


```python
#frac_ACE2bind_expr_pass_filter = (
#    escape_scores
#    .query('pass_pre_count_filter == True')
#    [['pre_sample', 'library', 'target', config['escape_score_group_by'],
#      'pre_count', 'pass_ACE2bind_expr_filter', 'variant_class']]
#    .drop_duplicates()
#    .groupby(['pre_sample', 'library', 'variant_class'], observed=True)
#    .aggregate(n_variants=pd.NamedAgg('pass_ACE2bind_expr_filter', 'count'),
#               n_pass_filter=pd.NamedAgg('pass_ACE2bind_expr_filter', 'sum')
#               )
#    .reset_index()
#    .assign(frac_pass_filter=lambda x: x['n_pass_filter'] / x['n_variants'],
#            pre_sample=lambda x: pd.Categorical(x['pre_sample'], x['pre_sample'].unique(), ordered=True).remove_unused_categories())
#    )

#p = (ggplot(frac_ACE2bind_expr_pass_filter) +
#     aes('variant_class', 'frac_pass_filter', fill='variant_class') +
#     geom_bar(stat='identity') +
#     facet_grid('library ~ pre_sample') +
#     theme(axis_text_x=element_text(angle=90),
#           figure_size=(3.3 * frac_ACE2bind_expr_pass_filter['pre_sample'].nunique(),
#                        2 * frac_ACE2bind_expr_pass_filter['library'].nunique()),
#           panel_grid_major_x=element_blank(),
#           ) +
#     scale_fill_manual(values=CBPALETTE[1:]) +
#     expand_limits(y=(0, 1))
#     )

#_ = p.draw()
```

## Examine and write escape scores
Plot the distribution of escape scores across variants of different classes **among those that pass both the pre-selection count filter and the ACE2-binding / expression filter**.
If things are working correctly, we don't expect escape in wildtype (or synonymous variants), but do expect escape for some small fraction of nonsynymous variants.
Also, we do not plot the scores for the stop codon variant class, as most stop-codon variants have already been filtered out so this category is largely noise:


```python
nfacets = len(escape_scores.groupby(['library', 'name']).nunique())
ncol = min(8, nfacets)
nrow = math.ceil(nfacets / ncol)

df = (escape_scores
      .query('(pass_pre_count_filter == True)')
      .query('variant_class != "stop"')
      )
     
p = (ggplot(df) +
     aes('variant_class', 'score', color='variant_class') +
     geom_boxplot(outlier_size=1.5, outlier_alpha=0.1) +
     facet_wrap('~ name + library', ncol=ncol) +
     theme(axis_text_x=element_text(angle=90),
           figure_size=(2.35 * ncol, 3 * nrow),
           panel_grid_major_x=element_blank(),
           ) +
     scale_fill_manual(values=CBPALETTE[1:]) +
     scale_color_manual(values=CBPALETTE[1:])
     )

_ = p.draw()
```


    
![png](counts_to_scores_ZIKV_files/counts_to_scores_ZIKV_60_0.png)
    


Also, we want to see how much the high escape scores are correlated with simple coverage.
To do this, we plot the correlation between escape score and pre-selection count just for the nonsynonymous variants (which are the ones that we expect to have true escape).
The plots below have a lot of overplotting, but are still sufficient to test of the score is simply correlated with the pre-selection counts or not.
The hoped for result is that the escape score doesn't appear to be strongly correlated with pre-selection counts:


```python
p = (ggplot(escape_scores
            .query('pass_pre_count_filter == True')
            .query('(pass_pre_count_filter == True)')
            .query('variant_class=="1 nonsynonymous"')
            ) +
     aes('pre_count', 'score') +
     geom_point(alpha=0.1, size=1) +
     facet_wrap('~ name + library', ncol=ncol) +
     theme(axis_text_x=element_text(angle=90),
           figure_size=(2.35 * ncol, 2.35 * nrow),
           ) +
     scale_fill_manual(values=CBPALETTE[1:]) +
     scale_color_manual(values=CBPALETTE[1:]) +
     scale_x_log10()
     )

_ = p.draw()
```


    
![png](counts_to_scores_ZIKV_files/counts_to_scores_ZIKV_62_0.png)
    


Write the escape scores to a file:


```python
print(f"Writing escape scores for ZIKV SSM to {config['escape_scores_ZIKV']}")
escape_scores.to_csv(config['escape_scores_ZIKV'], index=False, float_format='%.4g')
```

    Writing escape scores for ZIKV SSM to results/escape_scores/scores_ZIKV.csv


### Now we will also remove anything that did not pass all the filters above. 


```python
escape_scores_primary = (escape_scores
                         .query('(pass_pre_count_filter == True)')
                        )

display(HTML(escape_scores_primary.head().to_html()))
print(f"Read {len(escape_scores_primary)} scores.")
```


<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>name</th>
      <th>target</th>
      <th>library</th>
      <th>pre_sample</th>
      <th>post_sample</th>
      <th>barcode</th>
      <th>score</th>
      <th>score_var</th>
      <th>pre_count</th>
      <th>post_count</th>
      <th>codon_substitutions</th>
      <th>n_codon_substitutions</th>
      <th>aa_substitutions</th>
      <th>n_aa_substitutions</th>
      <th>variant_class</th>
      <th>pass_pre_count_filter</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>Z2029_1000</td>
      <td>ZIKV</td>
      <td>lib1</td>
      <td>exp1-Msc-1000-ref</td>
      <td>exp1-Z2029-1000-abneg</td>
      <td>GGGCATCGGCAAACCT</td>
      <td>0.038507</td>
      <td>2.261641e-10</td>
      <td>40297687</td>
      <td>6421676</td>
      <td></td>
      <td>0</td>
      <td></td>
      <td>0</td>
      <td>wildtype</td>
      <td>True</td>
    </tr>
    <tr>
      <th>1</th>
      <td>Z2029_1000</td>
      <td>ZIKV</td>
      <td>lib1</td>
      <td>exp1-Msc-1000-ref</td>
      <td>exp1-Z2029-1000-abneg</td>
      <td>GGGGCACTACGCTTCT</td>
      <td>0.048129</td>
      <td>4.385351e-10</td>
      <td>27677288</td>
      <td>5512578</td>
      <td></td>
      <td>0</td>
      <td></td>
      <td>0</td>
      <td>wildtype</td>
      <td>True</td>
    </tr>
    <tr>
      <th>2</th>
      <td>Z2029_1000</td>
      <td>ZIKV</td>
      <td>lib1</td>
      <td>exp1-Msc-1000-ref</td>
      <td>exp1-Z2029-1000-abneg</td>
      <td>ATTCATGAGTTTCAGT</td>
      <td>0.045075</td>
      <td>8.539195e-10</td>
      <td>14166710</td>
      <td>2642613</td>
      <td></td>
      <td>0</td>
      <td></td>
      <td>0</td>
      <td>wildtype</td>
      <td>True</td>
    </tr>
    <tr>
      <th>3</th>
      <td>Z2029_1000</td>
      <td>ZIKV</td>
      <td>lib1</td>
      <td>exp1-Msc-1000-ref</td>
      <td>exp1-Z2029-1000-abneg</td>
      <td>TACTGCACCCGTTGTG</td>
      <td>0.049546</td>
      <td>9.755161e-10</td>
      <td>13793392</td>
      <td>2828173</td>
      <td></td>
      <td>0</td>
      <td></td>
      <td>0</td>
      <td>wildtype</td>
      <td>True</td>
    </tr>
    <tr>
      <th>4</th>
      <td>Z2029_1000</td>
      <td>ZIKV</td>
      <td>lib1</td>
      <td>exp1-Msc-1000-ref</td>
      <td>exp1-Z2029-1000-abneg</td>
      <td>TAGTGCCGGGGCTGAT</td>
      <td>0.052769</td>
      <td>1.145835e-09</td>
      <td>12674369</td>
      <td>2767797</td>
      <td></td>
      <td>0</td>
      <td></td>
      <td>0</td>
      <td>wildtype</td>
      <td>True</td>
    </tr>
  </tbody>
</table>


    Read 9435 scores.


### Count number of barcodes per mutation and remove variants with >1 amino acid substitution
Also add the number of barocdes per mutation to the `escape_scores` dataframe and plot this. 
But first see how many variants there are with >1 mutation, and query the dataframe to look at them qualitatively. 


```python
p = (ggplot(escape_scores_primary) +
     aes('n_aa_substitutions', fill='variant_class') +
     geom_bar() +
     facet_wrap('~library + pre_sample', ncol=5) +
     theme(
#            figure_size=(3.3 * escape_scores_primary['pre_sample'].nunique(),
#                         2 * escape_scores_primary['library'].nunique()),
         figure_size=(12, 4),
           panel_grid_major_x=element_blank(),
           ) +
     scale_fill_manual(values=CBPALETTE[1:]) +
     expand_limits(y=(0, 1))
     )

_ = p.draw()
```


    
![png](counts_to_scores_ZIKV_files/counts_to_scores_ZIKV_68_0.png)
    


### Filter dataframe on single mutations that are present in at least `n` number of variants (specified in `config.yaml` file)
Now see how many `n_single_mut_measurements` there are for each variant:


```python
print(f'Remove anything with fewer than {config["escape_frac_min_single_mut_measurements"]} single mutant measurements (barcodes)')

raw_avg_single_mut_scores = (
    escape_scores_primary
    .query('n_aa_substitutions == 1')
    .rename(columns={'name': 'selection',
                     'aa_substitutions': 'mutation'})
    .groupby(['selection', 'library', 'mutation'])
    .aggregate(raw_single_mut_score=pd.NamedAgg('score', 'mean'),
               n_single_mut_measurements=pd.NamedAgg('barcode', 'count')
              )
    .assign(sufficient_measurements=lambda x: (
                (x['n_single_mut_measurements'] >= config['escape_frac_min_single_mut_measurements'])))
    .reset_index()
    )

# remove mutations with insufficient measurements
effects_df = (raw_avg_single_mut_scores
              .query('sufficient_measurements == True')
              .drop(columns='sufficient_measurements')
              )

# some error checks
assert len(effects_df) == len(effects_df.drop_duplicates()), 'duplicate rows in `effects_df`'
assert all(effects_df['raw_single_mut_score'].notnull() | (effects_df['n_single_mut_measurements'] == 0))

display(HTML(effects_df.head().to_html()))
```

    Remove anything with fewer than 1 single mutant measurements (barcodes)



<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>selection</th>
      <th>library</th>
      <th>mutation</th>
      <th>raw_single_mut_score</th>
      <th>n_single_mut_measurements</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>Z2029_1000</td>
      <td>lib1</td>
      <td>A12V</td>
      <td>0.007665</td>
      <td>1</td>
    </tr>
    <tr>
      <th>1</th>
      <td>Z2029_1000</td>
      <td>lib1</td>
      <td>A13L</td>
      <td>0.020803</td>
      <td>1</td>
    </tr>
    <tr>
      <th>2</th>
      <td>Z2029_1000</td>
      <td>lib1</td>
      <td>A13S</td>
      <td>0.056631</td>
      <td>1</td>
    </tr>
    <tr>
      <th>3</th>
      <td>Z2029_1000</td>
      <td>lib1</td>
      <td>A21E</td>
      <td>0.000443</td>
      <td>1</td>
    </tr>
    <tr>
      <th>4</th>
      <td>Z2029_1000</td>
      <td>lib1</td>
      <td>A21N</td>
      <td>0.000564</td>
      <td>1</td>
    </tr>
  </tbody>
</table>


### Should we exclude mutations for which the wildtype identity is a cysteine?
These would be mutations that break disulfide bonds (unless there is an unpaired cysteine in the protein that does not form a disulfide bond, of course). 

If we are excluding cysteines, do that now:


```python
# if we are excluding all cysteines to remove spurious mutations that break disulfide bonds:
if config['exclude_cysteines']:
    print(f'Excluding mutations where the wildtype identity is a cysteine')
    effects_df = effects_df.assign(pass_cysteine_filter=lambda x: x['mutation'].str[0] != "C",
                                  )
    disulfides_to_drop = effects_df.query('pass_cysteine_filter==False')['mutation'].unique()
    print(f'Specifically, excluding: {disulfides_to_drop}')
    effects_df=effects_df.query('pass_cysteine_filter').drop(columns='pass_cysteine_filter')

else:
    print(f'Retaining mutations where the wildtype identity is a cysteine')
```

    Excluding mutations where the wildtype identity is a cysteine
    Specifically, excluding: ['C10S' 'C41G' 'C41R']


We need to compute the escape scores (calculated as [here](https://jbloomlab.github.io/dms_variants/dms_variants.codonvarianttable.html?highlight=escape_scores#dms_variants.codonvarianttable.CodonVariantTable.escape_scores)) back to escape fractions. We define a function to do this depending on the score type:


```python
def score_to_frac(score):
    """Convert escape score to escape fraction."""
    if pd.isnull(score):
        return pd.NA
    floor = config['escape_score_floor_E']
    ceil = config['escape_score_ceil_E']
    if config['escape_score_type'] == 'frac_escape':
        return min(ceil, max(floor, score))  # just the score after applying ceiling and floor
    elif config['escape_score_type'] == 'log_escape':
        # convert score back to fraction, apply ceiling, then re-scale so floor is 0
        frac = 2**score
        frac = min(ceil, max(floor, frac))
        frac = (frac - floor) / (ceil - floor)
        return frac
    else:
        raise ValueError(f"invalid `escape_score_type` of {config['escape_score_type']}")

effects_df = (
    effects_df
    .assign(
            mut_escape_frac_single_mut=lambda x: x['raw_single_mut_score'].map(score_to_frac),
            )
    )
```

### Average the escape score across all barcodes of the same mutation, for each library, and the average of both libraries. 
Add rows that are the average of the two libraries for the fraction escape for all mutations that are present in both libraries (and if in just one library, the value in that or purge depending on config values printed here), the number of libraries each mutation is measured in, and the sum of the statistics giving the number of measurements:


```python
min_libs = config['escape_frac_avg_min_libraries']
min_single = config['escape_frac_avg_min_single']
print(f"Only taking average of mutations with escape fractions in >={min_libs} libraries "
      f"or with >={min_single} single-mutant measurements total.")

effects_df = (
    effects_df
    .query('library != "average"')  # git rid of averages if already there
    .assign(nlibs=1)
    .append(effects_df
            .query('library != "average"')
            .groupby(['selection', 'mutation'])
            .aggregate(nlibs=pd.NamedAgg('library', 'count'),
                       mut_escape_frac_single_mut=pd.NamedAgg('mut_escape_frac_single_mut',
                                                              lambda s: s.mean(skipna=True)),
                       n_single_mut_measurements=pd.NamedAgg('n_single_mut_measurements', 'sum'),
                       )
            .query('(nlibs >= @min_libs) or (n_single_mut_measurements >= @min_single)')
            .assign(library="average")
            .reset_index(),
            ignore_index=True, sort=False,
            )
    )

print(len(effects_df.query('nlibs>1')))
print(len(effects_df.query('nlibs==1')))
```

    Only taking average of mutations with escape fractions in >=1 libraries or with >=2 single-mutant measurements total.
    0
    8500


Plot the correlations of the escape fractions among the two libraries for all selections performed on both libraries. 


```python
#libraries = [lib for lib in effects_df['library'].unique() if lib != "average"]
#assert len(libraries) == 2, 'plot only makes sense if 2 libraries'

# wide data frame with each library's score in a different column
#effects_df_wide = (
#    effects_df
#    .query('library != "average"')
#    .query(f"n_single_mut_measurements >= 1")
    # just get selections with 2 libraries
#    .assign(nlibs=lambda x: x.groupby('selection')['library'].transform('nunique'))
#    .query('nlibs == 2')
    # now make columns for each library, only keep mutants with scores for both libs
 #   [['selection', 'mutation', 'library', 'mut_escape_frac_single_mut']]
 #   .pivot_table(index=['selection', 'mutation'],
#                 columns='library',
#                 values='mut_escape_frac_single_mut',
#                 aggfunc='first')
#    .reset_index()
#    .dropna(axis=0)
#    )

# correlations between libraries
#corrs = (
#    effects_df_wide
#    .groupby('selection')
#    [libraries]
#    .corr(method='pearson')
#    .reset_index()
#    .query('library == @libraries[0]')
#    .assign(correlation=lambda x: 'R=' + x[libraries[1]].round(2).astype(str))
#    [['selection', 'correlation']]
    # add number of mutations measured
#    .merge(effects_df_wide
#           .groupby('selection')
#           .size()
#           .rename('n')
#           .reset_index()
#           )
#    .assign(correlation=lambda x: x['correlation'] + ', N=' + x['n'].astype(str))
#    )

# plot correlations
#nfacets = effects_df_wide['selection'].nunique()
#ncol = min(nfacets, 5)
#nrow = math.ceil(nfacets / ncol)
#xmin = effects_df_wide[libraries[0]].min()
#xspan = effects_df_wide[libraries[0]].max() - xmin
#ymin = effects_df_wide[libraries[1]].min()
#yspan = effects_df_wide[libraries[1]].max() - ymin
#p = (ggplot(effects_df_wide) +
#     aes(libraries[0], libraries[1]) +
#     geom_point(alpha=0.2) +
#     geom_text(mapping=aes(label='correlation'),
#               data=corrs,
#               x=0.01 * xspan + xmin,
#               y=0.99 * yspan + ymin,
#               size=10,
#               ha='left',
#               va='top',
#               ) +
#     facet_wrap('~ selection', ncol=ncol) +
#     theme(figure_size=(2.5 * ncol, 2.5 * nrow),
#           plot_title=element_text(size=14)) +
#     ggtitle('Mutation-level escape fractions')
#     )

#_ = p.draw()
```

### Escape at site level
The above analysis estimates the effects of mutations. We also compute escape statistics at the site level. First, add sites to the data frame of mutational effects:


```python
effects_df = (
    effects_df
    .assign(site=lambda x: x['mutation'].str[1: -1].astype(int),
            wildtype=lambda x: x['mutation'].str[0],
            mutant=lambda x: x['mutation'].str[-1],
            )
    )
```

Now compute some site-level metrics. These are the average and total escape fraction at each site over all mutations at the site:


```python
site_effects_df = (
    effects_df
    .groupby(['selection', 'library', 'site'])
    .aggregate(
        site_avg_escape_frac_single_mut=pd.NamedAgg('mut_escape_frac_single_mut',
                                                    lambda s: s.mean(skipna=True)),
        site_total_escape_frac_single_mut=pd.NamedAgg('mut_escape_frac_single_mut', 'sum'),
        )
    .reset_index()
    )

display(HTML(site_effects_df.head().to_html(index=False)))
```


<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>selection</th>
      <th>library</th>
      <th>site</th>
      <th>site_avg_escape_frac_single_mut</th>
      <th>site_total_escape_frac_single_mut</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>Z2029_1000</td>
      <td>average</td>
      <td>1</td>
      <td>0.028309</td>
      <td>0.113236</td>
    </tr>
    <tr>
      <td>Z2029_1000</td>
      <td>average</td>
      <td>2</td>
      <td>0.038776</td>
      <td>0.116328</td>
    </tr>
    <tr>
      <td>Z2029_1000</td>
      <td>average</td>
      <td>4</td>
      <td>0.049676</td>
      <td>0.298056</td>
    </tr>
    <tr>
      <td>Z2029_1000</td>
      <td>average</td>
      <td>5</td>
      <td>0.049653</td>
      <td>0.049653</td>
    </tr>
    <tr>
      <td>Z2029_1000</td>
      <td>average</td>
      <td>6</td>
      <td>0.068192</td>
      <td>0.204575</td>
    </tr>
  </tbody>
</table>


Plot correlations between libraries of the same selection for the site-level statistics:


```python
#libraries = [lib for lib in effects_df['library'].unique() if lib != "average"]
#assert len(libraries) == 2, 'plot only makes sense if 2 libraries'

#for val in ['site_avg_escape_frac_single_mut', 'site_total_escape_frac_single_mut']:

    # wide data frame with each library's score in a different column
#    site_effects_df_wide = (
#        site_effects_df
#        .query('library != "average"')
        # just get selections with 2 libraries
#        .assign(nlibs=lambda x: x.groupby('selection')['library'].transform('nunique'))
#        .query('nlibs == 2')
        # now make columns for each library, only keep sites with scores for both libs
#        .pivot_table(index=['selection', 'site'],
#                     columns='library',
#                     values=val)
#        .reset_index()
#        .dropna(axis=0)
#        )

    # correlations between libraries
#    corrs = (
 #       site_effects_df_wide
 #       .groupby('selection')
 #       [libraries]
 #       .corr(method='pearson')
 #       .reset_index()
 #       .query('library == @libraries[0]')
 #       .assign(correlation=lambda x: 'R=' + x[libraries[1]].round(2).astype(str))
 #       [['selection', 'correlation']]
        # add number of mutations measured
#       .merge(site_effects_df_wide
#               .groupby('selection')
#               .size()
#               .rename('n')
#               .reset_index()
#               )
#        .assign(correlation=lambda x: x['correlation'] + ', N=' + x['n'].astype(str))
#        )

    # plot correlations
#    nfacets = site_effects_df_wide['selection'].nunique()
#    ncol = min(nfacets, 5)
#    nrow = math.ceil(nfacets / ncol)
#    xmin = site_effects_df_wide[libraries[0]].min()
#    xspan = site_effects_df_wide[libraries[0]].max() - xmin
#    ymin = site_effects_df_wide[libraries[1]].min()
#    yspan = site_effects_df_wide[libraries[1]].max() - ymin
#    p = (ggplot(site_effects_df_wide) +
#         aes(libraries[0], libraries[1]) +
#         geom_point(alpha=0.2) +
#         geom_text(mapping=aes(label='correlation'),
#                   data=corrs,
#                   x=0.01 * xspan + xmin,
#                   y=0.99 * yspan + ymin,
#                   size=10,
#                   ha='left',
#                   va='top',
#                   ) +
#         facet_wrap('~ selection', ncol=ncol) +
#         theme(figure_size=(2.5 * ncol, 2.5 * nrow),
#               plot_title=element_text(size=14)) +
#         ggtitle(val)
#         )

#    _ = p.draw()
```

## Write file with escape fractions at mutation and site levels
We write a files that has the mutation- and site-level escape fractions. This file has both the separate measurements for each library plus the average across libraries for all mutations measured in both libraries. We name the columns in such a way that this file can be used as [dms-view data file](https://dms-view.github.io/docs/dataupload):


```python
escape_fracs_to_write = (
    effects_df
    .merge(site_effects_df,
           how='left',
           validate='many_to_one',
           on=['selection', 'library', 'site'])
    .assign(protein_chain=config['escape_frac_protein_chain'],
            protein_site=lambda x: x['site'] + config['site_number_offset'],
            label_site=lambda x: x['protein_site'],
            condition=lambda x: x['selection'].where(x['library'] == "average", x['selection'] + '_' + x['library']),
            mutation=lambda x: x['mutant'],  # dms-view uses mutation to refer to mutant amino acid
            )
    [['selection', 'library', 'condition', 'site', 'label_site', 'wildtype', 'mutation',
      'protein_chain', 'protein_site', 'mut_escape_frac_single_mut', 'site_total_escape_frac_single_mut',
      'site_avg_escape_frac_single_mut', 'nlibs', 'n_single_mut_measurements',
      ]]
    .sort_values(['library', 'selection', 'site', 'mutation'])
    )

print('Here are the first few lines that will be written to the escape-fraction file:')
display(HTML(escape_fracs_to_write.head().to_html(index=False)))

print(f"\nWriting to {config['escape_fracs_ZIKV']}")
escape_fracs_to_write.to_csv(config['escape_fracs_ZIKV'], index=False, float_format='%.4g')

```

    Here are the first few lines that will be written to the escape-fraction file:



<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>selection</th>
      <th>library</th>
      <th>condition</th>
      <th>site</th>
      <th>label_site</th>
      <th>wildtype</th>
      <th>mutation</th>
      <th>protein_chain</th>
      <th>protein_site</th>
      <th>mut_escape_frac_single_mut</th>
      <th>site_total_escape_frac_single_mut</th>
      <th>site_avg_escape_frac_single_mut</th>
      <th>nlibs</th>
      <th>n_single_mut_measurements</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>Z2029_1000</td>
      <td>average</td>
      <td>Z2029_1000</td>
      <td>1</td>
      <td>299</td>
      <td>R</td>
      <td>C</td>
      <td>E</td>
      <td>299</td>
      <td>0.045825</td>
      <td>0.113236</td>
      <td>0.028309</td>
      <td>1</td>
      <td>1</td>
    </tr>
    <tr>
      <td>Z2029_1000</td>
      <td>average</td>
      <td>Z2029_1000</td>
      <td>1</td>
      <td>299</td>
      <td>R</td>
      <td>K</td>
      <td>E</td>
      <td>299</td>
      <td>0.016533</td>
      <td>0.113236</td>
      <td>0.028309</td>
      <td>1</td>
      <td>1</td>
    </tr>
    <tr>
      <td>Z2029_1000</td>
      <td>average</td>
      <td>Z2029_1000</td>
      <td>1</td>
      <td>299</td>
      <td>R</td>
      <td>L</td>
      <td>E</td>
      <td>299</td>
      <td>0.050612</td>
      <td>0.113236</td>
      <td>0.028309</td>
      <td>1</td>
      <td>1</td>
    </tr>
    <tr>
      <td>Z2029_1000</td>
      <td>average</td>
      <td>Z2029_1000</td>
      <td>1</td>
      <td>299</td>
      <td>R</td>
      <td>N</td>
      <td>E</td>
      <td>299</td>
      <td>0.000266</td>
      <td>0.113236</td>
      <td>0.028309</td>
      <td>1</td>
      <td>1</td>
    </tr>
    <tr>
      <td>Z2029_1000</td>
      <td>average</td>
      <td>Z2029_1000</td>
      <td>2</td>
      <td>300</td>
      <td>L</td>
      <td>C</td>
      <td>E</td>
      <td>300</td>
      <td>0.028572</td>
      <td>0.116328</td>
      <td>0.038776</td>
      <td>1</td>
      <td>1</td>
    </tr>
  </tbody>
</table>


    
    Writing to results/escape_scores/escape_fracs_ZIKV.csv



```python

```
