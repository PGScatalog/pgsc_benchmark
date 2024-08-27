# `pgsc_benchmark` User Manual
R implementation. Version 1.0 (Aug 2024)

For the Docker implementation please see `User-manual-Docker.md`.

## Software Requirements

This script was developed using R version 4.3.3. You can download the latest version of
R from _The R Project for Statistical Computing_ [website](<https://www.r-project.org>).

Additional software is required to run the PGS Catalog Calculator (e.g. Nextflow, Docker).
Please see the [online calculator documentation](<https://pgsc-calc.readthedocs.io/en/latest/>)
for details.

### Dependencies

This script requires the following R packages:

- boot
- data.table
- dplyr
- ggplot2
- parallel
- pROC
- r2redux
- R.utils
- survival

These packages will be automatically installed upon running the script for the first
time. If you would prefer to install these manually, please do so before attempting
to run the script.

Depending on your system setup and file permissions, you may have to specify the path
to the desired library using the `.libPaths()` function. A template for this is provided
on line 20 of the script (remember to uncomment this line).

## Running the PGS Catalog Calculator

The `pgsc_benchmark` tool is designed to run directly on the output files produced
by the PGS Catalog Calculator (`pgsc_calc`). A helpful step-by-step guide for running
the calculator can be found [here](<https://pgsc-calc.readthedocs.io/en/latest/>).

Since the benchmarking tool uses genetic ancestry data in all analyses, the calculator
must be run using a reference dataset. We recommend using the provided reference panel
comprising merged data from the Human Genome Diversity Project and 1000 Genomes Project
(HGDP+1kGP). Instructions for downloading and using this dataset can be found in the
[online documentation](<https://pgsc-calc.readthedocs.io/en/latest/getting-started.html>)
for the calculator.

The easiest and recommended method for including all relevant PGS IDs in a single run is by
specifying a trait experimental factor ontology (EFO) term using the <code>&#8209;&#8209;trait_efo</code>
option. This will automatically include all PGS IDs that are associated with the specified trait
and all child traits. A list of recommended EFO terms for several common traits is provided below.
In each case, the EFO term (or HP/MONDO term, where appropriate) has been chosen to maximise
the number of potentially relevant PGS IDs included in the analysis. Where multiple terms are
indicated for a single trait, specify all terms (as a comma-separated list) to ensure the
maximum number of relevant PGS IDs are included.

### Recommended EFO terms for common traits

#### Binary traits

<table>
  <thead>
    <tr>
      <th>Trait name</th>
      <th>Recommended EFO term <em>(trait)</em></th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th colspan=2><em>Cardiovascular diseases</em></th>
    </tr>
    <tr>
      <td>Atrial fibrillation</td>
      <td>EFO_0000275 <em>(atrial fibrillation)</em></td>
    </tr>
    <tr>
      <td>Coronary artery disease</td>
      <td>EFO_0001645 <em>(coronary artery disease)</em></td>
    </tr>
    <tr>
      <td>Hypertension</td>
      <td>EFO_0000537 <em>(hypertension)</em></td>
    </tr>
    <tr>
      <td>Ischemic stroke</td>
      <td>EFO_0000712 <em>(stroke)</em></td>
    </tr>
    <tr>
      <td>Stroke</td>
      <td>EFO_0000712 <em>(stroke)</em> &<br>EFO_0005669 <em>(intracerebral hemorrhage)</em></td>
    </tr>
    <tr>
      <td>Venous thromboembolism</td>
      <td>EFO_0004286 <em>(venous thromboembolism)</em></td>
    </tr>
    <tr>
      <th colspan=2><em>Cancers</em></th>
    </tr>
    <tr>
      <td>Basal cell carcinoma</td>
      <td>MONDO_0002898 <em>(skin cancer)</em></td>
    </tr>
    <tr>
      <td>Bladder cancer</td>
      <td>MONDO_0001187 <em>(urinary bladder cancer)</em></td>
    </tr>
    <tr>
      <td>Breast cancer</td>
      <td>EFO_0000305 <em>(breast carcinoma)</em></td>
    </tr>
    <tr>
      <td>Colorectal cancer</td>
      <td>MONDO_0005575 <em>(colorectal cancer)</em></td>
    </tr>
    <tr>
      <td>Lung cancer</td>
      <td>MONDO_0008903 <em>(lung cancer)</em></td>
    </tr>
    <tr>
      <td>Melanoma</td>
      <td>EFO_0000756 <em>(melanoma)</em> &<br>MONDO_0002898 <em>(skin cancer)</em></td>
    </tr>
    <tr>
      <td>Ovarian cancer</td>
      <td>EFO_0003893 <em>(ovarian neoplasm)</em></td>
    </tr>
    <tr>
      <td>Pancreatic cancer</td>
      <td>EFO_0003860 <em>(pancreatic neoplasm)</em></td>
    </tr>
    <tr>
      <td>Prostate cancer</td>
      <td>MONDO_0008315 <em>(prostate cancer)</em></td>
    </tr>
    <tr>
      <td>Testicular cancer</td>
      <td>EFO_0005088 <em>(testicular carcinoma)</em></td>
    </tr>
    <tr>
      <td>Thyroid cancer</td>
      <td>EFO_0002892 <em>(thyroid carcinoma)</em></td>
    </tr>
    <tr>
      <th colspan=2><em>Neurological diseases</em></th>
    </tr>
    <tr>
      <td>Alzheimer’s disease</td>
      <td>MONDO_0004975 <em>(Alzheimer disease)</em></td>
    </tr>
    <tr>
      <td>Major depressive disorder</td>
      <td>MONDO_0002050 <em>(depressive disorder)</em></td>
    </tr>
    <tr>
      <td>Parkinson’s disease</td>
      <td>MONDO_0005180 <em>(Parkinson disease)</em></td>
    </tr>
    <tr>
      <td>Schizophrenia</td>
      <td>MONDO_0005090 <em>(schizophrenia)</em></td>
    </tr>
    <tr>
      <th colspan=2><em>Autoimmune diseases</em></th>
    </tr>
    <tr>
      <td>Celiac disease</td>
      <td>EFO_0001060 <em>(celiac disease)</em></td>
    </tr>
    <tr>
      <td>Inflammatory bowel disease</td>
      <td>EFO_0003767 <em>(inflammatory bowel disease)</em></td>
    </tr>
    <tr>
      <td>Lupus erythematosus</td>
      <td>MONDO_0004670 <em>(lupus erythematosus)</em></td>
    </tr>
    <tr>
      <td>Rheumatoid arthritis</td>
      <td>EFO_0000685 <em>(rheumatoid arthritis)</em></td>
    </tr>
    <tr>
      <th colspan=2><em>Diabetes</em></th>
    </tr>
    <tr>
      <td>Type 1 diabetes</td>
      <td>MONDO_0005147 <em>(type 1 diabetes mellitus)</em></td>
    </tr>
    <tr>
      <td>Type 2 diabetes</td>
      <td>MONDO_0005148 <em>(type 2 diabetes mellitus)</em></td>
    </tr>
    <tr>
      <th colspan=2><em>Other diseases</em></th>
    </tr>
    <tr>
      <td>Asthma</td>
      <td>MONDO_0004979 <em>(asthma)</em></td>
    </tr>
    <tr>
      <td>Chronic kidney disease</td>
      <td>EFO_0003884 <em>(chronic kidney disease)</em></td>
    </tr>
    <tr>
      <td>Glaucoma</td>
      <td>MONDO_0005041 <em>(glaucoma)</em></td>
    </tr>
    <tr>
      <td>Gout</td>
      <td>EFO_0004274 <em>(gout)</em></td>
    </tr>
    <tr>
      <td>Hypercholesterolemia</td>
      <td>HP_0003124 <em>(hypercholesterolemia)</em></td>
    </tr>
    <tr>
      <td>Hypothyroidism</td>
      <td>EFO_0004705 <em>(hypothyroidism)</em></td>
    </tr>
    <tr>
      <td>Non-alcoholic fatty liver disease</td>
      <td>EFO_0003095 <em>(non-alcoholic fatty liver disease)</em></td>
    </tr>
  </tbody>
</table>

#### Quantitative traits

<table>
  <thead>
    <tr>
      <th>Trait name</th>
      <th>Recommended EFO term <em>(trait)</em></th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th colspan=2><em>Body measurements</em></th>
    </tr>
    <tr>
      <td>Body mass index</td>
      <td>EFO_0004340 <em>(body mass index)</em></td>
    </tr>
    <tr>
      <td>Height</td>
      <td>EFO_0004339 <em>(body height)</em>
    </tr>
    <tr>
      <th colspan=2><em>Cardiovascular measurements</em></th>
    </tr>
    <tr>
      <td>Diastolic blood pressure</td>
      <td>EFO_0006336 <em>(diastolic blood pressure)</em></td>
    </tr>
    <tr>
      <td>Resting heart rate</td>
      <td>EFO_0004326 <em>(heart rate)</em></td>
    </tr>
    <tr>
      <td>Systolic blood pressure</td>
      <td>EFO_0006335 <em>(systolic blood pressure)</em></td>
    </tr>
    <tr>
      <th colspan=2><em>Lipids</em></th>
    </tr>
    <tr>
      <td>HDL cholesterol</td>
      <td>EFO_0004612 <em>(HDL cholesterol measurement)</em></td>
    </tr>
    <tr>
      <td>LDL cholesterol</td>
      <td>EFO_0004611 <em>(LDL cholesterol measurement)</em></td>
    </tr>
    <tr>
      <td>Total cholesterol</td>
      <td>EFO_0004574 <em>(total cholesterol measurement)</em></td>
    </tr>
    <tr>
      <td>Triglycerides</td>
      <td>EFO_0004530 <em>(triglyceride measurement)</em></td>
    </tr>
    <tr>
      <th colspan=2><em>Other biomarkers</em></th>
    </tr>
    <tr>
      <td>C-reactive protein</td>
      <td>EFO_0004458 <em>(C-reactive protein measurement)</em></td>
    </tr>
    <tr>
      <td>Fasting glucose</td>
      <td>EFO_0004468 <em>(glucose measurement)</em></td>
    </tr>
    <tr>
      <td>Fasting insulin</td>
      <td>EFO_0004466 <em>(fasting blood insulin measurement)</em></td>
    </tr>
    <tr>
      <td>Glomerular filtration rate</td>
      <td>EFO_0005208 <em>(glomerular filtration rate)</em></td>
    </tr>
    <tr>
      <td>HbA1c</td>
      <td>EFO_0004541 <em>(HbA1c measurement)</em></td>
    </tr>
    <tr>
      <td>Lipoprotein A</td>
      <td>EFO_0006925 <em>(lipoprotein A measurement)</em></td>
    </tr>
    <tr>
      <td>Testosterone</td>
      <td>EFO_0004908 <em>(testosterone measurement)</em></td>
    </tr>
    <tr>
      <td>Vitamin D</td>
      <td>EFO_0004631 <em>(vitamin D measurement)</em></td>
    </tr>
  </tbody>
</table>

<table><tr><td><strong>Note:</strong> the files produced by the PGS Catalog Calculator
  can be very large, sometimes containing tens of millions of rows of data. Keep this
  in mind when deciding how much memory to allocate per core when running the benchmarking
  script, particularly for phenotypes with larger numbers of scores. E.g. we required
  ~3,500 MiB of RAM per core to analyse 70 scores in 500,000 individuals.</td></tr></table>

<table><tr><td><strong>Note:</strong> scoring files for some traits (e.g. Alzheimer’s
  disease or autoimmune diseases) may contain complex alleles (e.g. APOE or HLA alleles)
  that are not currently supported by the PGS Catalog Calculator. However, the calculator
  will still produce an output for these PGS (simply excluding such alleles from the final
  score), so care should be taken when interpreting these results. We recommend ignoring
  results from all scores including these complex allele types, as these alleles commonly
  have large effect sizes that would normally drive the overall score.</td></tr></table>

## Preparing your input demographic/phenotype file

In addition to the score and population similarity files generated by the PGS Catalog
Calculator, you must also provide a file containing the demographic and phenotype data for
each benchmarking trait with rows for each participant in your cohort. This will usually be
a CSV or TSV file with a single header row, followed by the data for a single participant on
each subsequent row. The tool is designed to be flexible with how this input data is
provided and only a few specific requirements must be met.

All demographic files must contain three columns for the identification number, age and sex
of each participant. The following default column names are recommended, however other names
can be specified using the <code>&#8209;&#8209;id</code>, <code>&#8209;&#8209;age</code> or
<code>&#8209;&#8209;sex</code> options when running the tool:

<table>
  <thead>
    <tr>
      <th>Column&nbsp;name</th>
      <th>Description</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>IID</td>
      <td>Identification numbers in this column must match those in the genotype files used by
        the PGS Catalog Calculator</td>
    </tr>
    <tr>
      <td>age</td>
      <td>Age at recruitment (in years)</td>
    </tr>
    <tr>
      <td>sex</td>
      <td>Coded as <em>Male/Female</em>, <em>M/F</em>, or <em>1/0</em>. Casing is not important. <strong>Note:</strong>
        this column is still needed even if you only have single-sex data</td>
    </tr>
  </tbody>
</table>

<table><tr><td><strong>Note:</strong> column names should begin with an alphabetic character
  and should not contain whitespace characters.</td></tr></table>

#### Binary traits

For a binary trait, the disease status of each participant must be coded as <em>1</em> (case) or
<em>0</em> (non-case). Up to three separate columns may be specified, representing overall disease
status (including both prevalent and incident cases) as well as prevalent and incident cases
separately. The column names used must be specified using the <code>&#8209;&#8209;any&#8209;case</code>,
<code>&#8209;&#8209;prev&#8209;case</code> and <code>&#8209;&#8209;inc&#8209;case</code> options,
respectively. Prevalent cases should be excluded from the incident case column, either using <em>NA</em>
or just leaving the entry blank.

A column containing the follow-up times (in years) for the incident cohort may be specified
using the <code>&#8209;&#8209;follow&#8209;time</code> option. This should represent the time to
event or censoring from recruitment and will be used as the underlying timescale for Cox proportional
hazards modelling.

#### Quantitative traits

For a quantitative trait, a column containing the trait measurement for each participant must
be provided. The column name containing this data must be specified using the
<code>&#8209;&#8209;trait&#8209;val</code> option. Values should represent the observed unscaled
measurements of a trait, using the standard units for that trait (which should be specified using
the <code>&#8209;&#8209;units</code> option).

#### Other columns

All regression analyses must include cohort defined covariates (age and sex) and will include
10 genetic principal components (extracted from the `pgsc_calc` population similarity file).
If there are additional variables that should be included as covariates, these can be specified
as a list of column names using the <code>&#8209;&#8209;covs&#8209;cat</code> option for categorical
variables, or the <code>&#8209;&#8209;covs&#8209;quant</code> option for quantitative variables.

<table><tr><td><strong>Note:</strong> if your input demographic file contains any other unessential
  columns, these will simply be ignored and do not need to be deleted.</td></tr></table>

## Running the PGS benchmarking tool

### Running the script from the command line

The script is designed to be run from the command line. Input data, relevant column names and
other parameters/run specifications must be specified as options. Most options take only a single
argument and can be specified as follows (no whitespace):

`--option=arg`

Some options (`--covs-cat`, `--covs-quant`, `--ancestry` and `--score-type`) can take more than
one argument. Use commas to separate multiple arguments (no whitespace):

`--option=arg1,arg2,arg3`

A full list of available options is provided in <em>Options</em> below.

#### Binary trait example

The following example code shows how to run the benchmarking script for a binary trait
(stroke) using participants in the UK Biobank stratified by ancestry. Separate columns
for any stroke, prevalent stroke and incident stroke are provided.

```
Rscript pgsc_benchmark.R \
   --sampleset=UKB \
    --phenotype=Stroke \
    --outdir=~/projects/results \
    --dem-file=~/UKB/phenotypes/stroke_data.csv \
    --pgs-file=~/pgs/stroke/UKB/score/UKB_pgs.txt.gz \
    --pop-file=~/pgs/stroke/UKB/score/UKB_popsimilarity.txt.gz \
    --any-case=stroke \
    --prev-case=prevalent_stroke \
    --inc-case=incident_stroke \
    --follow-time=incident_follow_up \
    --stratified=TRUE
```

#### Quantitative trait example

The following example code shows how to run the benchmarking script for a quantitative
trait (body mass index) using the subset of participants in the UK Biobank most
genetically similar to the labelled East Asian samples in the reference dataset used
by the calculator.

```
Rscript pgsc_benchmark.R \
    --sampleset=UKB \
    --phenotype=BMI \
    --outdir=~/projects/results \
    --dem-file=~/UKB/phenotypes/bmi_data.csv \
    --pgs-file=~/pgs/bmi/UKB/score/UKB_pgs.txt.gz \
    --pop-file=~/pgs/bmi/UKB/score/UKB_popsimilarity.txt.gz \
    --trait-val=bmi \
    --units=kg/m2 \
    --ancestry=EAS
```

### Running the script interactively

If you need to run the benchmarking script interactively (e.g. using RStudio) you
can specify options as variables within the script. All variables can be found
within the _‘Data, Parameters and Running Requirements’_ section starting from line
57. Remember to enclose character strings within quotation marks and use vector
notation to specify multiple arguments for a single option.

<table><tr><td><strong>Note:</strong> editing the R script directly will change
  the default value of each option, so take caution if you later try running the
  script from the command line.</td></tr></table>

### Options

<table>
  <thead>
    <tr>
      <th>Option</th>
      <th>Default</th>
      <th>Description</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th colspan=3><em>Data description</em></th>
    </tr>
    <tr>
      <td>&#8209;&#8209;sampleset</td>
      <td></td>
      <td><strong>Required.</strong> Specify the name of the cohort/study. Spaces
        are allowed. This text will be used in output tables/figures. E.g.
        <br>--sampleset=UKB</td>
    </tr>
    <tr>
      <td>&#8209;&#8209;phenotype</td>
      <td></td>
      <td><strong>Required.</strong> Specify the phenotype being examined. Spaces
        are allowed. This text will be used in output tables/figures. E.g.
        <br>--phenotype="LDL cholesterol"</td>
    </tr>
    <tr>
      <th colspan=3><em>Input and output</em></th>
    </tr>
    <tr>
      <td>&#8209;&#8209;dem&#8209;file</td>
      <td></td>
      <td><strong>Required.</strong> Specify the path to the file containing
        demographic/phenotype data for each participant. E.g.
        <br>--dem-filename=~/ukb/phenotypes/bmi.csv</td>
    </tr>
    <tr>
      <td>&#8209;&#8209;pgs&#8209;file</td>
      <td></td>
      <td><strong>Required.</strong> Specify the path to the PGS data generated
        by the PGS Catalog Calculator (compressed or uncompressed). E.g.
        <br>--pgs-filename=~/.../score/ukb_pgs.txt.gz</td>
    </tr>
    <tr>
      <td>&#8209;&#8209;pop&#8209;file</td>
      <td></td>
      <td><strong>Required.</strong> Specify the path to the population similarity
        data generated by the PGS Catalog Calculator (compressed or uncompressed).
        E.g.<br>--pop-filename=~/.../score/ukb_popsimilarity.txt.gz</td>
    </tr>
    <tr>
      <td>&#8209;&#8209;outdir</td>
      <td>NA</td>
      <td>Specify a path to where the results directory should be created. If
        <em>NA</em> (default), this will be set to the current working directory.
        You should use the same directory across all runs, cohorts, phenotypes,
        ancestries, etc. E.g.<br>--outdir=~/projects</td>
      </tr>
    <tr>
      <th colspan=3><em>Required columns</em></th>
    </tr>
    <tr>
      <td>--id</td>
      <td>IID</td>
      <td><strong>Required.</strong> Specify the column name in the input demographic
        data file containing participant ID numbers. These should match the IID numbers
        in the files generated by the PGS Catalog Calculator. E.g.<br>--id=idno</td>
    </tr>
    <tr>
      <td>--age</td>
      <td>age</td>
      <td><strong>Required.</strong> Specify the column name in the input demographic
        data file containing the age at recruitment (in years) for each participant.
        E.g.<br>--age=ages</td>
    </tr>
    <tr>
      <td>--sex</td>
      <td>sex</td>
      <td><strong>Required.</strong> Specify the column name in the input demographic
        data file containing the sex of each participant. This column should be coded
        as <em>Male/Female</em>, <em>M/F</em>, <em>1/0</em>, etc. E.g.
        <br>--sex=sexes</td>
    </tr>
    <tr>
      <th colspan=3><em>Binary phenotype data</em></th>
    </tr>
    <tr>
      <td>&#8209;&#8209;any&#8209;case</td>
      <td>NA</td>
      <td>Specify the column name in the input demographic data file containing the overall
        case status of each participant. This column should be coded as <em>1</em> (either
        prevalent or incident case) or <em>0</em> (non-case). If unspecified, comparable
        data will be created using the prevalent and incident case data, if these are both
        supplied. E.g.<br>--any-case=stroke</td>
    </tr>
    <tr>
      <td>&#8209;&#8209;prev&#8209;case</td>
      <td>NA</td>
      <td>Specify the column name in the input demographic data file containing the
        prevalent case status of each participant (i.e. disease status at baseline).
        This column should be coded as <em>1</em> (case) or <em>0</em> (control). E.g.
        <br>--prev-case=prevalent_stroke</td>
    </tr>
    <tr>
      <td>&#8209;&#8209;inc&#8209;case</td>
      <td>NA</td>
      <td>Specify the column name in the input demographic data file containing the incident
        case status of each participant (i.e. disease status after follow-up). This column
        should be coded as <em>1</em> (event) or <em>0</em> (non-event/censored). Prevalent
        cases should be coded as <em>NA</em> or left blank. E.g.
        <br>--inc-case=incident_stroke</td>
    </tr>
    <tr>
      <td>&#8209;&#8209;follow&#8209;time</td>
      <td>NA</td>
      <td>Specify the column name in the input demographic data file containing the
        follow-up time for each participant in the incident cohort. This should contain
        the time to event/ censoring from recruitment (in years). E.g.
        <br>--follow-time=follow_up</td>
    </tr>
    <tr>
      <th colspan=3><em>Quantitative phenotype data</em></th>
    </tr>
    <tr>
      <td>&#8209;&#8209;trait&#8209;val</td>
      <td>NA</td>
      <td>Specify the column name in the input demographic data file containing the
        trait measurement for each participant. These should be the observed (unscaled)
        values. E.g.<br>--trait-val=bmi</td>
    </tr>
    <tr>
      <td>--units</td>
      <td>NA</td>
      <td>Specify the units used for the phenotype measurement. E.g.
        <br>--units=kg/m2</td>
    </tr>
    <tr>
      <th colspan=3><em>Additional data</em></th>
    </tr>
    <tr>
      <td>&#8209;&#8209;covs&#8209;cat</td>
      <td>NA</td>
      <td>Specify the column names in the input demographic data file containing other
        <strong><em>categorical</em></strong> covariates that should be included in
        the regression analyses (sex is included by default). <strong>Note:</strong>
        inclusion of some covariates may lead to warning messages during Cox regression
        modelling, particularly where one or more factor levels show no events, so output
        should be used with caution. E.g.<br>--covs-cat=centre,genotype_batch</td>
    </tr>
    <tr>
      <td>&#8209;&#8209;covs&#8209;quant</td>
      <td>NA</td>
      <td>Specify the column names in the input demographic data file containing other
        <strong><em>quantitative</em></strong> covariates that should be included in the
        regression analyses (age and 10 PCs are included by default). E.g.
        <br>--covs-quant=recruit_year</td>
    </tr>
    <tr>
      <th colspan=3><em>Ancestry stratification</em></th>
    </tr>
    <tr>
      <td>&#8209;&#8209;ancestry</td>
      <td>all</td>
      <td>Specify the subset of ancestry codes to include in the analysis (based on genetic
        similarity with reference dataset). If <em>all</em> (default), all available ancestries
        will be included. For the HGDP+1kGP reference panel, ancestry codes include: AFR
        (African), AMR (Admixed American), CSA (Central and South Asian), EAS (East Asian),
        EUR (European), MID (Middle Eastern). E.g.<br>--ancestry=CSA,EAS</td>
    </tr>
    <tr>
      <td>&#8209;&#8209;stratified</td>
      <td>FALSE</td>
      <td>Specify whether the analysis should be stratified by the ancestries supplied to the
        <code>&#8209;&#8209;ancestry</code> option. If <em>FALSE</em> (default), all listed
        ancestries are analysed together in a combined analysis. E.g.<br>--stratified=TRUE</td>
    </tr>
    <tr>
      <th colspan=3><em>Run specifications</em></th>
    </tr>
    <tr>
      <td>&#8209;&#8209;score&#8209;type</td>
      <td>all</td>
      <td>Specify which score types should be used in the analysis (SUM, Z_MostSimilarPop,
        Z_norm1 or Z_norm2). If <em>all</em> (default), all score types will be included.
        The abbreviations ZMSP, ZN1 and ZN2 can be used. SUM values are centred and scaled
        within each stratum (see the <code>&#8209;&#8209;stratified</code> option above). E.g.
        <br>--score-type=SUM,ZN2</td>
    </tr>
    <tr>
      <td>&#8209;&#8209;min&#8209;cases</td>
      <td>50</td>
      <td>Specify the minimum number of cases required per analysis (for a binary trait). Demographic
        statistics, performance metrics and plots will not be generated for strata comprising fewer
        cases than this. Can also be used to specify the minimum number of total individuals required
        for a quantitative trait. E.g.<br>--min-cases=100</td>
    </tr>
    <tr>
      <td>--skip</td>
      <td>0</td>
      <td>Specify an integer number of PGS IDs to skip before starting analysis. Useful in combination
        with the <code>&#8209;&#8209;max&#8209;pgs</code> option if wanting to analyse a larger dataset
        in smaller groups of PGS IDs. E.g.<br>--skip=50</td>
    </tr>
    <tr>
      <td>&#8209;&#8209;max&#8209;pgs</td>
      <td>NA</td>
      <td>Specify an integer number for the maximum number of PGS IDs to analyse in a single run. If
        <em>NA</em> (default), all PGS IDs present in the input files will be analysed. E.g.
        <br>--max-pgs=50</td>
    </tr>
    <tr>
      <td>&#8209;&#8209;bootstraps</td>
      <td>200</td>
      <td>Specify an integer number of bootstrap replicates to use for generating 95% confidence
        intervals for <em>R</em><sup>2</sup> values (for binary traits). Set to <em>0</em> to
        turn off bootstrapping. E.g.<br>--bootstraps=1000</td>
    </tr>
    <tr>
      <td>--n-cores</td>
      <td>1</td>
      <td>Specify an integer number of cores to use (for parallel processing of PGS IDs). At least
        8 cores are recommended if bootstrapping is to be used. E.g.<br>--n-cores=8</td>
    </tr>
    <tr>
      <td>&#8209;&#8209;run&#8209;dems<br>&#8209;&#8209;run&#8209;metrics<br>
        &#8209;&#8209;run&#8209;plots</td>
      <td>TRUE<br>TRUE<br>TRUE</td>
      <td>Set any of these options to <em>FALSE</em> to prevent the specified analysis (demographics,
        performance metrics or plots) from being run. Can be used for debugging, recreating plots
        without re-running all analyses, etc. E.g.<br>--run-metrics=FALSE</td>
    </tr>
  </tbody>
</table>

## Output

All output files are written to a results directory (`[outdir]/pgsc_benchmark_results/`)
with the following structure. Additional subdirectories and files are created with each
subsequent run of the script, allowing for a full tree structure of results to be
produced (assuming the same output directory is used each time). A record of all
parameters and other options used in each run can be found in the `run_specifications.txt`
file. Depending on which analyses were run, some files or directories listed below may
not be present.

```
pgsc_benchmark_results
    └── [sampleset]
        └── [phenotype]
            └── [run]
                ├── run_specifications.txt
                ├── demographics.csv
                ├── metrics.csv
                └── Plots
                    └── [PGS]_[score type]_[ancestry]_[cohort].png
```

### Demographic data

Summary-level demographic data for the input cohort are reported in the `demographics.csv`
file, along with the mean and standard deviation of each PGS. Depending on the input data
and options selected, statistics will be stratified by ancestry, subcohort
(full/prevalent/incident) and disease status (case/control) (see <em>Subcohorts</em> below).

Each demographics file comprises the following columns:

<table>
  <thead>
    <tr>
      <th>Column&nbsp;name</th>
      <th>Description</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>sampleset</td>
      <td>Cohort/study name of the input dataset</td>
    </tr>
    <tr>
      <td>phenotype</td>
      <td>Disease or trait being examined</td>
    </tr>
    <tr>
      <td>ancestry</td>
      <td>List of ancestry codes used to subset the input data, based on genetic
        similarity to the labelled reference dataset used when running the PGS
        Catalog Calculator</td>
    </tr>
    <tr>
      <td>cohort</td>
      <td>Subset of participants included, based on study design and case/control
        definitions. For a binary trait, this will be one of <em>full</em>,
        <em>prevalent</em> or <em>incident</em>. For a quantitative trait, this
        will always be <em>full</em>. See <em>Subcohorts</em> below</td>
    </tr>
    <tr>
      <td>status</td>
      <td>Subset of participants included, based on disease status. For a binary trait,
        this will be one of <em>any</em>, <em>case</em> or <em>control</em> (or related
        term). This column can be ignored for quantitative traits. See <em>Subcohorts</em>
        below</td>
    </tr>
    <tr>
      <td>n.total</td>
      <td>Total number of participants</td>
    </tr>
    <tr>
      <td>n.males</td>
      <td>Total number of male participants</td>
    </tr>
    <tr>
      <td>age.recruit.mean</td>
      <td>Mean age of participants at recruitment (in years)</td>
    </tr>
    <tr>
      <td>age.recruit.sd</td>
      <td>Standard deviation of the age of participants at recruitment (in years)</td>
    </tr>
    <tr>
      <td>age.event.mean</td>
      <td>Mean age of participants at event or censoring (in years). Only reported for
        incident cohorts</td>
    </tr>
    <tr>
      <td>age.event.sd</td><td>Standard deviation of the age of participants at event
        or censoring (in years). Only reported for incident cohorts</td>
    </tr>
    <tr>
      <td>trait.value.mean</td>
      <td>Mean value of the trait measurement (in the units specified in
        <code>trait.value.units</code>). Only reported for quantitative traits</td>
    </tr>
    <tr>
      <td>trait.value.sd</td>
      <td>Standard deviation of the trait measurement (in the units specified in
        <code>trait.value.units</code>). Only reported for quantitative traits</td>
    </tr>
    <tr>
      <td>trait.value.units</td>
      <td>Units of measurement associated with the values reported in
        <code>trait.value.mean</code> and <code>trait.value.sd</code>. Only
        reported for quantitative traits</td>
    </tr>
    <tr>
      <td>[PGS ID].[score type].mean</td>
      <td>Mean value of the specified PGS. A separate column will be included for each
        PGS ID available in the input data. Abbreviations: ZMSP (Z_MostSimilarPop),
        ZN1 (Z_norm1), ZN2 (Z_norm2)</td>
    </tr>
    <tr>
      <td>[PGS ID].[score type].sd</td>
      <td>Standard deviation of the specified PGS. A separate column will be included
        for each PGS ID available in the input data. Abbreviations: ZMSP (Z_MostSimilarPop),
        ZN1 (Z_norm1), ZN2 (Z_norm2)</td>
    </tr>
  </tbody>
</table>

<table><tr><td><strong>Note:</strong> statistics are calculated after the input
  demographics file has been merged with the PGS and population similarity files,
  and all rows with missing data in essential columns (i.e. phenotype data and all
  covariates needed for regression) have been removed. Therefore, the reported
  number of participants may be lower than expected.</td></tr></table>

<table><tr><td><strong>Note:</strong> if you are seeing NA values where you were
  expecting data, check that the number of cases/total participants in the associated
  stratum is greater than or equal to the number supplied to the <code>--min-cases</code>
  option.</td></tr></table>

### Subcohorts

The following definitions are used to describe subsets of the input data (for a
binary trait):

<table>
  <thead>
    <tr>
      <th>Cohort</th>
      <th>Disease&nbsp;status</th>
      <th>Definition</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td rowspan=2>Full</td>
      <td>Case</td>
      <td>All prevalent and incident cases, regardless of whether the participant
        developed the disease before or after recruitment</td>
    </tr>
    <tr>
      <td>Non-case</td>
      <td>All participants who did not have the disease at recruitment, and who never
        developed the disease during follow up</td>
    </tr>
    <tr>
      <td rowspan=2>Prevalent<br>(baseline)</td>
      <td>Case</td>
      <td>All participants who developed the disease before recruitment (i.e. they
        already had the disease at baseline)</td>
    </tr>
    <tr>
      <td>Control</td>
      <td>All participants who did not have the disease at recruitment, regardless
        of whether they subsequently developed the disease during follow up or not</td>
    </tr>
    <tr><td rowspan=2>Incident<br>(follow&#8209;up)</td>
      <td>Event</td>
      <td>All participants who did not have the disease at recruitment, but who
        subsequently developed the disease during follow up</td>
    </tr>
    <tr>
      <td>Non-event</td>
      <td>All participants who did not have the disease at recruitment, and who
        never developed the disease during follow up</td>
    </tr>
  </tbody>
</table>

### Performance metrics data

The performance metrics for each PGS are reported in the `metrics.csv` file. The type
of metrics calculated will depend on the input data (see below) but will always include
a PGS effect size (e.g. odds ratio, hazard ratio or &beta; coefficient), as well as a
measurement for how well the full model (PGS and covariates) performs compared to a
baseline covariates-only model (e.g. &Delta;AUC, &Delta;C&#8209;index or
&Delta;<em>R</em><sup>2</sup>).

For all analyses, age, sex and 10 genetic principal components are automatically
included as covariates. Additional categorical and quantitative covariates can be
specified using the <code>&#8209;&#8209;covs&#8209;cat</code> or
<code>&#8209;&#8209;covs&#8209;quant</code> options respectively. For single-sex
input, sex will be automatically removed as a covariate.

Each metrics file comprises the following columns:

<table>
  <thead>
    <tr>
      <th>Column&nbsp;name</th>
      <th>Description</th>
    </tr>
  </thead>
  <tbody>  
    <tr>
      <td>sampleset</td>
      <td>Cohort/study name of the input dataset</td>
    </tr>
    <tr>
      <td>phenotype</td>
      <td>Disease or trait being examined</td>
    </tr>
    <tr><td>ancestry</td>
      <td>List of ancestry codes used to subset the input data, based on genetic
        similarity to the labelled reference dataset used when running the PGS
        Catalog Calculator</td>
    </tr>
    <tr>
      <td>pgs.id</td>
      <td>PGS identifier in the PGS Catalog</td>
    </tr>
    <tr>
      <td>score.type</td>
      <td>The score type used to calculate the performance metrics. This will
        be one of <em>SUM</em>, <em>ZMSP</em> (Z_MostSimilarPop), <em>ZN1</em>
        (Z_norm1) or <em>ZN2</em> (Z_norm2)</td>
    </tr>
    <tr>
      <td>cohort</td>
      <td>Subset of participants included, based on study design and case/control
        definitions. For a binary trait, this will be one of <em>full</em>,
        <em>prevalent</em> or <em>incident</em>. For a quantitative trait, this
        will always be <em>full</em>. See <em>Subcohorts</em> above</td>
    </tr>
    <tr>
      <td>model</td>
      <td>Model type used in the regression analysis. This will be one of
        <em>full</em> (the larger model, including the PGS and all covariates)
        or <em>baseline</em> (the reduced model, only including covariates)</td>
    </tr>
    <tr>
      <td>metric</td>
      <td>Type of statistic described by the next four columns. Abbreviations: AUC
        (area under the receiver operating characteristic curve), Beta (&beta; coefficient),
        HR (hazard ratio), OR (odds ratio), R2 (<em>R</em><sup>2</sup>), R2 liability
        (<em>R</em><sup>2</sup> on the liability scale)</td>
    </tr>
    <tr>
      <td>estimate</td>
      <td>Point estimate of the metric</td>
    </tr>
    <tr>
      <td>lower.95</td>
      <td>Lower bound of the 95% confidence interval for the estimate</td>
    </tr>
    <tr>
      <td>upper.95</td>
      <td>Upper bound of the 95% confidence interval for the estimate</td>
    </tr>
    <tr>
      <td>p.val</td>
      <td><em>P</em>-value for the estimate. Only reported for odds ratios, hazard ratios
        and &beta; coefficients. The <em>p</em>-values reported for all <em>delta</em>
        (&Delta;) metrics represent the <em>p</em>-values returned by the likelihood ratio
        Chi-squared test between the full and baseline models</td>
    </tr>
  </tbody>
</table>

#### Statistics reported for binary traits

(a) Full cohort/prevalent cohort:

- Odds ratio for the PGS, with 95% confidence interval and _p_\-value (logistic regression)
- <em>R</em><sup>2</sup> for the PGS on the liability scale (Robertson transformation), with bootstrapped 95% confidence interval
- AUC for the full and baseline logistic regression models, with 95% confidence intervals
- &Delta;AUC between the full and baseline logistic regression models, with 95% confidence interval (DeLong test) and _p_\-value (likelihood ratio Chi-squared test)

(b) Incident cohort (with time-to-event data):

- Hazard ratio for the PGS, with 95% confidence interval and _p_\-value (Cox proportional hazards regression, with time-on-study as timescale)
- C-index for the full and baseline Cox proportional hazards regression models, with 95% confidence intervals
- &Delta;C-index between the full and baseline Cox proportional hazards regression models, with 95% confidence interval and _p_\-value (likelihood ratio Chi-squared test)

#### Statistics reported for quantitative traits

- &beta; coefficient for the PGS, with 95% confidence interval and _p_\-value (linear regression)
- <em>R</em><sup>2</sup> for the full and baseline linear regression models, with 95% confidence intervals
- &Delta;<em>R</em><sup>2</sup> between the full and baseline linear regression models, with 95% confidence interval and _p_\-value (likelihood ratio Chi-squared test)

<table><tr><td><strong>Note:</strong> if time-to-event data is not provided,
  metrics for an incident cohort are calculated using the same statistics as
  (a) above.</td></tr></table>

<table><tr><td><strong>Note:</strong> if you are seeing NA values where you
  were expecting data, check that the number of cases/total participants in the
  associated stratum is greater than or equal to the number supplied to the
  <code>--min-cases</code> option.</td></tr></table>

### Plots

Figures showing the distribution of each PGS in the input dataset are stored as
PNG files in the `Plots/` directory. The type of figures produced depends on the
input data:

- For binary traits, figures will be smoothed kernel density estimates showing the score distributions in cases and controls. Separate figures will be produced for each subcohort (full/ prevalent/incident) in the input data.
- For quantitative traits, figures will be plots of the mean trait values (and 95% confidence intervals) stratified by PGS decile and overlaid with the regression lines from univariate linear regression.

<table><tr><td><strong>Note:</strong> if certain plots appear missing, check
  that the number of cases/total participants in the associated stratum is
  greater than or equal to the number supplied to the <code>--min-cases</code>
  option.</td></tr></table>

#### Binary trait example

<p align="center">
  <img width="75%" src="https://github.com/PGScatalog/pgsc_benchmark/blob/main/figures/Binary_trait.png">
</p>

#### Quantitative trait example

<p align="center">
<img width="75%" src="https://github.com/PGScatalog/pgsc_benchmark/blob/main/figures/Quantitative_trait.png">
</p>
