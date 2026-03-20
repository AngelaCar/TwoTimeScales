# Data from the chemotherapy for stace B/C colon cancer study

This dataset is a reduced version of the dataset colon from the package
survival (Therneau, T., 2023). Each observation is a transition from
recurrence of colon cancer to death or censoring. The time scales are
time from randomization to recurrence, time from randomization to death
or censoring and time from recurrence of the cancer to death or
censoring. Only observations about individuals with a recurrence of the
cancer are selected. Additionally, 7 individuals with exit times from
the risk set equal to entry times in the recurrence state (0 exposure
time) were dropped from the sample. In the original dataset, all times
of recurrence are known, so that after recurrence all individuals are
followed from entry in the state, without left truncation. To be able to
illustrate how to include left truncated times in the model, artificial
left truncated entry in the 'recurrence' state are randomly introduced
for 40 individuals.

## Usage

``` r
data(reccolon2ts)
```

## Format

### reccolon2ts A data.table with 461 rows and 25 columns:

- id:

  patient's id

- study:

  1 for all patients

- rx:

  Treatment - Obs(ervation), Lev(amisole), Lev(amisole)+5-FU

- sex:

  1=male, 0=female

- age:

  Age at transplant in years

- obstruct:

  obstruction of colon by tumor: 1=yes

- perfor:

  perforation of colon: 1=yes

- adhere:

  adherence to nearby organs: 1=yes

- nodes:

  number of lymph nodes with detectable cancer

- status:

  censoring status: 0=censoring, 1=event

- differ:

  differentiation of tumour: 1=well, 2=moderate, 3=poor

- extent:

  extent of local spread: 1=submucosa, 2=muscle, 3=serosa, 4=contigous
  structures

- surg:

  time from surgery to registration: 0=short, 1=long

- node4:

  more than 4 positive lymph nodes

- etype:

  2 for all patients (2=death)

- timedc:

  time in days from randomization to death or censoring

- timer:

  time in days from randomization to recurrence

- timesr:

  time in days from recurrence to death or censoring

- entrys:

  artificial entry time on the time since recurrence scale. For most of
  the individual this is 0 (no left truncation). For 40 individuals a
  random number between 1 and the exit time on the time since recurrence
  scale (timesr) is simulated.

- entryt:

  time in days from randomization to observation in the recurrence
  state. If the individual is observed from entry in the recurrence
  state this is equal to the time at recurrence. If the entry in the
  recurrence state is not observed from the beginning, left truncation
  is observed. This is not present in the original data, but has been
  here introduced artificially for 40 individuals. This is done by first
  increasing the time at recurrence by a random number between 1 and the
  exit time on the time since recurrence scale. Then, the time at
  recurrence is added to the artificial entry time.

- timedc_y:

  time in years from randomization to death or censoring

- timer_y:

  time in years from randomization to recurrence

- entrys_y:

  left truncated entry in the recurrence state measured in years since
  recurrence

- entryt_y:

  left truncated entry in the recurrence state measured in years since
  randomization

- timesr_y:

  time in years from recurrence to death or censoring

## Source

Therneau, T. (2023). A Package for Survival Analysis in R. R package
version 3.5-3, <https://CRAN.R-project.org/package=survival>

## References

Moertel, C.G, et al. (1995). Fluorouracil plus Levamisole as Effective
Adjuvant Theraphy after Resection of Stage III Colon Carcinoma: A Final
Report. Annals of Internal Medicine, 122:321-326

Moerel, C.G., et al. (1990). Levamisole and Fluorouracil for Adjvant
Theraphy of Resected Colon Carcinoma. The New England Journal of
Medicine, 322:352-8

## Examples

``` r
data(reccolon2ts)
rm(reccolon2ts)
#> Warning: object 'reccolon2ts' not found
```
