# Pierroton seedling experiment (2009–2017): plot-level extracts

Built by `lab/data/oak/scripts/pierroton_package_data.py` in the authors' workspace from the
sources below. Plot = 10 × block + unit plot. Treatments: L, myclobutanil; M, natural infection;
H, natural infection plus spores blown on (Desprez-Loustau et al. 2014; Barrès et al. 2024 pool M
and H as "natural").

- `plot_scores.csv`: per plot and mildew score, `n`, `mean`, `sd` and `zeros` of the seedlings'
  scores. `kind` is `percent` (leaf area infected, %), `presence` (0/1), `class` (0–3, 2014) or
  `percent_shoot1/2` (one shoot only). `date` is given (`date_known = 1`) only for the scores
  whose dates are published: 9 June and 15 October 2009, 17 May and 29 September 2010, 4 May
  and 27 July 2011 (Desprez-Loustau et al. 2014, Table 1). `order` is the position within the year.
- `plot_phenology.csv`: mean budburst score of the leading shoot (0 dormant, 1 leaves emerging,
  2 leaves and shoot developing, 3 leaves full size but tender, 4 fully developed, dark green)
  at each weekly scoring (from week 15; exact days unpublished).
- `weather.csv`: daily weather at 44.747 N, −0.767 E, 2008–2017 (local days): mean, min and max
  2 m temperature and mean relative humidity from ERA5-Land, precipitation from ERA5.

## Sources and licences

- Scores and budburst: Barrès B, Desprez-Loustau M-L et al. (2024), Realtime data and code,
  Zenodo, doi:10.5281/zenodo.13626006, **MIT licence**. Paper: Barrès et al. (2024) *Annals of
  Forest Science* 81:44, doi:10.1186/s13595-024-01259-2.
- Dates, design and scales: Desprez-Loustau M-L, Saint-Jean G, Barrès B, Dantec CF, Dutech C
  (2014) *Annals of Forest Science* 71:563–573, doi:10.1007/s13595-014-0364-6.
- Weather: **CC BY 4.0**. Weather data by Open-Meteo.com (https://open-meteo.com). Contains
  modified Copernicus Climate Change Service information (ERA5, ERA5-Land; Hersbach et al. 2020,
  Muñoz-Sabater et al. 2021).
