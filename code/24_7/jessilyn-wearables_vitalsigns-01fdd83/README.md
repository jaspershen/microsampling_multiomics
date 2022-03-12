# Wearable sensors enable personalized predictions of clinical laboratory measurements

Vital signs, including heart rate and body temperature, are useful in detecting or monitoring medical conditions, but are typically measured in the clinic and require follow-up laboratory testing for more definitive diagnoses. Here we examined whether vital signs as measured by consumer wearable devices (i.e. continuously monitored heart rate, body temperature, electrodermal activity, and movement) can predict clinical laboratory test results using machine learning models including random forests and LASSO.

See more details in our Nature Medicine paper.

## Scripts

All code for building and validating models as well as generation of figures is stored `PaperCode.R`.

## Data

Scripts require clinical data to run. They need to be stored in `SECURE_data` directory.

| File | Description |
|:------------- |:-------------|
| Basis2016_Cleaned_NotNorm0824_WeekPrior.csv | iPOP wearables/clinical combined data |
| vitals.csv | iPOP vitals |
| lab_results_20170717.csv | iPOP Labs |
| all_vitals.csv | SEHR vitals |
| all_labs.csv | SEHR labs |
| 20170905_Cleaned_joined_30k_labs_vitals.csv | cleaned SEHR file |
