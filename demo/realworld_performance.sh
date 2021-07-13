# determine realtime performance of the detector
for start_date in 2020-09-23 2020-09-30 2020-10-01 2020-10-08 2020-10-11 2020-10-15 2020-10-20 2020-10-21 2020-10-22 2020-10-23 2020-10-24 2020-10-25 2020-10-26 2020-06-03 2020-06-17 2020-07-02 2020-07-18 2020-07-21 2020-07-24 2020-07-25 2020-07-26 2020-07-29 2020-07-31 2020-08-01 2020-08-02 2020-08-03 2020-08-04 2020-08-05 2020-08-06 2020-08-07 2020-08-09 2020-12-11 2020-12-12 2020-12-14 2020-12-15 2020-12-16 2020-12-17 2020-12-18 2020-12-21 2020-12-22 2020-12-24 2020-12-26 2020-12-27 2020-12-28 2020-12-29 2020-12-31 2021-01-01 2021-01-03 2021-01-04;

do 
    echo $start_date
    pipenv run python3 fn4_pca.py demos/covid/covid_config_v3.json sqlite:////data/data/pca/rwd/${start_date}.sqlite --n_components 400 --focus_on_most_recent_n_days 90 --compute_slope_over 30 --analysis_dir /data/data/pca/rwd --analysis_date ${start_date} --remove_existing_data

done
