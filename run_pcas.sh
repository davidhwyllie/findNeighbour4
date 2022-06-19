
mkdir /data/data/pca/fn4_pca6/2021-11-01
pipenv run python3 fn4_pca.py demos/covid/atp.json sqlite:////data/data/pca/fn4_pca6/2021-11-01/2021-11-01.sqlite /data/logs/findNeighbour4/localcache/fndev_atptest/pca --n_components 100 \
--focus_on_most_recent_n_days 60 --compute_slope_over 30 --analysis_dir /data/data/pca/fn4_pca6 --analysis_date 2021-11-01 --analysis_window_start_date 2021-07-01 \
--remove_existing_data > pca6_2021-11-01.out 

mkdir /data/data/pca/fn4_pca6/2021-10-01
pipenv run python3 fn4_pca.py demos/covid/atp.json sqlite:////data/data/pca/fn4_pca6/2021-10-01/2021-10-01.sqlite /data/logs/findNeighbour4/localcache/fndev_atptest/pca --n_components 100 \
--focus_on_most_recent_n_days 60 --compute_slope_over 30 --analysis_dir /data/data/pca/fn4_pca6 --analysis_date 2021-10-01 --analysis_window_start_date 2021-06-01 \
--remove_existing_data > pca6_2021-10-01.out 

