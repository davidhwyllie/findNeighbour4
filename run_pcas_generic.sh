STARTDATE=${1}
ENDDATE=${2}
echo "Running analysis between ${STARTDATE} and ${ENDDATE}"
echo "Output is in directory ${ENDDATE}" 

mkdir /data/data/pca/fn4_pca6/${ENDDATE}
pipenv run python3 fn4_pca.py demos/covid/atp.json sqlite:////data/data/pca/fn4_pca6/${ENDDATE}/results.sqlite /data/logs/findNeighbour4/localcache/fndev_atptest/pca --n_components 100 \
--focus_on_most_recent_n_days 60 --compute_slope_over 30 --analysis_dir /data/data/pca/fn4_pca6/${ENDDATE} --analysis_date ${ENDDATE} --analysis_window_start_date ${STARTDATE} \
--min_variant_freq 5 --remove_existing_data > pca6_${ENDDATE}.out 
