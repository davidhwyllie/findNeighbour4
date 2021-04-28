# construct analysis name
cd /data/software/fn4dev
analysis_name=$(date +%F); echo $analysis_name
# run pca
pipenv run python3 pca/fn4_pca.py demos/covid/covid_config_v3.json --outputdir /data/data/pca/realtime_400 --analysis_name $analysis_name --n_components 400
