# methodscompare methods of clustering of transformed_coordinates
# resulting from pca of SARS-CoV-2 genomes.

# packages code from Brice Letcher, refactored to generate a 
# database backed application with a range of helpful functions
       
# consumes data in SQLite format exported by fn4pca.py, a component 
# of the findNeighbour4 server system

# note: if would be possible to code all this as part of a pure python
# web application, with the possible exception of some of the graphics produced
# by ggtree

rm(list=ls())

# functions
source('pca_depiction_functions.R')

## -----    example usage: loading, iteration  ------------------------
# set up parameters
# set up parameters
BASE_DIR <- "/data/data/pca/subsets_output"  # where the databases are
cogfile <- "/data/data/inputfasta/cog_metadata.csv"
interval <- 30
# find all sqlite dbs
glob_path <- paste0(BASE_DIR,"/0*.sqlite")
for (dbfile in sort(Sys.glob(glob_path), decreasing=FALSE)){
    
  db_stem <- substr(basename(dbfile), 3, 12) # used for output; the way these files are 
                                             # labelled, this is the date.
  PLOT_DIR <- paste0("pc_plots/strain_emergence") # where the output goes
  dir.create(PLOT_DIR, recursive=TRUE)
  
  db_conn <- connect_to_fn4sqlitedb(dbfile) # connect
  
  # ensure data is processed
  add_cog_metadata(db_conn, cogfile, db_stem)  # add cog data 
  fit_recent_trend_wrt_most_common_pcat_all_pcs(db_connection = db_conn, 
                                                analysis_family_id = db_stem,
                                                date_end = db_stem,
                                                interval_analysed = interval,
                                                max_pcs = NA,
                                                overwrite = FALSE)
  make_contingency_tables(db_conn, overwrite=FALSE)
  
  # ------- example plots -----
  
  # strains from https://www.gov.uk/government/publications/covid-19-variants-genomically-confirmed-case-numbers/variants-distribution-of-cases-data
  
  p1 <- plot_pc_cat_size_vs_change_marking_selected_lineages(
    db_connection =db_conn,
    date_end = db_stem,
    target_strains = c('B.1.1.7','B.1.351','P.1','P.1','B.1.1.318','P.3','B.1.617'),
    only_show_significant_trending = FALSE) 
  p1 <- p1 + ggtitle(db_stem) 
  to_tiff(to_plot=p1, 
            outputdir=PLOT_DIR, 
            db_name=db_stem, 
            plot_name='pc_cats_size_vs_change_all_marking_selected_lineages')
    
  dbDisconnect(db_conn)
  
  }
stop('FINISHED')
