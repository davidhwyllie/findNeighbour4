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

library(gridExtra)
# functions
source('pca_depiction_functions.R')

## -------------------------- example usage ------------------------
# set up parameters
BASE_DIR <- "/data/data/pca/subsets_output"  # where the databases are
cogfile <- "/data/data/inputfasta/cog_metadata.csv"
interval <- 30

dbfile <- '/data/data/pca/subsets_output/0-2021-03-29.sqlite' 
db_stem <- substr(basename(dbfile), 3, 12) # used for output; the way these files are 
                                            # labelled, this is the date.
PLOT_DIR <- paste0("/data/data/pca/pc_plots/",db_stem) # where the output goes
dir.create(PLOT_DIR, recursive=TRUE)  
db_conn <- connect_to_fn4sqlitedb(dbfile) # connect
add_cog_metadata(db_conn, cogfile)  # add cog data 

fit_recent_trend_wrt_most_common_pcat_all_pcs(db_connection = db_conn, 
                                            analysis_family_id = db_stem,
                                            date_end = db_stem,
                                            interval_analysed = interval,
                                            max_pcs = NA,
                                            overwrite= FALSE)

# compute contingency tables
make_contingency_tables(db_conn, overwrite=FALSE)

### example for all single components
for (pc in 1:number_of_pcs(db_conn)) {
    print(paste("Exporting trend for pc ",pc))
    p3 <- plot_single_pc_format_2(
      db_connection = db_conn, 
      pc = pc,
      date_end = db_stem) 
    to_tiff(to_plot=p3, 
            outputdir=PLOT_DIR, 
            db_name=db_stem, 
            plot_name=paste0('trend_pc_',pc))
  }

  dbDisconnect(db_conn)
  
}
print('FINISHED')

