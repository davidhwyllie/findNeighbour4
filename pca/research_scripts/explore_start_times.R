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
library(argparse)
# functions
source('pca_depiction_functions.R')

parser <- ArgumentParser()
parser$add_argument('-f', "--filepattern", default="*.sqlite", type="character",
        metavar="pattern to search for",
    help="file pattern to glob to look for sqlite files")

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults, 
args <- parser$parse_args()

## -----    example usage: loading, iteration  ------------------------
# set up parameters
# set up parameters
BASE_DIR <- "/data/data/pca/subsets_output_pca2_plus/"  # where the databases are
cogfile <- "/data/data/inputfasta/cog_metadata.csv"
interval <- 30

# find all sqlite dbs
glob_path <- paste0(BASE_DIR,args$filepattern)
print(paste0("Searching ",glob_path))
for (dbfile in sort(Sys.glob(glob_path), decreasing=FALSE)){
    
  db_stem <- substr(basename(dbfile), 3, 12) # used for output; the way these files are 
                                             # labelled, this is the date.
  PLOT_DIR <- paste0(BASE_DIR,"pc_plots/strain_emergence") # where the output goes
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

  trends <- load_trend_summary(db_conn)
  t1 <- subset(trends, param=='Rate ratio relative to ref at baseline')
  head(t1)

  print(t1[1:5,])
  print(table(t1$earliest_date))       # controls from a date window round the date - symmetric? will have > 20 controlsr 

  #get a list of similar 

  dbDisconnect(db_conn)
  stop()
}
stop('FINISHED')
