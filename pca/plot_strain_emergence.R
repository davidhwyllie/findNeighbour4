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
parser$add_argument('-d', "--searchdir", default="/data/data/pca/realtime_400/", type="character",
        metavar="directory in which to search for sqlite output",
        help="directory in which to search for sqlite files. Example: /data/pca/sqlite")

parser$add_argument('-f', "--filepattern", default="*.sqlite", type="character",
        metavar="pattern to search for in order to identify sqlite files",
        help="file pattern to glob to look for sqlite files. Example: 2021*.sqlite")

parser$add_argument('-m', "--cogfile", default="/data/data/inputfasta/cog_metadata.csv", type="character",
        metavar="metadata file; format used by cog-uk is currently expected",
        help="metadata file")

parser$add_argument('-i', "--interval", default=30, type="integer",
        metavar="number of days prior to today that the regression models examine",
        help="interval to compute trends over")

parser$add_argument('-e', "--highlight_min_estimate_IRR", default=3, type="double",
        metavar="highlight trending pc_cats with incidence rate ratio more than ",
        help="illustrate trending pc_cats above a rate ratio cutoff.  Suitable values 2-10")

parser$add_argument('-n', "--highlight_max_samples", default=300, type="double",
        metavar="highlight trending pcs_cats with fewer than max_samples cases ",
        help="highlight trending pc_cats with fewer than highlight_max_samples cutoff.  Suitable values 30 - 300")

parser$add_argument('-r', "--highlight_min_recent_proportion", default=0.5, type="double",
        metavar="highlight trending pcs_cats where > than highlight_min_recent_proportion cases are within --interval days",
        help="highlight trending pcs_cats where > than highlight_min_recent_proportion cases are within --interval days")

parser$add_argument('-l', "--lineage_association_or_cutoff", default=1e3, type="double",
        metavar="report associations if the presence vs. absence of a pc_cat modifies the odds ratio of a lineage being present vs. not being present by > this cutoff",
        help="report associations if the presence vs. absence of a pc_cat modifies the odds ratio of a lineage being present vs. not being present by > this cutoff")

parser$add_argument('-a', "--only_pc_cat_arising_after", default='', type="character",
        metavar="report associations if the earliest sample in the pc_cat is after an isoformat date (e.g. 2021-03-01).  Significantly speeds computation if set.",
        help="report associations if the earliest sample in the pc_cat is after an isoformat date (e.g. 2021-03-01).  Significantly speeds computation if set.")

parser$add_argument('-q', "--lineage_to_highlight", default="", type="character",
        metavar="highlight counts of, and associations with, lineage_to_depict",
        help="highlight counts of, and associations with, lineage_to_depict.  If omitted, no strains are highlighted")

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults, 
# example usage: 
# Rscript plot_strain_emergence_realtime.R --filepattern "2021-04-27.sqlite" --interval 30 --searchdir /data/data/pca/realtime_400/ --cogfile /data/data/inputfasta/cog_metadata.csv
# Rscript plot_strain_emergence.R --searchdir /data/data/pca/cp_subsets_output_400_plus/ --filepattern 2020*.sqlite --interval 30 -n 100 -r 0.5 -e 5


args <- parser$parse_args()

## -----    example usage: loading, iteration  ------------------------
# set up parameters 
BASE_DIR <- args$searchdir   # "/data/data/pca/realtime_400/"  # where the databases are
cogfile <- args$cogfile        #  "/data/data/inputfasta/cog_metadata.csv"
interval <- args$interval
min_estimate_IRR <- args$highlight_min_estimate_IRR
max_samples <- args$highlight_max_samples
min_recent_proportion <- args$highlight_min_recent_proportion  
filepattern <- args$filepattern
depiction_or_cutoff <- args$lineage_association_or_cutoff
lineage_to_highlight <- args$lineage_to_highlight
if (nchar(lineage_to_highlight)==0){
        print("Set lineage_to_highlight to NA")
        lineage_to_highlight <- NA
} else {
        print(paste0("Highlighting >",lineage_to_highlight,"<"))
}
only_pc_cat_arising_after <- NA
if (nchar(args$only_pc_cat_arising_after)>0){
        print(paste0("Will only report pc_cats associations if the pc_cat arose after ", args$only_pc_cat_arising_after))
        only_pc_cat_arising_after <- args$only_pc_cat_arising_after
} else {
        print("Will examine associations with lineage for all pca_cats.")
}

print("Subset selection criteria")
print(paste(min_estimate_IRR, max_samples,min_recent_proportion))

# find all sqlite dbs
glob_path <- paste0(BASE_DIR,filepattern)
print(paste0("Searching ",glob_path))
print(Sys.glob(glob_path))
for (dbfile in Sys.glob(glob_path)){
  print(dbfile)
  db_stem <- substr(basename(dbfile), 1, 10) # used for output; the way these files are 
                                             # labelled, this is the date.
  print(paste0("End date is ",db_stem))

  PLOT_DIR <- paste0(BASE_DIR,"pc_plots/strain_emergence") # where the output goes
  dir.create(PLOT_DIR, recursive=TRUE)
  
  FIT_DIR <- paste0(BASE_DIR,"pc_plots/fits_failed") # where the output goes
  dir.create(FIT_DIR, recursive=TRUE)
  
  db_conn <- connect_to_fn4sqlitedb(dbfile) # connect
  
  # ensure data is processed
  print("Adding cog data")
  add_cog_metadata(db_conn, cogfile, db_stem)  # add cog data 

  # recover counts to model
  print("Ensuring model is fitted")
  fit_recent_trend_in_counts(db_connection = db_conn, 
                                                analysis_family_id = db_stem,
                                                date_end = db_stem,
                                                interval_analysed = interval,
                                                fitting_failed_dir = FIT_DIR,
                                                )

  make_contingency_tables(db_conn, 
                          date_end = db_stem,
                          overwrite=FALSE,
                          only_pc_cat_arising_after= only_pc_cat_arising_after)          # debug
  # example strains from https://www.gov.uk/government/publications/covid-19-variants-genomically-confirmed-case-numbers/variants-distribution-of-cases-data
  # 'B.1.1.7','B.1.351','P.1','B.1.1.318',,'B.1.617','P.1'

  retVal <- depict_current_trends(
    db_connection = db_conn,
    db_stem = db_stem,
    PLOT_DIR= PLOT_DIR,
    lineage_to_depict = lineage_to_highlight,
    export_sequence_identifiers = TRUE,
    axis_start = '2020-06-01',
    axis_end = Sys.Date(),
    max_samples = max_samples
  )
  
  if (FALSE){
        for (pc in 1:number_of_pcs(db_conn)){
        print(paste0("Depicting ",pc))
        res <- plot_single_pc_format_3(db_conn, 
                                        pc=pc, 
                                        date_end=Sys.Date(),
                                        time_interval=interval)
        to_tiff(to_plot= res$p2, outputdir=PLOT_DIR, db_name=db_stem, plot_name=paste0("tiled_pc_",pc))
        
        p3 <- plot_pc_format_3(db_conn, pc, date_end=Sys.Date(), time_interval=interval)
        to_tiff(to_plot= p3, outputdir=PLOT_DIR, db_name=db_stem, plot_name=paste0("barchart_pc_",pc))
        
        }
  }
 dbDisconnect(db_conn)
}
print('FINISHED')
print("Warnings generated are as follows:")
warnings()