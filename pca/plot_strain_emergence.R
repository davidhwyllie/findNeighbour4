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
        metavar="number of days prior to today that the Poisson regression models examine",
    help="interval to compute trends over")
# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults, 
# example usage: 
# Rscript plot_strain_emergence_realtime.R --filepattern "2021-04-27.sqlite" --interval 30 --searchdir /data/data/pca/realtime_400/ --cogfile /data/data/inputfasta/cog_metadata.csv

args <- parser$parse_args()

## -----    example usage: loading, iteration  ------------------------
# set up parameters
# set up parameters
BASE_DIR <- args$searchdir   # "/data/data/pca/realtime_400/"  # where the databases are
cogfile <- args$cogfile #  "/data/data/inputfasta/cog_metadata.csv"
interval <- args$interval

# find all sqlite dbs
glob_path <- paste0(BASE_DIR,args$filepattern)
print(paste0("Searching ",glob_path))
for (dbfile in sort(Sys.glob(glob_path), decreasing=FALSE)){
    
  db_stem <- substr(basename(dbfile), 1, 10) # used for output; the way these files are 
                                             # labelled, this is the date.
  print(paste0("End date is ",db_stem))

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



  # strains from https://www.gov.uk/government/publications/covid-19-variants-genomically-confirmed-case-numbers/variants-distribution-of-cases-data
  # 'B.1.1.7','B.1.351','P.1','B.1.1.318',,'B.1.617','P.1'
  test_lineages <- c('B.1.1.7','B.1.351','B.1.1.318','B.1.617','P.1')
  #test_lineages <- c('B.1.351')
  for (lineage_to_depict in  test_lineages){

    ## start of future function
    # parameters controlling enhanced depiction options
    # lineage_to_depict
    min_estimate_IRR <- 3
    max_samples <- 300
    min_recent_proportion <- 0.5

    # function starts
    res <- plot_counts_per_lineage(db_connection = db_conn,
                                        target_lineage=lineage_to_depict, 
                                        date_end = db_stem,
                                        axis_end = Sys.Date(),
                                        axis_start = "2020-06-01"
    )

    p0 <- res[['p']]

    res <- plot_pc_cat_size_vs_change_marking_selected_lineages(
      db_connection =db_conn,
      date_end = db_stem,
      reporting_or_cutoff = 100,
      target_lineages = c(lineage_to_depict),
      only_show_significant_trending = TRUE) 
    
    p1 <- res[['p']]
    p1 <- p1 + geom_vline(xintercept = max_samples, lty =2)
    p1 <- p1 + geom_hline(yintercept = min_estimate_IRR, lty=2)

    trends <- res[['trends']]

    trends$recent_proportion <- 1-(trends$n_gt_30_days_before / trends$n_samples)
    trends <- subset(trends, 
        is_sig==1 & 
        natural_space_Estimate >= min_estimate_IRR & 
        recent_proportion >= min_recent_proportion & 
        n_samples <= max_samples)
      selection_criterion <- paste0(
              "* IRR estimate > ",
              min_estimate_IRR,
              " \n",
              "* prop seen in last 60d >= ",
              min_recent_proportion,
              " of total \n* < ",
              max_samples,
            " samples total"
            )

    trends$pc_cat <- factor(trends$pc_cat, levels=pc_cat_levels(trends$pc_cat))
    date_ends <- convert_date(db_stem)
    trends$marker_x <- date_ends[['dt']]-interval
    trends$marker_xend <- date_ends[['dt']]
    trends$marker_y <- 0
    trends$marker_yend <-0

    print(paste0("Trends highlighted = ", nrow(trends)))
    print(paste0("Highlighted pc cats: ", paste(trends$pc_cat, collapse=';')))
    print(paste0("Lineage associated pc cats: ", paste(res$associated_pc_cats, collapse=';')))

    # load data
    analyse_pc_cats <- unique(trends$pc_cat)
    read_pcs <- unique(trends$pc)
    cnts_overall <- list()

    print("Analysing these pc_cats")
    analyse_pc_cats_df <- data.frame(pc_cat = analyse_pc_cats)

    ############################ basis of an external function #################################

    ### basis of external function possibly without the option to highlight a single association

    # recover data about lineage associations
    cmd <- "
        select * 
        from feature_associations
       
        order by feature, pc_cat;"
          
      # get long dataset
    assoc_long <- database_query(
        db_connection = db_conn, 
        cmd = cmd) # recover everything, irrespective of p-value
    assoc_long <- subset(assoc_long, pc_cat %in% analyse_pc_cats)
    assoc_long$pc_cat <- factor(assoc_long$pc_cat, levels = pc_cat_levels(assoc_long$pc_cat))
    assoc_long$log_or <- log10(assoc_long$truncated_odds_ratio)

    assoc_long$proportion_feature_explained <- assoc_long$a/ assoc_long$a_c
    assoc_long$is_lineage <- ifelse(grepl("lineage:", assoc_long$feature), "Lineage", "Variant")
    assoc_long$feature <- gsub("lineage:","", assoc_long$feature, fixed=TRUE)
    assoc_long <- subset(assoc_long, !feature == "NA")

    # recover data about counts
    cnt_result <- count_per_pc_cats_in_time_interval(
        db_connection = db_conn,
        these_pc_cats = analyse_pc_cats,
        date_end = date_ends[['dt']], 
        interval_analysed= 3*interval
      )

    ## if there are counts
    if (cnt_result$success){
      # there are trending pc_cats
      cnts <- cnt_result$count_df
      cnts$pc_cat <- factor(cnts$pc_cat, levels = pc_cat_levels(cnts$pc_cat))
      print('plotting counts')
      
      p <- ggplot(cnts, aes(x=sample_date, y= n))+ theme_classic()
      p <- p + geom_col()
      p <- p + geom_hline(colour='black', lty=2,yintercept = c(5))
      p <- p + scale_x_date(limits = c(date_ends[['dt']]-4*interval, date_ends[['dt']])) 
      p <- p + geom_vline(colour='black', lty=1, xintercept = c( date_ends[['dt']]-4*interval))
      p <- p + geom_vline(colour='red', lty=1, xintercept = c(date_ends[['dt']]))

      p <- p + geom_hline(colour='black', lty=1,yintercept = 0)
      p <- p + geom_hline(colour='black', lty=2,yintercept = c(5))
      p <- p + facet_wrap(~pc_cat, scale='free_y')
      p <- p + theme(
          axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          axis.ticks.x=element_blank(),
          axis.ticks.y=element_blank(),
          legend.position="none",
          panel.background=element_blank(),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank()
      )

      p2 <- p + ggtitle(
        paste0("Highlighted count profiles for:\n",
        selection_criterion,
        "\nBar: last ",
        interval,
        ",days\n",
        "Gray/Red bar: not assoc/assoc with ",
        lineage_to_depict,
        "\n",
        "Horizontal dashes: 5 cases/day")
      )

      p2 <- p2 + geom_col(colour = 'black')
      p2 <- p2 + geom_segment(colour='gray', size= 3, lty=1, aes(x = marker_x, y=marker_y, xend=marker_xend, yend=marker_yend), data = trends)

      # and depict associations
      # check associations
      depiction_or_cutoff <- 1e2

      print("Checking associations with lineages")
      check_table_present(db_conn, 'feature_associations','run make_contingency_tables() first')
      depicted_associations <- subset(assoc_long,  truncated_odds_ratio > depiction_or_cutoff)

      trend_subset <- subset(trends, pc_cat %in% unique(depicted_associations$pc_cat))
      trend_subset <- subset(trend_subset, lineage == lineage_to_depict)
     
      p2 <- p2 + geom_segment(colour='red', size= 3, lty=1, aes(x = marker_x, y=marker_y, xend=marker_xend, yend=marker_yend), data = trend_subset)


      p3 <- ggplot(assoc_long, aes(y=pc_cat, x=feature))
      p3 <- p3 + geom_blank()   # fill the grid
      p3 <- p3 + geom_tile(aes(fill=log_or), data=subset(assoc_long, truncated_odds_ratio > depiction_or_cutoff))
      p3 <- p3 + geom_point(aes(size=proportion_feature_explained), colour='yellow', data= depicted_associations)
      p3 <- p3 + ggtitle("Sequence patterns\nassociations with lineage")
      p3 <- p3 + scale_y_discrete("Sequence pattern (pc_cat)")
      p3 <- p3 + scale_x_discrete("Lineage or mutation")
      p3 <- p3 + theme(axis.text.x=element_text(angle=90))
      p3 <- p3 + scale_fill_continuous(name = "Log OR favouring\nlineage given sequence\npattern")
      p3 <- p3 + theme(legend.position="top", legend.direction = 'vertical', legend.box='horizontal')
      p3 <- p3 + scale_size_continuous("Proportion of\nlineage comprising\nsequence pattern")
      p3 <- p3 + facet_grid(.~is_lineage, scales='free')

    } else {
        p2 <- ggplot() + theme_void()   # nothing
        p3 <- ggplot() + theme_void()   # nothing

    }
      
    print("Rendering to:")
    lineage_filename <- gsub('.','_',lineage_to_depict, fixed= TRUE)
    outfile <- paste0(PLOT_DIR,'/pc_cats_cf_', lineage_filename,"_",db_stem,".tiff")
    print(outfile)
    tiff(outfile,width=1600,height=1000, compression='lzw')
    gridExtra::grid.arrange(grobs = list(p0,p1,p2,p3), nrow=1, ncol=4, widths=c(1.5,1.5,2,3) )
    dev.off()
    print("Plot rendered")
  }

  
  # example code to dump trending sample lists
  print("Dumping contents of trending pc_cats to file")
  for (this_pc_cat in unique(assoc_long$pc_cat)){
    print(this_pc_cat)
    sqlcmd = "select distinct pc_cat, cm.sample_id, cm.sample_date from transformed_coordinate_categories tcc
              inner join clinical_metadata cm 
              on cm.sample_id = tcc.sample_id
              inner join sequence_metadata sm
              on sm.sample_id = tcc.sample_id
              where tcc.pc_cat = ?; "

    res <- database_query(db_conn, sqlcmd, params=list(this_pc_cat))
    outfile <- paste0(PLOT_DIR,"/members_of_",this_pc_cat,".tsv")
    write.table(res, file=outfile, sep="\t", row.names=FALSE)
  }

  dbDisconnect(db_conn)

}
print('FINISHED')
