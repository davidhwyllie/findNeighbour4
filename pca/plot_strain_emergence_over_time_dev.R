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

## -----    example usage: loading, iteration  ------------------------
# set up parameters
# set up parameters
BASE_DIR <- "/data/data/pca/subsets_output/"  # where the databases are
cogfile <- "/data/data/inputfasta/cog_metadata.csv"
interval <- 30

# find all sqlite dbs
glob_path <- paste0(BASE_DIR,"0-*.sqlite")
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



  # strains from https://www.gov.uk/government/publications/covid-19-variants-genomically-confirmed-case-numbers/variants-distribution-of-cases-data
  # 'B.1.1.7','B.1.351','P.1','B.1.1.318',,'B.1.617','P.1'
  for (lineage_to_depict in  c('P.1')){

    res <- plot_counts_per_lineage(db_connection = db_conn,
                                        target_lineage=lineage_to_depict, 
                                        date_end = db_stem,
                                        axis_end = "2021-04-01",
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
 
    trends <- res[['trends']]
    trends <- subset(trends, isSig==1 & tr_IRR> 3)

    ################################ TO ENCAPSULATE IN FUNCTION LATER ###########################################
    analyse_pc_cats <- unique(trends$pc_cat)
    read_pcs <- unique(trends$pc)
    cnts_overall <- list()
    date_ends <- convert_date(db_stem)
    print("Extracting counts for trending pcs")
    for (pc in read_pcs){
    
      cnts_overall[[pc]] <- count_per_pc_in_time_interval(
        db_connection = db_conn,
        this_pc = pc,
        date_end = date_ends[['dt']], 
        interval_analysed= 3*interval
        )
      cnts_overall[[pc]] <- subset(cnts_overall[[pc]]$count_df, pc_cat %in% analyse_pc_cats )
    }
    #print(names(cnts_overall))
    cnts <- do.call(rbind, cnts_overall)
    #print(cnts[1:10])

    #tmp <- count_per_pc_cats_in_time_interval(
    #    db_connection = db_conn,
    #    these_pc_cats = analyse_pc_cats,
    #    date_end = date_ends[['dt']], 
    #    interval_analysed= 30
    #    )
    cnts <- subset(cnts, pc_cat %in% analyse_pc_cats)
    
    # for each pc_cat, identify whether we are interested in pc_cats seen historically.  
    # let's consider 'historical' as interval.
    cnts$sample_date <- as.Date(cnts$sample_date)
    cnts$delta <- difftime( date_ends[['dt']] , cnts$sample_date, units = "days")
    cnts$older <- ifelse(cnts$delta > interval, "older", "recent")
    cnts_older <- reshape2::dcast(cnts, pc_cat ~ older , value.var= 'n', sum)
    cnts_older$ratio <- cnts_older$recent / (cnts_older$older + cnts_older$recent)
    
    recent_pc_cats <- cnts_older$pc_cat[cnts_older$ratio > 0.7]
    cnts <- subset(cnts, pc_cat %in% recent_pc_cats)
    p <- p + geom_hline(colour='black', lty=2,yintercept = c(5,10))
    p <- p + geom_hline(colour='black', lty=2,yintercept = c(5,10))
    label_horizontal_lines <- c(x= c(date_ends[['dt']]-interval, date_ends[['dt']]-interval), y=c(5,10), label=c(5,10))
    p <- ggplot(cnts, aes(x=sample_date, y= n))
    p <- p + geom_segment(colour='green', size= 5, lty=1, x = (date_ends[['dt']]-interval), xend=date_ends[['dt']],y=0, yend=0)
    p <- p + scale_x_date(limits = c(date_ends[['dt']]-4*interval, date_ends[['dt']])) 
    p <- p + geom_vline(colour='black', lty=1, xintercept = c( date_ends[['dt']]-4*interval))
    p <- p + geom_vline(colour='red', lty=1, xintercept = c(date_ends[['dt']]))
    p <- p + geom_hline(colour='black', lty=1,yintercept = 0)
    p <- p + geom_hline(colour='black', lty=2,yintercept = c(5,10))
    p <- p + geom_text(aes(x=x,y=y,label=label), data=label_horizontal_lines)
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
    p2 <- p + ggtitle(paste0("Count profiles over last 4 months\nFor trending pc_cats in last ",interval,"d (green bar)"))
    p2 <- p2 + geom_col()
    

    lineage_filename <- gsub('.','_',lineage_to_depict, fixed= TRUE)
    outfile <- paste0(PLOT_DIR,'/pc_cats_cf_', lineage_filename,"_",db_stem,".tiff")
    print(outfile)
    tiff(outfile,width=1200,height=1000, compression='lzw')
    #gridExtra::grid.arrange(p0,p1, nrow=1, ncol=2)
    gridExtra::grid.arrange(grobs = list(p0,p1,p2), nrow=1)
    dev.off()

  } 
  dbDisconnect(db_conn)
  
  }
stop('FINISHED')
