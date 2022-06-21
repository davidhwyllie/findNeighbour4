# methodscompare methods of clustering of transformed_coordinates
# resulting from pca of SARS-CoV-2 genomes.

# include code from Brice Letcher, refactored to generate a 
# database backed application with a range of helpful functions
       
# consumes data in SQLite format exported by fn4pca.py, a component 
# of the findNeighbour4 server system

# note: if would be easily possible to code all this as part of a pure python
# application
# but effort would need to be put into some of the more complex graphics


library(uuid)
library(DBI)
library(dplyr)
library(ggplot2)
library(parallel)
library(cowplot)
library(viridis)
library(data.table)
library(readr)
library(lubridate)
library(ggrepel)
library(tidyr)
library(DescTools)
library(gplots)
library(stringr)
library(svglite)
library(scales)
library(gridExtra)
library(rjson)
library(rms)
library(MASS)
# functions

# ---- helper functions for visualisation ----

pc_cat_levels <- function(x){
  # for pc_cats like 1_10, 1_1, 11_2, 2_3 return in the correct order (pc first, then cat)
  # useful for assigning levels to pc_cat fields
  if (length(x)==0){
    return(x)
  }
  x <- as.character(x)
  x <- unique(x)
  x <- x[!is.na(x)]
  bits <- strsplit(x,'_',fixed=TRUE)
  separated <- data.frame(t(data.frame(bits)))
  names(separated) <- c('pc','cat')
  separated$pc_cat <- x
  separated$pc <- as.integer(separated$pc)
  separated$cat <- as.integer(separated$cat)
  reorder <- order(separated$pc,separated$cat)
  separated <- separated[reorder,]
  separated$pc_cat
}

make_output_filename <- function(outputdir, db_name, plot_name, file_type){
  # helper function making a consistent file path for writing files
  #
  # parameters:
  #  outputdir: where the file goes
  #  db_name:  used as the first part of the filename
  #  plot_name: used as the second part of the filename
  #  file_type: e.g. tiff
  #
  # returns:
  #  filename
  return(
    paste0(paste0(outputdir,'/'), paste(db_name,plot_name,sep='-'), paste0('.',file_type))
  )
}
to_tiff <- function(to_plot, outputdir, db_name, plot_name, width=1200, height=1200) {
  # outputs a ggplot p to filename=filename
  # parameters
  #  to_plot : ggplot object
  #  filename: file to which it is to be written
  
  filename <- make_output_filename(outputdir=outputdir,
                                   db_name=db_name,
                                   plot_name=plot_name,
                                   file_type='tiff')
  
  print(paste("Writing file to ",filename))
  tiff(filename,width=width,height=height, compression='lzw')
  print(to_plot)
  dev.off()

}
to_svg <- function(to_plot, outputdir, db_name, plot_name) {
  # outputs a ggplot p to filename=filename
  # parameters
  #  to_plot : ggplot object
  #  filename: file to which it is to be written
  
  filename <- make_output_filename(outputdir=outputdir,
                                   db_name=db_name,
                                   plot_name=plot_name,
                                   file_type='svg')
  
  print(paste("Writing file to ",filename))
  dpi <- 144
  svglite(filename,width=1600/dpi,height=1200/dpi)
  print(to_plot)
  dev.off()

}
# ---- queries used multiply by visualisation and analysis ----
all_pc_cats <- function(db_connection){
  cmd <- "select distinct tcc.pc_cat pc_cat
                          from 
                          transformed_coordinate_categories tcc 
                          group by tcc.pc_cat;"
  database_query(db_connection, cmd, params=NA)

}
load_trend_summary <- function(db_connection, date_end=Sys.Date()) {
  # loads an integrated dataset including trends and the number of samples
  # 
  # parameters:
  #   db_connection: a database connection
  #   date_end:      if specified, discards any clinical data with sample_date > date_end
  #                   this is typically useful for removing malformatted dates (e.g. wrong year) from 
  #                   summary analyses
  

  # set today's date as a default date_end if we're not given one.
  if (is.na(date_end)){
    date_end <- Sys.Date()
  }
  
  date_ends <- convert_date(date_end) 
  
  # this query can have to collate large amounts of data.  the result is small.  we cache the result.
  target_table_name <- paste0('trend_summary')
  target_table_name <- gsub('-','_',target_table_name)
  if (!db_has_table(db_connection, target_table_name)) {
    print("Computing trend summary as no cached result exists")
    print(" Recovering model fits, currently for all models")
      # recover details of the Poisson model fitted
    cmd = "SELECT smm.analysis_id, analysis, analysis_type, date_end, interval_analysed_days, 
           pc_cat, param, param_desc, Estimate, Estimate_CI_low, Estimate_CI_high, p_value,
           is_reference, has_ci, Estimate2NaturalSpace
           FROM statistical_model_fits smf
    inner join statistical_model_metadata smm
    on smf.analysis_id = smm.analysis_id;"
    trend <- database_query(
      db_connection = db_connection, 
      cmd = cmd)
    print(paste0("Recovered ",nrow(trend)," rows"))
  
    check_table_present(db_connection, 'transformed_coordinate_categories','Table should be present, produced by fn4pca.py')
    check_table_present(db_connection, 'statistical_model_fits','Call fit_recent_trend_wrt_most_common_pcat_all_pcs() first')
    print('running query 1: properties of transformed coordinates')
    cmd = "
              select pc_cat, count(*) n_samples, 
              min(transformed_coordinate) min_trans_coord, 
              max(transformed_coordinate) max_trans_coord, 
              avg(transformed_coordinate) mean_trans_coord
              from transformed_coordinate_categories
              group by pc_cat;"

    summary_trends <- database_query(
      db_connection = db_connection, 
      cmd = cmd)
    print('query 1 complete')
    print(paste0("Recovered ",nrow(summary_trends)," rows"))
 

    if (!db_has_table(db_connection, 'clinical_metadata')){
      summary_trends$earliest_sample_date <- NA
      summary_trends$latest_sample_date <- NA
      summary_trends$n_gt_14_days_before <- NA
      summary_trends$n_gt_30_days_before <- NA
      summary_trends$n_gt_60_days_before <- NA
      summary_trends$n_gt_90_days_before <- NA
      summary_trends$n_gt_120_days_before <- NA
      summary_trends$n_days_observed <- NA

    } else {
      print('running query 2: date range of each transformed coordinate category')
      intervals <- c(120, 90, 60, 30, 14)

      print("Computing number of pc_cats")
      all_pcs <- all_pc_cats(db_connection)

      cmd <- "select tcc.pc_cat pc_cat, 
  sum(case when cm.sample_date < ? then 1 else 0 end) n_gt_14_days_before,
  sum(case when cm.sample_date < ? then 1 else 0 end) n_gt_30_days_before,
  sum(case when cm.sample_date < ? then 1 else 0 end) n_gt_60_days_before,
  sum(case when cm.sample_date < ? then 1 else 0 end) n_gt_90_days_before,
  sum(case when cm.sample_date < ? then 1 else 0 end) n_gt_120_days_before,
  min(sample_date) earliest_date, 
  max(sample_date) latest_date,
  count(distinct sample_date) n_days_observed
                          from 
                          transformed_coordinate_categories tcc 
                          LEFT JOIN
                          clinical_metadata cm
                          on cm.sample_id = tcc.sample_id
                          group by tcc.pc_cat;"
      print("Computing numbers of samples occurring > 14,30,60,90 and 120 days ago")
  
      datelist <- list()
      for (time_interval in intervals){
        cutoff_date <- date_ends[['dt']]-time_interval
        cutoff_date_iso <- lubridate::format_ISO8601(cutoff_date,usetz=FALSE, precision='ymd')
        datelist <- append(datelist, cutoff_date_iso)
      }
      historical_counts <- database_query(db_conn, cmd, params=datelist)
      all_pcs <- merge(all_pcs, historical_counts, all.x=TRUE, by='pc_cat')
      for (i in 2:6){
      all_pcs[,i][is.na(all_pcs[,i])] <- 0
      }

      print("query 2 complete")
      print(paste0("Recovered ",nrow(summary_trends)," rows"))
 
      summary_trends <- merge(summary_trends, all_pcs, all.x=TRUE, by=c('pc_cat'))
      print(paste0("merge 1 completed, leaving ",nrow(summary_trends)," rows"))
 
    }
    
    print("Merging")
    summary_trends <- merge(summary_trends, trend, by='pc_cat')
    print(paste0("merge 2 completed, leaving ",nrow(summary_trends)," rows"))
 
    # convert the date fields to dates
    # this is required by SQLite, which doesn't have a native date format and requires that
    # dates are stored as strings.  We do this in ymd iso format, and and cast as date using as.Date
    summary_trends$earliest_date <- as.Date(summary_trends$earliest_date)
    summary_trends$latest_date <- as.Date(summary_trends$latest_date)

    # cache result
    print(paste("Storing result as ", target_table_name))

    # due to SQLite limitations, store the dates in YMD iso format
    summary_trends$earliest_date <- lubridate::format_ISO8601(summary_trends$earliest_date,usetz=FALSE, precision='ymd')
    summary_trends$latest_date <- lubridate::format_ISO8601(summary_trends$latest_date,usetz=FALSE, precision='ymd')
    table_details <-  list()
    table_details[[target_table_name]] <- list(
        "data_table"=summary_trends, 
        "appropriate_indices"=c("pc_cat")
    )
    
    add_database_tables(db_connection, table_details, overwrite=FALSE)   # should not have to overwrite - should be once only
  }

  print(paste0("Using cached result: recovering from ",target_table_name))
  sqlcmd <- paste0("select * from ",target_table_name,";")
  summary_trends <- database_query(db_connection, sqlcmd)

  # due to SQLite limitations, we store the dates in YMD iso format.  Convert back to date
  summary_trends$earliest_date <- as.Date(summary_trends$earliest_date)
  summary_trends$latest_date <- as.Date(summary_trends$latest_date)
  print(paste0("Loaded from database table; recovered ",nrow(summary_trends)," rows"))
  
  # apply p-value cutoff
  summary_trends$is_sig <- ifelse(summary_trends$p_value <= 0.01, 1, 0)
  summary_trends$pc_cat <- factor(summary_trends$pc_cat, levels=pc_cat_levels(summary_trends$pc_cat))
  summary_trends     
}
depict_current_trends <- function(
    db_connection,
    db_stem,
    PLOT_DIR,
    lineage_to_depict = NA,
    export_sequence_identifiers = TRUE,
    axis_start = '2020-06-01',
    axis_end = Sys.Date(),
    max_samples = 1e6,
    reporting_OR_cutoff = 100,
    min_estimate_IRR =1,
    min_recent_proportion = 0.6
    ){
   
      # db_connection: database connection
      # db_stem: the part of the base name of the database, an isodate
      # PLOT_DIR: a writeable directory
      # lineage_to_depict: optionally, highlights a particular lineage & its association with trending pc_cats
      # export_sequence_identifiers:  writes sequence names belonging to trending pc_cats to file.
      # max_samples : do not highlight trends if > max_samples in the category
      # axis_start, axis_end - used for depiction.
      # reporting OR cutoff: report pc_cat/ lineage associations with OR > this
  
      # return value
      print("Plotting counts per lineage")
      retVal = list()
      res <- plot_counts_per_lineage(db_connection = db_connection,
                                            target_lineage=lineage_to_depict, 
                                            date_end = db_stem,
                                            axis_end = axis_end,
                                            axis_start = axis_start
      )
      
      p0 <- res[['p']]
      
      print("Marking selected lineages")
      res <- plot_pc_cat_size_vs_change_marking_selected_lineages(
        db_connection =db_connection,
        date_end = db_stem,
        reporting_or_cutoff = reporting_OR_cutoff,
        target_lineages = c(lineage_to_depict),
        only_show_significant_trending = TRUE) 
      
      p1 <- res[['p']]
  
      top_changers <- res[['top_changers']]
      trends <- res[['trends']]
     
      interval_analysed_days <- unique(trends$interval_analysed_days)
      
      trends$pc_cat <- as.character(trends$pc_cat)
      trends$recent_proportion <- 1-(trends$n_gt_60_days_before / trends$n_samples)
      print(paste0("Found ", nrow(trends), " trends in the database."))
  
      print("Selecting trending items")
      trends <- subset(trends, 
          is_sig==1 & 
          natural_space_Estimate >= min_estimate_IRR & 
          recent_proportion >= min_recent_proportion & 
          n_samples <= max_samples & 
          pc_cat %in% top_changers)
      print(paste0("Found ", nrow(trends), " trending pc_cats."))
  
      selection_criterion <- paste0(
              "* IRR estimate > ",
              min_estimate_IRR,
              "& significant;\n",
              "* prop seen in last 60d >= ",
              min_recent_proportion,
              " of total \n* < ",
              max_samples,
            " samples totals"
            )
  
     
        
      p2 <- ggplot() + theme_void()   # nothing
      p3 <- ggplot() + theme_void()   # nothing
     
      if (nrow(trends)>0){
        print("Trending pcs found.")
        date_ends <- convert_date(db_stem)
        trends$marker_x <- date_ends[['dt']]-interval_analysed_days
        trends$marker_xend <- date_ends[['dt']]
        trends$marker_y <- 0
        trends$marker_yend <-0
  
        print(paste0("Trends highlighted = ", nrow(trends)))
        print(paste0("Highlighted pc cats: ", paste(trends$pc_cat, collapse=';')))
        print(paste0("Lineage associated pc cats: ", paste(res$associated_pc_cats, collapse=';')))
  
        # load data
        analyse_pc_cats <- trends$pc_cat
        cnts_overall <- list()
  
        print("Depicting daily count data for these pc_cats:")
        print(analyse_pc_cats)
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
            db_connection = db_connection, 
            cmd = cmd) # recover everything, irrespective of p-value
        assoc_long <- subset(assoc_long, pc_cat %in% analyse_pc_cats)
        assoc_long$pc_cat <- factor(assoc_long$pc_cat, levels = pc_cat_levels(assoc_long$pc_cat))
        assoc_long$log_or <- log10(assoc_long$truncated_odds_ratio)
  
        assoc_long$proportion_feature_explained <- assoc_long$a/ assoc_long$a_c
        assoc_long$is_lineage <- ifelse(grepl("lineage:", assoc_long$feature), "Lineage", "Variant")
        assoc_long$is_lineage <- ifelse(grepl("mutations:", assoc_long$feature), "Mutation", assoc_long$is_lineage)
        assoc_long$feature <- gsub("lineage:","", assoc_long$feature, fixed=TRUE)
        assoc_long <- subset(assoc_long, !feature == "NA")
  
        # recover data about counts
        cnt_result <- count_per_pc_cats_in_time_interval(
            db_connection = db_connection,
            these_pc_cats = analyse_pc_cats,
            date_end = date_ends[['dt']], 
            interval_analysed= 3 * interval_analysed_days
          )
  
        ## if there are counts
        if (cnt_result$success){
          # there are trending pc_cats
          cnts <- cnt_result$count_df
          cnts$pc_cat <- factor(cnts$pc_cat, levels = pc_cat_levels(cnts$pc_cat))
          print('plotting counts')
          
          p <- ggplot(cnts, aes(x=sample_date, y= n+1))+ theme_classic()
          p <- p + geom_col()
        
          p <- p + geom_hline(colour='black', lty=2,yintercept = c(6,16,26))
          p <- p + scale_x_date(limits = c(date_ends[['dt']]-4*interval_analysed_days, date_ends[['dt']])) 
          p <- p + scale_y_log10()
          p <- p + geom_vline(colour='black', lty=1, xintercept = c( date_ends[['dt']]-4*interval_analysed_days))
          p <- p + geom_vline(colour='red', lty=1, xintercept = c(date_ends[['dt']]))
        
          p <- p + geom_hline(colour='black', lty=1,yintercept = 1)
          p <- p + facet_wrap(~pc_cat)
         
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
          if (is.na(lineage_to_depict)){
            p2_title <- paste0("Highlighted count profiles for:\n",
            selection_criterion,
            "\nBar: last ",
            interval_analysed_days,
            " days. \nHorizontal dashes: 5, 15, 25 cases/day")
  
          } else {
            p2_title <- paste0("Highlighted count profiles for:\n",
            selection_criterion,
            "\nBar: last ",
            interval_analysed_days,
            " days. Red bar = assoc with ",
            lineage_to_depict,
            "\n",
            "Horizontal dashes: 5, 15, 25 cases/day")
          }
          p2 <- p + ggtitle(
            p2_title
          )
          
          p2 <- p2 + geom_col(colour = 'black')
          p2 <- p2 + geom_segment(colour='gray', size= 3, lty=1, aes(x = marker_x, y=marker_y, xend=marker_xend, yend=marker_yend), data = trends)
          
          
          # and depict associations
          # check associations
          depiction_or_cutoff <- 1e3
  
          print("Checking associations with lineages")
          check_table_present(db_connection, 'feature_associations','run make_contingency_tables() first')
          depicted_associations <- subset(assoc_long,  truncated_odds_ratio > depiction_or_cutoff)
          
          trend_subset <- subset(trends, pc_cat %in% unique(depicted_associations$pc_cat))
          trend_subset <- subset(trend_subset, lineage == lineage_to_depict)
          if (nrow(trend_subset)>0){
            p2 <- p2 + geom_segment(colour='red', size= 3, lty=1, aes(x = marker_x, y=marker_y, xend=marker_xend, yend=marker_yend), data = trend_subset)
          }
  
          plotted <- subset(assoc_long, is_lineage %in% c('Lineage') )
          depicted_associations <- subset(depicted_associations, feature %in% unique(plotted$feature))
  
          trends$pc_cat <- as.character(trends$pc_cat)
  
          sig_assoc <- as.data.table(subset(assoc_long, truncated_odds_ratio > depiction_or_cutoff))
          collapse <- function(x){
            paste(x,collapse=';')
          }
          collapsen <- function(x){
            paste(x,collapse='\n')
          }
          
         
          if (nrow(sig_assoc)>0){
            label_df <- reshape2::dcast(pc_cat~., value.var='feature', fun.aggregate=collapsen, data =sig_assoc)
            label_df_n <- reshape2::dcast(pc_cat~., value.var='feature', fun.aggregate=length, data =sig_assoc)
          
            names(label_df)[2] <- 'associated_lineages'
            names(label_df_n)[2] <- 'n_pc_cats'
            label_df <- merge(label_df, label_df_n, by = 'pc_cat')
            label_df$pc_cat <- as.character(label_df$pc_cat)
            label_df$x <- date_ends[['dt']]-2*interval_analysed_days
            label_df$y <- 1000
  
            p2 <- p2 + geom_text(aes(x=x,y=y, label=associated_lineages), size=2, data=subset(label_df, n_pc_cats==1))
            output_summary <- merge(trends, label_df, by='pc_cat', all.x=TRUE)
            output_summary$associated_lineages[is.na(output_summary$associated_lineages)] <- 'No significant association found'  
          } else {
            output_summary <- trends
          }
         
          if (nrow(plotted)>0){
            plotted$or_hi <- ifelse(plotted$truncated_odds_ratio > depiction_or_cutoff, 1, 0)
            lineages_with_high_or <- reshape2::dcast(feature ~ ., value.var='or_hi', fun.aggregate=max, data=plotted)
            names(lineages_with_high_or)[2] <- 'n'
            lineages_with_high_or <- subset(lineages_with_high_or, n>0)
            plotted_no_empty_cols <- subset(plotted, feature %in% lineages_with_high_or$feature)
            
            
            p3 <- ggplot(plotted_no_empty_cols, aes(x=pc_cat, y=feature))
            p3 <- p3 + geom_blank()   # fill the grid
           
            p3 <- p3 + geom_tile(aes(fill=log_or), data=subset(plotted_no_empty_cols, truncated_odds_ratio > depiction_or_cutoff))
            
            p3 <- p3 + geom_point(aes(size=proportion_feature_explained), colour='yellow', data= depicted_associations)
            p3 <- p3 + ggtitle("Sequence patterns\nassociations with lineage")
            p3 <- p3 + scale_x_discrete("Sequence pattern (pc_cat)")
            p3 <- p3 + scale_y_discrete("Lineage or mutation")
            p3 <- p3 + theme(axis.text.x=element_text(angle=90))
            p3 <- p3 + scale_fill_continuous(name = "Log OR favouring\nlineage given sequence\npattern")
            p3 <- p3 + theme(legend.position="top", legend.direction = 'vertical', legend.box='horizontal')
            p3 <- p3 + scale_size_continuous("Proportion of\nlineage comprising\nsequence pattern")
            #p3 <- p3 + facet_grid(.~is_lineage, scales='free')
            
          # dump trending sample lists
          if (export_sequence_identifiers){
            date_ends <- convert_date(db_stem)
            cutoff_date <- as.character(date_ends[['dt']]-interval_analysed_days)
  
            print(paste0("Dumping contents of trending pc_cats to file, with sample_dates after ",cutoff_date))
            samples <- list()
            for (this_pc_cat in unique(plotted$pc_cat)){
              sqlcmd = "select distinct pc_cat, cm.sample_id, cm.sample_date from transformed_coordinate_categories tcc
                        inner join clinical_metadata cm 
                        on cm.sample_id = tcc.sample_id
                        inner join sequence_metadata sm
                        on sm.sample_id = tcc.sample_id
                        where tcc.pc_cat = ? and
                        cm.sample_date > ?; "
  
              samples[[this_pc_cat]] <- database_query(db_connection, sqlcmd, params=list(this_pc_cat, cutoff_date))
            }
            all_samples <- do.call(rbind, samples)
            retVal[['total_sequences']] <- length(unique(all_samples$sample_id))
            retVal[['sequences']] <- unique(all_samples$sample_id)
            retVal[['trending_lineages']] <- unique(plotted$feature)
            retVal[['trending_pc_cats']]  <- unique(plotted$pc_cat)
            outfile <- paste0(PLOT_DIR,'/',db_stem,"_seqlist.tsv")
            retVal[['sequences_list']] <- outfile
            write.table(all_samples, outfile, sep='\t', row.names=FALSE)
 
          }
            
        } 
         
        }
      }
  
      
    print("Starting output")

    lineage_filename <- gsub('.','_',lineage_to_depict, fixed= TRUE)
    outfile <- paste0(PLOT_DIR,'/',db_stem,'_pc_cats_highlighting_', lineage_filename,".tiff")

    retVal[['depiction_tiff']] <- outfile
    tiff(outfile,width=1600,height=1000, compression='lzw')
    gridExtra::grid.arrange(grobs = list(p0,p1,p2,p3), nrow=1, ncol=4, widths=c(1.5,1.5,2,3) )
    dev.off()
    print(paste0("Plot rendered to ",outfile))
  
    #outfile <- paste0(PLOT_DIR,'/',db_stem,"_data.tsv")
    #retVal[['depiction_data']] <- outfile
    #write.table(output_summary, outfile, sep='\t', row.names=FALSE)
    #print(paste0("Output summary written to ",outfile))
    
  retVal
}
count_per_pc_in_time_interval <- function(db_connection, 
                                          this_pc,
                                          date_end=Sys.Date(),
                                          interval_analysed=30){
  # construct a data.table comprising all date/category combinations
  # for the period interval_analysed days prior to date_end
  #
  # useful for regression modelling.
  #
  # parameters:
  #   db_connection: the database connection
  #   this_pc : the PC to analyse
  #   date_end: the last day analysed
  #   interval_analysed:  the number of days prior to date_end analysed
  # 
  
  # dates may be passed in iso format as strings.

  date_ends <- convert_date(date_end) 

  date_start <- date_ends[['dt']] - interval_analysed
  date_sequence <- seq(date_start, to = date_ends[['dt']], by =1)
  date_t <- as.numeric(as.integer(date_sequence-date_ends[['dt']]))/30 #per month
  last_week <- ifelse(date_t<= -7, 1, 0)
  dow_sequence <- lubridate::wday(date_sequence)
  time_series_meta <- data.table(
    sample_date = lubridate::format_ISO8601(date_sequence,usetz=FALSE, precision='ymd'),
    date_t = date_t,
    last_week = last_week,
    sample_dow = factor(dow_sequence)
  )
  
  # recover counts to model
  cmd <- "select ec.cat, cm.sample_date, count(*) n
  from transformed_coordinate_categories ec 
  INNER JOIN
  clinical_metadata cm
  ON cm.sample_id = ec.sample_id
  where cm.sample_date >= ? and cm.sample_date <= ? and pc= ?
  group by ec.cat, cm.sample_date;"
  
  params <- list(
    lubridate::format_ISO8601(date_start, usetz=FALSE, precision='ymd'),
    date_ends[['iso']],
    this_pc)

  # generate data frame for linear model, with zeros as appropriate
  count_df <- database_query(
    db_connection = db_connection, 
    cmd = cmd, 
    params = params)

  # check how many cats there are.  If there is only one, we give up
  # as we cannot fit relevant models.
  if (length(unique(count_df$cat)) < 2) {
    return(list('fitted' = FALSE, 'count_df' = count_df, 'reason' = "One or no categories"))
  }
  cats <- data.table(cat= count_df$cat)
  cats$unity <- 1
  time_series_meta$unity <- 1
  time_series_meta <- merge(time_series_meta, cats, allow.cartesian=TRUE, by='unity') # cross join
  
  count_df <- merge(time_series_meta, count_df, by = c('cat','sample_date'), all.x=TRUE)
  count_df$n[is.na(count_df$n)] <- 0  # missing date in the SQL query means count 0
  
  # convert cat to a factor, ordered such the the most common comes first
  count_df$pc_cat <- paste(this_pc, count_df$cat, sep = '_')
  total_counts <- count_df %>% group_by(pc_cat) %>% summarise(total_n=sum(n)) %>% arrange(total_n)
  count_df$cat <- factor(count_df$pc_cat,levels=rev(total_counts$pc_cat))
  return(list('fitted' = TRUE, 'count_df' = count_df, 'reason' = "2 or more categories"))

} 

count_per_pc_cats_in_time_interval <- function(db_connection, 
                                          these_pc_cats,
                                          date_end=Sys.Date(),
                                          interval_analysed=30){
  # construct a data.table comprising all date/category combinations
  # for the period interval_analysed days prior to date_end
  #
  # useful for regression modelling.
  #
  # parameters:
  #   db_connection: the database connection
  #   these_pc_cats : the pc_cats to analyse
  #   date_end: the last day analysed
  #   interval_analysed:  the number of days prior to date_end analysed
  # 
  ### NOT TESTEDF
  ### TODO Add pc_cat into transformed_coordinate_categories

  # dates may be passed in iso format as strings.

  print("Recovering counts by pc_cats")
  print("analysing")
  
  if (length(these_pc_cats)==0){
    return(list('success' = FALSE, 'count_df' = NA, 'reason' = "No pc_cats provided"))
  }

  date_ends <- convert_date(date_end) 
  date_start <- date_ends[['dt']] - interval_analysed
  date_sequence <- seq(date_start, to = date_ends[['dt']], by =1)
  date_t <- as.numeric(as.integer(date_sequence-date_ends[['dt']]))/30 #per month
  last_week <- ifelse(date_t<= -7, 1, 0)
  dow_sequence <- lubridate::wday(date_sequence)
  time_series_meta <- data.table(
    sample_date = lubridate::format_ISO8601(date_sequence,usetz=FALSE, precision='ymd'),
    date_t = date_t,
    last_week = last_week,
    sample_dow = factor(dow_sequence)
  )
  
  params_passed <- paste0(rep('?', length(these_pc_cats)),collapse=',')
  # recover counts to model
  cmd <- paste0("select ec.pc_cat, cm.sample_date, count(*) n
  from transformed_coordinate_categories ec 
  INNER JOIN
  clinical_metadata cm
  ON cm.sample_id = ec.sample_id
  where cm.sample_date >= ? and cm.sample_date < ? and pc_cat in (",
  params_passed,
  ") group by ec.pc_cat, cm.sample_date;"
  )
 
  params <- list(
      lubridate::format_ISO8601(date_start, usetz=FALSE, precision='ymd'),
      date_ends[['iso']]
  )
  element_number <- 2
  for (this_pc_cat in these_pc_cats){
    element_number <- element_number + 1
    params[[element_number]] <- this_pc_cat
  }

  # load count data
  count_df <- database_query(
    db_connection = db_connection, 
    cmd = cmd, 
    params = params)

  cats <- data.table(pc_cat=these_pc_cats)
  cats$unity <- 1
  time_series_meta$unity <- 1
  time_series_meta <- merge(time_series_meta, cats, allow.cartesian=TRUE, by='unity') # cross join
  
  count_df <- merge(time_series_meta, count_df, by = c('pc_cat','sample_date'), all.x=TRUE)
  count_df$n[is.na(count_df$n)] <- 0  # missing date in the SQL query means count 0
  count_df$sample_date <-  as.Date(count_df$sample_date)

  return(list('success' = TRUE, 'count_df' = count_df, 'reason' = "Success"))

} 

count_per_all_pc_cats_in_time_interval <- function(db_connection, 
                                          date_end=Sys.Date(),
                                          interval_analysed=30){
  # construct a data.table comprising all date/category combinations
  # for the period interval_analysed days prior to date_end
  #
  # useful for regression modelling.
  #
  # parameters:
  #   db_connection: the database connection
  #   date_end: the last day analysed
  #   interval_analysed:  the number of days prior to date_end analysed
  # 

  # dates may be passed in iso format as strings.

  print("Recovering counts by all pc_cats.")

  cats <- data.table(pc_cat=all_pc_cats(db_connection)$pc_cat)
  cats$unity <- 1

  date_ends <- convert_date(date_end) 
  date_start <- date_ends[['dt']] - interval_analysed
  date_sequence <- seq(date_start, to = date_ends[['dt']], by =1)
  date_t <- as.numeric(as.integer(date_sequence-date_ends[['dt']]))/30 #per month
  last_week <- ifelse(date_t<= -7, 1, 0)
  date_cat <- as.integer(date_sequence-date_ends[['dt']])
  dow_sequence <- lubridate::wday(date_sequence)
  time_series_meta <- data.table(
    sample_date = lubridate::format_ISO8601(date_sequence,usetz=FALSE, precision='ymd'),
    date_t = date_t,
    date_cat = factor(date_cat),
    last_week = last_week,
    sample_dow = factor(dow_sequence)
  )

 # recover counts to model
  cmd <- "select ec.pc_cat, cm.sample_date, count(*) n
  from transformed_coordinate_categories ec 
  INNER JOIN
  clinical_metadata cm
  ON cm.sample_id = ec.sample_id
  where cm.sample_date >= ? and cm.sample_date < ? 
  group by ec.pc_cat, cm.sample_date;"
   
  params <- list(
      lubridate::format_ISO8601(date_start, usetz=FALSE, precision='ymd'),
      date_ends[['iso']]
  )
  
  # load count data
  count_df <- database_query(
    db_connection = db_connection, 
    cmd = cmd, 
    params = params)

  time_series_meta$unity <- 1
  time_series_meta <- merge(time_series_meta, cats, allow.cartesian=TRUE, by='unity') # cross join
  
  count_df <- merge(time_series_meta, count_df, by = c('pc_cat','sample_date'), all.x=TRUE)
  count_df$n[is.na(count_df$n)] <- 0  # missing date in the SQL query means count 0
  count_df$sample_date <-  as.Date(count_df$sample_date)

  return(list('success' = TRUE, 'count_df' = count_df, 'reason' = "Success"))

} 

categorise_sample_dates <- function(summary_trends, date_end=Sys.Date()){
  # categorises the latest_date field relative to date_end; useful for highlighting recently
  # emerged samples
  # dates may be passed in iso format as strings.
  date_ends <- convert_date(date_end) 
  
  summary_trends$first_seen_days_ago <- as.integer(date_ends[['dt']] - summary_trends$earliest_date)
  summary_trends$first_seen_days_ago_cat <- cut(summary_trends$first_seen_days_ago, breaks=c(0, 7, 14, 28, 10000), include.lowest=TRUE)
  levels(summary_trends$first_seen_days_ago_cat) <- c('0 to 7 days before', '8 to 14 days before', '15 to 28 days ago', '>28 days before')
  summary_trends
}
uniquely_linked_to <- function(db_connection, 
                                 target_lineage,
                                 reporting_or_cutoff,
                                 date_end){
    # generates a list of pc_cats and whether they are linked to target_strain
    # parameters:
    #   db_connection: the connection to the database
    #   target_lineage: the lineage to look for unique links to
    #   reporting_or_cutoff:  report lineage - pc_cat association > reporting_or_cutoff
    #   date_end: the date of reporting.  Will not report associations to clones isolated in the future.
    #                                     If NA, ignores this setting.


    # determine whether there we are being asked to associated with a lineage which does not exist
    # before date_end.

    search_for_lineage <- TRUE

    if (!is.na(date_end)){
      date_ends <- convert_date(date_end)
      search_lineage <- paste0(target_lineage,'.%')
      cmd = "select count(*) n from clinical_metadata cm 
              INNER JOIN
              sequence_metadata sm
              WHERE
              cm.sample_id = sm.sample_id AND
              sm.variable = 'lineage' AND
              (sm.value LIKE ? or sm.value = ?) AND
              sample_date < ?;"

      count_previous_lineage <- database_query(
          db_connection = db_connection, 
          cmd = cmd, 
          params = list(search_lineage, target_lineage, date_ends[['iso']])
          )

      n_samples <- count_previous_lineage$n[1]
      if (n_samples == 0){
        # don't search for anything
        search_for_lineage <- FALSE
      }
    }
    
    #print(paste("Searching for associations of ",target_lineage))
    cmd <- "select substr(feature, 9) lineage, 
    pc_cat from 
    feature_associations
    where substr(feature,1,8)='lineage:' AND
    truncated_odds_ratio >= ?;"
    sa_map <- database_query(db_connection,
                                          cmd,
                                          params=list(reporting_or_cutoff))
    ## restrict to target strains.
    # contains a list of pcs which are targeted
    target_pattern1 <- paste0('^', target_lineage,"\\..*$")
    target_pattern2 <- paste0('^', target_lineage,"$")
    target_patterns <- c(target_pattern1,
                        target_pattern2)
  
    select_rows <- c()
    for (pattern in target_patterns){
      new_selection <- which(grepl(pattern, sa_map$lineage))
      select_rows <- union(select_rows, new_selection)
    }
  
    # mark rows which match
    sa_map$off_target <- 1
    sa_map$off_target[select_rows] <- 0
    
    # work out which PCs are ** only ** on-target
    sa_map_on_target <- reshape2::dcast(sa_map, pc_cat ~ ., value.var='off_target', fun.aggregate=sum)
    names(sa_map_on_target)[2] <- 'off_target'    # if a PC_cat maps to something other than the target, off_target > 0
    sa_map_on_target <- subset(sa_map_on_target, off_target == 0)
    
    # if we're told not to search, we just delete the records
    if (!search_for_lineage){
      sa_map_on_target <- sa_map_on_target[0,]
    }

    if (nrow(sa_map_on_target)==0){
      sa_map_on_target = data.frame(pc_cat=character(), off_target=integer(), lineage=character())
    } else {
      sa_map_on_target$lineage <- target_lineage
    }

    sa_map_on_target
    
}

# ---- specific visualisations ----
plot_pc_cat_size_vs_change_marking_selected_lineages <- function(
                                        db_connection,
                                        date_end,
                                        target_lineages = c(),
                                        reporting_or_cutoff = 1,
                                        mark_top = 40,
                                        only_show_significant_trending = FALSE){  
  # plots pc_categories relative to the most common category
  # over the time period analysed
  
  # parameters:
  #   database connection, as returned (for example) by get_db_connection
  #   date_end:  the end of the analysis period
  #   target_lineages: lineages to mark on the chart.  ex: c('B.1.1.380','B.1.1.97','B.1.9').  If c() or NA, no lineages are marked
  #   only_show_significant_trending : whether to only show pc_cats with significant trends
  # returns
  #   ggplot object
  # dates may be passed in iso format as strings.

  # SQL to recover counts

  print(paste("Marking trending lineages matching ",paste(target_lineages, collapse=';')))
  date_ends <- convert_date(date_end) 
  
  sa_maps = list()

  if (!is.na(target_lineages)){
    for (target_lineage in target_lineages){
      if (!is.na(target_lineage)){
        sa_maps[[target_lineage]] <- uniquely_linked_to(db_connection, 
                                      target_lineage,
                                      reporting_or_cutoff,
                                      date_end)
      }
    }
    sa_map <- unique(do.call(rbind, sa_maps))
    n_specific_pcs <- nrow(sa_map)
    #label_df <- data.frame(x=1,
    #                      y=110, label=paste0("Specific pc_cat=", n_specific_pcs))
    p1_title <- paste0("Change relative to total\nLineage specific pc_cats = ", n_specific_pcs)

  } else {  
    n_specific_pcs <- 0
    p1_title <- "Change relative to total"
  }
  # load trend summary 
  summary_trends <- load_trend_summary(
    db_connection = db_connection, 
    date_end = date_ends[['dt']])
  
  # add a date categorisation field
  summary_trends <- categorise_sample_dates(summary_trends, date_end = date_ends[['dt']])
  
  # there may be multiple analyses here.
  # these analyses are identified by the analysis_id, but the analysis_id is not human readable.
  
  # in this depiction, we don't show any reference categories
  summary_trends <- subset(summary_trends, is_reference == 0)

  # we generate incidence rate ratios.
  summary_trends$natural_space_Estimate <- ifelse(summary_trends$Estimate2NaturalSpace=='exp', 
                                                exp(summary_trends$Estimate), 
                                                summary_trends$Estimate
                                                )
  summary_trends$natural_space_Estimate_CI_high <- ifelse(summary_trends$Estimate2NaturalSpace=='exp', 
                                                exp(summary_trends$Estimate_CI_high), 
                                                summary_trends$Estimate_CI_high
                                                )
  summary_trends$natural_space_Estimate_CI_low <-  ifelse(summary_trends$Estimate2NaturalSpace=='exp',
                                                 exp(summary_trends$Estimate_CI_low), 
                                                 summary_trends$Estimate_CI_low
                                                )

  # generate an variable for plotting bounded by (0.1, 10)
  # in the models currently used, these are rate ratios
  summary_trends$plot_natural_space_Estimate <- summary_trends$natural_space_Estimate

  summary_trends$plot_natural_space_Estimate <- ifelse(
                                                  summary_trends$natural_space_Estimate> 3, 
                                                  3,
                                                  summary_trends$natural_space_Estimate)
  summary_trends$plot_natural_space_Estimate <- ifelse(
                                                  summary_trends$natural_space_Estimate< 0.3, 
                                                  0.3,
                                                  summary_trends$natural_space_Estimate)
  
  # incorporate any lineages linked in summary_trends
  if (!is.na(target_lineages)){
  summary_trends <- merge(summary_trends, sa_map, all.x=TRUE, by = 'pc_cat')
  } else {
    summary_trends$lineage <- NA
    summary_trends$off_target <- NA
  }

  if (only_show_significant_trending){
    summary_trends <- subset(summary_trends,  is_sig == 1)
  } 

  # bespoke for this analysis: need to edit if different statistical models are added
  summary_trends <- subset(summary_trends, param_desc == "Incidence rate ratio")
  summary_trends <- summary_trends[order(-summary_trends$natural_space_Estimate),]

  y_axis_breaks <- c(3,2,1.5,1,0.7, 0.5,0.3)
  p1<-ggplot(summary_trends,aes(x=n_samples,y=plot_natural_space_Estimate)) + 
    geom_point(alpha = 0.3, aes()) + theme(legend.position = "none") 

  p1 <- p1 +
    geom_hline(yintercept =1, lty =1) +
    geom_rug(alpha = 0.05)
  p1 <- p1 + scale_x_continuous("Total sample counts in category",trans='log10', limits=c(1,1e6)) 

  p1 <- p1 +
    scale_y_log10("30-day growth rate of samples, relative to total\n [Rate ratios > 10 plotted at 10]",
                  breaks=y_axis_breaks, labels=y_axis_breaks, limits = c(0.2, 3.3)) 
  
  p1 <- p1 + ggtitle(p1_title)
  marked <- subset(summary_trends, Estimate > 0 & !is.na(lineage))
  top_changers <- subset(summary_trends, Estimate>0)
  print(top_changers)
  if (mark_top > 0){
    if (mark_top>nrow(top_changers)){
      mark_top <- nrow(top_changers)
    }
    
    top_changers <- top_changers[1:mark_top,]
 
    p1 <- p1 + geom_text_repel(size=3,  aes(label = pc_cat),
                               data=subset(top_changers, !pc_cat %in% marked$pc_cat)) 
    
  }
  

  if (nrow(marked)>0){
    
    p1 <- p1 +  geom_point(size=2, aes(colour = lineage),
                           data=marked)
    
    p1 <- p1 + geom_point(size=4, shape=1, aes(colour = lineage),
                          data=marked)  
    p1 <- p1 + geom_text_repel(size=4, aes(label = pc_cat, colour = lineage),
                               data=marked) 
    
  }
  associated_pc_cats <- unique(summary_trends$pc_cat[!is.na(summary_trends$lineage)])
  summary_trends$pc_cat <- factor(summary_trends$pc_cat, levels=pc_cat_levels(summary_trends$pc_cat))
  res = list(trends = summary_trends, top_changers = as.character(top_changers$pc_cat), p= p1, associated_pc_cats = associated_pc_cats)
}
plot_counts_per_lineage <- function(db_connection,
                                        target_lineage, 
                                        date_end = Sys.Date(),
                                        axis_end = Sys.Date(),
                                        axis_start = "2020-06-01") {

  # plot counts per lineage over time.
  #
  # parameters:
  #   db_connection: DBI connection
  #   lineage:  lineage to plot, e.g. B.1.1.7. If NA, no lineage information is plotted.
  #   date_end: only report counts up this point.  Either ISO date string or Date()
  #   axis_start: LHS of the plot.  Ignore counts before this.  Either ISO date string or Date()
  #   axis_end: RHS of the plot.  Ignore counts after this.  Either ISO date string or Date()
  #
  # returns:
  #   ggplot object.

  # convert dates
  date_ends <- convert_date(date_end) 
  axis_ends <- convert_date(axis_end)
  axis_starts <- convert_date(axis_start)

# SQL to recover counts
  if (is.na(target_lineage)){
    target_lineage <- '-'
  
 } else {
    y_axis_title <- "Number of samples sequenced per day"

  }
  search_lineage <- paste0(target_lineage,'.%')     
  cmd = "
  select sample_date, value lineage_variant, count(*) n from clinical_metadata cm 
INNER JOIN
sequence_metadata sm
WHERE
cm.sample_id = sm.sample_id AND
sm.variable = 'lineage' AND
(sm.value = ? or sm.value LIKE ?)
and sample_date <= ?
and sample_date <= ?
and sample_date >= ?
group by sample_date, value
order by sample_date;"

  cnts <- database_query(db_connection,
                        cmd = cmd,
                        params = list(target_lineage,
                                      search_lineage,
                                      date_ends[['iso']],
                                      axis_ends[['iso']],
                                      axis_starts[['iso']]
                                       )
  )
  cnts$sample_date <- as.Date(cnts$sample_date)
  total_cnts <- sum(cnts$n)

  cmd = "
  select sample_date, count(*) n from 
  clinical_metadata cm 
  where
sample_date <= ?
and sample_date <= ?
and sample_date >= ?
group by sample_date
order by sample_date;"

  df_total_cnts <- database_query(db_connection,
                        cmd = cmd,
                        params = list(
                                      date_ends[['iso']],
                                      axis_ends[['iso']],
                                      axis_starts[['iso']]
                                       )
  )
  
  df_total_cnts$sample_date <- as.Date(df_total_cnts$sample_date)
  total_sequenced <- sum(df_total_cnts$n)
  label_df <- data.frame(
    x=axis_starts[['dt']],
    y=5000, label=paste("Gray bars are total sequenced numbers"))


  # make labels appropriate to whether we have a lineage to mark
  if ( !target_lineage == '-'){
    y_axis_title <-  paste0("Number samples sequenced per day (gray = total, coloured ", target_lineage,")")
    title_text <- paste0("Date: ", date_end,": ",target_lineage, "  n=", total_cnts)
  } else {
    y_axis_title <- "Number of samples sequenced per day"
    title_text <- ''
  }
  #title_text <- paste0(title_text, "\n", paste("Total sequences =", total_sequenced))
  breaks = c(1,3,10,30,100,300,1000,2000,3000,5000,15000,50000)
  p <- ggplot(cnts, aes(x=sample_date, y=n))
  p <- p + geom_col(data=df_total_cnts, colour='lightgray')
  p <- p + geom_col(aes(fill = factor(lineage_variant)), position='stack')
  p <- p + scale_x_date("Date of collection", 
  limits = c(axis_starts[['dt']], axis_ends[['dt']]),
        labels=date_format("%m-%Y"))
  p <- p + scale_fill_discrete("Lineage") 
  p <- p + scale_y_sqrt(name = y_axis_title, breaks=breaks, labels=breaks, limits=c(0,5000))
  p <- p + geom_vline(xintercept = date_ends[['dt']], colour='red')
  #p <- p + geom_label(aes(x=x,y=y,label=label), size= 3, data=label_df)
  p <- p + theme(legend.position = "none")
  p <- p + ggtitle(title_text)
  list(cnts = cnts, p = p)
}
plot_single_pc_format_1 <- function(db_connection, pc) {
  # plots a histogram illustrating transformed_coordinates(x) and the categories within it
  #
  # parameters:
  #   database connection, as returned (for example) by get_db_connection
  #   pc: the pc to plot
  #
  # returns
  #   ggplot object
  
  check_table_present(db_connection, 'transformed_coordinate_categories','regenerate sqlitedb using recent pca code which categorises transformed_coordinates')

  print("plot_single_pc_format_1: Loading data #1")
  cmd1 = "select transformed_coordinate ev, 
          cat from transformed_coordinate_categories 
          where pc = ?"
  df <- database_query(
    db_connection = db_connection, 
    cmd = cmd1, 
    params = list(pc))
  print("plot_single_pc_format_1: Loading data #2")
  cmd2 = "select pc, cat, count(*) n_samples, min(transformed_coordinate) min_trans_coord, max(transformed_coordinate) max_trans_coord, avg(transformed_coordinate) mean_trans_coord
        from transformed_coordinate_categories
        where pc= ? group by pc, cat;"
  summary_trends <- database_query(
    db_connection = db_connection, 
    cmd = cmd2, 
    params = list(pc))
  ll <- min(summary_trends$min_trans_coord)
  ul <- max(summary_trends$min_trans_coord)
  print("Plotting")  
  df$cat <- factor(df$cat)
  p1 <- ggplot(df, aes(x=ev))
  p1 <- p1 + geom_histogram(bins=500, aes(fill=cat))  
  p1 <- p1 + scale_x_continuous("transformed_coordinate",limits=c(ll,ul))
  p1 <- p1 + scale_fill_viridis_d("transformed_coordinate") + theme(legend.position = "none")
  p1 <- p1 + ggtitle(paste0("Principal component ",pc))
  br = c(10,20,30,50,100,250,500,1000,2500,5000,10000,20000,30000,50000,100000)
  p1 <- p1 + scale_y_sqrt("Sample count", labels=br, breaks=br)

  p1 <- p1 + geom_vline(aes(xintercept=min_trans_coord),data=summary_trends,alpha= 0.1)
  p1 <- p1 + geom_label(aes(x=mean_trans_coord,label=cat, y=2000), data=summary_trends)
  p1
}
plot_single_pc_format_2 <- function(db_connection, pc, date_end) {
  # plots a histogram illustrating transformed_coordinates(x) and the categories within it,
  # over time
  #
  # parameters:
  #   database connection, as returned (for example) by get_db_connection
  #   pc: the pc to plot
  #   date_end: the last date to plot
  
  # returns
  #   ggplot object
  
  check_table_present(db_connection, 'transformed_coordinate_categories','regenerate sqlitedb using recent pca code which categorises transformed_coordinates')
  check_table_present(db_connection, 'clinical_metadata','run add_cog_metadata() or equivalent')
  
  print("plot_single_pc_format_2: Loading data")
  cmd2 = "select pc, cat, cm.sample_date, count(*) n_samples from 
  transformed_coordinate_categories tcc
  INNER JOIN sample_id sm
  ON tcc.sample_id = sm.sample_id
  INNER JOIN clinical_metadata cm 
  on cm.sample_id = sm.sample_id
  where cm.sample_date <= ?
  and pc=?
  group by pc, cat, sample_date;"
  cnts <- database_query(
    db_connection = db_connection, 
    cmd = cmd2, 
    params = list(as.character(date_end), pc))
  print("Plotting")
  
  
  cnts$sample_date <- as.Date(cnts$sample_date)

  cnts$cat <- factor(cnts$cat)
  p1 <- ggplot(cnts, aes(x=sample_date))
  p1 <- p1 + geom_col(aes(y=n_samples, fill=cat),position='stack') 
  p1 <- p1 + scale_x_date("Specimen collection date")
  #p1 <- p1 + scale_y_continuous("Number of samples", limits=c(0,5000))
  p1 <- p1 + ggtitle(paste0("Principal component ",pc))
  br <- c(10,20,30,50,100,250,500,1000,2500,5000,10000,20000,30000,50000,100000)
  p1 <- p1 + scale_y_continuous("Number of samples", breaks= br, labels=br, limits = c(0,100000))
  p1
}
plot_single_pc_format_3 <- function(db_connection, pc, date_end,time_interval){

  cmd1 = "select round(transformed_coordinate,2) transformed_coordinate, cm.sample_date, count(*) n_samples from 
  transformed_coordinate_categories tcc
  INNER JOIN sample_id sm
  ON tcc.sample_id = sm.sample_id
  INNER JOIN clinical_metadata cm 
  on cm.sample_id = sm.sample_id
  where cm.sample_date <= ?
  and pc=?
  group by round(transformed_coordinate,2), sample_date;"
    cnts <- database_query(
      db_connection = db_connection, 
      cmd = cmd1, 
      params = list(as.character(date_end), pc))
    cmd2 = "select pc, cat, count(*) n_samples, min(transformed_coordinate) min_trans_coord, max(transformed_coordinate) max_trans_coord, avg(transformed_coordinate) mean_trans_coord
        from transformed_coordinate_categories
        where pc= ? group by pc, cat;"
    summary_trends <- database_query(
      db_connection = db_connection, 
      cmd = cmd2, 
      params = list(pc))
    summary_trends$date_end <- as.Date(date_end)+7
    print("Plotting")
    
    
    cnts$sample_date <- as.Date(cnts$sample_date)
    start_date <- as.Date(date_end)-time_interval
    p1 <- ggplot(cnts, aes(y=sample_date))+theme_classic()
    p1 <- p1 + geom_tile(aes(x=transformed_coordinate, fill=log10(1+n_samples)))
    p1 <- p1 + scale_y_date("Specimen collection date")
    p1 <- p1 + scale_x_continuous("Transformed coordinate")
    p1 <- p1 + ggtitle(paste0("Principal component ",pc))
    p1 <- p1 + scale_fill_viridis("Number of\nsamples\n(log10)")
    p1 <- p1 + geom_vline(aes(xintercept=min_trans_coord),data=summary_trends,alpha= 0.1)
    p1 <- p1 + geom_hline(yintercept=as.Date(date_end), color='red')
    p1 <- p1 + geom_text(aes(x=mean_trans_coord,label=cat, y=date_end), data=summary_trends)
    
    p2 <- ggplot(subset(cnts, sample_date > start_date), aes(y=sample_date))+theme_classic()
    p2 <- p2 + geom_tile(aes(x=transformed_coordinate, fill=log10(1+n_samples)))
    p1 <- p2 + scale_y_date("Specimen collection date")
    p2 <- p2 + scale_x_continuous("Transformed coordinate")
    p2 <- p2 + ggtitle(paste0("Principal component ",pc))
    p2 <- p2 + scale_fill_viridis("Number of\nsamples\n(log10)")
    p2 <- p2 + geom_vline(aes(xintercept=min_trans_coord),data=summary_trends,alpha= 0.1)
    p2 <- p2 + geom_hline(yintercept=as.Date(date_end), color='red')
    p2 <- p2 + geom_text(aes(x=mean_trans_coord,label=cat, y=date_end), data=summary_trends)
    
    list(p1=p1, p2=p2)
}
plot_all_pc_format_2 <- function(db_connection, date_end) {
  # plots a histogram illustrating transformed_coordinates(x) and the categories within it,
  # over time
  #
  # parameters:
  #   database connection, as returned (for example) by get_db_connection
  #   date_end: plots samples received up to date_end
  # returns
  #   ggplot object
  
  check_table_present(db_connection, 'transformed_coordinate_categories','regenerate sqlitedb using recent pca code which categorises transformed_coordinates')
  check_table_present(db_connection, 'clinical_metadata','run add_cog_metadata() or equivalent')
  
  print("plot_all_pcs_format_2: Loading data (note: may be slow)")
  cmd2 = "select pc, cat, substr(cm.sample_date,1,7) ym, count(*) n_samples from 
  transformed_coordinate_categories tcc
  INNER JOIN sample_id sm
  ON tcc.sample_id = sm.sample_id
  INNER JOIN clinical_metadata cm 
  on cm.sample_id = sm.sample_id
  where cm.sample_date <= ?
  group by pc, cat, substr(cm.sample_date,1,7);"
  cnts <- database_query(
    db_connection = db_connection, 
    cmd = cmd2, 
    params = list(date_end))
  print("Plotting")
  cnts$cat <- factor(cnts$cat)
  p1 <- ggplot(cnts, aes(x=ym))
  p1 <- p1 + geom_col(aes(y=n_samples, fill=cat),position='stack') 
  p1 <- p1 + scale_x_discrete("Specimen collection date (months)")
  p1 <- p1 + theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
  br <- c(100,1000,10000,30000,60000, 100000)
  p1 <- p1 + scale_y_sqrt("Number of samples", breaks=br, labels=br)
  p1 <- p1 + facet_wrap(~pc)
  p1
}
plot_pc_format_3 <- function(db_connection, pc, date_end, time_interval) {
  # plots counts per category in recent time: interval_days from date_end
  #
  # parameters:
  #   database connection, as returned (for example) by get_db_connection
  #   pc : the pc whose components to plot
  #   date_end: plots samples received up to date_end
  #   interval: plots most recent interval days samples
  # returns
  #   ggplot object
  
  check_table_present(db_connection, 'transformed_coordinate_categories','regenerate sqlitedb using recent pca code which categorises transformed_coordinates')
  check_table_present(db_connection, 'clinical_metadata','run add_cog_metadata() or equivalent')
  cutoff_date <- as.character(as.Date(date_end)-time_interval)
  print(paste(cutoff_date, date_end, pc))
  cmd3 = "
  select cat, 
  cm.sample_date, 
  count(*) n_samples from 
  transformed_coordinate_categories tcc
  INNER JOIN sample_id sm
  ON tcc.sample_id = sm.sample_id
  INNER JOIN clinical_metadata cm 
  on cm.sample_id = sm.sample_id
  where pc = ?
  and cm.sample_date >= ?
  and cm.sample_date <= ?
  group by cat, cm.sample_date;"
  
  #cm.sample_date <= ?
  #  and cm.sample_date >= ?
  cnts <- database_query(
    db_connection = db_connection, 
    cmd = cmd3, 
    params = list(pc, cutoff_date, as.character(date_end)))
  print("Plotting")
  cnts$cat <- factor(cnts$cat)
  cnts$sample_date <- as.Date(cnts$sample_date)
  p3 <- ggplot(cnts, aes(x=sample_date))
  p3 <- p3 + geom_col(aes(y=n_samples),position='stack') 
  p3 <- p3 + scale_x_date("Specimen collection date")
  br <- c(10,100,300,1000)
  p3 <- p3 + scale_y_sqrt("Number of samples", breaks=br, labels=br)
  p3 <- p3 + facet_wrap(~cat)
  p3 <- p3 + ggtitle(paste("PC", pc))
  p3
  
}
association_heatmap <- function(db_connection,
                                outputdir,
                                db_name,
                                plot_name,
                                pdf_size_in_inches){
  # displays a heatmap illustrating the associations of known features with PCs
  #
  # parameters:
  #  db_connection: connection to the database
  #  outputdir: where the file goes
  #  db_name:  used as the first part of the filename
  #  plot_name: used as the second part of the filename
  #
  # returns:
  #  filename written
  
  check_table_present(db_connection, 'feature_associations','run make_contingency_tables() first')
  
  # recover data
  cmd <- "
select feature, pc_cat, p_value, truncated_odds_ratio 
from feature_associations
where p_value < ?
order by feature, pc_cat;"
  
  # get long dataset
  assoc_long <- database_query(
    db_connection = db_connection, 
    cmd = cmd, 
    params = list(1)) # recover everything, irrespective of p-value
  assoc_long$log_or <- log10(assoc_long$truncated_odds_ratio)

  # build matrix
  m <- reshape2::dcast(pc_cat~feature, 
                       data=assoc_long, 
                       value.var='log_or', 
                       fun.aggregate=sum)
  m_num <- m[,2:ncol(m)]
  mat <- as.matrix(m_num)
  rownames(mat) <- m$pc_cat
  
  # construct and output heatmap
  filename <- make_output_filename(outputdir=outputdir,
                                   db_name=db_name,
                                   plot_name=plot_name,
                                   file_type='pdf')
  
  print(paste("Writing file to ",filename))
  pdf(filename,
      width=pdf_size_in_inches,
      height=pdf_size_in_inches)
  
  h <- gplots::heatmap.2(mat, 
                         col=colorRamps::matlab.like(200),
                         trace=c('none'),
                         key.title='Relative signal', key.xlab='Log Odds of feature\npresence given PC category', density.info=c('none')
  )
  dev.off()
  # end
  return(filename)
}
plot_sequenced_numbers_over_time <- function(db_connection, date_end = Sys.Date()) {
  # plots a histogram illustrating numbers of samples over time
  #
  # parameters:
  #   database connection, as returned (for example) by get_db_connection

  # returns
  #   ggplot object
  check_table_present(db_connection, 'clinical_metadata','Do you need to call add_cog_metadata()')
  print("plot_sequenced_numbers_over_time: Loading data #1")
  
  # dates may be passed in iso format as strings.
  date_ends <- convert_date(date_end) 
  
  
  cmd <- "select sample_date, adm1, count(*) n_samples
        from clinical_metadata
        where sample_date <= ?
        group by sample_date, adm1;"
 
  over_time <- database_query(
    db_connection = db_connection, 
    cmd = cmd, 
    params = list(date_ends[['iso']]))
  over_time$sample_date <- as.Date(over_time$sample_date)
  p <- ggplot(over_time, aes(x=sample_date,y=n_samples,fill=adm1)) + geom_col()
  p <- p + scale_fill_discrete(name = "Geography")
  p <- p + scale_x_date("Sample collection date")
  p <- p + scale_y_continuous("Number of sequenced isolates")
  p

}

# ---- get clinical and sequence data from COG .csv file ----
read_cog_metadata <- function(cogfile, date_end = Sys.Date()){
  
  date_ends <- convert_date(date_end)
  cog_metadata <- readr::read_csv(cogfile)
  
  # extract the country and sequenceId
  cog_metadata$sample_id <- sapply(strsplit(cog_metadata$sequence_name,"/"), `[`, 2)

  # we store dates as isoformat strings.  This is required by sqlite, and is compatible with searching
  # in other database too.
  cog_metadata$sample_date <- lubridate::format_ISO8601(
    cog_metadata$sample_date, 
    usetz = FALSE, 
    precision = 'ymd')
  cog_metadata <- subset(cog_metadata, sample_date <= date_ends[['iso']])
  
  # sequence features are everything except the fixed fields
  fixed_fields <- c("sequence_name",
                    "country",
                    "adm1",
                    "pillar_2",
                    "sample_date",
                    "lineage_support",
                    "lineages_version",
                    "sample_id",
                    "epi_week",
                    "variants" )
  sequence_features <- setdiff(names(cog_metadata), fixed_fields)
  
  cog_sequence_metadata <- subset(cog_metadata, select=c('sample_id', sequence_features))
  cog_sequence_metadata_long <- reshape2::melt(cog_sequence_metadata, id.vars=c('sample_id'))
  cog_sequence_metadata_long$feature <- paste(cog_sequence_metadata_long$variable,
                                              cog_sequence_metadata_long$value,
                                              sep = ':')
  
  cog_clinical_metadata_fields <- c( "sample_id",
                                     "country", 
                                     "adm1", 
                                      
                                     "sample_date"
                                     )
  cog_clinical_metadata <- subset(cog_metadata, select=cog_clinical_metadata_fields)

  list(
    clinical_metadata = list(
      "data_table"=cog_clinical_metadata, 
      "appropriate_indices"=c("sample_id","sample_date","adm1")
    ),
    
    sequence_metadata = list(
      "data_table"=cog_sequence_metadata_long, 
      "appropriate_indices"=c("sample_id", "variable","feature"))
  )    
}

# ---- database connection & mods ----
convert_date <- function(date_end){
  # sqlite stores dates in isoformat.
  # returns both the isoformat and Date versions of a date, passed in either format
  if (class(date_end) == 'character'){
    date_end_dt <- as.Date(date_end)
    date_end_iso <- date_end
  } else if (class(date_end)== 'Date') {
    date_end_dt <- date_end
    date_end_iso <- lubridate::format_ISO8601(date_end,usetz=FALSE, precision='ymd')
  } else {
    stop(paste("Date end is of the wrong class: allow date and character", date_end, class(date_end)))
  }
  list(dt = date_end_dt, iso = date_end_iso)
}
get_sqlite_connection  <-function(dbfile){
  # gets a database connection object for the sqllite db at dbfile
  # parameters:
  #    dbfile: the database file name
  # returns:
  #    RDBI connection object
  # raises:
  #    stops if file does not exist
  
  # note: to use a different database, you can make a new function
  # which returns an other valid DBI handle.
  
  if (!file.exists(dbfile)){
    stop("File does not exists")
  }
  dbConnect(RSQLite::SQLite(), dbfile)
}
database_query <- function(db_connection, cmd, params= NA) {
  # sends a query to the database, with safe parameter injection
  #
  # parameters:
  #   db_connection: a database connection, as from get_sqlite_connection()
  #   cmd: sql command
  #   params: a list of parameters
  
  # TODO: could cache queries, but 
  # need to encapsulate in an S4 class to make work.
  # can has the combination of db_connection@dbname, cmd, and params to act as a key
  
  if (class(params)=='list') {
    parameterised_query <- dbSendQuery(db_connection, cmd)
    dbBind(parameterised_query, params)
    res <- dbFetch(parameterised_query)
    dbClearResult(parameterised_query)
    return(res)
  } else if (is.na(params)){
    return(dbGetQuery(db_connection, cmd))
  } else {
    stop(paste("Invalid params passed: ", params))
  }
}
ensure_important_indices <- function(db_connection){
  # ensures that various important indexes are in place, particularly those on
  # the transformed_coordinate_categories tables (pc, cat, sample_id)
  # transformed_coordinates (pc, pos)
  important_indices = list(
    "eigenvectors" = c('pc','pos'),
    "transformed_coordinate_categories" = c('pc_cat','sample_id'),
    "sample_id" = c('sample_id')
  )
  
  for (table_name in names(important_indices)) {
    for (index_field in important_indices[[table_name]]){
      
      create_ix_sql = paste0("CREATE INDEX IF NOT EXISTS ", table_name,"_",index_field, " ON ", table_name, "(", index_field, ");")
      print(paste("Ensuring there is an index on table ", table_name, " field ", index_field))
      dbExecute(db_connection, create_ix_sql)

    }
  }
  # compound indexes
  if (db_has_table(db_connection,'transformed_coordinate_categories')){
    print("Building compound indices on pc_cat, sample_id")
    cmds = c("CREATE INDEX IF NOT EXISTS ix1 on transformed_coordinate_categories (pc_cat, sample_id );"
    )
    for (create_ix_sql in cmds){
    dbExecute(db_connection, create_ix_sql)
    }
  }
          
}
list_tables <- function(db_connection){
  # list tables in database to which you are connected
  #
  # parameters:
  #  db_connection: an database connection as generated by get_sqlite_connection
  table_list <- dbListTables(db_connection)
  table_list
}
number_of_pcs <- function(db_connection) {
  # returns the number of pcs in the database
  check_table_present(db_connection, 'explained_variance_ratio', 'Database does not contain PCA results.  ? connected to wrong database')
  cmd = "        select count(*) n_pcs from explained_variance_ratio;"
  res = database_query(db_connection = db_connection, 
                       cmd = cmd,
                       params = NA)

  as.integer(res[1])
}
build_date <- function(db_connection) {
  # returns the date and time the build was made
  check_table_present(db_connection, 'Metadata', 'Database does not contain PCA results.  ? connected to wrong database')
  cmd = "select value as build_time from Metadata where variable = 'build_time';"
  res = database_query(db_connection = db_connection, 
                       cmd = cmd,
                       params = NA)
  lubridate::as_date(res$build_time[1])
}

check_table_present <- function(db_connection, table_name, failure_message) {
    # checks that table_name is a table in the database at db_connection.
    # if it isn't, exits with failure_message
    if (!table_name %in% list_tables(db_connection)){
      msg <- paste("Table not present: ", table_name, "NOTE: ",failure_message)
      stop(msg)
    }
}
db_has_table <- function(db_connection, tbl_name){
  # check whether a table is present
  tbl_name %in% list_tables(db_connection) 
}
add_database_tables <- function(db_connection, list_of_tables, overwrite=TRUE){
  # adds new tables to the database connected on db_connection
  #
  # parameters:
  #  db_connection:  a database connection
  #  list_of_tables: a named list.  The names are to become the table name;
  #                                 the value of each item is a list
  #                                   data_table: the data
  #                                   appropriate_indices: a list of fields to index
  #                 such a list is generated by read_cog_metadata
  #  overwrite:  whether to overwrite the tables if they exist
  for (table_name in names(list_of_tables)) {
    print(paste("Writing ",table_name, " with ", nrow(list_of_tables[[table_name]]$data_table)," rows"))
    # Drop table if it already exists
    if (dbExistsTable(db_connection, table_name)) {
      if (overwrite) {
        dbRemoveTable(db_connection, table_name)
      } else {
        stop("Database table exists, but overwrite is false")
      }
    }
    # Write the data frame to the database
    dbWriteTable(db_connection, name = table_name, value = list_of_tables[[table_name]]$data_table, row.names = FALSE)
    for (index_name in list_of_tables[[table_name]]$appropriate_indices){
      create_ix_sql = paste0("CREATE INDEX ", table_name,"_",index_name, " ON ", table_name, "(", index_name, ");")
      print(paste("Building index", index_name))
      dbExecute(db_connection, create_ix_sql)
    }
  }
}
connect_to_fn4sqlitedb <- function(dbfile){
    # connects to a findNeighbour4 PCA output sqlite database,
    # performs checks to make sure it has the right tables, and
    # makes indexes if they are not there
  
    # parameters:
    #   dbfile: an sqlite database
    # returns:
    #   database connection
    print(paste0("Connecting to: ",dbfile))
    db_connection <- get_sqlite_connection(dbfile)  # open connection

    
    # set various pragmas which maybe improve performance  https://github.com/phiresky/blog/blob/master/posts/2020/sqlite-performance-tuning.md
    dbSendQuery(db_connection, "PRAGMA synchronous = normal;")
    dbSendQuery(db_connection, "PRAGMA temp_store = memory;")
    dbSendQuery(db_connection, "PRAGMA mmap_size = 20000000000;")  # 20G
    
    tbls <- list_tables(db_connection)  # list tables in  db
    print("Connected to database.  Tables are:")
    print(tbls)
    
    essential_tables <- c('transformed_coordinate_categories', 
                          'Metadata',
                          'pos_per_pc',
                          'contributing_pos',
                          'contributing_basepos',
                          'explained_variance_ratio',
                          'eigenvectors',
                          'suspect_quality_seqs',
                          'mix_quality_seqs')
    for (essential_table in essential_tables){
      check_table_present(db_connection, 'transformed_coordinate_categories', 'This should be present in fn4 pca output.  Check you have right database')
    }
    ensure_important_indices(db_connection)

    db_connection
}

add_cog_metadata <- function(db_connection, cogfile, date_end, overwrite=FALSE){
  # if cog meta data is not already in the database, then add it.
  # parameters:
  #   db_connection: the database connection
  #   cogfile: the cog-uk metadata file.
  #   date_end: don't add data after this date: can be either a Date or an isoformat string yyyy-mm-dd
  
  ## load metadata from cogfile if not present
  tbls <- list_tables(db_connection)
  if (!('sequence_metadata' %in% tbls & 'clinical_metadata' %in% tbls) | overwrite) {
    # no data loaded from cog metadata file
    print("Loading cogfile metadata")
    metadata_tables <- read_cog_metadata(cogfile, date_end)
    add_database_tables(db_conn, metadata_tables)
    # compound indexes
    print("Building compound indices on (sample_id, sample_date); (variable, sample_id)")
    cmds = c(
            "CREATE INDEX IF NOT EXISTS ix3 on clinical_metadata ( sample_id, sample_date );",
            "CREATE INDEX IF NOT EXISTS ix2 on sequence_metadata ( variable, sample_id );"
    )
    for (create_ix_sql in cmds){
    dbExecute(db_connection, create_ix_sql)
    }

           
  } else {
    print("sequence metadata is already present.")
  }
}

# ---- linear modelling of counts -----
# ---- linear modelling of counts -----
fit_recent_trends_after_emergence <- function(db_connection,
                                       date_end,
                                       interval_analysed,
                                       analysis_family_id= 'Not_provided',
                                       overwrite=FALSE, 
                                       fitting_failed_dir = '/tmp/',
                                       remove_first_n = 0
                                       )                                                          {
  # fits a Poisson regression model estimating the 
  # count for pc_cat category AFTER it emerges (not over a fixed period, incl.
  # time when it may not exist
  
  # Parameters
  #  db_connnection: the database connection
  #  this_pc:  the principal component to analyse
  #  interval_analysed: the number of days' data to analyse
  #  date_end: the last day of of the analysis period.  For today, set to 
  #     Sys.Date().  Note that date_end must be a Date, not an ISO format 
  #     string representing a date.
  #  analysis_family_id : optional, a reference number for this analysis.
  #                       can be used to link together analysis for different 
  #                       pcs in a database
  #  remove_first_n = 3 don't model the first 3 cases in each pc_cat
  
  
  # startup
  if (!db_has_table(db_connection, "statistical_model_fits") | overwrite){
    
    analysis_uuid <- as.character(uuid::UUIDgenerate())
    
    # recover counts in the relevant time period
    date_ends <- convert_date(date_end)
    date_start <- date_ends[['dt']] - interval_analysed
    date_sequence <- seq(date_start, to = date_ends[['dt']], by =1)
    date_cat <- as.integer(date_sequence-date_ends[['dt']])
    dow_sequence <- lubridate::wday(date_sequence)
    time_series_meta <- data.table(
      sample_date = lubridate::format_ISO8601(date_sequence,usetz=FALSE, precision='ymd'),
      date_cat = factor(date_cat),
      sample_dow = lubridate::wday(date_sequence),
      date_start = as.Date(date_start)
    )
    
    # recover all transformed values
    cmd <- "select * from transformed_coordinate_categories ec;"
    pca_output <- database_query(db_conn, cmd)
    pca_output$tip_label <- paste0('#', pca_output$sample_id)
    pca_output <- merge(pca_output, subset(metadata, select=c('sample_id','expanding_branch', 'sample_date')),
                        by = 'sample_id')
    
    # recover earliest sample per pc_cat
    cmd <- "select ec.pc_cat, min(cm.sample_date) earliest_sample_date
      from transformed_coordinate_categories ec 
      INNER JOIN
      clinical_metadata cm
      ON cm.sample_id = ec.sample_id
      group by ec.pc_cat;"
    
    # compute earliest date per pc
    earliest_per_pc_cat <- database_query(db_conn, cmd)
    one_per_date <- time_series_meta %>% full_join(earliest_per_pc_cat, by=character())
    
    cmd <- "select cm.sample_date, count(*) n
      from transformed_coordinate_categories ec 
      INNER JOIN
      clinical_metadata cm
      ON cm.sample_id = ec.sample_id
      where pc=0
      group by cm.sample_date;"
    perday_cnts <- as.data.table(database_query(db_conn, cmd))
    names(perday_cnts)[2] <- 'n_per_day'
    one_per_day <- one_per_date %>% left_join(perday_cnts, by='sample_date')
    one_per_day$n_per_day[is.na(one_per_day$n_per_day)] <- 0
    
    # now assemble count data per day
    # this data frame contains one entry for each pc_cat/day combination
    cmd <- "select ec.pc_cat, cm.sample_date, cm.sample_id
      from transformed_coordinate_categories ec 
      INNER JOIN
      clinical_metadata cm
      ON cm.sample_id = ec.sample_id
      order by ec.pc_cat, cm.sample_date, cm.sample_id;"
    all_samples <- as.data.table(database_query(db_conn, cmd))
    
    # mark the first n items from each pc_cat
    all_samples$rowid <- 1:nrow(all_samples)
    
    if (remove_first_n > 0) {
      for (i in 1:remove_first_n){
        print("Removing")
        remove_rows <- dcast(all_samples, pc_cat ~ ., value.var='rowid', min)
        names(remove_rows)[2]<- 'rowid'
        all_samples <- all_samples %>% filter(!rowid %in% remove_rows$rowid)
      }
    } else {
      print("All data analysed, including earliest points")
    }
    
    all_count_df <- dcast(all_samples ,pc_cat+sample_date ~ ., value.var='sample_id', length)
    names(all_count_df)[3] <- 'n_initial_removed'
    all_count_df <- merge(one_per_day, all_count_df, by=c('sample_date','pc_cat'), all.x=TRUE)
    all_count_df$n_initial_removed[is.na(all_count_df$n_initial_removed)] <- 0
    all_count_df$pc <- sapply(strsplit(all_count_df$pc_cat, '_', fixed=TRUE), `[`, 1)
    
    all_count_df <- all_count_df %>% mutate(earliest_sample_date = as.Date(earliest_sample_date,
                                                                           start_date = as.Date(start_date)))
    
    recent_pc_cats <- which(all_count_df$earliest_sample_date > all_count_df$date_start)
    all_count_df$date_start[recent_pc_cats] <- all_count_df$earliest_sample_date[recent_pc_cats]
    all_count_df$sample_date_dt <- as.Date(all_count_df$sample_date)
    all_count_df$date_t <- as.integer(all_count_df$sample_date_dt - all_count_df$date_start)
    # we drop any (zero) counts occurring before the first sample for each pc_cat
    all_count_df <- subset(all_count_df, earliest_sample_date <= sample_date)
    
    all_coeffs <- list()
    
    # fitting all the pcs at once is unnecessary and very slow.
    for (this_pc in unique(all_count_df$pc)){
      
      print(paste("Examining ",this_pc))
      count_df <- subset(all_count_df, pc == this_pc & n_per_day >0 ) 
      all_pc_cats <- unique(count_df$pc_cat)
      fit <- tryCatch(
        {
          # Just to highlight: if you want to use more than one 
          # R expression in the "try" part then you'll have to 
          # use curly brackets.
          # 'tryCatch()' will return the last evaluated expression 
          # in case the "try" part was completed successfully
          
          message("Fitting ..")
          
          glm.nb(n_initial_removed ~ date_t*pc_cat + sample_dow,  offset(log(n_per_day)), data=count_df)
          
        },
        error=function(cond) {
          message(paste("Fitting failed"))
          message("Original error message:")
          message(cond)
          
          # write failing data for investigation
          return(NA)
          export_to <- paste0(fitting_failed_dir, paste0(analysis_family_id,'_',this_pc,'.csv'))
          write.table(count_df, outfile, sep='\t', row.names=FALSE)
        },
        finally={
          message("Fit completed")
        }
      )  
      
      if (!is.na(fit)){
        # extract coefficients into a common framework, suitable for storage in a generic data model.
        coeffs<- tibble(data.frame(summary(fit)$coeff)) %>% 
          mutate(param=rownames(summary(fit)$coeff)) %>%
          mutate(Estimate=Estimate,Estimate_CI_low=Estimate-1.96*Std..Error, Estimate_CI_high=Estimate+1.96*Std..Error) %>%
          rename(p_value=Pr...z..) %>%
          rename(Std_Error=Std..Error) %>%
          filter(!grepl("sample_dow",param))
        
        coeffs <- subset(coeffs, select = c('param', 'Estimate', 'Std_Error', 'Estimate_CI_low', 'Estimate_CI_high', 'p_value'))
        
        coeffs$pc_cat <- gsub('pc_cat','',coeffs$param)
        coeffs$pc_cat <- gsub(':date_t','',coeffs$pc_cat)
        model_pc_cats <- unique(coeffs$pc_cat)
        reference_cat <- setdiff(all_pc_cats, model_pc_cats)[1]
        
        coeffs$param_desc <- ifelse(grepl(':date_t|date_t',coeffs$param), "Trend over time relative to ref category", "Rate ratio relative to ref at end of time period")
        coeffs$param_desc[coeffs$pc_cat=='(Intercept)'] <- 'Rate in reference pc_cat at end of time period'
        coeffs$param_desc[coeffs$pc_cat=='date_t'] <- 'Rate of change in reference pc_cat'
        coeffs$pc_cat <- gsub('(Intercept)',reference_cat, coeffs$pc_cat, fixed= TRUE)
        coeffs$pc_cat <- gsub('^date_t$', reference_cat, coeffs$pc_cat)
        coeffs$is_reference <- ifelse(coeffs$pc_cat == reference_cat,TRUE, FALSE)
        
        coeffs$comments <- "Estimates are ln IRRs"
        
        
        # now Var(X+Y) = Var(X) + Var(Y) + 2 Cov(X,Y)
        # but for Var(X*Y) we need to use Delta methods (essentially approximation by Taylor expansion).
        # https://stats.stackexchange.com/questions/62916/confidence-interval-for-the-product-of-two-parameters
        # now we compute rates over time for each component - not relative to the reference category
        # Var(XY) ~= MLE(Y)^2.Var(X) + MLE(X)^2.Var(Y) + 2 MLE(X)MLE(Y)COV(X,Y)
        # Background:  https://migariane.github.io/DeltaMethodEpiTutorial.nb.html
        
        # start with rate ratios
        vc <- vcov(fit)
        
        extract_element<- function(vc, row_name, col_name) {
          r_id <- which(rownames(vc) == row_name)
          c_id <- which(colnames(vc) == col_name)
          vc[r_id, c_id]
        }
        
        MLE_X <- as.numeric(subset(coeffs, param=="date_t")['Estimate'])
        VAR_X <- extract_element(vc, 'date_t','date_t')
        
        res <- list(
          param = paste0("IRR:", reference_cat),
          Estimate = MLE_X,
          Std_Error = as.numeric(subset(coeffs, param=="date_t")['Std_Error']),
          Estimate_CI_low = as.numeric(subset(coeffs, param=="date_t")['Estimate_CI_low']),
          Estimate_CI_high = as.numeric(subset(coeffs, param=="date_t")['Estimate_CI_high']),
          p_value = as.numeric(subset(coeffs, param=="date_t")['p_value']),
          pc_cat = reference_cat,
          param_desc = "Incidence rate ratio",
          is_reference = 0,
          comments = "Incidence ratio ratio computed by Delta method"
        )
        results <- list()
        results[[1]] <- res
        
        i <- 1
        for (this_param in coeffs$param[grepl("date_t:",coeffs$param)]){
          i <- i  +1
          this_pc_cat <- gsub("date_t:pc_cat","",this_param)
          MLE_Y <- as.numeric(subset(coeffs, param==this_param)['Estimate'])
          VAR_Y <- extract_element(vc, this_param,this_param)
          COV_XY <- extract_element(vc, 'date_t',this_param)
          VAR_XY <- (MLE_Y^2)*VAR_X + (MLE_X^2)*VAR_Y + 2*MLE_X*MLE_Y*COV_XY
          SE_XY <- sqrt(VAR_XY)
          MLE_XY <- MLE_X * MLE_Y
          Estimate_CI_low <- MLE_XY  - 1.96*SE_XY
          Estimate_CI_high <- MLE_XY + 1.96*SE_XY
          Z <- MLE_XY / SE_XY
          p_value <- pnorm(as.numeric(Z))
          
          res <- list(
            param = paste0("IRR:", this_pc_cat),
            Estimate = MLE_XY,
            Std_Error = SE_XY,
            Estimate_CI_low = Estimate_CI_low,
            Estimate_CI_high = Estimate_CI_high,
            p_value = p_value,
            pc_cat = this_pc_cat,
            param_desc = "Incidence rate ratio",
            is_reference = 0,
            comments = "Incidence ratio ratio computed by Delta method"
          )
          results[[i]] <- res
          
        }
        
        
        ## next we consider the rate estimates themselves
        MLE_X <- as.numeric(subset(coeffs, param=="(Intercept)")['Estimate'])
        VAR_X <- extract_element(vc, '(Intercept)','(Intercept)')
        
        res <- list(
          param = paste0("countPerDay:", reference_cat),
          Estimate = MLE_X,
          Std_Error = as.numeric(subset(coeffs, param=="(Intercept)")['Std_Error']),
          Estimate_CI_low = as.numeric(subset(coeffs, param=="(Intercept)")['Estimate_CI_low']),
          Estimate_CI_high = as.numeric(subset(coeffs, param=="(Intercept)")['Estimate_CI_high']),
          p_value = p_value,
          pc_cat = reference_cat,
          param_desc = "Estimated Counts per day",
          is_reference = 0,
          comments = "Estimated counts per day computed by Delta method"
        )
        i <- i + 1
        results[[i]] <- res
        
        
        for (this_param in coeffs$param[grepl("^pc_cat",coeffs$param)]){
          i <- i  +1
          this_pc_cat <- gsub("pc_cat","",this_param)
          MLE_Y <- as.numeric(subset(coeffs, param==this_param)['Estimate'])
          VAR_Y <- extract_element(vc, this_param,this_param)
          COV_XY <- extract_element(vc, '(Intercept)',this_param)
          VAR_XY <- (MLE_Y^2)*VAR_X + (MLE_X^2)*VAR_Y + 2*MLE_X*MLE_Y*COV_XY
          SE_XY <- sqrt(VAR_XY)
          MLE_XY <- MLE_X * MLE_Y
          Estimate_CI_low <- MLE_XY  - 1.96*SE_XY
          Estimate_CI_high <- MLE_XY + 1.96*SE_XY
          Z <- MLE_XY / SE_XY
          p_value <- pnorm(as.numeric(Z))
          
          res <- list(
            param = paste0("countPerDay:", this_pc_cat),
            Estimate = MLE_XY,
            Std_Error = SE_XY,
            Estimate_CI_low = Estimate_CI_low,
            Estimate_CI_high = Estimate_CI_high,
            p_value = p_value,
            pc_cat = this_pc_cat,
            param_desc = "Estimated Counts per day",
            is_reference = 0,
            comments = "Estimated counts per day computed by Delta method"
          )
          results[[i]] <- res
          
        }
        results_df <- data.frame(do.call(rbind.data.frame, results))
        coeffs <- data.frame(rbind(coeffs,results_df
        ))
        
        coeffs$analysis_id <- analysis_uuid
        coeffs$has_CI <- 1
        coeffs$Estimate2NaturalSpace <- "exp"
        all_coeffs[[as.character(this_pc)]] <- coeffs
      }
    }
    print("Fits done")
    
    all_coeffs_df <- data.frame(do.call(rbind.data.frame, all_coeffs))
    rownames(all_coeffs_df) <- 1:nrow(all_coeffs_df)
    
    metadata <- data.frame(list(
      'analysis_family_id' = analysis_family_id,
      'analysis_id' = analysis_uuid,
      'analysis' = 'fit_recent_trend_in_counts',
      'analysis_type' = 'Negative binomial regression',
      'readable_info' = 'Estimates rate of change of isolations of each pc_cat.  Rates relative to most common cat, and overall, are provided.  Assmues linear trend.  Controls for numbers of samples analysed each day and so reflects rates per sampled population',
      'other_parameters'='{}',
      'pc' = this_pc,
      'date_end' = date_ends[['iso']],
      'interval_analysed_days' = interval_analysed
    ))
    print("Packaged result with metadata")
    retVal <- list('fitted' = 1, 'reason' = 'Success','meta'= metadata, 'model_fit'= all_coeffs_df)
  
    
    add_database_tables(
      db_connection,
      list(
        statistical_model_metadata=list(
          data_table=metadata,
          appropriate_indices=list('analysis_id')
        ),
        statistical_model_fits=list(
          data_table=all_coeffs_df,
          appropriate_indices=list('analysis_id','pc_cat','param','p_value','Estimate')
        )
      ),
      overwrite=TRUE)
  } else {
    print("Using stored model")
  }
  
}

fit_recent_trend_in_counts <- function(db_connection,
                                                  date_end,
                                                  interval_analysed,
                                                  analysis_family_id= 'Not_provided',
                                                  overwrite=FALSE, 
                                                  fitting_failed_dir = '/tmp/')                                                          {
  # fits a Poisson regression model estimating the 
  # count for pc_cat category 
  
  # Parameters
  #  db_connnection: the database connection
  #  this_pc:  the principal component to analyse
  #  interval_analysed: the number of days' data to analyse
  #  date_end: the last day of of the analysis period.  For today, set to 
  #     Sys.Date().  Note that date_end must be a Date, not an ISO format 
  #     string representing a date.
  #  analysis_family_id : optional, a reference number for this analysis.
  #                       can be used to link together analysis for different 
  #                       pcs in a database
  

  # startup
  if (!db_has_table(db_connection, "statistical_model_fits")){

    # generate guid for this analysis
    analysis_uuid <- as.character(uuid::UUIDgenerate())

    date_ends <- convert_date(date_end)
    count_result <- count_per_all_pc_cats_in_time_interval(
        db_connection = db_connection,
        date_end = date_ends[['dt']], 
        interval_analysed= interval_analysed)
    
    if (!count_result$success==TRUE){
        return(count_result)  # failed, maybe no data
    } else {
        all_count_df <- count_result$count_df
    }


    all_count_df$pc <- sapply(strsplit(all_count_df$pc_cat, '_', fixed=TRUE), `[`, 1)

    all_coeffs <- list()
    # fitting all the pcs at once is unnecessary and very slow.
    for (this_pc in unique(all_count_df$pc)){
        print(paste("Examining ",this_pc))
        count_df <- subset(all_count_df, pc == this_pc)
        
        # compute counts per pc_cat.  We don't fit if < 3 counts in the interval studied, as we can't fit a regression line
        # meaningfully to less than 3 points
        
        total_cnts <- dcast(count_df, pc_cat ~ ., fun.aggregate=sum,  value.var='n')
        names(total_cnts)[2] <- 'total_n'
        total_cnts <- subset(total_cnts, total_n>=3)
        count_df <- merge(total_cnts, count_df, by='pc_cat')
        
        perday_cnts <- dcast(count_df, sample_date ~ ., fun.aggregate=sum,  value.var='n')
        names(perday_cnts)[2] <- 'per_day_n'
        perday_cnts <- subset(perday_cnts, per_day_n > 0) # no information if no observations
        count_df <- merge(perday_cnts, count_df, by='sample_date')
        all_pc_cats <- unique(count_df$pc_cat)
        # model secular time as a spline (there will be less data at the end).  
        # use negative binomial models, as poisson assumptions will be violated due to outbreaks (non-independence)
        # include day of the week

        # model pc_cat as a factor, and regard the most frequenct pc_cat as the reference
        #total_counts <- count_df %>% group_by(pc_cat) %>% summarise(total_n=sum(n)) %>% arrange(total_n)
        #count_df$pc_cat <- factor(count_df$pc_cat,levels=rev(total_counts$pc_cat))

        fit <- tryCatch(
                {
                    # Just to highlight: if you want to use more than one 
                    # R expression in the "try" part then you'll have to 
                    # use curly brackets.
                    # 'tryCatch()' will return the last evaluated expression 
                    # in case the "try" part was completed successfully

                    message("Fitting ..")

                    glm.nb(n ~ date_t*pc_cat + sample_dow,  offset(log(per_day_n)), data=count_df)

                },
                error=function(cond) {
                    message(paste("Fitting failed"))
                    message("Original error message:")
                    message(cond)
                    
                    # write failing data for investigation
                    return(NA)
                    export_to <- paste0(fitting_failed_dir, paste0(analysis_family_id,'_',this_pc,'.csv'))
                    write.table(count_df, outfile, sep='\t', row.names=FALSE)
                },
                finally={
                    message("Fit completed")
                }
            )              
        if (!is.na(fit)){
          # extract coefficients into a common framework, suitable for storage in a generic data model.
          coeffs<- tibble(data.frame(summary(fit)$coeff)) %>% 
            mutate(param=rownames(summary(fit)$coeff)) %>%
            mutate(Estimate=Estimate,Estimate_CI_low=Estimate-1.96*Std..Error, Estimate_CI_high=Estimate+1.96*Std..Error) %>%
            rename(p_value=Pr...z..) %>%
            rename(Std_Error=Std..Error) %>%
            filter(!grepl("sample_dow",param))

          coeffs <- subset(coeffs, select = c('param', 'Estimate', 'Std_Error', 'Estimate_CI_low', 'Estimate_CI_high', 'p_value'))
          
          coeffs$pc_cat <- gsub('pc_cat','',coeffs$param)
          coeffs$pc_cat <- gsub(':date_t','',coeffs$pc_cat)
          model_pc_cats <- unique(coeffs$pc_cat)
          reference_cat <- setdiff(all_pc_cats, model_pc_cats)[1]
          
          coeffs$param_desc <- ifelse(grepl(':date_t|date_t',coeffs$param), "Trend over time relative to ref category", "Rate ratio relative to ref at end of time period")
          coeffs$param_desc[coeffs$pc_cat=='(Intercept)'] <- 'Rate in reference pc_cat at end of time period'
          coeffs$param_desc[coeffs$pc_cat=='date_t'] <- 'Rate of change in reference pc_cat'
          coeffs$pc_cat <- gsub('(Intercept)',reference_cat, coeffs$pc_cat, fixed= TRUE)
          coeffs$pc_cat <- gsub('^date_t$', reference_cat, coeffs$pc_cat)
          coeffs$is_reference <- ifelse(coeffs$pc_cat == reference_cat,TRUE, FALSE)
          
          coeffs$comments <- "Estimates are ln IRRs"
    
          
          # now Var(X+Y) = Var(X) + Var(Y) + 2 Cov(X,Y)
          # but for Var(X*Y) we need to use Delta methods (essentially approximation by Taylor expansion).
          # https://stats.stackexchange.com/questions/62916/confidence-interval-for-the-product-of-two-parameters
          # now we compute rates over time for each component - not relative to the reference category
          # Var(XY) ~= MLE(Y)^2.Var(X) + MLE(X)^2.Var(Y) + 2 MLE(X)MLE(Y)COV(X,Y)
          # Background:  https://migariane.github.io/DeltaMethodEpiTutorial.nb.html
          
          # start with rate ratios
          vc <- vcov(fit)
          
          extract_element<- function(vc, row_name, col_name) {
              r_id <- which(rownames(vc) == row_name)
              c_id <- which(colnames(vc) == col_name)
              vc[r_id, c_id]
          }
          
          MLE_X <- as.numeric(subset(coeffs, param=="date_t")['Estimate'])
          VAR_X <- extract_element(vc, 'date_t','date_t')
          
          res <- list(
              param = paste0("IRR:", reference_cat),
              Estimate = MLE_X,
              Std_Error = as.numeric(subset(coeffs, param=="date_t")['Std_Error']),
              Estimate_CI_low = as.numeric(subset(coeffs, param=="date_t")['Estimate_CI_low']),
              Estimate_CI_high = as.numeric(subset(coeffs, param=="date_t")['Estimate_CI_high']),
              p_value = as.numeric(subset(coeffs, param=="date_t")['p_value']),
              pc_cat = reference_cat,
              param_desc = "Incidence rate ratio",
              is_reference = 0,
              comments = "Incidence ratio ratio computed by Delta method"
          )
          results <- list()
          results[[1]] <- res
          
          i <- 1
          for (this_param in coeffs$param[grepl("date_t:",coeffs$param)]){
              i <- i  +1
              this_pc_cat <- gsub("date_t:pc_cat","",this_param)
              MLE_Y <- as.numeric(subset(coeffs, param==this_param)['Estimate'])
              VAR_Y <- extract_element(vc, this_param,this_param)
              COV_XY <- extract_element(vc, 'date_t',this_param)
              VAR_XY <- (MLE_Y^2)*VAR_X + (MLE_X^2)*VAR_Y + 2*MLE_X*MLE_Y*COV_XY
              SE_XY <- sqrt(VAR_XY)
              MLE_XY <- MLE_X * MLE_Y
              Estimate_CI_low <- MLE_XY  - 1.96*SE_XY
              Estimate_CI_high <- MLE_XY + 1.96*SE_XY
              Z <- MLE_XY / SE_XY
              p_value <- pnorm(as.numeric(Z))
              
              res <- list(
              param = paste0("IRR:", this_pc_cat),
              Estimate = MLE_XY,
              Std_Error = SE_XY,
              Estimate_CI_low = Estimate_CI_low,
              Estimate_CI_high = Estimate_CI_high,
              p_value = p_value,
              pc_cat = this_pc_cat,
              param_desc = "Incidence rate ratio",
              is_reference = 0,
              comments = "Incidence ratio ratio computed by Delta method"
              )
              results[[i]] <- res
          
          }
          
          
          ## next we consider the rate estimates themselves
          MLE_X <- as.numeric(subset(coeffs, param=="(Intercept)")['Estimate'])
          VAR_X <- extract_element(vc, '(Intercept)','(Intercept)')
          
          res <- list(
              param = paste0("countPerDay:", reference_cat),
              Estimate = MLE_X,
              Std_Error = as.numeric(subset(coeffs, param=="(Intercept)")['Std_Error']),
              Estimate_CI_low = as.numeric(subset(coeffs, param=="(Intercept)")['Estimate_CI_low']),
              Estimate_CI_high = as.numeric(subset(coeffs, param=="(Intercept)")['Estimate_CI_high']),
              p_value = p_value,
              pc_cat = reference_cat,
              param_desc = "Estimated Counts per day",
              is_reference = 0,
              comments = "Estimated counts per day computed by Delta method"
          )
          i <- i + 1
          results[[i]] <- res
          
          
          for (this_param in coeffs$param[grepl("^pc_cat",coeffs$param)]){
              i <- i  +1
              this_pc_cat <- gsub("pc_cat","",this_param)
              MLE_Y <- as.numeric(subset(coeffs, param==this_param)['Estimate'])
              VAR_Y <- extract_element(vc, this_param,this_param)
              COV_XY <- extract_element(vc, '(Intercept)',this_param)
              VAR_XY <- (MLE_Y^2)*VAR_X + (MLE_X^2)*VAR_Y + 2*MLE_X*MLE_Y*COV_XY
              SE_XY <- sqrt(VAR_XY)
              MLE_XY <- MLE_X * MLE_Y
              Estimate_CI_low <- MLE_XY  - 1.96*SE_XY
              Estimate_CI_high <- MLE_XY + 1.96*SE_XY
              Z <- MLE_XY / SE_XY
              p_value <- pnorm(as.numeric(Z))
              
              res <- list(
              param = paste0("countPerDay:", this_pc_cat),
              Estimate = MLE_XY,
              Std_Error = SE_XY,
              Estimate_CI_low = Estimate_CI_low,
              Estimate_CI_high = Estimate_CI_high,
              p_value = p_value,
              pc_cat = this_pc_cat,
              param_desc = "Estimated Counts per day",
              is_reference = 0,
              comments = "Estimated counts per day computed by Delta method"
              )
              results[[i]] <- res
              
          }
          results_df <- data.frame(do.call(rbind.data.frame, results))
          coeffs <- data.frame(rbind(coeffs,results_df
                  ))

          coeffs$analysis_id <- analysis_uuid
          coeffs$has_CI <- 1
          coeffs$Estimate2NaturalSpace <- "exp"
          all_coeffs[[this_pc]] <- coeffs
          
        }
    }

    all_coeffs_df <- data.frame(do.call(rbind.data.frame, all_coeffs))
    rownames(all_coeffs_df) <- 1:nrow(all_coeffs_df)

    metadata <- data.frame(list(
        'analysis_family_id' = analysis_family_id,
        'analysis_id' = analysis_uuid,
        'analysis' = 'fit_recent_trend_in_counts',
        'analysis_type' = 'Negative binomial regression',
        'readable_info' = 'Estimates rate of change of isolations of each pc_cat.  Rates relative to most common cat, and overall, are provided.  Assmues linear trend.  Controls for numbers of samples analysed each day and so reflects rates per sampled population',
        'other_parameters'='{}',
        'pc' = this_pc,
        'date_end' = date_ends[['iso']],
        'interval_analysed_days' = interval_analysed
    ))
    retVal <- list('fitted' = 1, 'reason' = 'Success','meta'= metadata, 'model_fit'= all_coeffs_df)
    #saveRDS(retVal, file = "megbin_model.rds")
    #print("Wrote data to negbin_model.rds")


    add_database_tables(
    db_connection,
        list(
            statistical_model_metadata=list(
                                        data_table=metadata,
                                        appropriate_indices=list('analysis_id')
                                    ),
            statistical_model_fits=list(
                                        data_table=all_coeffs_df,
                                        appropriate_indices=list('analysis_id','pc_cat','param','p_value','Estimate')
                                    )
        ),
        overwrite=TRUE)
  } else {
      print("Using stored model")
  }

}


fit_recent_trend_wrt_most_common_pcat <- function(db_connection,
                                                  this_pc,
                                                  date_end,
                                                  interval_analysed,
                                                  analysis_family_id= 'Not provided')                                                          {
  # fits a Poisson regression model estimating the 
  # proportion of samples of each category of this_pc, 
  # relative to the most common category;
  # as well as the incidence rate ratio (relative growth)
  # over the past interval_analysed days to date_end
  
  # Parameters
  #  db_connnection: the database connection
  #  this_pc:  the principal component to analyse
  #  interval_analysed: the number of days' data to analyse
  #  date_end: the last day of of the analysis period.  For today, set to 
  #     Sys.Date().  Note that date_end must be a Date, not an ISO format 
  #     string representing a date.
  #  analysis_family_id : optional, a reference number for this analysis.
  #                       can be used to link together analysis for different 
  #                       pcs in a database
  

  ## DEPRECATED
  stop("Deprecated fucntion")
  # startup

  # generate guid for this analysis
  analysis_uuid <- uuid::UUIDgenerate()

  date_ends <- convert_date(date_end)
  count_result <- count_per_pc_in_time_interval(
    db_connection = db_connection,
    this_pc = this_pc,
    date_end = date_ends[['dt']], 
    interval_analysed= interval_analysed)
  
  if (!count_result$fitted==TRUE){
    return(count_result)  # failed, maybe no data
  } else {
    count_df <- count_result$count_df
  }
  total_counts <- count_df %>% group_by(pc_cat) %>% summarise(total_n=sum(n)) %>% arrange(total_n)
  total_counts$fit_cat <- paste0('cat', total_counts$pc_cat)

  reference_cat <- as.character(levels(count_df$cat)[1])

  # fit, contrasting to the most common category
  # sample_dow reflects day of week
  fit <- glm(n ~ cat * date_t + sample_dow, family="poisson", data=count_df)
  
  # extract coefficients into a common framework, suitable for storage in a generic data model.
  coeffs<- tibble(data.frame(summary(fit)$coeff)) %>% 
    mutate(param=rownames(summary(fit)$coeff)) %>%
    filter(!grepl("sample_dow",param)) %>%
    mutate(Estimate=Estimate,Estimate_CI_low=Estimate-1.96*Std..Error, Estimate_CI_high=Estimate+1.96*Std..Error) %>%
    rename(p_value=Pr...z..) %>%
    rename(Std_Error=Std..Error) %>%
    select(param, Estimate, Std_Error, Estimate_CI_low, Estimate_CI_high,p_value)

  coeffs$comments <- "Estimates are ln IRRs/rates"
  pc_classif<-sub("[^0-9]+([0-9]+)[^0-9]*","\\1",coeffs$param)
  baseline_locs <- grep("date_t|\\(Intercept\\)", pc_classif)
  is_baseline <- rep(FALSE, length(pc_classif))
  is_baseline[baseline_locs] <- TRUE
  pc_classif[baseline_locs] <- levels(count_df$cat)[1]
  #coeffs <- coeffs %>%
  #  mutate(baseline=is_baseline)

  coeffs$pc_cat <- gsub('cat','',coeffs$param)
  coeffs$pc_cat <- gsub(':date_t','',coeffs$pc_cat)
  coeffs$param <- ifelse(grepl(':date_t|date_t',coeffs$param), "Trend over time relative to ref category", "Rate ratio relative to ref at baseline")
  coeffs$param[coeffs$pc_cat=='(Intercept)'] <- 'Rate in most common pc_cat'
  coeffs$param[coeffs$pc_cat=='date_t'] <- 'Rate of change in most common pc_cat'
  coeffs$pc_cat <- gsub('(Intercept)',reference_cat,coeffs$pc_cat, fixed= TRUE)
  coeffs$pc_cat <- gsub('date_t',reference_cat,coeffs$pc_cat, fixed=TRUE)
 
  coeffs$is_reference <- ifelse(coeffs$pc_cat == reference_cat, TRUE, FALSE)
  coeffs <- merge(coeffs, subset(total_counts, select = c('pc_cat','total_n')), by = 'pc_cat')
  coeffs$analysis_id <- analysis_uuid
  coeffs$has_CI <- TRUE 
  coeffs$Estimate2NaturalSpace <- "exp"
  
  metadata <- data.frame(list(
    'analysis_family_id' = analysis_family_id,
    'analysis_id' = analysis_uuid,
    'analysis' = 'fit_recent_trend_wrt_most_common_pcat',
    'analysis_type' = 'Poisson regression',
    'readable_info' = 'Estimates relative preponderance of each PC category,\nrelative to the most common category.  Additionally, estimates linear trend in these',
    'other_parameters'='{}',
    'pc' = this_pc,
    'date_end' = date_ends[['iso']],
    'interval_analysed_days' = interval_analysed
  ))
  list('fitted' = TRUE, 'reason' = 'Success','meta'= metadata, 'model_fit'= coeffs)
}
# ------ association between PCs and known features  ---------------
compute_or <- function(rowid, a,b,c,d){
  # constructs a 2x2 contingency table
  #
  #  Feature    Present       Absent
  #  PC-CAT Y      a            c          a+c
  #         N      b            d          b+d
  #               a+b           c+d       a+b+c+d
  #
  # reports OR and G test as a list, with rowid in the list.
  #
  #
  mat <- matrix(c(a,b,c,d),nrow=2, ncol=2)
  
  gt <- DescTools::GTest(mat)
  res <- list(G=gt$statistic,
              df=gt$parameter,
              p_value = gt$p.value)

  # given pc_cat present, neither any row nor any col can be all zeros
  if (b>0){
    odds_feature_present <- a/b
  } else {
    odds_feature_present <- Inf
  }
  # given pc_cat absent
  if (d>0){
  odds_feature_absent <- c/d
  } else {
    odds_feature_absent <- Inf
  }
  
  # we know not all the cells are zero
  if (is.infinite(odds_feature_present)) {
    odds_ratio <- Inf
  } else if (is.infinite(odds_feature_absent)) {
    odds_ratio <- 0
  } else {
    odds_ratio <- odds_feature_present/ odds_feature_absent
  }
  #res['odds_feature_present'] <- odds_feature_present
  #res['odds_feature_absent']  <- odds_feature_absent
  #res['odds_ratio'] <- odds_ratio
  
  # generate a truncated OR for the purposes of plotting
  # truncate at 10-6 and 10-6
  res['truncated_odds_ratio'] <- odds_ratio
  if (is.infinite(odds_ratio)){
    res['truncated_odds_ratio'] <- 1e6
  } else if (odds_ratio > 1e6){
    res['truncated_odds_ratio'] <- 1e6
  } else if (odds_ratio < 1e-6){
    res['truncated_odds_ratio'] <- 1e-6
  }
  res['rowid'] <- rowid
  res
}
make_contingency_tables <- function(db_connection,
                                    date_end=Sys.Date(),
                                    overwrite = FALSE,
                                    only_pc_cat_arising_after = NA){
  # aims to construct a series of 2x2 contingency tables
  #           |  LINEAGE/OTHER FEATURE
  # ======================================================
  #           |   Present       Absent
  #  PC-CAT Y |     a            b          a+b
  #         N |     c            d          c+d
  #               a+c           b+d       a+b+c+d
  
  
  # Parameters:
  #  db_connection: the database connection
  #  date_end: do not use samples with collection dates later 
  #            than this.  Defaults to today.  ** CURRENTLY IGNORED **
  #  overwrite: whether to re-write data if present.
  
  # We use a series of SQL queries to measure
  # a+b+c+d : this is the number of samples from which the PCA was build, which
  #           can be obtained from the sample_id table.
  # a+b     : this can be obtained for all pc-categories
  # a+c     : this can be obtained for all features
  # a       : this can also be readily obtained, but the SQL is slower, and
  #           it is best run pc by pc.
 
  if (overwrite==FALSE & db_has_table(db_connection, 'feature_associations')) {
    # already computed
    print("Using stored associations.")
    return(0)
  }
  
  date_ends <- convert_date(date_end) 
  
  a_b_c_d_cmd = "select count(distinct cm.sample_id) n from 
  clinical_metadata cm
  inner join sequence_metadata sm
  on cm.sample_id = sm.sample_id
  where cm.sample_date <= ?;"
  
  # count by feature - restrict to lineages
  a_c_cmd = "select sm.feature, count(distinct cm.sample_id) n from 
  sequence_metadata sm
  INNER JOIN clinical_metadata cm 
  on cm.sample_id = sm.sample_id
  where cm.sample_date <= ?
  and  sm.variable = 'lineage'
  group by sm.feature;"
  
  a_b_cmd = "select pc_cat, count(distinct cm.sample_id) n from 
  transformed_coordinate_categories tcc
  INNER JOIN clinical_metadata cm 
  on cm.sample_id = tcc.sample_id
  where cm.sample_date <= ?
  group by pc_cat;"
  
  a_cmd= "select pc_cat,sm.feature, count(distinct cm.sample_id) n from 
  transformed_coordinate_categories tcc
  INNER JOIN 
  sequence_metadata sm
  on tcc.sample_id = sm.sample_id
  INNER JOIN clinical_metadata cm 
  on cm.sample_id = sm.sample_id
  where cm.sample_date <= ?
  and pc = ? 
  and  sm.variable = 'lineage'
  group by pc,cat, sm.feature;"
  
  a_pc_cat_cmd= "select pc_cat,sm.feature, count(distinct cm.sample_id) n from 
  transformed_coordinate_categories tcc
  INNER JOIN 
  sequence_metadata sm
  on tcc.sample_id = sm.sample_id
  INNER JOIN clinical_metadata cm 
  on cm.sample_id = sm.sample_id
  where cm.sample_date <= ?
  and tcc.pc_cat = ? 
  and  sm.variable = 'lineage'
  group by pc,cat, sm.feature;"

  # check the relevant data exists
  check_table_present(db_connection, 'transformed_coordinate_categories','Table should be present, produced by fn4pca.py')
  check_table_present(db_connection, 'sequence_metadata','Call add_cog_metadata() or equivalent first')
  
  print("Recovering counts")
  
  a_b_c_d <- database_query(
    db_connection = db_connection, 
    cmd = a_b_c_d_cmd, 
    params = list(date_ends[['iso']]))
  print(paste("number of sequences from which the model was built is ",a_b_c_d))

  print("Recovering marginal totals")
  a_c <- database_query(db_connection = db_connection, 
                        cmd = a_c_cmd, 
                        params = list(date_ends[['iso']]))
  names(a_c)[2] <- 'a_c'

  print(paste("number of features (restricted to lineages) is ",nrow(a_c)))
  a_b <- database_query(db_connection = db_connection, 
                        cmd = a_b_cmd, 
                        params = list(date_ends[['iso']]))

  names(a_b)[2] <- 'a_b'

  print(paste("number of pc_cats is ",nrow(a_b)))
  n_pc <- number_of_pcs(db_connection = db_connection)
  coexistence = list()
 
  paste(paste0("Considering whether to select pc_cats arising after ",only_pc_cat_arising_after))
  
  if (is.na(only_pc_cat_arising_after)){
    print(paste("Analysing all pc_cats.  number of pcs is ",n_pc))
    
    print("Recovering co-existence of pc_cats and features, pc by pc")
    for (this_pc in 1:n_pc){    
      print(this_pc)
      coexistence[[this_pc]] <- database_query(
        db_connection = db_connection,
        cmd = a_cmd, 
        params = list(date_ends[['iso']], this_pc)) 
    }
  } else {
    only_pc_cat_arising_after_dt <- convert_date(only_pc_cat_arising_after)

    print(paste0("Focusing on pc_cats observed after ",only_pc_cat_arising_after))
    print("Selecting pc_cats to study.")
    cmd = "select tcc.pc_cat pc_cat, 
  min(sample_date) earliest_date, 
  max(sample_date) latest_date,
  count(*) n_samples
                          from 
                          transformed_coordinate_categories tcc 
                          LEFT JOIN
                          clinical_metadata cm
                          on cm.sample_id = tcc.sample_id
                          group by tcc.pc_cat
                          having min(sample_date)>?;"
    recently_arisen <- database_query(
        db_connection = db_connection,
        cmd = cmd, 
        params = list(only_pc_cat_arising_after_dt[['iso']])
    )

    print(paste0("Selected ",nrow(recently_arisen)," pc_cats for analysis"))
    for (this_pc_cat in recently_arisen$pc_cat){
      print(this_pc_cat)
      coexistence[[this_pc_cat]] <- database_query(
        db_connection = db_connection,
        cmd = a_pc_cat_cmd, 
        params = list(date_ends[['iso']], this_pc_cat)) 
    }
  }

  coexistence_df <- do.call(rbind, coexistence)
 
  names(coexistence_df)[3] <- 'a'
 
  
  coexistence_df <- merge(coexistence_df, a_c, by = 'feature')

  coexistence_df <- merge(coexistence_df, subset(a_b, select=c('a_b','pc_cat')), by = 'pc_cat')

  # if there are no lineages to associated with, terminate
  if (nrow(coexistence_df)>0){
    coexistence_df$a_b_c_d <- as.integer(a_b_c_d)
  
    # now solve for b,c,d.
    coexistence_df$c <- coexistence_df$a_c - coexistence_df$a 
    coexistence_df$b <- coexistence_df$a_b - coexistence_df$a 
    coexistence_df$d <- coexistence_df$a_b_c_d - (coexistence_df$a +
      coexistence_df$b + 
      coexistence_df$c)
    
    # sanity check
    stopifnot(coexistence_df$a +
            coexistence_df$b + 
            coexistence_df$c + coexistence_df$d==coexistence_df$a_b_c_d)
  
    #print(subset(coexistence_df, b<0 | c<0 |a<0 | d<0))
    # compute OR for each row.
    print("Computing ORs")
    coexistence_df$rowid <-1:nrow(coexistence_df)
    ors <- list()
    for (rowid in 1:nrow(coexistence_df)){
      if (rowid %% 100000 ==0){
        print(rowid)
      }
      ors[[rowid]] <- compute_or(rowid,
                                 coexistence_df$a[rowid],
                                 coexistence_df$b[rowid],
                                 coexistence_df$c[rowid],
                                 coexistence_df$d[rowid])
    }
    ors_df <- plyr::ldply(ors, data.frame)
    coexistence_df <- merge(coexistence_df, ors_df, by='rowid')
  } else {
    warning("No 2x2 contingency tables could be constructed. Is there > 1 lineage defined?")
    return(NA)
  } 
  # write to database
  add_database_tables(
    db_connection,
    list(feature_associations=list(data_table=coexistence_df,
                                   appropriate_indices=c(
                                                         'pc_cat',
                                                         'feature',
                                                         'p_value')
                               )
    ),
    overwrite=TRUE)
    
  coexistence_df
}
