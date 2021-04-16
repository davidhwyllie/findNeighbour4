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
library(RSQLite)
library(reshape2)
library(dplyr)
library(uuid)
library(DBI)
library(dplyr)
library(ggplot2)
library(cowplot)
library(viridis)
library(data.table)
library(readr)
library(lubridate)
library(ggrepel)
library(tidyr)
library(DescTools)
library(gplots)
library(colorRamps)

# functions

# ---- helper functions for visualisation ----
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
to_tiff <- function(to_plot, outputdir, db_name, plot_name) {
  # outputs a ggplot p to filename=filename
  # parameters
  #  to_plot : ggplot object
  #  filename: file to which it is to be written
  
  filename <- make_output_filename(outputdir=outputdir,
                                   db_name=db_name,
                                   plot_name=plot_name,
                                   file_type='tiff')
  
  print(paste("Writing file to ",filename))
  tiff(filename,width=1200,height=1000, compression='lzw')
  print(to_plot)
  dev.off()

}
# ---- queries used multiply by visualisation and analysis ----
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
  
  print(paste("load_trend_summary; including clinical data up to ",date_ends[['iso']]))
  check_table_present(db_connection, 'transformed_coordinate_categories','Table should be present, produced by fn4pca.py')
  check_table_present(db_connection, 'trend_model_fits','Call fit_recent_trend_wrt_most_common_pcat_all_pcs() first')

  cmd = "
            select pc, cat, count(*) n_samples, 
            min(transformed_coordinate) min_trans_coord, 
            max(transformed_coordinate) max_trans_coord, 
            avg(transformed_coordinate) mean_trans_coord
            from transformed_coordinate_categories
            group by pc, cat;"

  summary_trends <- database_query(
    db_connection = db_connection, 
    cmd = cmd)

  if (!db_has_table(db_connection, 'clinical_metadata')){
    summary_trends$earliest_sample_date <- NA
    summary_trends$latest_sample_date <- NA
    summary_trends$n_available_sample_dates <- 0
  } else {
    
    cmd <- "select tcc.pc, tcc.cat, 
                      min(sample_date) earliest_date, 
                      max(sample_date) latest_date, 
                      count(sample_date) n_available_sample_dates 
                      from 
                      transformed_coordinate_categories tcc 
                      LEFT JOIN
                      clinical_metadata cm
                      on cm.sample_id = tcc.sample_id
                      where cm.sample_date <= ?
                      group by pc, cat;"
    extra_dates <- database_query(
      db_connection = db_connection,
      cmd = cmd, 
      params = list(date_ends[['iso']]))
    
    summary_trends <- merge(summary_trends, extra_dates,  all.x=TRUE, by=c('pc','cat'))
  }
  
  summary_trends$pc_cat <- paste0(summary_trends$pc,'_',summary_trends$cat)
  cmd = "SELECT pc_cat, tr_p_value, tr_IRR FROM trend_model_fits;"
  trend <- database_query(
    db_connection = db_connection, 
    cmd = cmd)
  ntest <- nrow(trend)
  
  summary_trends <- merge(summary_trends, trend, by='pc_cat')
  
  # convert the date fields to dates
  # this is required by SQLite, which doesn't have a native date format and requires that
  # dates are stored as strings.  We do this in ymd iso format, and and cast as date using as.Date
  summary_trends$earliest_date <- as.Date(summary_trends$earliest_date)
  summary_trends$latest_date <- as.Date(summary_trends$latest_date)

  # apply Bonferroni correction to trends
  adj_p <- 0.01 / ntest
  summary_trends$isSig <- ifelse(summary_trends$tr_p_value < adj_p, 1, 0)
  summary_trends$display_tr_IRR <- ifelse(summary_trends$isSig == 0, 1,  
                                          summary_trends$tr_p_value)  # plot as 1 if not signif diff
  summary_trends
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
  where cm.sample_date >= ? and cm.sample_date < ? and pc= ?
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
# ---- specific visualisations ----
plot_pc_cats_size_vs_change <- function(db_connection, 
                                        date_end = Sys.Date(),
                                        only_show_significant_trending = FALSE) {
  # plots pc_categories relative to the most common category
  # over the time period analysed
  
  # parameters:
  #   database connection, as returned (for example) by get_db_connection
  #   date_end:  the end of the analysis period
  #   only_show_significant_trending : whether to only show pc_cats with significant trends
  # returns
  #   ggplot object
  # dates may be passed in iso format as strings.
  date_ends <- convert_date(date_end) 
  
  
  summary_trends <- load_trend_summary(
    db_connection = db_connection, 
    date_end = date_ends[['dt']])
  
  summary_trends <- categorise_sample_dates(summary_trends, date_end = date_ends[['dt']])
  if (only_show_significant_trending){
    summary_trends <- subset(summary_trends,  isSig==1)
  }
  summary_trends$plot_tr_IRR <- summary_trends$tr_IRR
  summary_trends$plot_tr_IRR <- ifelse(summary_trends$tr_IRR> 10, 
                                 10, 
                                 summary_trends$tr_IRR)
  summary_trends$plot_tr_IRR <- ifelse(summary_trends$plot_tr_IRR< 0.1, 
                                 0.1, 
                                 summary_trends$plot_tr_IRR)
  
  y_axis_breaks <- c(10,3,2,1.5,1.25,1.1,1, 0.9, 0.8, 0.7, 0.5,0.3)
  p1<-ggplot(summary_trends,aes(x=n_samples,y=plot_tr_IRR)) + 
    geom_point(alpha = 0.3, aes(shape=first_seen_days_ago_cat)) + 
    geom_hline(yintercept =1, lty =1) +
    geom_hline(yintercept =2, lty =2) +
    geom_hline(yintercept =0.5, lty =2) +
    geom_rug(alpha = 0.05)+
        scale_x_continuous("Total sample counts in category",trans='log10') + 
    scale_y_log10("30-day growth rate of samples, relative to the most common \ncategory in each principal component [Rate ratios > 10 plotted at 10]",
                  breaks=y_axis_breaks, labels=y_axis_breaks) +
    scale_shape_discrete(paste('First seen w.r.t ',date_ends[['iso']]))
  p1
}
plot_pc_cats_vs_size <- function(db_connection, 
                                    date_end=Sys.Date(),
                                    only_show_significant_trending=FALSE) {
  # plots a bubble plot illustrating transformed_coordinates(x) vs. pcs (y)
  # illustrating significant rises and falls in recent rates
  
  # parameters:
  #   database connection, as returned (for example) by get_db_connection
  #   date_end:  the end of the analysis period
  #   only_show_significant_trending : whether to only show pc_cats with significant trends
  
  # returns
  #   ggplot object
  # dates may be passed in iso format as strings.
  date_ends <- convert_date(date_end) 
  
  summary_trends <- load_trend_summary(
    db_connection = db_connection, 
    date_end = date_ends[['dt']])
  
  summary_trends <- categorise_sample_dates(summary_trends, date_end = date_ends[['dt']]) 
  if (only_show_significant_trending==TRUE){
    summary_trends <- subset(summary_trends,  isSig==1)
  }
  y_axis_breaks <- c(10,3,2,1.5,1.25,1.1,1.05,1.02,1, 0.98, 0.95, 0.9, 0.8, 0.7, 0.5)

  print("plot_across_pcs: plotting")
  ll <- min(summary_trends$min_trans_coord)
  ul <- max(summary_trends$min_trans_coord)
  
  p2 <- ggplot(summary_trends, aes(x=mean_trans_coord))
  #p2 <- p2 + geom_errorbarh(aes(y=pc, xmin=min_trans_coord, xmax=min_trans_coord), alpha = 0.1)
  p2 <- p2 + geom_point(aes(y=pc, size= n_samples, colour=display_tr_IRR,
                            shape=first_seen_days_ago_cat), alpha = 0.2)
  p2 <- p2 + scale_y_continuous(name='Principal component')  # is is discrete, but plotting as continuous is much faster
  p2 <- p2 + scale_x_continuous("transformed_coordinate (mean, range per cluster)",limits=c(ll,ul))
  p2 <- p2 + scale_size_continuous("Number of samples", trans='log10')
  p2 <- p2 + scale_colour_gradient2(name="Trend in\ncategory", midpoint =1, low = 'blue', high ='red', mid = 'gray')
  p2 <- p2 +  scale_shape_discrete(paste('First seen w.r.t ',date_ends[['iso']]))
  p2
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
  p1 <- p1 + geom_text(aes(x=mean_trans_coord,label=cat, y=2000), data=summary_trends)
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
    params = list(date_end, pc))
  print("Plotting")
  cnts$sample_date <- as.Date(cnts$sample_date)
  cnts$cat <- factor(cnts$cat)
  p1 <- ggplot(cnts, aes(x=sample_date))
  p1 <- p1 + geom_col(aes(y=n_samples, fill=cat),position='stack') 
  p1 <- p1 + scale_x_date("Specimen collection date")
  p1 <- p1 + scale_y_continuous("Number of samples")
  p1 <- p1 + ggtitle(paste0("Principal component ",pc))
  br <- c(10,20,30,50,100,250,500,1000,2500,5000,10000,20000,30000,50000,100000)
  p1 <- p1 + scale_y_sqrt("Number of samples", breaks= br, labels=br)
  p1
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
read_cog_metadata <- function(cogfile){
  cog_metadata <- readr::read_csv(cogfile)
  # extract the country and sequenceId
  cog_metadata$sample_id <- sapply(strsplit(cog_metadata$sequence_name,"/"), `[`, 2)

  # we store dates as isoformat strings.  This is required by sqlite, and is compatible with searching
  # in other database too.
  cog_metadata$sample_date <- lubridate::format_ISO8601(
    cog_metadata$sample_date, 
    usetz = FALSE, 
    precision = 'ymd')

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
                                     "pillar_2", 
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
  # need to encapsulate in an S4 class
  
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
    "transformed_coordinate_categories" = c('pc','cat','sample_id'),
    "sample_id" = c('sample_id')
  )
  
  for (table_name in names(important_indices)) {
    for (index_field in important_indices[[table_name]]){
      
      create_ix_sql = paste0("CREATE INDEX IF NOT EXISTS ", table_name,"_",index_field, " ON ", table_name, "(", index_field, ");")
      print(paste("Ensuring there is an index on table ", table_name, " field ", index_field))
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
add_cog_metadata <- function(db_connection, cogfile){
  # if cog meta data is not already in the database, then add it.
  # parameters:
  #   db_connection: the database connection
  #   cogfile: the cog-uk metadata file.
  
  ## load metadata from cogfile if not present
  tbls <- list_tables(db_connection)
  if (!('sequence_metadata' %in% tbls & 'clinical_metadata' %in% tbls)) {
    # no data loaded from cog metadata file
    print("Loading cogfile metadata")
    metadata_tables <- read_cog_metadata(cogfile)
    add_database_tables(db_conn, metadata_tables)
    
  } else {
    print("sequence metadata is already present.")
  }
}

# ---- linear modelling of counts -----
fit_recent_trend_wrt_most_common_pcat_all_pcs <- function(db_connection,
                                                          interval_analysed=30,
                                                          date_end=NA, 
                                                          analysis_family_id= 'Not provided',
                                                          max_pcs = NA,
                                                          overwrite=FALSE
                                                          ) {
  # fits a Poisson regression model estimating the 
  # proportion of samples of each category of all pcs, 
  # relative to the most common category;
  # as well as the incidence rate ratio (relative growth)
  # over the past interval_analysed days to date_end
  
  # Parameters
  #  this_pc:  the principal component to analyse
  #  interval_analysed: the number of days' data to analyse
  #  date_end: the last day of of the analysis period.  For today, set to 
  #     Sys.Date().  Note that date_end can be a Date, or ISO format 
  #     string representing a date.  If NA, uses the PCA build time 
  #  analysis_family_id : optional, a reference number for this analysis.
  #                       can be used to link together analysis for different 
  #                       pcs in a database
  #  max_pcs: if set, only analyse the first max_pcs.  Useful mainly for debugging.
  #            If not set, all pcs will be analysed.
  #  overwrite: recompute if estimates already present
  #  returns:
  #           Nothing.  Output is written to database.

  if (overwrite==FALSE & db_has_table(db_connection, 'trend_model_metadata')) {
    # already computed
    print("Using stored models.")
    return(0)
  }
  if (is.na(date_end)){
    date_end <- build_date(db_connection)
    print(paste0("Set end of analysis period as PCA build date:", date_end))
  }
  
  date_ends <- convert_date(date_end)
  
  n_pcs <- number_of_pcs(db_connection)
  print(paste("There are", n_pcs, "principal components.  Fitting poisson models of trends."))
  
  # optionally restrict to a smaller number of pcs
  if (!is.na(max_pcs)){
    if (max_pcs < n_pcs) {
    n_pcs <- max_pcs
    print(paste("Restricting analysis to the first ",n_pcs," components."))
    }
  }

  ## fit models
  meta <- list()
  fits <- list()
  
  # iterate over each pc
  for (this_pc in 1:n_pcs) {
    print(this_pc)
    date_end <- lubridate::ymd(db_stem)  # as constructed, the database includes the latest sample date
    res <- fit_recent_trend_wrt_most_common_pcat(db_connection = db_connection,
                                                 this_pc = this_pc, 
                                                 date_end = date_ends[['dt']],
                                                 interval_analysed = interval_analysed,
                                                 analysis_family_id = analysis_family_id)
  
    if (res$fitted == TRUE){
      meta[[this_pc]] <- res$meta
      fits[[this_pc]] <- res$model_fit
    } else {
      warning(paste("Fitting could not occur: ", res$reason))
    }
    
  }
  
  ## store fits to db
  meta_all <- do.call(rbind, meta)
  meta_all$rowid <- 1:nrow(meta_all)
  fits_all <- do.call(rbind, fits)
  fits_all$rowid <- 1:nrow(fits_all)
  if (db_has_table(db_connection, 'trend_model_metadata')) {
    # drop & replace
    print("Dropping & replacing trend_model_metadata")
    dbRemoveTable(db_connection, 'trend_model_metadata')
  }
  dbCreateTable(db_connection, 'trend_model_metadata', meta_all)
  n_added <- dbAppendTable(db_conn, 'trend_model_metadata', meta_all)
  
  if (db_has_table(db_connection, 'trend_model_fits')) {
    # drop & replace
    print("Dropping & replacing trend_model_fits")
    dbRemoveTable(db_connection, 'trend_model_fits')
  }
  dbCreateTable(db_connection, 'trend_model_fits', fits_all)
  n_added <- dbAppendTable(db_connection, 'trend_model_fits', fits_all)
  print("Model fits written to database.")
  return(0)
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
  
  # fit, contrasting to the most common category
  # sample_dow reflects day of week
  # last_week is a fixed effect applied to samples in the last week relative to 
  fit <- glm(n ~ cat * date_t + sample_dow, family="poisson", data=count_df)
  
  # extract coefficients
  coeffs<- tibble(data.frame(summary(fit)$coeff)) %>% 
    mutate(param=rownames(summary(fit)$coeff)) %>%
    filter(!grepl("sample_dow",param)) %>%
    mutate(IRR=exp(Estimate),CI_low=exp(Estimate-1.96*Std..Error), CI_high=exp(Estimate+1.96*Std..Error)) %>%
    rename(p_value=Pr...z..) %>%
    select(param, IRR, CI_low, CI_high,p_value)
  
  coeffs$isRate <- as.integer(grepl(':date_t|date_t',coeffs$param))
  pc_classif<-sub("[^0-9]+([0-9]+)[^0-9]*","\\1",coeffs$param)
  baseline_locs <- grep("date_t|\\(Intercept\\)", pc_classif)
  is_baseline <- rep(FALSE, length(pc_classif))
  is_baseline[baseline_locs] <- TRUE
  pc_classif[baseline_locs] <- levels(count_df$cat)[1]
  coeffs <- coeffs %>%
    mutate(baseline=is_baseline)
  
  intercepts <- coeffs %>% filter(isRate == 0) %>%
    rename(inc_IRR = IRR, inc_CI_high=CI_high, inc_CI_low=CI_low,inc_p_value=p_value)
  slopes <- coeffs %>% filter(isRate == 1) %>%
    rename(tr_IRR = IRR, tr_CI_high=CI_high, tr_CI_low=CI_low,tr_p_value=p_value)
  slopes$fit_cat <- gsub(':date_t','',slopes$param, fixed=TRUE)
  intercepts$fit_cat <- gsub(':date_t','',intercepts$param, fixed=TRUE)
  
  joined_df <- inner_join(intercepts,slopes,by=c("fit_cat")) %>%
    select(-contains("param"),
           -contains("isRate"), 
           -contains('baseline.x'), 
           -contains('baseline.y'))
  
  joined_df <- inner_join(joined_df,total_counts,by=c("fit_cat"))
  joined_df$pc <- this_pc
  joined_df$analysis_id <- analysis_uuid
  
  metadata <- data.frame(list(
    'analysis_family_id' = analysis_family_id,
    'analysis_id' = analysis_uuid,
    'analysis' = 'fit_recent_trend_wrt_most_common_pcat',
    'analysis_type' = 'Poisson regression',
    'readable_info' = 'Estimates relative preponderance of each PC category,\nrelative to the most common category.  Additionally, estimates linear trend in these',
    'pc' = this_pc,
    'date_end' = date_ends[['iso']],
    'interval_analysed_days' = interval_analysed
  ))
  list('fitted' = TRUE, 'reason' = 'Success','meta'= metadata, 'model_fit'= joined_df)
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
                                    overwrite = FALSE){
  # aims to construct a series of 2x2 contingency tables
  #
  #  Feature    Present       Absent
  #  PC-CAT Y      a            c          a+c
  #         N      b            d          b+d
  #               a+b           c+d       a+b+c+d
  
  # Parameters:
  #  db_connection: the database connection
  #  date_end: do not use samples with collection dates later 
  #            than this.  Defaults to today.
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
  
  a_b_c_d_cmd = "select count(*) n from sample_id;"
  
  a_c_cmd = "select sm.feature, count(*) n from 
  sequence_metadata sm
  INNER JOIN sample_id
  ON sm.sample_id = sample_id.sample_id
  INNER JOIN clinical_metadata cm 
  on cm.sample_id = sm.sample_id
  where cm.sample_date <= ?
  group by sm.feature;"
  
  a_b_cmd = "select pc,cat, count(*) n from 
  transformed_coordinate_categories tcc
  INNER JOIN sample_id sm
  ON tcc.sample_id = sm.sample_id
  INNER JOIN clinical_metadata cm 
  on cm.sample_id = sm.sample_id
  where cm.sample_date <= ?
  group by pc, cat;"
  
  a_cmd= "select pc,cat,sm.feature, count(*) n from 
  transformed_coordinate_categories tcc
  INNER JOIN 
  sequence_metadata sm
  on tcc.sample_id = sm.sample_id
  INNER JOIN sample_id
  ON sm.sample_id = sample_id.sample_id
  INNER JOIN clinical_metadata cm 
  on cm.sample_id = sm.sample_id
  where cm.sample_date <= ?
  and pc = ? 
  group by pc,cat, sm.feature;"
  
  
  # check the relevant data exists
  check_table_present(db_connection, 'transformed_coordinate_categories','Table should be present, produced by fn4pca.py')
  check_table_present(db_connection, 'trend_model_fits','Call fit_recent_trend_wrt_most_common_pcat_all_pcs() first')
  check_table_present(db_connection, 'sequence_metadata','Call add_cog_metadata() or equivalent first')
  
  print("Recovering counts")
  a_b_c_d <- database_query(
    db_connection = db_connection, 
    cmd = a_b_c_d_cmd)
  print(paste("number of sequences from which the model was built is ",a_b_c_d))
  
  print("Recovering marginal totals for contingency tables")
  a_c <- database_query(db_connection = db_connection, 
                        cmd = a_c_cmd, 
                        params = list(date_ends[['iso']]))
  names(a_c)[2] <- 'a_c'
  
  print(paste("number of features is ",nrow(a_c)))
  a_b <- database_query(db_connection = db_connection, 
                        cmd = a_b_cmd, 
                        params = list(date_ends[['iso']]))
  names(a_b)[3] <- 'a_b'
  a_b$pc_cat <- paste(a_b$pc, a_b$cat, sep='_')
 
  print(paste("number of pc_cats is ",nrow(a_b)))
  n_pc <- number_of_pcs(db_connection = db_connection)
  print(paste("number of pcs is ",n_pc))
  print("Recovering co-existence of pc_cats and features")
  
  coexistence = list()
  for (this_pc in 1:n_pc){ #DEBUG 3

    print(this_pc)
    coexistence[[this_pc]] <- database_query(
      db_connection = db_connection,
      cmd = a_cmd, 
      params = list(date_ends[['iso']], this_pc)) 
  }
  
  coexistence_df <- do.call(rbind, coexistence)
  names(coexistence_df)[4] <- 'a'
  coexistence_df$pc_cat <- paste(coexistence_df$pc, coexistence_df$cat, sep='_')

  coexistence_df <- merge(coexistence_df, a_c, by = 'feature')

  coexistence_df <- merge(coexistence_df, subset(a_b, select=c('a_b','pc_cat')), by = 'pc_cat')
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
  
  # compute OR for each row.
  coexistence_df$rowid <-1:nrow(coexistence_df)
  ors <- list()
  for (rowid in 1:nrow(coexistence_df)){
    ors[[rowid]] <- compute_or(rowid,
                               coexistence_df$a[rowid],
                               coexistence_df$b[rowid],
                               coexistence_df$c[rowid],
                               coexistence_df$d[rowid])
  }
  ors_df <- plyr::ldply(ors, data.frame)
  coexistence_df <- merge(coexistence_df, ors_df, by='rowid')
  
  # write to database
  add_database_tables(
    db_connection,
    list(feature_associations=list(data_table=coexistence_df,
                                   appropriate_indices=c('pc',
                                                         'pc_cat',
                                                         'feature',
                                                         'p_value')
                               )
    ),
    overwrite=TRUE)
    
  coexistence_df
}

## -------------------------- example usage ------------------------
# set up parameters
BASE_DIR <- "/data/data/pca/subsets_output"  # where the databases are
cogfile <- "/data/data/inputfasta/cog_metadata.csv"
interval <- 30
# find all sqlite dbs
glob_path <- paste0(BASE_DIR,"/0*.sqlite")
for (dbfile in Sys.glob(glob_path)){
    
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
  
  ###  depictions
  # associations with known features
  output_file <- association_heatmap(
    db_connection = db_conn,
    outputdir = PLOT_DIR,
    db_name = db_stem,
    plot_name = 'association_heatmap',
    pdf_size_in_inches = 24)

  # specific PCs
  p1a <- plot_pc_cats_size_vs_change(db_conn)
  to_tiff(to_plot=p1a, 
          outputdir=PLOT_DIR, 
          db_name=db_stem, 
          plot_name='pc_cats_size_vs_change_all')
    
  p1b <- plot_pc_cats_size_vs_change(db_conn,
                                     date_end = db_stem,
                                     only_show_significant_trending=TRUE)
  to_tiff(to_plot=p1b, 
          outputdir=PLOT_DIR, 
          db_name=db_stem, 
          plot_name='pc_cats_size_vs_change_sig_trending')
  
  p1c <- plot_pc_cats_vs_size(db_connection = db_conn,
                              date_end = db_stem)

  to_tiff(to_plot=p1c, 
          outputdir=PLOT_DIR, 
          db_name=db_stem, 
          plot_name='pc_cats_vs_size_all')
  
  p1d <- plot_pc_cats_vs_size(db_connection = db_conn,
                              date_end = db_stem,
                              only_show_significant_trending=TRUE)
  to_tiff(to_plot=p1d, 
          outputdir=PLOT_DIR, 
          db_name=db_stem, 
          plot_name='pc_cats_vs_size_sig_trending')
  
  ### example for a single component
  pc <- 6
  p2 <- plot_single_pc_format_1(
    db_connection = db_conn, 
    pc = pc) 
  to_tiff(to_plot=p2, 
          outputdir=PLOT_DIR, 
          db_name=db_stem, 
          plot_name=paste0('categorised_pc_',pc))

  p3 <- plot_all_pc_format_2(db_conn, 
                             date_end = db_stem)
  to_tiff(to_plot=p3, 
          outputdir=PLOT_DIR, 
          db_name=db_stem, 
          plot_name='trend_all_pcs_scout'
  )
  
  ### over time
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

  
  # all samples, by region
  p4 <- plot_sequenced_numbers_over_time(
    db_connection = db_conn,
    date_end = db_stem) 
  to_tiff(to_plot=p4, 
          outputdir=PLOT_DIR, 
          db_name=db_stem, 
          plot_name='trend_time_region'
  )
  dbDisconnect(db_conn)
  
}
print('FINISHED')
