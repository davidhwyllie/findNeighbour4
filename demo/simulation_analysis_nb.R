# generate a depiction of the relationship between pc_cats and  tree
# as measured on SARS-CoV-2 genomes.

#library(ggplot2)
rm(list=ls())
library(treeio)
library(ggtree) 
library(ggplot2)
library(ggstance)
library(ape)
library(adephylo)
library(RSQLite)
library(grid)
library(gridExtra)
library(gtable)
library(data.table)
library(svglite)
library(lubridate)
library(trend)
library(argparse)
library(stats)

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
test_monophyletic <- function (phy, tips, reroot = !is.rooted(phy), plot = FALSE, 
                               ...) 
{
  # minor modification of ape is.monophyletic() function
  if (!inherits(phy, "phylo")) 
    stop("object 'phy' is not of class 'phylo'")
  n <- length(phy$tip.label)
  ROOT <- n + 1L
  if (is.numeric(tips)) {
    if (any(tips > n)) 
      stop("incorrect tip#: should not be greater than the number of tips")
    tips <- as.integer(tips)
  }
  if (is.character(tips)) {
    tips <- match(tips, phy$tip.label)
    if (anyNA(tips)) 
      stop("some tip label(s) not found in the tree")
  }
  tips <- sort(tips)
  if (length(tips) == 1L || length(tips) == n) 
    return(list(
      monophyly = TRUE,
      n_tips = length(tips),
      n_descendents = length(tips)
    ))
  if (reroot) {
    outgrp <- phy$tip.label[-tips][1]
    phy <- root(phy, outgroup = outgrp, resolve.root = TRUE)
    rerooted <- TRUE
  }
  else rerooted <- FALSE
  phy <- reorder(phy)
  seq.nod <- .Call(seq_root2tip, phy$edge, n, phy$Nnode)
  sn <- seq.nod[tips]
  newroot <- ROOT
  i <- 2
  repeat {
    x <- unique(unlist(lapply(sn, "[", i)))
    if (length(x) != 1) 
      break
    newroot <- x
    i <- i + 1
  }
  descendents <- which(unlist(lapply(seq.nod, function(x) any(x %in% 
                                                                newroot))))
  if (plot) {
    zoom(phy, tips, subtree = FALSE, ...)
    if (rerooted) 
      mtext("Input tree arbitrarily rerooted", side = 1, 
            cex = 0.9)
  }
  retVal <- list(
    monophyly = identical(tips, descendents),
    n_tips = length(tips),
    n_descendents = length(descendents)
  )
}

# start up ----
parser <- ArgumentParser()
parser$add_argument('-i', "--inputdir", default=file.path('sim_all','sqlite'), type="character",
                    metavar="directory in which to search for sqlite output",
                    help="directory in which to search for sqlite files. Example: /data/pca/sqlite")

parser$add_argument('-f', "--filepattern", default="sim0_*.sqlite", type="character",
                    metavar="pattern to search for in order to identify sqlite files",
                    help="file pattern to glob to look for sqlite files. Example: 2021*.sqlite")

parser$add_argument('-o', "--outputdir", default="output/monophyly", type="character",
                    metavar="output path",
                    help="directory to write output into")


parser$add_argument('-t', "--truthdir", default="truth", type="character",
                    metavar="where the tree files are",
                    help="where the tree files are")


args <- parser$parse_args()

inputdir <- args$inputdir
globpath <- file.path(inputdir,args$filepattern)
outputdir <- args$outputdir
truthdir <- args$truthdir
print("running with parameters:")
print(inputdir)
print(globpath)
print("------------------------")

if (TRUE) {
results <- list()
detections <- list()
ix <- 0
db_id <- 0
for (dbfile in sort(Sys.glob(globpath), decreasing=FALSE)){
  
  print(paste("Reading tree corresponding to ",dbfile))
  sim_id <- basename(dbfile)
  sim_id <- gsub('sim','',sim_id)
  sim_id <- gsub('.sqlite', '', sim_id, fixed=TRUE)
  sim_id <- strsplit(sim_id, '_')[[1]][1]
  tree_file <- file.path(truthdir,paste('truth',sim_id,'.nwk',sep=''))
  print(paste("Treefile is ",tree_file))
  original_t <- treeio::read.tree(tree_file)
  
  original_t$tip.label <- paste0('#',original_t$tip.label)
 
  intertip_file <- file.path(truthdir,paste('intertip_truth',sim_id,'.Rds',sep=''))
  if (!file.exists(intertip_file)){
    print("Computing intertip distances (once only)")
    intertip <- distTips(original_t,method='patristic',useC=TRUE)
    saveRDS(intertip, intertip_file)
  } else {
    print("Using stored intertip distances")
  }
  intertip_m <- as.matrix(readRDS(intertip_file))
  d2r_file <- file.path(truthdir,paste('distance2root_truth',sim_id,'.Rds',sep=''))
  if (!file.exists(d2r_file)){
    print("Computing distances to root (once only)")
    d2r <- distRoot(original_t,method='patristic')
    saveRDS(d2r, d2r_file)
  } else {
    print("Using stored  distance to root")
  }
  d2r <- readRDS(d2r_file)
  db_id <- db_id + 1

  #if (db_id > 150){
  #  break   ## DEBUG
  #}
  print(dbfile)
  analysis_date <- gsub('.sqlite','',basename(dbfile), fixed=TRUE)
  analysis_date <- gsub('sim._','',analysis_date, perl=TRUE)
  analysis_date <- gsub('sim.._','',analysis_date, perl=TRUE) 
  analysis_date <- gsub('sim..._','',analysis_date, perl=TRUE) # up to 999 sims will be coped with
  print(analysis_date)
  analysis_date <- as.Date(analysis_date)
  
  print(paste("Analysis date is ",analysis_date))

  sqlite <- dbDriver("SQLite")
  db_conn <- dbConnect(sqlite,dbfile)

  # load the sigtrends.txt (statistical tests)
  fp_sigtest_sql = "select ps.pop_int_id, ps.level_1_category_type, ps.level_1_category, ps.level_2_category_type, ps.level_2_category, pcas.pcas_int_id, pcas.pc_cat, pcas.earliest_date, pcas.latest_date,pcas.n, smf.statmodel_int_id, smf.estimate, 
  smf.estimate_ci_high, smf.estimate_ci_low, smf.p_value, smf.adj_p_value, sm.analysis_type from population_studied ps inner join pca_summary pcas on pcas.pop_int_id = ps.pop_int_id inner join statistical_model sm on sm.pcas_int_id = pcas.pcas_int_id inner join statistical_model_fit smf
  on smf.statmodel_int_id = sm.statmodel_int_id where smf.param='t' and smf.estimate >0  and level_2_category='lineage_0' and analysis_type = 'GLM:Poisson regression';"
  fp_tests <- database_query(db_conn, fp_sigtest_sql) 
  # load the sigtrends.txt (statistical tests)
  sigtest_sql = "select ps.pop_int_id, ps.level_1_category_type, ps.level_1_category, ps.level_2_category_type, ps.level_2_category, pcas.pcas_int_id, pcas.pc_cat, pcas.earliest_date, pcas.latest_date,pcas.n, smf.statmodel_int_id, smf.estimate, 
  smf.estimate_ci_high, smf.estimate_ci_low, smf.p_value, smf.adj_p_value, sm.analysis_type from population_studied ps inner join pca_summary pcas on pcas.pop_int_id = ps.pop_int_id inner join statistical_model sm on sm.pcas_int_id = pcas.pcas_int_id inner join statistical_model_fit smf
  on smf.statmodel_int_id = sm.statmodel_int_id where smf.param='t' and smf.estimate >0  and level_2_category='UKREGION';"
  tests <- database_query(db_conn, sigtest_sql)
  tests$pc_cat <- as.character(tests$pc_cat)
  nb <- subset(tests, analysis_type == 'GLM:NegativeBinomial', select = c('pc_cat','estimate','estimate_ci_high','estimate_ci_low','p_value'))
  names(nb) <- c('pc_cat','nb_estimate','nb_estimate_ci_high','nb_estimate_ci_low','nb_p_value')
  pois <- subset(tests, analysis_type == 'GLM:Poisson regression')
  both <- merge(pois,nb, all.x=TRUE, by='pc_cat')

  #p <- ggplot(both, aes(x= p_value, y=nb_p_value))
  #p <- p + geom_point()
  #p <- p + scale_y_log10("nb_p_value")
  #p <- p + scale_x_log10("pois_p_value")
  #p 
  #to_svg(p, outputdir, as.character(analysis_date), paste0(as.character(sim_id),'_bivar_p'))
        
 
  fp_tests$adj_p_value <- p.adjust(fp_tests$p_value, method = 'fdr')

  # Bonferroni if n small, otherwise FDR
  p_value_cutoff <- 0.02
  if (nrow(both)<20){
    both$adj_p_value  <- both$p_value
    p_value_cutoff <- p_value_cutoff / nrow(both)
  } else {
    both$adj_p_value <- p.adjust(both$p_value, method = 'fdr')    # poisson model
  }
  
  both <- subset(both, adj_p_value < p_value_cutoff)
  #print(both)
  # export any positive trends
  md_sql = "select * from modelled_data"
  md <- database_query(db_conn, md_sql)

  assocfile <- file.path(inputdir, gsub('.sqlite','.assocs.txt',basename(dbfile)))
 
  signal <- subset(both, level_2_category == 'UKREGION' & n >=3)
  assocs <- read.csv(assocfile, header= TRUE)
  lin1_assocs_10 <- subset(assocs, sequencefeature == 'pangolearn:lineage_1' & log_or > 1 & a > 0 ) # OR at least 10
  signal1 <- subset(signal, 
                      pc_cat %in% lin1_assocs_10$pc_cat ) 

  signal3 <- subset(fp_tests, adj_p_value <= p_value_cutoff) 
 
  #print(signal3)
  # BH method p_value control applied
  #print(nrow(both))
  # print(nrow(subset(both, adj_p_value < 0.01)))
  #print(nrow(subset(both, adj_p_value < 0.003)))
  #print(nrow(subset(both, adj_p_value < 0.001)))
  #print(nrow(signal))
  #print(nrow(signal1))
  #print(nrow(signal3))
  #if (nrow(signal)>0){
    #print(signal)
    #print(lin1_assocs_10)
  #}
  # get the categories

  #print(dbfile)
 
  cmd <- "select * from clinical_metadata;"
  clin_df <- database_query(db_conn, cmd)
  clin_df$sample_date_decimal <- decimal_date(as.Date(clin_df$sample_date))
  
  cmd <- "select * from sequence_feature where sequencefeature == 'pangolearn:lineage_1'";
  seq_df <- database_query(db_conn, cmd)
  
  
  cmd <- "select * from analysed_sample;"
  sample_df <- database_query(db_conn, cmd)
  sample_df$sample_date
  sample_df <- merge(sample_df, clin_df, by='sample_id')
  
  cmd <- "select pcas.* from pca_summary pcas inner join population_studied ps on ps.pop_int_id = pcas.pop_int_id where level_2_category = '--Any--';"
  pcas <- database_query(db_conn, cmd)
  #print(analysis_date)
  small_and_novel <- subset(pcas, n_days_observed < 10 & (analysis_date-30) <= latest_date )
 
  signal2 <- subset(small_and_novel, pc_cat %in% lin1_assocs_10$pc_cat)
 
  cmd <- "select * from transformed_coordinate_category tcc;"
  pca_output <- database_query(db_conn, cmd)
  pca_output <- merge(pca_output, sample_df, by = 'sample_int_id')
  
  pca_output$label <- paste0('#', pca_output$sample_id)
  pca_output$expanding_branch <- ifelse(as.integer(pca_output$sample_id) > 1000, 1,0)
 
  # remove any tips which we have not analysed.
  print("Dropping unused tips")
  t <- original_t
  to_remove <- setdiff(t$tip.label, pca_output$label)
  for (tip in to_remove){
    t <- drop.tip(t,tip)
  }

  # make a node2label lookup, which is required if a timescaled tree is to be constructed
  node2label <- t %>% as_tibble()
  #node2label <- merge(node2label, unique(subset(pca_output, select = c('sample_id','label','sample_date_decimal'))), by='label')
  #reorder_node2label <- data.frame(label = t$tip.label)
  #reorder_node2label <- plyr::join(reorder_node2label, node2label, by='label')
  
  # bactdating example [super slow] the depiction is a side effect
  #q <- roottotip(t, reorder_node2label$sample_date_decimal)
  #bd <- bactdate(unroot(t),reorder_node2label$sample_date_decimal,nbIts=10)
  max_estimate <- 0
  mean_estimate <- 0
  if (nrow(signal)>0){
    max_estimate = max(signal$estimate)
    mean_estimate = mean(signal$estimate)
  }
  res_detection = c("analysis_date_num" = decimal_date(as.Date(analysis_date)), 
                    "n_lineage_1" = nrow(seq_df),
                    "n_total" = nrow(sample_df),
                    "Lineage_1_associations_orgt10" = nrow(lin1_assocs_10),
                    "false_positives" = nrow(signal3),
                    "signal" = nrow(signal),
                    "true_positives_poisson" = nrow(signal1),
                    "true_positives_novelty" = nrow(signal2),
                    "max_estimate" = max_estimate,
                    "mean_estimate" = mean_estimate
  )
  #print(res_detection) 
  #print("-------------- signal 1 --------------------")   
  #print(signal1)
  #print(lin1_assocs_10)   
  detections[[db_id]] = res_detection
  
  for (this_pc in unique(pca_output$pc)){
    pca_output_subset <- subset(pca_output, pc == this_pc)
    #print(paste(analysis_date, "pc=", this_pc, "tree size=", t$Nnode))
    for (this_cat in unique(pca_output_subset$cat)){
      this_pc_cat <- paste(this_pc,this_cat, sep='_')
      pca_output_subset$is_cat <- 'No'
      pca_output_subset$is_cat <- ifelse(pca_output_subset$cat == this_cat, 'Yes', pca_output_subset$is_cat)
      pca_output_subset$is_cat <- ifelse(abs(pca_output_subset$cat-this_cat) == 1, 'Neighbour', pca_output_subset$is_cat)
      
      # find the node ids associated with these pcs
      tips <- pca_output_subset$label[pca_output_subset$cat==this_cat]
      ed <- pcas$earliest_date[pcas$pc_cat==this_pc_cat]
      ld <- pcas$latest_date[pcas$pc_cat==this_pc_cat]
      n <- pcas$n[pcas$pc_cat==this_pc_cat]
      nodes <- node2label$node[node2label$label %in% tips]
      #print(paste(this_pc_cat, 'with',length(nodes),'nodes'))
      monophyly_summary <-  test_monophyletic(t,nodes)
      monophyly <- monophyly_summary$monophyly
      prop_descendents_in_pc_cat <- monophyly_summary$n_tips / 
                                    monophyly_summary$n_descendents
     
      dist2root <- d2r[tips]
      # if includes root, we code min distance to root as zero
      dist2root <- dist2root[!is.na(dist2root)]
      if (length(dist2root)==0){
        dist2root = c(0)
      }
      intertip_subset <- intertip_m[tips,tips]
      ix <- ix +1
      results[[ix]] <- c(
                     analysis_id = db_id,
                     tree_size = t$Nnode,
                     pc=this_pc, 
                     cat=this_cat, 
                     n=length(nodes), 
                     dist2root_min = min(dist2root),
                     dist2root_mean = mean(dist2root),
                     dist2root_sd = sd(dist2root),
                     earliest_date_num = decimal_date(as.Date(ed)),
                     latest_date_num = decimal_date(as.Date(ld)),
                     intertip_mean = mean(intertip_subset),
                     intertip_median = median(intertip_subset),
                     intertip_sd = sd(intertip_subset),
                     monophyly = monophyly,
                     prop_descendents_in_pc_cat = prop_descendents_in_pc_cat
      
                     )
      # produce illustrations of phyly for one tree
      if (db_id == 1 & this_pc %in% c(3)){
        tip_meta <- data.frame(label = t$tip.label)  # 1st col
        tip_meta <- merge(
          tip_meta, 
          subset(
            pca_output,
            pc==this_pc,
            select = c('label','transformed_coordinate','sample_id','expanding_branch')
          ),
          by='label') # ensures 1st col tip label
        tip_meta$this_cat <- ifelse(tip_meta$label %in% tips, 1, 0)
        p <- ggtree(t, layout='rect', aes(), right=TRUE, ladderize=TRUE) %<+% tip_meta
        p <- p + geom_tippoint(aes(size= as.character(this_cat), colour=as.character(this_cat))) #expanding_branch
        p <- p + scale_colour_manual(values=c('0'='black','1'='red'))
        p <- p + scale_size_manual(values=c('0'=1,'1'=2))
        title_text <- paste(this_pc_cat, "(Observed",ed,"-",ld,'n=',n,') Prop. MRCA descendents in ',this_pc_cat,"=", as.integer(100*prop_descendents_in_pc_cat),'%')
        p <- p + ggtitle(title_text)
        p <- p + theme(legend.position = 'none')
        
        to_svg(p, outputdir, as.character(analysis_date), paste0(as.character(sim_id),'_monophyly_',db_id,'_',this_pc_cat))
        
      }
    }
  }
  print("Disconnected.")
  dbDisconnect(db_conn)
  
 
  # produce a depiction of cutting the tree for one tree
  if (db_id == 1) {
  
    res_df <- as.data.table(do.call(rbind, results))
    res_df <- subset(res_df, analysis_id == db_id)
    res_df$pc_cat <- paste(res_df$pc, res_df$cat, sep='_')
    res_df[,pc:=NULL]
    res_df[,cat:=NULL]  # drop columns
    res_df$prop_descendents_in_pc_cat_cat <- cut(res_df$prop_descendents_in_pc_cat, 
                                                 breaks=c(0,0.3,0.5,0.7,0.8,0.9,0.95,0.99,1))
    pcas2 <- merge(pcas, res_df, by=c('pc_cat'))
    
    first_pcs <- 4
    p <- ggplot(subset(pca_output, pc< first_pcs), aes(x=transformed_coordinate)) + theme_classic()
    p <- p + geom_vline(aes(xintercept= trans_coord_min), data= subset(pcas2, pc< first_pcs), alpha =0.3)
    
    p <- p + geom_histogram(binwidth = 0.05)
    p <- p + facet_grid(.~pc)
    p <- p + scale_x_continuous(name = 'transformed coordinate', limits = c(-2.2,2.2))
    p <- p + scale_y_continuous('Number of samples')
    p0 <- p + ggtitle(analysis_date)
    
    p <- ggplot(subset(pcas2, pc<first_pcs), aes(x=trans_coord_avg)) + theme_classic()
    p <- p + geom_segment(aes(y=as.Date(earliest_date), yend=as.Date(latest_date), x=trans_coord_avg, xend = trans_coord_avg))
    p <- p + facet_grid(.~pc)
    p <- p + scale_y_date(name = "Date range of samples observed")
    p <- p + scale_x_continuous(name = 'transformed coordinate', limits = c(-2.2,2.2))
    p <- p + theme(legend.position='none')
    p1 <- p
    p1
    
    col1 <- list(p0)
    col2 <- list(p1)
    gcol1 <- lapply(col1, ggplotGrob)
    gcol2 <- lapply(col2, ggplotGrob)
    g1 = do.call(cbind, c(gcol1, size="first"))
    g2 = do.call(cbind, c(gcol2, size="first"))
    g4 = do.call(rbind, c(list(g1,g2), size="first"))


    outputfile <- file.path(outputdir, paste0(as.character(sim_id),'Figure2_cut_tcc_example_nb_',as.character(analysis_date), "_", db_id,'.svg'))
    print(paste("Wrote Fig 2 output file to ",outputfile))
    svg(outputfile, width= 7, height = 7)
    grid.newpage()
    grid.draw(g4)    
    dev.off()
  }

  
}

res_df <- as.data.table(do.call(rbind, results))
res_df$tree_size_label = paste("Tree size = ", res_df$tree_size)
res_df$monophyly_cat <- ifelse(res_df$monophyly==0, 'No','Yes')
res_df$prop_descendents_in_pc_cat_cat <- cut(res_df$prop_descendents_in_pc_cat, 
                                             breaks=c(0,0.3,0.5,0.7,0.8,0.9,0.95,0.99,1))

# store the data frame
outputfile <- file.path(outputdir, paste0(sim_id, '_Figure3_TCC_characteristics_nb.Rds'))
print(paste("Wrote data file to ",outputfile))
saveRDS(res_df, outputfile)

outputfile <- file.path(outputdir, paste0(sim_id, '_Figure3_TCC_characteristics_nb.csv'))
write.csv(res_df, outputfile, row.names = FALSE)

outputfile <- file.path(outputdir, paste0(sim_id, '_Figure4_Detections_nb.Rds'))
detect_df <- as.data.table(do.call(rbind, detections))
print(paste("Wrote data file to ",outputfile))
saveRDS(detect_df, outputfile)


outputfile <- file.path(outputdir, paste0(sim_id, '_Figure4_Detections_nb.csv'))
write.csv(detect_df, outputfile, row.names = FALSE)

}

outputfile <- file.path(outputdir, paste0(sim_id, '_Figure3_TCC_characteristics_nb.Rds'))
res_df <- readRDS(outputfile)

# 251, 501, 750,2  , tree_size %in% c(931  )
max_tree_size = max(res_df$tree_size)

outputfile <- file.path(outputdir, paste0(sim_id, '_Figure4_Detections_nb.Rds'))
detect_df <- readRDS(outputfile)

# summary plots ----

p <- ggplot(subset(res_df, tree_size == max_tree_size), aes(x=earliest_date_num, y=intertip_median))+theme_classic()
p <- p + geom_point(aes(colour=prop_descendents_in_pc_cat_cat, size = n))
p <- p + facet_grid(.~tree_size_label)
p <- p + scale_colour_viridis_d(name="Proportion descended\nfrom MRCA in this\ntransformed coordinate category", 
                                direction=-1,
                                begin = 0, end = 1)
p <- p + scale_size_continuous(name="Number of samples\nin transformed coordinate category",
                               trans='log10', range = c(0.3, 3),
                               breaks=c(1,10,100,1000))
#p <- p  +scale_shape_discrete("Monophyletic group")
p <- p + scale_y_continuous("Median distance between members\nof the transformed coordinate category")
p0 <- p + scale_x_continuous("Earliest sample observed", limits=c(2020,2021.5))
p0

p <- ggplot(subset(res_df, tree_size == max_tree_size),
            aes(x=latest_date_num, y=intertip_median))+theme_classic()
p <- p + geom_point(aes(colour=prop_descendents_in_pc_cat_cat, 
                        size = n
)
)
p <- p + facet_grid(.~tree_size_label)
p <- p + scale_colour_viridis_d(name="Proportion descended\nfrom MRCA in this\ntransformed coordinate category", 
                                direction=-1,
                                begin = 0, end = 1)
p <- p + scale_size_continuous(name="Number of samples\nin transformed coordinate category",
                               trans='log10', range = c(0.3,3),
                               breaks=c(1,10,100,1000))
#p <- p  +scale_shape_discrete("Monophyletic group")
p <- p + scale_y_continuous("Median distance between members\n of the transformed coordinate category")
p1 <- p + scale_x_continuous("Latest sample observed", limits=c(2020,2021.5))
p1

col1 <- list(p0)
col2 <- list(p1)
gcol1 <- lapply(col1, ggplotGrob)
gcol2 <- lapply(col2, ggplotGrob)
g1 = do.call(cbind, c(gcol1, size="first"))
g2 = do.call(cbind, c(gcol2, size="first"))
g4 = do.call(rbind, c(list(g1,g2), size="first"))


outputfile <- file.path(outputdir, paste0(sim_id, '_Figure3_TCC_characteristics_nb.svg'))
print(paste("Wrote Fig 3 output file to ",outputfile))

svg(outputfile, width= 10, height = 10)
grid.newpage()
grid.draw(g4)    
dev.off()

detect_df$analysis_date <- date_decimal(detect_df$analysis_date_num)
detect_df$analysis_date_d <- as.Date(detect_df$analysis_date)

p <- ggplot(detect_df,
            aes(x=analysis_date_d))+theme_classic()
p <- p + geom_path(aes(y=n_total-n_lineage_1+1))
p <- p + geom_point(aes(y=n_total-n_lineage_1+1))

p <- p + geom_path(aes(y=n_lineage_1+1), colour = 'blue')
p <- p + geom_point(aes(y=n_lineage_1+1),colour = 'blue')
breaks = c(1,2,3,5,10, 20,30,50,75,100,200,300,500)
p <- p + scale_y_log10("Number of samples + 1", breaks = breaks, labels = breaks)
p0 <- p + scale_x_date("Date")
p0

p <- ggplot(detect_df,
            aes(x=analysis_date_d))+theme_classic()
p <- p + geom_vline(aes(xintercept = analysis_date_d), data = subset(detect_df, true_positives_novelty > 0), size=1, colour = 'yellow')
#p <- p + geom_vline(aes(xintercept = analysis_date_d), data = subset(detect_df, signal > 0), size=1, colour = 'red')
p <- p + geom_vline(aes(xintercept = analysis_date_d), data = subset(detect_df, true_positives_poisson > 0), size=1, colour = 'blue')

#p <- p + geom_vline(aes(xintercept = analysis_date_d), data = subset(detect_df, false_positives > 0), size=3, colour = 'lightgray')

p <- p + geom_path(aes(y=Lineage_1_associations_orgt10))
p <- p + geom_point(aes(y=Lineage_1_associations_orgt10))

p <- p + scale_y_continuous("Number of expanding clone associated pc_cats")

p1 <- p + scale_x_date("Date")
p1

col1 <- list(p0)
col2 <- list(p1)
gcol1 <- lapply(col1, ggplotGrob)
gcol2 <- lapply(col2, ggplotGrob)
g1 = do.call(cbind, c(gcol1, size="first"))
g2 = do.call(cbind, c(gcol2, size="first"))
g4 = do.call(rbind, c(list(g1,g2), size="first"))

outputfile <- file.path(outputdir, paste0(sim_id, '_Figure4_Detections_nb.svg'))
print(paste("Wrote Fig 4 output file to ",outputfile))

svg(outputfile, width= 10, height = 10)
grid.newpage()
grid.draw(g4)    
dev.off()

