# depict timing information from demos
rm(list=ls())
library(ggplot2)
library(quantreg)
library(argparse)

parser <- ArgumentParser()
parser$add_argument("inputdir", nargs=1, help="The directory containing the tsv files to read")

args = parser$parse_args()
inputdir <- args$inputdir

# scan the directory looking for tsv files
n_read <- 0
for (this_file in Sys.glob(file.path(inputdir,'*.tsv'))) {
  print(this_file)
  n_read <- n_read +1
  df <- read.table(this_file, sep='\t', header=  TRUE, stringsAsFactors = FALSE)
  if (n_read == 1) {
    timings <- df
  } else {
    timings <- rbind(timings, df)
  }
}

if (n_read == 0) {
  stop(paste("No files found at ",inputdir)) 
}  
# convert to native date formats, except for the first column
for (this_name in c('d_insert','d_read')) {
  timings[[this_name]]<- as.difftime(timings[[this_name]], format="%H:%M:%OS")
}
for (this_name in c('s_insert','e_insert','s_read','e_read')) {
  timings[[this_name]] <- strptime(timings[[this_name]], format="%Y-%m-%d %H:%M:%OS")
}


# quantile regression modelling
f_insert <- rq(d_insert ~ nSamples, data = timings)
timings$pred_d_insert <- predict(f_insert)


# time to add 100k samples
t_i_100 <- predict(f_insert, newdata = data.frame(nSamples = 100000))
# quantile regression modelling
f_read <- rq(d_read ~ nSamples, data = timings)
timings$pred_d_read <- predict(f_read)


p <- ggplot(subset(timings, nSamples>1), aes(x=nSamples))+ theme_classic()
p <- p + geom_point(aes(y=d_insert), size=0.5, colour='red', alpha=0.3)
p <- p + geom_point(aes(y=d_read), size=0.5, colour='blue', alpha=0.3)
p <- p + geom_line(aes(y=pred_d_insert), size=1, colour='red')
p <- p + geom_line(aes(y=pred_d_read), size=1, colour='blue')
p <- p + scale_y_continuous(name='Time to insert one sample(red) or \nread neighbours of 1 sample(blue)')
p <- p + scale_x_continuous(name='Number of samples added')
p <- p + geom_label(aes(x=x, y=y, label=label), data= data.frame(x=mean(timings$nSamples), y=2, label=paste('Est. insert time with 100k samples present:', as.integer(t_i_100), 's')))
ggsave(file=file.path(inputdir, 'insertion_timings.svg'))
