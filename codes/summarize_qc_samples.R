setwd('/home/i439h/projects/hipo_temp/results/cellranger/mouse/')
library(stringr)

sum_files = list.files(pattern = 'metrics_summary.csv',all.files = T,recursive = T,full.names = T)
sum_files= sum_files[grep('outs',sum_files)]
extracted_parts <- regmatches(sum_files, regexpr("K35R-[^/]+", sum_files))


dat.list = list()

for(i in seq_along(sum_files)){
  dat = read.csv(sum_files[i])
  dat$sample = extracted_parts[i]
  dat.list[[i]] = dat 
}

alldata = do.call('rbind',dat.list)
write.table(alldata,file = "summary_metrics_samples.tsv",quote = F)
