setwd("/home/i439h/projects/hipo_temp/")

library(stringr)

dir = "/home/i439h/projects/hipo_temp/results/cellranger/PDX/"

files = list.files(path = dir,pattern = 'metrics_summary.csv',full.names = T,recursive = T)
files = files[grep('outs',files)]

data_list = list()
# Loop over each file
for (file in files) {
  tryCatch({
    # Read the CSV file
    dat <- read.csv(file)
    
    # Extract sample information
    dat$sample <- str_extract(file, "K35R[^/]+")
    
    # Extract species information
    dat$species <- str_extract(file, "(human|mouse|PDX)")
    
    # Append the data frame to the list
    data_list[[length(data_list) + 1]] <- dat
  }, error = function(e) {
    # Print the error message and continue to the next file
    message("Error processing file: ", file)
    message("Error message: ", e$message)
    next
  })
}

# Combine all data frames in the list into a single data frame
if (length(data_list) > 0) {
  # Find the union of all column names
  all_columns <- unique(unlist(lapply(data_list, colnames)))
  
  # Ensure all data frames have the same columns
  data_list <- lapply(data_list, function(df) {
    df[setdiff(all_columns, colnames(df))] <- NA  # Add missing columns with NA values
    return(df[all_columns])  # Order columns uniformly
  })
  
  # Bind all data frames together
  alldata <- do.call(rbind, data_list)
} else {
  alldata <- data.frame()  # Empty data frame if no valid data
}

gsub(",","",alldata$Estimated.Number.of.Cells) %>% as.numeric()
# Print the combinedEstimated.Number.of.Cells# Print the combined data
print(alldata)

for( i in 1:4){
  alldata[,i] =readr::parse_number(alldata[,i])
}

## some plots 
(p1 = ggplot(alldata, aes(y = sample, x = Estimated.Number.of.Cells,fill=species))+
  geom_bar(stat = 'identity',color = 'black',width = 0.5)+
  facet_wrap(~species,scales = 'free')+
  geom_vline(xintercept = 200)+
  scale_x_log10()+
  scale_x_continuous(limits = c(0,20000))+
  scale_fill_brewer(palette = 'Accent')+
  ggthemes::theme_calc()+
  ggtitle("Estimated number of cells"))


(p2 = ggplot(alldata, aes(y = sample, x =Mean.Reads.per.Cell,fill=species))+
  geom_bar(stat = 'identity',color = 'black',width = 0.5)+
  facet_wrap(~species,scales = 'free')+
  # scale_x_continuous(limits = c(0,40000))+
  # scale_x_log10()+
  scale_fill_brewer(palette = 'Accent')+
  ggthemes::theme_calc()+
    ggtitle("Mean Reads per Cell"))

subdat = subset(alldata,species != 'PDX')
(p3 = ggplot(subdat, aes(y = sample, x =Median.Genes.per.Cell,fill=species))+
    geom_bar(stat = 'identity',color = 'black',width = 0.5)+
    facet_wrap(~species,scales = 'free')+
    # scale_x_continuous(limits = c(0,40000))+
    # scale_x_log10()+
    scale_fill_brewer(palette = 'Accent')+
    ggthemes::theme_calc()+
    ggtitle("Median Genes per Cell"))

write.table(alldata,"./results/cellranger/PDX/summary_metrics_samples.tsv",quote = F,col.names = T)
