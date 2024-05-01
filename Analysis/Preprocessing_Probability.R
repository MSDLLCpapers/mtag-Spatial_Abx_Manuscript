

# Masha Taguer
# Feb 22 2023

# Get the probability of two bacteria colocalizing with one another given their abundance in the image

# Bring in the raw per image data to get the probability of each pairwise

# Make a df that has each filename and its full path and then add in the metadata by parsing the string of the filename

## Bring data in by FOV
raw.files <- data_frame(filename = list.files('../data/Spatial_and_Abundance/'))
raw.file.paths <- raw.files  %>%
  mutate(filepath = paste0("../data/Spatial_and_Abundance/", filename))

read.csv.and.add.filename <- function(filepath){
  interim1 <- read_csv(filepath) %>%
    mutate(filepath=filepath)  # adds a column with the filepath
  interim1$Cell_Label <- str_c(interim1$Cell_Label, filepath, sep = "_") #add filepath to the cell_label so that rownames are unique
  interim1$Cell_Size <- as.numeric(interim1$Cell_Size) # make all cell sizes as doubles
  interim1$Centroid_X <- as.numeric(interim1$Centroid_X) # make all centroids as doubles
  interim1$Centroid_Y <- as.numeric(interim1$Centroid_Y) # make all centroids as doubles
  interim1$ASV_Cluster_ID <- as.character(interim1$ASV_Cluster_ID) # make all IDs as characters
  interim1$Intensity <- as.numeric(interim1$Intensity)
  as.data.frame(interim1)
}


## for each file, read the csv, and filename as column and to bacteria1, and get triangle

cell.info <- lapply(raw.file.paths$filepath, read.csv.and.add.filename) %>% bind_rows()

#turn bacteria 1 back into just the bacterial id
cell.info$Cell_Label <- sub("\\_.*", "", cell.info$Cell_Label) 


# Add metadata 

## add in column for mouse, location, tissue_replicate, and FOV
cell.info.wide <- cell.info %>% mutate(Mouse = str_sub(filepath, 37, 40)) %>% mutate(Tissue = str_sub(filepath, 42, 46)) %>%
  mutate(Section = str_sub(filepath, 65, 65)) %>% mutate(FoV = str_sub(filepath, 71, 71))

# Get the proportion of each bacteria in each image

cell.info.wide$Sample <- paste(cell.info.wide$Mouse, cell.info.wide$Tissue,
                             cell.info.wide$Section, cell.info.wide$FoV, sep = "_")

cell.info.summary <- cell.info.wide %>% group_by(Sample, ASV_Cluster_ID) %>% dplyr::summarize(Count = n())

cell.info.summary <- cell.info.summary %>% group_by(Sample) %>% mutate(Total = sum(Count), Prop = Count/Total)

## Now we're making the probability wide df to divide the actual values by to generate the score. 
# cell info summary has the proportion of each bacteria in each image
cell.info.summary <- cell.info.summary %>% select(Sample, ASV_Cluster_ID, Prop)

#breakdown by sample
split.cell.info.summary <- split(cell.info.summary, f = cell.info.summary$Sample) 

# if there is only one row, == 1 because that is what the probability is if there is only one type of Bug
not.diag <- map(split.cell.info.summary, ~ 
                if (.x %>% nrow() == 1) {
                  1 
                } else { 
                  (combn(.x$Prop, 2, prod)) 
                }
)

diag <- split.cell.info.summary %>% map(~ .$Prop^2) 
list.of.dfs <- map2(diag, not.diag, make.matrix) %>% lapply(as.data.frame) #now we have all our dfs
# get list of names
row.and.col.names <- split.cell.info.summary %>% map(~ .$ASV_Cluster_ID)

# get names into colnames and rownames
for (i in 1:length(list.of.dfs)){
  colnames(list.of.dfs[[i]]) <- row.and.col.names[[i]]
  row.names(list.of.dfs[[i]]) <- row.and.col.names[[i]]
}

list.long <- list.of.dfs %>% lapply(rownames_to_column, "B1") %>% 
  lapply(pivot_longer, !B1, names_to = "B2", values_to = "Probability")

for (i in 1:length(list.long)){
  list.long[[i]]$Coloc <- paste(list.long[[i]]$B1, list.long[[i]]$B2, sep = "_")
}

# make wide
list.wide <- list.long %>% lapply(select, c(Coloc, Probability)) %>% 
  lapply(pivot_wider, names_from = Coloc, values_from = Probability, names_prefix = "Coloc_")

#Add sample info
for (i in 1:length(list.wide)){
  list.wide[[i]]$Sample <- split.cell.info.summary[[i]]$Sample %>% unique()
}

Probability.wide.0 <- bind_rows(list.wide)
write.csv(cell.info.summary, "../data/cell_info_summary.csv", row.names = F)
write.csv(Probability.wide.0, "../data/Probability_wide_0.csv", row.names = F)

