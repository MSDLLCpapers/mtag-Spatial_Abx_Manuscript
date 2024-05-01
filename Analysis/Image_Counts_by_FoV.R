#Bring in image count data by FoV, not automatically grouped by mouse. 


# Need the replicates for full analysis, so let's bring in the data by FoV and not grouped by mouse

# Make a df that has each filename and its full path and then add in the metadata by parsing the string of the filename

## Bring data in by FOV
raw.files <- data_frame(filename = list.files('../data/Spatial_and_Abundance'))
raw.file.paths <- raw.files  %>%
  mutate(filepath = paste0("../data/Spatial_and_Abundance/", filename))

read.csv.and.add.filename <- function(filepath){
  interim1 <- read_csv(filepath) %>%
    mutate(filepath=filepath)  # adds a column with the filepath
  colnames(interim1)[1] <- "Cell_Label" #change the name of the first column
  interim1$Cell_Label <- str_c(interim1$Cell_Label, filepath, sep = "_") #add filepath to the bacterial name so that rownames are unique
  interim1 <- interim1 %>% column_to_rownames("Cell_Label")
  interim1$Cell_Size <- as.numeric(interim1$Cell_Size)
  interim1$Centroid_X <- as.numeric(interim1$Centroid_X)
  interim1$Centroid_Y <- as.numeric(interim1$Centroid_Y)
  interim1$Intensity <- as.numeric(interim1$Intensity)
  interim1$ASV_Cluster_ID <- as.numeric(interim1$ASV_Cluster_ID)
  interim1
}

## for each file, read the csv, and filename as column and to bacteria1, and get triangle

raw.data <- lapply(raw.file.paths$filepath, read.csv.and.add.filename) %>% bind_rows() %>%
  rownames_to_column("Cell_Label") 

#turn bacteria 1 back into just the bacterial id
raw.data$Cell_Label <- sub("\\_.*", "", raw.data$Cell_Label) 

## add in column for mouse, location, tissue_replicate, and FOV
raw.data.wide <- raw.data %>% mutate(Mouse = str_sub(filepath, 37, 40)) %>% mutate(Tissue = str_sub(filepath, 42, 46)) %>%
  mutate(Section = str_sub(filepath, 65, 65)) %>% mutate(FoV = str_sub(filepath, 71, 71))

# import tissue metadata
md <- read.csv("../data/Tissues_Metadata.csv")
#shorten metadata to only include samples we sent (3 mice per group)
mouse <- c("3402", "3405", "3406" , "3407", "3409", "3412", "3413", "3414", "3415", "3422", "3423", "3424", "3428", "3429", "3430", "3434", "3435", "3436", "3438", "3441", "3442", "3446", "3447", "3448", "3451", "3452", "3454")

# clean up the metadata to only include rows for the mice we submitted to kanvas
md.images <- md[,1:4] %>% filter (Mouse %in% mouse)
md.images$Mouse <- as.character(md.images$Mouse)

Counts.summary <- raw.data.wide %>% group_by(Mouse, Tissue, Section, FoV) %>% summarize(Total = n())

Counts.summary <- left_join(Counts.summary, md.images, by = "Mouse")
