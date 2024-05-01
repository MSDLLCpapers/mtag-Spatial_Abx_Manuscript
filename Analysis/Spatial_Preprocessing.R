# Preprocessing


# Make a df that has each filename and its full path and then add in the metadata by parsing the string of the filename

## Bring data in by FOV
raw.files <- data_frame(filename = list.files('../data/Milestone_2_Spatial_Association_Data_Package/Sample/'))
raw.file.paths <- raw.files  %>%
  mutate(filepath = paste0("../data/Milestone_2_Spatial_Association_Data_Package/Sample/", filename))

read.csv.and.add.filename <- function(filepath){
  interim1 <- read_csv(filepath) %>%
    mutate(filepath=filepath)  # adds a column with the filepath
  colnames(interim1)[1] <- "B1" #change the name of the first column
  interim1$B1 <- str_c(interim1$B1, filepath, sep = "_") #add filepath to the bacterial name so that rownames are unique
  interim1 <- interim1 %>% column_to_rownames("B1")
  interim1[lower.tri(interim1)] <- NA  #take just the triangle 
  as.data.frame(interim1)
}

## for each file, read the csv, and filename as column and to bacteria1, and get triangle

raw.data <- lapply(raw.file.paths$filepath, read.csv.and.add.filename) %>% bind_rows() %>%
  rownames_to_column("B1") 

#turn bacteria 1 back into just the bacterial id
raw.data$B1 <- sub("\\_.*", "", raw.data$B1) 

## add in column for mouse, location, tissue_replicate, and FOV
raw.data.wide <- raw.data %>% mutate(Mouse = str_sub(filepath, 73, 76)) %>% mutate(Tissue = str_sub(filepath, 78, 82)) %>%
  mutate(Section = str_sub(filepath, 101, 101)) %>% mutate(FoV = str_sub(filepath, 107, 107))


# make it long
## change first column to Bacteria 1 and as character
colnames(raw.data.wide)[1] <- "Bacteria1"
raw.data.wide$Bacteria1 <- as.character(raw.data.wide$Bacteria1)

data.long <- raw.data.wide %>% pivot_longer(where(is.numeric), names_to = "Bacteria2", values_to = "Count")

### sanity check
data.long$Mouse %>% unique()
data.long$Tissue %>% unique()
data.long$FoV %>% unique()
data.long$Section %>% unique()

# Add in the metadata

## import tissue metadata
md <- read.csv("../data/Tissues_Metadata.csv")
## shorten metadata to only include samples we sent (3 mice per group)
mouse <- c("3402", "3405", "3406" , "3407", "3409", "3412", "3413", "3414", "3415", "3422", "3423", "3424", "3428", "3429", "3430", "3434", "3435", "3436", "3438", "3441", "3442", "3446", "3447", "3448", "3451", "3452", "3454")
## clean up the metadata to only include rows for the mice we submitted to kanvas
md.images <- md[,1:4] %>% filter (Mouse %in% mouse)

## Add in metadata of the sample, ie treatment and timepoint
md.images$Mouse <- as.character(md.images$Mouse)
data.long.md <- left_join(data.long, md.images, by = "Mouse")

# Bring in the randomized data

# Make a df that has each filename and its full path and then add in the metadata by parsing the string of the filename

## Bring data in by FOV
raw.files <- data_frame(filename = list.files('../data/Milestone_2_Spatial_Association_Data_Package/Random/'))
raw.file.paths <- raw.files  %>%
  mutate(filepath = paste0("../data/Milestone_2_Spatial_Association_Data_Package/Random/", filename))

## bring in each csv and add filename and get just the triangle
raw.data.random <- lapply(raw.file.paths$filepath, read.csv.and.add.filename) %>% bind_rows() %>%
  rownames_to_column("B1") 

#turn bacteria 1 back into just the bacterial id
raw.data.random$B1 <- sub("\\_.*", "", raw.data.random$B1) 

## add in column for mouse, location, tissue_replicate, and FOV
raw.data.wide.random <- raw.data.random %>% mutate(Mouse = str_sub(filepath, 73, 76)) %>% mutate(Tissue = str_sub(filepath, 78, 82)) %>%
  mutate(Section = str_sub(filepath, 101, 101)) %>% mutate(FoV = str_sub(filepath, 107, 107))

# make it long
## change first column to Bacteria 1 and as character
colnames(raw.data.wide.random)[1] <- "Bacteria1"
raw.data.wide.random$Bacteria1 <- as.character(raw.data.wide.random$Bacteria1)

data.long.random <- raw.data.wide.random %>% pivot_longer(where(is.numeric), names_to = "Bacteria2", values_to = "Count")

### sanity check
data.long.random$Mouse %>% unique()
data.long.random$Tissue %>% unique()
data.long.random$FoV %>% unique()
data.long.random$Section %>% unique()

# Add in the metadata

## import tissue metadata
md <- read.csv("../data/Tissues_Metadata.csv")
## shorten metadata to only include samples we sent (3 mice per group)
mouse <- c("3402", "3405", "3406" , "3407", "3409", "3412", "3413", "3414", "3415", "3422", "3423", "3424", "3428", "3429", "3430", "3434", "3435", "3436", "3438", "3441", "3442", "3446", "3447", "3448", "3451", "3452", "3454")
## clean up the metadata to only include rows for the mice we submitted to kanvas
md.images <- md[,1:4] %>% filter (Mouse %in% mouse)

## Add in metadata of the sample, ie treatment and timepoint
md.images$Mouse <- as.character(md.images$Mouse)
random.long.md <- left_join(data.long.random, md.images, by = "Mouse")
write.csv(data.long.md, "../data/data_long_md.csv", row.names = F)
write.csv(random.long.md, "../data/random_long_md.csv", row.names = F)

