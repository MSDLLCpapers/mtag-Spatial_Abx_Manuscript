
# Metadata

# import tissue metadata
md <- read.csv("../data/Tissues_Metadata.csv")
#shorten metadata to only include samples we sent (3 mice per group)
mouse <- c("3402", "3405", "3406" , "3407", "3409", "3412", "3413", "3414", "3415", "3422", "3423", "3424", "3428", "3429", "3430", "3434", "3435", "3436", "3438", "3441", "3442", "3446", "3447", "3448", "3451", "3452", "3454")

# clean up the metadata to only include rows for the mice we submitted to kanvas
md.images <- md[,1:4] %>% filter (Mouse %in% mouse)


# Import image count data and combine csvs into one dataframe based on mouse and GI location

#CECUM
for (m in mouse) {
  ldf <- list.files(path = "../data/Spatial_and_Abundance", 
                       pattern = paste0(m, "_CECUM_."), full.names = TRUE) %>%
    lapply(read.csv)  #imports each csv
    ldf <- ldf[sapply(ldf, function(x) dim(x)[1]) > 0] #removes empty dfs that were empty csvs that had headers but not bacteria in the image
    cecum <- bind_rows(ldf) #combines df for each mouse into one df
  assign(paste0("cecum", m), cecum) #names that df based on the mouse
}

lcecum <- list(cecum3402, cecum3405, cecum3406 , cecum3407, cecum3409, cecum3412, cecum3413, cecum3414, cecum3415, cecum3422, cecum3423, cecum3424, cecum3428, cecum3429, cecum3430, cecum3434, cecum3435, cecum3436, cecum3438, cecum3441, cecum3442, cecum3446, cecum3447, cecum3448, cecum3451, cecum3452, cecum3454)

# COLON
for (m in mouse) {
  ldf <- list.files(path = "../data/Spatial_and_Abundance", 
                       pattern = paste0(m, "_COLON_."), full.names = TRUE) %>%
    lapply(read.csv)  
    ldf <- ldf[sapply(ldf, function(x) dim(x)[1]) > 0]
    colon <- bind_rows(ldf)
  assign(paste0("colon", m), colon)
} 

lcolon <- list(colon3402, colon3405, colon3406 , colon3407, colon3409, colon3412, colon3413, colon3414, colon3415, colon3422, colon3423, colon3424, colon3428, colon3429, colon3430, colon3434, colon3435, colon3436, colon3438, colon3441, colon3442, colon3446, colon3447, colon3448, colon3451, colon3452, colon3454)


#Get the count and relative abundance per bacterial ID
get_count_RA <- function(df){
  total <- dim(df)[1]
  df1 <- df %>% group_by(ASV_Cluster_ID) %>% summarize(count = n())
  df1 %>% mutate(RA = count/total *100)
}


# CECUM
# for all the cecum dfs, get the relative abundance for each ASV in the sample, and then add what mouse it came from 
lcecum.md <- list()
for (i in 1:length(lcecum)) {
  df.RA <- get_count_RA(lcecum[[i]])  #gets the relative abundance for each dataframe. 
  #Now we need to add a column called Mouse and fill it with the name of the mouse that has the same index as df index
  df.RA$Mouse <- mouse[i]
  lcecum.md[[i]] <- df.RA #add this df to the list of dfs
  lcecum.md.all <- bind_rows(lcecum.md)
} 



#Combine image dfs with metadata

#Combine the imaging dfs with the metadata now that we have the mouse information to join by
#sanity check
lcecum.md.all %>% group_by(Mouse) %>% summarize(sum(RA)) #good.

md.images$Mouse <- as.character(md.images$Mouse)
lcecum.md.all.md <- left_join(lcecum.md.all, md.images, by = "Mouse")

# We have three mice per group so the total relative abundance will be 300 and not 100

# COLON 
# for all the colon dfs, get the relative abundance for each ASV in the sample, and then add what mouse it came from 
lcolon.md <- list()
for (i in 1:length(lcolon)) {
  df.RA <- get_count_RA(lcolon[[i]])  #gets the relative abundance for each dataframe. 
  #Now we need to add a column called Mouse and fill it with the name of the mouse that has the same index as df index
  df.RA$Mouse <- mouse[i]
  lcolon.md[[i]] <- df.RA #add this df to the list of dfs
  lcolon.md.all <- bind_rows(lcolon.md)
}

#Combine the imaging dfs with the metadata now that we have the mouse information to join by
#sanity check
lcolon.md.all %>% group_by(Mouse) %>% summarize(sum(RA)) #good.

md.images$Mouse <- as.character(md.images$Mouse)
lcolon.md.all.md <- left_join(lcolon.md.all, md.images, by = "Mouse")

# We have three mice per group so the total relative abundance will be 300 and not 100
lcecum.md.all.md$Tissue <- "CECUM"
lcolon.md.all.md$Tissue <- "COLON"

summary.images <- rbind(lcecum.md.all.md, lcolon.md.all.md)
