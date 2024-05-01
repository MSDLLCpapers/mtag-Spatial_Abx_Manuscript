# Functions
read.csv.and.add.filename <- function(filepath){
  read_csv(filepath) %>%
    mutate(filepath=filepath) # adds a column with the filepath
}

read.csv.and.add.filename.CB <- function(filepath){
  interim1 <- read_csv(filepath) %>%
    mutate(filepath=filepath)  # adds a column with the filepath
  colnames(interim1)[1] <- "Center_Bug" #change the name of the first column
  interim1
}

KW.function <- function(df) {
  KW.result <- df %>% rstatix::kruskal_test(Score ~ Day) 
  KW.result$Coloc <- df$Coloc %>% unique()
  KW.result
}


KW.effect.function <- function(df) {
  KW.result <- df %>% rstatix::kruskal_effsize(Score ~ Day) 
  KW.result$Coloc <- df$Coloc %>% unique()
  KW.result
}

WC.incr <- function(df)  {
  Dayneg8 <- Score.van.cecum %>% filter(Day == -8, Coloc=="Coloc_1_1") %>% .$Score
  Day0 <- Score.van.cecum %>% filter(Day == 0, Coloc=="Coloc_1_1") %>% .$Score
  aa <- list(Day0 = Day0,Dayneg8 = Dayneg8)
  Day0neg8 <-  data.frame(lapply(aa, "length<-", max(lengths(aa))))
  KW.result <-  Day0neg8 %>% rstatix::wilcox_test(Day0 ~ 1, ref.group = Dayneg8)
  KW.result$Coloc <- df$Coloc %>% unique()
  #add score
  KW.result$ScoreDayneg8 <- mean(Dayneg8)
  KW.result$ScoreDay0 <- mean(Day0)
  KW.result
}

WC.effect.incr <- function(df) {
  Dayneg8 <- Score.van.cecum %>% filter(Day == -8, Coloc=="Coloc_1_1") %>% .$Score
  Day0 <- Score.van.cecum %>% filter(Day == 0, Coloc=="Coloc_1_1") %>% .$Score
  aa <- list(Day0 = Day0,Dayneg8 = Dayneg8)
  Day0neg8 <-  data.frame(lapply(aa, "length<-", max(lengths(aa))))
  KW.result <-  Day0neg8 %>% rstatix::wilcox_effsize(Day0 ~ 1, ref.group = Dayneg8)
  KW.result$Coloc <- df$Coloc %>% unique()
  #add score
  KW.result$ScoreDayneg8 <- mean(Dayneg8)
  KW.result$ScoreDay0 <- mean(Day0)
  KW.result
}

pnova <- function(dist) {
  in.sample <- dist  %>% as.matrix() %>% as.tibble()  %>% colnames()
  metadata.insample <- filter(metadata, Sample %in% in.sample)
  pnova.result <- adonis2(dist ~ Treatment, data = metadata.insample, permutations = 999)
  print(pnova.result)
}

distance.to.pcoa <- function(distance.matrix, negative.eigen, Group){
  # distance.matrix should be the output of vegdist
  # negative.eigen should be a TRUE or FALSE if you are expecting negative eigen values and need to correct
  # Group refers to how you want your points coloured
  pcoa <- cmdscale(distance.matrix, eig = TRUE, add = negative.eigen) 
  positions <- pcoa$points
  colnames(positions) <- c("PCoA1", "PCoA2")
  percent.explained <- pcoa$eig/sum(pcoa$eig) * 100
  pretty <- format(round(percent.explained[1:2], digits = 1), nsmall = 1)
  
  labs <- c(paste("PCoA 1", pretty[1], "%"), paste("PCoA 2", pretty[2], "%"))
  
  positions %>% as.tibble(rownames = "Sample") %>%
    left_join(md.plotting, by = "Sample") %>%
    full_join(md.images, by = "Mouse") %>%
    ggplot(aes_string("PCoA1", "PCoA2", colour = Group)) +
    geom_point() +
    labs(x = labs[1], y = labs[2])
}


# function to make matrix
make.matrix <- function(diag, not.diag) {
  m <- matrix(NA, ncol = length(diag), nrow = length(diag))
  m[lower.tri(m)] <- not.diag
  m[upper.tri(m)] <- t(m)[upper.tri(t(m))]
  diag(m) <- diag
  #take just the triangle
  m[lower.tri(m)] <- NA  #take just the triangle 
  m
}

# function for colocalizations between variables
KW.Day.Comparisons <- function(df, tissue, tx) {
  df <- df %>% filter(Tissue == tissue, Treatment == tx)
  df$Day <- as.factor(df$Day)
  
  # get rid of colocs that never happen and drop rows with NA in score
  Score.df.prob.no0 <- df %>% drop_na(Score) %>% group_by(Coloc) %>% 
    mutate(Total = sum(Score)) %>% filter(Total > 0) %>% ungroup() 
  
  # get rid of a day 0
  Score.df.prob.no0 <- Score.df.prob.no0 %>% filter(Day != 0)
  
  # break down into smaller dfs by Coloc, now it's in a list
  split_data.score <- split(Score.df.prob.no0, f = Score.df.prob.no0$Coloc) 
  
  # KW statistical test, adjusting p value with BH
  Score.df.prob.no0.results <- lapply(split_data.score, KW.function)  %>%
    bind_rows() %>% rstatix::adjust_pvalue(method = "BH") %>% 
    rstatix::add_significance(p.col = "p.adj", symbols = c("****", "***", "**", "*", " "))
  
  # add in tx and tissue information
  Score.df.prob.no0.results$Tissue <- tissue
  Score.df.prob.no0.results$Treatment <- tx
  total.possible <- Score.df.prob.no0.results %>% dim() %>% .[1]
  
  # make a smaller df of only significant
  Score.van.cecum.signif.01 <- Score.df.prob.no0.results %>% filter(p.adj < 0.1)
  signif.05 <- Score.df.prob.no0.results %>% filter(p.adj < 0.05) %>% dim() %>% .[1]
  print(signif.05)
  print(Score.van.cecum.signif.01)
}

# Kruskal wallis effect size function
effect.function <- function(df) {
  effect.size <- df %>% kruskal_effsize(Score ~ Day)
  effect.size$Coloc <- df$Coloc %>% unique()
  effect.size
}

# List of KW results
Effect.size <- function(df, tissue, tx, DayExclude) {
  df <- df %>% filter(Tissue == tissue, Treatment == tx)
  df$Day <- as.factor(df$Day)
  
  # get rid of colocs that never happen and drop rows with NA in score
  Score.df.prob.no0 <- df %>% drop_na(Score) %>% group_by(Coloc) %>% 
    mutate(Total = sum(Score)) %>% filter(Total > 0) %>% ungroup() 
  
  # get rid of a day, based on DayExclude
  Score.df.prob.no0 <- Score.df.prob.no0 %>% filter(Day != DayExclude)
  
  # break down into smaller dfs by Coloc, now it's in a list
  split_data.score <- split(Score.df.prob.no0, f = Score.df.prob.no0$Coloc) 
  
  # KW statistical test, adjusting p value with BH
  Score.df.prob.no0.results <- lapply(split_data.score, effect.function)  %>%
    bind_rows() 
  
  # add in tx and tissue information
  Score.df.prob.no0.results$Tissue <- tissue
  Score.df.prob.no0.results$Treatment <- tx
  Score.df.prob.no0.results
}

KW.Day.Comparisons <- function(df, tissue, tx) {
  df <- df %>% filter(Tissue == tissue, Treatment == tx)
  df$Day <- as.factor(df$Day)
  
  # get rid of colocs that never happen and drop rows with NA in score
  Score.df.prob.no0 <- df %>% drop_na(Score) %>% group_by(Coloc) %>% 
    mutate(Total = sum(Score)) %>% filter(Total > 0) %>% ungroup() 
  
  # get rid of a day 0
  Score.df.prob.no0 <- Score.df.prob.no0 %>% filter(Day != 35)
  
  # break down into smaller dfs by Coloc, now it's in a list
  split_data.score <- split(Score.df.prob.no0, f = Score.df.prob.no0$Coloc) 
  
  # KW statistical test, adjusting p value with BH
  Score.df.prob.no0.results <- lapply(split_data.score, KW.function)  %>%
    bind_rows() %>% rstatix::adjust_pvalue(method = "BH") %>% 
    rstatix::add_significance(p.col = "p.adj", symbols = c("****", "***", "**", "*", " "))
  
  # add in tx and tissue information
  Score.df.prob.no0.results$Tissue <- tissue
  Score.df.prob.no0.results$Treatment <- tx
  total.possible <- Score.df.prob.no0.results %>% dim() %>% .[1]
  
  # make a smaller df of only significant
  Score.van.cecum.signif.01 <- Score.df.prob.no0.results %>% filter(p.adj < 0.1)
  signif.05 <- Score.df.prob.no0.results %>% filter(p.adj < 0.05) %>% dim() %>% .[1]
  print(signif.05)
  print(Score.van.cecum.signif.01)
}

# UMAP PLOT
umap.plot <- function(df, md, Title) {
  # remove the metadata
  df <- df[,-(1:8)] 
  # remove those with count 0
  df.count <- df[, which(colSums(df) != 0)]
  
  umap.counts <- umap(df.count, n_neighbors = 80)
  
  umap_df <- umap.counts$layout %>%
    as.data.frame()%>%
    rename(UMAP1="V1",
           UMAP2="V2") %>%
    mutate(ID = row_number())%>%
    left_join(labels, by="ID")
  
  umap_df %>%
    ggplot(aes(x = UMAP1, y = UMAP2, color = as.character(Day) ))+
    geom_point(size = 2) +
    facet_grid(Treatment~Tissue) +
    stat_ellipse(type = "t") +
    labs(x = "UMAP1",
         y = "UMAP2",
         subtitle = Title) %>% print()
}

#Flatten Correlation Matrixes
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}

# annotating images for figure 4b
print.images <- function(day, Tx, Tis, ID1, ID2) {
  df <- Counts.df[Counts.df$Day == day & Counts.df$Treatment == Tx & Counts.df$Tissue == Tis,]
  df <- df[df$ASV_Cluster_ID == ID1 | df$ASV_Cluster_ID == ID2,] 
  df <- df %>% 
    mutate(Sample = paste(Mouse, Tissue, "tissue_replicate", Section, "fov", FoV, sep = "_")) %>%
    .$Sample
  
  Bacterial_ID <- c(ID1, ID2)
  
  for (sample in df) {
    coords <- read.csv(paste0(csv.filepath, sample, "_cell_information.csv"))
    coords.Am <- coords %>% filter(ASV_Cluster_ID == Bacterial_ID)
    img <- magick::image_read(paste0(image.filepath, sample, "_segmentation.png"))
    img <- image_rotate(img, 270)
    #image_flip(img)
    img.png <- image_convert(img, "png")
    
    theplot <-
      ggplot(coords.Am,
             aes(Centroid_X, Centroid_Y, colour = as.character(ASV_Cluster_ID))) +
      background_image(img.png) +
      geom_point(size = 2, shape = 1)  +
      ggtitle(sample) +
      coord_cartesian(xlim = c(0, 2000), ylim = c(0, 2000), expand = FALSE) +
      theme(
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank() ) +
      geom_circle(aes(x0 = 100, y0 = 1900, r = 71.428), 
                  inherit.aes = FALSE, 
                  color = "red" )
    print(theplot)
  }
}

KW.function.Tx <- function(df) {
  KW.result <- df %>% rstatix::kruskal_test(Score ~ Treatment) 
  KW.result$Coloc <- df$Coloc %>% unique()
  KW.result
}

KW.Tx.Comparisons.rev <- function(df, tissue, tx, day) {
  df <- df %>% filter(Tissue == tissue, Day == day, Treatment == tx | Treatment == "H2O")

  # get rid of colocs that never happen and drop rows with NA in score
  Score.df.prob.no0 <- df %>% drop_na(Score) %>% group_by(Coloc) %>% 
    mutate(Total = sum(Score)) %>% filter(Total > 0) %>% ungroup() 
  
  # break down into smaller dfs by Coloc, now it's in a list
  split_data.score <- split(Score.df.prob.no0, f = Score.df.prob.no0$Coloc) 
  
  # KW statistical test, adjusting p value with BH
  Score.df.prob.no0.results <- lapply(split_data.score, KW.function.Tx)  %>%
  bind_rows() %>% rstatix::adjust_pvalue(method = "BH") %>% 
    rstatix::add_significance(p.col = "p.adj", symbols = c("****", "***", "**", "*", " "))
  print(Score.df.prob.no0.results)

  # add in tx and tissue information
  Score.df.prob.no0.results$Tissue <- tissue
  Score.df.prob.no0.results$Treatment <- tx
  Score.df.prob.no0.results$Day <- day
  total.possible <- Score.df.prob.no0.results %>% dim() %>% .[1]
  
  # make a smaller df of only significant
  Score.van.cecum.signif.01 <- Score.df.prob.no0.results %>% filter(p.adj < 0.1)
  signif.05 <- Score.df.prob.no0.results %>% filter(p.adj < 0.05) %>% dim() %>% .[1]
  print(signif.05)
  print(Score.van.cecum.signif.01)
}

effect.function.tx <- function(df) {
  effect.size <- df %>% kruskal_effsize(Score ~ Treatment)
  effect.size$Coloc <- df$Coloc %>% unique()
  effect.size
}


Effect.size.rev <- function(df, tissue, tx, day) {
  df <- df %>% filter(Tissue == tissue, Day == day, Treatment == tx | Treatment == "H2O")
  df$Treatment <- as.factor(df$Treatment)
  print(df)
  # get rid of colocs that never happen and drop rows with NA in score
  Score.df.prob.no0 <- df %>% drop_na(Score) %>% group_by(Coloc) %>% 
    mutate(Total = sum(Score)) %>% filter(Total > 0) %>% ungroup() 
  
  # break down into smaller dfs by Coloc, now it's in a list
  split_data.score <- split(Score.df.prob.no0, f = Score.df.prob.no0$Coloc) 
  
  # KW statistical test, adjusting p value with BH
  Score.df.prob.no0.results <- lapply(split_data.score, effect.function.tx)  %>%
    bind_rows() 
  print(Score.df.prob.no0.results)
  # add in tx and tissue information
  Score.df.prob.no0.results$Tissue <- tissue
  Score.df.prob.no0.results$Treatment <- tx
  Score.df.prob.no0.results$Day <- day
  Score.df.prob.no0.results
}

get.matrix.rev <- function(df, BUG, Tissue, Abx) {
  # The ASV Cluster ID
  BUGID <- BUG 
  # The tissue group you want to compare
  Tis <- Tissue
  # The treatment group you want to compare
  Tx <- Abx
   
  # Filter  
  OneBug.wide <- df %>%
    filter(Center_Bug == BUGID) %>%
    filter(Tissue == Tis, Treatment == Tx | Treatment == "H2O") %>%
    filter(Day == 35) %>%
    dplyr::select(ASV_Cluster_Scientific_Name, Sample, Affinity_Score) %>%
    pivot_wider(names_from = ASV_Cluster_Scientific_Name, values_from = Affinity_Score, values_fill = 0) %>%
    column_to_rownames("Sample")
   #remove empty dfs
}

get.md.rev <- function(df, BUG, Tissue, Abx) { # get sample metadata 
  # The ASV Cluster ID
  BUGID <- BUG 
  # The tissue group you want to compare
  Tis <- Tissue
  # The treatment group you want to compare
  Tx <- Abx
  OneBug.wide.md <- df %>%
    filter(Center_Bug == BUGID) %>%
    filter(Tissue == Tis, Treatment == Tx | Treatment == "H2O") %>%
    filter(Day == 35) %>%
    dplyr::select(ASV_Cluster_Scientific_Name, Sample, Affinity_Score, Treatment, Tissue, Day, Center_Bug) %>%
    pivot_wider(names_from = ASV_Cluster_Scientific_Name, values_from = Affinity_Score) 
}