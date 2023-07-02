### Imputation paper analysis 

library(dplyr)
library(nationalparkcolors)
library(ggplot2)

################ Concordance Stats ####################
### Get file names --------------
files<-do.call(c, lapply(c("filtered", "unfiltered"), function(i){
    files_tmp <- list.files(path = paste0("/concord/", i),
                            pattern="*.concord.txt")
    return(paste0("/concord/", i, "/", files_tmp))
}))
files_GCsAF <- files[str_detect(files, "GCsAF")]
files_GCTs <- files[str_detect(files, "GCTs")]
#
### Read in GCsAFs ------------
afs_r = do.call(rbind, lapply(files_GCsAF, function(x){
  tmp<-read.delim(file = x, 
                  sep = '\t', header = T, skip = 1)
  # set a new column as the file name
  tmp$file_name <- gsub(x, pattern = ".concord.txt", replacement = "")
  # extract information from file name 
  tmp$eff_pop <- strsplit(tmp$file_name, split = "/")[[1]][1]
  tmp$status <- strsplit(tmp$file_name, split = "/")[[1]][2]
  tmp$animal_id <- strsplit(strsplit(tmp$file_name, split = "/")[[1]][3], "_")[[1]][2]
  tmp$coverage <- strsplit(strsplit(tmp$file_name, split = "/")[[1]][3], "_")[[1]][3]
  tmp$file_name <- NULL
  return(tmp)
}))
# Change column names
afs_r<-afs_r[,-2]
colnames(afs_r)[1:10] <- c("GCsAF","AF", "RR.Hom.matches", "RA.Het.matches", "AA.Hom.matches", "RR.Hom.mismatches", "RA.Het.mismatches", "AA.Hom.mismatches", "r.squared", "n_genotypes")
afs_r <- afs_r %>% 
  mutate(AF = as.factor(ifelse(AF==0.005, "0.001-1%", 
                               ifelse(AF==0.030, "1-5%", 
                                      ifelse(AF==0.075, "5-10%", "10-50%"))))) %>% 
  mutate(AF=fct_relevel(AF,c("0.001-1%","1-5%","5-10%", "10-50%"))) %>% 
  mutate(coverage = fct_relevel(coverage, "0.1x", "0.5x", "1x", "3x", "10x"))

# Make wide df 
afs_r_filt <- afs_r[afs_r$status=="filtered",-1]
colnames(afs_r_filt)[c(2:9)] <- paste0(colnames(afs_r_filt)[c(2:9)], "_filt")
afs_r_wide<-left_join(afs_r_filt, 
                    afs_r[afs_r$status=="unfiltered",-1], 
                    by = c("coverage","AF","animal_id","eff_pop"))
afs_r_wide$status.x<-NULL;afs_r_wide$status.y<-NULL
head(afs_r_wide,3)
#
### Number/percent of sites imputed to high-confidence -----------
afs_r_wide %>% 
  ggplot(aes(x = coverage, y = (n_genotypes_filt/n_genotypes)*100, color=AF)) + 
  geom_boxplot(outlier.colour = NA) +
  geom_point(size=2, shape = 1, position=position_jitterdodge()) + 
  theme_classic() + 
  scale_y_continuous(limits=c(82,100)) + 
  labs(y="percent genotypes pass filter (%)", color="MAF") + 
  scale_color_manual(values = brewer.pal(4,name = "Dark2")[c(4,1:3)]) 
#
### Factors that imputation accuracy ----------------------
afs_r %>% 
  filter(status == "filtered") %>%
  ggplot(aes(x = coverage, y = r.squared, color = AF)) + 
  geom_boxplot(outlier.colour = NA) + 
  geom_point(size=2, shape=1, position=position_jitterdodge()) + 
  theme_classic() + 
  labs(y = bquote(bold(r^2)), color="MAF") + 
  scale_color_manual(values = brewer.pal(4,name = "Dark2")[c(4,1:3)]) + 
  scale_y_continuous(limits=c(0.8,1),breaks=c(0.8,0.9,1))
#
### Read in GCTs ----------------------
gcts_ind_r = do.call(rbind, lapply(files_GCTs, function(x){
  # read files in
  tmp<-read.delim(file = x, sep = '\t', header = T, skip = 1)
  # set a new column as the file name
  tmp$file_name <- gsub(x, pattern = ".concord.txt", replacement = "")
  tmp$file_name <- gsub(tmp$file_name, pattern = "/concord/concord_", replacement = "")
  # extract information from file name 
  tmp$eff_pop <- strsplit(tmp$file_name, split = "/")[[1]][1]
  tmp$status <- strsplit(tmp$file_name, split = "/")[[1]][2]
  tmp$animal_id <- strsplit(strsplit(tmp$file_name, split = "/")[[1]][3], "_")[[1]][2]
  tmp$coverage <- strsplit(strsplit(tmp$file_name, split = "/")[[1]][3], "_")[[1]][3]
  tmp$file_name <- NULL
  tmp$p_concord_file <- sum(tmp[c(3,9,15)])/sum(tmp[c(3:5,8:10,13:15)])
  return(tmp)
}))

# Change column names
colnames(gcts_ind_r)[1:27] <- c("GCTs","sample", "RRHomRRHom", "RRHomRAHet", "RRHomAAHom", "RRHomAAHet", "RRHomMissing", "RAHetRRHom", "RAHetRAHet", "RAHetAAHom", "RAHetAAHet",	"RAHetMissing", "AAHomRRHom", "AAHomRAHet",	"AAHomAAHom",	"AAHomAAHet",	"AAHomMissing",	"AAHetRRHom", "AAHetRAHet",	"AAHetAAHom",	"AAHetAAHet", "AAHetMissing", "MissingRRHom", "MissingRAHet", "MissingAAHom", "MissingAAHet", "MissingMissing")
gcts_ind_r$coverage <- factor(gcts_ind_r$coverage, levels = c("0.1x", "0.5x", "1x", "3x","10x"))

for (i in c("RRHomAAHet", "RAHetAAHet", "AAHomAAHet", 
            "AAHetRRHom", "AAHetRAHet", "AAHetAAHom", "AAHetAAHet", "AAHetMissing", 
            "MissingRRHom", "MissingRAHet", "MissingAAHom", "MissingAAHet", "MissingMissing")) {
  print(range(gcts_ind_r[[i]]))
}

for (i in c("RRHomAAHet", "RAHetAAHet", "AAHomAAHet", 
            "AAHetRRHom", "AAHetRAHet", "AAHetAAHom", "AAHetAAHet", "AAHetMissing", 
            "MissingRRHom", "MissingRAHet", "MissingAAHom", "MissingAAHet", "MissingMissing")) {
  gcts_ind_r[[i]] <- NULL
}

# sum sites 
gcts_ind_r$total_sites = apply(gcts_ind_r, 1, function(x) {sum(as.numeric(x[c(3:5,7:9,11:13)]))})
#
### Heatmap error type ---------------------
error_heat_r <- gcts_ind_r %>% 
  mutate(perc_RRHomRRHom = RRHomRRHom/(RRHomRRHom+RRHomRAHet+RRHomAAHom),
         perc_RRHomRAHet = RRHomRAHet/(RRHomRRHom+RRHomRAHet+RRHomAAHom),
         perc_RRHomAAHom = RRHomAAHom/(RRHomRRHom+RRHomRAHet+RRHomAAHom),
         perc_RAHetRRHom = RAHetRRHom/(RAHetRRHom+RAHetRAHet+RAHetAAHom),
         perc_RAHetRAHet = RAHetRAHet/(RAHetRRHom+RAHetRAHet+RAHetAAHom),
         perc_RAHetAAHom = RAHetAAHom/(RAHetRRHom+RAHetRAHet+RAHetAAHom),
         perc_AAHomRRHom = AAHomRRHom/(AAHomRRHom+AAHomRAHet+AAHomAAHom),
         perc_AAHomRAHet = AAHomRAHet/(AAHomRRHom+AAHomRAHet+AAHomAAHom),
         perc_AAHomAAHom = AAHomAAHom/(AAHomRRHom+AAHomRAHet+AAHomAAHom), 
         perc_sum = (perc_RRHomRRHom+perc_RRHomRAHet+perc_RRHomAAHom+
                       perc_RAHetRRHom+perc_RAHetRAHet+perc_RAHetAAHom+
                       perc_AAHomRRHom+perc_AAHomRAHet+perc_AAHomAAHom)) %>% 
  filter(status=="filtered") %>% 
  filter(eff_pop==1000) %>% 
  group_by(coverage) %>% 
  summarise(
    # RR 
    mean_perc_RRHomRRHom = mean(perc_RRHomRRHom), 
    mean_perc_RRHomRAHet = mean(perc_RRHomRAHet), 
    mean_perc_RRHomAAHom = mean(perc_RRHomAAHom), 
    sd_RRHomRRHom=sd(perc_RRHomRRHom),
    sd_RRHomRAHet = sd(perc_RRHomRAHet),
    sd_RRHomAAHom = sd(perc_RRHomAAHom),
    # RA 
    mean_perc_RAHetRRHom = mean(perc_RAHetRRHom), 
    mean_perc_RAHetRAHet = mean(perc_RAHetRAHet), 
    mean_perc_RAHetAAHom = mean(perc_RAHetAAHom), 
    sd_RAHetRRHom = sd(perc_RAHetRRHom),
    sd_RAHetRAHet = sd(perc_RAHetRAHet),
    sd_RAHetAAHom = sd(perc_RAHetAAHom),
    # AA 
    mean_perc_AAHomRRHom = mean(perc_AAHomRRHom), 
    mean_perc_AAHomRAHet = mean(perc_AAHomRAHet), 
    mean_perc_AAHomAAHom = mean(perc_AAHomAAHom),
    sd_AAHomRRHom = sd(perc_AAHomRRHom),
    sd_AAHomRAHet = sd(perc_AAHomRAHet),
    sd_AAHomAAHom = sd(perc_AAHomAAHom))
error_heat_r_means <- error_heat_r %>% 
  dplyr::select("coverage","mean_perc_RRHomRRHom", "mean_perc_RRHomRAHet", "mean_perc_RRHomAAHom", 
                "mean_perc_RAHetRRHom", "mean_perc_RAHetRAHet", "mean_perc_RAHetAAHom",
                "mean_perc_AAHomRRHom", "mean_perc_AAHomRAHet", "mean_perc_AAHomAAHom") %>% 
  pivot_longer(cols = c("mean_perc_RRHomRRHom", "mean_perc_RRHomRAHet", "mean_perc_RRHomAAHom", 
                        "mean_perc_RAHetRRHom", "mean_perc_RAHetRAHet", "mean_perc_RAHetAAHom",
                        "mean_perc_AAHomRRHom", "mean_perc_AAHomRAHet", "mean_perc_AAHomAAHom"), 
               names_to="error_type", values_to="mean") %>% 
  mutate(true_GT = ifelse(error_type %in% c("mean_perc_RRHomRRHom", "mean_perc_RRHomRAHet", "mean_perc_RRHomAAHom"), "RR", 
                          ifelse(error_type %in% c("mean_perc_RAHetRRHom", "mean_perc_RAHetRAHet", "mean_perc_RAHetAAHom"), "RA", "AA"))) %>% 
  mutate(imputed_GT = ifelse(error_type %in% c("mean_perc_RRHomRRHom", "mean_perc_RAHetRRHom", "mean_perc_AAHomRRHom"), "RR", 
                             ifelse(error_type %in% c("mean_perc_RRHomRAHet", "mean_perc_RAHetRAHet", "mean_perc_AAHomRAHet"), "RA", "AA"))) 

error_heat_r_sds <- error_heat_r %>% 
  dplyr::select("coverage","sd_RRHomRRHom", "sd_RRHomRAHet", "sd_RRHomAAHom", 
                "sd_RAHetRRHom", "sd_RAHetRAHet", "sd_RAHetAAHom",
                "sd_AAHomRRHom", "sd_AAHomRAHet", "sd_AAHomAAHom") %>% 
  pivot_longer(cols = c("sd_RRHomRRHom", "sd_RRHomRAHet", "sd_RRHomAAHom", 
                        "sd_RAHetRRHom", "sd_RAHetRAHet", "sd_RAHetAAHom",
                        "sd_AAHomRRHom", "sd_AAHomRAHet", "sd_AAHomAAHom"), 
               names_to="error_type", values_to="sd") %>% 
  mutate(true_GT = ifelse(error_type %in% c("sd_RRHomRRHom", "sd_RRHomRAHet", "sd_RRHomAAHom"), "RR", 
                          ifelse(error_type %in% c("sd_RAHetRRHom", "sd_RAHetRAHet", "sd_RAHetAAHom"), "RA", "AA"))) %>% 
  mutate(imputed_GT = ifelse(error_type %in% c("sd_RRHomRRHom", "sd_RAHetRRHom", "sd_AAHomRRHom"), "RR", 
                             ifelse(error_type %in% c("sd_RRHomRAHet", "sd_RAHetRAHet", "sd_AAHomRAHet"), "RA", "AA"))) 

error_heat_r_comb <- merge(error_heat_r_means,error_heat_r_sds,by=c("coverage","true_GT", "imputed_GT"))

error_heat_r_comb %>% 
  mutate(text = paste0(round((mean)*100,1), "%", " (+/- ", round((sd)*100,2),")")) %>% 
  filter(coverage == "0.5x") %>% 
  ggplot(aes(x = true_GT, y = imputed_GT, fill = 100*(mean))) + 
  geom_tile() + 
  geom_text(aes(label = text), 
            fontface="bold", color="black", size=5)+ 
  scale_fill_gradient(low = "aliceblue",high="cornflowerblue") + 
  labs(fill="mean percent\nof sites", x = "genotype in truth dataset", y="genotype in imputed data") + 
  theme_classic() + 
  scale_x_discrete(limits=rev) + 
  coord_equal() + 
  theme(strip.background = element_rect(color = "white"),
        strip.placement = "inside", 
        strip.text = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank())
#
################ Genetic analyses ################
### Relatedness of close relatives ----------------------
# read in relatedness files 
po <- read.delim("/parent_offspring_pairs", sep = "", header=F)
colnames(po)<-c("offspring","parent")
half_sibs <- read.delim("/half_sib_pairs", sep = "", header=F)
colnames(half_sibs)<-c("sib1","sib2")
# read in VCFtools relatedness estimates
cayo_relate <- do.call(rbind, lapply(list.files(path = "/plink_analyses/", pattern="*.relatedness2"), 
                                     FUN = function(i){
                                       tmp <- read.delim(paste0("/plink_analyses/",i))
                                       tmp$coverage <- gsub(pattern = ".relatedness2", replacement = "", x = tmp$coverage)
                                       return(tmp)
                                     }))
cayo_relate[cayo_relate$coverage=="highcov",]$coverage<-"reference"
cayo_relate$coverage <- factor(cayo_relate$coverage, levels = c("0.1x","0.5x","1x","3x","10x","reference"))
# make new column of relationship
cayo_relate$relationship <- "unrelated"
cayo_relate[cayo_relate$INDV1 == cayo_relate$INDV2,]$relationship <- "self"
# for each in offspring/parent pair, change relationship
for (i in 1:nrow(po)) {
  cayo_relate[which(cayo_relate$INDV1 == po$parent[i] & 
                      cayo_relate$INDV2 == po$offspring[i] | 
                      cayo_relate$INDV1 == po$offspring[i] & 
                      cayo_relate$INDV2 == po$parent[i]),]$relationship <- "parent/offspring"
}
# for each in half pair, change relationship
for (i in 1:nrow(half_sibs)) {
  cayo_relate[which(cayo_relate$INDV1 == half_sibs$sib1[i] & 
                      cayo_relate$INDV2 == half_sibs$sib2[i] | 
                      cayo_relate$INDV1 == half_sibs$sib1[i] & 
                      cayo_relate$INDV2 == half_sibs$sib2[i]),]$relationship <- "half-sibling"
}
cayo_relate$relationship <- factor(cayo_relate$relationship, 
                                   levels=c("self","unrelated","half-sibling","parent/offspring"))

# Plot 
cayo_relate %>% 
  filter(!relationship=="self") %>% 
  filter(coverage %in% c("0.1x","0.5x","10x","reference")) %>% 
  ggplot(aes(x = relationship, y = RELATEDNESS_PHI, color=coverage)) + 
  geom_violin() + 
  geom_point(aes(group=coverage), shape=1, size = 1,
             position=position_jitterdodge(jitter.width = 0.5,dodge.width = 0.9)) +
  scale_y_continuous(limits=c(-0.05,0.3)) + 
  theme_classic() + 
  labs(x="relationship", y="relatedness Phi") + 
  scale_color_manual(values = c(park_palette("Badlands")[c(1)],
                                park_palette("Saguaro")[c(1)], 
                                park_palette("Everglades")[c(2)], 
                                park_palette("Badlands")[c(3)]))
#
## Supplemental analysis (SA) - downsample to n sites
cayo_relate_RAD <- read.delim("/scratch/mwatowic/impute/plink_analyses/inSilico_RADseq/cayo_highcov_RAD.relatedness2")
cayo_relate_RAD$coverage <- c("in silico RADseq")
# make new column of relationship
cayo_relate_RAD$relationship <- "unrelated"
cayo_relate_RAD[cayo_relate_RAD$INDV1 == cayo_relate_RAD$INDV2,]$relationship <- "self"
# for each in offspring/parent pair, change relationship
for (i in 1:nrow(po)) {
  cayo_relate_RAD[which(cayo_relate_RAD$INDV1 == po$parent[i] & 
                          cayo_relate_RAD$INDV2 == po$offspring[i] | 
                          cayo_relate_RAD$INDV1 == po$offspring[i] & 
                          cayo_relate_RAD$INDV2 == po$parent[i]),]$relationship <- "parent/offspring"
}
# for each in half pair, change relationship
for (i in 1:nrow(half_sibs)) {
  cayo_relate_RAD[which(cayo_relate_RAD$INDV1 == half_sibs$sib1[i] & 
                          cayo_relate_RAD$INDV2 == half_sibs$sib2[i] | 
                          cayo_relate_RAD$INDV1 == half_sibs$sib1[i] & 
                          cayo_relate_RAD$INDV2 == half_sibs$sib2[i]),]$relationship <- "half-sibling"
}
cayo_relate_RAD$relationship <- factor(cayo_relate_RAD$relationship, 
                                       levels=c("self","unrelated","half-sibling","parent/offspring"))

# Plot Supp Analyses
rhes_relat_vs_RAD <- rbind(cayo_relate, 
                           cayo_relate_RAD) %>% 
  filter(coverage %in% c("0.1x","0.5x","10x","reference","in silico RADseq")) %>% 
  mutate(coverage = factor(coverage, levels = c("0.1x","0.5x","10x","reference","in silico RADseq"))) %>% 
  filter(!relationship=="self") %>% 
  ggplot(aes(x = relationship, y = RELATEDNESS_PHI, color=coverage)) + 
  geom_violin() + 
  geom_point(aes(group=coverage), position=position_jitterdodge(jitter.width = 0.5,dodge.width = 0.9), shape=1) +
  scale_y_continuous(limits=c(-0.07,0.3)) + 
  theme_classic() + 
  labs(x="relationship", y="relatedness Phi") + 
  scale_color_manual(values = c(park_palette("Badlands")[c(1)],
                                park_palette("Saguaro")[c(1)], 
                                park_palette("Everglades")[c(2)], 
                                park_palette("Badlands")[c(3)],
                                park_palette("Badlands")[c(2)])) + 
  theme(strip.text = element_text(face="bold",size=16),
        axis.text.x = element_text(size=14))
# ggsave(rhes_relat_vs_RAD,
#        filename = "/scratch/mwatowic/impute/plots/rhes_relat_vs_RAD.pdf",
#        device = "pdf", width = 12, height = 10)
#
### PCA ------------------------
pca_gelada <- do.call(rbind, 
                      lapply(list.files(path = "/ped/results/", 
                                        pattern="*_pca.eigenvec"), 
                             FUN = function(i){
                               tmp <- read.table(paste0("/ped/results/",
                                                        i),header = F)[,1:12]
                               tmp$coverage <- gsub(pattern = "gelada_", replacement = "", x = i)
                               tmp$coverage <- gsub(pattern = "_pca.eigenvec", replacement = "", x = tmp$coverage)
                               return(tmp)
                             }))
pca_gelada[pca_gelada$coverage=="geladaRef",]$coverage <- "reference"
pca_gelada$coverage <- factor(pca_gelada$coverage, levels = c("reference","10x","3x","1x","0.5x","0.1x"))

pca_gelada %>% 
  filter(coverage %in% c("0.1x", "3x", "reference")) %>% 
  ggplot(aes(x = V3, y = V4, color = population)) + 
  geom_point(size = 2, shape = 1) + 
  theme_classic() + 
  facet_wrap(.~coverage, ncol = 1) + 
  scale_color_manual(values = park_palette("Everglades")[c(3,4,2)]) + #c(3,1,2,4)
  theme(strip.text = element_text(size = 16, face="bold")) + 
  labs(x="PC1", y="PC2") + 
  scale_x_continuous(limits=c(-0.3,0.3)) + 
  scale_y_continuous(limits=c(-0.47,0.4))
#
