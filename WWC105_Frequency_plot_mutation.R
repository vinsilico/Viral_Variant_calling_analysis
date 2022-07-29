library(data.table)
library(purrr)
library(dplyr)
library(ggplot2)
library(reshape2)
library(viridis)
library(stringr)
library(gtools)
library(stringr)
library(tidyr)
library(openxlsx)
library(wesanderson)

# freq_data <- fread("freq.txt")
# freq_data1 <- melt(freq_data)
# freq_data1
# freq_data1[, 'sample'] <- 1
# freq_data1
# write.xlsx(freq_data1, "freq_data1.xlsx")
# freq_data2 <- fread("freq_data2.txt")

freq_data <- read.table("All_VCF_depth_1.txt", sep='\t', strip.white = TRUE, fill = TRUE)
freq_data_group <- freq_data %>% mutate(group = cumsum(V1 =="Name"))
head(freq_data_group)

File_name_list <- freq_data_group %>% filter(V1=="Name")
File_name_list <- subset(File_name_list, select=-c(V1))


#names(File_name_list)[1] <- "File_Name"
#'freq_data_group[, 'File_Name'] <- NA
# temp2 <- freq_data_group %>% mutate(File_Name = if_else(group == File_name_list$group), File_name_list$File_Name)

freq_data_group_2 <- left_join(freq_data_group, File_name_list, by='group')

names(freq_data_group_2)[1] <- "POS"
names(freq_data_group_2)[2] <- "Depth"
names(freq_data_group_2)[4] <- "File_Name"

freq_data_group_3 <- freq_data_group_2 %>% filter(POS!="Name")
head(freq_data_group_3)

freq_data_group_3[, 'caller'] <- NA

temp2 <- freq_data_group_3 %>% mutate(caller = ifelse(str_detect(File_Name,"VS"), "VarScan", caller))
temp3 <- temp2 %>% mutate(caller = ifelse(str_detect(File_Name,"bcf"), "BCFtools", caller))
temp4 <- temp3 %>% mutate(caller = ifelse(str_detect(File_Name,"lfq"), "LoFreq", caller))
temp5 <- temp4 %>% mutate(caller = ifelse(str_detect(File_Name,"ivar"), "iVar", caller))
temp6 <- temp5 %>% mutate(caller = ifelse(str_detect(File_Name,"GATK"), "GATK", caller))
freq_data_group_4 <- temp6 %>% mutate(caller = ifelse(str_detect(File_Name,"fb"), "FreeBayes", caller))
#head(freq_data_group_4)
#write.csv(quasi_merged_caller, file="snpcall_benchmark_DELTA_merged_caller.txt")
#freq_data_group_4 %>% count(caller)
freq_data_group_4[, 'Replicate'] <- NA
temp7 <- freq_data_group_4 %>% mutate(Replicate = ifelse(str_detect(File_Name,"^01"), "Rep1", Replicate))
temp8 <- temp7 %>% mutate(Replicate = ifelse(str_detect(File_Name,"^02"), "Rep2", Replicate))
temp9 <- temp8 %>% mutate(Replicate = ifelse(str_detect(File_Name,"^03"), "Rep3", Replicate))
temp10 <- temp9 %>% mutate(Replicate = ifelse(str_detect(File_Name,"^04"), "Rep4", Replicate))
head(temp10)
temp10[, 'Rep_No'] <- NA
temp11 <- temp10 %>% mutate(Rep_No = ifelse(str_detect(File_Name,"^01"), "01", Rep_No))
temp12 <- temp11 %>% mutate(Rep_No = ifelse(str_detect(File_Name,"^02"), "02", Rep_No))
temp13 <- temp12 %>% mutate(Rep_No = ifelse(str_detect(File_Name,"^03"), "03", Rep_No))
temp14 <- temp13 %>% mutate(Rep_No = ifelse(str_detect(File_Name,"^04"), "04", Rep_No))
head(temp14)

temp14 %>% count(caller)
temp14[, 'Rep_Num'] <- NA
temp14[, 'Sample_No'] <- NA
temp15 <- temp14 %>% separate(File_Name, c("Rep_Num", "Sample_No"), extra='drop', remove=FALSE)
head(temp15)

temp16 <- temp15 %>% filter(Sample_No == "1" | Sample_No == "2")

alpha <- fread("alpha.txt", header=TRUE)
beta <- fread("beta.txt", header=TRUE)
delta <- fread("delta.txt", header=TRUE)

temp20 <- temp15 %>% filter(POS %in% alpha$POS)


temp21 <- temp15 %>% filter(POS %in% beta$POS)


temp22 <- temp15 %>% filter(POS %in% delta$POS)


jit <- position_jitter(seed=123)
point_plot <- ggplot(temp20, aes(x=factor(Sample_No, level =c(1,2,9,10,15,16,17,31,32,33,34,35,36,37,38,39,40,41,42)), y=factor(POS,level=c(241,913,3037,3267,5388,5986,6954,11288,14408,14676,15279,16176,21766,21994,23063,23271,23403,23604,23709,24506,24914,27972,28048,28111,28274,28280,28881,28977)),color=as.numeric(Depth))) + 
  geom_point(size=2) +
  #scale_color_viridis()+
  scale_colour_gradientn(colours=rainbow(7)) +
  #scale_color_gradient(low="red", high="blue") +
  #geom_jitter(aes(color=as.numeric(Depth)), size=1) +
  #scale_shape_manual(values=rep(15:22, len=6))  +
  #theme(legend.spacing.x = unit(0.2, 'cm'), text = element_text(size=8)) +
  #theme_bw(base_size=15) + theme(legend.position = "none") + 
  labs(x='Mix Samples', y='Alpha SNPs/Indels position')
point_plot + facet_wrap(~caller)

ggsave("Mix_Alpha_freq_plot-rainbow7-jitter_facet.png",units="in", width= 12, height = 6, device='png', dpi=300)


point_plot <- ggplot(temp21, aes(x=factor(Sample_No, level =c(1,2,9,10,15,16,17,31,32,33,34,35,36,37,38,39,40,41,42)), y=factor(POS,level=c(174,241,913,1059,2692,3037,5230,10323,11288,14408,21801,22206,22287,22813,23012,23063,23403,23664,25563,25904,26456,28253,28254,28887,29557)),color=as.numeric(Depth))) + 
  geom_point(size=2) +
  #scale_color_viridis()+
  scale_colour_gradientn(colours=rainbow(7)) +
  #scale_color_gradient(low="light blue", high="red") +
  #geom_jitter(aes(color=as.numeric(Depth), shape=caller), size=2, position=jit) +
  #scale_shape_manual(values=rep(0:5, len=6))  +
  #theme(legend.spacing.x = unit(0.2, 'cm'), text = element_text(size=8)) +
  #theme_bw(base_size=15) + theme(legend.position = "none") + 
  labs(x='Mix Samples', y='Beta SNPs/Indels position')
point_plot + facet_wrap(~caller)
ggsave("Mix_Beta_freq_plot-rainbow7-jitter_facet.png",units="in", width= 12, height = 6, device='png', dpi=300)

point_plot <- ggplot(temp22, aes(x=factor(Sample_No, level =c(1,2,9,10,15,16,17,31,32,33,34,35,36,37,38,39,40,41,42)), y=factor(POS,level=c(210,241,1191,1267,3037,5184,6539,9891,11418,12946,14408,15451,15603,16466,18486,20262,20320,21618,22029,22917,22995,23403,23604,24410,25469,26767,26985,27638,27739,27752,27916,28248,28274,28881,29402,29427,29742)),color=as.numeric(Depth))) + 
  geom_point(size=2) +
  #scale_color_viridis()+
  scale_colour_gradientn(colours=rainbow(7)) +
  #scale_color_gradient(low="light blue", high="red") +
  #geom_jitter(aes(color=as.numeric(Depth), shape=caller), size=2, position=jit) +
  #scale_shape_manual(values=rep(0:5, len=6))  +
  #theme(legend.spacing.x = unit(0.2, 'cm'), text = element_text(size=8)) +
  #theme_bw(base_size=15) + theme(legend.position = "none") + 
  labs(x='Mix Samples', y='Delta SNPs/Indels position')
point_plot
ggsave("Mix_Delta_freq_plot-rainbow7-jitter_facet.png",units="in", width= 12, height = 6, device='png', dpi=300)
