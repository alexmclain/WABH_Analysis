## Begin by setting the working drive to location to where the
## folders "Data Analysis", "Functions", and "Simulation Study" are.

delim <- "/"

require(caret)
library(dplyr)
library(ggplot2)

parlist <- as.matrix(data.frame(read.csv(paste0("Simulation study",delim,"Misspecified Setting Results",delim,"args_list.csv"),header = TRUE, fileEncoding="UTF-8-BOM")))

L <- nrow(parlist)
Sim_Result <- NULL
for(j in 1:L){
  args_list <- as.numeric(unlist(parlist[j,]))
  
  K <- as.numeric(args_list[1]) #number of false nulls
  eta <- as.numeric(args_list[2])    #theta value in paper draft
  C <- as.numeric(args_list[3]) #Controls the heterogeneity 0.5, 1.5, 3
  mis_sp <- as.numeric(args_list[4])*4 #Spatial clustering of misspecified covariate
  start <- 0
  
  res <- read.csv(paste0("Simulation study",delim,"Misspecified Setting Results",delim,"Summarized",delim,"Res_mis K=",K," eta=",eta," C=",C," mis_sp=",mis_sp," start=",start,".csv"))
  
  SetSim_Result <- cbind(rep(K,14),rep(eta,14),rep(C,14),rep(mis_sp,14),res,rep(paste0("K",K,"eta",eta),14))
  
  SetSim_Result <- data.frame(SetSim_Result)
  
  Sim_Result <- rbind(Sim_Result,SetSim_Result)
  
}

rownames(Sim_Result) <- NULL
colnames(Sim_Result) <- c("K","eta","C","mis_sp","method","TolRej","CorRej","FalRej","FDR","Group")

#change the eta to theta for data generating parameter not MMW eta
Sim_Result$mis_sp <- factor(Sim_Result$mis_sp, levels = c("0.04", "20", "40"))
Sim_Result$Group2 <- Sim_Result$Group 
Sim_Result$Group2[Sim_Result$Group2=="K100eta0.25"] <- "K = 100, \u03B8 = 0.25"
Sim_Result$Group2[Sim_Result$Group2=="K100eta0.5"] <- "K = 100, \u03B8 = 0.5"
Sim_Result$Group2[Sim_Result$Group2=="K100eta0.75"] <- "K = 100, \u03B8 = 0.75"
Sim_Result$Group2[Sim_Result$Group2=="K500eta0.25"] <- "K = 500, \u03B8 = 0.25"
Sim_Result$Group2[Sim_Result$Group2=="K500eta0.5"] <- "K = 500, \u03B8 = 0.5"
Sim_Result$Group2[Sim_Result$Group2=="K500eta0.75"] <- "K = 500, \u03B8 = 0.75"
Sim_Result$Group2 <- factor(Sim_Result$Group2, levels=c("K = 100, \u03B8 = 0.25","K = 100, \u03B8 = 0.5","K = 100, \u03B8 = 0.75","K = 500, \u03B8 = 0.25","K = 500, \u03B8 = 0.5","K = 500, \u03B8 = 0.75"))


Sim_Result$K_num <- Sim_Result$C 
Sim_Result$K_num[Sim_Result$Group =="K100eta0.25" | 
                   Sim_Result$Group == "K100eta0.5" | 
                   Sim_Result$Group == "K100eta0.75"] <- 100
Sim_Result$K_num[Sim_Result$Group =="K500eta0.25" | 
                   Sim_Result$Group == "K500eta0.5" | 
                   Sim_Result$Group == "K500eta0.75"] <- 500

Sim_Result$eta_num <- Sim_Result$C 
Sim_Result$eta_num[Sim_Result$Group =="K100eta0.25" | 
                     Sim_Result$Group =="K500eta0.25"] <- 0.25
Sim_Result$eta_num[Sim_Result$Group =="K100eta0.5" | 
                     Sim_Result$Group =="K500eta0.5"] <- 0.5
Sim_Result$eta_num[Sim_Result$Group =="K100eta0.75" | 
                     Sim_Result$Group =="K500eta0.75"] <- 0.75

Sim_Result$K <- paste0("K = ",Sim_Result$K_num)
Sim_Result$K <- factor(Sim_Result$K, levels=c("K = 100","K = 500"))

Sim_Result$C2 <- factor(Sim_Result$C, levels = c(0.5,1.5,3), 
                        labels =c("C = 0.5", "C = 1.5",
                                  "C = 3"))
Sim_Result$eta <- paste0("\u03B8 = ",Sim_Result$eta_num)
Sim_Result$eta <- factor(Sim_Result$eta, levels=c("\u03B8 = 0.25",
                                                  "\u03B8 = 0.5","\u03B8 = 0.75"))

Sim_Result$Group3 <- paste0("C = ",Sim_Result$C,", K = ",Sim_Result$K_num)
Sim_Result$Group3 <- as.factor(Sim_Result$Group3)

Sim_Result$Group <- factor(Sim_Result$Group, levels=c("K100eta0.25","K100eta0.5","K100eta0.75","K500eta0.25","K500eta0.5","K500eta0.75"))
Sim_Result$method[Sim_Result$method=="WABH MMW constant pm"] <- "WABH Constant"
Sim_Result$method[Sim_Result$method=="WABH MMW 0.5 camt pm"] <- "WABH 0.5 CAMT"
Sim_Result$method[Sim_Result$method=="WABH MMW 0.9 camt pm"] <- "WABH 0.9 CAMT"
Sim_Result$method[Sim_Result$method=="WABH MMW 0.9 adapt pm"] <- "WABH 0.9 ADAPT"

Sim_Result$Power <- Sim_Result$CorRej/Sim_Result$K_num
Sim_Result$Method <- Sim_Result$method

Sim_Res <- Sim_Result[Sim_Result$method %in% c("WABH Constant","WABH 0.5 CAMT", 
                                               "WABH 0.9 CAMT", "Ten Rule", 
                                               "Adaptive BH", "IHW",
                                               "ADAPT", "SWFDR", "CAMT"), ]

cols4methods <- c("blue","blueviolet","brown","black","red","darkgreen","darkgoldenrod1","chartreuse","darkgrey")


### Figure in supporting material


Sim_Res1 <- Sim_Res[Sim_Res$Method %in% c("ADAPT", "IHW",
                                          "SWFDR", "CAMT",
                                          "WABH 0.9 CAMT", 
                                          "Adaptive BH"),]

Sim_Res1$Method[Sim_Res1$Method=="WABH 0.9 CAMT"] <- "WABH"
Sim_Res1$Method[Sim_Res1$Method=="ADAPT"] <- "AdaPT"

grDevices::cairo_pdf(paste0("Simulation study",delim,"Figures",delim,"SupportInfoSimulationPlot3.pdf"),width=10,height=7)
ggplot(Sim_Res1, aes(x=mis_sp, y=FDR, color = Method, group = method)) +
  scale_color_manual(values=cols4methods)+
  geom_hline(yintercept=0.05, linetype = "dashed") +
  geom_line(position=position_dodge(0.05), size = 0.3) +
  geom_point(shape=19, size=2, position=position_dodge(width=0.05)) +
  xlab("s_w") + 
  ylab("FDR") + 
  coord_cartesian(ylim=c(0, NA), clip = "off") +
  facet_grid(eta~C2 + K, scales = "fixed", space = "fixed", switch = "y") +
  theme_bw(base_size = 12, base_family = "Helvetica") + 
  theme(strip.background = element_blank(),
        strip.placement = "outside",
        axis.text=element_text(color="#666666",size=10), 
        axis.title = element_text(color="#666666", 
                                  face="bold", size=22), 
        strip.text.x = element_text(size=12), 
        strip.text.y = element_text(size=12)) 
dev.off()

