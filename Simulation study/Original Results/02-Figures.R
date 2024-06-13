## Begin by setting the working drive to location to where the
## folders "Data Analysis", "Functions", and "Simulation Study" are.

delim <- "/"

require(caret)
library(dplyr)
library(ggplot2)

parlist <- as.matrix(data.frame(read.csv(paste0("Simulation study",delim,"Original Results",delim,"args_list.csv"),header = TRUE, fileEncoding="UTF-8-BOM")))

L <- nrow(parlist)
Sim_Result <- NULL
for(j in 1:L){
  args_list <- as.numeric(unlist(parlist[j,]))
  
  K <- as.numeric(args_list[1]) #number of false nulls
  eta <- as.numeric(args_list[2])    #theta value in paper draft
  C <- as.numeric(args_list[3]) #Controls the heterogeneity 0.5, 1.5, 3
  sig_sp <- as.numeric(args_list[4]) #Spatial clustering of signals
  start <- 0
  
  res <- read.csv(paste0("Simulation study",delim,"Original Results",delim,"Summarized",delim,"OriginalRes K=",K," eta=",eta," C=",C," sig_sp=",sig_sp," start=",start,".csv"))
  res_full <- read.csv(paste0("Simulation study",delim,"Original Results",delim,"By_iteration",delim,"OriginalRes_full K=",K," eta=",eta," C=",C," sig_sp=",sig_sp," start=",start,".csv"))
  
  disc_row <- seq(2,44,3)
  FDR_vals <- as.matrix(res_full[,disc_row+2]/res_full[,disc_row])
  FDR_vals[is.nan(FDR_vals)] <- 0
  
  summ_res <- data.frame(Method = res[,1], 
                     Total.Rej  = apply(res_full[,disc_row],2,mean),
                     Correct.Rej = apply(res_full[,disc_row+1],2,mean),
                     False.Rej = apply(res_full[,disc_row+2],2,mean), 
                 FDR = apply(FDR_vals,2,mean),
                FDR.sd = apply(FDR_vals,2,sd), row.names = NULL)
  
  SetSim_Result <- cbind(rep(K,15),rep(eta,15),rep(C,15),rep(sig_sp,15),summ_res,rep(paste0("K",K,"eta",eta),15))
  
  SetSim_Result <- data.frame(SetSim_Result)
  
  Sim_Result <- rbind(Sim_Result,SetSim_Result)
  
  }

rownames(Sim_Result) <- NULL
colnames(Sim_Result) <- c("K","eta","C","sig_sp","method","TolRej","CorRej","FalRej","FDR", "FDR.sd","Group")

#change the eta to theta for data generating parameter not MMW eta
Sim_Result$sig_sp <- factor(Sim_Result$sig_sp, levels = c("0.01", "5", "10"))
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
                                               "ADAPT", "SWFDR", "CAMT","LAWS"), ]

cols4methods <- c("blue","blueviolet","brown","black","darkgrey","red","darkgreen","darkgoldenrod1","chartreuse","cadetblue1")





Sim_Res1 <- Sim_Res[Sim_Res$Method %in% c("ADAPT", "IHW",
                                          "SWFDR", "CAMT",
                                          "WABH 0.9 CAMT", 
                                          "Adaptive BH","LAWS"),]

Sim_Res1$Method[Sim_Res1$Method=="WABH 0.9 CAMT"] <- "WABH"
Sim_Res1$Method[Sim_Res1$Method=="ADAPT"] <- "AdaPT"



### Figure 4 in the main document

grDevices::cairo_pdf(paste0("Simulation study",delim,"Figures",delim,"Power_goodFDR-only.pdf"),width=10,height=7)
ggplot(Sim_Res1 %>% filter(Method %in% c("AdaPT", "Adaptive BH", "IHW", "WABH", "LAWS")), 
       aes(x=sig_sp, y=Power, color = Method, group = method)) +
  scale_color_manual(values=cols4methods[c(1,2,4,5,7)])+
  geom_line(position=position_dodge(0.05), linewidth = 0.3) +
  geom_point(shape=19, size=2, position=position_dodge(width=0.05)) +
  xlab("s") + 
  ylab("Power") + 
  coord_cartesian(ylim=c(0, NA), clip = "off") +
  facet_grid(eta~C2 + K, scales = "free_y", space = "free_x", switch = "y") + 
  theme_bw(base_size = 12, base_family = "Helvetica") + 
  theme(strip.background = element_blank(),
        strip.placement = "outside",
        axis.text=element_text(color="#666666",size=10), 
        axis.title = element_text(color="#666666", 
                                  face="bold", size=22), 
        strip.text.x = element_text(size=12), 
        strip.text.y = element_text(size=12))  #+ theme_bw() 
dev.off()

### Figure 6 in supporting material

grDevices::cairo_pdf(paste0("Simulation study",delim,"Figures",delim,"Power_all-methods.pdf"),width=10,height=7)
ggplot(Sim_Res1, aes(x=sig_sp, y=Power, color = Method, group = method)) +
  scale_color_manual(values=cols4methods)+
  geom_line(position=position_dodge(0.05), linewidth = 0.3) +
  geom_point(shape=19, size=2, position=position_dodge(width=0.05)) +
  xlab("s") + 
  ylab("Power") + 
  coord_cartesian(ylim=c(0, NA), clip = "off") +
  facet_grid(eta~C2 + K, scales = "free_y", space = "free_x", switch = "y") + 
  theme_bw(base_size = 12, base_family = "Helvetica") + 
  theme(strip.background = element_blank(),
        strip.placement = "outside",
        axis.text=element_text(color="#666666",size=10), 
        axis.title = element_text(color="#666666", 
                                  face="bold", size=22), 
        strip.text.x = element_text(size=12), 
        strip.text.y = element_text(size=12))  #+ theme_bw() 
dev.off()


### Figure 3 in the main document

grDevices::cairo_pdf(paste0("Simulation study",delim,"Figures",delim,"FDR_all-methods.pdf"),width=10,height=7)
ggplot(Sim_Res1, aes(x=sig_sp, y=FDR, color = Method, group = method)) +
  scale_color_manual(values=cols4methods)+
  geom_hline(yintercept=0.05, linetype = "dashed") +
  geom_line(position=position_dodge(0.05), linewidth = 0.3) +
  geom_point(shape=19, size=2, position=position_dodge(width=0.05)) +
  xlab("s") + 
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


### Figure 7 in supporting material

grDevices::cairo_pdf(paste0("Simulation study",delim,"Figures",delim,"FDR-w-SE_all-methods.pdf"),width=10,height=7)
pos_wid <- 0.5
Sim_Res1 <- Sim_Res1 %>% 
  mutate(FDR.l  = FDR - FDR.sd, FDR.u  = FDR + FDR.sd) %>% 
  mutate(FDR.l = case_when(
    FDR.l < 0 ~ 0,
    TRUE ~ FDR.l
  ))
ggplot(Sim_Res1 , aes(x=sig_sp, y=FDR, color = Method, group = method)) +
  scale_color_manual(values=cols4methods)+
  geom_hline(yintercept=0.05, linetype = "dashed") +
  geom_point(shape=19, size=2, position=position_dodge(width=pos_wid)) +
  geom_errorbar(aes(ymin = FDR.l, ymax = FDR.u), 
                width = pos_wid, 
                position=position_dodge(width=pos_wid)) +  
  xlab("s") + 
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




cols4methods <- c("blue","blueviolet","brown","black","darkgrey","red","darkgreen","darkgoldenrod1","chartreuse")

Sim_Res2 <- Sim_Result[Sim_Result$method %in% c("Ten Rule","WABH 0.9 ADAPT",
                                                "WABH 0.9 CAMT", 
                                                "Adaptive BH", 
                                                "WABH Constant", 
                                                "Regular BH", 
                                                "WABH 0.5 CAMT"), ]

Sim_Res2$Method[Sim_Res2$Method=="WABH 0.9 ADAPT"] <- "WABH 0.9 AdaPT"

### Figure 5 in supporting material

grDevices::cairo_pdf(paste0("Simulation study",delim,"Figures",delim,"SupportInfoSimulationPlot5.pdf"),width=10,height=7)
ggplot(Sim_Res2, aes(x=sig_sp, y=Power, color = Method, group = method)) +
  scale_color_manual(values=cols4methods)+
  geom_line(position=position_dodge(0.05), linewidth = 0.3) +
  geom_point(shape=19, size=2, position=position_dodge(width=0.05)) +
  xlab("s") + 
  ylab("Power") + 
  coord_cartesian(ylim=c(0, NA), clip = "off") +
  facet_grid(eta~C2 + K, scales = "free_y", space = "free_x", switch = "y") + 
  theme_bw(base_size = 12, base_family = "Helvetica") + 
  theme(strip.background = element_blank(),
        strip.placement = "outside",
        axis.text=element_text(color="#666666",size=10), 
        axis.title = element_text(color="#666666", 
                                  face="bold", size=22), 
        strip.text.x = element_text(size=12), 
        strip.text.y = element_text(size=12))  #+ theme_bw() 
dev.off()


### Figure 4 in supporting material

grDevices::cairo_pdf(paste0("Simulation study",delim,"Figures",delim,"SupportInfoSimulationPlot4.pdf"),width=10,height=7)
ggplot(Sim_Res2, aes(x=sig_sp, y=FDR, color = Method, group = method)) +
  scale_color_manual(values=cols4methods)+
  geom_hline(yintercept=0.05, linetype = "dashed") +
  geom_line(position=position_dodge(0.05), linewidth = 0.3) +
  geom_point(shape=19, size=2, position=position_dodge(width=0.05)) +
  xlab("s") + 
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
