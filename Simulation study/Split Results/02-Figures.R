## Begin by setting the working drive to location to where the
## folders "Data Analysis", "Functions", and "Simulation Study" are.

delim <- "/"

require(caret)
library(dplyr)
library(ggplot2)

parlist <- as.matrix(data.frame(read.csv(paste0("Simulation study",delim,"Split Results",delim,"args_list.csv"),header = TRUE, fileEncoding="UTF-8-BOM")))

L <- nrow(parlist)
Sim_Result <- NULL
for(j in 1:L){
  args_list <- as.numeric(unlist(parlist[j,]))
  
  K <- as.numeric(args_list[1]) #number of false nulls
  eta <- as.numeric(args_list[2])    #theta value in paper draft
  C <- as.numeric(args_list[3]) #Controls the heterogeneity 0.5, 1.5, 3
  sig_sp <- as.numeric(args_list[4]) #Spatial clustering of signals
  start <- 0
  
  splitres_full <- read.csv(paste0("Simulation study",delim,"Split Results",delim,"By_iteration",delim,"SplitRes_full K=",K," eta=",eta," C=",C," sig_sp=",sig_sp," start=",start,".csv"))
  res_full <- read.csv(paste0("Simulation study",delim,"Original Results",delim,"By_iteration",delim,"OriginalRes_full K=",K," eta=",eta," C=",C," sig_sp=",sig_sp," start=",start,".csv"))
  
  res_full <- cbind(res_full[,which(colnames(res_full) == "WABH.MMW.0.9.camt.pm_Tot"):(which(colnames(res_full) == "WABH.MMW.0.9.camt.pm_Tot")+2)],
                    splitres_full[,-1])
  disc_row <- seq(1,10,3)
  FDR_vals <- as.matrix(res_full[,disc_row+2]/res_full[,disc_row])
  FDR_vals[is.nan(FDR_vals)] <- 0
  
  summ_res <- data.frame(Total.Rej  = apply(res_full[,disc_row],2,mean),
                         Correct.Rej = apply(res_full[,disc_row+1],2,mean),
                         False.Rej = apply(res_full[,disc_row+2],2,mean), 
                         FDR = apply(FDR_vals,2,mean),
                         FDR.sd = apply(FDR_vals,2,sd), row.names = NULL)
  
  
  SetSim_Result <- cbind(rep(K,4),rep(eta,4),rep(C,4),rep(sig_sp,4),c("Old Results","40/160 Split","100/100 Split","200/200 Split"),summ_res,rep(paste0("K",K,"eta",eta),4))
  
  SetSim_Result <- data.frame(SetSim_Result)
  
  Sim_Result <- rbind(Sim_Result,SetSim_Result)
  
}

rownames(Sim_Result) <- NULL
colnames(Sim_Result) <- c("K","eta","C","sig_sp","method","TolRej","CorRej","FalRej","FDR", "FDR.sd","Group")

#change the eta to theta for data generating parameter not MMW eta
Sim_Result$sig_sp <- factor(Sim_Result$sig_sp, levels = c("0.01", "5", "10"))
Sim_Result$method[Sim_Result$method == "Old Results"] <- "No Split"
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

Sim_Result$Power <- Sim_Result$CorRej/Sim_Result$K_num
Sim_Result$Method <- Sim_Result$method

cols4methods <- c("blue","black","red","darkgreen")

grDevices::cairo_pdf(paste0("Simulation study",delim,"Figures",delim,"SupportInfoSimulationPlot2.pdf"),width=10,height=7)
ggplot(Sim_Result, aes(x=sig_sp, y=Power, color = Method, group = method)) +
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


grDevices::cairo_pdf(paste0("Simulation study",delim,"Figures",delim,"SupportInfoSimulationPlot1.pdf"),width=10,height=7)
ggplot(Sim_Result, aes(x=sig_sp, y=FDR, color = Method, group = method)) +
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


grDevices::cairo_pdf(paste0("Simulation study",delim,"Figures",delim,"SupportInfoSimulationPlot1-SE.pdf"),width=10,height=7)
pos_wid <- 0.5
Sim_Res1 <- Sim_Result %>% 
  mutate(FDR.l  = FDR - FDR.sd, FDR.u  = FDR + FDR.sd) %>% 
  mutate(FDR.l = case_when(
    FDR.l < 0 ~ 0,
    TRUE ~ FDR.l
  ))
ggplot(Sim_Res1, aes(x=sig_sp, y=FDR, color = Method, group = method)) +
  scale_color_manual(values=cols4methods)+
  geom_hline(yintercept=0.05, linetype = "dashed") +
  geom_point(shape=19, size=2,
             position=position_dodge(width=pos_wid)) +
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

