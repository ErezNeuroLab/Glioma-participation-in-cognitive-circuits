library(lme4)
library(car)
library(ggplot2)
library(reshape)
source("~/Documents/Documents - Ayan’s MacBook Pro/Suckling_Lab/ecog/get_shape.R")


PSC_SCA_tbl <- read.csv("~/Documents/Documents - Ayan’s MacBook Pro/Suckling_Lab/ecog/PSC_SCA_tbl.txt")
View(PSC_SCA_tbl)

PSC_SCA_tumour_tbl = PSC_SCA_tbl[PSC_SCA_tbl$tumour==1,]

ggplot(PSC_SCA_tumour_tbl,aes(x=DAN_emp,y=cont_alt)) + geom_point(aes(colour=patient_id))


lmm <- lmer(DAN_emp ~ cont_alt + (1 | patient_id), data=PSC_SCA_tumour_tbl) # lmm model
Anova(lmm)



# now it's time to plot

PSC_SCA_tumour_melt<- melt(PSC_SCA_tumour_tbl,measure.vars = c("VN_emp", "SMN_emp","DAN_emp", "VAN_emp","LIM_emp", "FPN_emp","DMN_emp"))
PSC_SCA_tumour_melt<-melt(PSC_SCA_tumour_melt, measure.vars = c("cont_countF","cont_alt"))

colnames(PSC_SCA_tumour_melt)[16] <- "PSC"
colnames(PSC_SCA_tumour_melt)[15] <- "contrast"
colnames(PSC_SCA_tumour_melt)[14] <- "emp_val"
colnames(PSC_SCA_tumour_melt)[13] <- "network"

point_shape = c()

for (val in 1:length(rownames(PSC_SCA_tumour_melt)))
{
  
  if (PSC_SCA_tumour_melt$contrast[val]=="cont_countF") {
    if (PSC_SCA_tumour_melt$network[val]=="VN_emp") {
      p_n <- PSC_SCA_tumour_melt$VN_p[val]
      p_e <- PSC_SCA_tumour_melt$p_countF[val]
      point_shape[val] <- get_shape(p_n, p_e)
    }
    else if (PSC_SCA_tumour_melt$network[val]=="SMN_emp") {
      p_n <- PSC_SCA_tumour_melt$SMN_p[val]
      p_e <- PSC_SCA_tumour_melt$p_countF[val]
      point_shape[val] <- get_shape(p_n, p_e)
    }
    else if (PSC_SCA_tumour_melt$network[val]=="DAN_emp") {
      p_n <- PSC_SCA_tumour_melt$DAN_p[val]
      p_e <- PSC_SCA_tumour_melt$p_countF[val]
      point_shape[val] <- get_shape(p_n, p_e)
    }
    else if (PSC_SCA_tumour_melt$network[val]=="VAN_emp") {
      p_n <- PSC_SCA_tumour_melt$VAN_p[val]
      p_e <- PSC_SCA_tumour_melt$p_countF[val]
      point_shape[val] <- get_shape(p_n, p_e)
    }
    else if (PSC_SCA_tumour_melt$network[val]=="LIM_emp") {
      p_n <- PSC_SCA_tumour_melt$LIM_p[val]
      p_e <- PSC_SCA_tumour_melt$p_countF[val]
      point_shape[val] <- get_shape(p_n, p_e)
    }
    else if (PSC_SCA_tumour_melt$network[val]=="FPN_emp") {
      p_n <- PSC_SCA_tumour_melt$FPN_p[val]
      p_e <- PSC_SCA_tumour_melt$p_countF[val]
      point_shape[val] <- get_shape(p_n, p_e)
    }
    else if (PSC_SCA_tumour_melt$network[val]=="DMN_emp") {
      p_n <- PSC_SCA_tumour_melt$DMN_p[val]
      p_e <- PSC_SCA_tumour_melt$p_countF[val]
      point_shape[val] <- get_shape(p_n, p_e)
    }
  } 
  else if (PSC_SCA_tumour_melt$contrast[val] == "cont_alt"){
    if (PSC_SCA_tumour_melt$network[val]=="VN_emp") {
      p_n <- PSC_SCA_tumour_melt$VN_p[val]
      p_e <- PSC_SCA_tumour_melt$p_alt[val]
      point_shape[val] <- get_shape(p_n, p_e)
    }
    else if (PSC_SCA_tumour_melt$network[val]=="SMN_emp") {
      p_n <- PSC_SCA_tumour_melt$SMN_p[val]
      p_e <- PSC_SCA_tumour_melt$p_alt[val]
      point_shape[val] <- get_shape(p_n, p_e)
    }
    else if (PSC_SCA_tumour_melt$network[val]=="DAN_emp") {
      p_n <- PSC_SCA_tumour_melt$DAN_p[val]
      p_e <- PSC_SCA_tumour_melt$p_alt[val]
      point_shape[val] <- get_shape(p_n, p_e)
    }
    else if (PSC_SCA_tumour_melt$network[val]=="VAN_emp") {
      p_n <- PSC_SCA_tumour_melt$VAN_p[val]
      p_e <- PSC_SCA_tumour_melt$p_alt[val]
      point_shape[val] <- get_shape(p_n, p_e)
    }
    else if (PSC_SCA_tumour_melt$network[val]=="LIM_emp") {
      p_n <- PSC_SCA_tumour_melt$LIM_p[val]
      p_e <- PSC_SCA_tumour_melt$p_alt[val]
      point_shape[val] <- get_shape(p_n, p_e)
    }
    else if (PSC_SCA_tumour_melt$network[val]=="FPN_emp") {
      p_n <- PSC_SCA_tumour_melt$FPN_p[val]
      p_e <- PSC_SCA_tumour_melt$p_alt[val]
      point_shape[val] <- get_shape(p_n, p_e)
    }
    else if (PSC_SCA_tumour_melt$network[val]=="DMN_emp") {
      p_n <- PSC_SCA_tumour_melt$DMN_p[val]
      p_e <- PSC_SCA_tumour_melt$p_alt[val]
      point_shape[val] <- get_shape(p_n, p_e)
    }
    
  }
}


PSC_SCA_tumour_melt$contrast <- factor(PSC_SCA_tumour_melt$contrast, levels = c("cont_alt", "cont_countF"))
levels(PSC_SCA_tumour_melt$contrast) <- c("Hard > Easy", "Easy > Rest")
levels(PSC_SCA_tumour_melt$network) <- c("VN", "SMN", "DAN", "VAN", "LIM", "FPN", "DMN")

PSC_SCA_tumour_melt$patient_id <- replace(PSC_SCA_tumour_melt$patient_id, PSC_SCA_tumour_melt$patient_id=="2017_06", "Patient 1")
PSC_SCA_tumour_melt$patient_id <- replace(PSC_SCA_tumour_melt$patient_id, PSC_SCA_tumour_melt$patient_id=="2018_04", "Patient 2")
PSC_SCA_tumour_melt$patient_id <- replace(PSC_SCA_tumour_melt$patient_id, PSC_SCA_tumour_melt$patient_id=="2018_05", "Patient 3")
PSC_SCA_tumour_melt$patient_id <- replace(PSC_SCA_tumour_melt$patient_id, PSC_SCA_tumour_melt$patient_id=="2019_01", "Patient 4")

levels(PSC_SCA_tumour_melt$patient_id) <- c("Patient 1", "Patient 2", "Patient 3", "Patient 4")
PSC_SCA_tumour_melt$point_shape <- factor(point_shape, levels = c("Both significant", "Significant PSC", "Significant network", "Neither significant"))
shapes <- c("\u25C6", "\u25B2", "\u25B6", "\u25CF")

shapes <- c(23,24,22,21)
ggplot(PSC_SCA_tumour_melt, aes(x=emp_val,y=PSC,fill=patient_id)) + geom_point(aes(shape = point_shape), size = 2) +facet_grid(rows = vars(contrast), cols = vars(network)) + scale_shape_manual(values=shapes) + xlab("Functional Network Correspondence") + ylab("% Signal Change") + theme_bw() + geom_line(color = "black", size = 0.75) + coord_flip()

