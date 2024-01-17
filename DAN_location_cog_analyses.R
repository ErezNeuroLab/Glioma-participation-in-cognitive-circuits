# DAN analyses
library(ggplot2)
library(reshape)
library(lm.beta)

DAN_cog_data <- read.csv("~/Documents/Documents - Ayan’s MacBook Pro/Suckling_Lab/ecog/DAN_cog_data.txt")
View(DAN_cog_data)

# stats analyses

summary(lm(data=DAN_cog_data, DAN_connectivity ~ Location_simp + Hemisphere + overlap + Antiseizure + Age + Gender))
summary(lm.beta(lm(data=DAN_cog_data, DAN_connectivity ~ Location_simp + Hemisphere + overlap + Antiseizure + Age + Gender)))
ummary(lm(data=DAN_cog_data, latest_value_heart ~ preop_heart + DAN_connectivity + Hemisphere + Location_simp + latest_assessment_heart+Age + Gender))

summary(lm.beta(lm(data=DAN_cog_data, latest_value_heart ~ preop_heart + DAN_connectivity + Hemisphere + Location_simp + latest_assessment_heart+Age + Gender)))

summary(lm(data=DAN_cog_data, latest_value_salt_lighthouse ~ preop_salt_lighthouse + DAN_connectivity + Hemisphere + Location_simp + latest_assessment_salt_lighthouse+Age + Gender))
summary(lm.beta(lm(data=DAN_cog_data, latest_value_salt_lighthouse ~ preop_salt_lighthouse + DAN_connectivity + Hemisphere + Location_simp + latest_assessment_salt_lighthouse+Age + Gender)))

# visualizing

ggplot(data=DAN_cog_data, aes(x=Location_simp, y=DAN_connectivity)) + geom_boxplot()+ theme_classic()+xlab('Tumour Location') + ylab('Tumour-DAN Connectivity') + geom_point()

DAN_cog_data_melt <- melt(DAN_cog_data,measure.vars = c("preop_heart", "latest_value_heart", "delta_heart"))
ggplot(data=DAN_cog_data_melt, aes(x=DAN_connectivity, y=value)) + facet_grid(cols = vars(variable)) + geom_jitter(aes(color = Location_simp), size=2) + theme_bw() + xlab("Tumour-DAN Connectivity") + ylab("Overall Accuracy on Heart Cancellation Task") + geom_smooth(method = "glm", color = "black")

summary(lm(data=DAN_cog_data[DAN_cog_data$latest_assessment_heart!='Postop',], latest_value_heart ~ preop_heart + DAN_connectivity + Hemisphere + Location_simp + latest_assessment_heart + Age + Gender))
summary(lm.beta(lm(data=DAN_cog_data[DAN_cog_data$latest_assessment_heart!='Postop',], latest_value_heart ~ preop_heart + DAN_connectivity + Hemisphere + Location_simp + latest_assessment_heart + Age + Gender)))

#adding pathology subplots

ggplot(data=DAN_cog_data_melt, aes(x=DAN_connectivity, y=value)) + facet_grid(rows = vars(Pathology_simp), cols = vars(variable)) + geom_point(aes(color=Hemisphere, shape = Location_simp), size=2) + theme_bw() + xlab("Tumour-DAN Connectivity") + ylab("Overall Accuracy on Heart Cancellation Task") + geom_smooth(method = "lm", color = "black")


