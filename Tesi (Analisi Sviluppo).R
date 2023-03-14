# LIBRERIE ---------------------------------------------------------------------
library(survival)
library(survminer)
library(reshape)
library(reshape2)
library(ggplot2)
library(readxl)
library(MASS)
library(GGally)

# ATTENZIONE: NELLA PARTE DI ANALISI NON PARAMETRICA, AGE E NODES SONO STATE
# CODIFICATE PERMANENTEMENTE IN CLASSI DA QUEL PUNTO IN POI. PER TORNARE 
# AD ANALIZZARE LE DUE VARIABILI DI PARTENZA IN FORMA DISCRETA, E' NECESSARIO
# ELIMINARE QUELLE DUE VARIABILI E RILANCIARLE NELLA PRIMA PARTE DI COMANDI.

# UPLOAD E DATA CLEANING -------------------------------------------------------
# Il dataset è 182x21.
dati_D = read_excel("carvalho-prognostic-biomarkers-NSCLC.xlsx", 1)
View(dati_D)
dim(dati_D)


# Rename delle variabili con spazio (dataset di sviluppo)
colnames(dati_D)[colnames(dati_D) == "Lymph nodes"] = "Nodes"
colnames(dati_D)[colnames(dati_D) == "RT Protocol"] = "RTProtocol"
colnames(dati_D)[colnames(dati_D) == "Total dose (1st)"] = "TotalDose1st"
colnames(dati_D)[colnames(dati_D) == "Total Dose (2nd)"] = "TotalDose2nd"
colnames(dati_D)[colnames(dati_D) == "IL 6"] = "IL6"
colnames(dati_D)[colnames(dati_D) == "IL 8"] = "IL8"
colnames(dati_D)[colnames(dati_D) == "Cyfra 21-1"] = "Cyfra211"
colnames(dati_D)[colnames(dati_D) == "WHO-PS"] = "WHOPS"
colnames(dati_D)[colnames(dati_D) == "CA-9"] = "CA9"
colnames(dati_D)[colnames(dati_D) == "FEV1s%"] = "FEV1s"
names(dati_D)
View(dati_D)
dim(dati_D)
str(dati_D)


# Fattorizzazioni (dataset di sviluppo)
dati_D$histology = as.factor(dati_D$histology)
dati_D$stage = as.factor(dati_D$stage)
dati_D$Gender = as.factor(dati_D$Gender)
dati_D$RTProtocol = as.factor(dati_D$RTProtocol)
str(dati_D)


# Conversione dei valori "alive" e "dead" di dati_D di Status in 0 e 1
convert_status <- function(dati_D) {
  dati_D$Status <- ifelse(dati_D$Status == "alive", 0, 1)
  return(dati_D)
}
dati_D <- convert_status(dati_D)
View(dati_D)
str(dati_D)


# ANALISI ESPLORATIVA ----------------------------------------------------------
# Fattori
# Table per status (dataset di sviluppo)
table(dati_D$Status)

# Barplot per Gender (dataset di sviluppo)
ggplot(dati_D, aes(x = Gender)) + geom_bar(fill = c("steelblue", "pink")) + geom_text(stat = 'count', aes(label = ..count..), vjust = 2) +
  labs(title = "Distribuzione per genere", x = "Genere", y = "Frequenza")

# Barplot per histology (dataset di sviluppo)
ggplot(dati_D, aes(x = histology)) + geom_bar() + geom_text(stat = 'count', aes(label = ..count..), vjust = 2) +
  labs(title = "Distribuzione per istotipo", x = "Histology", y = "Frequenza")

# Barplot per stage (dataset di sviluppo)
ggplot(dati_D, aes(x = stage)) + geom_bar() + geom_text(stat = 'count', aes(label = ..count..), vjust = 2) +
  labs(title = "Distribuzione per stage tumorale", x = "Stage", y = "Frequenza")

# Barplot per RTProtocol (dataset di sviluppo)
ggplot(dati_D, aes(x = RTProtocol)) + geom_bar() + geom_text(stat = 'count', aes(label = ..count..), vjust = 2) +
  labs(title = "Distribuzione per protocollo RT", x = "Protocollo", y = "Frequenza")

# Barplot per classi d'età (dataset di sviluppo)
# NOTA: D'ora in poi questa variabile verrà considerata solo come suddivisa in classi
dati_D$age <- cut(dati_D$age, breaks = c(40, 50, 60, 70, 80, 90), labels = c("40-50", "50-60", "60-70", "70-80", "80-90"))
dati_D$age <- as.factor(dati_D$age)
table(dati_D$age)

ggplot(dati_D, aes(x = age)) + 
  geom_bar(stat = "count", fill = "steelblue") +
  labs(title = "Distribuzione per classid'età", x = "Età", y = "Frequenza")

# Barplot per Linfonodi (accorpati in due classi: 0-1 e 2-3-4) (dataset di sviluppo)
# NOTA: D'ora in poi questa variabile verrà considerata solo come suddivisa in classi
bins = c(0, 1, 4)
labels_lymph <- c("0-1", "2-3-4")
dati_D$Nodes <- cut(dati_D$Nodes, bins, labels = labels_lymph)

ggplot(dati_D, aes(x = Nodes)) + 
  geom_bar(stat = "count") +
  labs(title = "Distribuzione per linfonodi", x = "Linfonodi", y = "Frequenza")

# Variabili quantitative, dataset di sviluppo (uni- e bi-variate mediante GGpairs)
dati_D_quantitative = dati_D[,c(8,9,10,12,13,14,15,16,17,18,19,20,21)]
ggpairs(dati_D_quantitative)


# ANALISI NON PARAMETRICA ------------------------------------------------------
# Curva di sopravvivenza generale (Sviluppo)
sopr_D = survfit(Surv(Survival, Status) ~ 1, data = dati_D)
ggsurvplot(sopr_D, data = dati_D, conf.int = TRUE, conf.int.style = "step")

# Curva di sopravvivenza per Sesso (Sviluppo)
sopr_gender_D = survfit(Surv(Survival, Status) ~ Gender, data = dati_D)
ggsurvplot(sopr_gender_D, data = dati_D, conf.int = TRUE, pval = TRUE, pval.method = TRUE, conf.int.style = "step")

# Curva di sopravvivenza per Terapia (RTProtocol) (Sviluppo)
sopr_rtprotocol_D = survfit(Surv(Survival, Status) ~ RTProtocol, data = dati_D)
ggsurvplot(sopr_rtprotocol_D, data = dati_D, conf.int = TRUE, pval = TRUE, pval.method = TRUE, conf.int.style = "step")

# Curva di sopravvivenza per Istologico (Sviluppo)
sopr_hist_D = survfit(Surv(Survival, Status) ~ histology, data = dati_D)
ggsurvplot(sopr_hist_D, data = dati_D, conf.int = TRUE, pval = TRUE, pval.method = TRUE)

# Curva di sopravvivenza per Stage (Sviluppo)
sopr_stage_D = survfit(Surv(Survival, Status) ~ stage, data = dati_D)
ggsurvplot(sopr_stage_D, data = dati_D, conf.int = F, pval = TRUE, pval.method = TRUE)

# Curva di sopravvivenza per WHOPS (Sviluppo)
sopr_whops_D = survfit(Surv(Survival, Status) ~ WHOPS, data = dati_D)
ggsurvplot(sopr_whops_D, data = dati_D, conf.int = TRUE, pval = TRUE, pval.method = TRUE)

# Curva di sopravvivenza per Linfonodi suddiviso in classi (Sviluppo)
sopr_nodes_D_cl = survfit(Surv(Survival, Status) ~ Nodes, data = dati_D)
ggsurvplot(sopr_nodes_D_cl, data = dati_D, conf.int = TRUE, pval = TRUE, pval.method = TRUE, conf.int.style = "step")

# Curva di sopravvivenza per classi d'eta (Sviluppo)
sopr_age_D <- survfit(Surv(Survival, Status) ~ age, data = dati_D)
ggsurvplot(sopr_age_D, data = dati_D, conf.int = TRUE, pval = TRUE, pval.method = TRUE)


# CORRELAZIONI -----------------------------------------------------------------
# correlazione tra variabili fisiche e biomarker
daticor_D = na.omit(dati_D[,c(8,9,10,14,15,16,17,18,19,20,21)])
cormat <- round(cor(daticor_D),2)
melted_cormat <- melt(cormat)
ggplot(data = melted_cormat, aes(x=Var1, y=Var2, fill=value)) + geom_tile()


# MODELLI DI COX ---------------------------------------------------------------
# Modello di Cox standard
modcox_D_step = step(coxph(Surv(Survival, Status) ~ age + stage + histology + Gender + WHOPS + FEV1s + Nodes + RTProtocol + TotalDose1st + TotalDose2nd + GTV + OPN + CA9 + IL6 + IL8 + CRP + CEA + Cyfra211, data = na.omit(dati_D)))
modcox_D_final = coxph(Surv(Survival, Status) ~ CRP + CEA + IL6 + RTProtocol + WHOPS + histology + Gender + GTV, data = dati_D)
summary(modcox_D_final)
test.ph.D = cox.zph(modcox_D_final)
test.ph.D
ggcoxzph(test.ph.D)


# Modello di Cox con frailty() usata su RTProtocol
modcox_D_step_RTP_frailty = step(coxph(Surv(Survival, Status) ~ age + stage + histology + Gender + WHOPS + FEV1s + Nodes + frailty(RTProtocol) + TotalDose1st + TotalDose2nd + GTV + OPN + CA9 + IL6 + IL8 + CRP + CEA + Cyfra211, data = na.omit(dati_D)))


# Modello di Cox con cluster() usata su RTProtocol
modcox_D_step_RTP_cluster = step(coxph(Surv(Survival, Status) ~ age + stage + histology + Gender + WHOPS + FEV1s + Nodes + cluster(RTProtocol) + TotalDose1st + TotalDose2nd + GTV + OPN + CA9 + IL6 + IL8 + CRP + CEA + Cyfra211, data = na.omit(dati_D)))
modcox_D_final_RTP_cluster = coxph(Surv(Survival, Status) ~ cluster(RTProtocol) + TotalDose2nd + CRP + CEA + IL6 + histology + GTV + WHOPS + Gender , data = dati_D)
summary(modcox_D_final_RTP_cluster)
test.ph.D = cox.zph(modcox_D_final_RTP_cluster)
test.ph.D
ggcoxzph(test.ph.D)


m1 = coxph(Surv(Survival, Status) ~ age + histology*Gender + RTProtocol, data = na.omit(dati_D))
cox.zph(m1)
