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
# Il dataset è 181x25. Comprende 4 variabili in più rispetto al primo, queste
# sono sempre associate a dei biomarker.
dati_V = read_excel("carvalho-prognostic-biomarkers-NSCLC.xlsx", 2)
View(dati_V)
dim(dati_V)


# Rename delle variabili con spazio (dataset di validazione)
colnames(dati_V)[colnames(dati_V) == "Lymph nodes"] = "Nodes"
colnames(dati_V)[colnames(dati_V) == "RT Protocol"] = "RTProtocol"
colnames(dati_V)[colnames(dati_V) == "Total dose (1st)"] = "TotalDose1st"
colnames(dati_V)[colnames(dati_V) == "Total Dose (2nd)"] = "TotalDose2nd"
colnames(dati_V)[colnames(dati_V) == "IL 6"] = "IL6"
colnames(dati_V)[colnames(dati_V) == "IL 8"] = "IL8"
colnames(dati_V)[colnames(dati_V) == "Cyfra 21-1"] = "Cyfra211"
colnames(dati_V)[colnames(dati_V) == "WHO-PS"] = "WHOPS"
colnames(dati_V)[colnames(dati_V) == "CA-9"] = "CA9"
colnames(dati_V)[colnames(dati_V) == "FEV1s%"] = "FEV1s"
colnames(dati_V)[colnames(dati_V) == "TLR-4"] = "TLR4"
names(dati_V)
View(dati_V)
dim(dati_V)
str(dati_V)


# Fattorizzazioni (e numerizzazioni! Alcuni biomarker erano
# segnati come character, colpa di valori segnati come "< qualcosa")
dati_V$histology = as.factor(dati_V$histology)
dati_V$stage = as.factor(dati_V$stage)
dati_V$Gender = as.factor(dati_V$Gender)
dati_V$RTProtocol = as.factor(dati_V$RTProtocol)
dati_V$IL6 = as.numeric(dati_V$IL6)
dati_V$IL8 = as.numeric(dati_V$IL8)
dati_V$Cyfra211 = as.numeric(dati_V$Cyfra211)
str(dati_V)
View(dati_V)


# Conversione dei valori "alive" e "dead" di dati_V di Status in 0 e 1
convert_status <- function(dati_V) {
  dati_V$Status <- ifelse(dati_V$Status == "alive", 0, 1)
  return(dati_V)
}
dati_V <- convert_status(dati_V)
View(dati_V)


# ANALISI ESPLORATIVA PER DATASET DI VALIDAZIONE--------------------------------
# Fattori
# Table per status (dataset di validazione)
table(dati_V$Status)

# Barplot per Gender (dataset di validazione)
ggplot(dati_V, aes(x = Gender)) + geom_bar(fill = c("steelblue", "pink")) + geom_text(stat = 'count', aes(label = ..count..), vjust = 2) +
  labs(title = "Distribuzione per genere", x = "Genere", y = "Frequenza")

# Barplot per histology (dataset di validazione)
ggplot(dati_V, aes(x = histology)) + geom_bar() + geom_text(stat = 'count', aes(label = ..count..), vjust = 2) +
  labs(title = "Distribuzione per istotipo", x = "Histology", y = "Frequenza")

# Barplot per stage (dataset di validazione)
ggplot(dati_V, aes(x = stage)) + geom_bar() + geom_text(stat = 'count', aes(label = ..count..), vjust = 2) +
  labs(title = "Distribuzione per stage tumorale", x = "Stage", y = "Frequenza")

# Barplot per RTProtocol (dataset di Validazione)
ggplot(dati_V, aes(x = RTProtocol)) + geom_bar() + geom_text(stat = 'count', aes(label = ..count..), vjust = 2) +
  labs(title = "Distribuzione per protocollo RT", x = "Protocollo", y = "Frequenza")

# Barplot per classi d'età (dataset di validazione)
# NOTA: D'ora in poi questa variabile verrà considerata solo come suddivisa in classi
dati_V$age <- cut(dati_V$age, breaks = c(40, 50, 60, 70, 80, 90), labels = c("40-50", "50-60", "60-70", "70-80", "80-90"))
dati_V$age <- as.factor(dati_V$age)
table(dati_V$age)
ggplot(dati_V, aes(x = age)) + 
  geom_bar(stat = "count", fill = "steelblue") +
  labs(title = "Distribuzione per classid'età", x = "Età", y = "Frequenza")

# Barplot per Linfonodi (accorpati in due classi: 0-1 e 2-3-4) (dataset di validazione)
# NOTA: D'ora in poi questa variabile verrà considerata solo come suddivisa in classi
bins = c(0, 1, 4)
labels_lymph <- c("0-1", "2-3-4")
dati_V$Nodes <- cut(dati_V$Nodes, bins, labels = labels_lymph)

ggplot(dati_V, aes(x = Nodes)) + 
  geom_bar(stat = "count") +
  labs(title = "Distribuzione per linfonodi", x = "Linfonodi", y = "Frequenza")

# Variabili quantitative, dataset di validazione (uni- e bi-variate mediante GGpairs)
View(dati_V)
dati_V_quantitative = dati_V[,c(8,9,10,12,13,14,15,16,17,18,19,20,21,22,23,24,25)]
ggpairs(dati_V_quantitative)


# ANALISI NON PARAMETRICA ------------------------------------------------------
# Curva di sopravvivenza generale (Validazione)
sopr_V = survfit(Surv(Survival, Status) ~ 1, data = dati_V)
ggsurvplot(sopr_V, data = dati_V, conf.int = TRUE)

# Curva di sopravvivenza per Sesso (Validazione)
sopr_gender_V = survfit(Surv(Survival, Status) ~ Gender, data = dati_V)
ggsurvplot(sopr_gender_V, data = dati_V, conf.int = TRUE, pval = TRUE, pval.method = TRUE)

# Curva di sopravvivenza per Terapia (RTProtocol) (Validazione)
sopr_rtprotocol_V = survfit(Surv(Survival, Status) ~ RTProtocol, data = dati_V)
ggsurvplot(sopr_rtprotocol_V, data = dati_V, conf.int = TRUE, pval = TRUE, pval.method = TRUE)

# Curva di sopravvivenza per Istologico (Validazione)
sopr_hist_V = survfit(Surv(Survival, Status) ~ histology, data = dati_V)
ggsurvplot(sopr_hist_V, data = dati_V, conf.int = TRUE, pval = TRUE, pval.method = TRUE)

# Curva di sopravvivenza per Stage (Validazione)
sopr_stage_V = survfit(Surv(Survival, Status) ~ stage, data = dati_V)
ggsurvplot(sopr_stage_V, data = dati_V, conf.int = F, pval = TRUE, pval.method = TRUE)

# Curva di sopravvivenza per WHOPS (Validazione)
sopr_whops_V = survfit(Surv(Survival, Status) ~ WHOPS, data = dati_V)
ggsurvplot(sopr_whops_V, data = dati_V, conf.int = TRUE, pval = TRUE, pval.method = TRUE)

# Curva di sopravvivenza per Linfonodi (Validazione)
sopr_nodes_V = survfit(Surv(Survival, Status) ~ Nodes, data = dati_V)
ggsurvplot(sopr_nodes_V, data = dati_V, conf.int = TRUE, pval = TRUE, pval.method = TRUE)

# Curva di sopravvivenza per classi d'eta (Validazione)
sopr_age_V <- survfit(Surv(Survival, Status) ~ age, data = dati_V)
ggsurvplot(sopr_age_V, data = dati_V, conf.int = FALSE, pval = TRUE, pval.method = TRUE)


# CORRELAZIONI -----------------------------------------------------------------
# correlazione tra variabili fisiche e biomarker
daticor_V = na.omit(dati_V[,c(8,9,10,14,15,16,17,18,19,20,21,22,23,24,25)])
str(daticor_V)
cormat_V <- round(cor(daticor_V),2)
melted_cormat_V <- melt(cormat_V)
ggplot(data = melted_cormat_V, aes(x=Var1, y=Var2, fill=value)) + geom_tile()


# MODELLI DI COX ---------------------------------------------------------------
# Modello di Cox standard
modcox_V_step = step(coxph(Surv(Survival, Status) ~ age + stage + histology + Gender + WHOPS + FEV1s + Nodes + RTProtocol + TotalDose1st + TotalDose2nd + GTV + OPN + CA9 + IL6 + IL8 + CRP + CEA + Cyfra211 + α2M + sIL2R + TLR4 + VEGF, data = na.omit(dati_V)))
modcox_V_final = coxph(Surv(Survival, Status) ~ Gender + Cyfra211 + CEA + FEV1s + histology + a2M + CA9 + TotalDose2nd + WHOPS + Nodes + TotalDose1st + age + OPN, data = na.omit(dati_V))
?coxph
summary(modcox_V_final)
test.ph.V = cox.zph(modcox_V_final)
test.ph.V
ggcoxzph(test.ph.V)


# Modello di Cox con frailty() usata su RTProtocol
modcox_V_step_RTP_frailty = step(coxph(Surv(Survival, Status) ~ age + stage + histology + Gender + WHOPS + FEV1s + Nodes + frailty(RTProtocol) + TotalDose1st + TotalDose2nd + GTV + OPN + CA9 + IL6 + IL8 + CRP + CEA + Cyfra211 + α2M + sIL2R + TLR4 + VEGF, data = na.omit(dati_V)))


# Modello di Cox con cluster() usata su RTProtocol
modcox_V_step_RTP_cluster = step(coxph(Surv(Survival, Status) ~ age + stage + histology + Gender + WHOPS + FEV1s + Nodes + cluster(RTProtocol) + TotalDose1st + TotalDose2nd + GTV + OPN + CA9 + IL6 + IL8 + CRP + CEA + Cyfra211 + α2M + sIL2R + TLR4 + VEGF, data = na.omit(dati_V)))
modcox_V_final_RTP_cluster = coxph(Surv(Survival, Status) ~ Gender + Cyfra211 + CEA + FEV1s + histology + a2M + CA9 + TotalDose2nd + WHOPS + Nodes + TotalDose1st + age + OPN + cluster(RTProtocol), data = na.omit(dati_V))
cox.zph(modcox_V_final_RTP_cluster)