# LIBRERIE
library(survival)
library(survminer)
library(timereg)
library(reshape)
library(reshape2)
library(ggplot2)
library(readxl)
library(MASS)
library(GGally)
library(lmtest)
library(mice)
library(naniar)
library(glmnet)
library(ggfortify)
# UPLOAD
dati_V = read_excel("carvalho-prognostic-biomarkers-NSCLC_IMP_Mea_Med_MLRI.xlsx", 2)
View(dati_V)
dim(dati_V)
table(dati_V$RTProtocol)
str(dati_V)
# RENAMING DELLE VARIABILI CON IL CARATTERE DI SPAZIO
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
colnames(dati_V)[colnames(dati_V) == "α2M"] = "a2M"
names(dati_V)
View(dati_V)
dim(dati_V)
str(dati_V)
# IMPUTAZIONE DEI VALORI MANCANTI (MEDIANA, MODA E M.I.C.E.)
sum(is.na(dati_V$age))
sum(is.na(dati_V$CRP))
sum(is.na(dati_V$CEA))
58
3.3 Modelling: fitting, selection and comparisons
sum(is.na(dati_V$stage))
find_mode <- function(x) {
u <- unique(x)
tab <- tabulate(match(x, u))
u[tab == max(tab)]
}
moderes = find_mode(dati_V$stage)
ii = which(is.na(dati_V$stage))
dati_V$stage[ii] = moderes
sum(is.na(dati_V$TLR4))
shapiro.test(dati_V$TLR4)
boxplot(dati_V$TLR4)
medres = median(na.omit(dati_V$TLR4))
ii =which(is.na(dati_V$TLR4))
dati_V$TLR4[ii] = medres
sum(is.na(dati_V$VEGF))
shapiro.test(dati_V$VEGF)
boxplot(dati_V$VEGF)
median(na.omit(dati_V$VEGF))
sum(is.na(dati_V$Cyfra211))
shapiro.test(as.numeric(dati_V$Cyfra211))
boxplot(dati_V$as.numeric(Cyfra211))
median(na.omit(dati_V$Cyfra211))
sum(is.na(dati_V$OPN))
shapiro.test(dati_V$OPN)
boxplot(dati_V$OPN)
median(na.omit(dati_V$OPN))
sum(is.na(dati_V$CA9))
shapiro.test(dati_V$CA9)
boxplot(dati_V$CA9)
median(na.omit(dati_V$CA9))
sum(is.na(dati_V$GTV))
shapiro.test(dati_V$GTV)
boxplot(dati_V$GTV)
median(na.omit(dati_V$GTV))
sum(is.na(dati_V$a2M))
shapiro.test(dati_V$a2M)
boxplot(dati_V$a2M)
median(na.omit(dati_V$a2M))
sum(is.na(dati_V$histology))
59
3.3 Modelling: fitting, selection and comparisons
find_mode(dati_V$histology)
sum(is.na(as.numeric(dati_V$IL8))) # Post numericizzazione per conversione "<" in NA
shapiro.test(as.numeric(dati_V$IL8))
boxplot(as.numeric(dati_V$IL8))
# inserimento valori arbitrari compresi tra 0 e x per "valore < x"
sum(is.na(as.factor(dati_V$WHOPS)))
find_mode(na.omit(as.factor(dati_V$WHOPS)))
sum(is.na(as.numeric(dati_V$IL6))) # Post numericizzazione per conversione "<" in NA
shapiro.test(as.numeric(dati_V$IL6))
boxplot(as.numeric(dati_V$IL6))
# inserimento valori arbitrari compresi tra 0 e x per "valore < x"
sum(is.na(dati_V$TotalDose2nd)) # **: valori NA/dati eccedenti il 43%
mcar_test(dati_V[,-13])
sum(is.na(dati_V$FEV1s))
shapiro.test(dati_V$FEV1s)
boxplot(dati_V$FEV1s)
View(dati_V)
dati_V$Status = as.factor(dati_V$Status)
dati_V$stage = as.factor(dati_V$stage)
dati_V$histology = as.factor(dati_V$histology)
dati_V$Gender = as.factor(dati_V$Gender)
dati_V$WHOPS = as.factor(dati_V$WHOPS)
dati_V$RTProtocol = as.factor(dati_V$RTProtocol)
imp = mice(dati_V[, -c(1,13,22)], method = c(rep("",7), "midastouch",
rep("",14)), print=FALSE, seed=1234)
names(imp)
dim(imp$data)
imp$imp
imp$imp$FEV1s
imputed_values_FEV1s <- mice::complete(imp)$FEV1s
imputed_values_FEV1s[3]
imputed_values_FEV1s[5]
imputed_values_FEV1s[14]
dati_V$FEV1s = imputed_values_FEV1s
sum(is.na(dati_V$FEV1s))
# FATTORIZZAZIONI E NUMERICIZZAZIONI
dati_V$histology = as.factor(dati_V$histology)
dati_V$stage = as.factor(dati_V$stage)
dati_V$Gender = as.factor(dati_V$Gender)
dati_V$RTProtocol = as.factor(dati_V$RTProtocol)
dati_V$WHOPS = as.factor(dati_V$WHOPS)
60
3.3 Modelling: fitting, selection and comparisons
dati_V$IL6 = as.numeric(dati_V$IL6)
dati_V$IL8 = as.numeric(dati_V$IL8)
dati_V$Cyfra211 = as.numeric(dati_V$Cyfra211)
str(dati_V)
View(dati_V)
sum(is.na(dati_V$`IL 8`))
# CONVERSIONE
# (dei valori "alive" e "dead" di dati_V di Status in 0 e 1)
convert_status <- function(dati_V) {
dati_V$Status <- ifelse(dati_V$Status == "alive", 0, 1)
return(dati_V)
}
dati_V <- convert_status(dati_V)
View(dati_V)
# CREAZIONE DELLE VARIABILI AGE E NODES SUDDIVISE IN CLASSI
# Age_cat
dati_V$age_cat <- cut(dati_V$age, breaks = c(40, 50, 60, 70, 80, 90), labels = c("40-50",
"50-60", "60-70", "70-80", "80-90"))
dati_V$age_cat <- as.factor(dati_V$age_cat)
# Nodes_cat
bins = c(0, 1, 4)
labels_lymph <- c("0-1", "2-3-4")
dati_V$Nodes_cat <- cut(dati_V$Nodes, bins, labels = labels_lymph)
# ANALISI ESPLORATIVA
# Fattori
# Table per status (dataset di validazione)
table(dati_V$Status)
# Barplot per Sesso
ggplot(dati_V, aes(x = Gender)) + geom_bar(fill = c("steelblue", "pink"))
+ geom_text(stat = 'count', aes(label = ..count..), vjust = 2) +
labs(title = "Distribuzione per genere", x = "Genere", y = "Frequenza")
# Barplot per histology (dataset di validazione)
ggplot(dati_V, aes(x = histology)) + geom_bar() + geom_text(stat = 'count',
aes(label = ..count..), vjust = 2) +
labs(title = "Distribuzione per istotipo", x = "Histology", y = "Frequenza")
# Barplot per stage (dataset di validazione)
ggplot(dati_V, aes(x = stage)) + geom_bar() + geom_text(stat = 'count',
aes(label = ..count..), vjust = 2) +
labs(title = "Distribuzione per stage tumorale", x = "Stage", y = "Frequenza")
# Barplot per RTProtocol (dataset di Validazione)
ggplot(dati_V, aes(x = RTProtocol)) + geom_bar() + geom_text(stat = 'count',
61
3.3 Modelling: fitting, selection and comparisons
aes(label = ..count..),vjust = 2) +
labs(title = "Distribuzione per protocollo RT", x = "Protocollo", y = "Frequenza")
# Barplot per classi d'età (dataset di validazione)
table(dati_V$age_cat)
ggplot(dati_V, aes(x = age_cat)) +
geom_bar(stat = "count", fill = "steelblue") +
labs(title = "Distribuzione per classi d'età", x = "Età", y = "Frequenza")
# Barplot per Linfonodi (accorpati in due classi: 0-1 e 2-3-4) (dataset di validazione)
# NOTA: D'ora in poi questa variabile verrà considerata solo come suddivisa in classi
table(dati_V$Nodes_cat)
ggplot(dati_V, aes(x = Nodes_cat)) +
geom_bar(stat = "count") +
labs(title = "Distribuzione per linfonodi", x = "Linfonodi", y = "Frequenza")
# Variabili quantitative, dataset di validazione (uni- e bi-variate mediante GGpairs)
View(dati_V)
dati_V_quantitative = dati_V[,c(8,9,10,12,13,14,15,16,17,18,19,20,21,22,23,24,25)]
ggpairs(dati_V_quantitative)
# ANALISI NON PARAMETRICA
# Curva di sopravvivenza generale
sopr_V = survfit(Surv(Survival, Status) ~ 1, data = dati_V)
ggsurvplot(sopr_V, data = dati_V, conf.int = TRUE)
# Curva di sopravvivenza per Sesso
sopr_gender_V = survfit(Surv(Survival, Status) ~ Gender, data = dati_V)
ggsurvplot(sopr_gender_V, data = dati_V, conf.int = TRUE, pval = TRUE, pval.method = TRUE)
# Curva di sopravvivenza per Terapia (RTProtocol)
sopr_rtprotocol_V = survfit(Surv(Survival, Status) ~ RTProtocol, data = dati_V)
ggsurvplot(sopr_rtprotocol_V, data = dati_V, conf.int = TRUE, pval = TRUE, pval.method = TRUE)
# Curva di sopravvivenza per Istologico
sopr_hist_V = survfit(Surv(Survival, Status) ~ histology, data = dati_V)
ggsurvplot(sopr_hist_V, data = dati_V, conf.int = TRUE, pval = TRUE, pval.method = TRUE)
# Curva di sopravvivenza per Stage
sopr_stage_V = survfit(Surv(Survival, Status) ~ stage, data = dati_V)
ggsurvplot(sopr_stage_V, data = dati_V, conf.int = F, pval = TRUE, pval.method = TRUE)
# Curva di sopravvivenza per WHOPS (**)
sopr_whops_V = survfit(Surv(Survival, Status) ~ WHOPS, data = dati_V)
ggsurvplot(sopr_whops_V, data = dati_V, conf.int = TRUE, pval = TRUE, pval.method = TRUE)
# Curva di sopravvivenza per Linfonodi
sopr_nodes_V = survfit(Surv(Survival, Status) ~ Nodes, data = dati_V)
ggsurvplot(sopr_nodes_V, data = dati_V, conf.int = TRUE, pval = TRUE, pval.method = TRUE)
62
3.3 Modelling: fitting, selection and comparisons
# Curva di sopravvivenza per Linfonodi categoriale
sopr_nodes_V = survfit(Surv(Survival, Status) ~ Nodes_cat, data = dati_V)
ggsurvplot(sopr_nodes_V, data = dati_V, conf.int = TRUE, pval = TRUE, pval.method = TRUE)
# Curva di sopravvivenza per classi d'eta
sopr_age_V <- survfit(Surv(Survival, Status) ~ age_cat, data = dati_V)
ggsurvplot(sopr_age_V, data = dati_V, conf.int = FALSE, pval = TRUE, pval.method = TRUE)
# CURVE DI SOPRAVVIVENZA DEI BIOMARKER
# opn
dati_V$opn_cat <- cut(dati_V$OPN, breaks = c(0, 75,
150, 225, 300), labels = c("0-75", "75-150", "150-225", "225-300"))
dati_V$opn_cat <- as.factor(dati_V$opn_cat)
sopr_opn_cat = survfit(Surv(Survival, Status) ~ opn_cat, data = dati_V)
ggsurvplot(sopr_opn_cat, data = dati_V, conf.int = F, pval = TRUE, pval.method = TRUE)
# ca9
dati_V$ca9_cat <- cut(dati_V$CA9, breaks = c(0,
250, 500, 750, 1000), labels = c("0-250", "250-500", "500-750", "750-1000"))
dati_V$ca9_cat <- as.factor(dati_V$ca9_cat)
sopr_ca9_cat = survfit(Surv(Survival, Status) ~ ca9_cat, data = dati_V)
ggsurvplot(sopr_ca9_cat, data = dati_V, conf.int = F, pval = TRUE, pval.method = TRUE, xlab)
# il6
dati_V$il6_cat <- cut(dati_V$IL6, breaks = c(0, 10,
20, 30, 40), labels = c("0-10", "10-20", "20-30", "30-40"))
dati_V$il6_cat <- as.factor(dati_V$il6_cat)
sopr_il6_cat = survfit(Surv(Survival, Status) ~ il6_cat, data = dati_V)
ggsurvplot(sopr_il6_cat, data = dati_V, conf.int = F, pval = TRUE, pval.method = TRUE)
# il8
dati_V$il8_cat <- cut(dati_V$IL8, breaks = c(0, 10, 20, 30, 40), labels = c("0-10", "10-20",
"20-30", "30-40"))
dati_V$il8_cat <- as.factor(dati_V$il8_cat)
sopr_il8_cat = survfit(Surv(Survival, Status) ~ il8_cat, data = dati_V)
ggsurvplot(sopr_il8_cat, data = dati_V, conf.int = F, pval = TRUE, pval.method = TRUE)
# il8
dati_V$crp_cat <- cut(dati_V$CRP, breaks = c(0, 20, 40, 60, 80), labels = c("0-20",
"20-40", "40-60", "60-80"))
dati_V$crp_cat <- as.factor(dati_V$crp_cat)
sopr_crp_cat = survfit(Surv(Survival, Status) ~ crp_cat, data = dati_V)
ggsurvplot(sopr_crp_cat, data = dati_V, conf.int = F, pval = TRUE, pval.method = TRUE)
# cea
dati_V$cea_cat <- cut(dati_V$CEA, breaks = c(0, 1.5, 3, 4.5, 6), labels = c("0-1.5",
"1.5-3", "3-4.5", "4.5-6"))
dati_V$cea_cat <- as.factor(dati_V$cea_cat)
63
3.3 Modelling: fitting, selection and comparisons
sopr_cea_cat = survfit(Surv(Survival, Status) ~ cea_cat, data = dati_V)
ggsurvplot(sopr_cea_cat, data = dati_V, conf.int = F, pval = TRUE, pval.method = TRUE)
# cyfra
dati_V$cyfra_cat <- cut(dati_V$Cyfra211, breaks = c(0, 1, 2, 3, 4), labels = c("0-1",
"1-2", "2-3", "3-4"))
dati_V$cyfra_cat <- as.factor(dati_V$cyfra_cat)
sopr_cyfra_cat = survfit(Surv(Survival, Status) ~ cyfra_cat, data = dati_V)
ggsurvplot(sopr_cyfra_cat, data = dati_V, conf.int = F, pval = TRUE, pval.method = TRUE)
# cyfra
dati_V$a2m_cat <- cut(dati_V$a2M, breaks = c(0, 1, 2, 3, 4, 5), labels = c("0-1",
"1-2", "2-3", "3-4", "4-5"))
dati_V$a2m_cat <- as.factor(dati_V$a2m_cat)
sopr_a2m_cat = survfit(Surv(Survival, Status) ~ a2m_cat, data = dati_V)
ggsurvplot(sopr_a2m_cat, data = dati_V, conf.int = F, pval = TRUE, pval.method = TRUE)
# sil2r
dati_V$sil_cat <- cut(dati_V$sIL2R, breaks = c(0, 3125, 6250, 9375, 12500),
labels = c("0-3125", "3125-6250", "6250-9375", "9375-12500"))
dati_V$sil_cat <- as.factor(dati_V$sil_cat)
sopr_sil_cat = survfit(Surv(Survival, Status) ~ sil_cat, data = dati_V)
ggsurvplot(sopr_sil_cat, data = dati_V, conf.int = F, pval = TRUE, pval.method = TRUE)
# tlr4
dati_V$tlr_cat <- cut(dati_V$TLR4, breaks = c(0, 5, 10, 15), labels = c("0-5",
"5-10", "10-15"))
dati_V$tlr_cat <- as.factor(dati_V$tlr_cat)
sopr_tlr_cat = survfit(Surv(Survival, Status) ~ tlr_cat, data = dati_V)
ggsurvplot(sopr_tlr_cat, data = dati_V, conf.int = F, pval = TRUE, pval.method = TRUE)
# vegf
dati_V$vegf_cat <- cut(dati_V$VEGF, breaks = c(0, 75, 150, 225, 300), labels = c("0-75",
"75-150", "150-225", "225-300"))
dati_V$vegf_cat <- as.factor(dati_V$vegf_cat)
sopr_vegf_cat = survfit(Surv(Survival, Status) ~ vegf_cat, data = dati_V)
ggsurvplot(sopr_vegf_cat, data = dati_V, conf.int = F, pval = TRUE, pval.method = TRUE)
# CORRELAZIONI
# correlazione tra variabili fisiche e biomarker
str(dati_V)
daticor_V = na.omit(dati_V[,c(9,10,14,15,16,17,18,19,20,21,22,23,24,25)])
str(daticor_V)
cormat_V <- round(cor(daticor_V),2)
melted_cormat_V <- melt(cormat_V)
ggplot(data = melted_cormat_V, aes(x=Var1, y=Var2, fill=value)) + geom_tile()
# MODELLI DI COX (STANDARD)
64
3.3 Modelling: fitting, selection and comparisons
# MODELLO COMPRENSIVO DI TUTTE LE VARIABILI
M1step = coxph(Surv(Survival, Status) ~ age + stage + histology + Gender + WHOPS + Nodes
+ RTProtocol + TotalDose1st + GTV + OPN + CA9 + IL6 + IL8 + CRP + CEA + Cyfra211 + a2M + sIL2R
+ TLR4 + VEGF + FEV1s, data = dati_V)
summary(M1step)
cox.zph(M1step, transform = "log")
cox.zph(M1step, transform = "km")
cox.zph(M1step, transform = "identity")
cox.zph(M1step, transform = "rank")
cox.zph(M1step)
M2 = coxph(Surv(Survival, Status) ~ Gender + CEA + histology + WHOPS + Nodes
+ RTProtocol + GTV + OPN + a2M + sIL2R
+ TLR4, data = dati_V)
summary(M2)
cox.zph(M2, transform = "log")
cox.zph(M2, transform = "km")
cox.zph(M2, transform = "identity")
cox.zph(M2, transform = "rank")
table(dati_V$RTProtocol, dati_V$Status)
#simulazione under null hypothesis
fit.cox<-cox.aalen(Surv(dati_V$Survival, dati_V$Status)~prop(dati_V$Gender)+prop(dati_V$CEA)
+prop(dati_V$histology)+prop(dati_V$WHOPS)+prop(dati_V$Nodes)+prop(dati_V$RTProtocol)
+prop(dati_V$GTV)+prop(dati_V$OPN)+prop(dati_V$a2M)+prop(dati_V$sIL2R)+prop(dati_V$TLR4),
weighted.test=0,pbc);
plot(fit.cox,xlab="Time (years)",ylab="Test process",score=T, specific.comps=2)
plot(fit.cox,xlab="Time (years)",ylab="Test process",score=T, specific.comps=9)
plot(fit.cox,xlab="Time (years)",ylab="Test process",score=T, specific.comps=16)
# VALUTAZIONE DELLE FORME FUNZIONALI DEL MODELLO DI BASE
# Res. martingala vs variabile
# residui martingala sul modello finale M2
residui_martingala <- resid(M2, type="martingale")
plot(dati_V$Nodes, residui_martingala, xlab="Nodes", ylab="Residui Martingala")
plot(dati_V$CEA, residui_martingala, xlab="CEA", ylab="Residui Martingala") # **
plot(dati_V$OPN, residui_martingala, xlab="OPN", ylab="Residui Martingala")
plot(dati_V$age, residui_martingala, xlab="age", ylab="Residui Martingala")
plot(dati_V$sIL2R, residui_martingala, xlab="sil2r", ylab="Residui Martingala")
plot(dati_V$GTV, residui_martingala, xlab="GTV", ylab="Residui Martingala") # **
plot(dati_V$Cyfra211, residui_martingala, xlab="cyfra211", ylab="Residui Martingala") #**
plot(dati_V$TLR4, residui_martingala, xlab="TLR4", ylab="Residui Martingala") #**
plot(dati_V$a2M, residui_martingala, xlab="a2M", ylab="Residui Martingala")
plot(dati_V$TotalDose1st, residui_martingala, xlab="dose1st", ylab="Residui Martingala") #**
plot(dati_V$IL6, residui_martingala, xlab="il6", ylab="Residui Martingala")
plot(dati_V$IL8, residui_martingala, xlab="il8", ylab="Residui Martingala")
plot(dati_V$CRP, residui_martingala, xlab="crp", ylab="Residui Martingala") #**
plot(dati_V$VEGF, residui_martingala, xlab="vegf", ylab="Residui Martingala")
ggcoxdiagnostics(M2, type = "martingale", ggtheme = theme_bw())
65
3.3 Modelling: fitting, selection and comparisons
# Residui di Schoenfeld
# Residui schoenfeld CEA (*)
Mcea = coxph(Surv(Survival, Status) ~ CEA, data = dati_V)
ccea = cox.zph(Mcea, transform = "km")
plot(ccea)
ccea = cox.zph(Mcea, transform = "rank")
plot(ccea)
ccea = cox.zph(Mcea, transform = "identity")
plot(ccea)
ccea = cox.zph(Mcea, transform = "log")
plot(ccea)
# Residui schoenfeld OPN
Mopn = coxph(Surv(Survival, Status) ~ OPN, data = dati_V)
copn = cox.zph(Mopn, transform = "km")
plot(copn)
copn = cox.zph(Mopn, transform = "rank")
plot(copn)
copn = cox.zph(Mopn, transform = "identity")
plot(copn)
copn = cox.zph(Mopn, transform = "log")
plot(copn)
# Residui schoenfeld VEGF
Mvegf = coxph(Surv(Survival, Status) ~ VEGF, data = dati_V)
cvegf = cox.zph(Mvegf)
plot(cvegf)
# Residui schoenfeld Nodes (**)
Mnodes = coxph(Surv(Survival, Status) ~ Nodes, data = dati_V)
cnodes = cox.zph(Mnodes, transform = "km")
plot(cnodes)
cnodes = cox.zph(Mnodes, transform = "rank")
plot(cnodes)
cnodes = cox.zph(Mnodes, transform = "identity")
plot(cnodes)
cnodes = cox.zph(Mnodes, transform = "log")
plot(cnodes)
# Residui schoenfeld Cyfra211 (*)
Mcyf = coxph(Surv(Survival, Status) ~ Cyfra211, data = dati_V)
ccyf = cox.zph(Mcyf, transform = "km")
plot(ccyf)
ccyf = cox.zph(Mcyf, transform = "identity")
plot(ccyf)
ccyf = cox.zph(Mcyf, transform = "rank")
plot(ccyf)
ccyf = cox.zph(Mcyf, transform = "log")
66
3.3 Modelling: fitting, selection and comparisons
plot(ccyf)
# Residui schoenfeld WHOPS (*)
Mwhops = coxph(Surv(Survival, Status) ~ WHOPS, data = dati_V)
cwhops = cox.zph(Mwhops)
plot(cwhops)
# Residui schoenfeld Gender (*)
Mgen = coxph(Surv(Survival, Status) ~ Gender, data = dati_V)
cgen = cox.zph(Mgen)
plot(cgen)
# Residui schoenfeld TLR4
Mtlr4 = coxph(Surv(Survival, Status) ~ TLR4, data = dati_V)
ctlr4 = cox.zph(Mtlr4)
plot(ctlr4)
# Residui schoenfeld histology (*)
Mhist = coxph(Surv(Survival, Status) ~ histology, data = dati_V)
chist = cox.zph(Mhist)
plot(chist)
# Residui schoenfeld a2M
Ma2m = coxph(Surv(Survival, Status) ~ a2M, data = dati_V)
ca2M = cox.zph(Ma2m)
plot(ca2M)
# Residui schoenfeld RTProtocol (**)
Mrtp = coxph(Surv(Survival, Status) ~ RTProtocol, data = dati_V)
crtp = cox.zph(Mrtp)
plot(crtp)
?cox.zph
# Residui schoenfeld sIL2R
Msil = coxph(Surv(Survival, Status) ~ sIL2R, data = dati_V)
csil = cox.zph(Msil)
plot(csil)
# Residui schoenfeld GTV (*)
Mgtv = coxph(Surv(Survival, Status) ~ GTV, data = dati_V)
cgtv = cox.zph(Mgtv)
plot(cgtv)
# Residui schoenfeld Age
Mage = coxph(Surv(Survival, Status) ~ age, data = dati_V)
cage = cox.zph(Mage)
plot(cage)
# Residui schoenfeld TotalDose1st
67
3.3 Modelling: fitting, selection and comparisons
Mttd = coxph(Surv(Survival, Status) ~ TotalDose1st, data = dati_V)
cttd = cox.zph(Mttd)
plot(cttd)
# Residui schoenfeld IL6
Mil6 = coxph(Surv(Survival, Status) ~ IL6, data = dati_V)
cil6 = cox.zph(Mil6)
plot(cil6)
# Residui schoenfeld IL8
Mil8 = coxph(Surv(Survival, Status) ~ IL8, data = dati_V)
cil8 = cox.zph(Mil8)
plot(cil8)
# Residui schoenfeld CA9
Mca9 = coxph(Surv(Survival, Status) ~ CA9, data = dati_V)
cca9 = cox.zph(Mca9)
plot(cca9)
# Residui schoenfeld CRP
Mcrp = coxph(Surv(Survival, Status) ~ CRP, data = dati_V)
crp = cox.zph(Mcrp)
plot(crp)
# Residui schoenfeld FEV1s
Mfev = coxph(Surv(Survival, Status) ~ FEV1s, data = dati_V)
cfev = cox.zph(Mfev)
plot(cfev)
# curve di sopravvivenza stimate dai due modelli
soprM2 = survfit(M2, data = dati_V)
ggsurvplot(soprM2, data = dati_V, conf.int = TRUE, pval = TRUE, pval.method = TRUE)
sopr_V = survfit(Surv(Survival, Status) ~ 1, data = dati_V)
ggsurvplot(sopr_V, data = dati_V, conf.int = TRUE)
# MODELLO DI COX ESTESO
M2new = coxph(Surv(Survival, Status) ~ Gender + CEA + histology + WHOPS + Nodes
+ RTProtocol + pspline(GTV, df = 0) + OPN + a2M + sIL2R + TLR4, data = dati_V)
summary(M2new)
soprM2new = survfit(M2new, data = dati_V)
termplot(M2new, term=7, se=TRUE, col.term=1, col.se=1)
termplot(M2new, term=11, se=TRUE, col.term=1, col.se=1)
ptemp <- termplot(M2, se=TRUE, plot=FALSE)
attributes(ptemp)
termplot(M2new, term=11, se=TRUE, col.term=1, col.se=1)
ptemp <- termplot(M2, se=TRUE, plot=FALSE)
68
3.3 Modelling: fitting, selection and comparisons
attributes(ptemp)
?pspline
# VALUTAZIONE DELLE FORME FUNZIONALI DEL MODELLO A EFFETTI T.DIP.
res = resid(M2tt, "martingale")
ggcoxdiagnostics(M2tt, type = "martingale", ggtheme = theme_bw())
plot(M2new,xlab="Time (years)",ylab="Test process",score=T, specific.comps=10)
mm = timecox(Surv(Survival, Status) ~ Gender + CEA + histology + WHOPS + Nodes + RTProtocol
+ GTV + OPN + a2M, dati_V)
summary(mm)
names(dati_V)
# MODELLI FLESSIBILI: MODELLO A RISCHI ADDITIVI DI AALEN
library(timereg)
fit_aa = aalen(Surv(Survival, Status) ~ GTV + age + stage + histology + Gender
+ Nodes + RTProtocol
+ CRP + CEA + Cyfra211 + a2M + TLR4 + VEGF, data = dati_V)
summary(fit_aa)
# MODELLO DI MCKEAGUE E SASIENI
fit_ms = aalen(Surv(Survival, Status) ~ const(GTV) + const(age)
+ stage + const(histology) + const(Gender) + Nodes + RTProtocol + CRP + CEA + const(Cyfra211)
+ const(a2M) + const(TLR4) + const(VEGF), data = dati_V)
summary(fit_ms)
plot(fit_ms, what = "survival")
plot(fit_aa,xlab="Time (years)",ylab="Test process",score=T, specific.comps=10)
plot(fit_aa,xlab="Time (years)",ylab="Test process",score=T, specific.comps=13)
plot(fit_aa,xlab="Time (years)",ylab="Test process",score=T, specific.comps=14)
plot(fit_aa,xlab="Time (years)",ylab="Test process",score=T, specific.comps=4)
plot(fit_aa,xlab="Time (years)",ylab="Test process",score=T, specific.comps=12)
plot(fit_ms,score=T,xlab="Time (years)",ylab="Test process")
# MODELLO FINALE DI MCKEAGUE E SASIENI
fit_ms2 = aalen(Surv(Survival, Status) ~ const(GTV) + const(age) + const(stage)
+ const(histology) + const(Gender) + Nodes + const(RTProtocol) + const(CRP) + CEA
+ const(Cyfra211) + const(a2M) + const(TLR4) + const(VEGF),
data = dati_V)
summary(fit_ms2)
# MODELLO MOLTIPLICATIVO ADDITIVO DI COX AALEN
fit_ca = cox.aalen(Surv(Survival, Status) ~ prop(WHOPS) + prop(GTV) + prop(age) + prop(stage)
+ prop(histology) + prop(Gender) + Nodes + prop(RTProtocol) + prop(CRP)
+ CEA + prop(Cyfra211) + prop(a2M) + prop(TLR4) +
prop(VEGF), max.time = 8,Nit = 1000, dati_V)
summary(fit_ca)
cox.surv<-list(time=fit_ca$cum[,1],surv=exp(-fit_ca$cum[,2]))
lines(cox.surv$time,cox.surv$surv,type="s",lwd=2,lty=2)
69
3.3 Modelling: fitting, selection and comparisons
plot(fit_ca)
plot(fit_ca,score=T,xlab="Time (years)")
fit = aalen(Surv(Survival, Status) ~ const(GTV) + const(age) + const(stage) + const(histology)
+ const(Gender) + Nodes + const(RTProtocol) + const(CRP) + CEA + const(Cyfra211)
+ const(a2M) + const(TLR4) + const(VEGF), data = dati_V,max.time=8, resample.iid=1)
x0<-c(0,0,1); z0<-c(1,0,0);
delta<-matrix(0,length(fit$cum[,1]),181)
for (i in 1:181) {delta[,i]<-x0%*%t(fit$B.iid[[i]])+fit$cum[,1]*sum(z0*fit$gamma.iid[i,]);}
S0<-exp(- x0 %*% t(fit$cum[,-1])- fit$cum[,1]*sum(z0*fit$gamma))
se<-apply(delta^2,1,sum)^.5
plot(fit$cum[,1],S0,type="l",ylim=c(0,1),xlab="Time (years)",ylab="Survival")
fit_ca_s<-cox.aalen(Surv(Survival, Status) ~ prop(WHOPS) + prop(GTV) + prop(age) + prop(stage)
+ prop(histology) + prop(Gender) + Nodes + prop(RTProtocol) + prop(CRP) + CEA
+ prop(Cyfra211) + prop(a2M) + prop(TLR4) + prop(VEGF),
data = dati_V,max.time=8, resample.iid=1)
x0<-c(0,0,1); z0<-c(1,0,0);
delta<-matrix(0,length(fit_ca_s$cum[,1]),181)
for (i in 1:181) {delta[,i]<-x0%*%t(fit_ca_s$B.iid[[i]])+fit_ca_s$cum[,1]
*sum(z0*fit_ca_s$gamma.iid[i,]);}
S0<-exp(- x0 %*% t(fit_ca_s$cum[,-1])- fit_ca_s$cum[,1]*sum(z0*fit_ca_s$gamma))
se_ca_s<-apply(delta^2,1,sum)^.5
surv_ms2 = aalen(Surv(Survival, Status) ~ const(Gender) + const(RTProtocol)
+ const(WHOPS) + const(GTV) + const(stage) + const(histology) + const(Cyfra211)
+ const(a2M) + const(TLR4) + Nodes + CEA, resample.iid = 1, data = dati_V)
summary(surv_ms2)
sur_ms2 = predict.aalen(surv_ms2, dati_V, uniform = F, unif.bands = F)
plot(sur_ms2, col = "blue", ylab = "Survival", xlab = "Time (Years)")
surv_ca = cox.aalen(Surv(Survival, Status) ~ prop(Gender) + prop(RTProtocol)
+ prop(WHOPS) + prop(GTV) + prop(stage) + prop(histology) + prop(Cyfra211) + prop(a2M)
+ prop(TLR4) + Nodes + CEA, resample.iid = 1, data = dati_V)
summary(surv_ca)
sur_ca = predict.cox.aalen(surv_ca, dati_V, uniform = F, unif.bands = F)
plot(sur_ca, col = "red", ylab = "Survival", xlab = "Time (Years)")
?predict.aalen
M2 = coxph(Surv(Survival, Status) ~ Gender + RTProtocol + WHOPS + GTV
+ stage + histology + Cyfra211 + a2M + TLR4 + Nodes + CEA, data = dati_V)
soprM2 = survfit(M2, data = dati_V)
plot(soprM2)
lines(soprM2, conf.int = T, col = "black", lwd = 1)
lines(soprM2, conf.int = F, col = "black", lwd = 2)
lines(soprM2new, col = "orange", lwd = 1, conf.int = T)
lines(soprM2new, col = "orange", lwd = 2, conf.int = F)
