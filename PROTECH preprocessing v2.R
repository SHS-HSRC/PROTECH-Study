library(readxl)
library(comorbidity)
library(tidyverse)
library(caTools)
library(mltools)
library(caret)
library(ROCR) #Compute ROC-AUC
library(PRROC)
library(InformationValue) #For Plotting ROC-AUC
library(lubridate)
#library(MASS) #For robust regression but conflicts with dplyr
setwd("C:/Users/AI-Ball/Desktop/PROTECH")

#------------#
#---Cohort---#
#------------#
data <- as.data.frame(read_excel('Demographics_DCI__deid.xlsx'))
colnames(data)
data <- data[c('NRIC','DOB','Gender','Race','Age_At_Dx','FinalStage', 'A (Stage III/IV during study period)', 'B',
               'Recoded_ICD10','Hist_Grade', 'c_T','c_N','c_M','p_T','p_N','p_M','DeathDate','Diag_Date','Cons_Date')]

#Diagnosis Categorizing
data_diag <- read_excel('Demographics_DCI__deid.xlsx') %>% 
  filter(`A (Stage III/IV during study period)`=="A") %>% 
  select(NRIC, Recoded_ICD10, Dx) %>% 
  filter(NRIC %in% data$NRIC)
data_diag$Recoded_ICD10 <- as.character(data_diag$Recoded_ICD10)
data_diag$clean_ICD10 <- sapply(strsplit(data_diag$Recoded_ICD10,"\\."), `[`, 1)
data_diag$alphabet_ICD10 <- gsub("[^a-zA-Z]", "", data_diag$clean_ICD10)
data_diag$numeric_ICD10 <- as.numeric(gsub("[^0-9]", "", data_diag$clean_ICD10))

data_diag$category <- NA
data_diag[data_diag$numeric_ICD10<=14 & data_diag$alphabet_ICD10=="C",]$category <- "C00-C14"
data_diag[data_diag$numeric_ICD10>=15 & data_diag$numeric_ICD10<=21 & data_diag$alphabet_ICD10=="C",]$category <- "C15-C21"
data_diag[data_diag$numeric_ICD10>=22 & data_diag$numeric_ICD10<=26 & data_diag$alphabet_ICD10=="C",]$category <- "C22-C26"
data_diag[data_diag$numeric_ICD10>=30 & data_diag$numeric_ICD10<=39 & data_diag$alphabet_ICD10=="C",]$category <- "C30-C39"
data_diag[data_diag$numeric_ICD10>=40 & data_diag$numeric_ICD10<=41 & data_diag$alphabet_ICD10=="C",]$category <- "C40-C41"
data_diag[data_diag$numeric_ICD10>=43 & data_diag$numeric_ICD10<=44 & data_diag$alphabet_ICD10=="C",]$category <- "C43-C44"
data_diag[data_diag$numeric_ICD10>=45 & data_diag$numeric_ICD10<=49 & data_diag$alphabet_ICD10=="C",]$category <- "C45-C49"
data_diag[data_diag$numeric_ICD10==50 & data_diag$alphabet_ICD10=="C",]$category <- "C50"
data_diag[data_diag$numeric_ICD10>=51 & data_diag$numeric_ICD10<=58 & data_diag$alphabet_ICD10=="C",]$category <- "C51-C58"
data_diag[data_diag$numeric_ICD10>=60 & data_diag$numeric_ICD10<=63 & data_diag$alphabet_ICD10=="C",]$category <- "C60-C63"
data_diag[data_diag$numeric_ICD10>=64 & data_diag$numeric_ICD10<=68 & data_diag$alphabet_ICD10=="C",]$category <- "C64-C68"
data_diag[data_diag$numeric_ICD10>=69 & data_diag$numeric_ICD10<=72 & data_diag$alphabet_ICD10=="C",]$category <- "C69-C72"
data_diag[data_diag$numeric_ICD10>=73 & data_diag$numeric_ICD10<=75 & data_diag$alphabet_ICD10=="C",]$category <- "C73-C75"
data_diag[data_diag$numeric_ICD10>=76 & data_diag$numeric_ICD10<=80 & data_diag$alphabet_ICD10=="C",]$category <- "C76-C80"

data_diag <- data_diag %>% 
  filter(!is.na(category)) %>% 
  select(NRIC, category) %>%
  distinct()

dummy <- dummyVars("NRIC ~ category", data=data_diag)
final_df <- data.frame(predict(dummy, newdata=data_diag))
data_diag <- cbind(data_diag, final_df)
data_diag <- data_diag %>%  
  select(-category) %>% 
  group_by(NRIC) %>% 
  summarise(across(everything(), sum)) #Condense multiple entries by patient into 1 row
data <- left_join(data, data_diag, by='NRIC')

#Race
data$Race <- as.character(data$Race)
data$Race <- ifelse(data$Race!='CHINESE'&data$Race!='MALAY'&data$Race!='INDIAN', 'OTHERS', data$Race)

#Stage
stage3 <- c('III','IIIA','IIIA1','IIIA1i','IIIA1ii','IIIA2','IIIAS','IIIB','IIIC','IIIC1','IIIC2','IIICNOS','IIID','IIINOS','Stage III')
stage4 <- c('IV','IVA','IVA2','IVB','IVC','IVNOS','Stage IVA','Stage IVB')
data$stage <- NA
data[tolower(data$FinalStage) %in% tolower(stage3),]$stage <- 'Stage 3'
data[tolower(data$FinalStage) %in% tolower(stage4),]$stage <- 'Stage 4'

#Histology Grade
data$Hist_Grade <- ifelse(data$Hist_Grade=='GX', 9, data$Hist_Grade)
data$Hist_Grade <- gsub("[^0-9.-]", "", data$Hist_Grade)
data[!is.na(data$Hist_Grade),]$Hist_Grade <- paste0('Grade_',data[!is.na(data$Hist_Grade),]$Hist_Grade)

#TNM
convert_T <- function(data, col) {
  data[grepl('1',data[,col]),][,col] <- 'T1'
  data[grepl('2',data[,col]),][,col] <- 'T2'
  data[grepl('3',data[,col]),][,col] <- 'T3'
  data[grepl('4',data[,col]),][,col] <- 'T4'
  data[grepl('0',data[,col]),][,col] <- 'T0'
  data[grepl('X',toupper(data[,col])),][,col] <- 'TX'
  data[grepl('is',data[,col]),][,col] <- 'Tis'
  data[grepl('Ta',data[,col]),][,col] <- 'Tis'
  
  return(data[,col])
}

convert_N <- function(data, col) {
  data[grepl('1',data[,col]),][,col] <- 'N1'
  data[grepl('2',data[,col]),][,col] <- 'N2'
  data[grepl('3',data[,col]),][,col] <- 'N3'
  data[grepl('0',data[,col]),][,col] <- 'N0'
  data[grepl('X',toupper(data[,col])),][,col] <- 'NX'
  
  return(data[,col])
}

convert_M <- function(data, col) {
  data[grepl('1',data[,col]),][,col] <- 'M1'
  data[grepl('0',data[,col]),][,col] <- 'M0'
  data[grepl('X',toupper(data[,col])),][,col] <- 'MX'
  
  return(data[,col])
}

data$c_T <- convert_T(data, 'c_T')
data$c_N <- convert_N(data, 'c_N')
data$c_M <- convert_M(data, 'c_M')
data$p_T <- convert_T(data, 'p_T')
data$p_N <- convert_N(data, 'p_N')
data$p_M <- convert_M(data, 'p_M')

#Filtered by Dr QY on rows which duplicated same stage but different TNM rows are retained
filter_diag <- as.data.frame(read_excel('duplicated_stage_different_tnm (qy tagged).xlsx')) 
filter_diag <- filter_diag %>%
  filter(chosen_diagnosis==1) %>% 
  select(NRIC, Recoded_ICD10, c_T, c_N, c_M, p_T, p_N, p_M)
data1 <- data[do.call(paste0, data[c('NRIC','Recoded_ICD10', 'c_T', 'c_N', 'c_M', 'p_T', 'p_N', 'p_M')]) %in% do.call(paste0, filter_diag),]
data <- data[!(data$NRIC %in% data1$NRIC),]

#Select cohort A and remove duplicates
data <- data %>% 
  group_by(NRIC) %>% 
  filter(`A (Stage III/IV during study period)`=="A") %>% 
  arrange(desc(stage), .by_group = TRUE) %>%
  slice(1) %>% ungroup()
data <- rbind(data, data1)

#---------#
#---Lab---#
#---------#
cohort <- as.data.frame(read_excel('Demographics_DCI__deid.xlsx')) %>% 
  select(NRIC, Gender, DOB) %>% 
  distinct()

lab <- read.csv('PROTECH_Lab Results (Sensitive Data Excluded)_full_deid.csv')
lab <- left_join(lab, cohort, by=c("Patient.ID"="NRIC"))
lab$DOB <- as.Date(lab$DOB)
lab$Specimen.Received.Date <- as.Date(lab$Specimen.Received.Date)
lab$age_outpx <- floor(time_length(difftime(lab$Specimen.Received.Date, lab$DOB),"years"))
lab$Result.Value <- as.numeric(lab$Result.Value)
lab <- lab[!is.na(lab$Result.Value),] 

lab <- lab %>% 
  filter(grepl("^wbc$|white blood cell|white blood cells|haemoglobin|platelet count|^creatinine$|^albumin$|bilirubin total|bilirubin,total|^lymph absolute$|^neut absolute$|protein,total", 
               tolower(Lab.Resulted.Order.Test.Description)))
lab[grepl("^wbc$|white blood cell", tolower(lab$Lab.Resulted.Order.Test.Description)),]$Lab.Resulted.Order.Test.Description <- 'WBC'
lab[grepl("haemoglobin", tolower(lab$Lab.Resulted.Order.Test.Description)),]$Lab.Resulted.Order.Test.Description <- 'Haemoglobin'
lab[grepl("platelet count", tolower(lab$Lab.Resulted.Order.Test.Description)),]$Lab.Resulted.Order.Test.Description <- 'Platelet Count'
lab[grepl("^creatinine$", tolower(lab$Lab.Resulted.Order.Test.Description)),]$Lab.Resulted.Order.Test.Description <- 'Creatinine'
lab[grepl("^albumin$", tolower(lab$Lab.Resulted.Order.Test.Description)),]$Lab.Resulted.Order.Test.Description <- 'Albumin'
lab[grepl("bilirubin total|bilirubin,total", tolower(lab$Lab.Resulted.Order.Test.Description)),]$Lab.Resulted.Order.Test.Description <- 'Bilirubin Total'
lab[grepl("protein,total", tolower(lab$Lab.Resulted.Order.Test.Description)),]$Lab.Resulted.Order.Test.Description <- 'Protein Total'

#-- Transform Creatinine to Creatinine clearance --#
#Conversion factor 88.42 for umol/L to mg/dL
lab[lab$Lab.Resulted.Order.Test.Description=='Creatinine',]$Result.Value <- lab[lab$Lab.Resulted.Order.Test.Description=='Creatinine',]$Result.Value/88.42
#Female & Scr>0.7
lab[lab$Lab.Resulted.Order.Test.Description=='Creatinine'&lab$Gender=='Female'&lab$Result.Value>0.7,]$Result.Value <-  
  142*(lab[lab$Lab.Resulted.Order.Test.Description=='Creatinine'&lab$Gender=='Female'&lab$Result.Value>0.7,]$Result.Value/0.7)^(-1.2)*(0.9938^lab[lab$Lab.Resulted.Order.Test.Description=='Creatinine'&lab$Gender=='Female'&lab$Result.Value>0.7,]$age_outpx)*1.012
#Male & Scr>0.9
lab[lab$Lab.Resulted.Order.Test.Description=='Creatinine'&lab$Gender=='Male'&lab$Result.Value>0.9,]$Result.Value <-  
  142*(lab[lab$Lab.Resulted.Order.Test.Description=='Creatinine'&lab$Gender=='Male'&lab$Result.Value>0.9,]$Result.Value/0.9)^(-1.2)*(0.9938^lab[lab$Lab.Resulted.Order.Test.Description=='Creatinine'&lab$Gender=='Male'&lab$Result.Value>0.9,]$age_outpx)
#Female & Scr<=0.7
lab[lab$Lab.Resulted.Order.Test.Description=='Creatinine'&lab$Gender=='Female'&lab$Result.Value<=0.7,]$Result.Value <-  
  142*(lab[lab$Lab.Resulted.Order.Test.Description=='Creatinine'&lab$Gender=='Female'&lab$Result.Value<=0.7,]$Result.Value/0.7)^(-0.241)*(0.9938^lab[lab$Lab.Resulted.Order.Test.Description=='Creatinine'&lab$Gender=='Female'&lab$Result.Value<=0.7,]$age_outpx)*1.012
#Male & Scr<=0.9
lab[lab$Lab.Resulted.Order.Test.Description=='Creatinine'&lab$Gender=='Male'&lab$Result.Value<=0.9,]$Result.Value <- 
  142*(lab[lab$Lab.Resulted.Order.Test.Description=='Creatinine'&lab$Gender=='Male'&lab$Result.Value<=0.9,]$Result.Value/0.9)^(-0.302)*(0.9938^lab[lab$Lab.Resulted.Order.Test.Description=='Creatinine'&lab$Gender=='Male'&lab$Result.Value<=0.9,]$age_outpx)

NL_ratio <- lab %>%
  group_by(Patient.ID,Specimen.Received.Date) %>% 
  summarise(NL_ratio = Result.Value[Lab.Resulted.Order.Test.Description == "NEUT ABSOLUTE"] / 
              Result.Value[Lab.Resulted.Order.Test.Description == "LYMPH ABSOLUTE"]) %>% 
  ungroup()

AG_ratio <- lab %>%
  group_by(Patient.ID,Specimen.Received.Date) %>% 
  summarise(AG_ratio = Result.Value[Lab.Resulted.Order.Test.Description == "Albumin"] /
              (Result.Value[Lab.Resulted.Order.Test.Description == "Protein Total"]-Result.Value[Lab.Resulted.Order.Test.Description == "Albumin"])) %>% 
  ungroup()

table(lab$Lab.Resulted.Order.Test.Description)

#-----------#
#---OutPx---#
#-----------#
outpx <- as.data.frame(read_excel('PROTECH_Outpatient_NCC_Complete__deid.xlsx')) 
outpx <- outpx[outpx$`Patient ID` %in% data$NRIC,] 
outpx$`Visit Date` <- as.Date(outpx$`Visit Date`)
outpx <- outpx %>% 
  filter(grepl("Actualised|Walk-in", `Visit Status Description`)) %>%
  filter(grepl("Follow-up|New Case", `Visit Type Description`)) %>% 
  group_by(`Patient ID`, `Visit Date`) %>% 
  slice(1) %>% 
  arrange(`Visit Date`, .by_group = TRUE) %>% 
  ungroup() %>% 
  select(`Patient ID`, `Visit Date`)

## Albumin ##
for (i in 1:nrow(outpx)) {
  #1 year look back
  lab1 <- lab[lab$Patient.ID==outpx$`Patient ID`[i],]
  subset_1y <- subset(lab1, (Specimen.Received.Date>=outpx$`Visit Date`[i]-365 & Specimen.Received.Date<= (outpx$`Visit Date`[i]-1)) &
                        lab1$Lab.Resulted.Order.Test.Description=='Albumin')
  outpx$Albumin_med[i] <- median(subset_1y$Result.Value)
  outpx$Albumin_mean[i] <- mean(subset_1y$Result.Value)
  outpx$Albumin_max[i] = ifelse(is.infinite(max(subset_1y$Result.Value)), NA, max(subset_1y$Result.Value))
  outpx$Albumin_min[i] = ifelse(is.infinite(min(subset_1y$Result.Value)), NA, min(subset_1y$Result.Value))
  outpx$Albumin_sd[i] <- sd(subset_1y$Result.Value)
  outpx$Albumin_grad[i] <- tryCatch(rlm(Result.Value~Specimen.Received.Date, data=subset_1y)$coefficients['Specimen.Received.Date'], error=function(err) NA)
  outpx$Albumin_lat[i] <- (subset_1y %>% arrange(desc(Specimen.Received.Date)))$Result.Value[1]
  cat("Albumin current status:", i, '/', nrow(outpx), '\n')
}
save.image("1.RData")

#For index 141,626,748 - See difference between lm and rlm
#plot(subset_1y$Specimen.Received.Date, subset_1y$Result.Value)
#abline(lm(Result.Value~Specimen.Received.Date, data=subset_1y))
#abline(rlm(Result.Value~Specimen.Received.Date, data=subset_1y), col="blue")

## Bilirubin ##
for (i in 1:nrow(outpx)) {
  #1 year look back
  lab1 <- lab[lab$Patient.ID==outpx$`Patient ID`[i],]
  subset_1y <- subset(lab1, (Specimen.Received.Date>=outpx$`Visit Date`[i]-365 & Specimen.Received.Date<= (outpx$`Visit Date`[i]-1)) &
                        lab1$Lab.Resulted.Order.Test.Description=='Bilirubin Total')
  outpx$Bilirubin_med[i] <- median(subset_1y$Result.Value)
  outpx$Bilirubin_mean[i] <- mean(subset_1y$Result.Value)
  outpx$Bilirubin_max[i] = ifelse(is.infinite(max(subset_1y$Result.Value)), NA, max(subset_1y$Result.Value))
  outpx$Bilirubin_min[i] = ifelse(is.infinite(min(subset_1y$Result.Value)), NA, min(subset_1y$Result.Value))
  outpx$Bilirubin_sd[i] <- sd(subset_1y$Result.Value)
  outpx$Bilirubin_grad[i] <- tryCatch(rlm(Result.Value~Specimen.Received.Date, data=subset_1y)$coefficients['Specimen.Received.Date'], error=function(err) NA)
  outpx$Bilirubin_lat[i] <- (subset_1y %>% arrange(desc(Specimen.Received.Date)))$Result.Value[1]
  
  cat("Bilirubin current status:", i, '/', nrow(outpx), '\n')
}
save.image("2.RData")

## Creatinine Clearance ##
for (i in 1:nrow(outpx)) {
  #1 year look back
  lab1 <- lab[lab$Patient.ID==outpx$`Patient ID`[i],]
  subset_1y <- subset(lab1, (Specimen.Received.Date>=outpx$`Visit Date`[i]-365 & Specimen.Received.Date<= (outpx$`Visit Date`[i]-1)) &
                        lab1$Lab.Resulted.Order.Test.Description=='Creatinine')
  outpx$Creatinine_med[i] <- median(subset_1y$Result.Value)
  outpx$Creatinine_mean[i] <- mean(subset_1y$Result.Value)
  outpx$Creatinine_max[i] = ifelse(is.infinite(max(subset_1y$Result.Value)), NA, max(subset_1y$Result.Value))
  outpx$Creatinine_min[i] = ifelse(is.infinite(min(subset_1y$Result.Value)), NA, min(subset_1y$Result.Value))
  outpx$Creatinine_sd[i] <- sd(subset_1y$Result.Value)
  outpx$Creatinine_grad[i] <- tryCatch(rlm(Result.Value~Specimen.Received.Date, data=subset_1y)$coefficients['Specimen.Received.Date'], error=function(err) NA)
  outpx$Creatinine_lat[i] <- (subset_1y %>% arrange(desc(Specimen.Received.Date)))$Result.Value[1]
  
  cat("Creatinine current status:", i, '/', nrow(outpx), '\n')
}
save.image("3.RData")

## Haemoglobin ##
for (i in 1:nrow(outpx)) {
  #1 year look back
  lab1 <- lab[lab$Patient.ID==outpx$`Patient ID`[i],]
  subset_1y <- subset(lab1, (Specimen.Received.Date>=outpx$`Visit Date`[i]-365 & Specimen.Received.Date<= (outpx$`Visit Date`[i]-1)) &
                        lab1$Lab.Resulted.Order.Test.Description=='Haemoglobin')
  outpx$Haemoglobin_med[i] <- median(subset_1y$Result.Value)
  outpx$Haemoglobin_mean[i] <- mean(subset_1y$Result.Value)
  outpx$Haemoglobin_max[i] = ifelse(is.infinite(max(subset_1y$Result.Value)), NA, max(subset_1y$Result.Value))
  outpx$Haemoglobin_min[i] = ifelse(is.infinite(min(subset_1y$Result.Value)), NA, min(subset_1y$Result.Value))
  outpx$Haemoglobin_sd[i] <- sd(subset_1y$Result.Value)
  outpx$Haemoglobin_grad[i] <- tryCatch(rlm(Result.Value~Specimen.Received.Date, data=subset_1y)$coefficients['Specimen.Received.Date'], error=function(err) NA)
  outpx$Haemoglobin_lat[i] <- (subset_1y %>% arrange(desc(Specimen.Received.Date)))$Result.Value[1]
  
  cat("Haemoglobin current status:", i, '/', nrow(outpx), '\n')
}
save.image("4.RData")

## Platelet Count ##
for (i in 1:nrow(outpx)) {
  #1 year look back
  lab1 <- lab[lab$Patient.ID==outpx$`Patient ID`[i],]
  subset_1y <- subset(lab1, (Specimen.Received.Date>=outpx$`Visit Date`[i]-365 & Specimen.Received.Date<= (outpx$`Visit Date`[i]-1)) &
                        lab1$Lab.Resulted.Order.Test.Description=='Platelet Count')
  outpx$Platelet_med[i] <- median(subset_1y$Result.Value)
  outpx$Platelet_mean[i] <- mean(subset_1y$Result.Value)
  outpx$Platelet_max[i] = ifelse(is.infinite(max(subset_1y$Result.Value)), NA, max(subset_1y$Result.Value))
  outpx$Platelet_min[i] = ifelse(is.infinite(min(subset_1y$Result.Value)), NA, min(subset_1y$Result.Value))
  outpx$Platelet_sd[i] <- sd(subset_1y$Result.Value)
  outpx$Platelet_grad[i] <- tryCatch(rlm(Result.Value~Specimen.Received.Date, data=subset_1y)$coefficients['Specimen.Received.Date'], error=function(err) NA)
  outpx$Platelet_lat[i] <- (subset_1y %>% arrange(desc(Specimen.Received.Date)))$Result.Value[1]
  
  cat("Platelet current status:", i, '/', nrow(outpx), '\n')
}
save.image("5.RData")

## WBC ##
for (i in 1:nrow(outpx)) {
  #1 year look back
  lab1 <- lab[lab$Patient.ID==outpx$`Patient ID`[i],]
  subset_1y <- subset(lab1, (Specimen.Received.Date>=outpx$`Visit Date`[i]-365 & Specimen.Received.Date<= (outpx$`Visit Date`[i]-1)) &
                        lab1$Lab.Resulted.Order.Test.Description=='WBC')
  outpx$WBC_med[i] <- median(subset_1y$Result.Value)
  outpx$WBC_mean[i] <- mean(subset_1y$Result.Value)
  outpx$WBC_max[i] = ifelse(is.infinite(max(subset_1y$Result.Value)), NA, max(subset_1y$Result.Value))
  outpx$WBC_min[i] = ifelse(is.infinite(min(subset_1y$Result.Value)), NA, min(subset_1y$Result.Value))
  outpx$WBC_sd[i] <- sd(subset_1y$Result.Value)
  outpx$WBC_grad[i] <- tryCatch(rlm(Result.Value~Specimen.Received.Date, data=subset_1y)$coefficients['Specimen.Received.Date'], error=function(err) NA)
  outpx$WBC_lat[i] <- (subset_1y %>% arrange(desc(Specimen.Received.Date)))$Result.Value[1]
  
  cat("WBC current status:", i, '/', nrow(outpx), '\n')
}
save.image("6.RData")

## BMI ##
bmi <- read.csv('Weight_DCI_Updated_deid.csv')
bmi <- replace(bmi, bmi=='NULL', NA)
bmi <- bmi[complete.cases(bmi),]
bmi[,c('First_Height_Checked','Weight')] <- lapply(bmi[,c('First_Height_Checked','Weight')], function(x) as.numeric(x))
bmi[,c("First_Height_Dt","Obs_Date_Weight")] <- lapply(bmi[,c("First_Height_Dt","Obs_Date_Weight")], function(x) as.Date(x,format='%d/%m/%Y'))
bmi$First_Height_Checked <- ifelse(bmi$First_Height_Checked<=5, bmi$First_Height_Checked*100,bmi$First_Height_Checked)
bmi$bmi <- bmi$Weight /((bmi$First_Height_Checked/100)^2) 
bmi <- bmi[abs(bmi$First_Height_Checked-bmi$Weight)>=5 & bmi$bmi>=10,]

for (i in 1:nrow(outpx)) {
  #1 year look back
  bmi1 <- bmi[bmi$NRIC==outpx$`Patient ID`[i],]
  subset_1y <- subset(bmi1, (Obs_Date_Weight>=outpx$`Visit Date`[i]-365 & Obs_Date_Weight<= (outpx$`Visit Date`[i]-1)))
  subset_1y <- subset_1y %>% arrange(desc(subset_1y$Obs_Date_Weight))
  
  outpx$BMI_med[i] <- median(subset_1y$bmi)
  outpx$BMI_mean[i] <- mean(subset_1y$bmi)
  outpx$BMI_max[i] = ifelse(is.infinite(max(subset_1y$bmi)), NA, max(subset_1y$bmi))
  outpx$BMI_min[i] = ifelse(is.infinite(min(subset_1y$bmi)), NA, min(subset_1y$bmi))
  outpx$BMI_sd[i] <- sd(subset_1y$bmi)
  outpx$BMI_grad[i] <- tryCatch(rlm(bmi~Obs_Date_Weight, data=subset_1y)$coefficients['Obs_Date_Weight'], error=function(err) NA)
  outpx$BMI_lat[i] <- subset_1y$bmi[1]
  
  cat("BMI current status:", i, '/', nrow(outpx), '\n')
}
save.image("7.RData")

## Neut-Lymph Ratio - Latest within 1 yr ##
for (i in 1:nrow(outpx)) {
  #1 year look back
  NL_ratio1 <- NL_ratio[NL_ratio$Patient.ID==outpx$`Patient ID`[i],]
  subset_1y <- subset(NL_ratio1, (Specimen.Received.Date>=outpx$`Visit Date`[i]-365 & Specimen.Received.Date<= (outpx$`Visit Date`[i]-1)))
  subset_1y <- arrange(subset_1y, desc(Specimen.Received.Date))

  outpx$NL_ratio_med[i] <- tryCatch(median(subset_1y$NL_ratio), error=function(err) NA)
  outpx$NL_ratio_mean[i] <- tryCatch(mean(subset_1y$NL_ratio), error=function(err) NA)
  outpx$NL_ratio_max[i] <- ifelse(is.infinite(max(subset_1y$NL_ratio)), NA, max(subset_1y$NL_ratio))
  outpx$NL_ratio_min[i] <- ifelse(is.infinite(min(subset_1y$NL_ratio)), NA, min(subset_1y$NL_ratio))
  outpx$NL_ratio_sd[i] <- tryCatch(sd(subset_1y$NL_ratio), error=function(err) NA)
  outpx$NL_ratio_grad[i] <- tryCatch(rlm(NL_ratio~Specimen.Received.Date, data=subset_1y)$coefficients['Specimen.Received.Date'], error=function(err) NA)
  outpx$NL_ratio_lat[i] <- tryCatch(subset_1y$NL_ratio[1], error=function(err) NA)
  
  cat("NL ratio current status:", i, '/', nrow(outpx), '\n')
}
save.image("8.RData")

## Albumin-Gobulin Ratio - Latest within 1 yr ##
for (i in 1:nrow(outpx)) {
  #1 year look back
  AG_ratio1 <- AG_ratio[AG_ratio$Patient.ID==outpx$`Patient ID`[i],]
  subset_1y <- subset(AG_ratio1, (Specimen.Received.Date>=outpx$`Visit Date`[i]-365 & Specimen.Received.Date<= (outpx$`Visit Date`[i]-1)))
  subset_1y <- arrange(subset_1y, desc(Specimen.Received.Date))
  
  outpx$AG_ratio_med[i] <- tryCatch(median(subset_1y$AG_ratio), error=function(err) NA)
  outpx$AG_ratio_mean[i] <- tryCatch(mean(subset_1y$AG_ratio), error=function(err) NA)
  outpx$AG_ratio_max[i] <- ifelse(is.infinite(max(subset_1y$AG_ratio)), NA, max(subset_1y$AG_ratio))
  outpx$AG_ratio_min[i] <- ifelse(is.infinite(min(subset_1y$AG_ratio)), NA, min(subset_1y$AG_ratio))
  outpx$AG_ratio_sd[i] <- tryCatch(sd(subset_1y$AG_ratio), error=function(err) NA)
  outpx$AG_ratio_grad[i] <- tryCatch(rlm(AG_ratio~Specimen.Received.Date, data=subset_1y)$coefficients['Specimen.Received.Date'], error=function(err) NA)
  outpx$AG_ratio_lat[i] <- tryCatch(subset_1y$AG_ratio[1], error=function(err) NA)
  
  cat("AG ratio current status:", i, '/', nrow(outpx), '\n')
}
save.image("9.RData")

write.csv(outpx, 'outpx_lab-bmi_full.csv', row.names = FALSE)
#outpx <- read.csv('outpx_lab-bmi_full.csv', check.names = FALSE)
#outpx$`Visit Date` <- as.Date(outpx$`Visit Date`)

#-----------#
#---Drugs---#
#-----------#
px_diag_date <- data %>% dplyr::select(NRIC, Diag_Date)
px_diag_date$Diag_Date <- as.Date(px_diag_date$Diag_Date)

drugs <- read.csv('Cancer Drugs.csv', check.names = FALSE, fileEncoding="UTF-8-BOM")
drugs$Drug_Dispensed_Date_From <- as.Date(drugs$Drug_Dispensed_Date_From)
drugs$Combined_ATC_Codes <- ifelse(drugs$Trial_Drugs=='YES', 'TRIAL', drugs$Combined_ATC_Codes) #All trial drugs are grouped
drugs <- drugs %>% 
  select(Patient_ID, Drug_Dispensed_Date_From, Combined_ATC_Codes, L01A, L01B, L01C, L01D, L01E, L01F, L01X, L02A, L02B, Trial_Drugs) %>% 
  group_by(Patient_ID) %>% 
  arrange(Drug_Dispensed_Date_From, .by_group = TRUE) %>% 
  ungroup()
drugs <- inner_join(drugs, px_diag_date, by=c('Patient_ID'='NRIC'))
drugs <- select(drugs, Patient_ID, Diag_Date, everything())
drugs[,5:ncol(drugs)] <- ifelse(drugs[,5:ncol(drugs)]=='YES',1,0)

for (i in 1:nrow(outpx)) {
  #Drugs from diagnosis date to visit date
  drugs1 <- drugs[drugs$Patient_ID==outpx$`Patient ID`[i],]
  subset_drugs <- subset(drugs1, (Drug_Dispensed_Date_From>=(drugs1$Diag_Date[1]-14) & Drug_Dispensed_Date_From<= (outpx$`Visit Date`[i]-1)))
  subset_drugs <- arrange(subset_drugs, desc(Drug_Dispensed_Date_From))
  subset_drugs_sum <- subset_drugs[,5:ncol(drugs)] %>% summarise(across(everything(), sum))
  drug_outpx <- tryCatch(cbind(`Patient ID`=outpx$`Patient ID`[i],`Visit Date`=outpx$`Visit Date`[i], 'unique_drug_count'=length(unique(subset_drugs$Combined_ATC_Codes)), 
                               subset_drugs_sum), error=function(err) NA)
  
  if (i==1) {
    drug_df <- drug_outpx
  }
  else {
    drug_df <- tryCatch(rbind(drug_df, drug_outpx), error=function(err) NA)
  }

  cat("Drugs current status:", i, '/', nrow(outpx), '\n')
}       

drug_df[,4:ncol(drug_df)] <- ifelse(drug_df[,4:ncol(drug_df)]>=1, 1,0)
drug_df <- drug_df[,-c(1:2)]
outpx <- cbind(outpx, drug_df)

write.csv(outpx, 'outpx_lab-bmi-drugs.csv', row.names = FALSE)
#outpx <- read.csv('outpx_lab-bmi-drugs.csv', check.names = FALSE)
#outpx$`Visit Date` <- as.Date(outpx$`Visit Date`)

#---------------#
#---Diagnosis---#
#---------------#
diagnosis <- read_excel('PROTECH_Diagnosis_Complete__deid.xlsx')
diagnosis$Diagnosis_Date <- as.Date(diagnosis$Diagnosis_Date)
diagnosis <- diagnosis[diagnosis$Patient_ID %in% data$NRIC,] #Filter to just those in the cohort

diagnosis <- diagnosis %>% 
  filter(!grepl("^[[:digit:]]+", Diagnosis_Code_ICD10)) %>% #Removes Snowmed codes
  group_by(Patient_ID, Diagnosis_Date, Diagnosis_Code_ICD10) %>% 
  slice(1) %>% 
  arrange(Diagnosis_Date)

diagnosis$Patient_ID_date <- paste0(diagnosis$Patient_ID, ' ', diagnosis$Diagnosis_Date)

map_key <- data.frame('Patient_ID_date'=unique(diagnosis$Patient_ID_date))
map_key$Index <- 1:nrow(map_key)

diagnosis <- left_join(diagnosis, map_key, by = 'Patient_ID_date') %>%
  select(Patient_ID, Index, Patient_ID_date, everything())
elixhauser <- comorbidity(x = diagnosis, id = "Index", code = "Diagnosis_Code_ICD10", map = "elixhauser_icd10_quan", assign0 = FALSE)
elixhauser <- as.data.frame(elixhauser)

diag_elix <- left_join(map_key, elixhauser, by=c("Index"="Index"))
diag_elix <- diag_elix[, !(names(diag_elix) %in% 'Index')]
colnames(diag_elix)[2:32] <- paste0('elix_', colnames(diag_elix)[2:32])

diag_elix$Patient_ID <- sapply(strsplit(diag_elix$Patient_ID_date," "), `[`, 1)
diag_elix$Diagnosis_Date <- as.Date(sapply(strsplit(diag_elix$Patient_ID_date," "), `[`, 2))
diag_elix <- diag_elix %>% select(Patient_ID, Diagnosis_Date, everything()) %>% select(-Patient_ID_date)

for (i in 1:nrow(outpx)) {
  #All diagnosis from before visit date
  diag_elix1 <- diag_elix[diag_elix$Patient_ID==outpx$`Patient ID`[i],]
  subset_1y <- subset(diag_elix1, Diagnosis_Date<= (outpx$`Visit Date`[i]-1))
  subset_1y <- arrange(subset_1y, desc(Diagnosis_Date))
  subset_1y_sum <- subset_1y[,3:ncol(subset_1y)] %>% summarise(across(everything(), sum))
  diag_outpx <- tryCatch(cbind(`Patient ID`=outpx$`Patient ID`[i],`Visit Date`=outpx$`Visit Date`[i], subset_1y_sum), error=function(err) NA)
  
  if (i==1) {
    diag_df <- diag_outpx
  }
  else {
    diag_df <- tryCatch(rbind(diag_df, diag_outpx), error=function(err) NA)
  }
  
  cat("Diagnosis current status:", i, '/', nrow(outpx), '\n')
}       

diag_df[,3:ncol(diag_df)] <- ifelse(diag_df[,3:ncol(diag_df)]>=1, 1,0)
diag_df <- diag_df[,-c(1:2)]
outpx <- cbind(outpx, diag_df)

write.csv(outpx, 'outpx_lab-bmi-drugs_diag.csv', row.names = FALSE)

#-----------#
#---Visits--#
#-----------#
#Outpx visits / ED visits / Inpx visits - 1year/6m/3m lookback

#Outpx visits
for (i in 1:nrow(outpx)) {
  #All diagnosis from before visit date
  outpx1 <- outpx[outpx$`Patient ID`==outpx$`Patient ID`[i],]
  subset_1y <- subset(outpx1, `Visit Date`<=(outpx$`Visit Date`[i]-1) & `Visit Date`>=outpx$`Visit Date`[i]-365)
  subset_6m <- subset(outpx1, `Visit Date`<=(outpx$`Visit Date`[i]-1) & `Visit Date`>=outpx$`Visit Date`[i]-180)
  subset_3m <- subset(outpx1, `Visit Date`<=(outpx$`Visit Date`[i]-1) & `Visit Date`>=outpx$`Visit Date`[i]-90)
  
  outpx$outpx_visits_1y[i] <- nrow(subset_1y)
  outpx$outpx_visits_6m[i] <- nrow(subset_6m)
  outpx$outpx_visits_3m[i] <- nrow(subset_3m)

  cat("Outpx current status:", i, '/', nrow(outpx), '\n')
}       


#ED visits
ED_data <- read.csv('PROTECH_A&E Visits_Complete_deid.csv', check.names = FALSE)
ED_data$`Admit/Visit Date` <- as.Date(sapply(strsplit(ED_data$`Admit/Visit Date`," "), `[`, 1))
ED_data <- ED_data %>% 
  select(`Patient ID`, `Admit/Visit Date`) %>% 
  group_by(`Patient ID`, `Admit/Visit Date`) %>%
  slice(1) %>% 
  arrange(`Admit/Visit Date` , .by_group = TRUE)

for (i in 1:nrow(outpx)) {
  #All diagnosis from before visit date
  ED_data1 <- ED_data[ED_data$`Patient ID`==outpx$`Patient ID`[i],]
  subset_1y <- subset(ED_data1, `Admit/Visit Date`<=(outpx$`Visit Date`[i]-1) & `Admit/Visit Date`>=outpx$`Visit Date`[i]-365)
  subset_6m <- subset(ED_data1, `Admit/Visit Date`<=(outpx$`Visit Date`[i]-1) & `Admit/Visit Date`>=outpx$`Visit Date`[i]-180)
  subset_3m <- subset(ED_data1, `Admit/Visit Date`<=(outpx$`Visit Date`[i]-1) & `Admit/Visit Date`>=outpx$`Visit Date`[i]-90)
  
  outpx$ED_1y[i] <- nrow(subset_1y)
  outpx$ED_6m[i] <- nrow(subset_6m)
  outpx$ED_3m[i] <- nrow(subset_3m)
  
  cat("ED current status:", i, '/', nrow(outpx), '\n')
}       


#Inpx
inpx_data <- read.csv('PROTECH_Inpatient_Complete_deid.csv', check.names = FALSE)
inpx_data$`Admit/Visit Date` <- as.Date(sapply(strsplit(inpx_data$`Admit/Visit Date`," "), `[`, 1))
inpx_data <- inpx_data %>% 
  select(`Patient ID`, `Admit/Visit Date`) %>% 
  group_by(`Patient ID`, `Admit/Visit Date`) %>%
  slice(1) %>% 
  arrange(`Admit/Visit Date` , .by_group = TRUE)

for (i in 1:nrow(outpx)) {
  #All diagnosis from before visit date
  inpx_data1 <- inpx_data[inpx_data$`Patient ID`==outpx$`Patient ID`[i],]
  subset_1y <- subset(inpx_data1, `Admit/Visit Date`<=(outpx$`Visit Date`[i]-1) & `Admit/Visit Date`>=outpx$`Visit Date`[i]-365)
  subset_6m <- subset(inpx_data1, `Admit/Visit Date`<=(outpx$`Visit Date`[i]-1) & `Admit/Visit Date`>=outpx$`Visit Date`[i]-180)
  subset_3m <- subset(inpx_data1, `Admit/Visit Date`<=(outpx$`Visit Date`[i]-1) & `Admit/Visit Date`>=outpx$`Visit Date`[i]-90)
  
  outpx$inpx_1y[i] <- nrow(subset_1y)
  outpx$inpx_6m[i] <- nrow(subset_6m)
  outpx$inpx_3m[i] <- nrow(subset_3m)
  
  cat("Inpx current status:", i, '/', nrow(outpx), '\n')
}  

#------------------------#
#---Merge & Preprocess---#
#------------------------#
# load("C:/Users/AI-Ball/Desktop/PROTECH/PROTECH preprocessing (1yr death) - Before Merge.RData")
data <- left_join(data, outpx, by=c('NRIC'='Patient ID')) 
data[,c('Visit Date','Diag_Date','DeathDate','DOB')] <- lapply(data[,c('Visit Date','Diag_Date','DeathDate','DOB')] , date)
data$Age_At_visit <- floor(time_length(difftime(as.Date(data$`Visit Date`), as.Date(data$DOB)),"years"))

data <- data %>% 
  mutate(day_diff=`Visit Date`-Diag_Date) %>% 
  select(NRIC, Diag_Date, `Visit Date`, day_diff, everything()) %>% 
  filter(!is.na(Diag_Date)) %>%  #Remove NA diagnosis dates
  filter(`Visit Date`-Diag_Date >= 0) %>% #No visits before diagnosis dates
  group_by(NRIC, `Visit Date`) %>% 
  arrange(`Visit Date`) %>%
  slice(1) %>%  #Remove duplicate visit dates
  ungroup()
data <- data %>% group_by(`NRIC`) %>% slice(2:n()) #Second visit onwards

data$death1yr <- if_else(data$DeathDate-data$`Visit Date`<=365, 1, 0, 0)
table(data$death1yr)

#Down sample population
count_outpx <- data %>% 
  group_by(NRIC) %>%
  summarise(count = n()) %>% 
  ungroup() 
data <- left_join(data,count_outpx,by='NRIC')

#For those above median visits, take the last day of each month
median_outpx <- median((data %>% count(NRIC))$n)
data$Month_Yr <- format(as.Date(data$`Visit Date`), "%Y-%m")
data_a <- data %>%
  filter(count>median_outpx) %>% 
  group_by(NRIC,Month_Yr) %>% 
  slice(n()) %>% 
  ungroup()
data_b <- data %>% 
  filter(count<=median_outpx) %>% 
  ungroup()
data <- rbind(data_a, data_b)

#Filter the Dates
data$DeathDate <- as.character(data$DeathDate)
data$DeathDate <- ifelse(substr(data$DeathDate, 1,4)=="2022", NA, data$DeathDate)
data$DeathDate <- as.Date(data$DeathDate)
data <- data[data$`Visit Date`<="2020-12-31",] #censor date for collecting mortality was 31/12/2021, last prediction date is -1 year

length(unique(data$NRIC)) #Total no. of patients

#Save before removing variables
#save.image("C:/Users/AI-Ball/Desktop/PROTECH/PROTECH preprocessing (1yr death) descriptive.RData")
#load("C:/Users/AI-Ball/Desktop/PROTECH/PROTECH preprocessing (1yr death) descriptive.RData")

#Remove unnecessary variables (1st round)
data <- data %>% select(-c("Diag_Date","Visit Date","DeathDate","Cons_Date","FinalStage", "Recoded_ICD10",
                           "DOB","Age_At_Dx","A (Stage III/IV during study period)","B","Month_Yr", "count"))

#Save Data
save.image("C:/Users/AI-Ball/Desktop/PROTECH/PROTECH preprocessing (1yr death) descriptive.RData")
load("C:/Users/AI-Ball/Desktop/PROTECH/PROTECH preprocessing (1yr death).RData")
write.csv(data, 'PROTECH_preprocess.csv', row.names = FALSE)

