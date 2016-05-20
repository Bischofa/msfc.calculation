# Read in data ------------------------------------------------------------

raw_data = read.csv("data/epic_msfc_20160519.csv", stringsAsFactors = FALSE, sep = ";")

# Fill in missing scores that were due to MS ------------------------------

# For missing PASAT scores due to MS, fill in 0 (follow MSFC manual rules)
index = which(is.na(raw_data$PASAT3_Score) & raw_data$PASAT3_Incomplete_Value_Analytical_Status == "Missing data due to MS")
raw_data$TotalScore_missing_filled_in = raw_data$PASAT3_Score
raw_data$TotalScore_missing_filled_in[index] = 0

# For missing 9HPT scores due to MS, fill in 777 (follow MSFC manual rules)
index = which(is.na(raw_data$X9HPT_Dominant_Hand_Trial_1) & raw_data$X9HPT_Dominant_Hand_Trial_1_Incomplete_Value_Analytical_Status == "Missing data due to MS")
raw_data$X9HPT_Dominant_Trial1_missing_filled_in = raw_data$X9HPT_Dominant_Hand_Trial_1
raw_data$X9HPT_Dominant_Trial1_missing_filled_in[index] = 777

index = which(is.na(raw_data$X9HPT_Dominant_Hand_Trial_2) & raw_data$X9HPT_Dominant_Hand_Trial_2_Incomplete_Value_Analytical_Status == "Missing data due to MS")
raw_data$X9HPT_Dominant_Trial2_missing_filled_in = raw_data$X9HPT_Dominant_Hand_Trial_2
raw_data$X9HPT_Dominant_Trial2_missing_filled_in[index] = 777

index = which(is.na(raw_data$X9HPT_Nondominant_Hand_Trial_1) & raw_data$X9HPT_Nondominant_Hand_Trial_1_Incomplete_Value_Analytical_Status == "Missing data due to MS")
raw_data$X9HPT_NonDominant_Trial1_missing_filled_in = raw_data$X9HPT_Nondominant_Hand_Trial_1
raw_data$X9HPT_NonDominant_Trial1_missing_filled_in[index] = 777

index = which(is.na(raw_data$X9HPT_Nondominant_Hand_Trial_2) & raw_data$X9HPT_Nondominant_Hand_Trial_2_Incomplete_Value_Analytical_Status == "Missing data due to MS")
raw_data$X9HPT_NonDominant_Trial2_missing_filled_in = raw_data$X9HPT_Nondominant_Hand_Trial_2
raw_data$X9HPT_NonDominant_Trial2_missing_filled_in[index] = 777

# For missing T25W scores due to MS, fill in 165.7948 (follow MSFC manual rules)
# 165.7948 was obtained from the Task Force Database: T25W mean (9.5353), T25W sd (11.4058), and largest recorded T25W z-score (13.7)
t25wFill = 13.7*11.4058 + 9.5353 # = 165.7948
index = which(is.na(raw_data$X25TW_Trial_1) & raw_data$X25TW_Trial_1_Incomplete_Value_Analytical_Status == "Missing data due to MS")
raw_data$Trial1_missing_filled_in = raw_data$X25TW_Trial_1
raw_data$Trial1_missing_filled_in[index] = t25wFill

index = which(is.na(raw_data$X25TW_Trial_2) & raw_data$X25TW_Trial_2_Incomplete_Value_Analytical_Status == "Missing data due to MS")
raw_data$Trial2_missing_filled_in = raw_data$X25TW_Trial_2
raw_data$Trial2_missing_filled_in[index] = t25wFill

# Calculate mean values ---------------------------------------------------

# Add mean T25W (mean of two walk trials)
raw_data$MSFC_Walk_Trial_Mean = apply(cbind(raw_data$Trial1_missing_filled_in, raw_data$Trial2_missing_filled_in), 1, mean, na.rm = T)
raw_data$MSFC_Walk_Trial_Mean[is.nan(raw_data$MSFC_Walk_Trial_Mean)] = NA

# Add mean 9HPT: dominant hand (mean of two trials for dominant hand)
raw_data$MSFC_9HPT_DomHand_Trial_Mean = apply(cbind(raw_data$X9HPT_Dominant_Trial1_missing_filled_in, raw_data$X9HPT_Dominant_Trial2_missing_filled_in), 1, mean, na.rm = T)
raw_data$MSFC_9HPT_DomHand_Trial_Mean[is.nan(raw_data$MSFC_9HPT_DomHand_Trial_Mean)] = NA

# Add mean 9HPT: non-dominant hand (mean of two trials for non-dominant hand)
raw_data$MSFC_9HPT_NonDomHand_Trial_Mean = apply(cbind(raw_data$X9HPT_NonDominant_Trial1_missing_filled_in, raw_data$X9HPT_NonDominant_Trial2_missing_filled_in), 1, mean, na.rm = T)
raw_data$MSFC_9HPT_NonDomHand_Trial_Mean[is.nan(raw_data$MSFC_9HPT_NonDomHand_Trial_Mean)] = NA

# Create MSFC function ----------------------------------------------------

calculate_msfc_scores = function(nhpt_dominant_trial1_trial2_average, 
                                 nhpt_non_dominant_trial1_trial2_average,
                                 t25w_trial1_trial1_average, 
                                 pasat_score,
                                 baseline_nhpt_dominant_trial1_trial2_average, 
                                 baseline_nhpt_non_dominant_trial1_trial2_average,
                                 baseline_t25w_trial1_trial1_average, 
                                 baseline_pasat_score,
                                 between_studies = FALSE){
  
  
  nhpt_dominant_non_dominant_average = apply(data.frame(1/nhpt_dominant_trial1_trial2_average, 1/nhpt_non_dominant_trial1_trial2_average), 1, mean)
  baseline_nhpt_dominant_non_dominant_average = apply(data.frame(1/baseline_nhpt_dominant_trial1_trial2_average, 1/baseline_nhpt_non_dominant_trial1_trial2_average), 1, mean)
  
  if(between_studies){
    MSFC_9HPT_Zscore = (nhpt_dominant_non_dominant_average - .0439)/.0101 
    MSFC_Walk_Zscore = ((t25w_trial1_trial1_average - 9.5353)/11.4058) * -1
    MSFC_PASAT_ZScore = (pasat_score - 45.0311)/12.0771
  } else{
    
    # calculate mean of all baseline scores
    mean_baseline_nhpt = mean(baseline_nhpt_dominant_non_dominant_average, na.rm = TRUE)
    mean_baseline_t25w = mean(baseline_t25w_trial1_trial1_average, na.rm = TRUE)
    mean_baseline_pasat = mean(baseline_pasat_score, na.rm = TRUE)
    
    # calculate sd of all baseline scores
    sd_baseline_nhpt = sd(baseline_nhpt_dominant_non_dominant_average, na.rm = TRUE)
    sd_baseline_t25w = sd(baseline_t25w_trial1_trial1_average, na.rm = TRUE)
    sd_baseline_pasat = sd(baseline_pasat_score, na.rm = TRUE)
    
    # calculate z scores
    MSFC_9HPT_Zscore = (nhpt_dominant_non_dominant_average - mean_baseline_nhpt)/sd_baseline_nhpt
    MSFC_Walk_Zscore = ((t25w_trial1_trial1_average - mean_baseline_t25w)/sd_baseline_t25w) * -1
    MSFC_PASAT_ZScore = (pasat_score - mean_baseline_pasat)/sd_baseline_pasat
    
  }
  
  MSFC_Score = (MSFC_9HPT_Zscore + MSFC_Walk_Zscore + MSFC_PASAT_ZScore)/3
  return(data.frame(MSFC_9HPT_Zscore, MSFC_Walk_Zscore, MSFC_PASAT_ZScore, MSFC_Score))
  
}

# Calculate MSFC scores ---------------------------------------------------
baseline_nhpt_dominant_trial1_trial2_average = raw_data$MSFC_9HPT_DomHand_Trial_Mean[raw_data$VisitType == "Baseline"] 
baseline_nhpt_non_dominant_trial1_trial2_average = raw_data$MSFC_9HPT_NonDomHand_Trial_Mean[raw_data$VisitType == "Baseline"]
baseline_t25w_trial1_trial1_average = raw_data$MSFC_Walk_Trial_Mean[raw_data$VisitType == "Baseline"]
baseline_pasat_score = raw_data$TotalScore_missing_filled_in[raw_data$VisitType == "Baseline"]

msfc_score_info = calculate_msfc_scores(raw_data$MSFC_9HPT_DomHand_Trial_Mean, 
                                        raw_data$MSFC_9HPT_NonDomHand_Trial_Mean, 
                                        raw_data$MSFC_Walk_Trial_Mean, 
                                        raw_data$TotalScore_missing_filled_in, 
                                        baseline_nhpt_dominant_trial1_trial2_average,
                                        baseline_nhpt_non_dominant_trial1_trial2_average,
                                        baseline_t25w_trial1_trial1_average, 
                                        baseline_pasat_score, 
                                        between_studies = FALSE)

calculated_msfc_data = cbind(raw_data, msfc_score_info)

write.csv(calculated_msfc_data, "output/calculated_msfc_data.csv", row.names = FALSE)

