##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                                                                            --
##------------------ MATCHING PROCEDURE FOR PAIRED CONTROLS---------------------
##                          VERSION: JUN, 16TH 2023                           --
##                                                                            ~~
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



library("ccoptimalmatch")
load("base.RData") # can be any dataset with cases and controls


matching_variables<-c("sex", 
                      "I_age_1", 
                      "I_edu_VE",
                      "I_duur_klachten", 
                      "V_MMSE")



cognitive_unpaired<-base %>%  
  filter(diagnoseOms=="AD" & ( CSF_A=="A-"|PT_result_amyloid=="negatief"|PT_result_amyloid=="Negatief"))%>%
  select(I_ID, matching_variables)
cases$dx<-"AD" #you can skip this, I have more than two groups so I used to have them cleared
cases$case_control<-"case" ##### We going to need this variable  DONT ERASE THIS PART

cognitive_unpaired<-base %>%  
  filter(diagnoseOms=="Subjectieve klachten" & ( CSF_A=="A-"|PT_result_amyloid=="negatief"|PT_result_amyloid=="Negatief"))%>%
  select(I_ID, matching_variables)
cognitive_unpaired$dx<-"CU"
cognitive_unpaired$case_control<-"control"

rm(base)

# Step 1: Exact Matching on several variables

##We start by defining the subsets. In order to define the subsets, we have to decide wich variable we want exactly matched between groups (not posibility to have a range), useful for categorical variables
##we filter by the “cases”, take the distinct combination of exact variables

create_subset <- cases %>% 
  arrange(sex)%>%
  distinct(sex, #it can be a list of variables
           .keep_all = TRUE) %>%
  mutate(subset = 1:n()) %>%
  select( sex ,   subset)

##We merge the data that contains the “subset” variable with the data that contains the cases only
cases_subsetting <- cases %>% 
  full_join(create_subset, by = c("sex" ))

## We merge the data that contains the “subset” variable with the data that contains the controls only:
cognitive_unpaired_subsetting <- cognitive_unpaired %>% 
  full_join(create_subset, by = c("sex" )) %>% drop_na(subset,dx)


##Finally we bind the cases and the controls, which will have now the new variable “subset”:
not_processed <- rbind(cases_subsetting,cognitive_unpaired_subsetting)

# Step 2: Create deaviations variables and select the range of variables
##this step help us to pick a set of variables to be paired within a range (and the range)

bdd_controls <- not_processed[not_processed$dx=="CU",]
bdd_controls$cluster_case <- 0
bdd_cases <- not_processed[not_processed$dx=="AD",]
bdd_cases$cluster_case <- paste("case",1:nrow(bdd_cases),sep = "_")

not_processed <- rbind(bdd_cases,bdd_controls)
bdd_cases <- not_processed[not_processed$dx=="AD",]
bdd_control <- not_processed[not_processed$dx=="CU",]

bdd_temp_list <- list()
bdd_temp <- data.frame()
list_p <- unique(bdd_cases$cluster_case)

#This loops evaluates each of the controls (several per subject) with the subject 
#it is paired with, estimates the differences between the variables we want to pair 
#by rank, checks the rank on all of them and decides if it is acceptable to include or not.

for (i in 1:length(list_p)) {
  temp <- bdd_cases[bdd_cases$cluster_case == list_p[i], ]
  subset_identified <- temp$subset
  temp0 <- bdd_control[bdd_control$subset == temp$subset, ]
  temp_final <- rbind(temp, temp0)
  temp_final$cluster_case <- list_p[i]
  temp_final <- temp_final %>%
    group_by(cluster_case) %>%
    mutate(age_diff = abs(I_age_1 - I_age_1[dx == "bvAD"]),
           edu_diff = abs(I_edu_VE - I_edu_VE[dx == "bvAD"]))#can be as many as we want
  temp_final$matching_variable_criteria <- ifelse(temp_final$age_diff <= 5&temp_final$edu_diff<=2, "accept", "delete") #this is the condition of range that we want
  temp_final <- temp_final[temp_final$matching_variable_criteria == "accept", ]
  bdd_temp_list[[i]] <- temp_final
}

bdd_temp <- do.call(rbind, bdd_temp_list)

bdd_temp<-bdd_temp %>% drop_na()

# Step 3: Create the variables “total controls per case” and “frequency of controls”
bdd_temp = bdd_temp %>% group_by(cluster_case) %>% mutate(total_control_per_case = n()-1)
bdd_temp$case_ind <- ifelse(bdd_temp$dx=="bvAD",1,0)
bdd_temp <- subset(bdd_temp, select=c(cluster_case, I_ID, dx, case_ind,
                                      age_diff, edu_diff,sex, I_age_1,I_edu_VE, case_control,total_control_per_case))


bdd_temp = bdd_temp %>% group_by(I_ID) %>% mutate(freq_of_controls = n())

#this step is IMPORTANT, the order define the hierarchy of matching
bdd_temp<-bdd_temp[order(bdd_temp$cluster_case,bdd_temp$case_control, bdd_temp$age_diff, bdd_temp$edu_diff, bdd_temp$freq_of_controls),] 

test<-bdd_temp


#This function goes through the hierarchy of variables with ranks and choosing the closest ones
final_data <- optimal_matching(bdd_temp, 
                               n_con=5,#the max num of paired controls that we want to retrieve (if its possible) 
                               cluster_case, I_ID, 
                               total_control_per_case, case_control, with_replacement = FALSE)


final_data <- final_data %>% group_by(cluster_case) %>% mutate(total_control_matched = n()-1)

table(final_data$case_control,final_data$total_control_matched)


