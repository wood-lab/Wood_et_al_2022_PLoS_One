#---
title: "Permutation tests"
author: "Chelsea Wood"
date: "updated 21 October 2021"
#---


# Realizing that we had very low replication in one cell of our color * sex interaction 
# (n = 4 blue males), we sought to test whether the patterns revealed by the analysis above 
# could have arisen due to random chance. 

### First iteration of the permutation test - color randomly assigned only across males

nreps<-500
rand_estimate_interaction<-vector()
rand_estimate_color<-vector()
rand_estimate_sex<-vector()
rand_z_interaction<-vector()
rand_z_color<-vector()
rand_z_sex<-vector()
actual_z_interaction<--4.1193193
actual_z_color<-0.5945273
actual_z_sex<-5.4636742

males<-melted_data %>%
  filter(sex == "M")
females<-melted_data%>%
  filter(sex == "F")

females$rand_color<-females$color
females$rand_color<-gsub("BL","rand_blue",females$rand_color)
females$rand_color<-gsub("BR","rand_brown",females$rand_color)

for(i in 1:nreps){
  
  # First, randomly select your "blue" fish and make the remainder of fish "brown"
  
  male_sample_IDs<-unique(males$sample_ID)
  male_sample_IDs<-as.data.frame(male_sample_IDs)
  rand_blue_sample_IDs<-male_sample_IDs[sample(nrow(male_sample_IDs),4),]
  
  rand_blue_males<- males %>%
    filter(sample_ID == rand_blue_sample_IDs[1] | sample_ID == rand_blue_sample_IDs[2] |
             sample_ID == rand_blue_sample_IDs[3] | sample_ID == rand_blue_sample_IDs[4])
  
  rand_brown_males<- males %>%
    filter(sample_ID != rand_blue_sample_IDs[1] & sample_ID != rand_blue_sample_IDs[2] &
             sample_ID != rand_blue_sample_IDs[3] & sample_ID != rand_blue_sample_IDs[4])
  
  # Then, concatenate your "blue" and "brown" males and make sure that they are labelled in your new dataset.
  
  rand_data_males<-rbind(rand_blue_males,rand_brown_males)
  rand_blue_labels<-rep("rand_blue",length(rand_blue_males$sample_ID))
  rand_brown_labels<-rep("rand_brown",length(rand_brown_males$sample_ID))
  rand_color<-c(rand_blue_labels,rand_brown_labels)
  rand_data_males$rand_color<-rand_color
  
  # Put the males and females back together
  
  rand_data<-rbind(rand_data_males,females)
  
  # Then, run your analysis and store the result for the rand_color*sex interaction
  
  cat("starting model",i)
  nbmodel_rand<-glmer.nb(as.numeric(psite_count_sc)~depth_ft_sc+sex*rand_color+(1|sampling_region/sampling_location/sample_ID)+
                           (1 + sex*rand_color|psite_spp)+offset(log(TL_cm_sc)),data=rand_data,family="nbinom")
  
  fixed_effects<-fixef(nbmodel_rand)
  rand_estimate_interaction[i]<-fixed_effects[5]
  rand_estimate_color[i]<-fixed_effects[4]
  rand_estimate_sex[i]<-fixed_effects[3]
  rand_z_interaction[i]<-coef(summary(nbmodel_rand))[5,3]
  rand_z_color[i]<-coef(summary(nbmodel_rand))[4,3]
  rand_z_sex[i]<-coef(summary(nbmodel_rand))[3,3]
  cat("ending model",i)
  
}


hist(rand_z_interaction)

output<-cbind(rand_estimate_interaction,rand_estimate_color,rand_estimate_sex,rand_z_interaction,rand_z_color,rand_z_sex)
write.csv(output,"output_2020.08.13.csv")

first100<-read.csv("output_2020.07.17.csv",header=T,sep=",")
final408<-read.csv("output_2020.08.13.csv",header=T,sep=",")
length(final408$X)

final400<-final408[-c(401:408),]
length(final400$X)

final_500<-rbind(first100,final400)
length(final_500$X)

write.csv(final_500,"final_500_simulated_values_only_males_randomized.csv")

# Then, compare your actual value for the color*sex interaction to the distribution of randomly-created interaction values

prob_z_interaction<-length(final_500$rand_z_interaction[abs(final_500$rand_z_interaction)>=abs(actual_z_interaction)])/length(final_500$rand_estimate_interaction)

cat(" The probability value for the interaction is ",prob_z_interaction, "\n")



### Second iteration of the permutation test - color randomly assigned across both males and females

nreps<-500
rand_estimate_interaction<-vector()
rand_estimate_color<-vector()
rand_estimate_sex<-vector()
rand_z_interaction<-vector()
rand_z_color<-vector()
rand_z_sex<-vector()
actual_z_interaction<--4.1193193
actual_z_color<-0.5945273
actual_z_sex<-5.4636742

males<-melted_data %>%
  filter(sex == "M")
females<-melted_data%>%
  filter(sex == "F")


for(i in 1:nreps){
  
  # First, randomly select your "blue" fish and make the remainder of fish "brown"
  
  male_sample_IDs<-unique(males$sample_ID)
  male_sample_IDs<-as.data.frame(male_sample_IDs)
  rand_blue_sample_IDs<-male_sample_IDs[sample(nrow(male_sample_IDs),4),]
  
  rand_blue_males<- males %>%
    filter(sample_ID == rand_blue_sample_IDs[1] | sample_ID == rand_blue_sample_IDs[2] |
             sample_ID == rand_blue_sample_IDs[3] | sample_ID == rand_blue_sample_IDs[4])
  
  rand_brown_males<- males %>%
    filter(sample_ID != rand_blue_sample_IDs[1] & sample_ID != rand_blue_sample_IDs[2] &
             sample_ID != rand_blue_sample_IDs[3] & sample_ID != rand_blue_sample_IDs[4])
  
  # Then, concatenate your "blue" and "brown" males and make sure that they are labelled in your new dataset.
  
  rand_data_males<-rbind(rand_blue_males,rand_brown_males)
  rand_blue_labels<-rep("rand_blue",length(rand_blue_males$sample_ID))
  rand_brown_labels<-rep("rand_brown",length(rand_brown_males$sample_ID))
  rand_color<-c(rand_blue_labels,rand_brown_labels)
  rand_data_males$rand_color<-rand_color
  
  # Now do the same thing for the females, but select 11 (since that's the number from the actual dataset).
  
  female_sample_IDs<-unique(females$sample_ID)
  female_sample_IDs<-as.data.frame(female_sample_IDs)
  rand_blue_female_sample_IDs<-female_sample_IDs[sample(nrow(female_sample_IDs),11),]
  
  rand_blue_females<- females %>%
    filter(sample_ID == rand_blue_female_sample_IDs[1] | sample_ID == rand_blue_female_sample_IDs[2] |
             sample_ID == rand_blue_female_sample_IDs[3] | sample_ID == rand_blue_female_sample_IDs[4] | 
             sample_ID == rand_blue_female_sample_IDs[5] | sample_ID == rand_blue_female_sample_IDs[6] | 
             sample_ID == rand_blue_female_sample_IDs[7] | sample_ID == rand_blue_female_sample_IDs[8] | 
             sample_ID == rand_blue_female_sample_IDs[9] | sample_ID == rand_blue_female_sample_IDs[10] | 
             sample_ID == rand_blue_female_sample_IDs[11])
  
  rand_brown_females<- females %>%
    filter(sample_ID != rand_blue_female_sample_IDs[1] & sample_ID != rand_blue_female_sample_IDs[2] &
             sample_ID != rand_blue_female_sample_IDs[3] & sample_ID != rand_blue_female_sample_IDs[4] &
             sample_ID != rand_blue_female_sample_IDs[5] & sample_ID != rand_blue_female_sample_IDs[6] & 
             sample_ID != rand_blue_female_sample_IDs[7] & sample_ID != rand_blue_female_sample_IDs[8] &
             sample_ID != rand_blue_female_sample_IDs[9] & sample_ID != rand_blue_female_sample_IDs[10] &
             sample_ID != rand_blue_female_sample_IDs[11])
  
  # Then, concatenate your "blue" and "brown" males and make sure that they are labelled in your new dataset.
  
  rand_data_females<-rbind(rand_blue_females,rand_brown_females)
  rand_blue_labels<-rep("rand_blue",length(rand_blue_females$sample_ID))
  rand_brown_labels<-rep("rand_brown",length(rand_brown_females$sample_ID))
  rand_color<-c(rand_blue_labels,rand_brown_labels)
  rand_data_females$rand_color<-rand_color
  
  
  # Put the males and females back together
  
  rand_data<-rbind(rand_data_males,rand_data_females)
  
  # Then, run your analysis and store the result for the rand_color*sex interaction
  
  cat("starting model",i)
  nbmodel_rand<-glmer.nb(as.numeric(psite_count_sc)~depth_ft_sc+sex*rand_color+(1|sampling_region/sampling_location/sample_ID)+
                           (1 + sex*rand_color|psite_spp)+offset(log(TL_cm_sc)),data=rand_data,family="nbinom")
  
  fixed_effects<-fixef(nbmodel_rand)
  rand_estimate_interaction[i]<-fixed_effects[5]
  rand_estimate_color[i]<-fixed_effects[4]
  rand_estimate_sex[i]<-fixed_effects[3]
  rand_z_interaction[i]<-coef(summary(nbmodel_rand))[5,3]
  rand_z_color[i]<-coef(summary(nbmodel_rand))[4,3]
  rand_z_sex[i]<-coef(summary(nbmodel_rand))[3,3]
  cat("ending model",i)
  
}

# Then, compare your actual value for the color*sex interaction to the distribution of randomly-created interaction values

final_500<-cbind(rand_estimate_interaction,rand_estimate_color,rand_estimate_sex,rand_z_interaction,rand_z_color,rand_z_sex)
write.csv(final_500,"final_500_simulated_values_males_and_females.csv")

# Then, compare your actual value for the color*sex interaction to the distribution of randomly-created interaction values
str(final_500)
final_500<-as.data.frame(final_500)

final_500<-read.table("final_500_simulated_values_males_and_females.csv",header=T,sep=",")
prob_z_interaction<-length(final_500$rand_z_interaction[abs(final_500$rand_z_interaction)>=abs(actual_z_interaction)])/length(final_500$rand_estimate_interaction)

cat(" The probability value for the interaction is ",prob_z_interaction, "\n")



