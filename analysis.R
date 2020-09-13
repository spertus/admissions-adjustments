library(tidyverse)
# source functions
source("functions.R")

# Read in data
scores <- read_csv("../Data/applicants-genders.csv") %>%
  rename(applicant = `Applicant ID`, reviewer = Reviewer, score = Score, reviewer_sex = `Reviewer Sex`, applicant_gender = `Applicant Gender Identity`) %>%
  mutate(applicant = as_factor(applicant)) %>%
  mutate(applicant_gender = ifelse(applicant_gender == "Male", "Male", "NonMale"))

############## EDA ##############

# Visualize score distribution by reviewer
ggplot(scores, aes(x = reviewer, y = score)) + geom_boxplot() + 
  labs(x = "Reviewer", y = "Score") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 16)) 

# Summary statistics
reviewer_summary <- scores %>% 
  group_by(reviewer) %>%
  summarise(count = n())
mean(reviewer_summary$count)
median(reviewer_summary$count)
min(reviewer_summary$count)
max(reviewer_summary$count)

########## coin flip test ############

# Fisher's combining function -- or should really do nonparametric combination of tests here?
coin_scores <- run_coin_flip_test(scores)
chi_sq <- -2*sum(log(coin_scores$p.value))
combined_pvalue <- 1 - pchisq(chi_sq, nrow(coin_scores))

# Visualize p-values
ggplot(data = coin_scores, aes(reviewer, p.value)) + 
  geom_point() +
  labs(x = "Reviewer", y = "P-value") +
  geom_hline(yintercept = 0.05, col = "red") +
  geom_hline(yintercept = 0.05/22, col = "orange") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 16)) 

# holm bonferroni correction
coin_scores <- coin_scores %>% arrange(p.value) %>%
  mutate(rank = 1:nrow(coin_scores),
         reviewer = factor(reviewer, levels = reviewer),
         threshold = rank*0.05/nrow(coin_scores))

ggplot(data = coin_scores, aes(reviewer, p.value)) + 
  geom_line(data = coin_scores, aes(x = reviewer, y = threshold), col = "orange", group = 1, size = 1.05) +
  labs(x = "Reviewer", y = "P-value") +
  geom_hline(yintercept = 0.05, col = "red", size = 1.05, linetype = "dashed") +
  geom_point(size = 1.5) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 16))

############## estimate scores ##############
estimates <- estimate_all_scores(score_frame = scores) %>%
  rename("Observed Mean" = "observed_mean", "ANOVA" = "ANOVA_score", "LSE" = "LSE_score")

#score plots
ggplot(estimates, aes(x = `Observed Mean`, y = ANOVA)) + 
  geom_point(size = 1.5) +
  geom_abline(slope = 1, intercept = 0) +
  xlim(15,45) + 
  ylim(15,45) +
  theme_bw() +
  theme(text = element_text(size = 16))
ggplot(estimates, aes(x = `Observed Mean`, y = LSE)) + 
  geom_point(size = 1.5) +
  geom_abline(slope = 1, intercept = 0) +
  xlim(15,45) + 
  ylim(15,45) +
  theme_bw() +
  theme(text = element_text(size = 16))
ggplot(estimates, aes(x = ANOVA, y = LSE)) + 
  geom_point(size = 1.5) +
  geom_abline(slope = 1, intercept = 0) +
  xlim(15,45) + 
  ylim(15,45) +
  theme_bw() +
  theme(text = element_text(size = 16))

#rank plots
ranks <- compute_ranks(estimates)
ggplot(ranks, aes(x = `Observed Mean`, y = ANOVA)) + 
  geom_point() +
  geom_abline(slope = 1, intercept = 0) +
  geom_hline(yintercept = 100, linetype = 'dashed') +
  geom_vline(xintercept = 100, linetype = 'dashed') +
  theme_bw() +
  xlab("Observed Mean Ranks") +
  ylab("ANOVA Ranks") +
  theme(text = element_text(size = 14))
ggplot(ranks, aes(x = `Observed Mean`, y = LSE)) + 
  geom_point() +
  geom_abline(slope = 1, intercept = 0) +
  geom_hline(yintercept = 100, linetype = 'dashed') +
  geom_vline(xintercept = 100, linetype = 'dashed') +
  theme_bw() +
  xlab("Observed Mean Ranks") +
  ylab("LSE Ranks") +
  theme(text = element_text(size = 14))
ggplot(ranks, aes(x = ANOVA, y = LSE)) + 
  geom_point() +
  geom_abline(slope = 1, intercept = 0) +
  geom_hline(yintercept = 100, linetype = 'dashed') +
  geom_vline(xintercept = 100, linetype = 'dashed') +
  theme_bw() +
  xlab("ANOVA Ranks") +
  ylab("LSE Ranks") +
  theme(text = element_text(size = 14))

#how many applicants are above rank 200 in all 3 estimators?
top_applicants <- ranks %>%
  filter(ANOVA <= 100, `Observed Mean` <= 100, LSE <= 100)



########## test for gender bias in reviews #########
#reviewer S only reviewed Non Male applicants, drop reviewer S from this test
scores_without_S <- scores %>% filter(reviewer != "S")

reviewer_gender_bias_NAworstcase <- run_gender_test(score_frame = scores_without_S, na_handling = "worst-case")
reviewer_gender_bias_NAomit <- run_gender_test(score_frame = scores_without_S, na_handling = "omit")
