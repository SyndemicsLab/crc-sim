library(data.table)

sex <- fread("../Raw/OUDOrigin_Sex_12SEP2024.csv")
race <- fread("../Raw/OUDOrigin_Race_12SEP2024.csv")

sex_list <- c()
race_list <- c()
out_list <- c()
for (i in 1:100) {
  set.seed(2024 + i)
  for(j in unique(sex$final_sex)) {
    subDT <- Syndemics::crc(sex[final_sex == j & year == 2021, 
                                ][N_ID < 0, N_ID := sample(1:10, 1)],
                            "N_ID", c("Casemix", "APCD", "BSAS", "PMP", "Matris", "Death"))
    estimate <- subDT$estimate
    known <- as.numeric(sex[final_sex == j & year == 2021 & N_ID > 0,
                            ][, .(known = sum(N_ID))])
    sex_list[[j]] <- data.table(estimate = estimate, known = known, group = ifelse(j == 1, "Male", "Female"))
  } 
  sex_est <- rbindlist(sex_list)
  
  for(j in unique(race$final_re)) {
    subDT <- Syndemics::crc(race[final_re == j & year == 2021, 
                                ][N_ID < 0, N_ID := sample(1:10, 1)],
                            "N_ID", c("Casemix", "APCD", "BSAS", "PMP", "Matris", "Death"))
    estimate <- subDT$estimate
    known <- as.numeric(race[final_re == j & year == 2021 & N_ID > 0,
                            ][, .(known = sum(N_ID))])
    race_list[[j]] <- data.table(estimate = estimate, known = known, group = ifelse(j == 1, "White", 
                                                                                    ifelse(j == 2, "Black",
                                                                                           ifelse(j == 3, "Asian/PI",
                                                                                                  ifelse(j == 4, "Hispanic", "Other")))))
  } 
  race_est <- rbindlist(race_list)
  out_list[[i]] <- rbindlist(list(race_est, sex_est))
}
out_data <- rbindlist(out_list)
final <- out_data[, .(known = mean(known),
                      est = round(mean(estimate)),
                      se = sd(estimate) / sqrt(.N)), by = "group"
                  ][, `:=` (UCI = round(est + qt(0.975, est-1)*se),
                            LCI = round(est - qt(0.975, est-1)*se))]

total <- out_data[, total := estimate + known
                  ][group %in% c("White", "Black", "Hispanic", "Other", "Asian/PI"), gtype := "Race"
                    ][group %in% c("Male", "Female"), gtype := "Sex"
                      ][group %in% c("-20", "21-40", "41-60", "61-80", "81+"), gtype := "Age"
                        ][!group %in% c("Other")]

write.csv(final, "finalComp.csv", row.names = FALSE)
write.csv(out_data, "bootstrappeddata.csv", row.names = FALSE)

library(ggplot2)

ggplot(total, aes(x = group, y = total)) + 
  geom_boxplot() + 
  facet_grid(~gtype, scales = "free") + 
  theme_bw() + 
  labs(title = "Estimate of Total Prevalence for Opioid Use Disorder",
       subtitle = "State of Massachusetts, 2021")
ggsave("realEst.png", width = 12, height = 8)
