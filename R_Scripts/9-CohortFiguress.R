oneDose <- as_tibble(read.delim(paste0(dir_results,'Cohorts/one_dose_crude.txt'), sep = " "))
twoDose <- as_tibble(read.delim(paste0(dir_results, "Cohorts/two_dose_crude.txt"), sep = " "))


table <- oneDose
name <- "oneDose"
name1 <- "A) Seroconversion - One dose"

p1 <- table |> 
  select(days = Days_since_vaccine) |>
  mutate(days = as.numeric(days)) |>
  mutate(mean_days = mean(days)) |>
  ggplot(aes(x = days)) + 
  geom_bar(fill = "#123456") +
  labs(x = "Days between vaccination date and antibody test date", y = "Counts", title =  name1) +
  theme_bw() +
  scale_x_continuous(expand = c(0.01, 0.01)) +
  scale_y_continuous(expand = c(0.01, 0.05)) +
  geom_vline(aes(xintercept = mean_days), linetype = "dashed", colour = "red", size = 1.5)

ggsave(plot = p1, filename = paste0(dir_results,"Figures/", name, "_histogram.png"), width = 20, height = 15, units = "cm")
