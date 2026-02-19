library(tidyverse)
library(ggthemes)
library(RColorBrewer)

df <- read_csv("data/meta/all_data_oct_2025.csv")

greens_pal <- brewer.pal(9, "Greens")

vigor_palete <- c(greens_pal[4:9])

# data_manip --------------------------------------------------------------



time_df <- df %>% 
  pivot_longer(
    cols = starts_with("greenhouse_height_"),       
    names_to = "month",                  
    names_prefix = "greenhouse_height_",            
    values_to = "height"                 
  )%>%
  mutate(
    month = tolower(month),
    month = recode(month,
                   "febuary" = "february"),
    date = as.Date(paste0("2025-", match(month, tolower(month.name)), "-01")),
    month_num = case_when(
      month == "february" ~ 2,
      month == "march"    ~ 3,
      month == "april"    ~ 4,
      month == "august"   ~ 8,
      month == "october" ~10
    ),
    inocula_vigor = factor(inocula_vigor)
    )


# early figs height --------------------------------------------------------------

df_summary_time <- time_df %>%
  mutate(month = floor_date(date, "month")) %>%
  group_by(inocula_vigor, month) %>%
  summarise(
    mean_height = mean(height, na.rm = TRUE),
    se_height = sd(height, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  ) %>% 
  mutate(inocula_vigor = tidyr::replace_na(inocula_vigor, "Control"))

  
## time series figure ------------------------------------------------------

ggplot(df_summary_time, aes(x = month, y = mean_height, color = inocula_vigor)) +
  geom_point(size = 2) +
  geom_line()+
  geom_errorbar(aes(ymin =  mean_height- se_height,
                    ymax = mean_height + se_height), 
                width = 3
  )+
  scale_color_manual(values = vigor_palete)+
  scale_x_date(date_labels = "%b %Y", date_breaks = "1 month") +
  labs(
    x = "Date",
    y = "Plant height (cm)",
    color = "Donor Vigor"
  ) +
  theme_few()

a <- df_summary_time %>% 
  filter(month == "2025-08-01")

ggplot(a, aes(x = inocula_vigor, y = mean_height, fill = inocula_vigor))+
  geom_col()+
  geom_errorbar(aes(ymin =  mean_height- se_height,
                    ymax = mean_height + se_height), 
                width = .7
  )+
  scale_fill_manual( values = vigor_palete)+
  labs(
    x = "Inocula Vigor",
    y = "Plant height in October(cm)",
    fill = "Inocula Vigor"
  )+
  theme_few()


# volume  -----------------------------------------------------------------
vol_df <- df %>% 
  pivot_longer(
    cols = c(starts_with("greenhouse_height_"), starts_with("greenhouse_rcd_")),
    names_to = c(".value", "month"),
    names_pattern = "^greenhouse_(height|rcd)_(august|october)$"
  ) %>% 
  filter(!is.na(month)) %>% 
  mutate(
    month = tolower(month),
    month_num = match(month, tolower(month.name)),
    date = if_else(
      is.na(month_num),
      as.Date(NA),
      as.Date(sprintf("2025-%02d-01", month_num))
    ),
    inocula_vigor = factor(inocula_vigor)
  ) %>% 
  mutate(
    volume = (1/3)*(height)*((rcd/10)^2)
  )

df_summary_volume <- vol_df %>%
  mutate(month = floor_date(date, "month")) %>%
  group_by(inocula_vigor, month) %>%
  summarise(
    mean_volume = mean(volume, na.rm = TRUE),
    se_volume = sd(volume, na.rm = TRUE) / sqrt(n()),
    mean_rcd = mean(rcd, na.rm = TRUE),
    se_rcd = sd(rcd, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )


ggplot((df_summary_volume %>% filter(month == "2025-10-01")), 
       aes(x = inocula_vigor, y = mean_volume, fill = inocula_vigor))+
  geom_col()+
  geom_errorbar(aes(ymin =  mean_volume- se_volume,
                    ymax = mean_volume + se_volume), 
                width = .7
  )+
  scale_fill_manual( values = vigor_palete)+
  labs(
    x = "Inocula Vigor",
    y = "Plant volume in October(cm^2)",
    fill = "Inocula Vigor"
  )+
  theme_few()


ggplot(df_summary_volume, aes(x = month, y = mean_volume, color = inocula_vigor)) +
  geom_point(size = 2) +
  geom_line()+
  geom_errorbar(aes(ymin =  mean_volume - se_volume,
                    ymax = mean_volume + se_volume), 
                width = 3
  )+
  scale_color_manual(values = vigor_palete)+
  scale_x_date(date_labels = "%b %Y", date_breaks = "1 month") +
  labs(
    x = "Date",
    y = "Plant volume (cm)",
    color = "Donor Vigor"
  ) +
  theme_few()



  