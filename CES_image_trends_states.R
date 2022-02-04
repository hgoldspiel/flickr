# load image data with assigned themes from cluster analysis
load("data/cluster_analysis_output.RData")

# spatiotemporal summary stats for CES/non-CES images 
CES.themes <- c("scenery", "biota", "aquatics")

states.CES.hourly <- hourly.trends.fun(rural_photos_clustered[rural_photos_clustered$theme %in% CES.themes,])
states.nCES.hourly <- hourly.trends.fun(rural_photos_clustered[rural_photos_clustered$theme %notin% CES.themes,])
states.all.hourly <- hourly.trends.fun(rural_photos_clustered)

# bar plots separated by state
colors <- c("All images" = "grey", 
            "CES" = "forestgreen", 
            "Non-CES" = "black")

hourly.CES.trends.plot <-
  ggplot(states.CES.hourly, aes(x = Hour, y = Freq)) + 
  geom_col(data = states.all.hourly, inherit.aes = FALSE,
           aes(x = Hour, y = Freq, fill = "All images")) +
  geom_col(data = states.all.hourly, inherit.aes = FALSE,
           aes(x = Hour, y = -Freq, fill = "All images")) +
  geom_col(aes(fill = "CES"), width = 0.5) + 
  geom_col(data = states.nCES.hourly, width = 0.5,
           inherit.aes = FALSE, aes(x = Hour, y = -Freq, fill = "Non-CES")) +
  facet_grid(~State) +
  theme_bw() + mythemes + labs(x = NULL, y = NULL, fill = NULL) +
  scale_fill_manual(values = colors) +
  scale_x_continuous(breaks = 0:24, labels = c("0", as.character(1:24)), 
                     expand = c(.002,0)) +
  coord_polar(start = -.13, clip = "off") +
  theme(legend.position = "bottom", legend.text = element_text(size = 12),
        panel.border = element_blank(), strip.background = element_blank(), 
        axis.text.x = element_text(size = 11), axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(), panel.spacing = unit(1, "lines"))

## monthly image trends
states.CES.monthly <- monthly.trends.fun(rural_photos_clustered[rural_photos_clustered$theme %in% CES.themes,])
states.nCES.monthly <- monthly.trends.fun(rural_photos_clustered[rural_photos_clustered$theme %notin% CES.themes,])
states.all.monthly <- monthly.trends.fun(rural_photos_clustered)

monthly.CES.trends.plot <- 
  ggplot(states.CES.monthly, aes(x = factor(Month), y = Freq)) + 
  geom_col(data = states.all.monthly, inherit.aes = FALSE,
           aes(x = factor(Month), y = Freq, fill = "All images")) +
  geom_col(data = states.all.monthly, inherit.aes = FALSE,
           aes(x = factor(Month), y = -Freq, fill = "All images")) +
  geom_col(aes(fill = "CES"), width = 0.5) + 
  geom_col(data = states.nCES.monthly, width = 0.5,
           inherit.aes = FALSE, aes(x = Month, y = -Freq, fill = "Non-CES")) +
  facet_grid(~State) +
  scale_fill_manual(values = colors) +
  theme_bw() + mythemes + labs(x = NULL, y = NULL, fill = NULL) +
  scale_x_discrete(breaks = c(1:12), labels = month.abb) + 
  coord_polar(start = 6, clip = "off") +
  theme(legend.position = "bottom", legend.text = element_text(size = 12),
        panel.border = element_blank(), panel.spacing = unit(1, "lines"),
        strip.background = element_blank(), strip.text.x = element_blank(), 
        axis.text.x = element_text(size = 11), axis.text.y = element_blank(), 
        axis.ticks.y = element_blank())


## annual image trends
rural_photos_clustered$year <- year(rural_photos_clustered$date)
states.CES.yearly <- yearly.trends.fun(rural_photos_clustered[rural_photos_clustered$theme %in% CES.themes,])
states.nCES.yearly <- yearly.trends.fun(rural_photos_clustered[rural_photos_clustered$theme %notin% CES.themes,])
states.all.yearly <- yearly.trends.fun(rural_photos_clustered)

yearly.CES.trends.plot <- 
  ggplot(states.CES.yearly, aes(x = as.numeric(as.character(Year)), y = Freq)) + 
  geom_point(aes(color = "CES"), size = 2) + 
  geom_line(aes(color = "CES")) +
  geom_point(data = states.all.yearly, inherit.aes = FALSE, size = 2,
             aes(x = as.numeric(as.character(Year)), y = Freq, color = "All images")) +
  geom_line(data = states.all.yearly, inherit.aes = FALSE,
            aes(x = as.numeric(as.character(Year)), y = Freq, color = "All images")) +
  geom_point(data = states.nCES.yearly, inherit.aes = FALSE, size = 2,
             aes(x = as.numeric(as.character(Year)), y = Freq, color = "Non-CES")) +
  geom_line(data = states.nCES.yearly, inherit.aes = FALSE,
            aes(x = as.numeric(as.character(Year)), y = Freq, color = "Non-CES")) +
  facet_grid(~State) +
  scale_color_manual(values = colors) +
  scale_x_continuous(breaks = c(2012, 2014, 2016)) +
  facet_grid(~State) +
  theme_bw() + mythemes + labs(x = NULL, y = "Images (n)", fill = NULL) +
  theme(legend.position = "none", legend.text = element_text(size = 12), 
        strip.background = element_blank(), strip.text.x = element_blank(), 
        axis.text.x = element_text(size = 11))


ggarrange(hourly.CES.trends.plot, monthly.CES.trends.plot, yearly.CES.trends.plot,
          common.legend = TRUE, legend = "bottom",
          nrow = 3, labels = c("A", "B", "C"), align = "v",
          label.y = c(1,1,1.1))

ggsave("figures/rural_image_CES_trends.png", width = 8, height = 7, dpi = 600)