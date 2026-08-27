directory = "/path/10_Publication_figures/04_rg_mixer_gwaspw/"

rg_female = "/path/04_LDSC/SNPrg/Females/rg_femaleAD_correlation_table.txt"

rg_male = "/path/04_LDSC/SNPrg/Males/rg_maleAD_correlation_table.txt"


### LOAD PACKAGES ###

library(patchwork)
library(ggforce)      # ggplot circles
library(showtext)
showtext_auto()
library(tidyverse)



### LOAD DATA ###

rg_female_df <- read.table(rg_female, header = T, stringsAsFactors = F)

rg_male_df <- read.table(rg_male, header = T, stringsAsFactors = F)

rg_df <- bind_rows(rg_female_df, rg_male_df) %>%
  mutate(CI = 1.96 * se,
         p1 = gsub("_sumstats", "", p1),
         p2 = gsub("_sumstats", "", p2))


### FOREST PLOT: Rg WITH OTHER TRAITS ###
order <- c('Schizophrenia',
           'PTSD',
           'Bipolar',
           'Anxiety',
           'ADHD',
           'EA',
           'BMI',
           'Waist_to_hip_ratio',
           'Metabolic_syndrome',
           'Drinks_per_week',
           'Ever_smoked_regularly')
order <- rev(order)

plot_df_traits <- rg_df %>%
  filter(p2 == 'Schizophrenia' |
           p2 == 'PTSD' |
           p2 == 'Bipolar' |
           p2 == 'Anxiety' |
           p2 == 'ADHD' |
           p2 == 'EA' |
           p2 == 'BMI' |
           p2 == 'Waist_to_hip_ratio' |
           p2 == 'Metabolic_syndrome' |
           p2 == 'Drinks_per_week' |
           p2 == 'Ever_smoked_regularly') %>%
  mutate(p2 = factor(p2, levels = order),
         Sex = factor(p1, levels = c("Allcohorts_female_AD", "Allcohorts_male_AD")))



rg_plot_traits <- ggplot(plot_df_traits, aes(x = rg, y = p2, colour = p1)) +
  geom_vline(xintercept = 0,
             colour = "grey30",
             linetype = 'dashed') +
  geom_point(position = position_dodge(width = -0.8)) +
  geom_linerange(aes(xmin = rg + CI, xmax = rg - CI),
                 position = position_dodge(width = -0.8)) +
  scale_colour_viridis_d("Sex", begin = 1, end = 0,
                         labels = c("Females", "Males")) +
  scale_y_discrete("",
                   breaks = c('Schizophrenia',
                              'PTSD',
                              'Bipolar',
                              'Anxiety',
                              'ADHD',
                              'EA',
                              'BMI',
                              'Waist_to_hip_ratio',
                              'Metabolic_syndrome',
                              'Drinks_per_week',
                              'Ever_smoked_regularly'),
                   labels = c('Schizophrenia',
                              'PTSD',
                              'Bipolar',
                              'Anxiety',
                              'ADHD',
                              'Educational Attainment',
                              'BMI',
                              'Waist/hip ratio',
                              'Metabolic syndrome',
                              'Drinks/week',
                              'Regular smoker')) +
  scale_x_continuous("Genetic Correlation",
                     breaks = c(-0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1),
                     labels = c(-0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1),
                     limits = c(-0.41, 1.01)) +
  annotate("text", x = 0.7, y = "ADHD", label = "*", size = 7 / 2.845276) +
  annotate("text", x = 0.3, y = "BMI", label = "*", size = 7 / 2.845276) +
  annotate("text", x = 0.35, y = "Metabolic_syndrome", label = "*", size = 7 / 2.845276) +
  annotate("text", x = 0.45, y = "Ever_smoked_regularly", label = "*", size = 7 / 2.845276) +
  theme_classic() +
  theme(text = element_text(family = "Arial"),
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 7),
        axis.title.x = element_text(size = 7, colour = "black", margin = margin(10,0,0,0)),
        axis.title.y = element_text(size = 7, colour = "black", margin = margin(0,10,0,0)),
        legend.position = "top",
        legend.title = element_text(size = 7, colour = "black"),
        legend.text = element_text(size = 7, colour = "black"))

rg_plot_traits


### FOREST PLOT: Rg WITH BMI SEX-SPECIFIC ###

order <- c('BMI_female',
           'BMI_male')
order <- rev(order)

plot_df_bmi <- rg_df %>%
  filter(p2 == 'BMI_female' |
           p2 == 'BMI_male') %>%
  mutate(p2 = factor(p2, levels = order),
         Sex = factor(p1, levels = c("Allcohorts_female_AD", "Allcohorts_male_AD")))



rg_plot_BMI <- ggplot(plot_df_bmi, aes(x = rg, y = p2, colour = p1)) +
  geom_vline(xintercept = 0,
             colour = "grey30",
             linetype = 'dashed') +
  geom_point(position = position_dodge(width = -0.8)) +
  geom_linerange(aes(xmin = rg + CI, xmax = rg - CI),
                 position = position_dodge(width = -0.8)) +
  scale_colour_viridis_d("Sex", begin = 1, end = 0,
                         labels = c("Females", "Males")) +
  scale_y_discrete("",
                   breaks = c('BMI_female',
                              'BMI_male'),
                   labels = c('BMI: Females',
                              'BMI: Males')) +
  scale_x_continuous("Genetic Correlation",
                     breaks = c(0, 0.1, 0.2, 0.3),
                     labels = c(0, 0.1, 0.2, 0.3),
                     limits = c(0, 0.35)) +
  annotate("text", x = 0.32, y = "BMI_female", label = "*", size = 7 / 2.845276) +
  annotate("text", x = 0.32, y = "BMI_male", label = "*", size = 7 / 2.845276) +
  theme_classic() +
  theme(text = element_text(family = "Arial"),
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 7),
        axis.title.x = element_text(size = 7, colour = "black", margin = margin(10,0,0,0)),
        axis.title.y = element_text(size = 7, colour = "black", margin = margin(0,10,0,0)),
        legend.position = "none")

rg_plot_BMI


### Mixer Venn: female AD-BMI ###

# Data: sizes of each group
one <- 4695  # AD exclusive
two <- 407   # BMI exclusive
overlap <- 8549  # Shared between AD and BMI

# Calculate circle areas
area_one <- one + overlap
area_two <- two + overlap

# Convert areas to radii (since area = π * r^2, radius = sqrt(area / π))
r_1 <- sqrt(area_one / pi)
r_2 <- sqrt(area_two / pi)

# Calculate the proportional overlap distance
overlap_ratio <- overlap / min(area_one, area_two)
d <- (r_1 + r_2) * (1 - overlap_ratio / 2)  # Adjust overlap distance

# Create circle data with calculated positions
circle_data <- data.frame(
  x = c(0, d),  # x-coordinates for each circle
  y = c(0, 0),  # Keep y-coordinates the same for both circles
  radius = c(r_1, r_2),
  group = c("AD", "BMI")
)

# Plot the Euler diagram with numbers positioned correctly
mixer_venn_F_AD_BMI <- ggplot(circle_data) +
  geom_circle(aes(x0 = x, y0 = y, r = radius, fill = group), alpha = 0.6, color = "black") +
  scale_fill_manual(values = c("#FDE725FF", "#B8DE29FF")) +
  coord_fixed() +
  # Add labels inside each circle and in the overlapping section
  annotate("text",
           x = circle_data$x[1] - r_1 / 2,
           y = circle_data$y[1],
           label = one,
           size = 7 / 2.845276, color = "black", family = "Arial", fontface = "bold") +
  annotate("text",
           x = circle_data$x[2] + r_2 / 2,
           y = circle_data$y[2],
           label = two,
           size = 7 / 2.845276, color = "black", family = "Arial", fontface = "bold") +
  annotate("text",
           x = mean(circle_data$x) + 3,
           y = circle_data$y[1],
           label = overlap,
           size = 7 / 2.845276, color = "black", family = "Arial", fontface = "bold") +
  # Add labels above the circles
  annotate("text",
           x = circle_data$x[1],
           y = circle_data$y[1] + r_1 + 10,
           label = "AD",
           size = 7 / 2.845276, color = "black", family = "Arial", fontface = "bold") +
  annotate("text",
           x = circle_data$x[2],
           y = circle_data$y[2] + r_2 + 10,
           label = "BMI",
           size = 7 / 2.845276, color = "black", family = "Arial", fontface = "bold") +
  
  theme_classic() +
  theme(text = element_text(family = "Arial"),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.title = element_blank(),
        legend.position = "none")

mixer_venn_F_AD_BMI


### Mixer Venn: female AD-metS ###

# Data: sizes of each group
one <- 3780  # AD exclusive
two <- 958   # metS exclusive
overlap <- 9464  # Shared between AD and metS

# Calculate circle areas
area_one <- one + overlap
area_two <- two + overlap

# Convert areas to radii (since area = π * r^2, radius = sqrt(area / π))
r_1 <- sqrt(area_one / pi)
r_2 <- sqrt(area_two / pi)

# Calculate the proportional overlap distance
overlap_ratio <- overlap / min(area_one, area_two)
d <- (r_1 + r_2) * (1 - overlap_ratio / 2)  # Adjust overlap distance

# Create circle data with calculated positions
circle_data <- data.frame(
  x = c(0, d),  # x-coordinates for each circle
  y = c(0, 0),  # Keep y-coordinates the same for both circles
  radius = c(r_1, r_2),
  group = c("AD", "BMI")
)

# Plot the Euler diagram with numbers positioned correctly
mixer_venn_F_AD_metS <- ggplot(circle_data) +
  geom_circle(aes(x0 = x, y0 = y, r = radius, fill = group), alpha = 0.6, color = "black") +
  scale_fill_manual(values = c("#FDE725FF", "#B8DE29FF")) +
  coord_fixed() +
  # Add labels inside each circle and in the overlapping section
  annotate("text",
           x = circle_data$x[1] - r_1 / 2,
           y = circle_data$y[1],
           label = one,
           size = 7 / 2.845276, color = "black", family = "Arial", fontface = "bold") +
  annotate("text",
           x = circle_data$x[2] + r_2 / 2,
           y = circle_data$y[2],
           label = two,
           size = 7 / 2.845276, color = "black", family = "Arial", fontface = "bold") +
  annotate("text",
           x = mean(circle_data$x) + 3,
           y = circle_data$y[1],
           label = overlap,
           size = 7 / 2.845276, color = "black", family = "Arial", fontface = "bold") +
  # Add labels above the circles
  annotate("text",
           x = circle_data$x[1],
           y = circle_data$y[1] + r_1 + 10,
           label = "AD",
           size = 7 / 2.845276, color = "black", family = "Arial", fontface = "bold") +
  annotate("text",
           x = circle_data$x[2],
           y = circle_data$y[2] + r_2 + 10,
           label = "metS",
           size = 7 / 2.845276, color = "black", family = "Arial", fontface = "bold") +
  
  theme_classic() +
  theme(text = element_text(family = "Arial"),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.title = element_blank(),
        legend.position = "none")

mixer_venn_F_AD_metS

### Mixer Venn: male AD-BMI ###

# Data: sizes of each group
one <- 2376  # AD exclusive
two <- 3680   # BMI exclusive
overlap <- 4735  # Shared between AD and BMI

# Calculate circle areas
area_one <- one + overlap
area_two <- two + overlap

# Convert areas to radii (since area = π * r^2, radius = sqrt(area / π))
r_1 <- sqrt(area_one / pi)
r_2 <- sqrt(area_two / pi)

# Calculate the proportional overlap distance
overlap_ratio <- overlap / min(area_one, area_two)
d <- (r_1 + r_2) * (1 - overlap_ratio / 2)  # Adjust overlap distance

# Create circle data with calculated positions
circle_data <- data.frame(
  x = c(0, d),  # x-coordinates for each circle
  y = c(0, 0),  # Keep y-coordinates the same for both circles
  radius = c(r_1, r_2),
  group = c("AD", "BMI")
)

# Plot the Euler diagram with numbers positioned correctly
mixer_venn_M_AD_BMI <- ggplot(circle_data) +
  geom_circle(aes(x0 = x, y0 = y, r = radius, fill = group), alpha = 0.6, color = "black") +
  scale_fill_manual(values = c("#440154FF", "#453781FF")) +
  coord_fixed() +
  # Add labels inside each circle and in the overlapping section
  annotate("text",
           x = circle_data$x[1] - r_1 / 2,
           y = circle_data$y[1],
           label = one,
           size = 7 / 2.845276, color = "white", family = "Arial", fontface = "bold") +
  annotate("text",
           x = circle_data$x[2] + r_2 / 2,
           y = circle_data$y[2],
           label = two,
           size = 7 / 2.845276, color = "white", family = "Arial", fontface = "bold") +
  annotate("text",
           x = mean(circle_data$x),
           y = circle_data$y[1],
           label = overlap,
           size = 7 / 2.845276, color = "white", family = "Arial", fontface = "bold") +
  # Add labels above the circles
  annotate("text",
           x = circle_data$x[1],
           y = circle_data$y[1] + r_1 + 10,
           label = "AD",
           size = 7 / 2.845276, color = "black", family = "Arial", fontface = "bold") +
  annotate("text",
           x = circle_data$x[2],
           y = circle_data$y[2] + r_2 + 10,
           label = "BMI",
           size = 7 / 2.845276, color = "black", family = "Arial", fontface = "bold") +
  
  theme_classic() +
  theme(text = element_text(family = "Arial"),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.title = element_blank(),
        legend.position = "none")

mixer_venn_M_AD_BMI


### Mixer Venn: male AD-metS ###

# Data: sizes of each group
one <- 1033  # AD exclusive
two <- 4344   # metS exclusive
overlap <- 6078  # Shared between AD and metS

# Calculate circle areas
area_one <- one + overlap
area_two <- two + overlap

# Convert areas to radii (since area = π * r^2, radius = sqrt(area / π))
r_1 <- sqrt(area_one / pi)
r_2 <- sqrt(area_two / pi)

# Calculate the proportional overlap distance
overlap_ratio <- overlap / min(area_one, area_two)
d <- (r_1 + r_2) * (1 - overlap_ratio / 2)  # Adjust overlap distance

# Create circle data with calculated positions
circle_data <- data.frame(
  x = c(0, d),  # x-coordinates for each circle
  y = c(0, 0),  # Keep y-coordinates the same for both circles
  radius = c(r_1, r_2),
  group = c("AD", "BMI")
)

# Plot the Euler diagram with numbers positioned correctly
mixer_venn_M_AD_metS <- ggplot(circle_data) +
  geom_circle(aes(x0 = x, y0 = y, r = radius, fill = group), alpha = 0.6, color = "black") +
  scale_fill_manual(values = c("#440154FF", "#453781FF")) +
  coord_fixed() +
  # Add labels inside each circle and in the overlapping section
  annotate("text",
           x = circle_data$x[1] - r_1 / 2,
           y = circle_data$y[1],
           label = one,
           size = 7 / 2.845276, color = "white", family = "Arial", fontface = "bold") +
  annotate("text",
           x = circle_data$x[2] + r_2 / 2,
           y = circle_data$y[2],
           label = two,
           size = 7 / 2.845276, color = "white", family = "Arial", fontface = "bold") +
  annotate("text",
           x = mean(circle_data$x) - 5,
           y = circle_data$y[1],
           label = overlap,
           size = 7 / 2.845276, color = "white", family = "Arial", fontface = "bold") +
  # Add labels above the circles
  annotate("text",
           x = circle_data$x[1],
           y = circle_data$y[1] + r_1 + 10,
           label = "AD",
           size = 7 / 2.845276, color = "black", family = "Arial", fontface = "bold") +
  annotate("text",
           x = circle_data$x[2],
           y = circle_data$y[2] + r_2 + 10,
           label = "metS",
           size = 7 / 2.845276, color = "black", family = "Arial", fontface = "bold") +
  
  theme_classic() +
  theme(text = element_text(family = "Arial"),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.title = element_blank(),
        legend.position = "none")

mixer_venn_M_AD_metS


### Venn: gwas-pw female vs male for shared AD/BMI regions ###
# Data: sizes of each group
female <- 24
male <- 0
overlap <- 1

sizes <- c(Female = female + overlap,
           Male = male + overlap)

# Calculate circle radii proportional to the square root of sizes
r_female <- sqrt(sizes["Female"] / pi)
r_male <- sqrt(sizes["Male"] / pi)

# Calculate offset to shift Male circle to the right so its right edge touches Female's right edge
offset <- r_female - r_male

# Create circle data with adjusted position for Male circle
circle_data <- data.frame(
  x = c(1, 1 + offset),  # Offset the Male circle to the right
  y = c(1, 1),           # Keep y-coordinates the same for female and male circle
  radius = c(r_female, r_male),
  group = c("Female", "Male")
)

# Plot the Euler diagram
gwaspw_AD_BMI_venn <- ggplot(circle_data) +
  geom_circle(aes(x0 = x, y0 = y, r = radius, fill = group), alpha = 0.6) +
  scale_fill_manual(values = c("#FDE725FF", "#440154FF")) +
  coord_fixed(expand = TRUE) +
  annotate("text", x = circle_data$x[1], y = circle_data$y[1],
           label = female, size = 7 / 2.845276, color = "black", family = "Arial", fontface = "bold") +
  annotate("text", x = circle_data$x[2], y = circle_data$y[2],
           label = overlap, size = 7 / 2.845276, color = "white", family = "Arial", fontface = "bold") +
  # Add label above the circles
  # annotate("text",
  #          x = circle_data$x[1],
  #          y = circle_data$y[1] + circle_data$radius[1] + 0.7,
  #          label = "Shared regions \nAD/BMI",
  #          size = 7 / 2.845276, color = "black", family = "Arial", fontface = "bold") +
  ggtitle("Shared regions \nAD/BMI") +
  theme_classic() +
  theme(text = element_text(family = "Arial"),
        plot.title = element_text(size = 7, hjust = 0.5, face = "bold"),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.title = element_blank(),
        legend.position = "none")

gwaspw_AD_BMI_venn

### Venn: gwas-pw female vs male for shared AD/metS regions ###
# Data: sizes of each group
female <- 22  # female exclusive
male <- 4   # male exclusive
overlap <- 4  # Shared between female and male

# Calculate circle areas
sizes <- c(Female = female + overlap,
           Male = male + overlap)

# Convert areas to radii (since area = π * r^2, radius = sqrt(area / π))
r_female <- sqrt(sizes["Female"] / pi)
r_male <- sqrt(sizes["Male"] / pi)

# Calculate the proportional overlap distance
overlap_ratio <- overlap / min(sizes["Female"], sizes["Male"])
d <- (r_female + r_male) * (1 - overlap_ratio / 2)  # Adjust overlap distance

# Create circle data with calculated positions
circle_data <- data.frame(
  x = c(0, d),  # x-coordinates for each circle
  y = c(0, 0),  # Keep y-coordinates the same for both circles
  radius = c(r_female, r_male),
  group = c("Female", "Male")
)

# Plot the Euler diagram with numbers positioned correctly
gwaspw_AD_metS_venn <- ggplot(circle_data) +
  geom_circle(aes(x0 = x, y0 = y, r = radius, fill = group), alpha = 0.6, color = "black") +
  scale_fill_manual(values = c("#FDE725FF", "#440154FF")) +
  coord_fixed(expand = TRUE) +
  # Add labels inside each circle and in the overlapping section
  annotate("text",
           x = circle_data$x[1] - r_female / 2,
           y = circle_data$y[1],
           label = female,
           size = 7 / 2.845276, color = "black", family = "Arial", fontface = "bold") +
  annotate("text",
           x = circle_data$x[2] + r_male / 2,
           y = circle_data$y[2],
           label = male,
           size = 7 / 2.845276, color = "white", family = "Arial", fontface = "bold") +
  annotate("text",
           x = mean(circle_data$x) + 0.6,
           y = circle_data$y[1],
           label = overlap,
           size = 7 / 2.845276, color = "white", family = "Arial", fontface = "bold") +
  ggtitle("Shared regions \nAD/metS") +
  theme_classic() +
  theme(text = element_text(family = "Arial"),
        plot.title = element_text(size = 7, hjust = 0.5, face = "bold"),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.title = element_blank(),
        legend.position = "none")


gwaspw_AD_metS_venn




### PANELLED FIGURE ###

figure <- ((rg_plot_traits |
              (plot_spacer() / rg_plot_BMI / plot_spacer() +
                 plot_layout(heights = c(0.2, 1, 0.2)))) +
             plot_layout(widths = c(2, 1))) /
  (((mixer_venn_F_AD_BMI | mixer_venn_F_AD_metS) /
      (mixer_venn_M_AD_BMI | mixer_venn_M_AD_metS)) |
     (gwaspw_AD_BMI_venn + theme(plot.margin = unit(c(5,5,15,15), "pt")) |
        gwaspw_AD_metS_venn + theme(plot.margin = unit(c(5,5,15,15), "pt")))) +
  plot_annotation(tag_levels = list(c("A", "B", "C", "", "", "", "D"))) &
  theme(plot.tag = element_text(size = 7, face = "bold", family = "Arial"),
        plot.margin = unit(c(10, 10, 10, 10), "pt") )

figure

outfile <- paste(directory, "Rg_traits_BMI_sex_specific_mixer_gwaspw.eps", sep="")
ggsave(figure, width = 18, height = 18.5, units = "cm", file = outfile, device = cairo_ps, dpi = 300)
