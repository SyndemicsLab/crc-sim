library(ggplot2)
library(data.table)

ECNS <- fread("../Data/evenCaptures_noStrata.csv")
ECS <- fread("../Data/evenCaptures_Strata.csv")
UECS <- fread("../Data/unevenCaptures_fiveStrata.csv")
UECNS <- fread("../Data/unevencaptures_noStrata.csv")

makeFig <- function(data) {
  if(length(unique(data$Group)) != 1) stratified <- TRUE else stratified <- FALSE
  DT <- data[, diff := Estimate - GroundTruth]
  
  if(!stratified) {
    plabs <- labs(x = "")
    pticks <- theme(axis.text.x = element_blank(),
                    axis.ticks.x = element_blank())
  } else {
    plabs <- labs(x = "Group")
    pticks <- theme()
  }
  
  
  p <- ggplot(DT, aes(y = diff, x = Group, group = Group)) + 
    geom_hline(yintercept = 0, col = "red", linetype = "dashed") +
    geom_boxplot() + 
    guides(color = "none",
           size = "none") + 
    theme_bw() + 
    facet_grid(Model~Method) + 
    plabs + 
    pticks + 
    labs(y = "Estimate and Ground Truth\nDifference")
  
  return(p)
}

makeFig(UECNS) + labs(title = "Uneven Capture Probabilities",
                      subtitle = "Total Population Estimate")
ggsave("../Figures/unevenCaptures_noStrata.png")

makeFig(ECNS) + labs(title = "Even Capture Probabilities",
                     subtitle = "Total Population Estimate")
ggsave("../Figures/evenCaptures_noStrata.png")

makeFig(UECS) + labs(title = "Uneven Capture Probabilities",
                     subtitle = "By-group Estimation")
ggsave("../Figures/unevenCptures_strata.png", width = 6, height = 16)

makeFig(ECS) + labs(title = "Even Capture Probabilities",
                    subtitle = "Total Population Estimate")
ggsave("../Figures/evenCptures_strata.png")


avgnobs <- copy(UECS)[, .(gt = mean(GroundTruth)), by = c("Group")]
pct <- copy(UECS)[, pct := (Estimate - GroundTruth)/GroundTruth*100
                  ][, Group := factor(Group, labels = round(avgnobs$gt[order(avgnobs$gt, decreasing = FALSE)]), levels = c(5, 3, 2, 4, 1), ordered = TRUE)]
ggplot(pct[!Model %like% "Both"], aes(x = Group, y = pct, group = Group)) + 
  geom_hline(yintercept = 0, col = "red", linetype = "dashed") + 
  geom_boxplot() + 
  theme_bw() + 
  facet_grid(Model~Method) + 
  labs(y = "Percent Difference in\nEstimate Against Ground Truth",
       x = "Ground Truth Population",
       title = "Capture-Recapture Over Multiple Datasets",
       subtitle = "Stratified over 5 groups of varying sizes") +
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 22),
        strip.text = element_text(size = 15))

# Many data emulation
final <- fread("../Data/PHDEmulation.csv")

ground_truth <- setorderv(ground_truth, "N_ID", order = 1L)
plot_data <- copy(final)[, pct_diff := (estimates - gt) / gt * 100
                         ][, model := ifelse(model == "poisson", "Poisson", "NB")
                           ][, params := paste0(stringr::str_to_title(direction), "-", threshold)
                             ][, group := factor(group, levels = ground_truth$strata, labels = ground_truth$N_ID)]
library(ggplot2)
ggplot(plot_data[!direction %like% "both"], aes(x = group, y = pct_diff, group = group)) + 
  geom_hline(yintercept = 0, col = "red", linetype = "dashed") +
  geom_boxplot() + 
  facet_grid(params~model) + 
  labs(x = "Ground Truth Population",
       y = "Percent Difference in\nEstimate from Ground Truth",
       title = "Single Dataset Bootstrapping",
       subtitle = "Boostrapped Suppression Values ") +
  theme_bw() +
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 22),
        strip.text = element_text(size = 15))
ggsave("../Figures/PHDEmulation.png", width = 8, height = 8)



