# Cargar paquetes
library(tidyverse)
library(ggsignif)

# ---- Crear y limpiar el dataframe ----
# Control (algunos valores faltantes) 6.203174523, NA, NA, NA, 3.429225148, NA, NA, 11.24040789, 8.375072597  
# SI
df_raw <- data.frame(
  Sample = rep(c("Seed", "Shoot", "Root"), each = 9),
  Replicate = rep(c("R1.1", "R1.2", "R1.3", "R2.1", "R2.2", "R2.3", "R3.1", "R3.2", "R3.3"), times = 3),
  Value = c(
    1963.90987, 1198.15121, 914.6908556, 5444.940617, 7598.611202, 4668.43729, 3987.758766, 7080.157693, 2454.572822,  # Germinated Seed
    2620.210613, 2039.810602, 25865.20081, 7889.73164, 7547.703061, 6432.538435, 5621.671275, 7577.16525, 6010.782745,  # Shoot
    20475.94175, 27022.26093, 4206.089136, 24042.31852, 17951.63622, 18216.14443, 9328.259501, 9555.408674, 31396.5367  # Root
    )
)
#FP
df_raw <- data.frame(
  Sample = rep(c("Seed", "Shoot", "Root"), each = 9),
  Replicate = rep(c("R1.1", "R1.2", "R1.3", "R2.1", "R2.2", "R2.3", "R3.1", "R3.2", "R3.3"), times = 3),
  Value = c(
    3587.121624, 1865.726333, 6505.288015, 33223.39426, 1706.520855, 4583.529174, 8448.896695, 2147.662332, 596.0899504,
    2884.470506, 2874.305143, 2690.133513, 147.04525, 12.55353391, 21.71606835, 85.53817891, 35.65521011, 12.91546361,
    14349.04005, 10372.99559, 18881.09007, 31578.31845, 13850.74073, 311.0455554, 15910.76275, 713.0336075, 6438.728312
  )
)
df_raw <- df_raw %>%
  mutate(LogValue = log10(Value))
df_summary <- df_raw %>%
  group_by(Sample) %>%
  summarise(
    log_mean = mean(LogValue, na.rm = TRUE),
    log_se = sd(LogValue, na.rm = TRUE) / sqrt(n())
  )
anova_model <- aov(log10(Value) ~ Sample, data = df_raw)
tukey_df <- as.data.frame(TukeyHSD(anova_model)$Sample)
tukey_df$comparison <- rownames(tukey_df)
comparisons_list <- list(
  c("Seed", "Shoot"),
  c("Seed", "Root"),
  c("Shoot", "Root")
)
tukey_df$annotation <- case_when(
  tukey_df$`p adj` < 0.001 ~ "***",
  tukey_df$`p adj` < 0.01  ~ "**",
  tukey_df$`p adj` < 0.05  ~ "*",
  TRUE ~ "ns"
)
annotations <- tukey_df$annotation  # <--- esto faltaba
y_positions <- c(5.8, 6.2, 6.6)  # Ajusta según tu eje y
df_summary$Sample <- factor(df_summary$Sample, levels = c("Seed", "Shoot", "Root"))

p <-ggplot(df_summary, aes(x = Sample, y = log_mean)) +
  geom_bar(stat = "identity", fill = "skyblue", width = 0.6) +
  geom_errorbar(aes(ymin = log_mean - log_se, ymax = log_mean + log_se), width = 0.2) +
  geom_jitter(data = df_raw, aes(x = Sample, y = log10(Value)), 
              width = 0.1, alpha = 0.7, color = "darkblue")+
  geom_signif(
    comparisons = comparisons_list,
    annotations = annotations,
    y_position = y_positions,
    tip_length = 0.01,
    textsize = 4) +
  theme_minimal()


p

library(rvg)
library(officer)
# Crear el PowerPoint con dimensiones específicas para el gráfico
doc <- read_pptx() %>%
  add_slide() %>%
  ph_with(
    dml(ggobj = p), 
    location = ph_location(
      left = 1,      # Posición desde el borde izquierdo (en pulgadas)
      top = 1,       # Posición desde el borde superior (en pulgadas)
      width = 10,   # Ancho del gráfico (en pulgadas)
      height = 6   # Alto del gráfico (en pulgadas)
    )
  )

# Guardar el archivo PPTX
print(doc, target = "Wheat_FP_plot_noControl.pptx")

