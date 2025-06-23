library(phyloseq)
library(tidyverse)
library(biomformat)
library(vegan)
library(ggplot2)
library(ComplexUpset)
library(reshape2)
library(ggsci)
library(pheatmap)
library(FSA)
library(pairwiseAdonis)

#----Phyloseq Object----
# Leer biom
biom <- read_biom("asv-table-decont.biom")
otu_table <- otu_table(as(biom_data(biom), "matrix"), taxa_are_rows = TRUE)

# Leer taxonomía
taxa <- read.delim("taxonomy.tsv", sep = "\t", header = TRUE)
tax_split <- strsplit(as.character(taxa$Taxon), ";\\s*")
max_len <- max(lengths(tax_split))
tax_mat <- do.call(rbind, lapply(tax_split, function(x) { length(x) <- max_len; return(x) }))
rownames(tax_mat) <- taxa$Feature.ID
tax_table <- tax_table(as.matrix(tax_mat))

# Leer metadatos
metadata <- read.table("sample-metadata.tsv", sep = "\t", header = TRUE, row.names = 1)
sample_data <- sample_data(metadata)

# Crear objeto phyloseq
ps <- phyloseq(otu_table, tax_table, sample_data)


colnames(tax_table(ps)) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species") #Renombrar columnas taxonómicas
ps <- tax_glom(ps, taxrank = "Genus") # Agrupar a nivel de genero
# Sustituir g__Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium por g__Rhizobium
tax_table(ps)[tax_table(ps)[, "Genus"] == "g__Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium", "Genus"] <- "g__Rhizobium"
tax_table(ps)[tax_table(ps)[, "Genus"] == "g__Methylobacterium-Methylorubrum", "Genus"] <- "g__Methylorubrum"

#----SEED PCOA ----
ps_seed <- subset_samples(ps, Seed == "Yes") # Filtrar para quedarte solo con las muestras tipo "Seed"
ps_seed <- prune_taxa(taxa_sums(ps_seed) > 0, ps_seed) # Quitar OTUs sin abundancia después del filtrado
ordination_seed <- ordinate(ps_seed, method = "PCoA", distance = "bray") #Recalcular 
ordination_data_seed <- as.data.frame(ordination_seed$vectors) # Extraer coordenadas del PCoA
ordination_data_seed$Sample <- rownames(ordination_data_seed)
sample_meta <- data.frame(sample_data(ps_seed))# Extraer metadatos y unirlos
sample_meta$Sample <- rownames(sample_meta)
ordination_data_seed <- merge(ordination_data_seed, sample_meta, by = "Sample")# Merge para unir coordenadas y metadatos
ordination_data_seed <- subset(ordination_data_seed, !is.na(WheatSpecies)) # Quitar filas con NA en Wheat_specie

palette <- c(
  'A.tauschii' = "#fc2123",                 # Rojo fuerte
  'T.aestivum' = "#fc8d59",    # Naranja claro
  'T.aestivum_C' = "#fee08b",                 # Amarillo pálido
  'T.dicoccum' = "#91cf60",             # Verde medio
  'T.turgidum' = "darkgreen",             # Verde oscuro
  'T.durum' = "blue",          # Azul claro
  'T.monococcum' = "#762a83",               # Púrpura intenso
  'T.spelta' = "#f285c4"        # Magenta intermedio (heredado de Sphingomonas)
)
palette <- c(
  'America' = "#fc8d59",    # Naranja claro
  'Europe' = "#91cf60",             # Verde oscuro
  'Asia' = "blue",          # Azul claro
  'Oceania' = "#f285c4"        # Magenta intermedio (heredado de Sphingomonas)
)
p <- ggplot(ordination_data_seed, aes(x = Axis.1, y = Axis.2, color = WheatSpecies, shape = PloidyLevel)) +
  geom_point(size = 3) +
  #geom_text(aes(label = Sample), size = 2, vjust = -0.5) +
  theme_minimal() +
  scale_color_manual(values = palette) +
  labs(title = "PCoA - Bray-Curtis (Solo Semillas)", x = "PC1", y = "PC2", color = "WheatSpecies") +
  theme(legend.position = "right")
p

#PERMANOVA
dist_bray <- phyloseq::distance(ps_seed, method = "bray")
meta <- data.frame(sample_data(ps_seed))
#Full model
adonis_result <- vegan::adonis2(
  dist_bray ~ Domestication + PloidyLevel + WheatSpecies + Location,
  data = meta,
  permutations = 999
)
print(adonis_result)
#Individual model
adonis_result <- vegan::adonis2(
  dist_bray ~ Location,
  data = meta,
  permutations = 999
)
print(adonis_result)
# Obtener proporción de varianza explicada por los dos primeros ejes
variance_explained <- 100 * ordination_seed$values$Relative_eig[1:2]
names(variance_explained) <- c("PC1", "PC2")
print(variance_explained)
#----SEED TaxabarPlot----
ps_seed_rel <- transform_sample_counts(ps_seed, function(x) x / sum(x)) #Convertir a abundancia relativa
df <- psmelt(ps_seed_rel) #Pasar a data frame
df$Genus <- gsub("^g__", "", df$Genus)# Eliminar el prefijo 'g__' de los nombres de géneros
genus_mean <- df %>%
  group_by(Genus) %>%
  summarise(mean_abund = mean(Abundance)) %>%
  arrange(desc(mean_abund))
top_genera <- genus_mean %>%
  slice_max(mean_abund, n = 15) %>%
  pull(Genus) #Seleccionar los 10 géneros más abundantes 
df$Genus_filtered <- ifelse(df$Genus %in% top_genera, as.character(df$Genus), "Other") #Etiquetar los menos abundantes como "Other"
genus_order <- genus_mean %>%
  mutate(Genus = ifelse(Genus %in% top_genera, Genus, "Other")) %>%
  group_by(Genus) %>%
  summarise(mean_abund = mean(mean_abund)) %>%
  arrange(mean_abund) %>%
  pull(Genus)
df$Genus_filtered <- factor(df$Genus_filtered, levels = genus_order) #Reordenar Genus_filtered para que lo más abundante esté abajo en la leyenda y arriba en las barras
df$Sample <- reorder(df$Sample, -df$Abundance, sum) #Reordenar muestras por abundancia total
df$Sample <- sample_data(ps_seed_rel)$sample_id[match(df$Sample, rownames(sample_data(ps_seed_rel)))] #para poner los nombres del paper en el eje x

Genus.palette <- c(
  'Pantoea' = "#e41a1c",                       # Rojo fuerte (original)
  'Pseudomonas' = "#ffffbf",                   # Amarillo fuerte (original)
  'Staphylococcus' = "#d9ef8b",                # Verde amarillento (original)
  'Sphingomonas' = "#fee08b",                  # Amarillo pálido (original)
  'Bacteroides' = "#3288bd",                   # Azul claro (original)
  'Curtobacterium' = "#fc8d59",                # Naranja claro (original)
  'Massilia' = "#b695c0",                      # Lila (nuevo, similar a azul/rosa)
  'Haemophilus' = "#7FFFD4",                   # Aquamarine (nuevo, verde azulado)
  'Paenibacillus' = "#9e0142",                 # Magenta oscuro (original)
  'Hymenobacter' = "#B57EDC",                  # Lavanda (nuevo, violeta claro)
  'Streptococcus' = "#4CAF50",                 # Verde bosque (nuevo, similar a Enterobacter)
  'Enterobacter' = "#91cf60",                  # Verde medio (original)
  'Afipia' = "#9966CC",                        # Amatista (nuevo, violeta intermedio)
  'Rhizobium' = "darkgreen",                   # Verde pasto (original)
  'Methylobacterium-Methylorubrum' = "#de77ae", # Rosa claro (tomado de Rhodanobacter original)
  "Other" ="#E0E0E0"
)
p<-ggplot(df, aes(x = SampleName, y = Abundance, fill = Genus_filtered)) +
  geom_bar(stat = "identity") +
  facet_grid(~ PloidyLevel, scales = "free_x", space = "free_x") +
  theme_light() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 18),
    axis.text.y = element_text(size = 18),
    axis.title.x = element_text(size = 18),   
    axis.title.y = element_text(size = 18),
    strip.text = element_text(size = 18),
    legend.title = element_text(size = 18),
    legend.text = element_text(face = "italic", size = 18)
  ) +
  scale_fill_manual(values = Genus.palette) + 
  labs(
    x = "Sample",
    y = "Relative Abundance",
    fill = "Genus"
    #title = "Top 15 Genera present in Seeds"
  )
p

df_summary <- df %>%
  group_by(Genus_filtered) %>%
  summarise(total_abund = sum(Abundance)) %>%
  mutate(percent = round(total_abund / sum(total_abund) * 100, 2)) %>%
  arrange(desc(percent))
print(df_summary)

#----SEED Alfa_Diversidad_Wheat_species----
ps <- subset_samples(ps, !is.na(WheatSpecies) & WheatSpecies != "Unknown")#quitamos NAs y Unknown
sample_data(ps)$WheatSpecies <- factor(sample_data(ps)$WheatSpecies, 
                                       levels=c("A.tauschii", "T.monococcum", "T.turgidum", "T.dicoccum", "T.durum",
                                                "T.aestivum", "T.spelta"))
p <- plot_richness(ps, x = "WheatSpecies", measures = "Shannon", color = "WheatSpecies") +
  geom_boxplot(alpha = 0.6, width = 0.5) + 
  theme_minimal() +  # o theme_void() si quieres aún menos
  theme(
    legend.position = "none",
    panel.background = element_rect(fill = "transparent", colour = NA),
    plot.background = element_rect(fill = "transparent", colour = NA),
    axis.text.x = element_text(face = "italic", size = 9, family = "Arial"),
    axis.text.y = element_text(face = "plain", size = 9, family = "Arial"),
    axis.title.x = element_text(face = "plain", size = 10, family = "Arial"),
    axis.title.y = element_text(face = "plain", size = 10, family = "Arial"),
    strip.text = element_blank()
  ) +
  labs(
    x = "Wheat species", 
    y = "Alpha diversity (Shannon Index)"
  ) +
  scale_color_manual(values = c(
    "A.tauschii"="#260F99", "T.aestivum" = "#0F6B99", 
    "T.turgidum" = "#6B990F", "T.dicoccum" = "#F4D166", 
    "T.durum" = "#99540F", "T.monococcum"= "#990F0F", 
    "T.spelta"= "#CC5500"
  ))

p
# Analisis estadistico
shannon_df <- estimate_richness(ps_seed, measures = "Shannon")
shannon_df$Location <- sample_data(ps_seed)$Location
#Si no son normales se ejecuta el test de krustal wallis para ver si hay diferencias entre grupos. 
kruskal.test(Shannon ~ Location, data = shannon_df)

#----SEED Alfa_Diversidad_DomesticationStatus----
ps <- subset_samples(ps, !is.na(Domestication) & Domestication != "Unknown")#quitamos NAs y Unknown
sample_data(ps)$Domestication <- factor(sample_data(ps)$Domestication, 
                                       levels=c("Ancestral", "Commercial"))
p <- plot_richness(ps, x = "Domestication", measures = "Shannon", color = "Domestication") +
  geom_boxplot(alpha = 0.6, width = 0.5) + 
  theme_minimal() +  # o theme_void() si quieres aún menos
  theme(
    legend.position = "none",
    panel.background = element_rect(fill = "transparent", colour = NA),
    plot.background = element_rect(fill = "transparent", colour = NA),
    axis.text.x = element_text(face = "plain", size = 9, family = "Arial"),
    axis.text.y = element_text(face = "plain", size = 9, family = "Arial"),
    axis.title.x = element_text(face = "plain", size = 10, family = "Arial"),
    axis.title.y = element_text(face = "plain", size = 10, family = "Arial"),
    strip.text = element_blank()
  ) +
  labs(
    x = "Domestication", 
    y = "Alpha diversity (Shannon Index)"
  ) +
  scale_color_manual(values = c(
    "Ancestral"="#260F99", "Commercial"= "#990F0F"
  ))

p
#----SEED Alfa_Diversidad_PloidyLevel----
ps <- subset_samples(ps, !is.na(PloidyLevel) & PloidyLevel != "Unknown")#quitamos NAs y Unknown
sample_data(ps)$PloidyLevel <- factor(sample_data(ps)$PloidyLevel, 
                                        levels=c("2n", "4n", "6n"))
p <- plot_richness(ps, x = "PloidyLevel", measures = "Shannon", color = "PloidyLevel") +
  geom_boxplot(alpha = 0.6, width = 0.5) + 
  theme_minimal() +  # o theme_void() si quieres aún menos
  theme(
    legend.position = "none",
    panel.background = element_rect(fill = "transparent", colour = NA),
    plot.background = element_rect(fill = "transparent", colour = NA),
    axis.text.x = element_text(face = "plain", size = 9, family = "Arial"),
    axis.text.y = element_text(face = "plain", size = 9, family = "Arial"),
    axis.title.x = element_text(face = "plain", size = 10, family = "Arial"),
    axis.title.y = element_text(face = "plain", size = 10, family = "Arial"),
    strip.text = element_blank()
  ) +
  labs(
    x = "PloidyLevel", 
    y = "Alpha diversity (Shannon Index)"
  ) +
  scale_color_manual(values = c(
    "2n"="#260F99", "4n"= "#990F0F", "6n" = "#F4D166"
  ))

p
#----SEED Alfa_Diversidad_Location----
ps <- subset_samples(ps, !is.na(Location) & Location != "Unknown")#quitamos NAs y Unknown
sample_data(ps)$Location <- factor(sample_data(ps)$Location, 
                                       levels=c("Europe", "America","Asia", "Oceania"))
p <- plot_richness(ps, x = "Location", measures = "Shannon", color = "Location") +
  geom_boxplot(alpha = 0.6, width = 0.2) + 
  theme_bw() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(face = "plain", size = 9, family = "Arial"),  # Cambia el tamaño del texto del eje x
    axis.text.y = element_text(face = "plain", size = 9, family = "Arial"),             # Cambia el tamaño del texto del eje y
    axis.title.x = element_text(face = "plain", size = 10, family = "Arial"),             # Cambia el tamaño del título del eje x
    axis.title.y = element_text(face = "plain", size = 10, family = "Arial"),              # Cambia el tamaño del título del eje y
    strip.text = element_blank()
  ) +
  labs(
    x = "Location", 
    y = "Alpha diversity (Shannon Index)"
  ) +
  scale_color_manual(values = c("Europe"="#260F99","Asia" = "#6B990F", "Oceania" = "#F4D166", 
                                 "America"= "#990F0F"))

p
#----SEED Abundance-Occupancy----
otu <- otu_table(ps_seed)

otu_PA <- 1*((otu>0)==1)                                               # presence-absence data
otu_occ <- rowSums(otu_PA)/ncol(otu_PA)                                # occupancy calculation
otu_rel <- apply(decostand(otu, method="total", MARGIN=2),1, mean)     # mean relative abundance
occ_abun <- data.frame(otu_occ=otu_occ, otu_rel=otu_rel) %>%           # combining occupancy and abundance data frame
  rownames_to_column('otu')

p<-ggplot(data=occ_abun, aes(x=log10(otu_rel), y=otu_occ, color = otu_occ)) +
  geom_point(alpha = 0.95) +
  geom_text(aes(label = otu), size = 2, vjust = -0.5) +
  xlim (-5,0) +
  scale_color_gradient(low = "grey", high = "red") +
  labs(x="log10(mean relative abundance)", y="Occupancy",
       color = "Abundance")+
  theme(
    axis.title = element_text(size = 9),      # tamaño títulos ejes
    axis.text = element_text(size = 8),       # tamaño etiquetas ejes
    legend.title = element_text(size = 9),    # tamaño título leyenda
    legend.text = element_text(size = 8)      # tamaño texto leyenda
  )+
  theme_light()

p #"AbundanceOccupancy_Seed"
ggsave("Ocupancy_Seed_text.pdf", plot = p, width = 15, height = 8, units = "cm")

features_a_buscar <- c(
  "fc5860decf6565201a343d865b30951f",
  ) #para buscar los puntos q aparecen como abundancia-ocupancy distribution
resultados <- taxa[taxa$Feature.ID %in% features_a_buscar, c("Feature.ID", "Taxon")] # Filtrar el data.frame para quedarte solo con esos Feature.ID
print(resultados) # Imprimir resultados

#----SEED UpsetR NO Core----
otu_matrix <- as(otu_table(ps_seed), "matrix")
if (!identical(colnames(otu_matrix), sample_names(ps_seed))) {
  otu_matrix <- t(otu_matrix)  # Asegura que OTUs están en filas, muestras en columnas
}
muestras_turgidum  <- sample_names(ps_seed)[sample_data(ps_seed)$WheatSpecies == "T.turgidum"]
muestras_spelta    <- sample_names(ps_seed)[sample_data(ps_seed)$WheatSpecies == "T.spelta"]
muestras_durum     <- sample_names(ps_seed)[sample_data(ps_seed)$WheatSpecies == "T.durum"]
muestras_aestivum  <- sample_names(ps_seed)[sample_data(ps_seed)$WheatSpecies == "T.aestivum"]
muestras_tauschii  <- sample_names(ps_seed)[sample_data(ps_seed)$WheatSpecies == "A.tauschii"]
muestras_dicoccum    <- sample_names(ps_seed)[sample_data(ps_seed)$WheatSpecies == "T.dicoccum"]
muestras_mono     <- sample_names(ps_seed)[sample_data(ps_seed)$WheatSpecies == "T.monococcum"]

otus_turgidum  <- unique(rownames(otu_matrix)[rowSums(otu_matrix[, muestras_turgidum, drop = FALSE]) > 0])
otus_spelta <- unique(rownames(otu_matrix)[rowSums(otu_matrix[, muestras_spelta, drop = FALSE]) > 0])
otus_durum <- unique(rownames(otu_matrix)[rowSums(otu_matrix[, muestras_durum, drop = FALSE]) > 0])
otus_aestivum  <- unique(rownames(otu_matrix)[rowSums(otu_matrix[, muestras_aestivum, drop = FALSE]) > 0])
otus_tauschii <- unique(rownames(otu_matrix)[rowSums(otu_matrix[, muestras_tauschii, drop = FALSE]) > 0])
otus_dicoccum <- unique(rownames(otu_matrix)[rowSums(otu_matrix[, muestras_dicoccum, drop = FALSE]) > 0])
otus_mono  <- unique(rownames(otu_matrix)[rowSums(otu_matrix[, muestras_mono, drop = FALSE]) > 0])

# Crear una tabla binaria de presencia/ausencia
otu_presence <- data.frame(
  T.turgidum  = taxa_names(ps_seed) %in% otus_turgidum,
  T.durum = taxa_names(ps_seed) %in% otus_durum,
  T.spelta = taxa_names(ps_seed) %in% otus_spelta,
  T.aestivum  = taxa_names(ps_seed) %in% otus_aestivum,
  A.tauschii = taxa_names(ps_seed) %in% otus_tauschii,
  T.dicoccum = taxa_names(ps_seed) %in% otus_dicoccum,
  T.monococcum  = taxa_names(ps_seed) %in% otus_mono
)
rownames(otu_presence) <- taxa_names(ps_seed)# Asignar nombres de OTUs

otu_presence <- as.data.frame(lapply(otu_presence, as.integer))  # Convertir a 0/1
# Convertir 0/1 a valores lógicos TRUE/FALSE
otu_presence <- otu_presence %>%
  mutate(across(everything(), ~ . > 0))  
otu_presence_clean <- otu_presence %>%
  filter(rowSums(.) > 0) #eliminar los de suelo

p <- upset(
  otu_presence_clean, 
  intersect = c("T.aestivum", "T.spelta", "T.durum", "T.turgidum", "T.monococcum", "T.dicoccum", "A.tauschii"),
  base_annotations = list(
    "Intersection Size" = intersection_size()
  ),
  keep_empty_groups = FALSE  # Eliminar "Outside of Known Sets"
) + 
  theme_classic() +  
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14)
  ) + 
  labs(title = "Wheat Related Species ASVs distribution", 
       x = "Genotype", 
       y = "ASVs number")

p
otus_comunes <- Reduce(intersect, list(otus_spelta, otus_tauschii))
print(otus_comunes) #conocer los taxones comunes
features_a_buscar <- c(
#  "4711420ea5416c1b6935904216a58a50", Spingo
 # "f8ef89b0ef91e703f5a73299b597887b", Pseudomonas
  #"fc5860decf6565201a343d865b30951f", Pantoea
  #"1355705d5f5537fdd6063e8aef08650d" Stafi
  #"ab589cea8df6e4d57063760719baaacc" Bacteroides
  #"264be41ad2f4137d7cd5fefe393e1b92" Haemophilus
  "8bd283fe27d86af4e15aeebeae7b2a25"
)

resultados <- taxa[taxa$Feature.ID %in% features_a_buscar, c("Feature.ID", "Taxon")] # Filtrar el data.frame para quedarte solo con esos Feature.ID
print(resultados) # Imprimir resultados


#----GNOTO PCOA----
ps_gno <- subset_samples(ps, Gnotobiotic == "Yes") # Filtrar para quedarte solo con las muestras tipo "Gnotobiotioc"
ps_gno <- prune_taxa(taxa_sums(ps_gno) > 0, ps_gno) # Quitar OTUs sin abundancia después del filtrado
ordination_gno <- ordinate(ps_gno, method = "PCoA", distance = "bray") #Recalcular 
ordination_data_gno <- as.data.frame(ordination_gno$vectors) # Extraer coordenadas del PCoA
ordination_data_gno$Sample <- rownames(ordination_data_gno)
sample_meta <- data.frame(sample_data(ps_gno))# Extraer metadatos y unirlos
sample_meta$Sample <- rownames(sample_meta)
ordination_data_gno <- merge(ordination_data_gno, sample_meta, by = "Sample")# Merge para unir coordenadas y metadatos
ordination_data_gno <- subset(ordination_data_gno, !is.na(PlantPart)) # Quitar filas con NA en Wheat_specie

p <- ggplot(ordination_data_gno, aes(x = Axis.1, y = Axis.2, color = PlantPart)) +
  geom_point(size = 3) +
  stat_ellipse()+
  #geom_text(aes(label = Sample), size = 2, vjust = -0.5) +
  theme_minimal() +
  labs(title = "PCoA - Bray-Curtis Gnotobiotico", x = "PC1", y = "PC2", color = "Plant_part") +
  scale_color_manual(values = c(
    "Seed" = "#5C3A21",   # Marrón
    "Shoot" = "#4CAF50",   # Verde bosque
    "Root" = "#D4B85A"   # Verde bosque
  )) +
  theme(legend.position = "right")
p

#PERMANOVA
dist_bray <- phyloseq::distance(ps_gno, method = "bray")
meta <- data.frame(sample_data(ps_gno))
adonis_result <- vegan::adonis2(
  dist_bray ~ PlantPart,
  data = meta,
  permutations = 999
)
print(adonis_result)
pairwise.adonis2(dist_bray ~ PlantPart, data = meta)

# Obtener proporción de varianza explicada por los dos primeros ejes
variance_explained <- 100 * ordination_gno$values$Relative_eig[1:2]
names(variance_explained) <- c("PC1", "PC2")
print(variance_explained)

#----GNOTO TaxabarPlot----
ps_gno_rel <- transform_sample_counts(ps_gno, function(x) x / sum(x))
sample_meta <- as.data.frame(sample_data(ps_gno_rel))
plant_order <- c("Seed", "Root", "Shoot")
sample_meta$PlantPart <- factor(sample_meta$PlantPart, levels = plant_order)
df <- psmelt(ps_gno_rel)
df$Genus <- gsub("^g__", "", df$Genus)
df$SampleName <- sample_meta$SampleName[match(df$Sample, rownames(sample_meta))]
df$PlantPart <- sample_meta$PlantPart[match(df$Sample, rownames(sample_meta))]

genus_mean <- df %>%
  group_by(Genus) %>%
  summarise(mean_abund = mean(Abundance)) %>%
  arrange(desc(mean_abund))
top_genera <- genus_mean %>%
  slice_max(mean_abund, n = 15) %>%
  pull(Genus)
df$Genus_filtered <- ifelse(df$Genus %in% top_genera, as.character(df$Genus), "Other")
genus_order <- genus_mean %>%
  mutate(Genus = ifelse(Genus %in% top_genera, Genus, "Other")) %>%
  group_by(Genus) %>%
  summarise(mean_abund = mean(mean_abund)) %>%
  arrange(mean_abund) %>%
  pull(Genus)
df$Genus_filtered <- factor(df$Genus_filtered, levels = genus_order)

# Ordenar muestras primero por PlantPart, luego por abundancia total
sample_order <- df %>%
  group_by(SampleName, PlantPart) %>%
  summarise(total_abundance = sum(Abundance), .groups = "drop") %>%
  arrange(PlantPart, desc(total_abundance)) %>%
  pull(SampleName)
df$SampleName <- factor(df$SampleName, levels = sample_order)

Genus.palette <- c(
  'Pantoea' = "#e41a1c",                      # Rojo fuerte (original)
  'Pseudomonas' = "#ffffbf",                  # Amarillo fuerte (original)
  'Staphylococcus' = "#d9ef8b",               # Verde amarillento (original)
  'Sphingomonas' = "#fee08b",                 # Amarillo pálido (original)
  'Mycobacterium' = "#4a0c6b",                  # Azul púrpura oscuro (original)
  'Curtobacterium' = "#fc8d59",               # Naranja claro (original)
  'Massilia' = "#b695c0",                     # Verde oscuro (original)
  'Cloacibacterium' = "#5288bd",                  # Azul claro (original)
  'Paenibacillus' = "#9e0142",                # Magenta oscuro (original)
  'Hymenobacter' = "#B57EDC",                 # Lavanda (nuevo, violeta claro)
  'Streptococcus' = "#de77ae",                # Rosa claro (original)
  'Enterobacter' = "#91cf60",                 # Azul violáceo (original)
  'Afipia' = "#1966CC",                       # Amatista (nuevo, violeta intermedio)
  'Rhizobium' = "darkgreen",                  # Verde pasto (original)
  'Methylorubrum' = "#20B2AA",                # Lila (nuevo, violeta suave)
  "Other" = "#E0E0E0" 
)

p <- ggplot(df, aes(x = SampleName, y = Abundance, fill = Genus_filtered)) +
  geom_bar(stat = "identity") +
  facet_grid(~ PlantPart, scales = "free_x", space = "free_x") +
  theme_light() +
  theme(
    axis.text.x = element_text(hjust = 1, angle=45, size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 12),   
    axis.title.y = element_text(size = 12),
    strip.text = element_text(size = 12),
    legend.title = element_text(size = 12),
    legend.text = element_text(face = "italic", size =11)
  ) +
  scale_fill_manual(values = Genus.palette) + 
  labs(
    x = "Sample",
    y = "Relative Abundance",
    fill = "Genus"
  )

p
df_summary <- df %>%
  group_by(PlantPart, Genus_filtered) %>%
  summarise(total_abund = sum(Abundance), .groups = "drop") %>%
  group_by(PlantPart) %>%
  mutate(percent = round(total_abund / sum(total_abund) * 100, 2)) %>%
  arrange(PlantPart, desc(percent))
print(df_summary)
#write.csv(df_summary, "GNO_abundancia_por_genero_y_parte.csv", row.names = FALSE)

#----GNOTO Alfa_Diversidad----
# Asegurémonos de que la columna "Wheat_specie" esté presente en los metadatos
sample_data(ps_gno)$PlantPart <- factor(sample_data(ps_gno)$PlantPart, 
                                       levels=c("Seed", "Root", "Shoot"))
p <- plot_richness(ps_gno, x = "PlantPart", measures = "Shannon", color = "PlantPart") +
  geom_boxplot(alpha = 0.6, width = 0.5) + 
  theme_bw() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(face = "plain", size = 7, family = "Arial"),  # Cambia el tamaño del texto del eje x
    axis.text.y = element_text(face = "plain", size = 7, family = "Arial"),             # Cambia el tamaño del texto del eje y
    axis.title.x = element_text(face = "plain", size = 8, family = "Arial"),             # Cambia el tamaño del título del eje x
    axis.title.y = element_text(face = "plain", size = 8, family = "Arial"),              # Cambia el tamaño del título del eje y
    strip.text = element_blank()
  ) +
  labs(
    x = "Tissue", 
    y = "Alpha diversity (Shannon Index)"
  ) +
  scale_color_manual(values = c("Seed" = "#5C3A21", "Shoot" = "#4CAF50", "Root" = "#D4B85A"))

p
# Analisis estadistico
shannon_df <- estimate_richness(ps_gno, measures = "Shannon")
shannon_df$PlantPart <- sample_data(ps_gno)$PlantPart
by(shannon_df$Shannon, shannon_df$PlantPart, shapiro.test) 
#Comprobamos la homogeneidad de las varianzas
bartlett.test(Shannon ~ PlantPart, data = shannon_df)
#Si son normales, se realiza un Anova
anova_result <- aov(Shannon ~ PlantPart, data = shannon_df)
summary(anova_result)
#Si sale significativo se hace un test post hoc
TukeyHSD(anova_result)

#----FIELD PCOA----
ps_field <- subset_samples(ps, Field == "Yes") # Filtrar para quedarte solo con las muestras tipo "Gnotobiotioc"
ps_field <- prune_taxa(taxa_sums(ps_field) > 0, ps_field) # Quitar OTUs sin abundancia después del filtrado
ordination_field <- ordinate(ps_field, method = "PCoA", distance = "bray") #Recalcular 
ordination_data_field <- as.data.frame(ordination_field$vectors) # Extraer coordenadas del PCoA
ordination_data_field$Sample <- rownames(ordination_data_field)
sample_meta <- data.frame(sample_data(ps_field))# Extraer metadatos y unirlos
sample_meta$Sample <- rownames(sample_meta)
ordination_data_field <- merge(ordination_data_field, sample_meta, by = "Sample")# Merge para unir coordenadas y metadatos
ordination_data_field <- subset(ordination_data_field, !is.na(PlantPart)) # Quitar filas con NA en Wheat_specie

p <- ggplot(ordination_data_field, aes(x = Axis.1, y = Axis.2, color = PlantPart, shape=WheatSpecies)) +
  geom_point(size = 3) +
  stat_ellipse()+
  theme_minimal() +
  labs(title = "PCoA - Bray-Curtis Field", x = "PC1", y = "PC2", color = "PlantPart") +
  scale_color_manual(values = c("Root"="#FFA500", "Soil" = "#0F6B99", "Shoot" = "#6B990F", 
                                "Seed" = "#99540F", "Spike"= "#FFDB58"
  )) +
  theme(legend.position = "right")
p
#PERMANOVA
dist_bray <- phyloseq::distance(ps_field, method = "bray")
meta <- data.frame(sample_data(ps_field))
adonis_result <- vegan::adonis2(
  dist_bray ~ PlantPart,
  data = meta,
  permutations = 999
)
print(adonis_result)
pairwise.adonis2(dist_bray ~ PlantPart, data = meta)

# Obtener proporción de varianza explicada por los dos primeros ejes
variance_explained <- 100 * ordination_field$values$Relative_eig[1:2]
names(variance_explained) <- c("PC1", "PC2")
print(variance_explained)

#----Field TaxabarPlot----
ps_filtered_rel <- transform_sample_counts(ps_field, function(x) x / sum(x))
df <- psmelt(ps_filtered_rel)
df$Genus <- gsub("^g__", "", df$Genus)
genus_mean <- df %>%
  group_by(Genus) %>%
  summarise(mean_abund = mean(Abundance)) %>%
  arrange(desc(mean_abund))
top_genera <- genus_mean %>%
  slice_max(mean_abund, n = 15) %>%
  pull(Genus)
df$Genus_filtered <- ifelse(df$Genus %in% top_genera, as.character(df$Genus), "Other")
genus_order <- genus_mean %>%
  mutate(Genus = ifelse(Genus %in% top_genera, Genus, "Other")) %>%
  group_by(Genus) %>%
  summarise(mean_abund = mean(mean_abund)) %>%
  arrange(mean_abund) %>%
  pull(Genus)
df_summary <- df %>%
  group_by(WheatSpecies, Genus_filtered, PlantPart) %>%
  summarise(mean_abund = mean(Abundance), .groups = "drop")
genus_order <- df_summary %>%
  group_by(Genus_filtered) %>%
  summarise(mean_abund = mean(mean_abund)) %>%
  arrange(mean_abund) %>%
  pull(Genus_filtered)
df_summary$Genus_filtered <- factor(df_summary$Genus_filtered, levels = genus_order)
df_summary$PlantPart <- factor(df_summary$PlantPart, levels = c("Soil", "Root", "Shoot", "Spike", "Seed"))
df_summary$mean_abund <- as.numeric(df_summary$mean_abund)

Genus.palette <- c(
  'Pantoea' = "#e41a1c",            # Del seed y gno
  'Pseudomonas' = "#ffffbf",        # Del seed y gno
  'Sphingomonas' = "#fee08b",       # Del seed y gno
  'Enterobacter' = "#91cf60",       # Del seed (puedes cambiar a gno: #5e4fa2 si prefieres)
  'Rhizobium' = "darkgreen",        # Del seed y gno
  'Curtobacterium' = "#fc8d59",     # Del seed y gno
  'Massilia' = "#b695c0",           # Del seed (puedes cambiar a gno: #1a9850 si prefieres)
  'Staphylococcus' = "#d9ef8b",     # Del seed y gno
  'Nocardioides' = "#8DA0CB",       # Nuevo, color por definir
  'Pseudarthrobacter' = "#FFD92F",  # Nuevo, color por definir
  'Bacteroides' = "#3288bd",        # Del seed
  'Mycobacterium' = "#4a0c6b",      # Del gno
  'Hymenobacter' = "#B57EDC",       # Del seed y gno
  'Leifsonia' = "#A6D854",          # Nuevo, color por definir
  'Bacillus' = "#E78AC3"      ,     # Nuevo, color por definir
  "Other" ="grey"
)
p<-ggplot(df_summary, aes(x = PlantPart, y = mean_abund, fill = Genus_filtered)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_grid(~WheatSpecies) +  # Añadir facet_wrap para separar por especie de trigo
  theme_light() +
  theme(
    axis.text.x = element_text(hjust = 1, angle=45, size = 18),
    axis.text.y = element_text(size = 18),
    axis.title.x = element_text(size = 18),   
    axis.title.y = element_text(size = 18),
    strip.text = element_text(size = 18, face = "italic"),
    legend.title = element_text(size = 18),
    legend.text = element_text(face = "italic", size =18)
      ) +
  scale_fill_manual(values = Genus.palette) +
  labs(
    x = "Plant Tissue",
    y = "Mean Relative Abundance",
    fill = "Genus",
    #title = "Mean Relative Abundance of Top Genera per Tissue and Wheat Species"
  )
p

df_summary <- df %>%
  group_by(PlantPart, Genus_filtered) %>%
  summarise(total_abund = sum(Abundance), .groups = "drop") %>%
  group_by(PlantPart) %>%
  mutate(percent = round(total_abund / sum(total_abund) * 100, 2)) %>%
  arrange(PlantPart, desc(percent))
print(df_summary)
write.csv(df_summary, "ro_y_parte.csv", row.names = FALSE)


#----Field Alfa_Diversidad----
# Asegurémonos de que la columna "Wheat_specie" esté presente en los metadatos
sample_data(ps_field)$PlantPart <- factor(sample_data(ps_field)$PlantPart, 
                                         levels=c("Soil" ,"Root", "Shoot", 
                                                  "Spike", "Seed"))
p <- plot_richness(ps_field, x = "PlantPart", measures = "Shannon", color = "PlantPart") +
  geom_boxplot(alpha = 0.6, width = 0.5) + 
  theme_bw() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(face = "plain", size = 16, family = "Arial"),  # Cambia el tamaño del texto del eje x
    axis.text.y = element_text(face = "plain", size = 16, family = "Arial"),             # Cambia el tamaño del texto del eje y
    axis.title.x = element_text(face = "plain", size = 17, family = "Arial"),             # Cambia el tamaño del título del eje x
    axis.title.y = element_text(face = "plain", size = 17, family = "Arial"),              # Cambia el tamaño del título del eje y
    strip.text = element_blank()
  ) +
  labs(
    x = "PlantPart", 
    y = "Alpha diversity (Shannon Index)"
  ) +
  scale_color_manual(values = c("Soil" = "#0F6B99", 
                                "Root"="#FFA500", 
                                "Shoot" = "#6B990F", 
                                "Spike"= "#FFDB58",
                                "Seed" = "#99540F"))

p
# Analisis estadistico
shannon_df <- estimate_richness(ps_field, measures = "Shannon")
shannon_df$PlantPart <- sample_data(ps_field)$PlantPart
by(shannon_df$Shannon, shannon_df$PlantPart, shapiro.test) 
#Comprobamos la homogeneidad de las varianzas
bartlett.test(Shannon ~ PlantPart, data = shannon_df)
#Si son normales, se realiza un Anova
anova_result <- aov(Shannon ~ PlantPart, data = shannon_df)
summary(anova_result)
#Si sale significativo se hace un test post hoc
TukeyHSD(anova_result)

#----FIELD Abundance-Occupancy----
otu <- otu_table(ps_field)

otu_PA <- 1*((otu>0)==1)                                               # presence-absence data
otu_occ <- rowSums(otu_PA)/ncol(otu_PA)                                # occupancy calculation
otu_rel <- apply(decostand(otu, method="total", MARGIN=2),1, mean)     # mean relative abundance
occ_abun <- data.frame(otu_occ=otu_occ, otu_rel=otu_rel) %>%           # combining occupancy and abundance data frame
  rownames_to_column('otu')

p<-ggplot(data=occ_abun, aes(x=log10(otu_rel), y=otu_occ, color = otu_occ)) +
  geom_point(alpha = 0.95) +
  #geom_text(aes(label = otu), size = 2, vjust = -0.5) +
  xlim (-5,0) +
  scale_color_gradient(low = "grey", high = "red") +
  labs(title = "Abundance-Occupancy Distribution_Field samples",
       x="log10(mean relative abundance)", y="Occupancy",
       color = "Abundance")+
  theme_light()

p #"AbundanceOccupancy_Field"

occ_abun %>%
  arrange(desc(otu_occ)) %>%
  head(10)

features_a_buscar <- c(
  "fc5860decf6565201a343d865b30951f",
  "f8ef89b0ef91e703f5a73299b597887b",
  "4711420ea5416c1b6935904216a58a50",
  "18f23f300ca2237253c4703cf55ed601",
  "183c5b10e7b4afead83fc4e006886f37",
  "5251ba579faaf24f00548cbf01c04376",
  "92a4da6a12eb5275bc9880149003e7a6",
  "edfaf9bebbf2d56f902832bf67c724ec",
  "c3d6bd8829e50e91207303ffc8e67530",
  "fec74a96ae35bad77273b2d97898925c"
) #para buscar los puntos q aparecen como abundancia-ocupancy distribution
resultados <- taxa[taxa$Feature.ID %in% features_a_buscar, c("Feature.ID", "Taxon")] # Filtrar el data.frame para quedarte solo con esos Feature.ID
print(resultados) # Imprimir resultados



#----Field UpsetR per wheat specie 70%Core----
otu_matrix <- as(otu_table(ps_field), "matrix")
if (!identical(colnames(otu_matrix), sample_names(ps_field))) {
  otu_matrix <- t(otu_matrix)  # Asegura que OTUs están en filas, muestras en columnas
}
muestras_turgidum  <- sample_names(ps_field)[sample_data(ps_field)$WheatSpecies == "T.turgidum"]
muestras_spelta    <- sample_names(ps_field)[sample_data(ps_field)$WheatSpecies == "T.spelta"]
muestras_durum     <- sample_names(ps_field)[sample_data(ps_field)$WheatSpecies == "T.durum"]
muestras_aestivum  <- sample_names(ps_field)[sample_data(ps_field)$WheatSpecies == "T.aestivum"]

#Función para filtrar OTUs presentes en >=70 % de las muestras
get_otus_70p <- function(muestras) {
  sub_matrix <- otu_matrix[, muestras, drop = FALSE]
  min_muestras <- ceiling(length(muestras) * 0.7)
  otus <- rownames(sub_matrix)[rowSums(sub_matrix > 0) >= min_muestras]
  return(otus)
}
#Obtener OTUs por especie con >=70% de presencia
otus_turgidum  <- get_otus_70p(muestras_turgidum)
otus_spelta    <- get_otus_70p(muestras_spelta)
otus_durum     <- get_otus_70p(muestras_durum)
otus_aestivum  <- get_otus_70p(muestras_aestivum)

# Crear una tabla binaria de presencia/ausencia
otu_presence <- data.frame(
  T.turgidum  = taxa_names(ps_field) %in% otus_turgidum,
  T.durum = taxa_names(ps_field) %in% otus_durum,
  T.spelta = taxa_names(ps_field) %in% otus_spelta,
  T.aestivum  = taxa_names(ps_field) %in% otus_aestivum
)
rownames(otu_presence) <- taxa_names(ps_field)# Asignar nombres de OTUs

otu_presence <- as.data.frame(lapply(otu_presence, as.integer))  # Convertir a 0/1
# Convertir 0/1 a valores lógicos TRUE/FALSE
otu_presence <- otu_presence %>%
  mutate(across(everything(), ~ . > 0))  
otu_presence_clean <- otu_presence %>%
  filter(rowSums(.) > 0) #eliminar los de suelo

p <- upset(
  otu_presence_clean, 
  intersect = c("T.aestivum", "T.spelta", "T.durum", "T.turgidum"),
  base_annotations = list(
    "Intersection Size" = intersection_size()
  ),
  keep_empty_groups = FALSE  # Eliminar "Outside of Known Sets"
) + 
  theme_classic() +  
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14)
  ) + 
  labs(title = "Wheat Related Species ASVs distribution", 
       x = "Genotype", 
       y = "ASVs number")

p
otus_comunes <- Reduce(intersect, list(otus_turgidum, otus_spelta, otus_durum, otus_aestivum))
print(otus_comunes) #conocer los taxones comunes
features_a_buscar <- c(
  "fc5860decf6565201a343d865b30951f",
  "f8ef89b0ef91e703f5a73299b597887b"
)
resultados <- taxa[taxa$Feature.ID %in% features_a_buscar, c("Feature.ID", "Taxon")] # Filtrar el data.frame para quedarte solo con esos Feature.ID
print(resultados) # Imprimir resultados

#----Field UpsetR per tissue 70%Core----
ps_field <- subset_samples(ps_field, Field == "Yes") #filtramos las muestras de campo

otu_matrix <- as(otu_table(ps_field), "matrix")
if (!identical(colnames(otu_matrix), sample_names(ps_field))) {
  otu_matrix <- t(otu_matrix)  # Asegura que OTUs están en filas, muestras en columnas
}
muestras_seed  <- sample_names(ps_field)[sample_data(ps_field)$PlantPart == "Seed"] #Identificar las muestras por tipo de tejido
muestras_spike <- sample_names(ps_field)[sample_data(ps_field)$PlantPart == "Spike"]
muestras_shoot <- sample_names(ps_field)[sample_data(ps_field)$PlantPart == "Shoot"]
muestras_root  <- sample_names(ps_field)[sample_data(ps_field)$PlantPart == "Root"]

#Función para filtrar OTUs presentes en >=70 % de las muestras
get_otus_70p <- function(muestras) {
  sub_matrix <- otu_matrix[, muestras, drop = FALSE]
  min_muestras <- ceiling(length(muestras) * 0.7)
  otus <- rownames(sub_matrix)[rowSums(sub_matrix > 0) >= min_muestras]
  return(otus)
}
#Obtener OTUs por especie con >=70% de presencia
otus_seed  <- get_otus_70p(muestras_seed)
otus_spike    <- get_otus_70p(muestras_spike)
otus_shoot     <- get_otus_70p(muestras_shoot)
otus_root  <- get_otus_70p(muestras_root)

# Crear una tabla binaria de presencia/ausencia
otu_presence <- data.frame(
  Seed  = taxa_names(ps_field) %in% otus_seed,
  Spike = taxa_names(ps_field) %in% otus_spike,
  Shoot = taxa_names(ps_field) %in% otus_shoot,
  Root  = taxa_names(ps_field) %in% otus_root
)
rownames(otu_presence) <- taxa_names(ps_field)# Asignar nombres de OTUs

otu_presence <- as.data.frame(lapply(otu_presence, as.integer))  # Convertir a 0/1
# Convertir 0/1 a valores lógicos TRUE/FALSE
otu_presence <- otu_presence %>%
  mutate(across(everything(), ~ . > 0))  
otu_presence_clean <- otu_presence %>%
  filter(rowSums(.) > 0) #eliminar los de suelo

p <- upset(
  otu_presence_clean, 
  intersect = c("Seed", "Spike", "Shoot", "Root"),
  base_annotations = list(
    "Intersection Size" = intersection_size()
  ),
  keep_empty_groups = FALSE  # Eliminar "Outside of Known Sets"
) + 
  theme_classic() +  
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14)
  ) + 
  labs(title = "Tissue Related Species ASVs distribution", 
       x = "Plant part", 
       y = "ASVs number")

p

otus_comunes <- Reduce(intersect, list(otus_shoot, otus_spike, otus_root))
print(otus_comunes) #conocer los taxones comunes
features_a_buscar <- c(
  "fc5860decf6565201a343d865b30951f","5251ba579faaf24f00548cbf01c04376", "4711420ea5416c1b6935904216a58a50",
  "18f23f300ca2237253c4703cf55ed601", "f8ef89b0ef91e703f5a73299b597887b")
resultados <- taxa[taxa$Feature.ID %in% features_a_buscar, c("Feature.ID", "Taxon")] # Filtrar el data.frame para quedarte solo con esos Feature.ID
print(resultados) # Imprimir resultados

#----Field enrichment----
selected_genera <- c("Pantoea", "Sphingomonas", "Pseudomonas",
                     "Massilia", "Methylorubrum")

tissue_order <- c("Seed", "Spike", "Shoot", "Root", "Soil")

df_selected <- df %>%
  filter(Genus %in% selected_genera)

df_summary <- df_selected %>%
  group_by(Genus, Sample, PlantPart) %>%  # Agrupar por muestra también
  summarise(Abundance = sum(Abundance), .groups = "drop") %>%
  group_by(Genus, PlantPart) %>%
  summarise(
    mean_abund = mean(Abundance),
    se = sd(Abundance) / sqrt(n()),
    .groups = "drop"
  ) %>%
  mutate(
    Genus = factor(Genus, levels = selected_genera),
    PlantPart = factor(PlantPart, levels = tissue_order)
  )
df_summary <- df_summary %>%
  mutate(
    mean_abund = mean_abund * 100,
    se = se * 100
  )

tissue_palette <- c("Soil" = "#8B8970",   # Gris oscuro para Soil
                    "Root" = "#40E0D0",   # Turquesa para Root
                    "Shoot" = "#228B22",  # Verde para Shoot
                    "Spike" = "#FFD700",  # Amarillo para Spike
                    "Seed" = "#8B4513")   # Marrón para Seed

p <- ggplot(df_summary, aes(x = mean_abund, y = PlantPart, fill = PlantPart)) +
  geom_bar(stat = "identity", width = 0.6) +
  geom_errorbar(
    aes(xmin = mean_abund - se, xmax = mean_abund + se),
    width = 0.2,
    color = "black"
  ) +
  facet_wrap(~Genus, scales = "free_x", nrow = 1) +
  scale_fill_manual(values = tissue_palette) +
  scale_y_discrete(limits = rev(tissue_order)) +
  #scale_x_continuous(limits = c(0, 20)) +
  labs(
    x = "Relative Abundance (%)", 
    y = "Tissue"
  ) +
  theme_light() +
  theme(
    axis.text.x = element_text(size = 13),
    axis.text.y = element_text(size = 13),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    legend.position = "none",
    strip.text = element_text(size = 14, face = "italic"),
    strip.background = element_rect(fill = "lightgray", color = "gray"),
    panel.border = element_rect(color = "gray", size = 0.5),
    panel.spacing.x = unit(1, "lines"),
    plot.margin = margin(10, 10, 10, 10)
  )
p

#ANALISIS ESTADISTICO
df_selected <- df %>%
  filter(Genus %in% selected_genera)
df_stats <- df_selected %>%
  group_by(Genus, Sample, PlantPart) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop") %>%
  mutate(Abundance = Abundance * 100)  # Transformar a porcentaje si quieres
df_genus <- df_stats %>%
  filter(Genus == "Pantoea")
anova_result <- aov(Abundance ~ PlantPart, data = df_genus)
summary(anova_result)
TukeyHSD(anova_result)

#----PPT_Figure----
library(rvg)
library(officer)

doc <- read_pptx() %>%
  add_slide() %>%
  ph_with(
    dml(ggobj = p), 
    location = ph_location(
      left = 1,      # Posición desde el borde izquierdo (en pulgadas)
      top = 1,       # Posición desde el borde superior (en pulgadas)
      width = 15,   # Ancho del gráfico (en pulgadas)15 cm
      height = 7   # Alto del gráfico (en pulgadas)8 cm
    )
  )

# Guardar el archivo PPTX
print(doc, target = "Field_Upset_Genus_Tissue.pptx")

