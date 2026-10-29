library(sf)
library(dplyr)
library(broom)
library(tidyverse) 
library(readr)

# Importation des données et préparation

TRI <- st_read("Données/Floods_Spatial_RiskZone_Template_-_4326_MultiPolygon.gpkg")

Caracteristiques_Tout <- read_delim("Données/Caractéristiques communes.csv", 
    delim = ";", escape_double = FALSE, locale = locale(date_names = "fr"), 
    trim_ws = TRUE)
Caracteristiques_Commune <- Caracteristiques_Tout %>% filter(GEO_OBJECT == "COM")

Caracteristiques_Commune <- Caracteristiques_Commune %>% filter(TAB_MEASURE %in% c("DWELLINGS_OCS_DW_MAIN_TPH_MF21", "MED_SL", "POP_EMPSTA_ENQ_1T2_AGE_Y15T64", "PR_MD60")) # Suppression des données inutiles pour l'analyse

# Les informations sur la situation dans les communes comprennent des résultats de plusieurs enquêtes. Seule la plus récente est retenue
Caracteristiques_Commune <- Caracteristiques_Commune[
 ave(as.numeric(as.character(Caracteristiques_Commune$TIME_PERIOD)),
 Caracteristiques_Commune$GEO,
 Caracteristiques_Commune$TAB_MEASURE,
 FUN = function(x) rank(-x, ties.method = "first")) == 1,
]

Contours_Communes <- st_read("Données/Contours communes/Communes.shp") # Contours des communes

# Alignement des représentations spatiales
crs_TRI <- st_crs(TRI)
crs_Contours <- st_crs(Contours_Communes)
TRI <- st_transform(TRI, crs_Contours)

# Corriger les erreurs dans les géométries des contours des communes
Contours_Communes <- Contours_Communes %>%
  mutate(geometry = st_make_valid(geometry))

# Vérifier la validité des géométries après les tentatives de corrections
sum(!st_is_valid(Contours_Communes))

## Calcul, pour chaque commune, la part du territoire exposé aux inondations. Le résultat sera entre 0 et 1. ##

# Calculer l'aire de chaque commune (en m² ou km² selon le CRS)
Contours_Communes <- Contours_Communes %>% mutate(Aire_Totale = st_area(.))

# Identifier les parties du territoire exposées au risque d'inondation
intersections <- st_intersection(Contours_Communes, TRI)

# Calculer l'aire de chaque intersection
intersections <- intersections %>%
 mutate(Aire_Exposée = st_area(.))

# Agrégation par commune
Aire_Exposée_par_Commune <- intersections %>% group_by(DCOE_C_COD) %>% summarise(Aire_Exposée_Totale = sum(Aire_Exposée, na.rm = TRUE)) %>% ungroup() 

# A l’aide de la surface du territoire de chaque commune et de l’étendu des surfaces inondables pour chaque commune, calculer un coefficient d’exposition aux inondations pour chaque commune. Ce coefficient sera stocké dans l’objet « Risque_Commune ».
Risque_Commune <- Contours_Communes %>%
  select(DCOE_C_COD, Aire_Totale) %>%
  left_join(st_drop_geometry(Aire_Exposée_par_Commune), by = "DCOE_C_COD") %>%
  mutate(
    Aire_Totale = as.numeric(Aire_Totale),
    Aire_Exposée_Totale = as.numeric(Aire_Exposée_Totale),
    Aire_Exposée_Totale = replace_na(Aire_Exposée_Totale, 0),
    Part_Territoire_Exposé = Aire_Exposée_Totale / Aire_Totale
  ) %>%
  select(DCOE_C_COD, Part_Territoire_Exposé)

## Comparaison de l'exposition des communes au risque d'exposition le revenu médian et le taux de pauvreté à 60 % du revenu médian ##
# Harmoniser l'appellation des colonnes des données sur les indicateurs socio-économiques et des risques d'inondation
Risque_Commune <- Risque_Commune %>% rename(GEO = DCOE_C_COD)

# Extraire les données de chaque indicateur
Revenu_Median <- Caracteristiques_Commune %>% filter(TAB_MEASURE == "MED_SL") %>% select(GEO, Revenu_Median = OBS_VALUE)
Pauvrete <- Caracteristiques_Commune %>% filter(TAB_MEASURE == "PR_MD60") %>% select(GEO, Pauvrete = OBS_VALUE)

# Fusionner toutes les données avec Risque_Commune
Vulnerabilite <- Risque_Commune %>%
 left_join(Revenu_Median, by = "GEO") %>%
 left_join(Pauvrete, by = "GEO") %>%
 left_join(Chomage, by = "GEO") %>%
 left_join(Menages_Monoparentaux, by = "GEO")

# Modèles de régression
Revenu médian ~ Part_Territoire_Exposé
modele_revenu <- lm(Revenu_Median ~ Part_Territoire_Exposé, data = Vulnerabilite, na.action = na.exclude)
summary(modele_revenu)

# Pauvreté ~ Part_Territoire_Exposé
modele_pauvrete <- lm(Pauvrete ~ Part_Territoire_Exposé, data = Vulnerabilite, na.action = na.exclude)
summary(modele_pauvrete)


## Comparaison du taux de pauvreté entre les communes exposées au risque d'inondation et celles qui en restent épargnées. Comparant deux groupes (part du territoire exposé < 10 % et > 90 %) ##
Communes_Non_Exposées <- Vulnerabilite %>%
 filter(Part_Territoire_Exposé < 0.1) %>%
 pull(Pauvrete)

Communes_Exposées <- Vulnerabilite %>%
 filter(Part_Territoire_Exposé > 0.9) %>%
 pull(Pauvrete)
Moyenne_Non_Exposées <- mean(Communes_Non_Exposées, na.rm = TRUE)
Moyenne_Exposées <- mean(Communes_Exposées, na.rm = TRUE)
print(Moyenne_Non_Exposées)
print(Moyenne_Exposées)

# Test de comparaison (test t de Student pour vérifier si la différence est significative)
resultat_test <- t.test(Communes_Non_Exposées, Communes_Exposées)
print(resultat_test)

### Comparaison avec les quartiers prioritaires de la politique de la ville ##

QPV <- st_read("Données/QPV/QP2024_France_Hexagonale_Outre_Mer_WGS84.gpkg")
QPV <- QPV %>% mutate(the_geom = st_make_valid(the_geom))

# Alignement de la projection des données géographiques
QPV <- st_transform(QPV, 2154)
TRI <- st_transform(TRI, 2154)

QPV <- QPV %>%
  mutate(Superficie_QPV = as.numeric(st_area(.) / 1000000))
# Calculer, pour chaque QPV, la part de la superficie en TRI
intersections_qpv_tri <- st_intersection(QPV, TRI) %>%
  mutate(Aire_Exposée = as.numeric(st_area(.) / 1000000))  # en km²

# Agrégation par QPV (en utilisant code_qp)
QPV_inondable <- intersections_qpv_tri %>%
  group_by(code_qp) %>%
  summarise(Aire_Exposée_Totale = sum(Aire_Exposée, na.rm = TRUE)) %>%
  ungroup() %>%
  st_drop_geometry()

# Fusionner et calculer la part exposée
QPV_Exposition <- QPV %>%
  select(code_qp, Superficie_QPV) %>%
  left_join(QPV_inondable, by = "code_qp") %>%
  mutate(
    Aire_Exposée_Totale = as.numeric(Aire_Exposée_Totale),
    Aire_Exposée_Totale = replace_na(Aire_Exposée_Totale, 0),
    Part_QPV_Exposé = Aire_Exposée_Totale / Superficie_QPV
  ) %>%
  select(code_qp, Part_QPV_Exposé)