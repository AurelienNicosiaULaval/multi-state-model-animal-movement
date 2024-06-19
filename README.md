# Code R et Ensemble de Données pour la Thèse "Modèles Multi-États Généraux pour l’Analyse du Mouvement Animalier"

Bienvenue dans le dépôt GitHub associé à ma thèse intitulée "Modèles multi-états généraux pour l’analyse du mouvement animalier". Ce dépôt contient tout le code R et les ensembles de données utilisés dans mes chapitres.

## Structure du Dépôt

Le dépôt est structuré de la manière suivante, pour chapitre Chapitre (1 à 4) :

- **/data** : Contient tous les ensembles de données utilisés dans les analyses.
- **/scripts** : Contient les scripts R pour les analyses effectuées dans le chapitre.

## Utilisation

### Prérequis

Pour reproduire les analyses, vous aurez besoin des packages R suivants :

- `tidyverse`
- `geosphere`
- `ggplot2`
- `dplyr`
- `lubridate`
- `survival`
- `...` (ajoutez tous les autres packages nécessaires)

### Installation des Packages

Vous pouvez installer les packages nécessaires en exécutant le script `install_packages.R` :

```r
source("scripts/install_packages.R")
