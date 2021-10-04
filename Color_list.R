##############################
#
# Color list
#
# Written by Beáta Szeitz
#
##############################


##############################
# This script is for specifying the column annotation colors on the
# heatmaps.
# 
# colorlist will be used as an input for ComplexHeatmap::HeatmapAnnotation()
#
# Packages required: none
#
colorlist <- list(
  "Type of Samples" = 
    c("Prim"="white",
      "Met"="black"),
  "Organ of the Samples" = 
    c("Cut"="white",
      "LN"="black"),
  "Tumor content (%).Cat" = 
    c("[0,10]"="white",
      "(10,25]"="#FED976",
      "(25,55]" ="#FD8D3C",
      "(55,100]"="#E31A1C"),
  "Therapy groups" = 
    c("No therapy"="white",
      "Immunotherapy"="red",
      "Other therapy"="black"),
  "Age at sample collectin.Cat"=
    c("[38,60]"="white",
      "(60,65]" ="#D9F0A3",
      "(65,70]"="#78C679",
      "(70,91]"="#238443"),
  "Sex" = 
    c("male"="white",
      "female"="black"),
  "DFS (m).Cat" = 
    c("[0,10]"="white",
      "(10,30]" = "#C7E9B4",
      "(30,60]"="#41B6C4",
      "(60,150]"="#225EA8"),
  "PFS (m).Cat" = 
    c("[0,10]"="white",
      "(10,30]" = "#C7E9B4",
      "(30,60]"="#41B6C4",
      "(60,100]"="#225EA8",
      "(100,205]" = "black"),
  "OS (m).Cat" = 
    c("[0,10]"="white",
      "(10,30]" = "#C7E9B4",
      "(30,60]"="#41B6C4",
      "(60,100]"="#225EA8",
      "(100,205]" = "black"),
  "Live" = 
    c("0"="white",
      "1"="black"),
  "Loc_p" = 
    c("0"="white",
      "1" = "#B3CDE3",
      "2"="#8C96C6",
      "3"="#8856A7",
      "4" = "#810F7C",
      "6" = "black"),
  "AJCC8" = 
    c("1A"="white",
      "1B" = "#1B9E77",
      "2A"="#D95F02",
      "2B"="#7570B3",
      "3A" = "#E7298A",
      "3B" = "#66A61E",
      "3C"="#E6AB02",
      "3D" = "#A6761D",
      "4" = "black"),
  "Type.Simplified" = 
    c("ALM"="white",
      "ALM with SSM" = "#1B9E77",
      "Blue naevus with vertical growth"="#D95F02",
      "LMM"="#7570B3",
      "MUP" = "#E7298A",
      "NM" = "#66A61E",
      "Spitzoid melanoma"="#E6AB02",
      "SSM" = "#A6761D",
      "SSM with vertical growth" = "#3d087b",
      "Uveal melanoma" = "black"),
  "T" = 
    c("1"="white",
      "2"="#D9F0A3",
      "3" ="#78C679",
      "4"="#238443"),
  "U" = 
    c("0"="white",
      "1"="black"),
  "Breslow (mm).Cat" = 
    c("[0,2]"="white",
      "(2,4]"="#D9F0A3",
      "(4,8]" ="#78C679",
      "(8,31]"="#238443"),
  "Regres" = 
    c("0"="white",
      "1"="black"),
  "BRAFstate_simple" = 
    c("NO"="white",
      "NO (cKIT)"="black",
      "BRAF"="#865c47"),
  "BRAFstate" = 
    c("NO"="white",
      "NO (cKIT)"="black",
      "BRAF D587G"="#c1cd77",
      "BRAFV600E"="#9e1b27",
      "BRAFV600K"="#865c47"),
  "MSOrder.Cat" = 
    c("[1,24]"="white",
      "(24,46]"="#ffe3d6",
      "(46,68]" ="#e300c3",
      "(68,90]"="#a70098"),
  "Year.of.collection.Cat" = 
    c("[2005,2015]"="white",
      "(2015,2017]" = "#ffe3d6",
      "(2017,2018]"="#e300c3",
      "(2018,2020]"="#a70098"),
  "No.MVs.Cat" = 
    c("[300,400]"="white",
      "(400,450]"="#ffe3d6",
      "(450,550]"="#e300c3",
      "(550,750]" ="#a70098",
      "(750,1100]"="#3f004b"),
  "Stage" = 
    c("1"="white",
      "2"="#ffe3d6",
      "3" ="#e300c3",
      "4"="#a70098"),
  "Progression.time.months.Cat" = 
    c("[0,3]"="white",
      "(3,7]"="#FED976",
      "(7,12]" ="#FD8D3C",
      "(12,60]"="#E31A1C"),
  "Progression.event" = 
    c("0"="white",
      "1"="black"),
  "Progression.therapyType" = 
    c("Target"="white",
      "Immun"="black"),
  "SampleCluster" = 
    c("0"="white",
      "1" = "#1B9E77",
      "2"="#D95F02",
      "3"="#7570B3",
      "4" = "#E7298A",
      "5" = "#66A61E",
      "6"="black")
)