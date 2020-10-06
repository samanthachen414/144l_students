# 144l-students
# Samantha Chen (04/10/2020)
library(tidyverse)
library(readxl)
excel_sheets("Input_Data/week1/2013_2019_CALFIRE_Redbook.xlsx")
calfire.metadata <- read_excel("Input_Data/week1/2013_2019_CALFIRE_Redbook.xlsx", sheet = "Metadata")
calfire.data <- read_excel("Input_Data/week1/2013_2019_CALFIRE_Redbook.xlsx", sheet = 2)
#####

This repository includes data, code, and supplementary information to be used by students for the remotely instructed EEMB144L during the fall quarter 2020 (@) UCSB. 