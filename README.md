# 144l-students
# Samantha Chen (04/10/2020)

library(tidyverse)

library(readxl)

excel_sheets("Input_Data/week1/2013_2019_CALFIRE_Redbook.xlsx")

calfire.metadata <- read_excel("Input_Data/week1/2013_2019_CALFIRE_Redbook.xlsx", sheet = "Metadata")

calfire.data <- read_excel("Input_Data/week1/2013_2019_CALFIRE_Redbook.xlsx", sheet = 2)

##### Initial Exploration #####
Samantha Chen (10/10/2020)

names(calfire.data) #shows variable (column) names

dim(calfire.data) #dimensions of data set (row, columns)

class(calfire.data) #shows data class

head(calfire.data) #shows first six observations

tail(calfire.data)  #shows last six observations

# ?? brings up every package that might have a certain function

#Single columns can be referred to using a "$"
county <- calfire.data$County_Unit

max_acres <- max(calfire.data$Total_Acres_Burned, na.rm = T) #na.rm = T ignores any "NA"

max(calfire.data$Structures_Destroyed)
max(calfire.data$Structures_Destroyed, na.rm = T)

##### Basic Data Wrangling (dplyr functions) #####

df1 <- select(calfire.data, County_Unit:Controlled_Date, Total_Acres_Burned:Civil_Fatalities) #Pulls out variables

view(df1)

df2 <- filter(df1, County_Unit %in% c("SANTA BARBARA", "VENTURA", "LOS ANGELES", "ORANGE", "SAN DIEGO") & Total_Acres_Burned >= 500 | Fire_Name == "THOMAS") #Pulls out rows

View(df2)

df3 <- arrange(df2, desc(Start_Date), Total_Acres_Burned) #Can order data of choosing

df4 <- mutate_at(df3, vars(Structures_Destroyed:Civil_Fatalities), replace_na, 0) #replaces NA with something else
  
df5 <- mutate(df4, Fatalities = Fire_Fatalities + Civil_Fatalities) #Can create new rows like totals

#Mess with Time
library(lubridate)

df6 <- mutate(df5, 
            interv = interval(Start_Date, Controlled_Date), 
            dur = as.duration(interv), 
            days = as.numeric(dur, "days")
        
# We used 15 lines to do all of that. Now we have 5 dataframes which is a little inefficient. So, the better way to do this is called "piping"

#### Intro to Piping ####

#We want to restrict our data to the Socal coast, exclude fires that burned less than 500 acres, add column that sums the number of fatalities, change NA's to 0's, arrange data

#Pipe Operator: %>% (control + shift + m)

socal.fires <- calfire.data %>%
  select(calfire.data, County_Unit:Controlled_Date, Total_Acres_Burned:Civil_Fatalities) %>% 
  filter(County_Unit %in% c("SANTA BARBARA", "VENTURA", "LOS ANGELES", "ORANGE", "SAN DIEGO") & Total_Acres_Burned >= 500 | Fire_Name == "THOMAS") %>% 
  arrange(desc(Start_Date), Total_Acres_Burned) %>% 
  mutate_at(vars(Structures_Destroyed:Civil_Fatalities), replace_na, 0) %>% 
  mutate(Fatalities = Fire_Fatalities + Civil_Fatalities) %>%   mutate(interv = interval(Start_Date, Controlled_Date), dur = as.duration(interv), days = as.numeric(dur, "days"), County_Unit = ifelse(County_Unit == "VENTURA/SANTA BARBARA", "VENTURA", County_Unit)

#### Our First graphs in ggplot ####

# Three things you must tell R to make a graph in ggplot
# (1) That you're using ggplot
# (2) What data you're using (including whay should be x and y)
# (3) Type of graph that you want to create

socal.plot <- socal_fires %>%
  rename(start = Start_Date,
          acres = Total_Acres_Burned) %>% 
  ggplot(aes(x = start, y = acres)) + geom_point(aes(color = county)) +
  ggtitle("California SoCal Major Fires 2013 - 2018") +
    xlab("Date") +
    ylab("Acres Burned") +
    theme(panel.grid.major = element_blank())
    
# Separate into multiple plots

socal.plot + facet_grid(~county)
  
  


This repository includes data, code, and supplementary information to be used by students for the remotely instructed EEMB144L during the fall quarter 2020 (@) UCSB. 

