excel_sheets("Input_Data/week1/2013_2019_CALFIRE_Redbook.xlsx")
calfire.metadata <- read_excel("Input_Data/week1/2013_2019_CALFIRE_Redbook.xlsx", sheet = "Metadata")
calfire.data <- read_excel("Input_Data/week1/2013_2019_CALFIRE_Redbook.xlsx", sheet = 2)

socal.fires <- calfire.data %>%
select(County_Unit:Controlled_Date, Total_Acres_Burned:Civil_Fatalities) %>% 
filter(County_Unit %in% c("SANTA BARBARA", "VENTURA", "LOS ANGELES", "ORANGE", "SAN DIEGO") & Total_Acres_Burned >= 500 | Fire_Name == "THOMAS") %>% 
arrange(desc(Start_Date), Total_Acres_Burned) %>% 
mutate_at(vars(Structures_Destroyed:Civil_Fatalities), replace_na, 0) %>% 
mutate(Fatalities = Fire_Fatalities + Civil_Fatalities) %>%
mutate(interv = 
           interval(Start_Date, Controlled_Date), 
         dur = as.duration(interv), 
         days = as.numeric(dur, "days"), 
         County_Unit = ifelse(County_Unit == "VENTURA/SANTA BARBARA", "VENTURA", County_Unit))

ggplot(socal.fires, aes(x = Start_Date, y = Total_Acres_Burned)) + 
  geom_point(aes(color = County_Unit)) + 
  ggtitle("CA South Coast Major Fires \n2014 - 2018") +
  labs(x = "", y = "Total Acres Burned", color = "County") + theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  facet_grid(rows = "County_Unit", scales = "free")

plot.data <- socal.fires %>% 
  rename(county = County_Unit,
         acres = Total_Acres_Burned,
         start = Start_Date,
         end = Controlled_Date) %>% 
  mutate(year = year(start), 
         county = ifelse(county == "VENTURA/SANTA BARBARA", 
                         "VENTURA", county))

incidents <- plot.data %>% 
  group_by(county, year) %>% 
  tally() %>% 
  ungroup()
  
incidents.plot <- incidents %>% 
  ggplot(aes(x = year, y = n)) +
  geom_point(aes(color = county)) +
  geom_line(aes(color = county)) +
  labs(title = "CA South Coast Major Fire Incidents \n 2013 - 2018",
       x = "", y = "Incidents", color = "County") +
  theme_bw() +
  facet_grid(rows = "county", scales = "free")

