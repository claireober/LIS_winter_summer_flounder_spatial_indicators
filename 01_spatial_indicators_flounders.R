
# ============================================================
# LIS Spatial Indicator Analysis (CT DEEP trawl survey)
# Supporting code for: "Evaluating spatial and temporal dynamics of cold- and warm-adapted fish in a changing Long Island Sound"
#
# NOTE ON DATA ACCESS:
# Raw CT DEEP data are not publicly available. This script assumes
# the user has obtained permission and has local access to the files.
# Abundance index estimates were also provided by CT DEEP for these analyes 
# ============================================================

library(mapdata)
library(maps)
library(RGeostats)
library(dplyr)
library(rstatix)
library(ggplot2)
library(gridExtra)
library(grid)
library(ggpubr)
library(viridis)
library(rworldxtra); data("countriesHigh")

### DATA PREP ###

# OPTIONAL: length-based filtering
# Not used in the final manuscript analyses
# If your species is using length-specific data, use the following code and modify the length cutoffs for your specific species

Station_Info=read.csv("C:/Users/oberc/OneDrive/Desktop/FVCOM_LIS/BC_LOB/Tow_data.csv")
Num_Len=read.csv("C:/Users/oberc/OneDrive/Desktop/FVCOM_LIS/BC_LOB/Fish_len.csv")
Num_Len=Num_Len[Num_Len$ExpCatchNum > 0 ,]

# Make a data set that contains the total number of fish caught at each length at each station
Num_Len_agg=aggregate(ExpCatchNum ~ Sample_Number + Year + Season + Common_name + Length_mm, Num_Len, sum)
Num_Len_agg$Length_mm=as.numeric(Num_Len_agg$Length_mm)

df_len=Num_Len_agg[Num_Len_agg$Common_name== "summer flounder",] # enter the species name exactly as it appears in the Num_Len_agg file

# Remove length frequencies that are greater than, or smaller than the length cutoff, but always include the cutoff value
df_len=df_len[df_len$Length_mm >= 300 ,] # The value will not be in quotations


# Count up the total number of fish at each station that is larger or smaller than our length cutoff
df_len_total=aggregate(ExpCatchNum ~ Sample_Number + Year + Season, df_len, sum)
names(df_len_total)[4] = "SUMMERFL_NUM"

# Associate the number caught with the station information
Station_Info_catch=left_join(Station_Info,df_len_total)

# Replace NAs with 0 only during years when length data is available

catch_new = Station_Info_catch %>% 
  mutate(across(SUMMERFL_NUM, ~ replace(., is.na(.), 0)))

# Remove the years where there is no data
catch_new=catch_new[!is.na(catch_new$SUMMERFL_NUM),]



### FORLOOP ###

start_y=40.8
end_y=41.4
start_x=-73.72
end_x=-72

start_y_nm=-20
end_y_nm=15
start_x_nm=-50
end_x_nm=40

polygon_df_LIS=read.csv("C:/Users/oberc/OneDrive/Desktop/FVCOM_LIS/BC_LOB/polygon_df_LIS.csv")
polygon_df_LIS=polygon_df_LIS[,-c(1)]
poly.data= polygon.create(polygon_df_LIS)

camargo_eveness=function(n_spec,include_zeros=TRUE){
  if(is.vector(n_spec)==FALSE){stop("\n n_spec must be a vector of abundance of species \n")}  
  if(include_zeros){n<-n_spec}else{n<-n_spec[n_spec>0]}
  S<-length(n)
  camar=matrix(nrow=length(n), ncol=length(n))
  for(i in 1:S){
    for(j in 1:S){
      p_i=n[i]/sum(n)
      p_j=n[j]/sum(n)
      camar[i,j]=((abs(p_i-p_j))/S)
    }
  }
  sum.camar=abs(sum(as.dist(camar,diag=FALSE,upper=FALSE)))
  return(1-sum.camar)
}


catch_new$Year=as.character(catch_new$Year)
catch_new= catch_new[order(catch_new$Year, decreasing = FALSE),]


# Diagnostic section where we can identify if there are any years and seasons where there was zero catch

catch_new_diag=aggregate(SUMMERFL_NUM ~ Year + Season, data=catch_new, sum)


season=c("FA", "SP")[1:2]
yr_u=unique(catch_new$Year)
yr=c(unique(catch_new$Year))[1:length(yr_u)]

model_setting=expand.grid(yr, season)
colnames(model_setting)=c("Year", "Season")

model_setting = model_setting %>%
  filter(!(Year=="2010" & Season=="FA"))

print(model_setting)


res.table= cbind(model_setting, CG_Lon=NA, CG_Lat=NA, Inertia=NA, Isotropy=NA, AI=NA, Pos_Area=NA, Equiv_Area=NA, Spread_Area=NA, Evenness=NA)


for (i in 1:nrow(model_setting)){
  
  temp_setting=model_setting[i,]
  Year=temp_setting$Year
  Season=temp_setting$Season
  
  
  ##Data prep, extracting the season, year data for each iteration
  Data_Set=catch_new[catch_new$Season== Season,]
  Data_Set=Data_Set[Data_Set$Year== Year,]
  
  
  ##Creating our db object
  #This should represent the longitude, latitude and abundance of your specific species. Column 3 is SUMMERFL
  Data_Set=Data_Set[,c("Longitude", "Latitude", "SUMMERFL_NUM")]
  colnames(Data_Set)[3]= "CATCH_NUM"
  
  
  #Add a title that changes with the iterations, add map, add a way to automatically save it to a specific file. This is making the map using lat and long coords
  projec.toggle(0)
  db.data=new("db", flag.grid=FALSE, ndim=2, x0=0, dx=0, nx=0, angles=0, flag.rotation=0, locators= c("x1",   "x2",   "z1"), items=Data_Set)
  
  db.data <- db.delete(db=db.data,names=6)
  db.data <- infl(db.data,nodes=c(400,400),origin=c(-74,40.5),extend=c(2,1),
                  dmax=100,polygon=poly.data,plot=T,asp=1)
  
  
  
  mypath <- file.path("C:/Users/oberc/OneDrive/Desktop/Flounders/Summer/",paste(temp_setting$Year, temp_setting$Season,"- SUMMERFL", ".jpg", sep = ""))
  
  jpeg(file=mypath)
  
  mytitle = paste("SUMMERFL Density, CG and Inertia", temp_setting$Year, "-", temp_setting$Season)
  
  plot(db.data,title=mytitle ,ylim=c(start_y, end_y), xlim=c(start_x, end_x),
       asp=1/cos(mean(db.extract(db=db.data,names="x2"))*pi/180),inches=4, pch = ifelse(db.data@items$CATCH_NUM > 0, 1, 4))
  
  legend("bottomright",legend=c(0 , max(db.data@items$CATCH_NUM)), pch = c(4,1),
         col  = c("blue", "blue"), pt.cex = c(0.5,4), cex=1.1)
  
  # Center of gravity, Inertia and Isotropy  
  
  cgi = SI.cgi(db.data,flag.plot=T,flag.inertia=T,col=2)
  map(database = "worldHires", ylim=c(start_y, end_y), xlim=c(start_x, end_x),col = "gray90", fill = TRUE, add = TRUE, lwd=0.1)
  box()
  dev.off() 
  
  
  ## Now switch to nautical miles 
  
  projec.define(projection="mean",db=db.data)
  db.data=new("db", flag.grid=FALSE, ndim=2, x0=0, dx=0, nx=0, angles=0, flag.rotation=0, locators= c("x1",   "x2",   "z1"), items=Data_Set)
  
  db.data <- db.delete(db=db.data,names=6)
  db.data <- infl(db.data,nodes=c(400,400),origin=c(-74,40.5),extend=c(2,1),
                  dmax=100,polygon=poly.data,plot=T,asp=1)
  
  
  
  plot(db.data,title=paste0("SUMMERFL Density, CG and Inertia - ", temp_setting$Year, "-", temp_setting$Season) ,ylim=c(start_y_nm, end_y_nm), xlim=c(start_x_nm, end_x_nm),
       asp=1/cos(mean(db.extract(db=db.data,names="x2"))*pi/180),inches=4, pch = ifelse(db.data@items$CATCH_NUM > 0, 1, 4))
  
  
  legend("bottomright",legend=c(0 , max(db.data@items$CATCH_NUM)), pch = c(4,1),
         col  = c("blue", "blue"), pt.cex = c(0.5,4), cex=1.1)
  
  
  cgi = SI.cgi(db.data,flag.plot=T,flag.inertia=T,col=2)
  plot(poly.data,col=NA,add=TRUE)
  
  
  #Calculate Results and put them in a table
  res = projec.invert(cgi$center[1],cgi$center[2])
  stats = SI.stats(db.data,flag.plot=F)
  
  res.table[i,3]=res$x
  res.table[i,4]=res$y
  res.table[i,5]=cgi$inertia
  res.table[i,6]=cgi$iso
  res.table[i,7]=stats$totab
  res.table[i,8]=stats$parea
  res.table[i,9]=stats$eqarea
  res.table[i,10]=stats$sparea
  res.table[i,11]=camargo_eveness(Data_Set$CATCH_NUM, FALSE)
  
}


write.csv(res.table, "C:/Users/oberc/OneDrive/Desktop/Flounders/Summer/res_table.csv")



res.table=read.csv("C:/Users/oberc/OneDrive/Desktop/Flounders/Summer/res_table.csv")
abundance_ind=read.csv("C:/Users/oberc/OneDrive/Desktop/LIS_HSI/Data/CTDEEP_AI_All_Lenth.csv")
ai_SUMMERFL=abundance_ind[abundance_ind$Species=="SFL",]
ai_SUMMERFL=ai_SUMMERFL[ai_SUMMERFL$Year<2022,]
res_all=left_join(res.table, ai_SUMMERFL)
res.table=res_all

### Plots


##Spring

res_table_S=res.table[res.table$Season=="SP",]
res_table_S$Year=as.numeric(as.character(res_table_S$Year))

even=ggplot(res_table_S, aes(x=Year, y=Evenness)) + 
  geom_line() +
  labs(title = "Evenness", y = "Evenness", x ="Year")   + 
  theme(panel.background = element_blank()) +
  geom_smooth(method = 'loess',se = TRUE)


lon=ggplot(res_table_S, aes(x=Year, y=CG_Lon)) + 
  geom_line() +
  labs(title = "CG Longitude", y = "CG Lon", x ="Year")   + 
  theme(panel.background = element_blank()) +
  geom_smooth(method = 'loess',se = TRUE)

lat=ggplot(res_table_S, aes(x=Year, y=CG_Lat)) + 
  geom_line() +
  labs(title = "CG Latitiude", y = "CG Lat", x ="Year")   + 
  theme(panel.background = element_blank()) +
  geom_smooth(method = 'loess',se = TRUE)

inertia=ggplot(res_table_S, aes(x=Year, y=Inertia)) + 
  geom_line() +
  labs(title = "Inertia", y = "Inertia (nautical miles^2)", x ="Year")   + 
  theme(panel.background = element_blank()) +
  geom_smooth(method = 'loess',se = TRUE)

ai=ggplot(res_table_S, aes(x=Year, y=GmMn)) +
  geom_line() +
  labs(title = "Abundance Index", y = "Abundance Index", x ="Year")   +
  theme(panel.background = element_blank()) +
  geom_smooth(method = 'loess',se = TRUE)

pos_area=ggplot(res_table_S, aes(x=Year, y=Pos_Area)) + 
  geom_line() +
  labs(title = "Positive Area", y = "Positive Area (nautical miles^2)", x ="Year")   + 
  theme(panel.background = element_blank()) +
  geom_smooth(method = 'loess',se = TRUE)

plot_S=arrangeGrob(top=textGrob("Summer Flounder - Spring Spatial Indicators", gp=gpar(fontsize=20)),lon, lat, inertia, ai, pos_area, even, nrow = 2)

ggsave("C:/Users/oberc/OneDrive/Desktop/Flounders/Summer/SUMMERFL_S_plot.png", plot_S, width = 30, height = 20, units = "cm")






## Fall

res_table_F=res.table[res.table$Season=="FA",]
res_table_F$Year=as.numeric(as.character(res_table_F$Year))

even=ggplot(res_table_F, aes(x=Year, y=Evenness)) + 
  geom_line() +
  labs(title = "Evenness", y = "Evenness", x ="Year")   + 
  theme(panel.background = element_blank()) +
  geom_smooth(method = 'loess',se = TRUE)


lon=ggplot(res_table_F, aes(x=Year, y=CG_Lon)) + 
  geom_line() +
  labs(title = "CG Longitude", y = "CG Lon", x ="Year")   + 
  theme(panel.background = element_blank()) +
  geom_smooth(method = 'loess',se = TRUE)

lat=ggplot(res_table_F, aes(x=Year, y=CG_Lat)) + 
  geom_line() +
  labs(title = "CG Latitiude", y = "CG Lat", x ="Year")   + 
  theme(panel.background = element_blank()) +
  geom_smooth(method = 'loess',se = TRUE)

inertia=ggplot(res_table_F, aes(x=Year, y=Inertia)) + 
  geom_line() +
  labs(title = "Inertia", y = "Inertia (nautical miles^2)", x ="Year")   + 
  theme(panel.background = element_blank()) +
  geom_smooth(method = 'loess',se = TRUE)

ai=ggplot(res_table_F, aes(x=Year, y=GmMn)) +
  geom_line() +
  labs(title = "Abundance Index", y = "Abundance Index", x ="Year")   +
  theme(panel.background = element_blank()) +
  geom_smooth(method = 'loess',se = TRUE)

pos_area=ggplot(res_table_F, aes(x=Year, y=Pos_Area)) + 
  geom_line() +
  labs(title = "Positive Area", y = "Positive Area (nautical miles^2)", x ="Year")   + 
  theme(panel.background = element_blank()) +
  geom_smooth(method = 'loess',se = TRUE)

plot_F=arrangeGrob(top=textGrob("Summer Flounder- Fall Spatial Indicators", gp=gpar(fontsize=20)),lon, lat, inertia, ai, pos_area, even, nrow = 2)

ggsave("C:/Users/oberc/OneDrive/Desktop/Flounders/Summer/SUMMERFL_F_plot.png", plot_F, width = 30, height = 20, units = "cm")




## Center of Gravity Maps


res.table$Season=recode(res.table$Season, "SP" = "Spring", "FA"="Fall")
res.table$Year=as.numeric(as.character(res.table$Year))

##CG plot
#Changing ggplots 
theme_update(plot.title = element_text(hjust = 0.5))

CG_all=ggplot(res.table, aes(CG_Lon, CG_Lat))+
  geom_point(aes(color = Year), size=2) +
  scale_color_viridis(discrete= FALSE, direction= -1, option = "mako") +
  ggtitle("SUMMERFL Center of Gravity") +
  theme_bw()+theme(strip.text = element_text(size=10, face="bold"))+
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(strip.background = element_rect(fill="white"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  coord_quickmap(xlim = c(start_x, end_x), ylim=c(start_y, end_y)) +
  borders(database = "usa", regions = ".",fill = "gray90",xlim = NULL,ylim = NULL)+
  facet_wrap(~factor(Season, c("Spring", "Fall")), nrow=1, ncol=2)+
  labs(y = "Latitude", x ="Longitude")



##Save plot as a pdf
ggsave("C:/Users/oberc/OneDrive/Desktop/Flounders/Summer/SUMMERFL_CG_plot.png", CG_all, width = 30, height = 10, units = "cm") 



############ Nice plot for flounders paper. Combined spring and fall spatial indicator plots




library(ggplot2)
library(gridExtra)
library(viridis)
library(dplyr)


res.table$Year <- as.numeric(as.character(res.table$Year))

# Split the data by season
res_table_S <- res.table[res.table$Season == "SP", ]
res_table_F <- res.table[res.table$Season == "FA", ]


res.table$Season <- recode(res.table$Season, "SP" = "Spring", "FA" = "Fall")

#  indicator labels
indicator_titles <- c(
  "Evenness" = "F) Evenness",
  "CG_Lon" = "A) Center of Gravity Longitude",
  "CG_Lat" = "B) Center of Gravity Latitude",
  "Inertia" = "C) Inertia (nautical miles²)",
  "GmMn" = "D) Abundance Index",
  "Pos_Area" = "E) Positive Area (nautical miles²)"
)

# Function to create combined plots
create_combined_plot <- function(data, y_var, label, y_label) {
  ggplot(data, aes(x = Year, y = !!sym(y_var), color = Season)) +
    geom_line(size = 1.2) +  # Thicker lines
    geom_smooth(method = 'loess', se = TRUE, size = 1.2) +  # Thicker smooth line
    scale_color_manual(values = c("Spring" = "blue", "Fall" = "red")) +
    labs(title = label, y = y_label, x = "Year") +
    theme_minimal() +
    theme(panel.background = element_blank(),
          plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
}

# Create combined plots for each indicator with 
even_plot <- create_combined_plot(res.table, "Evenness", indicator_titles["Evenness"], "Evenness")
lon_plot <- create_combined_plot(res.table, "CG_Lon", indicator_titles["CG_Lon"], "CG Longitude")
lat_plot <- create_combined_plot(res.table, "CG_Lat", indicator_titles["CG_Lat"], "CG Latitude")
inertia_plot <- create_combined_plot(res.table, "Inertia", indicator_titles["Inertia"], "Inertia (nautical miles²)")
ai_plot <- create_combined_plot(res.table, "GmMn", indicator_titles["GmMn"], "Abundance Index")
pos_area_plot <- create_combined_plot(res.table, "Pos_Area", indicator_titles["Pos_Area"], "Positive Area (nautical miles²)")


combined_plot <- arrangeGrob(
  top = textGrob("Summer Flounder Spring and Fall Spatial Indicators", gp = gpar(fontsize = 20, fontface = "bold")),
  lon_plot, lat_plot, inertia_plot, ai_plot, pos_area_plot, even_plot,
  nrow = 2
)

# Save
ggsave("C:/Users/oberc/OneDrive/Desktop/Flounders/Summer/SUMMERFL_Combined_Spring_Fall_Plot2.png",
       combined_plot, width = 30, height = 20, units = "cm")



######################### new COG maps for paper 


res.table$Season <- recode(res.table$Season, "SP" = "Spring", "FA" = "Fall")
res.table$Year <- as.numeric(as.character(res.table$Year))


res.table <- res.table %>%
  group_by(Season) %>%
  mutate(
    Symbol = case_when(
      Year == min(Year) ~ "First Year",
      Year == max(Year) ~ "Last Year",
      TRUE ~ "Intermediate Year"
    )
  )


symbol_shapes <- c("First Year" = 17, "Intermediate Year" = 16, "Last Year" = 15)
symbol_colors <- c("First Year" = "red", "Intermediate Year" = "black", "Last Year" = "red")

# Updated COG plot
CG_all <- ggplot(res.table, aes(CG_Lon, CG_Lat)) +
  geom_point(
    data = res.table %>% filter(Symbol == "Intermediate Year"),
    aes(color = Year, shape = Symbol),
    size = 2, show.legend = TRUE
  ) +
  # First and Last years: fixed color and give it with distinct shapes
  geom_point(
    data = res.table %>% filter(Symbol %in% c("First Year", "Last Year")),
    aes(shape = Symbol),
    color = "red", size = 3, show.legend = TRUE
  ) +
  scale_color_viridis(discrete = FALSE, direction = -1, option = "mako", name = "Year") +
  scale_shape_manual(values = symbol_shapes, name = "Year Type") +
  guides(shape = guide_legend(override.aes = list(color = c("red", "black", "red")))) +
  ggtitle("Summer Flounder Center of Gravity") +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    strip.text = element_text(size = 10, face = "bold"),
    strip.background = element_rect(fill = "white"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  coord_quickmap(xlim = c(start_x, end_x), ylim = c(start_y, end_y)) +
  borders(database = "usa", regions = ".", fill = "gray90", xlim = NULL, ylim = NULL) +
  facet_wrap(~factor(Season, c("Spring", "Fall")), nrow = 1, ncol = 2) +
  labs(y = "Latitude", x = "Longitude")


plot(CG_all)

# Save 
ggsave("C:/Users/oberc/OneDrive/Desktop/Flounders/Summer/SUMMERFL_CG_plot_updated.png", 
       CG_all, width = 30, height = 10, units = "cm")



##########################################################################################################################################

##########################################################################################################################################

############### winter flounder below ####################################################################################################

##########################################################################################################################################


### DATA PREP ###


# OPTIONAL: length-based filtering
# Not used in the final manuscript analyses
# If your species is using length-specific data, use the following code and modify the length cutoffs for your specific species

Station_Info=read.csv("C:/Users/oberc/OneDrive/Desktop/FVCOM_LIS/BC_LOB/Tow_data.csv")
Num_Len=read.csv("C:/Users/Oberc/OneDrive/Desktop/FVCOM_LIS/BC_LOB/Fish_len.csv")
Num_Len=Num_Len[Num_Len$ExpCatchNum > 0 ,]

# Make a data set that contains the total number of fish caught at each length at each station
Num_Len_agg=aggregate(ExpCatchNum ~ Sample_Number + Year + Season + Common_name + Length_mm, Num_Len, sum)
Num_Len_agg$Length_mm=as.numeric(Num_Len_agg$Length_mm)

df_len=Num_Len_agg[Num_Len_agg$Common_name== "winter flounder",] # enter the species name exactly as it appears in the Num_Len_agg file

# Remove length frequencies that are greater than, or smaller than the length cutoff, but always include the cutoff value
df_len=df_len[df_len$Length_mm >= 300 ,] # The value will not be in quotations


# Count up the total number of fish at each station that is larger or smaller than our length cutoff
df_len_total=aggregate(ExpCatchNum ~ Sample_Number + Year + Season, df_len, sum)
names(df_len_total)[4] = "WINTERFL_NUM"

# Associate the number caught with the station information
Station_Info_catch=left_join(Station_Info,df_len_total)

# Replace NAs with 0 only during years when length data is available
catch_new = Station_Info_catch %>% 
  mutate(across(WINTERFL_NUM, ~ replace(., is.na(.), 0)))

# Remove the years where there is no data
catch_new=catch_new[!is.na(catch_new$WINTERFL_NUM),]



### FORLOOP ###

start_y=40.8
end_y=41.4
start_x=-73.72
end_x=-72

start_y_nm=-20
end_y_nm=15
start_x_nm=-50
end_x_nm=40

polygon_df_LIS=read.csv("C:/Users/oberc/OneDrive/Desktop/FVCOM_LIS/BC_LOB/polygon_df_LIS.csv")
polygon_df_LIS=polygon_df_LIS[,-c(1)]
poly.data= polygon.create(polygon_df_LIS)

camargo_eveness=function(n_spec,include_zeros=TRUE){
  if(is.vector(n_spec)==FALSE){stop("\n n_spec must be a vector of abundance of species \n")}  
  if(include_zeros){n<-n_spec}else{n<-n_spec[n_spec>0]}
  S<-length(n)
  camar=matrix(nrow=length(n), ncol=length(n))
  for(i in 1:S){
    for(j in 1:S){
      p_i=n[i]/sum(n)
      p_j=n[j]/sum(n)
      camar[i,j]=((abs(p_i-p_j))/S)
    }
  }
  sum.camar=abs(sum(as.dist(camar,diag=FALSE,upper=FALSE)))
  return(1-sum.camar)
}


catch_new$Year=as.character(catch_new$Year)
catch_new= catch_new[order(catch_new$Year, decreasing = FALSE),]


# Diagnostic section where we can identify if there are any years and seasons where there was zero catch

catch_new_diag=aggregate(WINTERFL_NUM ~ Year + Season, data=catch_new, sum)


season=c("FA", "SP")[1:2]
yr_u=unique(catch_new$Year)
yr=c(unique(catch_new$Year))[1:length(yr_u)]

model_setting=expand.grid(yr, season)
colnames(model_setting)=c("Year", "Season")

model_setting = model_setting %>%
  filter(!(Year=="2010" & Season=="FA"))

print(model_setting)


res.table= cbind(model_setting, CG_Lon=NA, CG_Lat=NA, Inertia=NA, Isotropy=NA, AI=NA, Pos_Area=NA, Equiv_Area=NA, Spread_Area=NA, Evenness=NA)


for (i in 1:nrow(model_setting)){
  
  temp_setting=model_setting[i,]
  Year=temp_setting$Year
  Season=temp_setting$Season
  
  
  ##Data prep, extracting the season, year data for each iteration
  Data_Set=catch_new[catch_new$Season== Season,]
  Data_Set=Data_Set[Data_Set$Year== Year,]
  
  
  ##Creating our db object
  #This should represent the longitude, latitude and abundance of your specific species. 
  Data_Set=Data_Set[,c("Longitude", "Latitude", "WINTERFL_NUM")]
  colnames(Data_Set)[3]= "CATCH_NUM"
  
  
  #Add a title that changes with the iterations, add map, add a way to automatically save it to a specific file. This is making the map using lat and long coords
  projec.toggle(0)
  db.data=new("db", flag.grid=FALSE, ndim=2, x0=0, dx=0, nx=0, angles=0, flag.rotation=0, locators= c("x1",   "x2",   "z1"), items=Data_Set)
  
  db.data <- db.delete(db=db.data,names=6)
  db.data <- infl(db.data,nodes=c(400,400),origin=c(-74,40.5),extend=c(2,1),
                  dmax=100,polygon=poly.data,plot=T,asp=1)
  
  
  
  mypath <- file.path("C:/Users/oberc/OneDrive/Desktop/Flounders/Winter/",paste(temp_setting$Year, temp_setting$Season,"- WINTERFL", ".jpg", sep = ""))
  
  jpeg(file=mypath)
  
  mytitle = paste("WINTERFL Density, CG and Inertia", temp_setting$Year, "-", temp_setting$Season)
  
  plot(db.data,title=mytitle ,ylim=c(start_y, end_y), xlim=c(start_x, end_x),
       asp=1/cos(mean(db.extract(db=db.data,names="x2"))*pi/180),inches=4, pch = ifelse(db.data@items$CATCH_NUM > 0, 1, 4))
  
  legend("bottomright",legend=c(0 , max(db.data@items$CATCH_NUM)), pch = c(4,1),
         col  = c("blue", "blue"), pt.cex = c(0.5,4), cex=1.1)
  
  # Center of gravity, Inertia and Isotropy 
  
  cgi = SI.cgi(db.data,flag.plot=T,flag.inertia=T,col=2)
  map(database = "worldHires", ylim=c(start_y, end_y), xlim=c(start_x, end_x),col = "gray90", fill = TRUE, add = TRUE, lwd=0.1)
  box()
  dev.off() 
  
  
  ## Now switch to nautical miles 
  
  projec.define(projection="mean",db=db.data)
  db.data=new("db", flag.grid=FALSE, ndim=2, x0=0, dx=0, nx=0, angles=0, flag.rotation=0, locators= c("x1",   "x2",   "z1"), items=Data_Set)
  
  db.data <- db.delete(db=db.data,names=6)
  db.data <- infl(db.data,nodes=c(400,400),origin=c(-74,40.5),extend=c(2,1),
                  dmax=100,polygon=poly.data,plot=T,asp=1)
  
  
  
  plot(db.data,title=paste0("WINTERFL Density, CG and Inertia - ", temp_setting$Year, "-", temp_setting$Season) ,ylim=c(start_y_nm, end_y_nm), xlim=c(start_x_nm, end_x_nm),
       asp=1/cos(mean(db.extract(db=db.data,names="x2"))*pi/180),inches=4, pch = ifelse(db.data@items$CATCH_NUM > 0, 1, 4))
  
  
  legend("bottomright",legend=c(0 , max(db.data@items$CATCH_NUM)), pch = c(4,1),
         col  = c("blue", "blue"), pt.cex = c(0.5,4), cex=1.1)
  
  
  cgi = SI.cgi(db.data,flag.plot=T,flag.inertia=T,col=2)
  plot(poly.data,col=NA,add=TRUE)
  
  
  #Calculate Results and put them in a table
  res = projec.invert(cgi$center[1],cgi$center[2])
  stats = SI.stats(db.data,flag.plot=F)
  
  res.table[i,3]=res$x
  res.table[i,4]=res$y
  res.table[i,5]=cgi$inertia
  res.table[i,6]=cgi$iso
  res.table[i,7]=stats$totab
  res.table[i,8]=stats$parea
  res.table[i,9]=stats$eqarea
  res.table[i,10]=stats$sparea
  res.table[i,11]=camargo_eveness(Data_Set$CATCH_NUM, FALSE)
  
}


write.csv(res.table, "C:/Users/oberc/OneDrive/Desktop/Flounders/Winter/res_table.csv")



res.table=read.csv("C:/Users/oberc/OneDrive/Desktop/Flounders/Winter/res_table.csv")
abundance_ind=read.csv("C:/Users/oberc/OneDrive/Desktop/LIS_HSI/Data/CTDEEP_AI_All_Lenth.csv")
ai_WINTERFL=abundance_ind[abundance_ind$Species=="WFL",]
ai_WINTERFL=ai_WINTERFL[ai_WINTERFL$Year<2022,]
res_all=left_join(res.table, ai_WINTERFL)
res.table=res_all

### Plots


##Spring

res_table_S=res.table[res.table$Season=="SP",]
res_table_S$Year=as.numeric(as.character(res_table_S$Year))

even=ggplot(res_table_S, aes(x=Year, y=Evenness)) + 
  geom_line() +
  labs(title = "Evenness", y = "Evenness", x ="Year")   + 
  theme(panel.background = element_blank()) +
  geom_smooth(method = 'loess',se = TRUE)


lon=ggplot(res_table_S, aes(x=Year, y=CG_Lon)) + 
  geom_line() +
  labs(title = "CG Longitude", y = "CG Lon", x ="Year")   + 
  theme(panel.background = element_blank()) +
  geom_smooth(method = 'loess',se = TRUE)

lat=ggplot(res_table_S, aes(x=Year, y=CG_Lat)) + 
  geom_line() +
  labs(title = "CG Latitiude", y = "CG Lat", x ="Year")   + 
  theme(panel.background = element_blank()) +
  geom_smooth(method = 'loess',se = TRUE)

inertia=ggplot(res_table_S, aes(x=Year, y=Inertia)) + 
  geom_line() +
  labs(title = "Inertia", y = "Inertia (nautical miles^2)", x ="Year")   + 
  theme(panel.background = element_blank()) +
  geom_smooth(method = 'loess',se = TRUE)

ai=ggplot(res_table_S, aes(x=Year, y=GmMn)) +
  geom_line() +
  labs(title = "Abundance Index", y = "Abundance Index", x ="Year")   +
  theme(panel.background = element_blank()) +
  geom_smooth(method = 'loess',se = TRUE)

pos_area=ggplot(res_table_S, aes(x=Year, y=Pos_Area)) + 
  geom_line() +
  labs(title = "Positive Area", y = "Positive Area (nautical miles^2)", x ="Year")   + 
  theme(panel.background = element_blank()) +
  geom_smooth(method = 'loess',se = TRUE)

plot_S=arrangeGrob(top=textGrob("Winter Flounder - Spring Spatial Indicators", gp=gpar(fontsize=20)),lon, lat, inertia, ai, pos_area, even, nrow = 2)

ggsave("C:/Users/oberc/OneDrive/Desktop/Flounders/Winter/WINTERFL_S_plot.png", plot_S, width = 30, height = 20, units = "cm")




## Fall

res_table_F=res.table[res.table$Season=="FA",]
res_table_F$Year=as.numeric(as.character(res_table_F$Year))

even=ggplot(res_table_F, aes(x=Year, y=Evenness)) + 
  geom_line() +
  labs(title = "Evenness", y = "Evenness", x ="Year")   + 
  theme(panel.background = element_blank()) +
  geom_smooth(method = 'loess',se = TRUE)


lon=ggplot(res_table_F, aes(x=Year, y=CG_Lon)) + 
  geom_line() +
  labs(title = "CG Longitude", y = "CG Lon", x ="Year")   + 
  theme(panel.background = element_blank()) +
  geom_smooth(method = 'loess',se = TRUE)

lat=ggplot(res_table_F, aes(x=Year, y=CG_Lat)) + 
  geom_line() +
  labs(title = "CG Latitiude", y = "CG Lat", x ="Year")   + 
  theme(panel.background = element_blank()) +
  geom_smooth(method = 'loess',se = TRUE)

inertia=ggplot(res_table_F, aes(x=Year, y=Inertia)) + 
  geom_line() +
  labs(title = "Inertia", y = "Inertia (nautical miles^2)", x ="Year")   + 
  theme(panel.background = element_blank()) +
  geom_smooth(method = 'loess',se = TRUE)

ai=ggplot(res_table_F, aes(x=Year, y=GmMn)) +
  geom_line() +
  labs(title = "Abundance Index", y = "Abundance Index", x ="Year")   +
  theme(panel.background = element_blank()) +
  geom_smooth(method = 'loess',se = TRUE)

pos_area=ggplot(res_table_F, aes(x=Year, y=Pos_Area)) + 
  geom_line() +
  labs(title = "Positive Area", y = "Positive Area (nautical miles^2)", x ="Year")   + 
  theme(panel.background = element_blank()) +
  geom_smooth(method = 'loess',se = TRUE)

plot_F=arrangeGrob(top=textGrob("Winter Flounder- Fall Spatial Indicators", gp=gpar(fontsize=20)),lon, lat, inertia, ai, pos_area, even, nrow = 2)

ggsave("C:/Users/oberc/OneDrive/Desktop/Flounders/Winter/WINTERFL_F_plot.png", plot_F, width = 30, height = 20, units = "cm")


## Center of Gravity Maps


res.table$Season=recode(res.table$Season, "SP" = "Spring", "FA"="Fall")
res.table$Year=as.numeric(as.character(res.table$Year))

##CG plot
#Changing ggplots 
theme_update(plot.title = element_text(hjust = 0.5))

CG_all=ggplot(res.table, aes(CG_Lon, CG_Lat))+
  geom_point(aes(color = Year), size=2) +
  scale_color_viridis(discrete= FALSE, direction= -1, option = "mako") +
  ggtitle("Winter Flounder Center of Gravity") +
  theme_bw()+theme(strip.text = element_text(size=10, face="bold"))+
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(strip.background = element_rect(fill="white"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  coord_quickmap(xlim = c(start_x, end_x), ylim=c(start_y, end_y)) +
  borders(database = "usa", regions = ".",fill = "gray90",xlim = NULL,ylim = NULL)+
  facet_wrap(~factor(Season, c("Spring", "Fall")), nrow=1, ncol=2)+
  labs(y = "Latitude", x ="Longitude")



##Save plot as a pdf
ggsave("C:/Users/oberc/OneDrive/Desktop/Flounders/Winter/WINTERFL_CG_plot.png", CG_all, width = 30, height = 10, units = "cm") 






####### Nice new spatial indicator plots for paper


res.table$Year <- as.numeric(as.character(res.table$Year))


res_table_S <- res.table[res.table$Season == "SP", ]
res_table_F <- res.table[res.table$Season == "FA", ]


res.table$Season <- recode(res.table$Season, "SP" = "Spring", "FA" = "Fall")

# Define indicator labels 
indicator_titles <- c(
  "Evenness" = "F) Evenness",
  "CG_Lon" = "A) Center of Gravity Longitude",
  "CG_Lat" = "B) Center of Gravity Latitude",
  "Inertia" = "C) Inertia",
  "GmMn" = "D) Abundance Index",
  "Pos_Area" = "E) Positive Area"
)


create_combined_plot <- function(data, y_var, label, y_label) {
  ggplot(data, aes(x = Year, y = !!sym(y_var), color = Season)) +
    geom_line(size = 1.2) +  # Thicker lines
    geom_smooth(method = 'loess', se = TRUE, size = 1.2) +  # Thicker smooth line
    scale_color_manual(values = c("Spring" = "blue", "Fall" = "red")) +
    labs(title = label, y = y_label, x = "Year") +
    theme_minimal() +
    theme(panel.background = element_blank(),
          plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
}

# Create combined plots for each indicator
even_plot <- create_combined_plot(res.table, "Evenness", indicator_titles["Evenness"], "Evenness")
lon_plot <- create_combined_plot(res.table, "CG_Lon", indicator_titles["CG_Lon"], "CG Longitude")
lat_plot <- create_combined_plot(res.table, "CG_Lat", indicator_titles["CG_Lat"], "CG Latitude")
inertia_plot <- create_combined_plot(res.table, "Inertia", indicator_titles["Inertia"], "Inertia (nautical miles²)")
ai_plot <- create_combined_plot(res.table, "GmMn", indicator_titles["GmMn"], "Abundance Index")
pos_area_plot <- create_combined_plot(res.table, "Pos_Area", indicator_titles["Pos_Area"], "Positive Area (nautical miles²)")

# Arrange 
combined_plot <- arrangeGrob(
  top = textGrob("Winter Flounder Spring and Fall Spatial Indicators", gp = gpar(fontsize = 20, fontface = "bold")),
  lon_plot, lat_plot, inertia_plot, ai_plot, pos_area_plot, even_plot,
  nrow = 2
)

plot(combined_plot)

# Save
ggsave("C:/Users/oberc/OneDrive/Desktop/Flounders/Winter/winter_Combined_Spring_Fall_Plot_new.png",
       combined_plot, width = 30, height = 20, units = "cm")





#########updated cog plot for paper. winter flounder


res.table$Season <- recode(res.table$Season, "SP" = "Spring", "FA" = "Fall")
res.table$Year <- as.numeric(as.character(res.table$Year))

# Identify first and last year for each season
res.table <- res.table %>%
  group_by(Season) %>%
  mutate(
    Symbol = case_when(
      Year == min(Year) ~ "First Year",
      Year == max(Year) ~ "Last Year",
      TRUE ~ "Intermediate Year"
    )
  )

#  shapes for First and Last Year only
symbol_shapes <- c("First Year" = 17, "Last Year" = 15)

# Updated Center of Gravity plot
CG_all <- ggplot(res.table, aes(CG_Lon, CG_Lat)) +
  geom_point(
    data = res.table %>% filter(Symbol == "Intermediate Year"),
    aes(color = Year),
    shape = 16, size = 2, show.legend = TRUE
  ) +
  geom_point(
    data = res.table %>% filter(Symbol %in% c("First Year", "Last Year")),
    aes(shape = Symbol),
    color = "red", size = 3, show.legend = TRUE
  ) +
  scale_color_viridis(discrete = FALSE, direction = -1, option = "mako", name = "Year") +
  scale_shape_manual(values = symbol_shapes, name = "") +  # Blank legend title for shapes
  guides(
    shape = guide_legend(override.aes = list(color = "red")),  # Shape legend for First/Last Year
    color = guide_colorbar(title = "Year")  # Color legend for intermediate years
  ) +
  ggtitle("Winter Flounder Center of Gravity") +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    strip.text = element_text(size = 10, face = "bold"),
    strip.background = element_rect(fill = "white"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  coord_quickmap(xlim = c(start_x, end_x), ylim = c(start_y, end_y)) +
  borders(database = "usa", regions = ".", fill = "gray90", xlim = NULL, ylim = NULL) +
  facet_wrap(~factor(Season, c("Spring", "Fall")), nrow = 1, ncol = 2) +
  labs(y = "Latitude", x = "Longitude")


plot(CG_all)


ggsave("C:/Users/oberc/OneDrive/Desktop/Flounders/Winter/WinterFL_CG_plot_updated.png", 
       CG_all, width = 30, height = 10, units = "cm")

