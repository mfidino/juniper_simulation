library(ggmap)
library(rgdal)

coyote <- creds[creds$species == "Coyote", ]
points(coyote[,5] ~ coyote[,1], type = "l", lwd = 2, col = "darkviolet")
points(coyote[,5] ~ coyote[,2], type = "l", lwd = 4, col = "darkviolet")
points(coyote[,5] ~ coyote[,3], type = "p", cex = 1, pch = 19, col = "darkviolet")

covs <- read.table("C:/Users/mfidino/Documents/R/bayes/bayes/cov_dat_1k_hous.txt", header = TRUE, sep = "\t")


cov_data <- read.csv("C:/Users/MFidino/Documents/R/Transect/SummaryFigures/RawMasterSheetforR.csv",
                     header = TRUE)
colnames(cov_data) <- tolower(colnames(cov_data))

to_keep <- c("transect", "section", "site", "stationid", "easting", "northing" )

cov_data  <-  cov_data[1:118,colnames(cov_data) %in% to_keep]
colnames(cov_data)[3] <- "site" 


utmcoor<-SpatialPoints(cbind(cov_data$easting,cov_data$northing),
                       proj4string=CRS("+proj=utm +zone=16N"))


longlatcoor<-spTransform(utmcoor,CRS("+proj=longlat"))

cov_data$easting <- longlatcoor@coords[,1]
cov_data$northing <- longlatcoor@coords[,2]

for_mapping <- data.frame(est = sp_rich$mode,
                          season = gl(10, 118),
                          lon = rep(cov_data$easting, 10),
                          lat = rep(cov_data$northing, 10))

chicago_map<-get_map(location=c(lon = -87.89529, lat =    41.89542),
                     maptype="toner-2011", zoom = 10)

windows(6,6)

for_mapping$est <- factor(for_mapping$est, ordered = TRUE)

setwd("./gif_stop/")
png(file = "map_plot%02d.png", height = 10, width = 10, units = "in", res = 400)


my_map <- ggmap(chicago_map)+
  geom_point(mapping = aes(x = lon, y = lat, fill = est, shape = est), size = 4, data = for_mapping[for_mapping$season == 1,])+
  scale_fill_discrete(h = c(0, 200),l = 65, c = 300, limits = c(0,1, 2, 3, 4, 5), name = "Predicted\nSpecies\nRichness")+
  xlab("Longitude")+
  ylab("Latitude")+
theme(text = element_text(size = 20),
      plot.margin = unit(c(1, 1, 1, 1), "cm"),
      legend.position = c(.93,.13))+
  scale_shape_manual(values = c(21, 21, 21, 21, 24, 24), name = "Predicted\nSpecies\nRichness" , limits = c(0,1, 2, 3, 4, 5))+
  annotate("text", x = -87.51, y = 42.20, label = "Spring", color = "White", size = 8)+
  annotate("text", x = -87.53, y = 42.17, label = "Summer", color = "gray50", size = 8)+
  annotate("text", x = - 87.49, y = 42.14, label = "Fall", color = "gray50", size = 8)+
  annotate("text", x = - 87.513, y = 42.11, label = "Winter", color = "gray50", size = 8)+
  annotate("text", x = - 87.50, y = 42.02, label = "2010", color = "White", size = 8)+
  annotate("text", x = - 87.50, y = 41.99, label = "2011", color = "gray50", size = 8)+
  annotate("text", x = - 87.50, y = 41.96, label = "2012", color = "gray50", size = 8)
# ggsave("my_map01.png", height = 10, width = 10, units = "in", dpi = 200)

my_sub <- ggplot(for_mapping[for_mapping$season == 1,], aes(x = est))+
  geom_histogram(aes(fill = est))+
  scale_fill_discrete(h = c(0, 200),l = 65, c = 300, limits = c(0,1, 2, 3, 4, 5))+
  xlab("Species Richness")+
  ylab("Count")+
  scale_x_discrete(limits = c("0", "1", "2", "3", "4", "5"))+
  scale_y_continuous(limits = c(0, 80))+
  theme(legend.position = "none")

windows(10, 10)

my_map
print(my_sub, vp = viewport(.384, .86, .17, .2))

my_map <- ggmap(chicago_map)+
  geom_point(mapping = aes(x = lon, y = lat, fill = est, shape = est), size = 4, data = for_mapping[for_mapping$season == 2,])+
  scale_fill_discrete(h = c(0, 200),l = 65, c = 300, limits = c(0,1, 2, 3, 4, 5), name = "Predicted\nSpecies\nRichness")+
  scale_shape_manual(values = c(21, 21, 21, 21, 24, 24), name = "Predicted\nSpecies\nRichness" , limits = c(0,1, 2, 3, 4, 5))+
  xlab("Longitude")+
  ylab("Latitude")+
  theme(text = element_text(size = 20),
        plot.margin = unit(c(1, 1, 1, 1), "cm"),
        legend.position = c(.93,.13))+
  annotate("text", x = -87.51, y = 42.20, label = "Spring", color = "gray50", size = 8)+
  annotate("text", x = -87.53, y = 42.17, label = "Summer", color = "white", size = 8)+
  annotate("text", x = - 87.49, y = 42.14, label = "Fall", color = "gray50", size = 8)+
  annotate("text", x = - 87.513, y = 42.11, label = "Winter", color = "gray50", size = 8)+
  annotate("text", x = - 87.50, y = 42.02, label = "2010", color = "White", size = 8)+
  annotate("text", x = - 87.50, y = 41.99, label = "2011", color = "gray50", size = 8)+
  annotate("text", x = - 87.50, y = 41.96, label = "2012", color = "gray50", size = 8)
# ggsave("my_map02.png", height = 10, width = 10, units = "in", dpi = 200)

my_sub <- ggplot(for_mapping[for_mapping$season == 2,], aes(x = est))+
  geom_histogram(aes(fill = est))+
  scale_fill_discrete(h = c(0, 200),l = 65, c = 300, limits = c(0,1, 2, 3, 4, 5))+
  xlab("Species Richness")+
  ylab("Count")+
  scale_x_discrete(limits = c("0", "1", "2", "3", "4", "5"))+
  scale_y_continuous(limits = c(0, 80))+
  theme(legend.position = "none")

windows(10, 10)

my_map
print(my_sub, vp = viewport(.384, .86, .17, .2))

my_map <- ggmap(chicago_map)+
  geom_point(mapping = aes(x = lon, y = lat, fill = est, shape = est), size = 4, data = for_mapping[for_mapping$season == 3,])+
  scale_fill_discrete(h = c(0, 200),l = 65, c = 300, limits = c(0,1, 2, 3, 4, 5), name = "Predicted\nSpecies\nRichness")+
  scale_shape_manual(values = c(21, 21, 21, 21, 24, 24), name = "Predicted\nSpecies\nRichness" , limits = c(0,1, 2, 3, 4, 5))+
  xlab("Longitude")+
  ylab("Latitude")+
  theme(text = element_text(size = 20),
        plot.margin = unit(c(1, 1, 1, 1), "cm"),
        legend.position = c(.93,.13))+
  annotate("text", x = -87.51, y = 42.20, label = "Spring", color = "gray50", size = 8)+
  annotate("text", x = -87.53, y = 42.17, label = "Summer", color = "gray50", size = 8)+
  annotate("text", x = - 87.49, y = 42.14, label = "Fall", color = "white", size = 8)+
  annotate("text", x = - 87.513, y = 42.11, label = "Winter", color = "gray50", size = 8)+
  annotate("text", x = - 87.50, y = 42.02, label = "2010", color = "White", size = 8)+
  annotate("text", x = - 87.50, y = 41.99, label = "2011", color = "gray50", size = 8)+
  annotate("text", x = - 87.50, y = 41.96, label = "2012", color = "gray50", size = 8)
# ggsave("my_map03.png", height = 10, width = 10, units = "in", dpi = 200)

my_sub <- ggplot(for_mapping[for_mapping$season == 3,], aes(x = est))+
  geom_histogram(aes(fill = est))+
  scale_fill_discrete(h = c(0, 200),l = 65, c = 300, limits = c(0,1, 2, 3, 4, 5))+
  xlab("Species Richness")+
  ylab("Count")+
  scale_x_discrete(limits = c("0", "1", "2", "3", "4", "5"))+
  scale_y_continuous(limits = c(0, 80))+
  theme(legend.position = "none")

windows(10, 10)

my_map
print(my_sub, vp = viewport(.384, .86, .17, .2))

my_map <- ggmap(chicago_map)+
  geom_point(mapping = aes(x = lon, y = lat, fill = est, shape = est), size = 4, data = for_mapping[for_mapping$season == 4,])+
  scale_fill_discrete(h = c(0, 200),l = 65, c = 300, limits = c(0,1, 2, 3, 4, 5), name = "Predicted\nSpecies\nRichness")+
  scale_shape_manual(values = c(21, 21, 21, 21, 24, 24), name = "Predicted\nSpecies\nRichness" , limits = c(0,1, 2, 3, 4, 5))+
  xlab("Longitude")+
  ylab("Latitude")+
  theme(text = element_text(size = 20),
        plot.margin = unit(c(1, 1, 1, 1), "cm"),
        legend.position = c(.93,.13))+
  annotate("text", x = -87.51, y = 42.20, label = "Spring", color = "gray50", size = 8)+
  annotate("text", x = -87.53, y = 42.17, label = "Summer", color = "gray50", size = 8)+
  annotate("text", x = - 87.49, y = 42.14, label = "Fall", color = "gray50", size = 8)+
  annotate("text", x = - 87.513, y = 42.11, label = "Winter", color = "white", size = 8)+
  annotate("text", x = - 87.50, y = 42.02, label = "2010", color = "gray50", size = 8)+
  annotate("text", x = - 87.50, y = 41.99, label = "2011", color = "white", size = 8)+
  annotate("text", x = - 87.50, y = 41.96, label = "2012", color = "gray50", size = 8)
# ggsave("my_map04.png", height = 10, width = 10, units = "in", dpi = 200)

my_sub <- ggplot(for_mapping[for_mapping$season == 4,], aes(x = est))+
  geom_histogram(aes(fill = est))+
  scale_fill_discrete(h = c(0, 200),l = 65, c = 300, limits = c(0,1, 2, 3, 4, 5))+
  xlab("Species Richness")+
  ylab("Count")+
  scale_x_discrete(limits = c("0", "1", "2", "3", "4", "5"))+
  scale_y_continuous(limits = c(0, 80))+
  theme(legend.position = "none")

windows(10, 10)

my_map
print(my_sub, vp = viewport(.384, .86, .17, .2))

my_map <- ggmap(chicago_map)+
  geom_point(mapping = aes(x = lon, y = lat, fill = est, shape = est), size = 4, data = for_mapping[for_mapping$season == 5,])+
  scale_fill_discrete(h = c(0, 200),l = 65, c = 300, limits = c(0,1, 2, 3, 4, 5), name = "Predicted\nSpecies\nRichness")+
  scale_shape_manual(values = c(21, 21, 21, 21, 24, 24), name = "Predicted\nSpecies\nRichness" , limits = c(0,1, 2, 3, 4, 5))+
  xlab("Longitude")+
  ylab("Latitude")+
  theme(text = element_text(size = 20),
        plot.margin = unit(c(1, 1, 1, 1), "cm"),
        legend.position = c(.93,.13))+
  annotate("text", x = -87.51, y = 42.20, label = "Spring", color = "white", size = 8)+
  annotate("text", x = -87.53, y = 42.17, label = "Summer", color = "gray50", size = 8)+
  annotate("text", x = - 87.49, y = 42.14, label = "Fall", color = "gray50", size = 8)+
  annotate("text", x = - 87.513, y = 42.11, label = "Winter", color = "gray50", size = 8)+
  annotate("text", x = - 87.50, y = 42.02, label = "2010", color = "gray50", size = 8)+
  annotate("text", x = - 87.50, y = 41.99, label = "2011", color = "white", size = 8)+
  annotate("text", x = - 87.50, y = 41.96, label = "2012", color = "gray50", size = 8)
# ggsave("my_map05.png", height = 10, width = 10, units = "in", dpi = 200)

my_sub <- ggplot(for_mapping[for_mapping$season == 9,], aes(x = est))+
  geom_histogram(aes(fill = est))+
  scale_fill_discrete(h = c(0, 200),l = 65, c = 300, limits = c(0,1, 2, 3, 4, 5))+
  xlab("Species Richness")+
  ylab("Count")+
  scale_x_discrete(limits = c("0", "1", "2", "3", "4", "5"))+
  scale_y_continuous(limits = c(0, 80))+
  theme(legend.position = "none")

windows(10, 10)

my_map
print(my_sub, vp = viewport(.384, .86, .17, .2))

my_map <- ggmap(chicago_map)+
  geom_point(mapping = aes(x = lon, y = lat, fill = est, shape = est), size = 4, data = for_mapping[for_mapping$season == 6,])+
  scale_fill_discrete(h = c(0, 200),l = 65, c = 300, limits = c(0,1, 2, 3, 4, 5), name = "Predicted\nSpecies\nRichness")+
  scale_shape_manual(values = c(21, 21, 21, 21, 24, 24), name = "Predicted\nSpecies\nRichness" , limits = c(0,1, 2, 3, 4, 5))+
  xlab("Longitude")+
  ylab("Latitude")+
  theme(text = element_text(size = 20),
        plot.margin = unit(c(1, 1, 1, 1), "cm"),
        legend.position = c(.93,.13))+
  annotate("text", x = -87.51, y = 42.20, label = "Spring", color = "gray50", size = 8)+
  annotate("text", x = -87.53, y = 42.17, label = "Summer", color = "white", size = 8)+
  annotate("text", x = - 87.49, y = 42.14, label = "Fall", color = "gray50", size = 8)+
  annotate("text", x = - 87.513, y = 42.11, label = "Winter", color = "gray50", size = 8)+
  annotate("text", x = - 87.50, y = 42.02, label = "2010", color = "gray50", size = 8)+
  annotate("text", x = - 87.50, y = 41.99, label = "2011", color = "white", size = 8)+
  annotate("text", x = - 87.50, y = 41.96, label = "2012", color = "gray50", size = 8)
# ggsave("my_map06.png", height = 10, width = 10, units = "in", dpi = 200)

my_sub <- ggplot(for_mapping[for_mapping$season == 6,], aes(x = est))+
  geom_histogram(aes(fill = est))+
  scale_fill_discrete(h = c(0, 200),l = 65, c = 300, limits = c(0,1, 2, 3, 4, 5))+
  xlab("Species Richness")+
  ylab("Count")+
  scale_x_discrete(limits = c("0", "1", "2", "3", "4", "5"))+
  scale_y_continuous(limits = c(0, 80))+
  theme(legend.position = "none")

windows(10, 10)

my_map
print(my_sub, vp = viewport(.384, .86, .17, .2))

my_map <- ggmap(chicago_map)+
  geom_point(mapping = aes(x = lon, y = lat, fill = est, shape = est), size = 4, data = for_mapping[for_mapping$season == 7,])+
  scale_fill_discrete(h = c(0, 200),l = 65, c = 300, limits = c(0,1, 2, 3, 4, 5), name = "Predicted\nSpecies\nRichness")+
  scale_shape_manual(values = c(21, 21, 21, 21, 24, 24), name = "Predicted\nSpecies\nRichness" , limits = c(0,1, 2, 3, 4, 5))+
  xlab("Longitude")+
  ylab("Latitude")+
  theme(text = element_text(size = 20),
        plot.margin = unit(c(1, 1, 1, 1), "cm"),
        legend.position = c(.93,.13))+
  annotate("text", x = -87.51, y = 42.20, label = "Spring", color = "gray50", size = 8)+
  annotate("text", x = -87.53, y = 42.17, label = "Summer", color = "gray50", size = 8)+
  annotate("text", x = - 87.49, y = 42.14, label = "Fall", color = "white", size = 8)+
  annotate("text", x = - 87.513, y = 42.11, label = "Winter", color = "gray50", size = 8)+
  annotate("text", x = - 87.50, y = 42.02, label = "2010", color = "gray50", size = 8)+
  annotate("text", x = - 87.50, y = 41.99, label = "2011", color = "white", size = 8)+
  annotate("text", x = - 87.50, y = 41.96, label = "2012", color = "gray50", size = 8)
# ggsave("my_map07.png", height = 10, width = 10, units = "in", dpi = 200)

my_sub <- ggplot(for_mapping[for_mapping$season == 7,], aes(x = est))+
  geom_histogram(aes(fill = est))+
  scale_fill_discrete(h = c(0, 200),l = 65, c = 300, limits = c(0,1, 2, 3, 4, 5))+
  xlab("Species Richness")+
  ylab("Count")+
  scale_x_discrete(limits = c("0", "1", "2", "3", "4", "5"))+
  scale_y_continuous(limits = c(0, 80))+
  theme(legend.position = "none")

windows(10, 10)

my_map
print(my_sub, vp = viewport(.384, .86, .17, .2))

my_map <- ggmap(chicago_map)+
  geom_point(mapping = aes(x = lon, y = lat, fill = est, shape = est), size = 4, data = for_mapping[for_mapping$season == 8,])+
  scale_fill_discrete(h = c(0, 200),l = 65, c = 300, limits = c(0,1, 2, 3, 4, 5), name = "Predicted\nSpecies\nRichness")+
  scale_shape_manual(values = c(21, 21, 21, 21, 24, 24), name = "Predicted\nSpecies\nRichness" , limits = c(0,1, 2, 3, 4, 5))+
  xlab("Longitude")+
  ylab("Latitude")+
  theme(text = element_text(size = 20),
        plot.margin = unit(c(1, 1, 1, 1), "cm"),
        legend.position = c(.93,.13))+
  annotate("text", x = -87.51, y = 42.20, label = "Spring", color = "gray50", size = 8)+
  annotate("text", x = -87.53, y = 42.17, label = "Summer", color = "gray50", size = 8)+
  annotate("text", x = - 87.49, y = 42.14, label = "Fall", color = "gray50", size = 8)+
  annotate("text", x = - 87.513, y = 42.11, label = "Winter", color = "white", size = 8)+
  annotate("text", x = - 87.50, y = 42.02, label = "2010", color = "gray50", size = 8)+
  annotate("text", x = - 87.50, y = 41.99, label = "2011", color = "gray50", size = 8)+
  annotate("text", x = - 87.50, y = 41.96, label = "2012", color = "white", size = 8)
# ggsave("my_map08.png", height = 10, width = 10, units = "in", dpi = 200)

my_sub <- ggplot(for_mapping[for_mapping$season == 8,], aes(x = est))+
  geom_histogram(aes(fill = est))+
  scale_fill_discrete(h = c(0, 200),l = 65, c = 300, limits = c(0,1, 2, 3, 4, 5))+
  xlab("Species Richness")+
  ylab("Count")+
  scale_x_discrete(limits = c("0", "1", "2", "3", "4", "5"))+
  scale_y_continuous(limits = c(0, 80))+
  theme(legend.position = "none")

windows(10, 10)

my_map
print(my_sub, vp = viewport(.384, .86, .17, .2))

my_map <- ggmap(chicago_map)+
  geom_point(mapping = aes(x = lon, y = lat, fill = est, shape = est), size = 4, data = for_mapping[for_mapping$season == 9,])+
  scale_fill_discrete(h = c(0, 200),l = 65, c = 300, limits = c(0,1, 2, 3, 4, 5), name = "Predicted\nSpecies\nRichness")+
  scale_shape_manual(values = c(21, 21, 21, 21, 24, 24), name = "Predicted\nSpecies\nRichness" , limits = c(0,1, 2, 3, 4, 5))+
  xlab("Longitude")+
  ylab("Latitude")+
  theme(text = element_text(size = 20),
        plot.margin = unit(c(1, 1, 1, 1), "cm"),
        legend.position = c(.93,.13))+
  annotate("text", x = -87.51, y = 42.20, label = "Spring", color = "white", size = 8)+
  annotate("text", x = -87.53, y = 42.17, label = "Summer", color = "gray50", size = 8)+
  annotate("text", x = - 87.49, y = 42.14, label = "Fall", color = "gray50", size = 8)+
  annotate("text", x = - 87.513, y = 42.11, label = "Winter", color = "gray50", size = 8)+
  annotate("text", x = - 87.50, y = 42.02, label = "2010", color = "gray50", size = 8)+
  annotate("text", x = - 87.50, y = 41.99, label = "2011", color = "gray50", size = 8)+
  annotate("text", x = - 87.50, y = 41.96, label = "2012", color = "white", size = 8)
# ggsave("my_map09.png", height = 10, width = 10, units = "in", dpi = 200)

my_sub <- ggplot(for_mapping[for_mapping$season == 9,], aes(x = est))+
  geom_histogram(aes(fill = est))+
  scale_fill_discrete(h = c(0, 200),l = 65, c = 300, limits = c(0,1, 2, 3, 4, 5))+
  xlab("Species Richness")+
  ylab("Count")+
  scale_x_discrete(limits = c("0", "1", "2", "3", "4", "5"))+
  scale_y_continuous(limits = c(0, 80))+
  theme(legend.position = "none")

windows(10, 10)

my_map
print(my_sub, vp = viewport(.384, .86, .17, .2))


my_map <- ggmap(chicago_map)+
  geom_point(mapping = aes(x = lon, y = lat, fill = est, shape = est), size = 4, data = for_mapping[for_mapping$season == 10,])+
  scale_fill_discrete(h = c(0, 200),l = 65, c = 300, limits = c(0,1, 2, 3, 4, 5), name = "Predicted\nSpecies\nRichness")+
  scale_shape_manual(values = c(21, 21, 21, 21, 24, 24), name = "Predicted\nSpecies\nRichness" , limits = c(0,1, 2, 3, 4, 5))+
  xlab("Longitude")+
  ylab("Latitude")+
  theme(text = element_text(size = 20),
        plot.margin = unit(c(1, 1, 1, 1), "cm"),
        legend.position = c(.93,.13))+
  annotate("text", x = -87.51, y = 42.20, label = "Spring", color = "gray50", size = 8)+
  annotate("text", x = -87.53, y = 42.17, label = "Summer", color = "white", size = 8)+
  annotate("text", x = - 87.49, y = 42.14, label = "Fall", color = "gray50", size = 8)+
  annotate("text", x = - 87.513, y = 42.11, label = "Winter", color = "gray50", size = 8)+
  annotate("text", x = - 87.50, y = 42.02, label = "2010", color = "gray50", size = 8)+
  annotate("text", x = - 87.50, y = 41.99, label = "2011", color = "gray50", size = 8)+
  annotate("text", x = - 87.50, y = 41.96, label = "2012", color = "white", size = 8)
# ggsave("my_map10.png", height = 10, width = 10, units = "in", dpi = 200)
  
my_sub <- ggplot(for_mapping[for_mapping$season == 10,], aes(x = est))+
  geom_histogram(aes(fill = est))+
  scale_fill_discrete(h = c(0, 200),l = 65, c = 300, limits = c(0,1, 2, 3, 4, 5))+
  xlab("Species Richness")+
  ylab("Count")+
  scale_x_discrete(limits = c("0", "1", "2", "3", "4", "5"))+
  scale_y_continuous(limits = c(0, 80))+
  theme(legend.position = "none")

windows(10, 10)

my_map
print(my_sub, vp = viewport(.384, .86, .17, .2))

test <- arrangeGrob(my_map, main = textGrob("Spring", hjust = 4.25, vjust = 5, gp = gpar(fontsize = 25, face = "bold", col = "green")))
test <- arrangeGrob(test, main = textGrob("Summer", hjust = 2.25, vjust = 8, gp = gpar(fontsize = 20, col = "black")))

dev.off()
system('convert -delay 100 *.png butts.gif')

