

covar_est <- m_sum[grep("beta_est", row.names(m_sum)),-1]
covar_est$species <- factor(rep(c("coyote", "opossum", "raccoon", "red_fox", "skunk"), 100))


# scaled covariate for prediction
housing_data <- seq(-1, 1, length.out = 100)
# unscaled
housing_unscale <- seq(0, 7655.838, length.out = length(housing_data))

covar_est$housing <- rep(housing_unscale, each = 5)



library(reshape2)
library(ggplot2)
library(gridExtra)
library(RColorBrewer)


windows(height = 12, width = 12)



plot_theme <-   theme(axis.line = element_line(colour = "black"),
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      panel.border = element_blank(),
                      panel.background = element_blank(),
                      text = element_text(size = 20, color = "black"),
                      axis.ticks = element_line(color = "black"))


p_color = "darkviolet"
coyote <- ggplot(covar_est[covar_est$species =="coyote",], 
                 aes(x = housing, y = mode))+
  geom_line(aes(x = housing, y = uci),linetype = "dotdash", color = p_color)+
  geom_line(aes(x = housing, y = lci), linetype = "dotdash", color = p_color)+
  geom_line(size = 1.25, color = p_color)+
  ylab("Probability of colonization")+
  scale_y_continuous(limits = c(0,.8))+
  xlab("Housing Density")+
  ggtitle("Coyote")+
  plot_theme

windows(height = 24, width = 40)
grid.arrange(coyote, raccoon, opossum, fox, skunk, my_cat, ncol = 3)

windows(height = 12, width = 12)
caterplot(mod_mcmc, regex = "beta[^_]", reorder = FALSE,
          labels = c("Coyote", "Opossum", "Raccoon", "Red Fox", "Striped Skunk"),
          style = "plain", cex = 1.25, labels.loc = "above")
abline(v = 0, col = "red")
title(xlab = "logit scale beta estimates with 95% credible intervals", cex.lab = 1.25)


# make a big model matrix
mod_mat <- as.matrix(mod_mcmc, chain = TRUE)
# just grab the betas
mod_mat <- mod_mat[, 23:27]

# get credible intervals
# get quantiles
creds <- apply(mod_mat, 2, quantile, probs = c(0.025, 0.975, 0.25,  0.75, 0.50))
# change column anmes
colnames(creds) <- c("Coyote", "Opossum", "Raccoon", "Red Fox", "Striped Skunk")
# make long format
creds <- melt(creds)
# make wide format
creds <- dcast(creds, Var2~ Var1)
# change column names
colnames(creds) <- c("species", "ci", "uci", "ci2", "uci2", "mode")
# add some more mode for the matrix
creds$mode2 <- creds$mode
# make it a matrix without species
creds <- data.frame(matrix(creds, ncol = 3 ))
# change coumn names
colnames(creds) <- c("ci", "ci2", "mode")
# add a species column
creds$species <- factor(rep(c("Coyote", "Opossum", "Raccoon", "Red Fox", "Striped Skunk"), 2),
                        levels = c("Striped Skunk", "Red Fox", "Opossum", 
                                   "Raccoon", "Coyote"))



# plot it out
my_cat <- ggplot(creds, aes(x = ci, y = species, color = species))+
  geom_line(size = 1.25)+
  geom_line(aes(x = ci2, y = species), size = 1.75)+
  geom_point(aes(x = mode, y = species, size = 4))+
  plot_theme+
  xlab("logit scale beta estimates")+
  ylab("")+
  geom_vline(x = 0, color = "red")+
  theme(legend.position = "none")+
  scale_color_manual(values = c("black", "olivedrab4", 
                                "dodgerblue4", "firebrick2", "darkviolet"))


# plot out known vs unknown

sp_rich <- m_sum[grep("sp_rich", m_sum$parameter),]
sp_rich$season <- factor(rep(c(1:10), each = 118))
known_sp <- apply(zinit, c(2,3),sum,  na.rm = TRUE)
days_samp <- apply(jmat, c(2,3), max)
sp_rich$known <- as.vector(known_sp)
sp_rich$site <- factor(rep(c(1:118), 10))
sp_rich$day_samp <- factor(floor(1 + (3 * as.vector(days_samp) / max(days_samp))))
levels(sp_rich$day_samp) <- c("Sample: 0-9 days", "Sample: 10-18 days", 
                              "Sample: 19-27 days", "Sample: 28 days")
sp_rich$day_samp <- factor(sp_rich$day_samp, levels = c("Sample: 28 days",
                                                        "Sample: 19-27 days",
                                                        "Sample: 10-18 days",
                                                        "Sample: 0-9 days"))
sp_rich$day_samp_num <- as.vector(days_samp)

limits <- aes(ymax = uci, ymin = lci)

create_plot <- function(species = NULL, my_col = NULL){
  
  simpleCap <- function(x) {
    s <- strsplit(x, " ")[[1]]
    paste(toupper(substring(s, 1, 1)), substring(s, 2),
          sep = "", collapse = " ")
  }

ym<-covar_est[covar_est$species==species,3]
yl<-covar_est[covar_est$species==species,2]
yu<-covar_est[covar_est$species==species,4]
xm<-covar_est[covar_est$species==species,6]


tiff(paste(species, ".tif", sep = ""), height = 4.5, width = 5.35, units="in", res=400)

par(mar=c(4.5,4.5,1.5,1.5))
plot(1~1, type='n', xlim=c(0,8000), ylim=c(0,1), xlab="", ylab="", xaxt="n", yaxt="n")

axis(1, at=seq(0,8000, 2000), labels=F, tck=-.035)
axis(1, at=seq(0,8000, 1000), labels=F, tck=-.02)
axis(1, at=seq(0,8000, 500), labels=F, tck=-.01)
mtext(text=seq(0,8000, 2000), 1, line=1, at=seq(0,8000, 2000), cex=1.25)

axis(2, at=seq(0,1, .2), labels=F, tck=-.035)
axis(2, at=seq(0,1, .1), labels=F, tck=-.02)
axis(2, at=seq(0,1, .05), labels=F, tck=-.01)
mtext(text=seq(0,1, .2), 2, line=.75, at=seq(0,1, .2), cex=1.25, las=1)

mtext("Housing Density", 1, line=3, cex=1.5)
mtext(expression(paste("Probability of Colonization (", gamma, ")", sep = "")), 2, line=3, cex=1.5)

points(ym~xm, type='l', lwd=4, col=my_col)
points(yl~xm, type='l', lwd=2, lty=3, col=my_col)
points(yu~xm, type='l', lwd=2, lty=3, col=my_col)

if(species == "red_fox"){
  mtext("Red Fox", 3, line = 0.2, cex = 1.5)
}
if(species == "skunk"){
  mtext("Skunk", 3, line = 0.2, cex = 1.5)
}

if(species %in% c("coyote", "opossum", "raccoon")){
mtext(simpleCap(species), 3, line = 0.2, cex = 1.5)
}
dev.off()
}

create_plot("raccoon", "firebrick2")
create_plot("coyote", "darkviolet")
create_plot("opossum", "dodgerblue4")
create_plot("red_fox", "olivedrab4")
create_plot("skunk", "black")



windows(4.5, 4.5, )


mod_mat <- mod_matrix[, 23:27]

# get credible intervals
# get quantiles
creds <- apply(mod_mat, 2, quantile, probs = c(0.025, 0.975, 0.25,  0.75, 0.50))
# change column anmes
colnames(creds) <- c("Coyote", "Opossum", "Raccoon", "Red Fox", "Striped Skunk")
# make long format
creds <- melt(creds)
# make wide format
creds <- dcast(creds, Var2~ Var1)
# change column names
colnames(creds) <- c("species", "ci", "uci", "ci2", "uci2", "mode")
# add some more mode for the matrix
creds$mode2 <- creds$mode
# remove species
creds <- as.matrix(creds[,-1])
# make it a matrix without species
creds <- data.frame(matrix(creds, ncol = 3 ))
# change coumn names
colnames(creds) <- c("ci", "ci2", "mode")
# add a species column
creds$species <- factor(rep(c("Coyote", "Opossum", "Raccoon", "Red Fox", "Striped Skunk"), 2),
                        levels = c("Striped Skunk", "Red Fox", "Opossum", 
                                   "Raccoon", "Coyote"))

creds$sp <- as.numeric(creds$species)

tiff("housing_density_slope.tif", height = 3.5, width = 6.3, units = "in", res = 400)
par(mar=c(3.5,7.0,2.5,1))
plot(1~1, type='n', xlim=c(-2,.5),ylim = c(.8,5.2), xlab="", ylab="", xaxt="n", yaxt="n")
axis(1, at = seq(-2, .5, .5), labels=F, tck=-.035)
axis(1, at=seq(-2, .5, .25), labels=F, tck=-.02)
axis(1, at=seq(-2, .5, .125), labels=F, tck=-.01)
mtext(text=seq(-2, .5, .5), 1, line=1, at=seq(-2, .5, .5), cex=1.25)

axis(2, at = seq(1, 5, 1), labels=F, tck=-.035)
mtext(text = c("Skunk", "Red Fox", "Opossum", "Raccoon", "Coyote"),
      2, line = .75, at = seq(1, 5, 1), las = 1, cex = 1.5)

skunk <- creds[creds$species == "Striped Skunk", ]
points(skunk[,4] ~ skunk[,1], type = "l", lwd = 3.5)
points(skunk[,4] ~ skunk[,2], type = "l", lwd = 6.5)
points(skunk[,4] ~ skunk[,3], type = "p", cex = 1.75, pch = 19, col=0)
points(skunk[,4] ~ skunk[,3], type = "p", cex = 1.75, pch = 10, lwd=2)

fox <- creds[creds$species == "Red Fox", ]
points(fox[,4] ~ fox[,1], type = "l", lwd = 3.5, col = "olivedrab4")
points(fox[,4] ~ fox[,2], type = "l", lwd = 6.5, col = "olivedrab4")
points(fox[,4] ~ fox[,3], type = "p", cex = 1.75, pch = 19, col =0)
points(fox[,4] ~ fox[,3], type = "p", cex = 1.75, pch = 10, col ="olivedrab4", lwd = 2)

possum <- creds[creds$species == "Opossum", ]
points(possum[,4] ~ possum[,1], type = "l", lwd = 3.5, col = "dodgerblue4")
points(possum[,4] ~ possum[,2], type = "l", lwd = 6.5, col = "dodgerblue4")
points(possum[,4] ~ possum[,3], type = "p", cex = 1.75, pch = 19, col = 0)
points(possum[,4] ~ possum[,3], type = "p", cex = 1.75, pch = 10, col = "dodgerblue4", lwd = 2)

raccoon <- creds[creds$species == "Raccoon", ]
points(raccoon[,4] ~ raccoon[,1], type = "l", lwd = 3.5, col = "firebrick2")
points(raccoon[,4] ~ raccoon[,2], type = "l", lwd = 6.5, col = "firebrick2")
points(raccoon[,4] ~ raccoon[,3], type = "p", cex = 1.75, pch = 19, col = 0)
points(raccoon[,4] ~ raccoon[,3], type = "p", cex = 1.75, pch = 10, col = "firebrick2", lwd = 2)

coyote <- creds[creds$species == "Coyote", ]
points(coyote[,4] ~ coyote[,1], type = "l", lwd = 3.5, col = "darkviolet")
points(coyote[,4] ~ coyote[,2], type = "l", lwd = 6.5, col = "darkviolet")
points(coyote[,4] ~ coyote[,3], type = "p", cex = 1.75, pch = 19, col = 0)
points(coyote[,4] ~ coyote[,3], type = "p", cex = 1.75, pch = 10, col = "darkviolet", lwd = 2)

abline(v = 0, lty = 3)
mtext(expression(paste("Impact of Housing Density on Colonization (", gamma, ")", sep = "")), 3, line = 0.4, cex = 1.4)

dev.off()


mod_matrix <- as.matrix(mod_mcmc, chain = TRUE)
# just grab the betas
mod_mat <- mod_matrix[, 2:6]

# get credible intervals
# get quantiles
creds <- apply(mod_mat, 2, quantile, probs = c(0.025, 0.975, 0.25,  0.75, 0.50))
# change column anmes
colnames(creds) <- c("Coyote", "Opossum", "Raccoon", "Red Fox", "Striped Skunk")
# make long format
creds <- melt(creds)
# make wide format
creds <- dcast(creds, Var2~ Var1)
# change column names
colnames(creds) <- c("species", "ci", "uci", "ci2", "uci2", "mode")
# add some more mode for the matrix
creds$mode2 <- creds$mode
# remove species
creds <- as.matrix(creds[,-1])
# make it a matrix without species
creds <- data.frame(matrix(creds, ncol = 3 ))
# change coumn names
colnames(creds) <- c("ci", "ci2", "mode")
# add a species column
creds$species <- factor(rep(c("Coyote", "Opossum", "Raccoon", "Red Fox", "Striped Skunk"), 2),
                        levels = c("Striped Skunk", "Red Fox", "Opossum", 
                                   "Raccoon", "Coyote"))

creds$sp <- as.numeric(creds$species)

tiff("colonization_cater.tif", height = 3.5, width = 6.3, units = "in", res = 400)
par(mar=c(3.5,7.0,2.5,1))
plot(1~1, type='n', xlim=c(-3.5,0),ylim = c(.8,5.2), xlab="", ylab="", xaxt="n", yaxt="n")
axis(1, at = seq(-3.5, 0, .5), labels=F, tck=-.035)
axis(1, at=seq(-3.5, 0, .25), labels=F, tck=-.02)
axis(1, at=seq(-3.5, 0, .125), labels=F, tck=-.01)
mtext(text=seq(-3.5, 0, .5), 1, line=1, at=seq(-3.5, 0, .5), cex=1.25)

axis(2, at = seq(1, 5, 1), labels=F, tck=-.035)
mtext(text = c("Skunk", "Red Fox", "Opossum", "Raccoon", "Coyote"),
      2, line = .75, at = seq(1, 5, 1), las = 1, cex = 1.5)

skunk <- creds[creds$species == "Striped Skunk", ]
points(skunk[,5] ~ skunk[,1], type = "l", lwd = 3.5)
points(skunk[,5] ~ skunk[,2], type = "l", lwd = 6.5)
points(skunk[,5] ~ skunk[,3], type = "p", cex = 1.75, pch = 19, col=0)
points(skunk[,5] ~ skunk[,3], type = "p", cex = 1.75, pch = 10, lwd=2)

fox <- creds[creds$species == "Red Fox", ]
points(fox[,5] ~ fox[,1], type = "l", lwd = 3.5, col = "olivedrab4")
points(fox[,5] ~ fox[,2], type = "l", lwd = 6.5, col = "olivedrab4")
points(fox[,5] ~ fox[,3], type = "p", cex = 1.75, pch = 19, col =0)
points(fox[,5] ~ fox[,3], type = "p", cex = 1.75, pch = 10, col ="olivedrab4", lwd = 2)

possum <- creds[creds$species == "Opossum", ]
points(possum[,5] ~ possum[,1], type = "l", lwd = 3.5, col = "dodgerblue4")
points(possum[,5] ~ possum[,2], type = "l", lwd = 6.5, col = "dodgerblue4")
points(possum[,5] ~ possum[,3], type = "p", cex = 1.75, pch = 19, col = 0)
points(possum[,5] ~ possum[,3], type = "p", cex = 1.75, pch = 10, col = "dodgerblue4", lwd = 2)

raccoon <- creds[creds$species == "Raccoon", ]
points(raccoon[,5] ~ raccoon[,1], type = "l", lwd = 3.5, col = "firebrick2")
points(raccoon[,5] ~ raccoon[,2], type = "l", lwd = 6.5, col = "firebrick2")
points(raccoon[,5] ~ raccoon[,3], type = "p", cex = 1.75, pch = 19, col = 0)
points(raccoon[,5] ~ raccoon[,3], type = "p", cex = 1.75, pch = 10, col = "firebrick2", lwd = 2)

coyote <- creds[creds$species == "Coyote", ]
points(coyote[,5] ~ coyote[,1], type = "l", lwd = 3.5, col = "darkviolet")
points(coyote[,5] ~ coyote[,2], type = "l", lwd = 6.5, col = "darkviolet")
points(coyote[,5] ~ coyote[,3], type = "p", cex = 1.75, pch = 19, col = 0)
points(coyote[,5] ~ coyote[,3], type = "p", cex = 1.75, pch = 10, col = "darkviolet", lwd = 2)

abline(v = 0, lty = 3)
mtext(expression(paste("Colonization (", gamma, ") Intercept", sep = "")), 3, line = 0.4, cex = 1.4)

dev.off()




mod_mat <- mod_matrix[, 7:11]

# get credible intervals
# get quantiles
creds <- apply(mod_mat, 2, quantile, probs = c(0.025, 0.975, 0.25,  0.75, 0.50))
# change column anmes
colnames(creds) <- c("Coyote", "Opossum", "Raccoon", "Red Fox", "Striped Skunk")
# make long format
creds <- melt(creds)
# make wide format
creds <- dcast(creds, Var2~ Var1)
# change column names
colnames(creds) <- c("species", "ci", "uci", "ci2", "uci2", "mode")
# add some more mode for the matrix
creds$mode2 <- creds$mode
# remove species
creds <- as.matrix(creds[,-1])
# make it a matrix without species
creds <- data.frame(matrix(creds, ncol = 3 ))
# change coumn names
colnames(creds) <- c("ci", "ci2", "mode")
# add a species column
creds$species <- factor(rep(c("Coyote", "Opossum", "Raccoon", "Red Fox", "Striped Skunk"), 2),
                        levels = c("Striped Skunk", "Red Fox", "Opossum", 
                                   "Raccoon", "Coyote"))

creds$sp <- as.numeric(creds$species)

tiff("persistence_cater.tif", height = 3.5, width = 6.3, units = "in", res = 400)
par(mar=c(3.5,7.0,2.5,1))
plot(1~1, type='n', xlim=c(-2,1.5),ylim = c(.8,5.2), xlab="", ylab="", xaxt="n", yaxt="n")
axis(1, at = seq(-2, 1.5, .5), labels=F, tck=-.035)
axis(1, at=seq(-2, 1.5, .25), labels=F, tck=-.02)
axis(1, at=seq(-2, 1.5, .125), labels=F, tck=-.01)
mtext(text=seq(-2, 1.5, .5), 1, line=1, at=seq(-2, 1.5, .5), cex=1.25)

axis(2, at = seq(1, 5, 1), labels=F, tck=-.035)
mtext(text = c("Skunk", "Red Fox", "Opossum", "Raccoon", "Coyote"),
      2, line = .75, at = seq(1, 5, 1), las = 1, cex = 1.5)

skunk <- creds[creds$species == "Striped Skunk", ]
points(skunk[,5] ~ skunk[,1], type = "l", lwd = 3.5)
points(skunk[,5] ~ skunk[,2], type = "l", lwd = 6.5)
points(skunk[,5] ~ skunk[,3], type = "p", cex = 1.75, pch = 19, col=0)
points(skunk[,5] ~ skunk[,3], type = "p", cex = 1.75, pch = 10, lwd=2)

fox <- creds[creds$species == "Red Fox", ]
points(fox[,5] ~ fox[,1], type = "l", lwd = 3.5, col = "olivedrab4")
points(fox[,5] ~ fox[,2], type = "l", lwd = 6.5, col = "olivedrab4")
points(fox[,5] ~ fox[,3], type = "p", cex = 1.75, pch = 19, col =0)
points(fox[,5] ~ fox[,3], type = "p", cex = 1.75, pch = 10, col ="olivedrab4", lwd = 2)

possum <- creds[creds$species == "Opossum", ]
points(possum[,5] ~ possum[,1], type = "l", lwd = 3.5, col = "dodgerblue4")
points(possum[,5] ~ possum[,2], type = "l", lwd = 6.5, col = "dodgerblue4")
points(possum[,5] ~ possum[,3], type = "p", cex = 1.75, pch = 19, col = 0)
points(possum[,5] ~ possum[,3], type = "p", cex = 1.75, pch = 10, col = "dodgerblue4", lwd = 2)

raccoon <- creds[creds$species == "Raccoon", ]
points(raccoon[,5] ~ raccoon[,1], type = "l", lwd = 3.5, col = "firebrick2")
points(raccoon[,5] ~ raccoon[,2], type = "l", lwd = 6.5, col = "firebrick2")
points(raccoon[,5] ~ raccoon[,3], type = "p", cex = 1.75, pch = 19, col = 0)
points(raccoon[,5] ~ raccoon[,3], type = "p", cex = 1.75, pch = 10, col = "firebrick2", lwd = 2)

coyote <- creds[creds$species == "Coyote", ]
points(coyote[,5] ~ coyote[,1], type = "l", lwd = 3.5, col = "darkviolet")
points(coyote[,5] ~ coyote[,2], type = "l", lwd = 6.5, col = "darkviolet")
points(coyote[,5] ~ coyote[,3], type = "p", cex = 1.75, pch = 19, col = 0)
points(coyote[,5] ~ coyote[,3], type = "p", cex = 1.75, pch = 10, col = "darkviolet", lwd = 2)

abline(v = 0, lty = 3)
mtext(expression(paste("Persistence (", Phi, ") Intercept", sep = "")), 3, line = 0.4, cex = 1.4)

dev.off()

mod_mat <- mod_matrix[, 30:34]

# get credible intervals
# get quantiles
creds <- apply(mod_mat, 2, quantile, probs = c(0.025, 0.975, 0.25,  0.75, 0.50))
# change column anmes
colnames(creds) <- c("Coyote", "Opossum", "Raccoon", "Red Fox", "Striped Skunk")
# make long format
creds <- melt(creds)
# make wide format
creds <- dcast(creds, Var2~ Var1)
# change column names
colnames(creds) <- c("species", "ci", "uci", "ci2", "uci2", "mode")
# add some more mode for the matrix
creds$mode2 <- creds$mode
# remove species
creds <- as.matrix(creds[,-1])
# make it a matrix without species
creds <- data.frame(matrix(creds, ncol = 3 ))
# change coumn names
colnames(creds) <- c("ci", "ci2", "mode")
# add a species column
creds$species <- factor(rep(c("Coyote", "Opossum", "Raccoon", "Red Fox", "Striped Skunk"), 2),
                        levels = c("Striped Skunk", "Red Fox", "Opossum", 
                                   "Raccoon", "Coyote"))



creds$sp <- as.numeric(creds$species)

tiff("detection_cater.tif", height = 3.5, width = 6.3, units = "in", res = 400)
par(mar=c(3.5,7.0,2.5,1))
plot(1~1, type='n', xlim=c(0,0.25),ylim = c(.8,5.2), xlab="", ylab="", xaxt="n", yaxt="n")
axis(1, at = seq(0, 0.25, .05), labels=F, tck=-.035)
axis(1, at=seq(0, 0.25, .025), labels=F, tck=-.02)
axis(1, at=seq(0, 0.25, .0125), labels=F, tck=-.01)
mtext(text=seq(0, 0.25, .05), 1, line=1, at=seq(0, 0.25, .05), cex=1.25)

axis(2, at = seq(1, 5, 1), labels=F, tck=-.035)
mtext(text = c("Skunk", "Red Fox", "Opossum", "Raccoon", "Coyote"),
      2, line = .75, at = seq(1, 5, 1), las = 1, cex = 1.5)

skunk <- creds[creds$species == "Striped Skunk", ]
points(skunk[,5] ~ skunk[,1], type = "l", lwd = 3.5)
points(skunk[,5] ~ skunk[,2], type = "l", lwd = 6.5)
points(skunk[,5] ~ skunk[,3], type = "p", cex = 1.75, pch = 19, col=0)
points(skunk[,5] ~ skunk[,3], type = "p", cex = 1.75, pch = 10, lwd=2)

fox <- creds[creds$species == "Red Fox", ]
points(fox[,5] ~ fox[,1], type = "l", lwd = 3.5, col = "olivedrab4")
points(fox[,5] ~ fox[,2], type = "l", lwd = 6.5, col = "olivedrab4")
points(fox[,5] ~ fox[,3], type = "p", cex = 1.75, pch = 19, col =0)
points(fox[,5] ~ fox[,3], type = "p", cex = 1.75, pch = 10, col ="olivedrab4", lwd = 2)

possum <- creds[creds$species == "Opossum", ]
points(possum[,5] ~ possum[,1], type = "l", lwd = 3.5, col = "dodgerblue4")
points(possum[,5] ~ possum[,2], type = "l", lwd = 6.5, col = "dodgerblue4")
points(possum[,5] ~ possum[,3], type = "p", cex = 1.75, pch = 19, col = 0)
points(possum[,5] ~ possum[,3], type = "p", cex = 1.75, pch = 10, col = "dodgerblue4", lwd = 2)

raccoon <- creds[creds$species == "Raccoon", ]
points(raccoon[,5] ~ raccoon[,1], type = "l", lwd = 3.5, col = "firebrick2")
points(raccoon[,5] ~ raccoon[,2], type = "l", lwd = 6.5, col = "firebrick2")
points(raccoon[,5] ~ raccoon[,3], type = "p", cex = 1.75, pch = 19, col = 0)
points(raccoon[,5] ~ raccoon[,3], type = "p", cex = 1.75, pch = 10, col = "firebrick2", lwd = 2)

coyote <- creds[creds$species == "Coyote", ]
points(coyote[,5] ~ coyote[,1], type = "l", lwd = 3.5, col = "darkviolet")
points(coyote[,5] ~ coyote[,2], type = "l", lwd = 6.5, col = "darkviolet")
points(coyote[,5] ~ coyote[,3], type = "p", cex = 1.75, pch = 19, col = 0)
points(coyote[,5] ~ coyote[,3], type = "p", cex = 1.75, pch = 10, col = "darkviolet", lwd = 2)


mtext(substitute(paste("Daily Probability (", italic("p"), ") of Detection", sep = "")), 3, line = 0.4, cex = 1.4)

dev.off()









 windows(height = 30, width = 30)
par(fig=c(.8,1,0.5,1))
plot(rnorm(100))
par(fig=c(0,1,0.15,.61), new=T)
plot(rnorm(100))

ggplot(sp_rich, aes(x = known, y = mode, fill = day_samp, shape = day_samp, size = day_samp))+
  geom_point(position = "jitter", color = "black", alpha = 0.5)+
  scale_shape_manual(values = c(21:25))+
  scale_size_manual(values = c(2, 2.25, 3, 3.5))
  geom_point(aes(x = known, y = uci), alpha = 0.5, color = "black", position = "jitter",
             size = 2, shape = 1)

ggplot(sp_rich, aes(x = day_samp_num, y = mode-known))+
  geom_point(position = "jitter")+
  scale_x_continuous(breaks = seq(1, 28, 1))+
  stat_smooth(se = FALSE)

sp_rich <- sp_rich[-which(sp_rich$day_samp_num == 0),]


tiff("est_vs_known.tif", height = 4.0 , width = 6.3, units = "in", res = 400)


par(mar=c(4.5,4.5,1.5,1.5))
plot(1~1, type='n', xlim=c(.5,28.5),ylim = c(-0.2,3.2), xlab="", ylab="", xaxt="n", yaxt="n")
axis(1, at = c(1,seq(5, 25, 5)), labels=F, tck=-.035)
axis(1, at = seq(1, 28, 1), labels = F, tck = -0.02)
mtext(text=c(1,seq(5, 25, 5)), 1, line=1, at=c(1,seq(5, 25, 5)), cex=1.25)
mtext(text = "Days Sampled", 1, line = 3, cex = 1.5)
axis(2, at = seq(0, 3, 1), label = F, tck = -0.035)
mtext(text = "Predicted - Observed", 2, line = 3, cex = 1.5)
mtext(text = seq(0, 3, 1), 2, line = 0.75, at = seq(0,3,1), cex = 1.25, las = 1)
mtext(text = "Species Richness", 3, line = 0.2, cex = 1.5)
points(jitter(mode - known) ~ jitter(day_samp_num, amount = .5), data = sp_rich, col=alpha("dodgerblue4", 0.4))
lines(smooth.spline((sp_rich$mode - sp_rich$known)~ sp_rich$day_samp_num, spar = .6), lwd = 2)
dev.off()






