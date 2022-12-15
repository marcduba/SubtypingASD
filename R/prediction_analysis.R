### HYPOTHESIS 3
### GLM ANALYSIS: PREDICTION OF DIAGNOSIS

library(dplyr)
library(tidyr)
library(stringr)
library(jtools)
library(caret)
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(pROC)
library(cvAUC)
source('~/Documents/Psychologie/MA/Code/R/GeneralParameters.R')
load('~/Documents/Psychologie/MA/Code/R/SignificantRegions.Rda')

# Load and prepare data
setwd('~/Desktop/Clustering2/Data/')
allData <- read.csv('allData.csv')

# separating data into variables that are needed for analysis
VarNames <- colnames(allData[2:152])
VarNames <- c(VarNames, params$BehavMeasures[5], params$BehavMeasures[6])
braindata <- allData[, VarNames]
colnames(braindata) <- colnames(braindata) %>%
  str_replace_all('[.]', '-') %>%
  str_replace_all('Basic_Demos-', '')

significantRegions$CT
regionnamesCT <- significantRegions$CT$ASD1v2
regionnamesCT <- paste(regionnamesCT, '_CT', sep = '')

significantRegions$SA
regionnamesSA <- significantRegions$SA$ASD1v2
regionnamesSA <- paste(regionnamesSA, '_SA', sep = '')

# filtering data sets for regions and clusters
datGLM1 <- braindata %>%
  dplyr::select(c(diag, clust, Sex, Age, matches(regionnamesCT), matches(regionnamesSA))) %>%
  filter(clust == 1 | clust == 0)
datGLM1 <- datGLM1 %>%
  dplyr::select(-c(clust))
datGLM1$diag <- as.factor(datGLM1$diag)
datGLM1$Sex <- as.factor(datGLM1$Sex)
nyholt1 <- getMeff(datGLM1[,3:6])

datGLM2 <- braindata %>%
  dplyr::select(c(diag, clust, Sex, Age, matches(regionnamesCT), matches(regionnamesSA))) %>%
  filter(clust == 2 | clust == 0)
datGLM2 <- datGLM2 %>%
  dplyr::select(-c(clust))
datGLM2$diag <- as.factor(datGLM2$diag)
datGLM2$Sex <- as.factor(datGLM2$Sex)
nyholt2 <- getMeff(datGLM2[,3:6])

datGLM3 <- braindata %>%
  dplyr::select(c(diag, clust, Sex, Age, matches(regionnamesCT), matches(regionnamesSA))) %>%
  filter(clust == 1 | clust == 2 | clust == 0)
datGLM3 <- datGLM3 %>%
  dplyr::select(-c(clust))
datGLM3$diag <- as.factor(datGLM3$diag)
datGLM3$Sex <- as.factor(datGLM3$Sex)
nyholt3 <- getMeff(datGLM3[,3:6])

# spliting data into train (80%) and test (20%) set 
set.seed(123)
index1 <- sample(nrow(datGLM1), nrow(datGLM1)*0.8)
dat1_train <- datGLM1[index1, ]
dat1_test <- datGLM1[-index1, ]

index2 <- sample(nrow(datGLM2), nrow(datGLM2)*0.8)
dat2_train <- datGLM2[index2, ]
dat2_test <- datGLM2[-index2, ]

index3 <- sample(nrow(datGLM3), nrow(datGLM3)*0.8)
dat3_train <- datGLM3[index3, ]
dat3_test <- datGLM3[-index3, ]

## GLM analysis
# Cluster 1 v controls
set.seed(456)
train.control <- trainControl(method = 'cv', number = 10)
model1 <- train(diag ~ . + Age:. + Sex:., 
                data = dat1_train, 
                method = 'glm', 
                trControl = train.control, 
                family = binomial())

model2 <- train(diag ~ . + Age:. + Sex:., 
                data = dat2_train, 
                method = 'glm', 
                trControl = train.control, 
                family = binomial())

model3 <- train(diag ~ . + Age:. + Sex:., 
                data = dat3_train, 
                method = 'glm', 
                trControl = train.control, 
                family = binomial())

summary(model1)
summary(model2)
summary(model3)

# function to display glm results
regression_table <- function(glm, nyholt) {
  coefficients <- data.frame(summary(glm)$coefficients)
  coefficients$p.cor <- coefficients$Pr...z.. * nyholt
  significant <- data.frame(coefficients[which(coefficients$p.cor <= 0.05),])
  
  # significance stars
  sig.stars <- ifelse(coefficients$p.cor < 0.0001, '***', 
                      ifelse(coefficients$p.cor < 0.001, '*** ', 
                             ifelse(coefficients$p.cor < 0.01, '**  ', 
                                    ifelse(coefficients$p.cor < 0.05, '*   ', '    '))))
  # display coefficients with stars
  coefficients.stars <- coefficients
  coefficients.stars[1:3] <- format(round(coefficients.stars[1:3], digits = 2), nsmall = 2)
  coefficients.stars[4:5] <- format(round(coefficients.stars[4:5], digits = 3), nsmall = 3)
  coefficients.stars$p.cor <- paste(coefficients.stars$p.cor, sig.stars, sep = '')
  colnames(coefficients.stars) <- c('Estimate', 'SE', 'z.value', 'p', 'p corr')
  rownames(coefficients.stars) <- rownames(coefficients.stars) %>%
    str_replace_all('`', '') %>%
    str_replace_all(fixed('\\'), '') %>%
    str_replace_all('[:]', ' * ') %>%
    str_replace_all('1', '(fem)')
  
  results <- list(coef = coefficients, 
                  sig = significant, 
                  coef.stars = coefficients.stars, 
                  null = summary(glm)$null, 
                  dev = summary(glm)$deviance, 
                  aic = glm$finalModel$aic, 
                  acc = glm$results$Accuracy,
                  accSD = glm$results$AccuracySD,
                  kappa = glm$results$Kappa,
                  kappaSD = glm$results$KappaSD)
  return(results)
}

# model results
regTab_m1 <- regression_table(model1, nyholt1)
print(regTab_m1$coef.stars)
print(regTab_m1$aic)
print(regTab_m1$acc)
print(regTab_m1$accSD)
print(regTab_m1$kappa)
print(regTab_m1$kappaSD)
print(regTab_m1$null)
print(regTab_m1$dev)

regTab_m2 <- regression_table(model2, nyholt2)
print(regTab_m2$coef.stars)
print(regTab_m2$aic)
print(regTab_m2$acc)
print(regTab_m2$accSD)
print(regTab_m2$kappa)
print(regTab_m2$kappaSD)
print(regTab_m2$null)
print(regTab_m2$dev)

regTab_m3 <- regression_table(model3, nyholt3)
print(regTab_m3$coef.stars)
print(regTab_m3$aic)
print(regTab_m3$acc)
print(regTab_m3$accSD)
print(regTab_m3$kappa)
print(regTab_m3$kappaSD)
print(regTab_m3$null)
print(regTab_m3$dev)


# saving model output
write.table(regTab_m1$coef.stars, file = '~/Documents/Psychologie/MA/Results/Prediction/model1.txt', row.names = T, col.names = T, quote = F, sep = '\t', dec = '.')
write.table(regTab_m2$coef.stars, file = '~/Documents/Psychologie/MA/Results/Prediction/model2.txt', row.names = T, col.names = T, quote = F, sep = '\t', dec = '.')
write.table(regTab_m3$coef.stars, file = '~/Documents/Psychologie/MA/Results/Prediction/model3.txt', row.names = T, col.names = T, quote = F, sep = '\t', dec = '.')

# table for model comparison
models_table <- data.frame(Measure = c('Accuracy', 'SD', 'Kappa', 'SD', 'AIC'), 
                           Model_1 = c(round(regTab_m1$acc, digits = 4), round(regTab_m1$accSD, digits = 2), round(regTab_m1$kappa, digits = 4), round(regTab_m1$kappaSD, digits = 2), round(regTab_m1$aic, digits = 2)), 
                           Model_2 = c(round(regTab_m2$acc, digits = 4), round(regTab_m2$accSD, digits = 2), round(regTab_m2$kappa, digits = 4), round(regTab_m2$kappaSD, digits = 2), round(regTab_m2$aic, digits = 2)), 
                           Model_3 = c(round(regTab_m3$acc, digits = 4), round(regTab_m3$accSD, digits = 2), round(regTab_m3$kappa, digits = 4), round(regTab_m3$kappaSD, digits = 2), round(regTab_m3$aic, digits = 2)))

colnames(models_table) <- colnames(models_table) %>%
  str_replace_all('_', ' ')

write.table(models_table, file = '~/Documents/Psychologie/MA/Results/Prediction/summary_table.txt', row.names = F, col.names = T, quote = F, sep = '\t', dec = '.')


## Comparing model performance

# Prediction in test data
predict1 <- predict(model1, newdata = dat1_test, type = 'prob')
predict2 <- predict(model2, newdata = dat2_test, type = 'prob')
predict3 <- predict(model3, newdata = dat3_test, type = 'prob')

# ROC
roc1 <- roc(response = dat1_test$diag, predictor = predict1$ASD, auc = T, ci = T, plot = T)
roc2 <- roc(response = dat2_test$diag, predictor = predict2$ASD, auc = T, ci = T, plot = T)
roc3 <- roc(response = dat3_test$diag, predictor = predict3$ASD, auc = T, ci = T, plot = T)

# AUC
auc1 <- auc(roc1)
auc2 <- auc(roc2)
auc3 <- auc(roc3)

ci.auc1 <- ci.auc(roc1, conf.level = .95, method = 'bootstrap', boot.n = 5000)
ci.auc2 <- ci.auc(roc2, conf.level = .95, method = 'bootstrap', boot.n = 5000)
ci.auc3 <- ci.auc(roc3, conf.level = .95, method = 'bootstrap', boot.n = 5000)

ci.se1 <- ci.se(roc1, boot.n = 5000)
ci.se2 <- ci.se(roc2, boot.n = 5000)
ci.se3 <- ci.se(roc3, boot.n = 5000)

cidat1 <- data.frame(x = rev(as.numeric(rownames(ci.se1))), lower = ci.se1[,1], upper = ci.se1[,3])
cidat2 <- data.frame(x = rev(as.numeric(rownames(ci.se2))), lower = ci.se2[,1], upper = ci.se2[,3])
cidat3 <- data.frame(x = rev(as.numeric(rownames(ci.se3))), lower = ci.se3[,1], upper = ci.se3[,3])

rocplot_function <- function(title, roc, auc, ciauc, cidat, color) {
  label1 <- paste('AUC:', format(round(auc, 4), nsmall = 4), sep = ' ')
  label2 <- paste('[', format(round(ciauc[1], 4), nsmall = 4), '-', format(round(ciauc[3], 4), nsmall = 4), ']', sep = '')
  
  plot <- ggroc(roc, legacy.axes = T, col = color, linetype = 1, size = .75) +
    geom_abline(intercept = 0, slope = 1, col = 'black', lty = 'dashed', alpha = .75) +
    geom_ribbon(data = cidat, aes(x = x, ymin = lower, ymax = upper), fill = color, alpha = .1) + 
    theme_apa() + 
    ggtitle(label = '', subtitle = title) +
    annotate('text', x = .8, y = .2, label = label1) + 
    annotate('text', x = .8, y = .1, label = label2)
  
  return(plot)
}

plot1 <- rocplot_function('Model 1', roc1, auc1, ci.auc1, cidat1, params$paletteASD[1])
plot2 <- rocplot_function('Model 2', roc2, auc2, ci.auc2, cidat2, params$paletteASD[3])
plot3 <- rocplot_function('Model 3', roc3, auc3, ci.auc3, cidat3, params$paletteASD[2])

print(plot1)
print(plot2)
print(plot3)

ggsave(plot = plot1, filename = '~/Documents/Psychologie/MA/Results/Prediction/Figures/ROC_1.png', width = 600, height = 500, units = 'px', scale = 2.5)
ggsave(plot = plot2, filename = '~/Documents/Psychologie/MA/Results/Prediction/Figures/ROC_2.png', width = 600, height = 500, units = 'px', scale = 2.5)
ggsave(plot = plot3, filename = '~/Documents/Psychologie/MA/Results/Prediction/Figures/ROC_3.png', width = 600, height = 500, units = 'px', scale = 2.5)

allROC <- grid.arrange(plot1, plot2, plot3, nrow = 2, ncol = 2)

ggsave(allROC, filename = '~/Documents/Psychologie/MA/Results/Prediction/Figures/allROC.png', width = 700, height = 600, units = 'px', scale = 3)

# testing for difference in AUCs
roc.test(roc1, roc2, method = 'bootstrap', boot.n = 5000, conf.level = .95)
roc.test(roc1, roc3, method = 'bootstrap', boot.n = 5000, conf.level = .95)
roc.test(roc2, roc3, method = 'bootstrap', boot.n = 5000, conf.level = .95)

