library(rJava)
library(UpSetR)
library(tidyverse)
library(venneuler)
library(grid)

dat <- read.table('/media/hamreni/19ff7812-0544-4d55-bb1c-75ed0293e108/home/hamreni/Documents/snorna_project/predictors_summary.csv')
dat
dat %>% print(n = nrow(dat))
subsets <- dat$combination
subsets

haca <- c(snoReport = 47, snoGPS = 563, cmsearch = 232, `snoReport&snoGPS` = 39, `snoReport&cmsearch` = 17, 
            `snoGPS&cmsearch` = 151, `snoReport&snoGPS&cmsearch` = 40)
cd <- c(snoReport = 16, snoSCAN = 30, cmsearch = 811, `snoReport&snoSCAN` = 3, `snoReport&cmsearch` = 38, 
          `snoSCAN&cmsearch` = 18, `snoReport&snoSCAN&cmsearch` = 28)
  
# Create an UpsetR Plot
upset(fromExpression(haca), order.by = "freq", 
      #sets.bar.color=c("maroon","blue","orange"), 
      point.size=5,
      queries = list(
        list(
          query = intersects,
          params = list("snoReport", "snoGPS", "cmsearch"), 
          color = "#Df5286", 
          active = T,
          query.name = "Predicted by all"
        )
      ),
      mainbar.y.label = "Number of predicted HACA snoRNAs"
      )
upset(fromExpression(cd), order.by = "freq", 
      #sets.bar.color=c("maroon","blue","orange"), 
      point.size=5,
      queries = list(
        list(
          query = intersects,
          params = list("snoReport", "snoSCAN", "cmsearch"), 
          color = "#Df5286", 
          active = T,
          query.name = "Predicted by all"
        )
      ),
      mainbar.y.label = "Number of predicted CD snoRNAs"
)
