library(RSQLite)
setwd("/home/rytis/Downloads/downloads_birth/preprocessed")


# Fetch aggregate tables from the SQLite database
con <- dbConnect(SQLite(), "data.db")

rs <- dbSendQuery(con, "SELECT COUNT(*),age,year FROM data GROUP BY age,year ORDER BY year")
total_births <- dbFetch(rs, -1)

rs2 <- dbSendQuery(con, "SELECT COUNT(*),age,year FROM data WHERE children > 1 GROUP BY age,year ORDER BY year")
twin_births <- dbFetch(rs2, -1)

dbClearResult(rs)
dbClearResult(rs2)

dbDisconnect(con)

# Data treatment
library(plyr)
library(ggplot2)
library(gridExtra)

names(total_births) <- names(twin_births) <- c("count", "age", "year")

total_births <- ddply(total_births, .(year), mutate, density=count/sum(count))
twin_births <- ddply(twin_births, .(year), mutate, density=count/sum(count))



plotf <- function(x) {
  p <- ggplot(x, aes(x=age, y=density, col=as.factor(year))) +
          guides(colour=FALSE) + geom_line() +
          xlab("Age") + ylab("Fertility density") +
          theme(text = element_text(size=10))
  return(p)
}

plot1 <- plotf(total_births)
plot2 <- plotf(twin_births)
grid.arrange(plot1, plot2, ncol=2)
