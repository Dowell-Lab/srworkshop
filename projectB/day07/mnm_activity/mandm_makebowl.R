library(tidyverse)

outdir = "/Users/maryallen/srworkshop/projectB/day07/mnm_activity/"
#you set these
totalMandMs = 38*32
skittlestotal = (27*16)/2
mike_ike = 16*4 
redhot = 25*10 

mike_ikepercent = 1
redhotpercent = 1

#percent M and M 
mmbluepercent = 0.24
mmorangepercent = 0.2
mmgreenpercent = 0.16
mmyellowpercent = 0.14
mmredpercent = 0.13
mmbrownpercent = 0.13

#percent skittles
sktlorangepecent = 0.20
sktlgreenpecent = 0.20
sktlredpecent = 0.20
sktlyellowpecent = 0.20
sktlpurplepecent = 0.20

colors = c("mmblue", "mmorange", "mmgreen", "mmyellow",
           "mmred", "mmbrown", "mike_ikepercent", "redhotpercent",
           "sktlorange", "sktlgreen", "sktlred", 
           "sktlyellow", "sktlpurple")

percents_item = rbind(mmbluepercent,mmorangepercent,mmgreenpercent,mmyellowpercent, 
              mmredpercent, mmbrownpercent, mike_ikepercent, redhotpercent,
                sktlorangepecent, sktlgreenpecent, sktlredpecent, 
                sktlyellowpecent, sktlpurplepecent)
totals = c(rep(totalMandMs, 6),rep(mike_ike, 1),rep(redhot, 1) , rep(skittlestotal, 5))

df = cbind(colors, percents_item, totals)

df <- data.frame(df)

df$percents_item <- df$V2
df$percents_item <- as.numeric(df$percents_item)
df$totals <- as.numeric(df$totals)


createbowl <- function(df, extraname, extranumber){
  df2 <- df
  df2$bowlnumber <- df2$totals*df2$percents_item
  df2$bowlnumber <- as.integer(df2$bowlnumber)
  df2$name <- df2$colors
  #df2$frac_of_bowl = df2$bowlnumber/sum(df2$bowlnumber) two years ago all I did was point out that it's already written so if it's that tversions too old so we had them W get an older version and thou it wor
  return (df2)
}

add_to_bowl <- function(df, extraname, extranumber) {
  df$bowlnumber[df$name == extraname] <- df$bowlnumber[df$name == extraname] + extranumber
  return(df)
}

equalbowl = createbowl(df)
redbowl = add_to_bowl(equalbowl, extraname="mmblue", extranumber=0)
greenbowl = add_to_bowl(equalbowl, extraname="mmblue", extranumber=230)

justcountsred = redbowl %>% select(name, bowlnumber)
colnames(justcountsred) <- c("name", "redbowl")
justcountsgreen = greenbowl %>% select(name, bowlnumber)
colnames(justcountsgreen) <- c("name", "greenbowl")

bowlsdf= merge(justcountsred, justcountsgreen, by="name")

write.csv(bowlsdf, paste(outdir,"bowldf.csv", sep=""))

