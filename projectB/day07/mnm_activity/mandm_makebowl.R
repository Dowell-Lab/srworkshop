
#you set these
totalMandMs = 10000
skittlestotal = 100

#percent M and M spike ins do you want
bigmarsh = 0.05
tinymarsh = 0.25

mmbluepercent = 0.24
mmorangepercent = 0.2
mmgreenpercent = 0.16
mmyellowpercent = 0.14
mmredpercent = 0.13
mmbrownpercent = 0.13

sktlorangepecent = 0.20
sktlgreenpecent = 0.20
sktlredpecent = 0.20
sktlyellowpecent = 0.20
sktlpurplepecent = 0.20

colors = c("mmblue", "mmorange", "mmgreen", "mmyellow",
           "mmgreen", "mmyellow", "mmbrown", "bigmarsh", "tinymarsh",
           "sktlorangepecent", "sktlgreenpecent", "sktlredpecent", 
           "sktlyellowpecent", "sktlpurplepecent")

percents_item = rbind(mmbluepercent,mmorangepercent,mmgreenpercent,mmyellowpercent, 
                 mmgreenpercent, mmredpercent, mmbrownpercent, bigmarsh, tinymarsh,
                sktlorangepecent, sktlgreenpecent, sktlredpecent, 
                sktlyellowpecent, sktlpurplepecent)
totals = c(rep(totalMandMs, 9), rep(skittlestotal, 5))

df = cbind(colors, percents_item, totals)

df <- data.frame(df)

df$percents_item <- df$V2