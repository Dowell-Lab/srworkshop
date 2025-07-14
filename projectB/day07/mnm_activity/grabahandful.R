
#indir="/Users/maryallen/srworkshop/projectB/day07/mnm_activity/"
indir=""
infile = "bowldf.csv"
howmanyreads = 100
whichbowl = "redbowl" #other option is greenbowl


#don't edit stuff below

grabahandful <- function(whichbowl, howmanyreads){
  bowls = read.csv(paste(indir, infile, sep=""))
  if (whichbowl=="redbowl"){
    samp_idx <- sample(seq_len(nrow(bowls)), howmanyreads, prob = df$redbowl, replace = TRUE)
    
  } else {
    samp_idx <- sample(seq_len(nrow(bowls)), howmanyreads, prob = df$greenbowl, replace = TRUE)
  }
  new_data <- bowls[samp_idx, ] 
  handful = table(new_data$name)
  return(handful)}

handful <- grabahandful(whichbowl, howmanyreads)
print(handful)

