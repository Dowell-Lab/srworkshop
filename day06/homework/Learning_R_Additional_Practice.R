####### Learning R Additional Practice
####### Author: Taylor Jones (2022), Georgia Barone (2023)
####### Here we will practice what we learned today in class. This will be less commented than the work in class.

####### You will be using a provided dataframe called: pokemon_data.csv
####### pokemon_data.csv is on the sr2023 GitHub under day06/data.
####### Download this data to your local machine before beginning this exercise. 
workdir <- '/path/to/your/working/dir/'
setwd(workdir)
data <- read_csv('/path/to/pokemon_data.csv')

####### 1. Summary Stats
#1. do a quick summary and head to get a sense for the dataset
#2. what are the column names?
#3. how many of each "Type 1" powers are there?
#4. how many legendary pokemon?
#5. how many generations are there?=
#6. average HP (health)? attack? defense?
#7. stdev of these things?
#8. what are the dimensions of this dataframe?

####### 2. Generate plots- basic
#1. plot attack vs defense
#2. plot various stats vs generation
#3. do boxplots for all stats (hint: make a list of columns you care about, parse the dataframe then run boxplot(data))
#4. plot histograms of speed, attack, and defense
#5. what is the highest HP (health) pokemon? highest attack? highest defense?
#6. what other trends do you notice?

####### 3. Manipulating dataframes- basic
#1. column names with spaces can cause issues when manipulating dataframes.
#rename the columns so there are no spaces to make life easier.
#2. create a dataframe with only pokemons that have Fire powers in the Type 1 column.
#3. what is the max value in the Defense column for this new dataframe?
#4. add a new column to the original dataframe that is the ratio of sp.atk vs speed
#5. what it the average ratio of sp. atk vs speed?
#6. now delete your new ratio column to restore the original dataframe
#7. create a new df that contains only the Names, Type 1, and Type 2 columns

####### 4. Manipulating dataframes- advanced (remember to load in tidyverse)
#1. How are type 1 and type 2 related? What is the most common combination?
#2. Filter for the top quartile attack OR top quartile defense pokemon
#3. Filter for the top quartile Speed AND generation 3 or greater
#4. how pokemon are generation 4? How many are NOT generation 4?
#5. how pokemon are have a name that start with "B"?
#6. how many different types are there and which is the most common? Which is the rarest?

####### 5. Generate plots- advanced
#1. find a pair of stats that are interesting and make a beautiful pokemon plot!
#note** consider all the possible additions through ggplot that we practiced on iris and mtcars

####### 6. Install R packages rsubread and DESeq2 (we will be using these packages tomorrow)
# Both packages can be found on bioconductor
# rsubread (https://bioconductor.org/packages/release/bioc/html/Rsubread.html)
# DESeq2 (https://bioconductor.org/packages/release/bioc/html/DESeq2.html)
