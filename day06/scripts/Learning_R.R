####### Learning R
####### Author: Taylor Jones (2022), Georgia Barone (2023)
####### Here we will learn a little more about R. R is built as a statistical program,
####### so we will play with some dataframes, statistics and plotting to get used to R.
####### The top portion IN THIS PANEL is the script. This is a great place to record what functions you used.
####### The bottom panel is called your Console. That is where your output will be located.
####### To the right is your variable list under "Environment" and your output for plots under "Plots".
####### The pound/hashtag sign used here allows you to comment out code so it won't run.
####### It is excellent practice to write out what you are doing in your code with a "#" and note what you're doing.
####### To run a line you can either hit Ctrl+Enter OR the "Run" button to the top right of this cell.

##############################################################################################################################
# Let's begin by setting a working directory, where the magic happens.
# The command for this is setwd()
# **Side note: make sure to edit the path below to match your desired working directory path 
setwd('/path/to/your/working/dir/')
getwd()

# We can also do this through a variable. This can be helpful if we want to output something specifically to our working directory
workdir <- '/path/to/your/working/dir/'
setwd(workdir)
getwd()

# We are going to pull in data that is part of base R as practice.
#So first we grab the data
data(iris) #grabbing the data

####### NOTE: To pull in a csv file you set the data as a variable: data<-read.csv(file='FILE.csv')
####### You can think of a dataframe almost like a table in excel.

#first let's do a quick check on what is in iris
names(iris) #names gives us column names inside of iris.
#like before we can do sum, mean and sd, but this time on a dataframe!
#to do this we need to call the dataframe then the column like so: dataframe_name$column_name
sum(iris$Sepal.Length) #sum
mean(iris$Sepal.Length) #mean
sd(iris$Sepal.Length) #stdev

# We can also look at more general information about iris by using the following commands, summary, names, dim
head(iris) #see the first few rows for iris
summary(iris) #summary gives min, 1st quart, mean, 3rd quart, max for each column
dim(iris) #dimensions of dataframe


####### 1. Get general information about mtcars.
#NOW try the same thing for a new dataset called mtcars. Check out the sum, mean and standard deviation for mpg.
#grab mtcars
#Note: now mtcars should be pulled into this environment. You should see it listed under the "Environment" tab to the right.

#what columns are in mtcars?

#sum of mpg in mtcars

#mean of mpg in mtcars

#stdev of mpg in mtcars

#check the first few lines of mtcars

#check the dimensions of mtcars

#summarize the data in mtcars

##############################################################################################################################
#Now that we have checked out our data we can plot it!
#The base plots aren't beautiful but are quick and easy ways to assess a large amount of data.

#######R for Data Science says it well: Data exploration is the art of looking at your data, 
#######rapidly generating hypotheses, 
#######quickly testing them, 
#######then repeating again and again and again.
#######The goal of data exploration is to generate many promising leads that you can later explore in more depth.

#Here is a quick scatter plot
plot(x=iris$Petal.Length, y=iris$Petal.Width)
#This is correlated such that the longer the petal, the wider the petal

#Here is a quick histogram
hist(iris$Petal.Length)
#This looks like it has a bimodal distribution. The second part of the distribution looks normal.

#Now lets make a boxplot
boxplot(iris)

#This is a nice box plot. Let's save it.
png(paste0(workdir, 'iris_box_plot.png')) #here we are "pasting" our working directory and our plot name together for our output file
boxplot(iris) #here we generate the plot
dev.off() #here we tell R "ok, I'm done"
#Now check your folder and see if it's there! A common error is missing the final "/" in the workdir

#Perhaps you want to make a special output for your scatter plot comparisons
#We want to go back to our petal length vs width and make a special outdir for it.

condx<-'Petal.Length' #we can turn our column names into variables of x and y for easy swapping if desired
condy<-'Petal.Width'

outdir <- paste(workdir, condy, '_vs_', condx, '/', sep='') ##Generic code that will allow us to create a SPECIFIC folder
dir.create(outdir) ###creating the directory

png(paste0(outdir,'iris_scatter.png'))
plot(x=iris[[condx]], y=iris[[condy]])  # Since condx and condy haven't changed we are SURE that the data matches the folder
#The square brackets allow us to call a column name as a string
dev.off()

####### 2. Generate a plots for mtcars.
#plot mpg vs weight (wt)

#what trends do you notice?

#plot a histogram for mpg

#what do you notice?

#Make a boxplot for mtcars and save it to your working directory

##############################################################################################################################
# Now let's practice working with and manipulating the dataframes.

#Let's start by looking at only long petal flowers
long_iris <- subset(iris, Petal.Length > 4.5) #you can see subset is a slightly more complicated command. It takes two inputs separated by a comma
dim(long_iris) #check dimensions
hist(long_iris$Petal.Length) #our new histogram
#notice it is now just the tail of the normal distribution

#What if we want to know the difference between petal length and width?
iris$petal_len_minus_width = iris$Petal.Length - iris$Petal.Width #this will subtract width from length across a given row
head(iris) #we now have a new column
hist(iris$petal_len_minus_width)
#you can see the histogram looks very similar to the petal length

#We decide to get rid of our new column since it isn't very informative.
names(iris) #we see that the new column is the 6th column of iris
iris <- iris[-c(6)] #this deletes the 6th element of the table.
head(iris) #now when we check the dataframe that column is gone.

#Something that is VERY useful is being able to rename columns
#Your count table from featureCounts may have some snarly names
#You'll want some easier names to call before you go into differential expression analysis
names(iris) <- c("sepal_len", "sepal_wid", "petal_len",
                 "petal_wid", "species")
#Check out the new names
head(iris)

#we can make a list of column names using the names command also
names_iris <- names(iris)
#we can then SELECT names we care about from this list, say we only care about petals
names_iris
names_iris[c(3:4)]
#we then can filter our dataframe based on this list
iris_petals <- subset(iris, select = names_iris[c(3:4)])
head(iris_petals)

#transpose iris
tiris <- t(iris)
head(tiris)

####### 3. Manipulating mtcars.
#check only cars with good mpg and create a new table called enviro_friendly

#you want to look at how mpg and cylinders may be related. You decide to make a new column that takes the RATIO of mpg to cylinder

#look at the data and see if there are any interesting trends (hint: checking a scatter between mpg and mpg/cyl might be interesting)

#now delete your new ratio column to restore mtcars to its original state

#rename the columns of mtcars to be more descriptive (like cyl is cylinder for example)

#select mpg and cylinder columns only

#transpose mtcars

##############################################################################################################################
# Package management

# you can see all installed packages (or called libraries in R) under the package tab to the right (Files/Plots/Packages)
# Take a second to look at what is installed for base R!
# now that we are fresh we can download our first package
# it is good code to save but hash out your code you use to install a package. Once it is installed you don't need to install it again!

# most packages you can install with this syntax: install.packages("PACKAGE"). Such as this package:
#install.packages("tidyverse") #--great package, uncomment this line and run it to install tidyverse on your local machine

library(tidyverse) #This loads tidyverse into our environment

??tidyverse #this calls for help of a package. This is useful if you have no idea where to start

##############################################################################################################################
# "Advanced" dataframe manipulation
# tidyverse contains dplyr which is a package that permits more advanced dataframe manipulations
??dplyr

##### Let's look at some logical operators (filter).
#Understanding how to use logical operators will help you assess your data and generate meaningful plots.
#When you get a differential expression table what if you want to look at significant genes?
#How about genes that are significant AND upregulated?
# we can filter with "and" using "&"

#We want sepal length greater than (or equal to) 6 AND petal length greater than (or equal to) 2 using dplyr
filter(iris, sepal_len >= 6 & petal_len >= 2)

#We only care about the versicolor irises from this list
filter(iris, sepal_len >= 6 & petal_len >= 2 & species == "versicolor")

#what if we want to make a plot looking at this new filtered data?
big_versicolor <- filter(iris, sepal_len >= 6 & petal_len >= 2 & species == "versicolor")
plot(x=big_versicolor$sepal_len, y=big_versicolor$petal_len)

#What if we want significant genes that have a log2(fold change) of MORE than 2 OR LESS than -2 (ie greater or less than a 2 fold change)
# we can filter with "or" using "|"

#We want petal length greater than 6 OR petal length less than 2
filter(iris, petal_len > 6 | petal_len < 2)

#what if we want to make a plot looking at this new filtered data?
sep_petal_len <- filter(iris, petal_len > 6 | petal_len < 2)
hist(sep_petal_len$petal_len)

#Say we want to ONLY look at the species "virginica"
filter(iris, species == "virginica")

#Say we want to DON'T want look at the species "virginica"
filter(iris, species != "virginica")

#Say we want species starting with v- we can use str_detect for more general filtering
summary(filter(iris, str_detect(species, "v")))

##### So a recap:
# & means and
# | means or
# < means less than (> greater)
# <= means less than or equal to (>= greater)
# == means equal to
# != means NOT equal to

##### Before we filtered rows, now we can select columns (select).

#select species column and list just the unique species
unique(select(iris, 'species'))

#select all but species- this drops the species column
select(iris, -'species')

##### We can add new variables (mutate)
#adding petal ratio as a column
iris <- mutate(iris, petal_ratio = petal_len/petal_wid)

#adding sepal ratio as a column
iris <- mutate(iris, sepal_ratio = sepal_len/sepal_wid)
head(iris)

#adding an id column
iris$index <- rownames(iris)
iris <- mutate(iris,id=as.character(paste0(index, '_',species, '_',petal_ratio, '_', sepal_ratio)))

#changing the index to the id
rownames(iris) <- iris$id
head(iris)

#changing it back because it is sort of long
rownames(iris) <- iris$index

##### We can order rows (arrange)
#arrange from smallest to largest petal ratio
head(arrange(iris, petal_ratio))
tail(arrange(iris, petal_ratio))

#to flip and go largest to smallest add "desc"
head(arrange(iris, desc(petal_ratio)))

##### We can use (group_by) to look at certain groups at a time
#this is useful for looking at differences in iris species
ir <- group_by(iris, species)
group_keys(ir) #listing how it was grouped

#we can then split this dataframe
group_split(ir)

####### 4. Tidy-ing mtcars
#plot a scatter plot of cars with a weight greater than (or equal to) 3.6 and greater than (or equal to) 15 mpg (wt vs mpg)

#plot a histogram of weights for cars with a weigh less than 1.5 or greater than 3.6

#plot a scatter plot for mpg for this filtered dataframe

#how many cars have exactly 4 cylinders? How many do NOT have 4 cylinders? NOTE: you can use count() to count how many rows are in a dataframe.

#how many cars have the same number of gears and carborators?

#make mtcars car names into a column (hint: rownames)

#filter for all Toyta and Mazdas

#group mtcars by cylinder and make a boxplot for gears

##############################################################################################################################
# "Advanced" Plotting
#quick R plots are really useful for exploratory data analysis. They are not good for your publications!
#So now we will generate some pretty plots! We can add dimensions like color to get more information from our dataframe
#let's start with sep_petal_len that we generated before
hist(sep_petal_len$petal_len)
plot(x=sep_petal_len$petal_len, y=sep_petal_len$petal_wid)

#we can use ggplot to improve this scatter plot
g<-ggplot(sep_petal_len, aes(petal_len, petal_wid))+geom_point()
g #run g to view the plot

#let's color the points based on the species- this will give us additional information
g<-ggplot(sep_petal_len, aes(petal_len, petal_wid, color=species))+geom_point()
g

#what if I want exact colors?
g<-ggplot(sep_petal_len, aes(petal_len, petal_wid, color=factor(species)))+
  scale_color_manual(values=c("#2540D8", "red")) + geom_point(shape = 16, size = 1.4)  #you can use color names such as red or hex codes to be more precise
g

#let's change the labels for the axes
#NOTE: g already exists so we can just add to it!
g<-g+ggtitle('Short and Long Petals')+
  labs(x='Petal Length', y='Petal Width')
g

#let's change change our font themes.
g <- g+theme(plot.title= element_text(size=20, face='bold', margin = margin(10, 0, 10, 0))) +
  theme(axis.text = element_text(size=12, face='bold')) +
  theme(axis.title = element_text(size=12, face = 'bold'))
g

#let's make out background white
g <- g+ theme(panel.background = element_rect(fill = 'white'),
              panel.grid.major = element_line(colour = "white",linewidth =0),
              panel.grid.minor = element_line(colour = "white", linewidth =0))
g

### or
g <- g+ theme_minimal()
g

#let's set the x and y axis limits
g <- g+xlim(c(1,7))+ylim(c(0,3))
g

#let's add a verticle line and a horizontal line that are the middle between the species
summary(sep_petal_len) #looking at our df
vir <- filter(sep_petal_len, species=='virginica') #pulling species 1
summary(vir) #looking at stats for species 1
set <- filter(sep_petal_len, species=='setosa') #pulling species 2
summary(set) #looking at stats for species 2
v_line <- (mean(vir$petal_len) + mean(set$petal_len))/2 #calculating where to put the vertical line
v_line #looking at that location
h_line <- (mean(vir$petal_wid) + mean(set$petal_wid))/2 #calculating where to put the horizontal line
h_line #looking at that location

#now we can actually plot those lines!
g <- g+ 
  geom_hline(yintercept=h_line, linetype="dashed", color = "black") + 
  geom_vline(xintercept=v_line, linetype="dashed", color = "black")
g

#let's add in species names ONTO the plot.
g <- g + 
  geom_text(x=5.8, y=1.5, label="virginica", color='red') +
  geom_text(x=2.1, y=0.8, label="setosa", color='#2540D8')

g

####Here is all we did in one big plot input!
g_combined <- ggplot(sep_petal_len, aes(petal_len, petal_wid, color=factor(species)))+
  scale_color_manual(values=c("#2540D8", "red")) + geom_point(shape = 16, size = 1.4)+
  ggtitle('Short and Long Petals')+
  labs(x='Petal Length', y='Petal Width')+
  theme(plot.title= element_text(size=20, face='bold', margin = margin(10, 0, 10, 0))) +
  theme(axis.text = element_text(size=12, face='bold')) +
  theme(axis.title = element_text(size=12, face = 'bold')) +
  theme_minimal() +
  xlim(c(1,7))+ylim(c(0,3))+ 
  geom_hline(yintercept=h_line, linetype="dashed", color = "black") + 
  geom_vline(xintercept=v_line, linetype="dashed", color = "black") +
  geom_text(x=5.8, y=1.5, label="virginica", color='red') +
  geom_text(x=2.1, y=0.8, label="setosa", color='#2540D8')
g_combined
  
#finally, let's save all our work on g!!
ggsave("pretty_short_vs_long_petals.png", plot = last_plot(), path = workdir,
       scale = 1, width = NA, height = NA, units = c("in", "cm", "mm"),
       dpi = 600, limitsize = TRUE)


####### 5. Pretty plotting mtcars
#Create a pretty plot for mtcars mpg vs weight
plot(mtcars$wt, mtcars$mpg)

#color dots by cylinder number
#change the dots to another shape (hint: google geom_point shape)
#change the theme to white/simple
#add titles and change x/y axis titles,
#make font bold
#set the x and y limits
#add a horizontal line for where "efficient" cars begin (judgement call- decide what makes a car efficient)
#label this new line on the plot "Efficient cars"
#save your plot

####Advanced: Google to figure out how to add a regression line to your plot!

####Homework is finishing this. If you want additional practice go to Learning_R_Additional_Practice.R

