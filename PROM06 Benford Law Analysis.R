# USING BENFORD’S LAW TO DETECT INTERESTING PATTERNS IN DATA RESEARCH PROJECT
# PROM06 Assignment 2 (University of Sunderland)
# Developed by Ik Chin Fong (BI02SY)
#*************** Initial Commands ***************
cat("\014") #Clean console screen
rm(list=ls()) #clean all variables on global environment
setwd("C:/Users/Swift 3/Desktop/Benford's 2 digit law")

#*************** Packages Required ***************
library(dplyr) #loading dplyr package for data manipulation
library(ggplot2) #plot graphs
library(stringr) # tidy up data
library(moments) # compute skewness
library(pracma) #used for calculate cube root
library(tidyverse) # for arranging tables

#library(patchwork) #to combine plots
#library(benford.analysis) #validates a dataset using Benford’s Law
#*************** Reading Data ***************
#var <- "Global YouTube Statistics.csv"
#var <- "BolsaFamilia_Municipalities_Jan_2018.csv"
#var <- "Malaysia Election 2022 GE15.csv"
#var <- "Malaysia Election 2018 GE14.csv"
#var <- "Malaysia Elections GE11 to GE15.csv"
#var <- "governors_county_candidate.csv"
var <- "US and Malaysia General Election.csv"

df <- read.csv(file=var, header=TRUE) #load csv file

# remove .csv from file name to save name of data
dname <- str_replace(var, ".csv", ".") 

#*************** Data Matrix ***************
df_detail <- function(dataset) {
  print(paste('Number of columns of the matrix: ', length(dataset))) #print number of columns
  summary(dataset)
} 

df_detail(df)
head(df) # Print the records

#*************** Prepare Data *************** 
numdata <- df %>% 
  select_if(is.numeric) %>% #select column with numeric
  abs()%>% #remove -ve signs
  print() 

char_del <- setdiff(colnames(df), colnames(numdata)) #Identify removed column due to string

#loop and print the deleted columns
for (x in char_del) {
  print(paste("omitted: ", x, "values are strings.")) 
}

# count the number of columns with numeric data
count <- length(numdata)
print(paste("No. of columns left after removing columns with strings:", count))

#create empty lists to store omitted consecutive column names
seq_del <- c() 

# identify & remove columns with consecutive values
i <- 1
while (i <= count) {
  vec <- as.vector(numdata[[i]]) # change all values to vector format
  vec <- na.omit(vec) # NA of each columns will be omitted
  result <- rle(diff(vec)) #using rle function to check the difference between each values
  consec_values <- any(result$lengths>=2 & result$values==1) #consecutive values
  diff_values <- max(vec) - min(vec) #range of values
  
  #add column to seq_del list when range or values are too tight
  if(consec_values == TRUE | diff_values <= 500) { 
    seq_del <- append(seq_del, colnames(numdata[i])) 
  }
  i <- i + 1
}

#remove columns with consecutive values
if(!is.null(seq_del)){
  numdata <- numdata[, -which(names(numdata) %in% seq_del)] 
} 

#print out deleted columns due to sequence
for (x in seq_del) {
  print(paste("omitted: ", x, "due to sequence."))
}

numdata[numdata<=10] <- NaN #change all values < 10 to NaN for final omission

count <- length(numdata) #count the number of columns left
print(paste("No. of columns left:", count)) 

col_names <- names(numdata) #column names of the data

data <- list() #create empty list to store new NA omitted values 

i <- 1 # loop through the column with numeric data
while (i <= count) {
  data[i] <- na.omit(numdata[i]) %>% 
    format(scientific = TRUE, digits = 5) %>% #get values in significant digits
    print()
  
  i <- i+1
}

dtrim <- list() #create empty list to store data with white space removed.

i <- 1 # loop through the column with numeric data
while (i <= count){
  dtrim[[i]] <- str_trim(data[[i]], "left") %>% #trim white space to prevent "" values
    as.numeric() %>%
    print() # check results
  i <- i+1
}

dtrim

options(scipen=999)
#extract first digit by taking only first values available
fdigit <- lapply(dtrim, substr, 1, 1) %>% 
  print()
#extract second digit by taking only second values available
sdigit <- lapply(dtrim, substr, 2, 2) %>% 
  print()
#extract first and second digit by taking only first and second values available
fsdigit <- lapply(dtrim, substr, 1, 2) %>% 
  print()

#*************** Plotting Graphs *************** 

create_hist <- function(dataset, title){
  i <- 1 #loop over the columns
  while ( i <= count) {
    value <- as.numeric(dataset[[i]]) #change data type to number
    col_data <- as.data.frame(value) #set the data as data frame
    
    t_row <- nrow(col_data)
    
    #plot labelling
    p_title <- paste("A Histogram of", title, "N =", t_row) #using database name as plot title
    cname <- col_names[i]
    p_col <- paste(cname) #using column names as x-labels
    
    first_quart <- as.numeric(quantile(value, c(0.25))) #calculating 1st quartile
    third_quart <- as.numeric(quantile(value, c(0.75))) #calculating 4th quartile 
    max_val <- as.numeric(max(value)) #calculating greatest value of the column
    min_val <- as.numeric(min(value)) # calculating lowest value of the column
    
    #Freedman-Diaconis rule to calculate bins and binwidth. 
    #ref:https://medium.datadriveninvestor.com/how-to-decide-on-the-number-of-bins-of-a-histogram-3c36dc5b1cd8
    b_width <- (2 * (third_quart - first_quart) ) / nthroot(t_row,3)
    b_count <- ceiling((max_val - min_val)/b_width)
    
    #create Histogram
    print(ggplot(col_data, aes(x=value)) + 
            geom_histogram(color="darkgreen", fill="lightgreen", bins=b_count ) + 
            ggtitle(p_title) +
            theme_minimal() +
            theme(
              plot.title=element_text(family='', face='bold', colour='darkblue', size=18)
            ) +
            xlab(p_col))
    
    i <- i+1
  }
}



create_graphs <- function(ndigit, title){
  i <- 1
  while ( i <= count) {
    #creating frequency table
    tdigit <- table(ndigit[i]) # change individual column into individual frequency table
    pdigit <- prop.table(tdigit, margin = NULL)*100 #calculate the probability of each frequency table
    pdigit<- round(pdigit, digit=2) # rounding off each probability into 1 decimal place.
    prob_table <- as.data.frame(pdigit) #data frame for plotting bar graph and histogram using ggplot
    
    
    #leading digits is saved for title lableling
    ldigit <- ""
    
    #expected Benford's values 
    xaxis <- list()
    yaxis <- list()

    expected <- function(input) {
      if (0 %in% input$Var1) {
        # creating reference bar graph for 2nd digits
        xaxis <<- as.numeric(as.character(prob_table$Var1)) #x-axis with 2nd leading digits 
        ldigit <<- "Second Digit"
        # y axis with benford's law formula 
        for (x in xaxis) {
          yaxis <<- append(yaxis, sum(log10(1+(1/(10*1:9+x))))) #calculation of expected distribution
        }
        yaxis <<- round(as.numeric(yaxis)*100, digit=3)
      } else {
        # creating reference bar graph for 1st & 12 digits
        xaxis <<- as.numeric(as.character(prob_table$Var1)) #x-axis with leading digits
        expected <- log10(1+1/xaxis)*100 # y axis with benford's law formula 
        yaxis <<- round(expected, digit=2)
        if(10 %in% xaxis){
          ldigit <<- "First-Two Digits"
        } else {
          ldigit <<- "First Digit"
        }
      }
    }
    expected(prob_table)
    
    #plot labelling
    p_title <- paste(ldigit, "of", title) #using database name as plot title
    cname <- col_names[i]
    p_col <- paste(cname) #using column names as x-labels
    
    print(ggplot(prob_table) +
            geom_bar( aes(x=Var1, y=Freq), stat="identity", colour="darkblue", fill="lightblue") +
            geom_crossbar( aes(x=Var1, y=yaxis, ymin=yaxis, ymax=yaxis), 
                           colour="darkblue", 
                           linewidth=0.5, 
                           width=0.4, 
                           alpha= 0.5) + 
            ggtitle(p_title) +
            theme_minimal()  +
            theme(
              plot.title=element_text(family='', face='bold', colour='darkblue', size=18),
              axis.text.x = element_text(angle = 90,
                                         vjust = 0.5,
                                         hjust = 0.5)
            ) +
            xlab(p_col) +
            ylab("Probability [%]"))
    i <- i+1
  }
}

create_hist(data, dname)

create_graphs(fdigit, dname)
create_graphs(sdigit, dname)
create_graphs(fsdigit, dname)

#*************** Statistical Test *************** 
stat_results <- list()
stat_calc <- function(ndigit) {
  i <- 1
  while ( i <= count) {
    #creating observed frequency table
    tdigit <- table(ndigit[i]) # change individual column into individual frequency table
    pdigit <- prop.table(tdigit, margin = NULL) #calculate the probability of each frequency table
    observed<- round(pdigit, digit=2) # rounding off each probability into 1 decimal place.
    obs_table <<- as.data.frame(pdigit) #data frame for plotting bar graph and histogram using ggplot
    
    value <- as.numeric(data[[i]]) 
    data_col <- as.data.frame(value)
    t_row <- nrow(data_col) #total number of rows
    
    cname <- col_names[i] #column names
    
    #expected Benford's values 
    xaxis <- list()
    yaxis <- list()
    
    if (0 %in% obs_table$Var1) {
      # creating expected values of 2nd digit
      xaxis <- as.numeric(as.character(obs_table$Var1)) #x-axis with 2nd leading digits 
      # y axis with benford's law formula 
      for (x in xaxis) {
        yaxis <- append(yaxis, sum(log10(1+(1/(10*1:9+x))))) #calculation of expected distribution
      }
      yaxis <- round(as.numeric(yaxis), digit=10)
      
    } else {
      # creating expected values of 1st or 12 digits
      xaxis <- as.numeric(as.character(obs_table$Var1)) #x-axis with leading digits
      expected <- log10(1+(1/xaxis)) # y axis with benford's law formula 
      yaxis <- round(expected, digit=10)
    }
    
    exp_table <- as.data.frame(yaxis) # create table for expected Benford
    
    #Create function to display results
    output <- function(test_name, result){
      if (0 %in% obs_table$Var1) {
        stat_results$test <<- append(stat_results$test, test_name)
        stat_results$digit <<- append(stat_results$digit, "Second")
        stat_results$col <<- append(stat_results$col, cname)
        stat_results$val <<- append(stat_results$val, result)
        
      } else if (10 %in% obs_table$Var1) {
        stat_results$test <<- append(stat_results$test, test_name)
        stat_results$digit <<- append(stat_results$digit, "First-Two")
        stat_results$col <<- append(stat_results$col, cname)
        stat_results$val <<- append(stat_results$val, result)
        
      } else {
        stat_results$test <<- append(stat_results$test, test_name)
        stat_results$digit <<- append(stat_results$digit, "First")
        stat_results$col <<- append(stat_results$col, cname)
        stat_results$val <<- append(stat_results$val, result)
        
      }
    }
    
    #CHI SQUARE TEST 
    chisq_list <- list() #create list to store individual chisq values
    
    #function to calculate Chi-square individual values
    chisq_test <- function(obs, exp) {
      j <- 1
      while (j <= nrow(obs)) {
        chi_sq = (obs$Freq[j] - exp[j])**2/exp[j]
        chisq_list <<- append(chisq_list, chi_sq)
        j <- j + 1
      }
    }
    
    chisq_test(obs_table, yaxis) #call chi-square function
    
    chisq <- sum(as.numeric(chisq_list))*t_row #completion of chisq calculations
    chisq <- round(chisq, digit=3) #rounding the values to 3 decimal place
    
    output("chi_square", chisq) # call output function to print results)
    
    #SSD TEST 
    ssd_list <- list() #create list to store individual SSD values
    
    #Create function to calculate individual SSD
    ssd_test <- function(obs, exp) {
      l <- 1
      while (l <= nrow(obs)) {
        ssd <- (obs$Freq[l] - exp[l])**2
        ssd_list <<- append(ssd_list, ssd)
        l <- l + 1
      }
    }
    
    ssd_test(obs_table, yaxis) #Calling function to calculate SSD
    
    ssd <- sum(as.numeric(ssd_list)) * 100 *100 #Final SSD calculation
    ssd <- round(ssd, digit=3) #rounding to 3 decimal place
    
    output("SSD", ssd)# call output function to print results
    
    #MAD TEST 
    mad_list <- list() #create list to store individual MAD values
    
    #create function to calculating MAD (Mean Absolute Deviation)
    mad_test <- function(obs, exp) {
      t_digits <- nrow(obs)
      k <- 1
      while (k <= t_digits) {
        mad <- abs(obs$Freq[k] - exp[k])/t_digits
        mad_list <<- append(mad_list, mad)
        k <- k + 1
      }
    }
    mad_test(obs_table, yaxis)
    
    mad_total <- sum(as.numeric(mad_list))
    mad_total <- round(mad_total, digit = 4)
    
    output("MAD", mad_total) # call output function to print results

    i <- i+1
  }
}
stat_calc(fdigit) 
stat_calc(sdigit) 
stat_calc(fsdigit)
stat_results  <- as.data.frame(stat_results)
stat_results

test_result <- list()

calc <- function(dataset) {
  i <- 1
  while(i <= count){
    value <- as.numeric(dataset[[i]]) 
    obs_data <- as.data.frame(value)
    cname <- col_names[i]
    
    
    t_row <- nrow(obs_data)
    P99 <- as.numeric(quantile(value, c(0.99))) #calculating 99st quartile
    P1 <- as.numeric(quantile(value, c(0.01))) #calculating 1st quartile 
    max_val <- as.numeric(max(value)) #calculating greatest value of the column
    min_val <- as.numeric(min(value)) # calculating lowest value of the column
    
    #Order of Magnitude (OOM)
    oom <- log10(max_val/min_val)
    oom <- round(oom, digit=3)
    test_result$testname <<- append(test_result$testname, "OOM")
    test_result$colname <<- append(test_result$colname, cname)
    test_result$N <<- append(test_result$N, t_row)
    test_result$result <<- append(test_result$result, oom)
    
    
    #Robust Order of Magnitude (ROM)
    rom <- log10(P99/P1)
    rom <- round(rom, digit=3)
    test_result$testname <<- append(test_result$testname, "ROM")
    test_result$colname <<- append(test_result$colname, cname)
    test_result$N <<- append(test_result$N, t_row)
    test_result$result <<- append(test_result$result, rom)
    
    #skewness of graph
    sk <- skewness(value)
    sk <- round(sk, digit = 3)
    test_result$testname <<- append(test_result$testname, "Skewness")
    test_result$colname <<- append(test_result$colname, cname)
    test_result$N <<- append(test_result$N, t_row)
    test_result$result <<- append(test_result$result, sk)
    
    
    
    i <- i+1
  }
}

#calculate ROM for all columns. 
calc(data) # criteria, ROM > 2.5 to fit Benford
test_result <- as.data.frame(test_result)
test_result

#*************** order of magnitude Grading ***************
oom_gr <- filter(test_result, testname == "OOM")
rom_gr <- filter(test_result, testname == "ROM")


test_grade <- function(testr, testname) {
  testr$Threshold<-cut(testr$result,breaks=c(-Inf,2.5,3,Inf),labels=c("<2.5 (Non-Benford)","<=3 (Marginal)",">3 (Benford)"))
  print(testr)
  laby <- paste(testname, "Values")
  ptitle <- paste("A Scatter Plot of General Elections (X-Axis) against", testname, "Value (Y-Axis)")
  ggplot(testr,aes(as.factor(colname),result,color=Threshold))+
    geom_point()+
    ggtitle(ptitle)  +
    theme_minimal()  +
    xlab("General Elections") +
    ylab(laby) 
}
test_grade(oom_gr, "OOM")
test_grade(rom_gr, "ROM")

#*************** Scatter plots ***************

scat_plot <- function(table, testname, digit) {
  ptitle <- paste("A Scatter Plot of General Elections (X-Axis) against",testname, "Test Value (Y-Axis) of", digit, "digit")
  laby <- paste(testname, "Test")
  ggplot(table,aes(as.factor(col),val,color=Threshold))+
    geom_point()+
    ggtitle(ptitle)  +
    theme_minimal()  +
    xlab("General Elections") +
    ylab(laby) 
}

#*Chisq Test First-digit 
chisq_f <- filter(stat_results,  digit == "First" & test == "chi_square")
chisq_f$Threshold<-cut(chisq_f$val,breaks=c(-Inf,15.5,Inf),labels=c("<=15.5 (Benford)",">15.5 (Non-Benford)"))
scat_plot(chisq_f, "Chi-Square", "First")
chisq_f

#*Chisq Test Second-digit 
chisq_s <- filter(stat_results,  digit == "Second" & test == "chi_square")
chisq_s$Threshold<-cut(chisq_s$val,breaks=c(-Inf,16.9,Inf),labels=c("<=16.9 (Benford)",">16.9 (Non-Benford)"))
scat_plot(chisq_s, "chi_square", "Second")
chisq_s

#*Chisq Test First-Two-digits 
chisq_fs <- filter(stat_results,  digit == "First-Two" & test == "chi_square")
chisq_fs$Threshold<-cut(chisq_fs$val,breaks=c(-Inf,112,Inf),labels=c("<=112 (Benford)",">112 (Non-Benford)"))
scat_plot(chisq_fs, "chi_square", "First-Two")
chisq_fs

#*MAD Test First-digits
mad_f <- filter(stat_results,  digit == "First" & test == "MAD")
mad_f$Threshold<-cut(mad_f$val,breaks=c(-Inf,0.004, 0.008, 0.012,Inf),labels=c("<=0.004 (Benford)","<0.008 (Acceptable)","<0.012 (Marginal)",">0.012 (Non-Benford)"))
scat_plot(mad_f, "MAD", "First")
mad_f

#*MAD Test Second-digits 
mad_s <- filter(stat_results,  digit == "Second" & test == "MAD")
mad_s$Threshold<-cut(mad_s$val,breaks=c(-Inf,0.008, 0.012, 0.016,Inf),labels=c("<=0.008 (Benford)","<0.012 (Acceptable)","<0.016 (Marginal)",">0.016 (Non-Benford)"))
scat_plot(mad_s, "MAD", "Second")
mad_s

#*MAD Test First-Two-digits
mad_fs <- filter(stat_results,  digit == "First-Two" & test == "MAD")
mad_fs$Threshold<-cut(mad_fs$val,breaks=c(-Inf,0.0006, 0.0012, 0.0018,Inf),labels=c("<=0.0006 (Benford)","<0.0012 (Acceptable)","<0.0018 (Marginal)",">0.0018 (Non-Benford)"))
scat_plot(mad_fs, "MAD", "First-Two")
mad_fs

#*SSD Test First-digits
ssd_f <- filter(stat_results,  digit == "First" & test == "SSD")
ssd_f$Threshold<-cut(ssd_f$val,breaks=c(-Inf,2, 25, 100,Inf),labels=c("<=2 (Benford)","<25 (Acceptable)","<100 (Marginal)",">100 (Non-Benford)"))
scat_plot(ssd_f, "SSD", "First")
ssd_f

#*SSD Test Second-digits 
ssd_s <- filter(stat_results,  digit == "Second" & test == "SSD")
ssd_s$Threshold<-cut(ssd_s$val,breaks=c(-Inf,2, 10, 50,Inf),labels=c("<=2 (Benford)","<10 (Acceptable)","<50 (Marginal)",">50 (Non-Benford)"))
scat_plot(ssd_s, "SSD", "Second")
ssd_s

#*SSD Test First-Two-digits 
ssd_fs <- filter(stat_results,  digit == "First-Two" & test == "SSD")
ssd_fs$Threshold<-cut(ssd_fs$val,breaks=c(-Inf,2, 10, 50,Inf),labels=c("<=2 (Benford)","<10 (Acceptable)","<50 (Marginal)",">50 (Non-Benford)"))
scat_plot(ssd_fs, "SSD", "First-Two")
ssd_fs

#*************** Compiling all data ***************

all_gr <- filter(test_result, testname == "Skewness")
names(all_gr)[names(all_gr) == 'result'] <- 'Skewness'
all_gr <- subset(all_gr, select = -1)
all_gr
 
all_gr$OOM <- oom_gr$result
all_gr$ROM <- rom_gr$result

com_data <- function (col, colname) {
  all_gr$val <<- col$val
  #all_gr$Threshold <<- col$Threshold
  
  cname <- paste(colname)
  #tname <- paste(colname, "_T")
  
  names(all_gr)[names(all_gr) == 'val'] <<- cname
  #names(all_gr)[names(all_gr) == 'Threshold'] <<- tname
}
com_data(chisq_f, "XSQ_F")
com_data(chisq_s, "XSQ_S")
com_data(chisq_fs, "XSQ_FS")
com_data(mad_f, "MAD_F")
com_data(mad_s, "MAD_S")
com_data(mad_fs, "MAD_FS")
com_data(ssd_f, "SSD_F")
com_data(ssd_s, "SSD_S")
com_data(ssd_fs, "SSD_FS")

all_gr

ggplot(all_gr,aes(ROM,OOM))+
  geom_point()+
  ggtitle("A scatter plot of ROM values (x-Axis) against Skewness (Y-Axis)")  +
  theme_minimal()  +
  xlab("ROM values") +
  ylab("Skewness")


corr <- subset(all_gr, select = -1)

library(corrplot)
library(RColorBrewer)
M <-cor(corr)
corrplot(M, type="upper", order="hclust",
         col=brewer.pal(n=8, name="RdYlBu"))

library(PerformanceAnalytics)

chart.Correlation(M, histogram = TRUE, method = "pearson")

panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...) {
  usr <- par("usr")
  on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  Cor <- abs(cor(x, y)) # Remove abs function if desired
  txt <- paste0(prefix, format(c(Cor, 0.123456789), digits = digits)[1])
  if(missing(cex.cor)) {
    cex.cor <- 0.4 / strwidth(txt)
  }
  text(0.5, 0.5, txt,
       cex = 1 + cex.cor * Cor) # Resize the text by level of correlation
}

# Plotting the correlation matrix
pairs(M,
      upper.panel = panel.cor,    # Correlation panel
      lower.panel = panel.smooth) # Smoothed regression lines

#*************** Statistical Test Double Check ***************
#first digits
i <- 1 
while (i <= count) {
  cname <- col_names[i]
  value <- as.numeric(data[[i]]) 
  obs_data <- as.data.frame(value)
  bfordResult <- benford(obs_data$value, 1, discrete=TRUE, sign="both")
  print(cname)
  print(bfordResult)
  i <- i + 1
}

#First-Two digits
i <- 1 
while (i <= count) {
  cname <- col_names[i]
  value <- as.numeric(data[[i]]) 
  obs_data <- as.data.frame(value)
  bfordResult <- benford(obs_data$value, 2, discrete=TRUE, sign="both")
  print(cname)
  print(bfordResult)
  i <- i + 1
}




#*************** Generating Random Data - superseded ***************
count <- 1
random <- list()
random[[1]] <- sample(10:10000000, 10000, replace=TRUE)


samplefd <- lapply(random, substr, 1, 1) %>% #extract second digit by taking only second values available
  print()
samplesd <- lapply(random, substr, 2, 2) %>% #extract second digit by taking only second values available
  print()
samplefsd <- lapply(random, substr, 1, 2) %>% #extract second digit by taking only second values available
  print()

create_hist(random, "Normal Distribution")

create_graphs(samplefd, "Normal Distribution")
create_graphs(samplesd, "Normal Distribution")
create_graphs(samplefsd, "Normal Distribution")
#*************** creating graph showing ranges - superseded *************** 

i <- 1
value <- as.numeric(dataset[[i]]) 
    col_data <- as.data.frame(value)
    print(col_data)while (i <= count) {
  print(ggplot(numdata[i], aes(x=col_names[i], y=numdata[[i]], fill=col_names[i])) +
          geom_boxplot(alpha=0.7) +
          stat_summary(fun.y=mean, geom="point", shape=23, size=5, color="slateblue", fill="slateblue", alpha=0.9) +
          theme_minimal() +
          theme(legend.position="none") +
          scale_fill_brewer(palette="Set1") +
          xlab(dname) +
          ylab(""))
  i <- i + 1
  
}


labels <- col_names
dmin <- list()
dmax <- list()
gap <- list()


i <- 1
while (i <= count) {
  v <- na.omit(numdata[i])
  cname <- colnames(numdata[i])
  
  cat(paste("\n", cname,
            "\nskewness: ", skewness(v),
            "\nsample size: ", nrow(v)))
  
  maxv <- ceiling(as.numeric(max(v)))
  minv <- ceiling(as.numeric(min(v)))
  
  dmax <- append(dmax, maxv)
  dmin <- append(dmin, minv)
  gap <- append(gap, maxv-minv)
  i <- i + 1
}

min_val <- as.double(dmin)
max_val <- as.double(dmax)
gap <- as.double(gap)
max <- max_val
span <- tibble (labels, min_val, max_val, gap, max)

# make into long format for easier plotting  
df_long=span %>% 
  pivot_longer(
    c(min_val,max_val)
  )

df_long

nudge_value=.6

p_main=
  df_long %>% 
  
  # the following 3 lines of code are the same
  ggplot(aes(x=value,y=labels)) +
  geom_line(aes(group=labels), color="#E7E7E7", linewidth=3.5) +
  geom_point(aes(color=name), size=3) +
  
  # but we want geom_text for the data callouts and the legend
  
  # data callout
  geom_text(aes(label=value, color=name),
            size=3.25,
            nudge_x=if_else(
              df_long$value==df_long$max, # if it's the larger value...
              nudge_value,   # move it to the right of the point
              -nudge_value), # otherwise, move it to the left of the point
            hjust=if_else(
              df_long$value==df_long$max, #if it's the larger value
              0, # left justify
              1),# otherwise, right justify      
  )+
  
  # legend
  geom_text(aes(label=name, color=name), 
            data=. %>% filter(gap==max(gap)),
            nudge_y =.5, 
            fontface="bold",
            size=3.25)+  
  
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.y = element_text(color="black"),
        axis.text.x = element_text(color="#989898"),
        axis.title = element_blank(),
        panel.grid = element_blank()
  ) +
  labs(x="%",y=NULL) +
  scale_color_manual(values=c("#436685", "#BF2F24")) +
  
  #extend the y-axis otherwise the legend is cut off
  coord_cartesian(ylim=c(1, 7.5)) +
  
  #display percentages with % appended
  scale_x_continuous(labels = scales::percent_format(scale = 1)) 

p_main

p_gap=
  df_long %>% 
  ggplot(aes(x=gap,y=labels)) +
  geom_text(aes(x=0, label=gap),
            fontface="bold",
            size=3.25) +
  
  geom_text(aes(x=0, y=7), # 7 because that's the # of y-axis values
            label="Diff",
            nudge_y =.5, # match the nudge value of the main plot legend    
            fontface="bold",
            size=3.25) +
  
  theme_void() +
  coord_cartesian(xlim = c(-.05, 0.05), 
                  ylim=c(1,7.5) # needs to match main plot
  )+
  theme(
    plot.margin = margin(l=0, r=0, b=0, t=0), #otherwise it adds too much space
    panel.background = element_rect(fill="#EFEFE3", color="#EFEFE3"),
    legend.position = "none"
  )+
  scale_color_manual(values=c("#436685", "#BF2F24"))

p_gap

p_whole=
  # syntax from `patchwork`
  p_main + p_gap + plot_layout(design=
                                 c(
                                   area(l=0,  r=45, t=0, b=1), # defines the main figure area
                                   area(l=46, r=52, t=0, b=1)  # defines the gap figure area
                                 )) 

p_whole













i <- 1
while ( i<= count) {
  tdigit <- table(fdigit[i])
  pdigit <- prop.table(tdigit)*100 
  pdigit<- round(pdigit, digit=1)
  df <- as.data.frame(pdigit)
  
  #op <- par(mfrow=c(1, 2)) # to put graphs side by side
  #histdata <- as.numeric(data[[i]])
  #h <- hist(histdata, 
  #main=paste("Histogram of", col_names[i]),
  #col = "#a1e9f0", border = "black")
  
  
  #par(op) # restore
  
  
  i <- i+1
}

#*************** Statistic tests - superseded ***************
stat_results <- list()
stat_calc <- function(ndigit) {
  i <- 1
  while ( i <= count) {
    #creating observed frequency table
    tdigit <- table(ndigit[i]) # change individual column into individual frequency table
    pdigit <- prop.table(tdigit, margin = NULL) #calculate the probability of each frequency table
    observed<- round(pdigit, digit=2) # rounding off each probability into 1 decimal place.
    obs_table <<- as.data.frame(pdigit) #data frame for plotting bar graph and histogram using ggplot
    
    value <- as.numeric(data[[i]]) 
    data_col <- as.data.frame(value)
    t_row <- nrow(data_col) #total number of rows
    
    cname <- col_names[i] #column names
    
    #expected Benford's values 
    xaxis <- list()
    yaxis <- list()
    
    if (0 %in% obs_table$Var1) {
      # creating expected values of 2nd digit
      xaxis <- as.numeric(as.character(obs_table$Var1)) #x-axis with 2nd leading digits 
      # y axis with benford's law formula 
      for (x in xaxis) {
        yaxis <- append(yaxis, sum(log10(1+(1/(10*1:9+x))))) #calculation of expected distribution
      }
      yaxis <- round(as.numeric(yaxis), digit=10)
      
    } else {
      # creating expected values of 1st or 12 digits
      xaxis <- as.numeric(as.character(obs_table$Var1)) #x-axis with leading digits
      expected <- log10(1+(1/xaxis)) # y axis with benford's law formula 
      yaxis <- round(expected, digit=10)
    }
    
    exp_table <- as.data.frame(yaxis) # create table for expected Benford
    
    #Create function to display results
    output <- function(test_name, result){
      if (0 %in% obs_table$Var1) {
        stat_results$test <<- append(stat_results$test, test_name)
        stat_results$digit <<- append(stat_results$digit, "Second")
        stat_results$col <<- append(stat_results$col, cname)
        stat_results$val <<- append(stat_results$val, result)
        print(paste(as.character(test_name), "of second digit of", cname,"is:", result))
      } else if (10 %in% obs_table$Var1) {
        stat_results$test <<- append(stat_results$test, test_name)
        stat_results$digit <<- append(stat_results$digit, "First-Two")
        stat_results$col <<- append(stat_results$col, cname)
        stat_results$val <<- append(stat_results$val, result)
        print(paste(as.character(test_name), "of first-two digits of", cname,"is:", result))
      } else {
        stat_results$test <<- append(stat_results$test, test_name)
        stat_results$digit <<- append(stat_results$digit, "First")
        stat_results$col <<- append(stat_results$col, cname)
        stat_results$val <<- append(stat_results$val, result)
        print(paste(as.character(test_name), "of first digit of", cname,"is:", result))
      }
    }
    
    #CHI SQUARE TEST 
    chisq_list <- list() #create list to store individual chisq values
    
    #function to calculate Chi-square individual values
    chisq_test <- function(obs, exp) {
      j <- 1
      while (j <= nrow(obs)) {
        chi_sq = (obs$Freq[j] - exp[j])**2/exp[j]
        chisq_list <<- append(chisq_list, chi_sq)
        j <- j + 1
      }
    }
    
    chisq_test(obs_table, yaxis) #call chi-square function
    
    chisq <- sum(as.numeric(chisq_list))*t_row #completion of chisq calculations
    chisq <- round(chisq, digit=3) #rounding the values to 3 decimal place
    
    output("chi_square", chisq) # call output function to print results)
    
    #SSD TEST 
    ssd_list <- list() #create list to store individual SSD values
    
    #Create function to calculate individual SSD
    ssd_test <- function(obs, exp) {
      l <- 1
      while (l <= nrow(obs)) {
        ssd <- (obs$Freq[l] - exp[l])**2
        ssd_list <<- append(ssd_list, ssd)
        l <- l + 1
      }
    }
    
    ssd_test(obs_table, yaxis) #Calling function to calculate SSD
    
    ssd <- sum(as.numeric(ssd_list)) * 100 *100 #Final SSD calculation
    ssd <- round(ssd, digit=3) #rounding to 3 decimal place
    
    output("SSD", ssd)# call output function to print results
    
    #MAD TEST 
    mad_list <- list() #create list to store individual MAD values
    
    #create function to calculating MAD (Mean Absolute Deviation)
    mad_test <- function(obs, exp) {
      t_digits <- nrow(obs)
      k <- 1
      while (k <= t_digits) {
        mad <- abs(obs$Freq[k] - exp[k])/t_digits
        mad_list <<- append(mad_list, mad)
        k <- k + 1
      }
    }
    mad_test(obs_table, yaxis)
    
    mad_total <- sum(as.numeric(mad_list))
    mad_total <- round(mad_total, digit = 4)
    
    output("MAD", mad_total) # call output function to print results
    
    i <- i+1
  }
}
stat_calc(fdigit) 
stat_calc(sdigit) 
stat_calc(fsdigit)
stat_results  <- as.data.frame(stat_results)
stat_results

test_result <- list()

calc <- function(dataset) {
  i <- 1
  while(i <= count){
    value <- as.numeric(dataset[[i]]) 
    obs_data <- as.data.frame(value)
    cname <- col_names[i]
    
    t_row <- nrow(obs_data)
    P99 <- as.numeric(quantile(value, c(0.99))) #calculating 1st quartile
    P1 <- as.numeric(quantile(value, c(0.01))) #calculating 4th quartile 
    max_val <- as.numeric(max(value)) #calculating greatest value of the column
    min_val <- as.numeric(min(value)) # calculating lowest value of the column
    
    #Order of Magnitude (OOM)
    oom <- log10(max_val/min_val)
    oom <- round(oom, digit=3)
    print(paste("Order of Magnitude (OOM) of", cname, ":", oom))
    test_result$testname <<- append(test_result$testname, "OOM")
    test_result$colname <<- append(test_result$colname, cname)
    test_result$result <<- append(test_result$result, oom)
    
    #Robust Order of Magnitude (ROM)
    rom <- log10(P99/P1)
    rom <- round(rom, digit=3)
    print(paste("Robust Order of Magnitude (ROM) of", cname, ":", rom))
    test_result$testname <<- append(test_result$testname, "ROM")
    test_result$colname <<- append(test_result$colname, cname)
    test_result$result <<- append(test_result$result, rom)
    
    #skewness of graph
    sk <- skewness(value)
    sk <- round(sk, digit = 3)
    print(paste("Skewness of", cname, ":", sk))
    test_result$testname <<- append(test_result$testname, "Skewness")
    test_result$colname <<- append(test_result$colname, cname)
    test_result$result <<- append(test_result$result, sk)
    
    i <- i+1
  }
}

#calculate ROM for all columns. 
calc(data) # criteria, ROM > 2.5 to fit Benford
as.data.frame(test_result)

#*************** grading of Statistical tests - superseded ***************

#Statistical grading
#########################First Digit - Chi - square#####################
chisq_f <- filter(stat_results,  digit == "First" & test == "chi_square")

analysis <- list()

for(x in chisq_f$val) {
  if (x >=15.5){
    analysis <- append(analysis, "Non-Benford")
  } else if (x < 15.5) {
    analysis <- append(analysis, "Benford")
  }
}

chisq_f$analysis <- analysis

##################Second Digits - Chi-Square#######################

chisq_s <- filter(stat_results,  digit == "Second" & test == "chi_square")

analysis <- list()

for(x in chisq_s$val) {
  if (x >=16.9){
    analysis <- append(analysis, "Non-Benford")
  } else if (x < 16.9) {
    analysis <- append(analysis, "Benford")
  }
}

chisq_s$analysis <- analysis

#########################First two Digit - Chi - square#####################
chisq_fs <- filter(stat_results,  digit == "First-Two" & test == "chi_square")

analysis <- list()

for(x in chisq_fs$val) {
  if (x >=112){
    analysis <- append(analysis, "Non-Benford")
  } else if (x < 112) {
    analysis <- append(analysis, "Benford")
  }
}

chisq_fs$analysis <- analysis

#########################First Digit - MAD#####################
mad_f <- filter(stat_results,  digit == "First" & test == "MAD")

analysis <- list()

for(x in mad_f$val) {
  if (x >0.012){
    analysis <- append(analysis, "Non-Benford")
  } else if (x > 0.008) {
    analysis <- append(analysis, "MArginal")
  } else if (x >= 0.004) {
    analysis <- append(analysis, "Acceptable")
  } else if (x < 0.004) {
    analysis <- append(analysis, "Benford")
  }
}

mad_f$analysis <- analysis

##################Second Digits - MAD#######################

mad_s <- filter(stat_results,  digit == "Second" & test == "MAD")

analysis <- list()

for(x in mad_s$val) {
  if (x >0.016){
    analysis <- append(analysis, "Non-Benford")
  } else if (x > 0.012) {
    analysis <- append(analysis, "MArginal")
  } else if (x >= 0.008) {
    analysis <- append(analysis, "Acceptable")
  } else if (x < 0.008) {
    analysis <- append(analysis, "Benford")
  }
}

mad_s$analysis <- analysis

#########################First Two Digit - MAD#####################
mad_fs <- filter(stat_results,  digit == "First-Two" & test == "MAD")

analysis <- list()

for(x in mad_fs$val) {
  if (x >0.0018){
    analysis <- append(analysis, "Non-Benford")
  } else if (x > 0.012) {
    analysis <- append(analysis, "MArginal")
  } else if (x >= 0.006) {
    analysis <- append(analysis, "Acceptable")
  } else if (x < 0.006) {
    analysis <- append(analysis, "Benford")
  }
}

mad_fs$analysis <- analysis



#########################First Digit - SSD#####################
ssd_f <- filter(stat_results,  digit == "First" & test == "SSD")

analysis <- list()

for(x in ssd_f$val) {
  if (x >100){
    analysis <- append(analysis, "Non-Benford")
  } else if (x > 25) {
    analysis <- append(analysis, "MArginal")
  } else if (x >= 2) {
    analysis <- append(analysis, "Acceptable")
  } else if (x < 2) {
    analysis <- append(analysis, "Benford")
  }
}

ssd_f$analysis <- analysis

##################Second Digits - SSD#######################

ssd_s <- filter(stat_results,  digit == "Second" & test == "SSD")

analysis <- list()

for(x in ssd_s$val) {
  if (x >50){
    analysis <- append(analysis, "Non-Benford")
  } else if (x > 10) {
    analysis <- append(analysis, "MArginal")
  } else if (x >= 2) {
    analysis <- append(analysis, "Acceptable")
  } else if (x < 2) {
    analysis <- append(analysis, "Benford")
  }
}

ssd_s$analysis <- analysis



#########################First Two Digit - SSD#####################
ssd_fs <- filter(stat_results,  digit == "First-Two" & test == "SSD")

analysis <- list()

for(x in ssd_fs$val) {
  if (x >50){
    analysis <- append(analysis, "Non-Benford")
  } else if (x > 10) {
    analysis <- append(analysis, "MArginal")
  } else if (x >= 2) {
    analysis <- append(analysis, "Acceptable")
  } else if (x < 2) {
    analysis <- append(analysis, "Benford")
  }
}

ssd_fs$analysis <- analysis

chisq_f
chisq_s
chisq_fs

mad_f
mad_s
mad_fs

ssd_s
ssd_s
ssd_fs
oom_gr <- filter(test_result, testname == "OOM")
analysis <- list()

for(x in oom_gr$result) {
  if (x >3){
    analysis <- append(analysis, "Benford")
  } else if (x >= 2.5) {
    analysis <- append(analysis, "MArginal")
  } else if (x < 2.5) {
    analysis <- append(analysis, "Non-Benford")
  }
}

oom_gr$analysis <- analysis
oom_gr

rom_gr <- filter(test_result, testname == "ROM")

analysis <- list()

for(x in rom_gr$result) {
  if (x >3){
    analysis <- append(analysis, "Benford")
  } else if (x >= 2.5) {
    analysis <- append(analysis, "MArginal")
  } else if (x < 2.5) {
    analysis <- append(analysis, "Non-Benford")
  }
}

rom_gr$analysis <- analysis
rom_gr







  



