### Transform raw count data into appropriate matrix for use in survival and dispersal analyses
### Created by Joseph Scavetta --- (thank you!!!) 
### used in R version 3.6.3

setwd("xxx")

# Read in CSV of count table
data <- read.csv('raw_survival_2019_2020_example.csv')

# Rename features
colnames(data) <- c("bird_id", "host_species", "mite_species")

# Empty data frame to populate survival data
surv_data <- data.frame()

# The column indices of each day in the count table (increases by 2)
day_indices <- seq(4, 38, 2) # this is 17 indices for the 17 days (38-4 = 17*2 = 34)

# Loop through each row in the count table
for(row_index in 1:nrow(data))
{
    row <- data[row_index, ]
    
    # Loop through each day (ignore day 0 [index == 1], end at day 12 [index == 13]) - original data
    # Loop through each day (ignore day 0 [index == 1], end at day 17 [index == 18]) - new 2019_2020 data
    for(index in 2:18)
    {
        # Number of deaths on this day = number alive on previous day - (number alive this day + number that went missing)
        deaths <- as.numeric(row[day_indices[index-1]]) - (as.numeric(row[day_indices[index]]) + as.numeric(row[day_indices[index]+1]))
        if(deaths > 0)
        {
            # Create a death event for each death on this day
            for(i in 1:deaths)
            {
                surv_data <- rbind(surv_data, data.frame(day=index-1, event=1, row[1:3]))
            }
        }
        
        missing <- as.numeric(row[day_indices[index-1]+1])
        if(missing > 0)
        {
            # Create a censor event for each missing mite on this day
            for(i in 1:missing)
            {
                surv_data <- rbind(surv_data, data.frame(day=index-1, event=0, row[1:3]))
            }
        }
    }
}

# Total event (mites) for each bird just as a sanity check. Should equal number of mites on day 0 for each bird id
data.frame(table(surv_data$bird_id))

# Write survival data to a csv
write.csv(surv_data, 'mite_survival_2019_2020_example.csv', row.names = FALSE)

