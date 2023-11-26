# change to development_data

development_data[development_data$colony=="SD18-11"&development_data$sang_prop>0.1,"remarks"] <- "queenless colony"
