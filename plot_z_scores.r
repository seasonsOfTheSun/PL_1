library(ggplot2)
library(tidyr)
library(tibble)
library(hrbrthemes)
library(dplyr)

# Volcano dataset
#volcano



data <- read.csv("data/intermediate/disease_factor_z_score/rwd_factors_jensen_disease_db")


# Heatmap 
# Heatmap 
data_df = (data %>% as_tibble() %>% pivot_longer(col=-c(1:1), names_to="disease",values_to="Z_Score"))


plot <- (data_df %>% ggplot(aes(x=disease, y=X, fill= Z_Score)) + 
          geom_tile())
          
          
          
#           +
#          theme_ipsum() +
#          theme(legend.position="none"))

# rowid_to_column(var="factor_index")

#ndwjfkkjdsnakln
    
#  %>% rowid_to_column(var="Disease")

#  data %>% 
#    # Data wrangling
#    as_tibble() %>%
#    rowid_to_column(var="X") %>%
#    gather(key="Y", value="Z", -1) %>%

#    # Change Y to numeric
#    mutate(Y=as.numeric(gsub("V","",Y))))
    
    

#  dldfl %>%

    # Viz
#    ggplot(aes(X, Y, fill= Z)) + 
#      geom_tile() +
#      theme_ipsum() +
#      theme(legend.position="none")

#  ggsave(
#    "test.png",
#    plot)
