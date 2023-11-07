p <- p100
p <- p[p$year %in% c(2011:2019),] ## Legacy code. p100$year is already within this range
df <- data.frame(m100=p$m100,  year=p$year)
df
duplicated(df$m100)


df2 <- df %>% group_by(m100)
df2 %>% tally ## shows the number of rows in each group

df4 <- df %>% group_by(m100) %>% filter(year >= min(year)) ## when filtering on grouped tibble, the filter values are calculated relative to other rows in teh group
df4 %>% tally ## shows the number of rows in each group. Hasn't worked.


df2 %>% slice(1) ## this works in that it selects only one row per group
df2 %>% arrange(year) %>% slice(1) ## can arrange years within group and then pick the first or the last
df2 %>% arrange(desc(year)) %>% slice(1)

df2 %>% sample(year) %>% slice(1) ## fails

df2 %>% sample_n(size=1) ## works - randomly samples one row. Fumnction superceded by slice_sample

df2 %>% slice_sample(n=1, weight_by(sample(a$year, size=3))) ## Don't know how to set the size for each group in the sample function

##### This works! #####
## Try this with a function instead
subdf <- split(df, df$m100)[[1]]
nrow(subdf)
subdf
x <- sample(a$year, size=nrow(subdf), replace=T) ## create a distribution of absence years the length needed to sample the presence years
sample(subdf$year, size=1, prob=x)
sample(subdf$year, size=10, prob=x, replace=T) ## Demonstrate that sample can return different years each time

df2 <- lapply(split(df, df$m100),
              function(subdf){
                
                # test <- split(df, df$m100)
                # for(i in 1:620) {
                #   subdf <- test[[i]]
                if(nrow(subdf) >1) {   
                  x <- sample(a$year, size=nrow(subdf), replace=T) ## create a distribution of absence years the length needed to sample the presence years
                  y <- sample(subdf$year, size=1, prob=x)
                } else { y <- subdf$year}
                
                out <- cbind(subdf$m100[1], y)
                # subdf[sample(1:nrow(subdf), 1),]
              }
              
)

df2 <- do.call('rbind', df2)
colnames(df2) <- c("m100", "year")


