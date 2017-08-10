library(ggplot2)
library(dplyr)
library(isomiRs)
data(mirData)
ids<-mirData
column <- "group"
rawList <- metadata(ids)[[1]]
metadata <- as.data.frame(colData(ids))
metadata$sample <- as.character(row.names(metadata))
metadata$condition <- metadata[,column]
.l <- lapply(row.names(colData(ids)), function(sample){
    .d <- rawList[[sample]] %>%
        mutate(mismtag=ifelse(mism != "0", "Yes", "No")) %>%
        mutate(addtag=ifelse(add != "0", "Yes", "No")) %>%
        mutate(t5tag=ifelse(t5 != "0", "Yes", "No")) %>%
        mutate(t3tag=ifelse(t3 != "0", "Yes", "No"))
    .s <- data.frame(sample=sample,
                     type=c("mm","add","t5","t3","ref"),
                     method="freq",
                     freq=c(sum(.d$freq[.d$mismtag=="Yes"], na.rm=TRUE),
                            sum(.d$freq[.d$addtag=="Yes"], na.rm=TRUE),
                            sum(.d$freq[.d$t5tag=="Yes"], na.rm=TRUE),
                            sum(.d$freq[.d$t3tag=="Yes"], na.rm=TRUE),
                            sum(.d$freq[.d$mismtag!="Yes" & .d$addtag!="Yes"  & .d$t5tag!="Yes" & .d$t3tag!="Yes"])),
                     stringsAsFactors = FALSE
    )
    .s$freq <- .s$freq/sum(.s$freq) * 100
    .u <- data.frame(sample=sample,
                     type=c("mm","add","t5","t3","ref"),
                     method="unique",
                     freq=c(sum(.d$mismtag=="Yes"),
                            sum(.d$addtag=="Yes"),
                            sum(.d$t5tag=="Yes"),
                            sum(.d$t3tag=="Yes"),
                            sum(.d$mismtag!="Yes" & .d$addtag!="Yes"  & .d$t5tag!="Yes" & .d$t3tag!="Yes")),
                     stringsAsFactors = FALSE
    )
    .u$freq <- .u$freq/sum(.u$freq) * 100
    rbind(.s, .u)
})
p<-do.call(rbind, .l) %>% arrange(type) %>%
    left_join(metadata, by="sample") %>%
    ggplot(aes_string(x="type", y="freq", group="sample", color="condition")) +
    geom_polygon(fill=NA, size=1.5) +
    coord_polar(start=-pi) +
    scale_color_manual(values = c("steelblue", "steelblue")) +
    theme_bw() + ylab("") + xlab("") + guides(color=FALSE) +
    theme(panel.grid.major = element_line(size = 1)) +
    theme(axis.title.x=element_blank(),
          axis.ticks=element_blank(),
          panel.border = element_blank(),
          axis.text.y = element_blank(),
          plot.background = element_rect(fill = "transparent",colour = NA),
          panel.background = element_rect(fill = "transparent",colour = NA),
          axis.text.x = element_blank()) +
    scale_y_continuous(expand = c(0,0))

hexSticker::sticker(p,
        package="isomiRs", 
        s_x=0.93,
        s_y=0.8, 
        s_width=1.5, 
        s_height=1.5,
        p_x = 1,
        p_y = 1.55,
        h_color = "orange2",
        h_fill = rgb(0, 0, 0, alpha = 0.7),
        h_size = 1,
        p_color = "orange2",
        p_size = 24,
        filename="inst/stickers/isomirs.png")
