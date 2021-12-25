library(tidyverse)


Cell_number<-table(test.integrated@meta.data$seurat_clusters)
names(Cell_number)<-paste(rep("C",length(names(Cell_number))),names(Cell_number),sep = "")

cir_bar_data<-data.frame(cluster=names(Cell_number),number=as.numeric(Cell_number),id=1:length(names(Cell_number)))
# 构建输入数据


# 绘制基础环状条形图

label_data <- cir_bar_data

# calculate the ANGLE of the labels
number_of_bar <- nrow(label_data)
angle <-  90 - 360 * (label_data$id-0.5) /number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)

# calculate the alignment of labels: right or left
# If I am on the left part of the plot, my labels have currently an angle < -90
label_data$hjust<-ifelse( angle < -90, 1, 0)

# flip angle BY to make them readable
label_data$angle<-ifelse(angle < -90, angle+180, angle)
#查看标签数据
head(label_data)


# 绘制带标签的环状条形图
p <- ggplot(label_data, aes(x=as.factor(id), y=number)) +       # Note that id is a factor. If x is numeric, there is some space between the first bar
  
  # This add the bars with a blue color
  geom_bar(stat="identity", fill=alpha("skyblue", 0.7)) +
  
  # Limits of the plot = very important. The negative value controls the size of the inner circle, the positive one is useful to add size over each bar
  ylim(-2500,8500) +
  
  # Custom the theme: no axis title and no cartesian grid
  theme_minimal() +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(rep(-1,4), "cm")      # Adjust the margin to make in sort labels are not truncated!
  ) +
  
  # This makes the coordinate polar instead of cartesian.
  coord_polar(start = 0) +
  
  # Add the labels, using the label_data dataframe that we have created before
  geom_text(data=label_data, aes(x=id, y=number+100, label=cluster, hjust=hjust), 
            color="black", fontface="bold",alpha=0.6, size=2.5, 
            angle= label_data$angle, inherit.aes = FALSE ) 


# 绘制分组颜色的环状条形图

p <- ggplot(label_data, aes(x=as.factor(id), y=number, fill=cluster)) +       # Note that id is a factor. If x is numeric, there is some space between the first bar
  geom_bar(stat="identity", alpha=0.5) +
  ylim(-2500,8500) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(rep(-1,4), "cm") 
  ) +
  coord_polar() + 
  geom_text(data=label_data, aes(x=id, y=number+100, label=cluster, hjust=hjust), color="black", fontface="bold",alpha=0.6, size=2.5, angle= label_data$angle, inherit.aes = FALSE ) 


# 绘制分组+坐标+间隔的环状条形图

# 设置添加分组间隔
empty_bar <- 1
to_add <- data.frame( matrix(NA, empty_bar*nlevels(as.factor(label_data$cluster)), ncol(label_data)) )
colnames(to_add) <- colnames(label_data)
to_add$cluster <- rep(levels(as.factor(label_data$cluster)), each=empty_bar)
data <- rbind(label_data, to_add)
data <- data %>% arrange(group)
data$id <- seq(1, nrow(data))
head(data)

# 设置添加label标签信息

label_data <- data
number_of_bar <- nrow(label_data)
angle <- 90 - 360 * (label_data$id-0.5) /number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
label_data$hjust <- ifelse( angle < -90, 1, 0)
label_data$angle <- ifelse(angle < -90, angle+180, angle)
head(label_data)

#画图
p <- ggplot(label_data, aes(x=as.factor(id), y=number, fill=cluster)) +       # Note that id is a factor. If x is numeric, there is some space between the first bar
  geom_bar(stat="identity", alpha=0.5) +
  ylim(-2500,8500) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(rep(-1,4), "cm") 
  ) +
  coord_polar() + 
  geom_text(data=label_data, aes(x=id, y=number+100, label=cluster, hjust=hjust), color="black", fontface="bold",alpha=0.6, size=2.5, angle= label_data$angle, inherit.aes = FALSE ) 

p+ggplot2::annotate("text", x = rep(max(data$id),5), y = c(0, 2000, 4000, 6000, 8000), label = c("0", "2000", "4000", "6000", "8000") , color="grey", size=3 , angle=0, fontface="bold", hjust=1) 

