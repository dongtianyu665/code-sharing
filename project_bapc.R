#################加载包
{library(Rgraphviz)
  library(Epi)
  library(BAPC)
  library(INLA)
  library(reshape2)
}
data(whostandard)
sum(whostandard$whoStandard[1:7])#58.38
sum(whostandard$whoStandard[8:9])#13.74
sum(whostandard$whoStandard[10:18])#27.91
sum(whostandard$whoStandard[19:21])#100.03
# 修改第7-9行的数值列
whostandard$agegroup <- factor(whostandard$agegroup, 
                               levels = c(levels(whostandard$agegroup), "<34", "34-", "45-"))
whostandard[19:21, "agegroup"] <- c("<34", "34-", "45-")
whostandard[19:21, "whoStandard"] <- c(7.76, 32.76, 59.47)
#################设置路径
setwd('D:/qs_paper')

#############################宏函数##################################
Prevalence_process <- function(location_expr, output_expr) {
  expr <- substitute({
    ################## 处理counts数据 ########################
    counts.uf <- read.csv("p_data/project.csv", fileEncoding = "UTF8")
    agegroups <- c("<34", "34-", "45-")
    x_list <- vector("list", length(agegroups))
    
    for (i in seq_along(agegroups)) {
      filtered_data <- subset(counts.uf,
                                tech ==LOCATION_PH &  # 占位符
                                agegroup == agegroups[i])
      x_list[[i]] <- matrix(filtered_data$percent)
    }
    
    names(x_list) <- paste0("x", 1:3)
    df <- as.data.frame(do.call(cbind, x_list))
    
    years <- 2016:2050
    counts_df <- data.frame(Year = years)
    for (i in 1:3) counts_df[[paste0("x", i)]] <- NA
    counts_df[1:9, 2:4] <- df
    
    write.table(counts_df, OUTPUT_FILE_PH, sep = "\t", row.names = F, col.names = F, quote = F)
    counts_data <- read.table(OUTPUT_FILE_PH, fileEncoding = "UTF8", row.names = 1, header = F)
    counts_data <- round(counts_data)
    
    write.table(counts_df, OUTPUT_FILE_PH, sep = "\t", row.names = F, col.names = F, quote = F)
    counts_data <- read.table(OUTPUT_FILE_PH, fileEncoding = "UTF8", row.names = 1, header = F)
    counts_data <- round(counts_data)
    ##################################### fit APC model #############################
    apc_obj <- APCList(counts_data, pop_100, gf = 5, agelab = agegroups)
    
    ######################## fit BAPC model ########################
    bapc_model <- BAPC(apc_obj, predict = list(npredict = 26, retro = FALSE),
                       model = list(
                         age = list(model = "rw2", prior = "loggamma", param = c(1, 5e-5)),
                         period = list(include = TRUE, model = "rw2", prior = "loggamma", param = c(1, 5e-5)),
                         cohort = list(include = TRUE, model = "rw2", prior = "loggamma", param = c(1, 5e-5)),
                         overdis = list(include = TRUE, model = "iid", prior = "loggamma", param = c(1, 0.005))
                       ),
                       secondDiff = F, stdweight = whostandard[19:21, 2], verbose = F)
    
    ############################## 标准化率表格 ###########################
    agestd_rate_df <- as.data.frame(bapc_model@agestd.rate)
    agestd_rate_df$mean_per_100 <- round(agestd_rate_df$mean * 1e2, 2)
    agestd_rate_df$sd_per_100 <- round(agestd_rate_df$sd * 1e2, 2)
    
    ########################################## 绘图 ################################
    jpeg(file = PLOT_FILE_PH, width = 2000, height = 1800, units = "px", res = 400)
    plotBAPC(bapc_model, scale = 1e2, type = "ageStdRate", obs.lwd = 0, obs.cex = 0.8,
             probs = seq(0.05, 0.95, by = 0.01))
    dev.off()
  })
  
  # 替换占位符为实际值
  expr <- do.call(substitute, list(expr, list(
    LOCATION_PH = location_expr,
    OUTPUT_FILE_PH = paste0("p_data/counts_p_", output_expr, ".txt"),
    PLOT_FILE_PH = paste0("p_data/counts_p_", output_expr, ".jpg")
  )))
  
  eval(expr, envir = parent.frame())  # 在调用环境中执行
}
########### 调用宏处理不同地区 ###########
pop_100 <- read.table("p_data/pop_100.txt", fileEncoding = "UTF8", row.names = 1, header = F)
Prevalence_process("AM", "AM")
Prevalence_process("HM", "HM")
Prevalence_process("LM", "LM")
Prevalence_process("VM", "VM")
Prevalence_process("LH", "LH")
Prevalence_process("OH", "OH")
Prevalence_process("FUAS", "FUAS")
Prevalence_process("UAE", "UAE")
Prevalence_process("RFA", "RFA")



Prevalence_process <- function(location_expr, output_expr) {
  expr <- substitute({
    ################## 处理counts数据 ########################
    counts.uf <- read.csv("p_data/project2.csv", fileEncoding = "UTF8")
    agegroups <- c("<34", "34-", "45-")
    x_list <- vector("list", length(agegroups))
    
    for (i in seq_along(agegroups)) {
      filtered_data <- subset(counts.uf,
                              tech ==LOCATION_PH &  # 占位符
                                agegroup == agegroups[i])
      x_list[[i]] <- matrix(filtered_data$pop_freq1)
    }
    
    names(x_list) <- paste0("x", 1:3)
    df <- as.data.frame(do.call(cbind, x_list))
    
    years <- 2016:2030
    counts_df <- data.frame(Year = years)
    for (i in 1:3) counts_df[[paste0("x", i)]] <- NA
    counts_df[1:9, 2:4] <- df
    
    write.table(counts_df, OUTPUT_FILE_PH, sep = "\t", row.names = F, col.names = F, quote = F)
    counts_data <- read.table(OUTPUT_FILE_PH, fileEncoding = "UTF8", row.names = 1, header = F)
    counts_data <- round(counts_data)
    
    write.table(counts_df, OUTPUT_FILE_PH, sep = "\t", row.names = F, col.names = F, quote = F)
    counts_data <- read.table(OUTPUT_FILE_PH, fileEncoding = "UTF8", row.names = 1, header = F)
    counts_data <- round(counts_data)
    ##################################### fit APC model #############################
    apc_obj <- APCList(counts_data, pop_100000, gf = 1, agelab = agegroups)
    
    ######################## fit BAPC model ########################
    bapc_model <- BAPC(apc_obj, predict = list(npredict = 6.5, retro = FALSE),
                       model = list(
                         age = list(model = "rw2", prior = "loggamma", param = c(1, 5e-5)),
                         period = list(include = TRUE, model = "rw2", prior = "loggamma", param = c(1, 5e-5)),
                         cohort = list(include = TRUE, model = "rw2", prior = "loggamma", param = c(1, 5e-5)),
                         overdis = list(include = TRUE, model = "iid", prior = "loggamma", param = c(1, 0.005))
                       ),
                       secondDiff = F, stdweight = whostandard[19:21, 2], verbose = F)
    
    ############################## 标准化率表格 ###########################
    agestd_rate_df <- as.data.frame(bapc_model@agestd.rate)
    agestd_rate_df$mean_per_100 <- round(agestd_rate_df$mean * 1e5, 2)
    agestd_rate_df$sd_per_100 <- round(agestd_rate_df$sd * 1e5, 2)
    
    ########################################## 绘图 ################################
    jpeg(file = PLOT_FILE_PH, width = 2000, height = 1800, units = "px", res = 400)
    plotBAPC(bapc_model, scale = 1e5, type = "ageStdRate", obs.lwd = 0, obs.cex = 0.8,
             probs = seq(0.05, 0.95, by = 0.01)
#注释             
         , ylim = c(0, 400),
           yaxt = "n")
    y_ticks <- c(0,100,200,300,400) # 刻度位置
    y_labels <- c("0","100","200","300",'400') # 刻度标签
    axis(2, at = y_ticks, labels = y_labels)
#注释 
    dev.off()
  })
  
  # 替换占位符为实际值
  expr <- do.call(substitute, list(expr, list(
    LOCATION_PH = location_expr,
    OUTPUT_FILE_PH = paste0("p_data/counts_p_", output_expr, ".txt"),
    PLOT_FILE_PH = paste0("p_data/counts_p_", output_expr, ".jpg")
  )))
  
  eval(expr, envir = parent.frame())  # 在调用环境中执行
}
########### 调用宏处理不同地区 ###########
pop_100000 <- read.table("p_data/pop_100k_2030.txt", fileEncoding = "UTF8", row.names = 1, header = F)

Prevalence_process("HM", "HM")
Prevalence_process("LM", "LM")
Prevalence_process("FUAS", "FUAS")

Prevalence_process("VM", "VM")
Prevalence_process("RFA", "RFA")
Prevalence_process("AM", "AM")
#Prevalence_process("LH", "LH")
#Prevalence_process("OH", "OH")
#Prevalence_process("UAE", "UAE")

