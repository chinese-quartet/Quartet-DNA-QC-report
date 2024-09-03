# #' Get query data and calculate score for generating DNA QC report
# #' @export
# #' @importFrom data.table fread
# #' @importFrom dplyr mutate 
# #' @importFrom dplyr across 
# #' @importFrom dplyr case_when 
# #' @importFrom dplyr group_by 
# #' @importFrom dplyr summarize 
# #' @importFrom dplyr ungroup 
# #' @importFrom dplyr arrange 
# #' @importFrom dplyr pull
# #' @importFrom dplyr select
# #' @importFrom dplyr filter 
# #' @importFrom dplyr left_join
# #' @importFrom stringr str_detect
# #' @importFrom tidyr pivot_longer 
# #' @importFrom tidyr unnest_wider
# #' @importFrom stats quantile


# Dnaseqc <- function(variant_qc_file, mendelian_qc_file, data_type){
  
#   ## 检查数据形式
  
#   ## 确认测序类型WGS或者WES
#   if (data_type %in% c("WGS","wgs")){
#     WGS_ref_metric <- system.file("extdata","WGS_historical_qc_metrics.csv",package = "dnaseqc")
#     ref_metric_dt <- fread(WGS_ref_metric,sep = '\t')
#   }
#   else {
#     WES_ref_metric <- system.file("extdata","WES_historical_qc_metrics.csv",package = "dnaseqc")
#     ref_metric_dt <- fread(WES_ref_metric,sep = '\t')
#   }
  
#   variant_qc_dt <- fread(variant_qc_file,sep='\t')
#   mendelian_qc_dt <- fread(mendelian_qc_file,sep='\t')
  
#   ### 处理variant_qc数据
#   variant_qc_dt <- variant_qc_dt %>%
#     mutate(across(c(`SNV precision`, `INDEL precision`, `SNV recall`, `INDEL recall`), 
#                   ~round(.x / 100, 4)))
#   ## 添加query数据名称
#   one_set <- c('Queried_Data_Set', 'Queried', 'Queried_Data')
  
#   ## 对相同名称的样本取平均值
#   metrics <- variant_qc_dt %>%
#     mutate(Sample_Group = case_when(
#       str_detect(Sample, "LCL5|D5") ~ "D5",
#       str_detect(Sample, "LCL6|D6") ~ "D6",
#       str_detect(Sample, "LCL7|F7") ~ "F7",
#       str_detect(Sample, "LCL8|M8") ~ "M8"
#     )) %>%
#     group_by(Sample_Group) %>%
#     summarize(
#       snv_precision = mean(`SNV precision`, na.rm = TRUE), # 计算平均值，忽略NA
#       snv_recall = mean(`SNV recall`, na.rm = TRUE), # 计算平均值，忽略NA
#       indel_precision = mean(`INDEL precision`, na.rm = TRUE), # 计算平均值，忽略NA
#       indel_recall = mean(`INDEL recall`, na.rm = TRUE) # 计算平均值，忽略NA
#     ) %>%
#     ungroup() %>%
#     pivot_longer(-Sample_Group, names_to = "Metric", values_to = "Value") %>%
#     arrange(Sample_Group, Metric)
  
#   # 添加metrics值到one_set
#   one_set <- c(one_set, metrics$Value)
#   # 定义列名
#   column_names <- c('sample', 'group', 'batch', 'snv_D5-precision', 'snv_D5-recall', 'snv_D6-precision', 'snv_D6-recall', 'snv_F7-precision', 'snv_F7-recall', 'snv_M8-precision', 'snv_M8-recall', 'indel_D5-precision', 'indel_D5-recall', 'indel_D6-precision', 'indel_D6-recall', 'indel_F7-precision', 'indel_F7-recall', 'indel_M8-precision', 'indel_M8-recall')
#   variant_metric_df <- data.frame(t(one_set))
#   names(variant_metric_df) <- column_names
  
#   # 初始化列表
#   one_set <- c("Queried_Data_Set", "Queried", "Queried_Data") 
  
#   ## 计算孟德尔遗传率平均值
#   snv_mendelian <- mendelian_qc_dt %>% 
#     filter(str_detect(Family, "SNV$")) %>% 
#     summarize(Mean_Mendelian_Concordance_Rate = mean(Mendelian_Concordance_Rate, na.rm = TRUE)) %>%
#     pull(Mean_Mendelian_Concordance_Rate)
  
#   indel_mendelian <- mendelian_qc_dt %>% 
#     filter(str_detect(Family, "INDEL$")) %>% 
#     summarize(Mean_Mendelian_Concordance_Rate = mean(Mendelian_Concordance_Rate, na.rm = TRUE)) %>%
#     pull(Mean_Mendelian_Concordance_Rate)
  
#   ## 添加到one_set
#   one_set <- c(one_set, snv_mendelian, indel_mendelian)
  
#   ## 转换成DataFrame
#   mendelian_df <- data.frame(matrix(one_set, ncol=5, byrow=TRUE),
#                              stringsAsFactors=FALSE)
#   # 设置列名
#   names(mendelian_df) <- c('sample', 'group', 'batch', 'snv_mendelian', 'indel_mendelian')
  
#   ##合并variant 和 mendelian数据
#   queried_performance <- variant_metric_df %>%
#     left_join(mendelian_df, by = c("sample", "group", "batch"))
  
#   ## 计算qc 指标
#   quality_metrics_list <- list()
#   batches <- unique(queried_performance$batch)
  
#   for (bat in batches) {
#     # 使用dplyr过滤特定批次的数据
#     df_batch <- queried_performance %>% filter(batch == bat)
    
#     ## 转换数据类型
#     df_batch <- df_batch %>%
#       mutate(across(-c(batch,group,sample), as.numeric))
#     # 计算指标
#     precision_snv <- mean(df_batch %>% dplyr::select(matches('snv_.*precision')) %>% unlist(), na.rm = TRUE)
#     precision_indel <- mean(df_batch %>% dplyr::select(matches('indel_.*precision')) %>% unlist(), na.rm = TRUE)
#     recall_snv <- mean(df_batch %>% dplyr::select(matches('snv_.*recall')) %>% unlist(), na.rm = TRUE)
#     recall_indel <- mean(df_batch %>% dplyr::select(matches('indel_.*recall')) %>% unlist(), na.rm = TRUE)
#     mendelian_snv <- mean(df_batch %>% dplyr::select(contains('snv_mendelian')) %>% unlist(), na.rm = TRUE)
#     mendelian_indel <- mean(df_batch %>% dplyr::select(contains('indel_mendelian')) %>% unlist(), na.rm = TRUE)
    
#     # 计算总分
#     beta2 <- 0.25
#     precision_beta <- (1 + beta2) * precision_snv * precision_indel / (beta2 * precision_snv + precision_indel)
#     recall_beta <- (1 + beta2) * recall_snv * recall_indel / (beta2 * recall_snv + recall_indel)
#     mendelian_beta <- (1 + beta2) * mendelian_snv * mendelian_indel / (beta2 * mendelian_snv + mendelian_indel)
#     total <- (precision_beta + recall_beta + mendelian_beta) / 3
    
#     # 添加到列表
#     quality_metrics_list[[length(quality_metrics_list) + 1]] <- list(bat, round(precision_snv, 5), round(precision_indel, 5), round(recall_snv, 5), round(recall_indel, 5), round(mendelian_snv, 5), round(mendelian_indel, 5), round(total, 5))
#   }
  
#   quality_metrics_df <- do.call(rbind, lapply(quality_metrics_list, function(x) as.data.frame(t(x), stringsAsFactors = FALSE)))
#   names(quality_metrics_df) <- c('batch', 'precision_snv', 'precision_indel', 'recall_snv', 'recall_indel', 'mendelian_snv', 'mendelian_indel', 'total')
  
#   ## 添加表现列
#   quality_metrics_df$precision_snv_performance <- NA
#   quality_metrics_df$precision_indel_performance <- NA
#   quality_metrics_df$recall_snv_performance <- NA
#   quality_metrics_df$recall_indel_performance <- NA
#   quality_metrics_df$mendelian_snv_performance <- NA
#   quality_metrics_df$mendelian_indel_performance <- NA
#   quality_metrics_df$total_performance <- NA
  
  
#   quantile_results <- data.frame()
#   ## 从历史数据中获取百分位数值，并进行评分
#   columns_to_process <- names(ref_metric_dt)[!names(ref_metric_dt) %in% c("batch", grep("*_performance", names(ref_metric_dt), value = TRUE))]
  
#   for (col_name in columns_to_process) {
#     # 对每个列展平列表并转换为数值型
#     column_data <- unlist(ref_metric_dt[[col_name]])
#     column_data_numeric <- as.numeric(column_data)
    
#     # 计算分位数
#     Q1 <- quantile(column_data_numeric, probs = 0.2, na.rm = TRUE)
#     Q2 <- quantile(column_data_numeric, probs = 0.5, na.rm = TRUE)
#     Q3 <- quantile(column_data_numeric, probs = 0.8, na.rm = TRUE)
    
#     # 构造临时数据框以存储结果
#     temp_df <- data.frame(column_name = col_name,
#                           Q1 = Q1, Q2 = Q2, Q3 = Q3)
    
#     # 合并结果到最终数据框
#     quantile_results <- rbind(quantile_results, temp_df)
#   }
  
#   # 转换结果数据框以便于分析
#   quantile_df <- pivot_longer(quantile_results, cols = c("Q1", "Q2", "Q3"), names_to = "Quantile", values_to = "Value")
  
#   # 重新整理数据框，确保每个分位数的值与对应的列名和分位数级别相匹配
#   quantile_df <- quantile_df %>%
#     mutate(Quantile = factor(Quantile, levels = c("Q1", "Q2", "Q3"))) %>%
#     arrange(column_name, Quantile)
  
#   ## 定义评级函数
#   assign_performance <- function(metric_value, q1, q2, q3) {
#     if (is.na(metric_value)) return(NA)
#     if (metric_value < q1) {
#       return("Bad")
#     } else if (metric_value < q2) {
#       return("Fair")
#     } else if (metric_value < q3) {
#       return("Good")
#     } else {
#       return("Great")
#     }
#   }
  
#   # 为每个指标应用性能等级
#   metrics <- c("precision_snv", "precision_indel", "recall_snv", "recall_indel", "mendelian_snv", "mendelian_indel", "total")
#   for (metric in metrics) {
#     q1 <- quantile_df %>% filter(column_name == metric, Quantile == "Q1") %>% pull(Value)
#     q2 <- quantile_df %>% filter(column_name == metric, Quantile == "Q2") %>% pull(Value)
#     q3 <- quantile_df %>% filter(column_name == metric, Quantile == "Q3") %>% pull(Value)
    
#     quality_metrics_df[[paste0(metric, "_performance")]] <- mapply(assign_performance, quality_metrics_df[[metric]], q1, q2, q3)
#   }
  
#   ## 整合历史数据与Query数据
#   Merge_data <- rbind(ref_metric_dt,quality_metrics_df)
  
#   ## 获取数据排名
#   full_name <- list(
#     precision_snv = 'Precision (SNV)', 
#     precision_indel = 'Precision (INDEL)', 
#     recall_snv = 'Recall (SNV)', 
#     recall_indel = 'Recall (INDEL)',
#     mendelian_snv = 'Mendelian Concordance Rate (SNV)', 
#     mendelian_indel = 'Mendelian Concordance Rate (INDEL)',
#     total = 'Total Score'
#   )
  
#   # 初始化一个空的数据框用于存储排名结果
#   rank_df <- data.frame(batch = unlist(Merge_data$batch))
#   # 需要计算排名的列名
#   metrics_columns <- names(full_name)
#   # 对每个需要计算排名的列进行循环
#   for (col_name in metrics_columns) {
#     # 展平列表并转换为数值向量
#     numeric_vector <- sapply(Merge_data[[col_name]], function(x) as.numeric(unlist(x)))
#     # 计算排名
#     rank_column <- rank(-numeric_vector, ties.method = "min", na.last = "keep")
#     # 将排名结果添加到rank_df中
#     rank_df[[col_name]] <- rank_column
#   }
  
  
#   ## scale total score
#   score_raw_ref = Merge_data[Merge_data$batch != 'Queried_Data',"total"]
#   score_raw_ref <- score_raw_ref %>% unlist()
#   score_raw_ref <- as.numeric(score_raw_ref)
#   score_raw = Merge_data$total
#   score_raw <- as.numeric(score_raw)
  
#   k <- 9/(max((score_raw_ref)-min(score_raw_ref)))
#   score_norm <- round(1 + k * (score_raw - min(score_raw_ref)), 2)
  
#   # Q1_total <- quantile_df %>% filter(column_name == "total", Quantile == "Q1") %>% pull(Value)
#   # Q2_total <- quantile_df %>% filter(column_name == "total", Quantile == "Q2") %>% pull(Value)
#   # Q3_total <- quantile_df %>% filter(column_name == "total", Quantile == "Q3") %>% pull(Value)
  
#   Merge_data$total_norm <- score_norm
  
#   # Q1_v = round(1 + k * (Q1_total - min(score_raw_ref)), 2)
#   # Q2_v = round(1 + k * (Q2_total - min(score_raw_ref)), 2)
#   # Q3_v = round(1 + k * (Q3_total - min(score_raw_ref)), 2)
  
#   if (Merge_data[Merge_data$batch == "Queried_Data","total_norm"] <= 1){
#     Merge_data[Merge_data$batch == "Queried_Data","total_norm"] <- 1
#   } else if (Merge_data[Merge_data$batch == "Queried_Data","total_norm"] >= 10){
#     Merge_data[Merge_data$batch == "Queried_Data","total_norm"] <- 10
#   }
  
#   Merge_data
#   ## 获取评估结果
#   # 初始化结果列表
  
#   evaluation_metrics <- list()
#   for (metric in names(full_name)) {
#     queried_value <- Merge_data[Merge_data$batch == "Queried_Data", ][[metric]]
#     mean_std <- sprintf('%.3f ± %.3f', mean(ref_metric_dt[[metric]], na.rm = TRUE), sd(ref_metric_dt[[metric]], na.rm = TRUE))
#     rank_info <- sprintf('%.0f / %.0f', rank_df[rank_df$batch == "Queried_Data", ][[metric]], nrow(rank_df))
#     performance <- Merge_data[Merge_data$batch == "Queried_Data", ][[paste0(metric, "_performance")]]
    
#     # 添加到列表
#     evaluation_metrics[[length(evaluation_metrics) + 1]] <- list(
#       Name = full_name[[metric]],
#       Queried_Value = ifelse(length(queried_value) > 0, queried_value[1], NA),
#       Mean_Std = mean_std,
#       Rank_Info = rank_info,
#       Performance = ifelse(length(performance) > 0, performance[1], NA)
#     )
#   }
  
#   evaluation_metrics_df <- do.call(rbind, lapply(evaluation_metrics, function(x) as.data.frame(t(unlist(x)), stringsAsFactors = FALSE)))
#   names(evaluation_metrics_df) <- c("Metric", "Queried Value", "Mean ± SD", "Rank", "Performance")
  
#   evaluation_metrics_df$`Queried Value` <- round(as.numeric(evaluation_metrics_df$`Queried Value`),3)
  
#   Merge_data <- Merge_data %>%
#     mutate(across(everything(), ~unlist(.)))
  
#   ### 计算F1 score(
#   Merge_data$snv_f1 = 2*(unlist(Merge_data$precision_snv)*unlist(Merge_data$recall_snv))/(unlist(Merge_data$precision_snv)+unlist(Merge_data$recall_snv))
#   Merge_data$indel_f1 = 2*(Merge_data$precision_indel*Merge_data$recall_indel)/(Merge_data$precision_indel+Merge_data$recall_indel)
  

  
  
#   ## 添加绘图标签
#   Merge_data$Type <- "Reference"
#   Merge_data[Merge_data$batch == "Queried_Data",]$Type <- "Query"
  
#   pt_snv_mendelian_f1 <- plot_scatter_box(Merge_data, var_x = 'mendelian_snv', var_y = 'snv_f1', 
#                                           col_g = 'Type', xlab = 'Mendelian Concordence Rate', ylab = 'F1 score', 
#                                           title_lab = 'SNV Performance')
  
#   pt_indel_mendelian_f1 <- plot_scatter_box(Merge_data, var_x = 'mendelian_indel', var_y = 'indel_f1', 
#                                             col_g = 'Type', xlab = 'Mendelian Concordence Rate', ylab = 'F1 score', 
#                                             title_lab = 'INDEL Performance')
  
#   return(list(conclusion=evaluation_metrics_df,
#               rank_table=Merge_data,
#               vcf_table=variant_qc_dt,
#               mendelian_table=mendelian_df,
#               p_mendelian_f1_snv=pt_snv_mendelian_f1,
#               p_mendelian_f1_indel=pt_indel_mendelian_f1))
  
  
# }
