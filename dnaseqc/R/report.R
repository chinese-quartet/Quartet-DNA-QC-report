# --------------------------------------------------------------------------- #
#' Generate Quartet Genomics report
#'
#' @description
#' Use calculated DNA result to generate QC report
#'
#' @param qc_result list
#' @param report_template character
#' @param report_dir character
#' @param report_name character
#'
#' @importFrom dplyr %>%
#' @importFrom flextable flextable
#' @importFrom flextable theme_vanilla
#' @importFrom flextable theme_box
#' @importFrom flextable color
#' @importFrom flextable set_caption
#' @importFrom flextable align
#' @importFrom flextable width
#' @importFrom flextable bold
#' @importFrom flextable nrow_part
#' @importFrom flextable bg
#' @importFrom flextable fontsize
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 theme_minimal
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 element_blank
#' @importFrom ggplot2 element_text
#' @importFrom ggplot2 geom_rect
#' @importFrom ggplot2 scale_x_continuous
#' @importFrom ggplot2 theme
#' @importFrom officer body_add_par
#' @importFrom flextable body_add_flextable
#' @importFrom officer body_add_gg
#' @importFrom officer body_add_break
#' @importFrom officer read_docx
#' @importFrom officer fp_text
#' @importFrom officer fpar
#' @importFrom officer ftext
#' @importFrom officer body_add_fpar
#' @importFrom stats quantile
#'
#' @export
generate_dna_report <- function(qc_result,
                                report_template,
                                report_dir = NULL,
                                report_name = NULL) {
  if (is.null(qc_result) || is.null(report_template)) {
    stop("All arguments (qc_result, report_template) are required.")
  }

  if (is.null(report_dir)) {
    path <- getwd()
    sub_dir <- "output"
    dir.create(file.path(path, sub_dir), showWarnings = FALSE)
    report_dir <- file.path(path, "output")
  }
  if (is.null(report_name)) {
    report_name <- "Quartet_DNA_Report.docx"
  }
  output_file <- file.path(report_dir, report_name)

  # --- 1. 定义中文文本内容 ---
  
  # 摘要
  text_sum_intro <- "本报告基于多项组学关键质量控制指标，总结了 Quartet DNA 参考物质所生成数据的质量情况。质量控制流程从变异检测文件（variant calling file, VCF）开始，分别计算了SNV和INDEL的精确度(Precision）、灵敏度(Recall) 及整体质量判别。"
  
  # 定义
  text_def_title <- "质量控制指标"
  
  text_prec_title <- "精确度 (Precision)"
  text_prec_desc <- "采用 hap.py工具 (https://github.com/Illumina/hap.py) 将测试数据集中的变异与基准变异集进行比较。精确度定义为测试数据集中被判定为真实变异的比例。"
  
  text_rec_title <- "灵敏度 (Recall)"
  text_rec_desc <- "灵敏度定义为所有真实变异中，在测试数据集中被成功检测到的比例。"
  
  # 参考文献
  text_ref_title <- "参考文献"
  text_ref_1 <- "1. Zheng Y, Liu Y, Yang J, et al. Multi-omics data integration using ratio-based quantitative profiling with Quartet reference materials. Nature Biotechnology, 2024."
  text_ref_2 <- "2. Ren L, Duan X, Dong L, et al. Quartet DNA reference materials and datasets for comprehensively evaluating germline variant calling performance. Genome Biology, 2023."
  text_ref_3 <- "3. GB/T 45214-2025《人全基因组高通量测序数据质量评价方法》"
  text_ref_4 <- "4. 上海临床队列组学检测工作指引（征求意见稿）, 2025/11/26"
  
  # 免责声明
  text_disc_title <- "免责声明"
  text_disc_content <- "本数据质量报告仅针对所评估的特定数据集提供分析结果，仅供信息参考之用。尽管已尽最大努力确保分析结果的准确性和可靠性，但本报告按“现状（AS IS）”提供，不附带任何形式的明示或暗示担保。报告作者及发布方不对基于本报告内容所采取的任何行动承担责任。本报告中的结论不应被视为对任何产品或流程质量的最终判定，也不应用于关键应用场景、商业决策或法规合规用途，除非经过专业核查和独立验证。对于分析结果的正确性、准确性、可靠性或适用性，不作任何明示或暗示的保证。"

  
  # --- 2. 创建 VCF 质量控制表格 ---
  
  # 获取原始数据 (从 qc_result$vcf_table 中)
  # 假设 vcf_table 包含列: Sample, SNV, INDEL, SNV precision, INDEL precision, SNV recall, INDEL recall
  # 注意: qc.R 中已经将 precision/recall 除以 100 转换为 0-1 小数
  raw_df <- qc_result$vcf_table
  
  # 辅助函数: 格式化数值并检查阈值
  # 如果低于阈值，添加 " ↓"
  fmt_val <- function(val, threshold) {
    s <- sprintf("%.3f", val)
    if (!is.na(val) && val < threshold) {
      return(paste0(s, " ↓"))
    }
    return(s)
  }
  
  # 准备数据行
  data_rows <- list()
  
  for (i in 1:nrow(raw_df)) {
    # 获取数值
    # 注意列名可能包含空格，需根据实际 qc.R 输出调整，这里使用标准名称
    v_sp <- raw_df$`SNV precision`[i]
    v_ip <- raw_df$`INDEL precision`[i]
    v_sr <- raw_df$`SNV recall`[i]
    v_ir <- raw_df$`INDEL recall`[i]
    
    # 获取计数值 (处理可能的不同列名情况，优先使用 'SNV number')
    # 添加逗号分隔符 (big.mark)
    c_snv <- if("SNV number" %in% names(raw_df)) {
      format(as.numeric(raw_df$`SNV number`[i]), big.mark=",")
    } else if ("SNV" %in% names(raw_df)) {
      format(as.numeric(raw_df$SNV[i]), big.mark=",")
    } else {
      "-"
    }
    
    c_indel <- if("INDEL number" %in% names(raw_df)) {
      format(as.numeric(raw_df$`INDEL number`[i]), big.mark=",")
    } else if ("INDEL" %in% names(raw_df)) {
      format(as.numeric(raw_df$INDEL[i]), big.mark=",")
    } else {
      "-"
    }
    
    # 判断是否全部通过
    # 标准: SNV P >= 0.99, INDEL P >= 0.90, SNV R >= 0.98, INDEL R >= 0.90
    is_pass <- (v_sp >= 0.99 && v_ip >= 0.90 && v_sr >= 0.98 && v_ir >= 0.90)
    overall_q <- ifelse(is_pass, "Yes", "No")
    
    # 构建行向量
    row_vec <- c(
      raw_df$Sample[i],
      c_snv,
      c_indel,
      fmt_val(v_sp, 0.99),
      fmt_val(v_ip, 0.90),
      fmt_val(v_sr, 0.98),
      fmt_val(v_ir, 0.90),
      overall_q
    )
    data_rows[[i]] <- row_vec
  }
  
  # 构建完整表格数据框
  # 第一行: 推荐质量标准
  rec_row <- c("推荐质量标准", "-", "-", "≥0.99", "≥0.90", "≥0.98", "≥0.90", "-")
  
  tbl_matrix <- do.call(rbind, data_rows)
  final_df <- rbind(rec_row, tbl_matrix)
  
  # 设置列名
  final_df <- as.data.frame(final_df, stringsAsFactors = FALSE)
  colnames(final_df) <- c("样本", "#SNV", "#INDEL", "SNV精确度", "INDEL精确度", "SNV 灵敏度", "INDEL灵敏度", "是否通过")
  
  # 生成 Flextable
  ft <- flextable(final_df) %>%
    theme_box() %>%
    align(align = "center", part = "all") %>%
    # 将整个表格设置为 Times New Roman。
    flextable::font(part = "all", fontname = "Times New Roman") %>%
    # --- 修改列宽逻辑 ---
    # 第1列(样本名)给 2 英寸 (足够长)
    width(j = 1, width = 1.8) %>% 
    # 第2-8列(数值)给 0.6 英寸 (紧凑)
    width(j = 2, width = 0.6) %>%
    # width(j = 3, width = 0.7) %>%
    # width(j = 4, width = 0.6) %>%
    # width(j = 5, width = 0.7) %>%
    width(j = 3:7, width = 0.7) %>%
    # 第8列(数值)给 0.5 英寸 (紧凑)
    width(j = 8, width = 0.8) %>%
    # ------------------
  bold(part = "header") %>%
    bold(i = 1, part = "body") %>% # 加粗第一行推荐标准
    bg(i = 1, bg = "#EFEFEF", part = "body") %>% # 第一行背景灰
    # 调整字号，防止表格过宽换行
    fontsize(part = "all", size = 10) %>%
    # 动态标红: 是否通过为 No
    color(j = "是否通过", i = ~ `是否通过` == "No", color = "#B80D0D") %>%
    # 动态标红: 数值带有 ↓ 符号
    color(j = "SNV精确度", i = ~ grepl("↓", `SNV精确度`), color = "#B80D0D") %>%
    color(j = "INDEL精确度", i = ~ grepl("↓", `INDEL精确度`), color = "#B80D0D") %>%
    color(j = "SNV 灵敏度", i = ~ grepl("↓", `SNV 灵敏度`), color = "#B80D0D") %>%
    color(j = "INDEL灵敏度", i = ~ grepl("↓", `INDEL灵敏度`), color = "#B80D0D")

  
  # --- 3. 生成报告文档 ---
  read_docx(report_template) %>%
    # 标题
    body_add_par(value = "Quartet基因组质量报告", style = "heading 1") %>%
    
    # 摘要
    body_add_par(value = "摘要", style = "heading 2") %>%
    body_add_par(value = text_sum_intro, style = "Normal") %>%
    body_add_par(value = " ", style = "Normal") %>%
    
    # 插入表格
    body_add_flextable(ft) %>%
    # body_add_break() %>% # 根据需要分页
    
    # 定义部分
    body_add_par(value = text_def_title, style = "heading 2") %>%
    
    body_add_par(value = text_prec_title, style = "heading 3") %>%
    body_add_par(value = text_prec_desc, style = "Normal") %>%
    
    body_add_par(value = text_rec_title, style = "heading 3") %>%
    body_add_par(value = text_rec_desc, style = "Normal") %>%
    
    # 参考文献
    body_add_par(value = text_ref_title, style = "heading 2") %>%
    body_add_par(value = text_ref_1, style = "Normal") %>%
    body_add_par(value = text_ref_2, style = "Normal") %>%
    body_add_par(value = text_ref_3, style = "Normal") %>%
    body_add_par(value = text_ref_4, style = "Normal") %>%
    body_add_par(value = " ", style = "Normal") %>%
    
    # 免责声明
    body_add_par(value = text_disc_title, style = "heading 3") %>%
    body_add_par(value = text_disc_content, style = "Normal") %>%
    
    # 输出
    print(target = output_file)
}
