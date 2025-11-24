# --------------------------------------------------------------------------- #
#' @title Generate Quartet Genomics report
#'
#' @description Use calculated Met result to generate report
#'
#' @param dna_result list
#' @param temp_doc_path character
#' @param output_path character
#'
#' @return word file
#'
#' @importFrom dplyr %>%
#' @importFrom flextable flextable
#' @importFrom flextable theme_vanilla
#' @importFrom flextable color
#' @importFrom flextable set_caption
#' @importFrom flextable align
#' @importFrom flextable width
#' @importFrom flextable bold
#' @importFrom flextable nrow_part
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
#'

generate_dna_report <- function(dna_result = NULL,
                                doc_file_path = NULL,
                                output_path = NULL) {
  if (is.null(dna_result) || is.null(doc_file_path)) {
    stop("All arguments (dna_result, doc_file_path) are required.")
  }

  if (is.null(output_path)) {
    path <- getwd()
    sub_dir <- "output"
    dir.create(file.path(path, sub_dir), showWarnings = FALSE)
    output_path <- file.path(path, "output")
  }
  output_file <- file.path(output_path, "Quartet_DNA_Report.docx")

  ### 创建Evaluate Metrics 表格

  ## all metrics table
  summary_ft <- flextable(dna_result$conclusion)
  summary_ft <- summary_ft %>%
    color(~ Performance == "Bad", color = "#B80D0D", ~Performance) %>%
    color(~ Performance == "Fair", color = "#D97C11", ~Performance) %>%
    color(~ Performance == "Good", color = "#70C404", ~Performance) %>%
    color(~ Performance == "Great", color = "#0F9115", ~Performance) %>%
    width(width = 1.25) %>%
    align(align = "center", part = "all") %>%
    bold(i = nrow_part(., part = "body"), part = "body") %>%
    bold(i = 1, part = "header")

  ## vcf table
  vcf_ft <- flextable(dna_result$vcf_table)
  vcf_ft <- vcf_ft %>%
    bold(i = 1, part = "header") %>%
    width(j = 1, width = 2) %>%
    align(align = "center", part = "all")

  mendelian_ft <- dna_result$mendelian_table
  if (!is.null(mendelian_ft)) {
    mendelian_ft <- flextable(mendelian_ft)
    mendelian_ft <- mendelian_ft %>%
      bold(i = 1, part = "header") %>%
      width(j = 1, width = 2) %>%
      align(align = "center", part = "all")
  }

  ### 绘制Total score 历史分数排名散点图
  historical_rank <- dna_result$rank_table[, c("batch", "total_norm")]
  historical_rank_ref <- historical_rank[historical_rank$batch != "Queried_Data", "total_norm"]

  q1 <- quantile(historical_rank_ref, probs = 0.2, na.rm = TRUE)
  q2 <- quantile(historical_rank_ref, probs = 0.5, na.rm = TRUE)
  q3 <- quantile(historical_rank_ref, probs = 0.8, na.rm = TRUE)

  p_rank_scatter_plot <- ggplot(data = historical_rank) +
    # 添加四个区域
    # geom_rect(aes(xmin = 1, xmax = q1, ymin = -Inf, ymax = Inf), fill = "#B80D0D", alpha = 0.08) +
    # geom_rect(aes(xmin = q1, xmax = q2, ymin = -Inf, ymax = Inf), fill = "#D97C11", alpha = 0.08) +
    # geom_rect(aes(xmin = q2, xmax = q3, ymin = -Inf, ymax = Inf), fill = "#70C404", alpha = 0.08) +
    # geom_rect(aes(xmin = q3, xmax = 10, ymin = -Inf, ymax = Inf), fill = "#0F9115", alpha = 0.08) +
    geom_rect(
      xmin = 1,
      xmax = q1,
      ymin = -Inf,
      ymax = Inf,
      fill = "#B80D0D",
      alpha = 0.08
    ) +
    geom_rect(
      xmin = q1,
      xmax = q2,
      ymin = -Inf,
      ymax = Inf,
      fill = "#D97C11",
      alpha = 0.08
    ) +
    geom_rect(
      xmin = q2,
      xmax = q3,
      ymin = -Inf,
      ymax = Inf,
      fill = "#70C404",
      alpha = 0.08
    ) +
    geom_rect(
      xmin = q3,
      xmax = 10,
      ymin = -Inf,
      ymax = Inf,
      fill = "#0F9115",
      alpha = 0.08
    ) +
    # 添加基础点图层
    geom_point(aes(x = total_norm, y = reorder(batch, total_norm))) +
    # 突出显示 "QUERIED DATA" 对应的点
    geom_point(
      data = subset(historical_rank, batch == "Queried_Data"),
      aes(x = total_norm, y = reorder(batch, total_norm)),
      color = "orange",
      size = 3
    ) +
    # 自定义x轴刻度
    scale_x_continuous(breaks = c(
      1,
      round(as.numeric(q1), 2),
      round(as.numeric(q2), 2),
      round(as.numeric(q3), 2),
      10
    )) +
    theme_minimal() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text = element_text(face = "bold"),
      plot.title = element_text(
        hjust = 0.5,
        size = 16,
        face = "bold"
      )
    ) +
    labs(x = " ", y = " ", title = "Total Score")


  #### 设置输出文本

  ###### 第一部分
  ### Summary

  text_sum_intro <- "This report summarizes the quality of the data generated from Quartet DNA reference materials based on several key quality control (QC) metrics. Each metric is accompanied by its current value, historical average, ranking among all datasets evaluated, and the corresponding performance grade."

  text_1 <- "The submitted data will be graded as Bad, Fair, Good, or Great, depending on how the total score stands against historical data. Total score = (1+0.5^2) x SNV_score x INDEL_score / (0.5^2 x SNV_score + INDEL_score). SNV_score and INDEL_score are obtained by calculating the mean values of Precision, Recall, and MCR, respectively."
  ### Four levels of performance
  text_1_sup_1 <- "Based on the scaled total score, the submitted data will be ranked together with all Quartet historical datasets. The higher the score, the higher the ranking. After this, the performance levels will be assigned based on their ranking ranges."
  text_1_sup_2 <- fpar(
    ftext("· Bad: ", fp_text(bold = TRUE)),
    ftext("Lowest quintile (0-20th percentile).", fp_text())
  )
  text_1_sup_3 <- fpar(
    ftext("· Fair: ", fp_text(bold = TRUE)),
    ftext("Lower middle quartile (21st-50th percentile).", fp_text())
  )
  text_1_sup_4 <- fpar(
    ftext("· Good: ", fp_text(bold = TRUE)),
    ftext("Upper middle quartile (51st-80th percentile).", fp_text())
  )
  text_1_sup_5 <- fpar(
    ftext("· Great: ", fp_text(bold = TRUE)),
    ftext("Highest quintile (81st-100th percentile).", fp_text())
  )


  #### 第二部分 Quality control metric

  ### Performance Score
  text_2 <- "Scores of evaluation metrics for the current batch and all historical batches assessed. For better comparison and presentation, the total score was scaled to the interval [1, 10], with the worst dataset being 1 and the best dataset scoring 10. Please note that the results shown here are scaled values for all batches in each metric. The name of your data is Queried_Data."
  ### Signal-to-Noise Ratio
  text_3 <- "The SNV performance of evaluated data compared to the Quartet historical batches is shown in this section. Each data point represents a set of Quartet samples, i.e., one each of D5, D6, F7, and M8."
  ### Correlation with Reference Datasets
  text_4 <- "The Indel performance of evaluated data compared to the Quartet historical batches is shown in this section. Each data point represents a set of Quartet samples, i.e., one each of D5, D6, F7, and M8."
  ###
  text_5 <- " "
  ###
  text_6 <- "Each row represents a set of Quartet samples, i.e. one each of D5, D6, F7 and M8. When multiple sets of technical replicates are measured, the performance of each set will be represented by row."


  ### Method
  supplementary_info_1 <- "The QC pipeline starts from the variant calling file, enabling the calculation of the following three metrics. The total score is an F0.5-measure of the SNV score and the INDEL score, which are the mean values of Precision, Recall, and MCR, respectively."
  supplementary_info_1_1 <- "Tested call sets were compared with benchmark small variants using hap.py (https://github.com/Illumina/hap.py). Precision is the fraction of called variants in the test dataset that are true."
  supplementary_info_1_2 <- "Recall is the fraction of true variants are called in the test dataset."
  supplementary_info_1_3 <- "Mendelian concordance rate (MCR) is the number of variants following Mendelian inheritance laws divided by the total number of variants called among the four Quartet samples. Mendelian concordant variants are the variants shared by the twins (D5 and D6) and following Mendelian inheritance laws with parents (Father: F7 and Mother M8). Mendelian analysis was performed using VBT (https://github.com/sbg/VBT-TrioAnalysis). When calculating Mendelian concordance rate of small variants, variants on large deletions were not included, because VBT takes these variants as Mendelian violations."
  # supplementary_info_1_4 = "Pearson correlation coefficient reflects the overall reproducibility within replicates. We calculate correlation coefficients between each two replicates within each biological sample, and take the median as the definitive absolute correlation."
  # supplementary_info_1_5 = "SNR is established to characterize the capability of a platform or lab or batch to differentiate the intrinsic differences among distinct biological sample groups ('signal') from variations in technical replicates of the same sample group ('noise')."
  # supplementary_info_1_6 = "RC is employed to evaluate the quantitative agreement with the reference dataset on a relative scale. Since protein quantities are reassembled by peptide intensities, the reference dataset is established at ratio-scale quantified peptide level across sample pairs (D5/D6, F7/D6, M8/D6), based on historical datasets. To ensure a robust assessment, we only consider peptides that are shared with reference datasets and whose log2FCs reached statistical significance (using limma R package; FDR adjusted P < 0.05). The RC value is subsequently determined by the Pearson correlation coefficient between the test dataset and the reference dataset."

  ### Reference
  supplementary_info_ref1 <- "1. Zheng Y, Liu Y, Yang J, et al. Multi-omics data integration using ratio-based quantitative profiling with Quartet reference materials. Nat Biotechnol. Published online September 7, 2023."
  supplementary_info_ref2 <- "2. Ren L, Duan X, Dong L, et al. Quartet DNA reference materials and datasets for comprehensively evaluating germline variant calling performance. Genome Biol. 2023; 24(1): 270."

  ### Contact us
  # supplementary_info_2_1 = "Fudan University Pharmacogenomics Research Center"
  # supplementary_info_2_2 = "Project manager: Quartet Team"
  # supplementary_info_2_3 = "Email: quartet@fudan.edu.cn"

  ### Disclaimer
  supplementary_info_3 <- 'This Data Quality Report is provided as an analysis of the specific dataset evaluated and is intended for informational purposes only. While every effort has been made to ensure the accuracy and reliability of the analysis, the information is presented "AS IS" without warranty of any kind, either express or implied. The authors and distributors of this report shall not be held liable for any actions taken in reliance thereon. Users are advised that the findings within this report are not to be used as definitive statements on the quality of any product or process beyond the scope of the dataset assessed. This report is not intended for use in critical applications, commercial decision-making, or for regulatory compliance without professional verification and independent validation. No guarantee, either expressed or implied, is made regarding the use or results of the analysis, including without limitation, the correctness, accuracy, reliability, or applicability of the findings.'

  read_docx(doc_file_path) %>%
    ## 添加报告标题
    body_add_par(value = "Quartet Report for Genomics", style = "heading 1") %>%
    ## 第一部分，Assessment Summary
    body_add_par(value = "Summary", style = "heading 2") %>%
    body_add_par(value = " ", style = "Normal") %>%
    body_add_par(value = text_sum_intro, style = "Normal") %>%
    body_add_par(value = " ", style = "Normal") %>%
    body_add_flextable(summary_ft) %>%
    body_add_break() %>%
    ### 第二部分 Quality control metric
    body_add_par(value = "QC Metrics", style = "heading 2") %>%
    body_add_par(value = supplementary_info_1, style = "Normal") %>%
    body_add_par(value = "Precision:", style = "heading 3") %>%
    body_add_par(value = supplementary_info_1_1, style = "Normal") %>%
    body_add_par(value = "Recall:", style = "heading 3") %>%
    body_add_par(value = supplementary_info_1_2, style = "Normal") %>%
    body_add_par(value = "Mendelian Concordance Rate:", style = "heading 3") %>%
    body_add_par(value = supplementary_info_1_3, style = "Normal") %>%
    # body_add_par(value = "Absolute correlation:",style = "heading 3") %>%
    # body_add_par(value = supplementary_info_1_4,style = "Normal") %>%
    # body_add_par(value = "Signal-to-Noise Ratio (SNR):",style = "heading 3") %>%
    # body_add_par(value = supplementary_info_1_5,style = "Normal") %>%
    # body_add_par(value = "Relative Correlation with the Reference Dataset (RC):",style = "heading 3") %>%
    # body_add_par(value = supplementary_info_1_6,style = "Normal") %>%

    body_add_par(value = "Total Score:", style = "heading 3") %>%
    body_add_par(value = text_1, style = "Normal") %>%
    ## 分页
    # body_add_break()%>%

    body_add_par(value = "Performance Category:", style = "heading 3") %>%
    body_add_par(value = text_1_sup_1, style = "Normal") %>%
    body_add_fpar(value = text_1_sup_2, style = "Normal") %>%
    body_add_fpar(value = text_1_sup_3, style = "Normal") %>%
    body_add_fpar(value = text_1_sup_4, style = "Normal") %>%
    body_add_fpar(value = text_1_sup_5, style = "Normal") %>%
    ### 排名散点图
    body_add_par(value = "Performance Grade", style = "heading 2") %>%
    body_add_gg(value = p_rank_scatter_plot, style = "centered") %>%
    body_add_par(value = text_2, style = "Normal") %>%
    ## 分页
    body_add_break() %>%
    {
      if (!is.null(mendelian_ft)) {
        ## SNV 散点图
        body_add_par(x = ., value = "Performance of SNV", style = "heading 2") %>%
          body_add_gg(dna_result$p_mendelian_f1_snv, style = "centered") %>%
          body_add_par(value = text_3, style = "Normal") %>%
          ## 分页
          body_add_break() %>%
          ## indel 散点图
          body_add_par(value = "Performance of Indel", style = "heading 2") %>%
          body_add_gg(dna_result$p_mendelian_f1_indel, style = "centered") %>%
          body_add_par(value = text_4, style = "Normal") %>%
          ## 分页
          body_add_break()
      } else {
        .
      }
    } %>%
    ## vcf qc
    body_add_par(value = "Variant Calling Quality Control", style = "heading 2") %>%
    body_add_par(value = "Details based on reference datasets", style = "heading 3") %>%
    body_add_flextable(vcf_ft) %>%
    {
      if (!is.null(mendelian_ft)) {
        body_add_par(x = ., value = text_5, style = "Normal") %>%
          ## mendelian qc
          body_add_par(value = "Details based on Quartet genetic built-in truth", style = "heading 3") %>%
          body_add_flextable(mendelian_ft)
      } else {
        .
      }
    } %>%
    body_add_par(value = text_6, style = "Normal") %>%
    ## 分页
    body_add_break() %>%
    ### 附加信息
    body_add_par(value = "Supplementary Information", style = "heading 2") %>%
    body_add_par(value = "Reference", style = "heading 3") %>%
    body_add_par(value = supplementary_info_ref1, style = "Normal") %>%
    body_add_par(value = supplementary_info_ref2, style = "Normal") %>%
    # body_add_par(value = supplementary_info_ref3, style = "Normal") %>%
    # body_add_par(value = supplementary_info_ref4, style = "Normal") %>%
    # body_add_par(value = "Contact us", style = "heading 3") %>%
    # body_add_par(value = supplementary_info_2_1, style = "Normal") %>%
    # body_add_par(value = supplementary_info_2_2, style = "Normal") %>%
    # body_add_par(value = supplementary_info_2_3, style = "Normal") %>%
    body_add_par(value = "Disclaimer", style = "heading 3") %>%
    body_add_par(value = supplementary_info_3, style = "Normal") %>%
    ## 输出文件
    print(target = output_file)
}
