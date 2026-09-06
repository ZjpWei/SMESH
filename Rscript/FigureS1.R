# =============================================================================
#   Figure S1 - simulation: signature-type-specific detection, benchmark methods
# =============================================================================
#
#  Reads ./Simulation/Sim_CRC_stage1/, writes ./Simulation/Figure/SFig1_*.png
#  Run from the project root: Rscript Rscript/FigureS1.R
# =============================================================================

  library('ggplot2')
  library("tidyverse")
  library("latex2exp")
  library("ggtext")
  
  rm(list = ls())

  ## Aggregate results ----
  data.loc.L10 <- "./Simulation/Sim_CRC_stage1/res_Ka"
  
  PRC_all <- NULL
  for(Ka in c(0.1, 0.2, 0.3, 0.4)){
    for(scenario in 1:5){
      PRC_tmp <- NULL
      for(s in 1:50){
        if(file.exists(paste0(data.loc.L10, Ka, "_pos0.6_u0_scenario", scenario, "_s", s, ".Rdata"))){
          load(paste0(data.loc.L10, Ka, "_pos0.6_u0_scenario", scenario, "_s", s, ".Rdata"))
          PRC_tmp <- rbind(PRC_tmp, result_mat)
        }
      }
      PRC_tmp$L <- paste0("L = 10, ", PRC_tmp$data)
      PRC_tmp$scenario <- paste0("Scenario ", scenario)
      PRC_tmp$Ka <- as.character(Ka)
      PRC_all <- rbind(PRC_all, PRC_tmp)
    }
  }
  PRC_all$Method <- factor(PRC_all$method, levels =  c("SMESH", "True-cluster FE", "SKM + FE", "SHC + FE", "Melody"))
  
  PRC_sum <- PRC_all %>%
    group_by(scenario, Ka, Method) %>%
    summarise(
      mean_P = mean(Precision_hom, na.rm = TRUE),
      sd_P   = sd(Precision_hom,  na.rm = TRUE),
      mean_R = mean(Recall_hom, na.rm = TRUE),
      sd_R   = sd(Recall_hom,  na.rm = TRUE),
      .groups = "drop"
    ) 
  
  ## Generate figure result for all signature (FigS1 a) ----
  p.precison <- PRC_sum %>%
    ggplot(aes(x = Ka, y = mean_P, color = Method, group = Method)) +
    geom_line(linewidth = 0.8, position = position_dodge(width = 0.4)) +
    geom_point(size = 1.5, position = position_dodge(width = 0.4)) +
    geom_errorbar(
      aes(ymin = pmax(mean_P - sd_P, 0), ymax = pmin(mean_P + sd_P, 1)),
      width = 0.15,
      position = position_dodge(width = 0.4)
    ) +
    facet_grid(~scenario) +
    ylim(0, 1) +
    xlab("Proportion of true signatures") +
    ylab("Precision") +
    scale_color_manual(
      name = "Summary statistics",
      breaks = c("SMESH", "True-cluster FE", "SKM + FE", "SHC + FE", "Melody"),
      values = c("red", "orange", "skyblue", "blue","#4dac26")
    ) +
    theme_bw() +
    theme(
      plot.title = element_blank(),
      axis.title.y = element_text(size = 18),
      axis.title.x = element_blank(),
      axis.text.y = element_text(size = 18),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.line = element_line(colour = "black"),
      
      panel.grid.major = element_line(colour = "white", linetype = "dotted"),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(colour = "black", fill = NA),
      panel.background = element_rect(fill = "white"),
      
      legend.position = "bottom",
      legend.box = "vertical",
      legend.title = element_blank(),
      legend.text = element_text(size = 18),
      strip.text = element_markdown(size = 18),
      legend.key.width = unit(1.5, "cm")
    ) +
    guides(color = guide_legend(order = 1, nrow = 1))
  
  p.recall <- PRC_sum %>%
    ggplot(aes(x = Ka, y = mean_R, color = Method, group = Method)) +
    geom_line(linewidth = 0.8, position = position_dodge(width = 0.4)) +
    geom_point(size = 1.5, position = position_dodge(width = 0.4)) +
    geom_errorbar(
      aes(ymin = pmax(mean_R - sd_R, 0), ymax = pmin(mean_R + sd_R, 1)),
      width = 0.15,
      position = position_dodge(width = 0.4)
    ) +
    facet_grid(~scenario) +
    ylim(0, 1) +
    xlab("Proportion of true signatures") +
    ylab("Recall") +
    scale_color_manual(
      name = "Summary statistics",
      breaks = c("SMESH", "True-cluster FE", "SKM + FE", "SHC + FE", "Melody"),
      values = c("red", "orange", "skyblue", "blue","#4dac26")
    ) +
    theme_bw() +
    theme(
      plot.title = element_blank(),
      axis.title.y = element_text(size = 18),
      axis.title.x = element_text(size = 18),
      axis.text = element_text(size = 18),
      axis.line = element_line(colour = "black"),
      
      panel.border = element_rect(colour = "black", fill = NA),
      panel.background = element_rect(fill = "white"),
      panel.grid.major = element_line(colour = "white", linetype = "dotted"),
      panel.grid.minor = element_blank(),
      
      legend.position = "bottom",
      legend.box = "vertical",
      legend.title = element_blank(),
      legend.text = element_text(size = 18),
      strip.text = element_blank(),
      legend.key.width = unit(1.5, "cm")
    ) +
    guides(color = guide_legend(order = 1, nrow = 1))
  
  g_A <- ggpubr::ggarrange(p.precison, p.recall, ncol = 1, common.legend = TRUE,  
                           heights = c(1, 1.1), legend = "none")
  
  ggsave(
    filename = "./Simulation/Figure/SFig1_A.png",
    plot = g_A,
    width = 350,
    height = 150,
    units = "mm",
    dpi = 300
  )
  
  ## Supplementary Figure1 B
  PRC_sum <- PRC_all %>%
    group_by(scenario, Ka, Method) %>%
    summarise(
      mean_P = mean(Precision_het, na.rm = TRUE),
      sd_P   = sd(Precision_het,  na.rm = TRUE),
      mean_R = mean(Recall_het, na.rm = TRUE),
      sd_R   = sd(Recall_het,  na.rm = TRUE),
      .groups = "drop"
    ) %>% dplyr::filter(scenario != "Scenario 1")
  
  ## Generate figure result for all signature
  p.precison <- PRC_sum %>%
    ggplot(aes(x = Ka, y = mean_P, color = Method, group = Method)) +
    geom_line(linewidth = 0.8, position = position_dodge(width = 0.4)) +
    geom_point(size = 1.5, position = position_dodge(width = 0.4)) +
    geom_errorbar(
      aes(ymin = pmax(mean_P - sd_P, 0), ymax = pmin(mean_P + sd_P, 1)),
      width = 0.15,
      position = position_dodge(width = 0.4)
    ) +
    facet_grid(~scenario) +
    ylim(0, 1) +
    xlab("Proportion of true signatures") +
    ylab("Precision") +
    scale_color_manual(
      name = "Summary statistics",
      breaks = c("SMESH", "True-cluster FE", "SKM + FE", "SHC + FE", "Melody"),
      values = c("red", "orange", "skyblue", "blue","#4dac26")
    ) +
    theme_bw() +
    theme(
      plot.title = element_blank(),
      axis.title.y = element_text(size = 18),
      axis.title.x = element_blank(),
      axis.text.y = element_text(size = 18),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.line = element_line(colour = "black"),
      
      panel.grid.major = element_line(colour = "white", linetype = "dotted"),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(colour = "black", fill = NA),
      panel.background = element_rect(fill = "white"),
      
      legend.position = "bottom",
      legend.box = "vertical",
      legend.title = element_blank(),
      legend.text = element_text(size = 18),
      strip.text = element_markdown(size = 18),
      legend.key.width = unit(1.5, "cm")
    ) +
    guides(color = guide_legend(order = 1, nrow = 1))
  
  p.recall <- PRC_sum %>%
    ggplot(aes(x = Ka, y = mean_R, color = Method, group = Method)) +
    geom_line(linewidth = 0.8, position = position_dodge(width = 0.4)) +
    geom_point(size = 1.5, position = position_dodge(width = 0.4)) +
    geom_errorbar(
      aes(ymin = pmax(mean_R - sd_R, 0), ymax = pmin(mean_R + sd_R, 1)),
      width = 0.15,
      position = position_dodge(width = 0.4)
    ) +
    facet_grid(~scenario) +
    ylim(0, 1) +
    xlab("Proportion of true signatures") +
    ylab("Recall") +
    scale_color_manual(
      name = "Summary statistics",
      breaks = c("SMESH", "True-cluster FE", "SKM + FE", "SHC + FE", "Melody"),
      values = c("red", "orange", "skyblue", "blue","#4dac26")
    ) +
    theme_bw() +
    theme(
      plot.title = element_blank(),
      axis.title.y = element_text(size = 18),
      axis.title.x = element_text(size = 18),
      axis.text = element_text(size = 18),
      axis.line = element_line(colour = "black"),
      
      panel.border = element_rect(colour = "black", fill = NA),
      panel.background = element_rect(fill = "white"),
      panel.grid.major = element_line(colour = "white", linetype = "dotted"),
      panel.grid.minor = element_blank(),
      
      legend.position = "bottom",
      legend.box = "vertical",
      legend.title = element_blank(),
      legend.text = element_text(size = 18),
      strip.text = element_blank(),
      legend.key.width = unit(1.5, "cm")
    ) +
    guides(color = guide_legend(order = 1, nrow = 1))
  
  g_B <- ggpubr::ggarrange(p.precison, p.recall, ncol = 1, common.legend = TRUE,  
                           heights = c(1, 1.1), legend = "none")
  
  ggsave(
    filename = "./Simulation/Figure/SFig1_B.png",
    plot = g_B,
    width = 280,
    height = 150,
    units = "mm",
    dpi = 300
  )
  
  