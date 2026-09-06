# =============================================================================
#   Figure 3 - simulation: effect of the input summary statistics
# =============================================================================
#
#  Reads ./Simulation/Sim_CRC_stage2/, writes ./Simulation/Figure/Fig3_*.png
#  Run from the project root: Rscript Rscript/Figure3.R
# =============================================================================

  # Packages ----
  library('ggplot2')
  library("tidyverse")
  library("latex2exp")
  library("ggtext")
  
  rm(list = ls())
  
  ########################### Stage 1 simulation results ###########################
  ## Proportion of true signatures
  data.loc.L10 <- "./Simulation/Sim_CRC_stage2/res_Ka"
  
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
  PRC_all$Method <- factor(PRC_all$method, levels = c("SMESH-PALM", "SMESH-ANCOMBC2", "SMESH-LinDA", "SMESH-MaAsLin3"), 
                           labels = c("PALM", "ANCOMBC2", "LinDA", "MaAsLin3"))
  
  ## Figure A
  PRC_sum <- PRC_all %>%
    group_by(scenario, Ka, Method) %>%
    summarise(
      mean_G = mean(G_hat, na.rm = TRUE),
      sd_G   = sd(G_hat,  na.rm = TRUE),
      mean_P = mean(Precision, na.rm = TRUE),
      sd_P   = sd(Precision,  na.rm = TRUE),
      mean_R = mean(Recall, na.rm = TRUE),
      sd_R   = sd(Recall,  na.rm = TRUE),
      mean_F1 = mean(F1, na.rm = TRUE),
      sd_F1   = sd(F1,  na.rm = TRUE),
      mean_ARI = mean(ARI, na.rm = TRUE),
      sd_ARI   = sd(ARI,  na.rm = TRUE),
      .groups = "drop"
    )

  true_G_df <- data.frame(
    scenario = paste0("Scenario ", 1:5),
    true_G = c(1, 2, 3, 2, 3)
  )
  
  ## Figure A ----
  ## Generate figure result for all signature
  p_ARI <- PRC_sum %>%
    ggplot(aes(x = Ka, y = mean_ARI, color = Method, group = Method)) +
    geom_line(linewidth = 0.8, position = position_dodge(width = 0.4)) +
    geom_point(size = 1.5, position = position_dodge(width = 0.4)) +
    geom_errorbar(
      aes(ymin = pmax(mean_ARI - sd_ARI, 0), ymax = pmin(mean_ARI + sd_ARI, 1)),
      width = 0.15,
      position = position_dodge(width = 0.4)
    ) +
    facet_grid(~scenario) +
    ylim(0, 1) +
    xlab("Proportion of true signatures") +
    ylab("ARI") +
    scale_color_manual(
      name = "SMESH input summary statistics",
      breaks = c("PALM", "ANCOMBC2", "LinDA", "MaAsLin3"),
      values = c("red", "brown", "pink", "#C77CFF")
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
      
      legend.position = "none",
      legend.box = "vertical",
      legend.title = element_blank(),
      legend.text = element_text(size = 18),
      strip.text = element_blank(),
      legend.key.width = unit(1.5, "cm")
    ) +
    guides(color = guide_legend(order = 1, nrow = 1))
    
  p_G <- PRC_sum %>%
    ggplot(aes(x = Ka, y = mean_G, color = Method, group = Method)) +
    geom_line(linewidth = 0.8, position = position_dodge(width = 0.4)) +
    geom_point(size = 1.5, position = position_dodge(width = 0.4)) +
    geom_errorbar(
      aes(ymin = pmax(mean_G - sd_G, 0), ymax = pmin(mean_G + sd_G, 5)),
      width = 0.15,
      position = position_dodge(width = 0.4)
    ) +
    geom_hline(
      data = true_G_df,
      aes(yintercept = true_G),
      inherit.aes = FALSE,
      linetype = "dashed",
      color = "grey",
      linewidth = 0.5
    ) +
    facet_grid(~scenario) +
    ylim(0, 5) +
    xlab("Proportion of true signatures") +
    ylab(expression("Selected " * italic(G))) +
    scale_color_manual(
      name = "SMESH input summary statistics",
      breaks = c("PALM", "ANCOMBC2", "LinDA", "MaAsLin3"),
      values = c("red", "brown", "pink", "#C77CFF")
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
      
      legend.position = "none",
      legend.box = "vertical",
      legend.title = element_text(size = 18),
      legend.text = element_text(size = 18),
      strip.text = element_markdown(size = 18),
      legend.key.width = unit(1.5, "cm")
    ) +
    guides(color = guide_legend(order = 1, nrow = 1))
  
  ggsave(
    filename = "./Simulation/Figure/Fig3_A1.png",
    plot = p_ARI,
    width = 350,
    height = 90,
    units = "mm",
    dpi = 300
  )
  
  ggsave(
    filename = "./Simulation/Figure/Fig3_A2.png",
    plot = p_G,
    width = 350,
    height = 80,
    units = "mm",
    dpi = 300
  )
  
  ## Figure B ----
  ## Generate figure result for all signature
  p.precison <- PRC_sum %>%
    ggplot(aes(x = Ka, y = mean_P, color = Method, group = Method)) +
    geom_line(linewidth = 0.8, position = position_dodge(width = 0.4)) +
    geom_point(size = 1.5, position = position_dodge(width = 0.4)) +
    geom_errorbar(
      aes(ymin = pmax(mean_P - sd_P, 0.4), ymax = pmin(mean_P + sd_P, 1)),
      width = 0.15,
      position = position_dodge(width = 0.4)
    ) +
    facet_grid(~scenario) +
    scale_y_continuous(
      limits = c(0.4, 1),
      breaks = c(0.4, 0.6, 0.8, 1.0),
      labels = scales::number_format(accuracy = 0.01)
    ) +
    xlab("Proportion of true signatures") +
    ylab("Precision") +
    scale_color_manual(
      name = "SMESH input summary statistics",
      breaks = c("PALM", "ANCOMBC2", "LinDA", "MaAsLin3"),
      values = c("red", "brown", "pink", "#C77CFF")
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
      
      legend.position = "none",
      legend.box = "vertical",
      legend.title = element_text(size = 18),
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
      aes(ymin = pmax(mean_R - sd_R, 0.4), ymax = pmin(mean_R + sd_R, 1)),
      width = 0.15,
      position = position_dodge(width = 0.4)
    ) +
    facet_grid(~scenario) +
    scale_y_continuous(
      limits = c(0.4, 1),
      breaks = c(0.4, 0.6, 0.8, 1.0),
      labels = scales::number_format(accuracy = 0.01)
    ) +
    xlab("Proportion of true signatures") +
    ylab("Recall") +
    scale_color_manual(
      name = "SMESH input summary statistics",
      breaks = c("PALM", "ANCOMBC2", "LinDA", "MaAsLin3"),
      values = c("red", "brown", "pink", "#C77CFF")
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
                           legend = "none")
  
  ggsave(
    filename = "./Simulation/Figure/Fig3_B.png",
    plot = g_B,
    width = 350,
    height = 150,
    units = "mm",
    dpi = 300
  )
  