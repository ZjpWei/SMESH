# =============================================================================
#   Figure 2 - simulation: benchmark-method comparison
# =============================================================================
#
#  Reads ./Simulation/Sim_CRC_stage1/, writes ./Simulation/Figure/Fig2_*.png
#  Run from the project root: Rscript Rscript/Figure2.R
# =============================================================================

  library('ggplot2')
  library("tidyverse")
  library("latex2exp")
  library("ggtext")
  
  rm(list = ls())
  
  ## Aggregate simulation results ----
  data.loc.L10 <- "./Simulation/Sim_CRC_stage1/res_Ka"
  
  PRC_all <- NULL
  for(Ka in c(0.1, 0.2, 0.3, 0.4)){
    for(scenario in 1:5){
      PRC_tmp <- NULL
      for(s in 1:50){
        if(file.exists(paste0(data.loc.L10, Ka, "_pos0.6_u0_scenario", scenario, "_s", s, ".Rdata"))){
          load(paste0(data.loc.L10, Ka, "_pos0.6_u0_scenario", scenario, "_s", s, ".Rdata"))
          result_mat$linetype <- "median"
          PRC_tmp <- rbind(PRC_tmp, result_mat)
        }
      }
      PRC_tmp$L <- paste0("L = 10, ", PRC_tmp$data)
      PRC_tmp$scenario <- paste0("Scenario ", scenario)
      PRC_tmp$Ka <- as.character(Ka)
      PRC_all <- rbind(PRC_all, PRC_tmp)
    }
  }
  PRC_all$Method <- factor(PRC_all$method, levels = c("SMESH", "True-cluster FE", "SKM + FE", "SHC + FE", "Melody"))
  
  PRC_sum <- PRC_all %>%
    group_by(scenario, Ka, Method, linetype) %>%
    summarise(
      mean_P = mean(Precision, na.rm = TRUE),
      sd_P   = sd(Precision,  na.rm = TRUE),
      mean_R = mean(Recall, na.rm = TRUE),
      sd_R   = sd(Recall,  na.rm = TRUE),
      mean_F1 = mean(F1, na.rm = TRUE),
      sd_F1   = sd(F1,  na.rm = TRUE),
      mean_ari = mean(ARI, na.rm = TRUE),
      sd_ari   = sd(ARI,  na.rm = TRUE),
      .groups = "drop"
    )
  
  ## Generate ARI (Fig2 B) ----
  p.ari <- PRC_sum %>% filter(scenario != "Scenario 1") %>%
    ggplot(aes(x = Ka, y = mean_ari, color = Method, group = Method)) +
    geom_line(linewidth = 0.8, position = position_dodge(width = 0.4)) +
    geom_point(size = 1.5, position = position_dodge(width = 0.4)) +
    geom_errorbar(
      aes(ymin = mean_ari - sd_ari, ymax = pmin(mean_ari + sd_ari, 1.2)),
      width = 0.15,
      position = position_dodge(width = 0.4)
    ) +
    facet_grid(~scenario) +
    ylim(-1.2, 1.2) +
    xlab("Proportion of true signatures") +
    ylab("ARI") +
    scale_color_manual(
      name = "Summary statistics",
      breaks = c("SMESH", "True-cluster FE", "SKM + FE", "SHC + FE", "Melody"),
      values = c("red", "orange", "skyblue", "blue","#4dac26")
    ) +
    theme_bw() +
    theme(
      plot.title = element_blank(),
      axis.title.y = element_text(size = 18),
      axis.title.x = element_text(size = 18),,
      axis.text.y = element_text(size = 18),
      axis.text.x = element_text(size = 18),
      axis.line = element_line(colour = "black"),
      
      panel.grid.major = element_line(colour = "white", linetype = "dotted"),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(colour = "black", fill = NA),
      panel.background = element_rect(fill = "white"),
      
      legend.position = "none",
      legend.box = "vertical",
      legend.title = element_blank(),
      legend.text = element_text(size = 18),
      strip.text = element_markdown(size = 18),
      legend.key.width = unit(1.5, "cm")
    ) +  guides(color = guide_legend(order = 1, nrow = 1))
  
  ggsave(
    filename = "./Simulation/Figure/Fig2_b.png",
    plot = p.ari,
    width = 280,
    height = 90,
    units = "mm",
    dpi = 300
  )
  
  ## Overall signature detection (Fig2 c) ----
  p.precison <- PRC_sum %>%
    ggplot(aes(x = Ka, y = mean_P, color = Method, linetype = linetype,
               group = interaction(Method, linetype))) +
    geom_line(linewidth = 0.8, position = position_dodge(width = 0.4)) +
    geom_point(size = 1.5, position = position_dodge(width = 0.4)) +
    geom_errorbar(
      aes(ymin = pmax(mean_P - sd_P, 0.5), ymax = pmin(mean_P + sd_P, 1)),
      width = 0.15,
      position = position_dodge(width = 0.4)
    ) +
    facet_grid(~scenario) +
    ylim(0.5, 1) +
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
      
      panel.border = element_rect(colour = "black", fill = NA),
      panel.background = element_rect(fill = "white"),
      panel.grid.major = element_line(colour = "white", linetype = "dotted"),
      panel.grid.minor = element_blank(),
      
      legend.position = "bottom",
      legend.box = "vertical",
      legend.title = element_blank(),
      legend.text = element_text(size = 18),
      strip.text = element_markdown(size = 18),
      legend.key.width = unit(1.5, "cm")
    ) +
    guides(color = guide_legend(order = 1, nrow = 1))
  
  p.recall <- PRC_sum %>%
    ggplot(aes(x = Ka, y = mean_R, color = Method, linetype = linetype,
               group = interaction(Method, linetype))) +
    geom_line(linewidth = 0.8, position = position_dodge(width = 0.4)) +
    geom_point(size = 1.5, position = position_dodge(width = 0.4)) +
    geom_errorbar(
      aes(ymin = pmax(mean_R - sd_R, 0.5), ymax = pmin(mean_R + sd_R, 1)),
      width = 0.15,
      position = position_dodge(width = 0.4)
    ) +
    facet_grid(~scenario) +
    ylim(0.5, 1) +
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
  
  ## Generate figure result for all signature
  g_c <-  ggpubr::ggarrange(p.precison, p.recall, ncol = 1, common.legend = TRUE,
                            legend = "none", align = "v", axis = "lr", heights = c(1, 1.1))
  
  ggsave(
    filename = "./Simulation/Figure/Fig2_c.png",
    plot = g_c,
    width = 350,
    height = 240,
    units = "mm",
    dpi = 300
  )
  
  ## Signature detection by type for scenario 3 (Fig d) ----
  data.loc.L10 <- "./Simulation/Sim_CRC_stage1/res_Ka"
  
  PRC_all <- NULL
  for(Ka in c(0.1, 0.2, 0.3, 0.4)){
    PRC_tmp <- NULL
    for(s in 1:50){
      if(file.exists(paste0(data.loc.L10, Ka, "_pos0.6_u0_scenario3_s", s, ".Rdata"))){
        load(paste0(data.loc.L10, Ka, "_pos0.6_u0_scenario3_s", s, ".Rdata"))
        PRC_tmp <- rbind(PRC_tmp, result_mat)
      }
    }
    PRC_tmp$L <- paste0("L = 10, ", PRC_tmp$data)
    PRC_tmp$scenario <- "Scenario 3"
    PRC_tmp$Ka <- as.character(Ka)
    PRC_all <- rbind(PRC_all, PRC_tmp)
  }
  PRC_all$Method <- factor(PRC_all$method, levels = c("SMESH", "True-cluster FE", "SKM + FE", "SHC + FE", "Melody"))
  
  PRC_sum <- rbind(
    PRC_all %>%
      group_by(scenario, Ka, Method) %>%
      summarise(
        mean = mean(Precision_hom, na.rm = TRUE),
        sd   = sd(Precision_hom,  na.rm = TRUE),
        setting = "Shared signatures",
        value = "Precision",
        .groups = "drop"
      ),
    PRC_all %>%
      group_by(scenario, Ka, Method) %>%
      summarise(
        mean = mean(Recall_hom, na.rm = TRUE),
        sd   = sd(Recall_hom,  na.rm = TRUE),
        setting = "Shared signatures",
        value = "Recall",
        .groups = "drop"
      ),
    PRC_all %>%
      group_by(scenario, Ka, Method) %>%
      summarise(
        mean = mean(pmin(Precision_het,1), na.rm = TRUE),
        sd   = sd(Precision_het,  na.rm = TRUE),
        setting = "Cluster-specific",
        value = "Precision",
        .groups = "drop"
      ),
    PRC_all %>%
      group_by(scenario, Ka, Method) %>%
      summarise(
        mean = mean(pmin(Recall_het, 1), na.rm = TRUE),
        sd   = sd(Recall_het,  na.rm = TRUE),
        setting = "Cluster-specific",
        value = "Recall",
        .groups = "drop"
      )
  ) %>% dplyr::mutate(setting = factor(setting, levels = c("Shared signatures", "Cluster-specific")),
                      value   = factor(value, levels = c("Precision", "Recall")))
  
  ## Generate figure result for all signature
  p.1 <- PRC_sum %>% filter(setting == "Shared signatures") %>%
    ggplot(aes(x = Ka, y = mean, color = Method, group = Method )) +
    geom_line(linewidth = 0.8, position = position_dodge(width = 0.4)) +
    geom_point(size = 1.5, position = position_dodge(width = 0.4)) +
    geom_errorbar(
      aes(ymin = pmax(mean - sd, 0), ymax = pmin(mean + sd, 1)),
      width = 0.15,
      position = position_dodge(width = 0.4)
    ) +
    facet_grid(~value) +
    ylim(0, 1) +
    xlab("Proportion of true signatures") +
    ylab("") +
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
      axis.text = element_text(size = 18),
      axis.line = element_line(colour = "black"),
      
      panel.grid.major = element_line(colour = "white", linetype = "dotted"),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(colour = "black", fill = NA),
      panel.background = element_rect(fill = "white"),
      
      legend.position = "none",
      legend.box = "vertical",
      legend.title = element_blank(),
      legend.text = element_text(size = 18),
      strip.text = element_markdown(size = 18),
      legend.key.width = unit(1.5, "cm")
    ) +
    guides(color = guide_legend(order = 1, nrow = 1))
  
  p.2 <- PRC_sum %>% filter(setting == "Cluster-specific") %>%
    ggplot(aes(x = Ka, y = mean, color = Method, group = Method )) +
    geom_line(linewidth = 0.8, position = position_dodge(width = 0.4)) +
    geom_point(size = 1.5, position = position_dodge(width = 0.4)) +
    geom_errorbar(
      aes(ymin = pmax(mean - sd, 0), ymax = pmin(mean + sd, 1)),
      width = 0.15,
      position = position_dodge(width = 0.4)
    ) +
    facet_grid(~value) +
    ylim(0, 1) +
    xlab("Proportion of true signatures") +
    ylab("") +
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
      axis.text = element_text(size = 18),
      axis.line = element_line(colour = "black"),
      
      panel.grid.major = element_line(colour = "white", linetype = "dotted"),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(colour = "black", fill = NA),
      panel.background = element_rect(fill = "white"),
      
      legend.position = "none",
      legend.box = "vertical",
      legend.title = element_blank(),
      legend.text = element_text(size = 18),
      strip.text = element_markdown(size = 18),
      legend.key.width = unit(1.5, "cm")
    ) +
    guides(color = guide_legend(order = 1, nrow = 1))
  
  g_d <- ggpubr::ggarrange(p.1, p.2, ncol = 2, align = "v", common.legend = TRUE, legend = "none")
  
  ggsave(
    filename = "./Simulation/Figure/Fig2_d.png",
    plot = g_d,
    width = 350,
    height = 90,
    units = "mm",
    dpi = 300
  )
  