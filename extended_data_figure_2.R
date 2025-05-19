#code to visualize extended data figure 2

library(ggplot2)
library(dplyr)
library(patchwork)
df <- read.csv("path_to_folder_with_metabolic_cage_data/tnf_metabolic_cage.csv")
df$genotype <- factor(df$genotype, 
                      levels = c("WT","TNF", "uMTTNF"), 
                      labels = c("WT", "TNF\u0394ARE/+", "\u03bcMT-TNF\u0394ARE/+"))

p1 <- ggplot(df, aes(total_mass, average_heat_light_unnorm, color = genotype, shape = genotype)) + 
      geom_point(size = 1, stroke = 1) + 
      geom_smooth(method = lm) + 
      labs(title = "Energy Expenditure\nLight") +
      xlab("Body weight (g)") +
      ylab("Heat (kcal/hr)") +
      coord_cartesian(xlim = c(16, 23), ylim = c(0.25, 0.55)) +
      scale_x_continuous(expand = c(0, 0), n.breaks = 8) +
      scale_y_continuous(expand = c(0, 0), n.breaks = 9) + 
      theme(plot.title = element_text(hjust = 0.5)) +
      scale_color_manual(values = c("#4DAF4A", "#377EB8", "#E41A1C")) +
      scale_shape_manual(values = c(1, 15, 17)) +
      theme(text = element_text(size = 12),
            panel.background = element_blank(),
            panel.border = element_rect(fill = NA, color = "black"),
            panel.grid = element_line(color = "gray"),
            axis.line = element_line(color = "black"),
            legend.position = "none")
p2 <- ggplot(df, aes(total_mass, average_heat_dark_unnorm, color = genotype, shape = genotype)) + 
      geom_point(size = 1, stroke = 1) + 
      geom_smooth(method = lm) + 
      labs(title = "Energy Expenditure\nDark") +
      xlab("Body weight (g)") +
      coord_cartesian(xlim = c(16, 23), ylim = c(0.25, 0.55)) +
      scale_x_continuous(expand = c(0, 0), n.breaks = 8) +
      scale_y_continuous(expand = c(0, 0), n.breaks = 9) + 
      theme(plot.title = element_text(hjust = 0.5)) +
      scale_color_manual(values = c("#4DAF4A", "#377EB8", "#E41A1C")) +
      scale_shape_manual(values = c(1, 15, 17)) +
      theme(text = element_text(size = 12),
            panel.background = element_blank(),
            panel.border = element_rect(fill = NA, color = "black"),
            panel.grid = element_line(color = "gray"),
            axis.line = element_line(color = "black"),
            axis.title.y = element_blank(),
            legend.position = "none")
p3 <- ggplot(df, aes(total_mass, food_light, color = genotype, shape = genotype)) + 
      geom_point(size = 1, stroke = 1) + 
      geom_smooth(method = lm) + 
      labs(title = "Food Intake\nLight") +
      xlab("Body weight (g)") +
      ylab("Food intake (g)") +
      coord_cartesian(xlim = c(16, 23), ylim = c(0, 5)) +
      scale_x_continuous(expand = c(0, 0), n.breaks = 8) +
      scale_y_continuous(expand = c(0, 0), n.breaks = 7) + 
      theme(plot.title = element_text(hjust = 0.5)) +
      scale_color_manual(values = c("#4DAF4A", "#377EB8", "#E41A1C")) +
      scale_shape_manual(values = c(1, 15, 17)) +
      theme(text = element_text(size = 12),
            panel.background = element_blank(),
            panel.border = element_rect(fill = NA, color = "black"),
            panel.grid = element_line(color = "gray"),
            axis.line = element_line(color = "black"),
            legend.position = "none")
p4 <- ggplot(df, aes(total_mass, food_dark, color = genotype, shape = genotype)) + 
      geom_point(size = 1, stroke = 1) + 
      geom_smooth(method = lm) + 
      labs(title = "Food Intake\nDark") +
      xlab("Body weight (g)") +
      coord_cartesian(xlim = c(16, 23), ylim = c(0, 5)) +
      scale_x_continuous(expand = c(0, 0), n.breaks = 8) +
      scale_y_continuous(expand = c(0, 0), n.breaks = 7) + 
      theme(plot.title = element_text(hjust = 0.5)) +
      scale_color_manual(values = c("#4DAF4A", "#377EB8", "#E41A1C")) +
      scale_shape_manual(values = c(1, 15, 17)) +
      theme(text = element_text(size = 12),
            panel.background = element_blank(),
            panel.border = element_rect(fill = NA, color = "black"),
            panel.grid = element_line(color = "gray"),
            axis.line = element_line(color = "black"),
            axis.title.y = element_blank(),
            legend.position = "none")

panel_b <- p1 | p2
panel_c <- p3 | p4

ggsave("path_to_output/extended_fig2b_R_graph.pdf", plot = panel_b, height = 3, units = "in")
ggsave("path_to_output/extended_fig2c_R_graph.pdf", plot = panel_c, height = 3, units = "in")
