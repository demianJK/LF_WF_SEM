## Figures of Results

# load packages
library(ggplot2)
library(dplyr) 
library(patchwork)
library(cowplot)

# set paths and load data
pathO <- "objects/"
pathG <- "graphics/"
dat <- readRDS(paste0(pathO, "dat.rds"))

# data prep
dat$ICC <- factor(dat$ICC, ordered=TRUE, levels=unique(dat$ICC))
dat$n <- paste0("n = ", dat$n)
dat$n <- factor(dat$n, ordered=TRUE, levels=unique(dat$n))
dat$p <- paste0("p = ", dat$p)
dat$p <- factor(dat$p, ordered=TRUE, levels=unique(dat$p))

which <- filter(dat, approach == "WF", conv == 100)$cond
datWFconv <- filter(dat, cond %in% which)

# cols:rows
# LF-B: p / g
# WF-T: (p*n) / g
limCR <- c(0, 1) # truncated
breaksCR <- sort(unique(dat$ratio[dat$ratio <= 1]))
labelsCR <- c(0.004, rep("", 7), 0.1, rep("", 3), 0.25, "", 0.5, "", 0.8, 1)

# rows
limR <- c(0, max(dat$g))
breaksR <- sort(unique(dat$g))
labelsR <- c(2, "", "", "", "", "", "", 50, 100, 200, 500)

# note that there is a bug in ggplot when scale_color_discrete() is used, the color changes
# sometimes we need the function for renaming the legend etc, thus we use the function throughout, even if we do not print the legend so that the colors match in all figures


## The Effect of Data Format

# Fig. 4 Convergence Aggregated By Sample Size at Level-2

ggplot(dat, aes(x=g, y=conv, col=ICC)) + 
  stat_summary(alpha=0.5) +
  facet_grid(cols = vars(approach)) +
  scale_y_continuous(name="Convergence Rate (%)", expand=c(0.01,0.01)) +
  scale_x_continuous(limits=limR, breaks=breaksR, labels=labelsR) +
  theme_minimal() + 
  scale_color_manual(values=c("#6b3772", "#289491", "#dac72c"), labels=c("0.05", "0.25", "0.50")) +
  guides(col = guide_legend(override.aes=list(linetype = c(0, 0, 0)))) +
  coord_cartesian(ylim = c(0, 100)) +
  theme(text = element_text(family="serif"), panel.grid.minor = element_blank(), 
        axis.text.x = element_text(angle=45, margin=margin(10,0,0,0)),
        legend.pos = "bottom", panel.border = element_rect(color = "grey",  fill = NA, linewidth = 0.5))

ggsave(
  filename = "MS_fig4.jpeg",
  path = pathG,
  device = "jpeg",
  plot = last_plot(),
  width = 500*3, 
  height = 300*3,
  units = "px",
  dpi = 300
)


# Fig. 5 Estimation Accuracy of Between-Group Parameters Aggregated By Sample Size at Level-2

relRMSE_B_R <- 
  ggplot(dat, aes(x=g, y=relRMSE_B, col=ICC)) + stat_summary(alpha=0.5) +
  facet_grid(cols = vars(approach)) + 
  ylab("Relative RMSE (%)") +
  scale_x_continuous(limits=limR, breaks=breaksR, labels=NULL) +
  theme_minimal() + 
  scale_color_manual(values=c("#6b3772", "#289491", "#dac72c"), labels=c("0.05", "0.25", "0.50")) +
  guides(col = guide_legend(override.aes=list(linetype = c(0, 0, 0)))) +
  theme(text = element_text(family="serif"), panel.grid.minor = element_blank(), 
        axis.title.x= element_blank(),
        legend.pos = "none", panel.border = element_rect(color = "grey",  fill = NA, linewidth = 0.5)) +
  labs(title="A")

relBias_B_R <- 
  ggplot(dat, aes(x=g, y=relBias_B, col=ICC)) + stat_summary(alpha=0.5) +
  facet_grid(cols = vars(approach)) + 
  ylab("Relative Bias (%)") +
  scale_x_continuous(limits=limR, breaks=breaksR, labels=NULL) +
  theme_minimal() + 
  scale_color_manual(values=c("#6b3772", "#289491", "#dac72c"), labels=c("0.05", "0.25", "0.50")) +
  guides(col = guide_legend(override.aes=list(linetype = c(0, 0, 0)))) +
  theme(text = element_text(family="serif"), panel.grid.minor = element_blank(), 
        axis.title.x= element_blank(), strip.text.x = element_blank(),
        legend.pos = "none", panel.border = element_rect(color = "grey",  fill = NA, linewidth = 0.5)) +
  labs(title="B")

relVar_B_R <- 
  ggplot(dat, aes(x=g, y=relVar_B, col=ICC)) + stat_summary(alpha=0.5) +
  facet_grid(cols = vars(approach)) + 
  ylab("Relative Variance (%)") +
  scale_x_continuous(limits=limR, breaks=breaksR, labels=labelsR) +
  theme_minimal() + 
  scale_color_manual(values=c("#6b3772", "#289491", "#dac72c"), labels=c("0.05", "0.25", "0.50")) +
  guides(col = guide_legend(override.aes=list(linetype = c(0, 0, 0)))) +
  theme(text = element_text(family="serif"), panel.grid.minor = element_blank(), 
        strip.text.x = element_blank(),
        axis.text.x = element_text(angle=45, margin=margin(10,0,0,0)),
        legend.pos = "bottom", panel.border = element_rect(color = "grey",  fill = NA, linewidth = 0.5)) +
  labs(title="C")

relRMSE_B_R + relBias_B_R + relVar_B_R + plot_layout(ncol=1)

ggsave(
  filename = "MS_fig5.jpeg",
  path = pathG,
  device = "jpeg",
  plot = last_plot(),
  width = 500*3, 
  height = 700*3,
  units = "px",
  dpi = 300
)


# Fig. 6 Estimation Accuracy of Within-Group Parameters Aggregated By Sample Size at Level-2

relRMSE_W_R <- 
  ggplot(dat, aes(x=g, y=relRMSE_W, col=ICC)) + stat_summary(alpha=0.5) +
  facet_grid(cols = vars(approach)) + 
  ylab("Relative RMSE (%)") +
  scale_y_continuous(limits=c(0,80), breaks=seq(0,75, 25), expand=c(0,0)) +
  scale_x_continuous(limits=limR, breaks=breaksR, labels=NULL) +
  theme_minimal() + 
  scale_color_manual(values=c("#6b3772", "#289491", "#dac72c"), labels=c("0.05", "0.25", "0.50")) +
  guides(col = guide_legend(override.aes=list(linetype = c(0, 0, 0)))) +
  theme(text = element_text(family="serif"), panel.grid.minor = element_blank(), 
        axis.title.x= element_blank(),
        legend.pos = "none", panel.border = element_rect(color = "grey",  fill = NA, linewidth = 0.5)) +
  labs(title="A")

relBias_W_R <- 
  ggplot(dat, aes(x=g, y=relBias_W, col=ICC)) + stat_summary(alpha=0.5) +
  facet_grid(cols = vars(approach)) + 
  ylab("Relative Bias (%)") +
  scale_y_continuous(limits=c(-2,4), breaks=seq(-2,4, 1), expand=c(0,0)) +
  scale_x_continuous(limits=limR, breaks=breaksR, labels=NULL) +
  theme_minimal() + 
  scale_color_manual(values=c("#6b3772", "#289491", "#dac72c"), labels=c("0.05", "0.25", "0.50")) +
  guides(col = guide_legend(override.aes=list(linetype = c(0, 0, 0)))) +
  theme(text = element_text(family="serif"), panel.grid.minor = element_blank(), 
        axis.title.x= element_blank(), strip.text.x = element_blank(),
        legend.pos = "none", panel.border = element_rect(color = "grey",  fill = NA, linewidth = 0.5)) +
  labs(title="B")

relVar_W_R <- 
  ggplot(dat, aes(x=g, y=relVar_W, col=ICC)) + stat_summary(alpha=0.5) +
  facet_grid(cols = vars(approach)) + 
  ylab("Relative Variance (%)") +
  scale_y_continuous(limits=c(0,30), breaks=seq(0,30, 5), expand=c(0,0)) +
  scale_x_continuous(limits=limR, breaks=breaksR, labels=labelsR) +
  theme_minimal() + 
  scale_color_manual(values=c("#6b3772", "#289491", "#dac72c"), labels=c("0.05", "0.25", "0.50")) +
  guides(col = guide_legend(override.aes=list(linetype = c(0, 0, 0)))) +
  theme(text = element_text(family="serif"), panel.grid.minor = element_blank(), 
        strip.text.x = element_blank(),
        axis.text.x = element_text(angle=45, margin=margin(10,0,0,0)),
        legend.pos = "bottom", panel.border = element_rect(color = "grey",  fill = NA, linewidth = 0.5)) +
  labs(title="C")

relRMSE_W_R + relBias_W_R + relVar_W_R + plot_layout(ncol=1)

ggsave(
  filename = "MS_fig6.jpeg",
  path = pathG,
  device = "jpeg",
  plot = last_plot(),
  width = 500*3,
  height = 700*3,
  units = "px",
  dpi = 300
)



## The Effect of Cols:Rows in Each Data Format

# Fig. 7 Convergence Aggregated by Cols:Rows

ggplot(dat, aes(x=ratio, y=conv, col=ICC)) + stat_summary(alpha=0.5) +
  facet_grid(cols = vars(approach), labeller=label_parsed) +
  scale_y_continuous(name="Convergence Rate (%)", expand=c(0.01,0.01)) +
  scale_x_continuous(limits=limCR, breaks=breaksCR, labels=labelsCR, name=expression(cols:rows)) +
  theme_minimal() +
  scale_color_manual(values=c("#6b3772", "#289491", "#dac72c"), labels=c("0.05", "0.25", "0.50")) +
  guides(col = guide_legend(override.aes=list(linetype = c(0, 0, 0)))) +
  theme(text = element_text(family="serif"), panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle=45, margin=margin(10,0,0,0)),
        legend.pos = "bottom", panel.border = element_rect(color = "grey",  fill = NA, linewidth = 0.5))

ggsave(
  filename = "MS_fig7.jpeg",
  path = pathG,
  device = "jpeg",
  plot = last_plot(),
  width = 500*3,
  height = 300*3,
  units = "px",
  dpi = 300
)

# Fig. 8 Estimation Accuracy of Between-Group Parameters Aggregated by Cols:Rows

relRMSE_B_CR <- 
  ggplot(datWFconv, aes(x=ratio, y=relRMSE_B, col=ICC)) + stat_summary(alpha=0.5) +
  facet_grid(cols = vars(approach), labeller=label_parsed) + 
  ylab("Relative RMSE (%)") +
  scale_x_continuous(limits=limCR, breaks=breaksCR, labels=NULL) +
  theme_minimal() + 
  theme(text = element_text(family="serif"), panel.grid.minor = element_blank(), 
        axis.title.x= element_blank(),
        legend.pos = "none", panel.border = element_rect(color = "grey",  fill = NA, linewidth = 0.5)) +
  labs(title="A")

relBias_B_CR <- 
  ggplot(datWFconv, aes(x=ratio, y=relBias_B, col=ICC)) + stat_summary(alpha=0.5) +
  facet_grid(cols = vars(approach), labeller=label_parsed) + 
  ylab("Relative Bias (%)") +
  scale_x_continuous(limits=limCR, breaks=breaksCR, labels=NULL) +
  theme_minimal() + 
  theme(text = element_text(family="serif"), panel.grid.minor = element_blank(), 
        axis.title.x= element_blank(), strip.text.x = element_blank(),
        legend.pos = "none", panel.border = element_rect(color = "grey",  fill = NA, linewidth = 0.5)) +
  labs(title="B")

relVar_B_CR <- 
  ggplot(datWFconv, aes(x=ratio, y=relVar_B, col=ICC)) + stat_summary(alpha=0.5) +
  facet_grid(cols = vars(approach), labeller=label_parsed) + 
  ylab("Relative Variance (%)") +
  scale_x_continuous(limits=limCR, breaks=breaksCR, labels=labelsCR, name=expression(cols:rows)) +
  theme_minimal() + 
  theme(text = element_text(family="serif"), panel.grid.minor = element_blank(), 
        strip.text.x = element_blank(), axis.text.x = element_text(angle=45, margin=margin(10,0,0,0)),
        legend.pos = "bottom", panel.border = element_rect(color = "grey",  fill = NA, linewidth = 0.5)) +
  labs(title="C")

relRMSE_B_CR + relBias_B_CR + relVar_B_CR + plot_layout(ncol=1)

ggsave(
  filename = "MS_fig8.jpeg",
  path = pathG,
  device = "jpeg",
  plot = last_plot(),
  width = 500*3, 
  height = 700*3,
  units = "px",
  dpi = 300
)


# Fig. 9 Estimation Accuracy ofWithin-Group Parameters Aggregated by Cols:Rows

relRMSE_W_CR <- 
  ggplot(datWFconv, aes(x=ratio, y=relRMSE_W, col=ICC)) + stat_summary(alpha=0.5) +
  facet_grid(cols = vars(approach)) + 
  ylab("Relative RMSE (%)") +
  scale_y_continuous(limits=c(0,80), breaks=seq(0,75, 25), expand=c(0,0)) +
  scale_x_continuous(limits=limCR, breaks=breaksCR, labels=NULL) +
  theme_minimal() + 
  theme(text = element_text(family="serif"), panel.grid.minor = element_blank(), 
        axis.title.x= element_blank(),
        legend.pos = "none", panel.border = element_rect(color = "grey",  fill = NA, linewidth = 0.5)) +
  labs(title="A")

relBias_W_CR <- 
  ggplot(datWFconv, aes(x=ratio, y=relBias_W, col=ICC)) + stat_summary(alpha=0.5) +
  facet_grid(cols = vars(approach)) + 
  ylab("Relative Bias (%)") +
  scale_y_continuous(limits=c(-2,4), breaks=seq(-2,4, 1), expand=c(0,0)) +
  scale_x_continuous(limits=limCR, breaks=breaksCR, labels=NULL) +
  theme_minimal() + 
  theme(text = element_text(family="serif"), panel.grid.minor = element_blank(), 
        axis.title.x= element_blank(), strip.text.x = element_blank(),
        legend.pos = "none", panel.border = element_rect(color = "grey",  fill = NA, linewidth = 0.5)) +
  labs(title="B")

relVar_W_CR <- 
  ggplot(datWFconv, aes(x=ratio, y=relVar_W, col=ICC)) + stat_summary(alpha=0.5) +
  facet_grid(cols = vars(approach)) + 
  ylab("Relative Variance (%)") +
  scale_y_continuous(limits=c(0,30), breaks=seq(0,30, 5), expand=c(0,0)) +
  scale_x_continuous(limits=limCR, breaks=breaksCR, labels=labelsCR, name=expression(cols:rows)) +
  theme_minimal() + 
  theme(text = element_text(family="serif"), panel.grid.minor = element_blank(), 
        strip.text.x = element_blank(),
        axis.text.x = element_text(angle=45, margin=margin(10,0,0,0)),
        legend.pos = "bottom", panel.border = element_rect(color = "grey",  fill = NA, linewidth = 0.5)) +
  labs(title="C")

relRMSE_W_CR + relBias_W_CR + relVar_W_CR + plot_layout(ncol=1)

ggsave(
  filename = "MS_fig9.jpeg",
  path = pathG,
  device = "jpeg",
  plot = last_plot(),
  width = 500*3, 
  height = 700*3,
  units = "px",
  dpi = 300
)


# Online Supplemental 

# Fig. 1 Convergence and Estimation Accuracy Aggregated by All Simulation Factors

# kappa removed because information not well presented anymore in figure

# How many conditions with only infinite kappas?
infK_LF <- table(is.na(dat$kappa_B[dat$approach == "LF"])) #
round( (infK_LF[2]/(infK_LF[1]+infK_LF[2])*100), 2 ) # 3.75%
infK_WF <- table(is.na(dat$kappa_B[dat$approach == "WF"])) #
round( (infK_WF[2]/(infK_WF[1]+infK_WF[2])*100), 2 ) # NA = 0%

# Mean over all conditions of infinite kappas?
round(mean(dat$Infkappa_B[dat$approach == "LF"]), 2) # 32.65%
round(mean(dat$Infkappa_B[dat$approach == "WF"]), 2) # 0.03%

ns <- ggplot(dat) + geom_point(aes(x=g, y=ns_B, col=approach, shape=ICC), show.legend = FALSE, alpha=0.5) +
  facet_grid(rows = vars(p), cols = vars(n), drop=TRUE, scales="free_x") +
  ylab("Non-Singular (%)") +
  scale_x_continuous(name=NULL, labels=NULL) +
  scale_color_discrete(name="(Sample) Covariance Matrix", labels=c(expression(hat(Sigma)[LF-B]), expression(hat(Sigma)[WF-T]))) +
  theme_minimal() +
  theme(text = element_text(family="serif"), panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "grey",  fill = NA, linewidth = 0.5)) +
  guides(col = guide_legend(order = 1), shape = guide_legend(order = 2)) +
  labs(title="A")

pd <- ggplot(dat) + geom_point(aes(x=g, y=pd_B, col=approach, shape=ICC), alpha=0.5) +
  facet_grid(rows = vars(p), cols = vars(n), drop=TRUE, scales="free_x") +
  ylab("Positive Definite (%)") +
  theme_minimal() +
  scale_color_discrete(name="(Sample) Covariance Matrix", labels=c(expression(hat(Sigma)[LF-B]), expression(hat(Sigma)[WF-T]))) +
  scale_shape_discrete(labels=c("0.05", "0.25", "0.50")) +
  theme(text = element_text(family="serif"), panel.grid.minor = element_blank(), axis.text.x = element_text(angle=45),
        strip.text.x = element_blank(),
        panel.border = element_rect(color = "grey",  fill = NA, linewidth = 0.5)) +
  guides(col = guide_legend(order = 1), shape = guide_legend(order = 2)) +
  labs(title="B")

ns + pd + plot_layout(ncol = 1, guides='collect') & theme(text = element_text(family="serif"), legend.position = "bottom")

ggsave(
  filename = "OS_fig1.jpeg",
  path = pathG,
  device = "jpeg",
  plot = last_plot(),
  width = 600*3,
  height = 750*3,
  units = "px",
  dpi = 300
)


# Fig. 2 Convergence and Estimation Accuracy Aggregated by All Simulation Factors

conv <- ggplot(dat) + geom_point(aes(x=g, y=conv, col=approach, shape=ICC), alpha=0.5) +
  facet_grid(rows = vars(p), cols = vars(n), drop=TRUE, scales="free_x") +
  scale_y_continuous(name="Convergence Rate (%)", expand=c(0.01,0.01)) +
  scale_x_continuous(name=NULL, labels=NULL) +
  scale_color_discrete(name="Data Format") +
  theme_minimal() +
  theme(text = element_text(family="serif"), panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "grey",  fill = NA, linewidth = 0.5)) +
  labs(title="A")

relRMSE_B <-
  ggplot(dat) + geom_point(aes(x=g, y=relRMSE_B, col=approach, shape=ICC), show.legend = FALSE, alpha=0.5) +
  facet_grid(rows = vars(p), cols = vars(n), drop=TRUE, scales="free_x") +
  ylab("Relative RMSE (%)") +
  scale_x_continuous(name=NULL, labels=NULL) +
  scale_color_discrete(name="Data Format") +
  geom_text(data = subset(subset(dat, n == "n = 100"), p == "p = 2"), aes(x = Inf, y= Inf, label = "hat(theta)[between]"), col=c(rep("white", 11), "black"), parse=TRUE, hjust=1, vjust=1) +
  theme_minimal() +
  theme(text = element_text(family="serif"), panel.grid.minor = element_blank(),
        strip.text.x = element_blank(), strip.text.y = element_blank(),
        panel.border = element_rect(color = "grey",  fill = NA, linewidth = 0.5)) +
  labs(title="B")

relRMSE_W <- ggplot(dat) + geom_point(aes(x=g, y=relRMSE_W, col=approach, shape=ICC), show.legend = FALSE, alpha=0.5) +
  facet_grid(rows = vars(p), cols = vars(n), drop=TRUE, scales="free_x") +
  ylab("Relative RMSE (%)") +
  ylim(c(0, 80)) +
  scale_x_continuous(name=NULL, labels=NULL) +
  scale_color_discrete(name="Data Format") +
  geom_text(data = subset(subset(dat, n == "n = 100"), p == "p = 2"), aes(x = Inf, y= Inf, label = "hat(theta)[within]"), col=c(rep("white", 11), "black"), parse=TRUE, hjust=1, vjust=1) +
  theme_minimal() +
  theme(text = element_text(family="serif"), panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "grey",  fill = NA, linewidth = 0.5)) +
  labs(title="C")

relBias_B <- ggplot(dat) + geom_point(aes(x=g, y=relBias_B, col=approach, shape=ICC), show.legend = FALSE, alpha=0.5) +
  facet_grid(rows = vars(p), cols = vars(n), drop=TRUE, scales="free_x") +
  ylab("Relative Bias (%)") +
  scale_x_continuous(name=NULL, labels=NULL) +
  scale_color_discrete(name="Data Format") +
  geom_text(data = subset(subset(dat, n == "n = 100"), p == "p = 2"), aes(x = Inf, y= Inf, label = "hat(theta)[between]"), col=c(rep("white", 11), "black"), parse=TRUE, hjust=1, vjust=1) +
  theme_minimal() +
  theme(text = element_text(family="serif"), panel.grid.minor = element_blank(), axis.text.x = element_text(angle=45),
        strip.text.x = element_blank(), strip.text.y = element_blank(),
        panel.border = element_rect(color = "grey",  fill = NA, linewidth = 0.5)) +
  labs(title="D")

relBias_W <- ggplot(dat) + geom_point(aes(x=g, y=relBias_W, col=approach, shape=ICC), show.legend = FALSE, alpha=0.5) +
  facet_grid(rows = vars(p), cols = vars(n), drop=TRUE, scales="free_x") +
  ylab("Relative Bias (%)") +
  scale_x_continuous(name=NULL, labels=NULL) +
  scale_color_discrete(name="Data Format") +
  geom_text(data = subset(subset(dat, n == "n = 100"), p == "p = 2"), aes(x = Inf, y= Inf, label = "hat(theta)[within]"), col=c(rep("white", 11), "black"), parse=TRUE, hjust=1, vjust=1) +
  theme_minimal() +
  theme(text = element_text(family="serif"), panel.grid.minor = element_blank(),
        strip.text.x = element_blank(),
        panel.border = element_rect(color = "grey",  fill = NA, linewidth = 0.5)) +
  labs(title="E")

relVar_B <- ggplot(dat) + geom_point(aes(x=g, y=relVar_B, col=approach, shape=ICC), alpha=0.5) +
  facet_grid(rows = vars(p), cols = vars(n), drop=TRUE, scales="free_x") +
  ylab("Relative Variance (%)") +
  geom_text(data = subset(subset(dat, n == "n = 100"), p == "p = 2"), aes(x = Inf, y= Inf, label = "hat(theta)[between]"), col=c(rep("white", 11), "black"), parse=TRUE, hjust=1, vjust=1) +
  scale_color_discrete(name="Data Format") +
  scale_shape_discrete(labels=c("0.05", "0.25", "0.50")) +
  theme_minimal() +
  theme(text = element_text(family="serif"), panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle=45), legend.position = "bottom", legend.direction = "vertical", legend.box = "horizontal",
        strip.text.x = element_blank(), strip.text.y = element_blank(),
        panel.border = element_rect(color = "grey",  fill = NA, linewidth = 0.5)) +
  guides(col = guide_legend(order = 1), shape = guide_legend(order = 2))+
  labs(title="F")

legAll <- get_legend(relVar_B)

relVar_W <- ggplot(dat) + geom_point(aes(x=g, y=relVar_W, col=approach, shape=ICC), show.legend = FALSE, alpha=0.5) +
  facet_grid(rows = vars(p), cols = vars(n), drop=TRUE, scales="free_x") +
  ylab("Relative Variance (%)") +
  scale_color_discrete(name="Data Format") +
  geom_text(data = subset(subset(dat, n == "n = 100"), p == "p = 2"), aes(x = Inf, y= Inf, label = "hat(theta)[within]"), col=c(rep("white", 11), "black"), parse=TRUE, hjust=1, vjust=1) +
  theme_minimal() +
  theme(text = element_text(family="serif"), panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle=45),
        strip.text.x = element_blank(),
        panel.border = element_rect(color = "grey",  fill = NA, linewidth = 0.5))+
  labs(title="G")

conv + legAll + relRMSE_B + relRMSE_W + relBias_B + relBias_W + relVar_B + relVar_W + plot_layout(ncol=2, guides='collect') & theme(legend.position = "none")

ggsave(
  filename = "OS_fig2.jpeg",
  path = pathG,
  device = "jpeg",
  plot = last_plot(),
  width = 800*3,
  height = 1000*3,
  units = "px",
  dpi = 300
)


# Fig. 3 Matrix Properties Aggregated by Sample Size at Level-2

ns_R <- 
  ggplot(dat, aes(x=g, y=ns_B, col=ICC)) + stat_summary(alpha=0.5) +
  facet_grid(cols = vars(approach), labeller=label_parsed) + 
  ylab("Non-Singular (%)") +
  scale_y_continuous(limits=c(0,100), breaks=seq(0,100, 25), expand=c(0.01,0.01)) +
  scale_x_continuous(limits=limR, breaks=breaksR, labels=NULL) +
  theme_minimal() + 
  scale_color_manual(values=c("#6b3772", "#289491", "#dac72c"), labels=c("0.05", "0.25", "0.50")) +
  guides(col = guide_legend(override.aes=list(linetype = c(0, 0, 0)))) +
  theme(text = element_text(family="serif"), panel.grid.minor = element_blank(), 
        axis.title.x= element_blank(),
        legend.pos = "none", panel.border = element_rect(color = "grey",  fill = NA, linewidth = 0.5)) +
  labs(title="A")

pd_R <- 
  ggplot(dat, aes(x=g, y=pd_B, col=ICC)) + stat_summary(alpha=0.5) +
  facet_grid(cols = vars(approach), labeller=label_parsed) + 
  ylab("Positive Definite (%)") +
  scale_y_continuous(limits=c(0,100), breaks=seq(0,100, 25), expand=c(0.01,0.01)) +
  scale_x_continuous(limits=limR, breaks=breaksR, labels=labelsR) +
  theme_minimal() + 
  scale_color_manual(values=c("#6b3772", "#289491", "#dac72c"), labels=c("0.05", "0.25", "0.50")) +
  guides(col = guide_legend(override.aes=list(linetype = c(0, 0, 0)))) +
  theme(text = element_text(family="serif"), panel.grid.minor = element_blank(), 
        strip.text.x = element_blank(), axis.text.x = element_text(angle=45, margin=margin(10,0,0,0)),
        legend.pos = "none", panel.border = element_rect(color = "grey",  fill = NA, linewidth = 0.5)) +
  labs(title="B")


ns_R + pd_R + plot_layout(ncol = 1, guides='collect') & theme(legend.position = "bottom")


ggsave(
  filename = "OS_fig3.jpeg",
  path = pathG,
  device = "jpeg",
  plot = last_plot(),
  width = 500*3, 
  height = 500*3,
  units = "px",
  dpi = 300
)


# Fig. 4 Matrix Properties Aggregated by Cols:Rows

# How many conditions with only infinite kappas?
infK_LF <- table(is.na(dat$kappa_B[dat$ratio <= 1 & dat$approach == "LF"])) # 0%
infK_WF <- table(is.na(dat$kappa_B[dat$ratio <= 1 & dat$approach == "WF"]))  
round( (infK_WF[2]/(infK_WF[1]+infK_WF[2])*100), 2 ) # 0%

# Mean over all conditions of infinite kappas?
round(mean(dat$Infkappa_B[dat$ratio <= 1 & dat$approach == "LF"]), 2) # 22.98%
round(mean(dat$Infkappa_B[dat$ratio <= 1 & dat$approach == "WF"]), 2) # 0%

ns_CR <- 
  ggplot(dat, aes(x=ratio, y=ns_B, col=ICC)) + 
  stat_summary(alpha=0.5) + # If no aggregation functions are supplied, will default to mean_se().
  facet_grid(cols = vars(approach), labeller=label_parsed) + 
  ylab("Non-Singular (%)") +
  scale_y_continuous(limits=c(0,100), breaks=seq(0,100, 25), expand=c(0.01,0.01)) +
  scale_x_continuous(limits=limCR, breaks=breaksCR, labels=NULL) +
  scale_color_manual(values=c("#6b3772", "#289491", "#dac72c"), labels=c("0.05", "0.25", "0.50")) +
  guides(col = guide_legend(override.aes=list(linetype = c(0, 0, 0)))) +
  theme_minimal() + 
  theme(text = element_text(family="serif"), panel.grid.minor = element_blank(), 
        axis.title.x= element_blank(),
        legend.pos = "none", panel.border = element_rect(color = "grey",  fill = NA, linewidth = 0.5)) +
  labs(title="A")

pd_CR <- 
  ggplot(dat, aes(x=ratio, y=pd_B, col=ICC)) + 
  stat_summary(alpha=0.5) +
  facet_grid(cols = vars(approach), labeller=label_parsed) + 
  ylab("Positive Definite (%)") +
  scale_y_continuous(limits=c(0,100), breaks=seq(0,100, 25), expand=c(0.01,0.01)) +
  scale_x_continuous(limits=limCR, breaks=breaksCR, labels=labelsCR, name=expression(cols:rows)) +
  scale_color_manual(values=c("#6b3772", "#289491", "#dac72c"), labels=c("0.05", "0.25", "0.50")) +
  guides(col = guide_legend(override.aes=list(linetype = c(0, 0, 0)))) +
  theme_minimal() + 
  theme(text = element_text(family="serif"), panel.grid.minor = element_blank(), 
        strip.text.x = element_blank(), axis.text.x = element_text(angle=45, margin=margin(10,0,0,0)),
        legend.pos = "none", panel.border = element_rect(color = "grey",  fill = NA, linewidth = 0.5)) +
  labs(title="B")

ns_CR + pd_CR + plot_layout(ncol = 1, guides='collect') & theme(legend.position = "bottom")


ggsave(
  filename = "OS_fig4.jpeg",
  path = pathG,
  device = "jpeg",
  plot = last_plot(),
  width = 500*3, 
  height = 500*3,
  units = "px",
  dpi = 300
)
