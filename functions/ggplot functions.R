gg_text <- function(
    size_minor = 8, size_major = 9, size_title = 10, x_axis_rm = F,
    x_angle = 0, ystrip_angle = 270, family = "sans", 
    title_hjust = 0.5, subtitle_hjust = 0.5, caption_hjust = 1) 
{
  theme(text = element_text(family = family, colour = "black", size = size_minor)) +
    
    theme(axis.text.x = element_text(size = size_minor, colour = "black", angle = x_angle)) + 
    theme(axis.text.y = element_text(size = size_minor, colour = "black")) + 
    theme(axis.ticks = element_line(size = 0.5, colour = "black")) + 
    
    theme(axis.title.x = element_text(size = size_major, colour = "black")) + 
    theme(axis.title.y = element_text(size = size_major, colour = "black")) +
    
    theme(strip.text.x = element_text(size = size_major, colour = "black")) +
    theme(strip.text.y = element_text(size = size_major, colour = "black", angle = ystrip_angle)) +
    theme(legend.text = element_text(size = size_minor, colour = "black")) +
    theme(legend.title = element_text(size = size_major, colour = "black",)) +
    
    theme(plot.title = element_text(size = size_title, colour = "black", 
                                    hjust = title_hjust)) +
    theme(plot.subtitle = element_text(size=size_major, colour = "black", 
                                       hjust = subtitle_hjust)) +
    theme(plot.caption = element_text(size = size_minor, colour = "grey45", 
                                      hjust = caption_hjust)) +
    if(x_axis_rm == T){
      theme(axis.title.x=element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank())
    }
  
}


gg_aesthetic <- function(
    yint = 0, xint = NULL, theme = 1, 
    yint_alpha = 1,
    xlab_vec = NULL, alpha_vec = NULL, colour_vec = NULL, linetype_vec = NULL, shape_vec = NULL,  size_vec = NULL, 
    panel.spacing.x = 0.075, panel.spacing.y = 0.075, facet = F)
{
  list(
    geom_hline(yintercept = yint, lty = 5, linewidth = 0.5, alpha = yint_alpha),
    geom_vline(xintercept = xint, lty = 5, linewidth = 0.5),
    theme(panel.spacing.y = unit(panel.spacing.y, "inch"),
          panel.spacing.x = unit(panel.spacing.x, "inch") ),
    {
      if (theme == 1){
        theme_bw()
      } else {
        list(
          theme_minimal(), # gets of the boxes around facet titles.
          annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, size=1), 
          theme(plot.background = element_blank(),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.border = element_blank() )
        )
      }
    },
    {
      if(is.null(xlab_vec)==F){
        scale_x_discrete(labels = xlab_vec) #set x axis labels manually
      }
    },
    
    {
      if(is.null(alpha_vec)==F){
        scale_alpha_manual(values = alpha_vec) #set alpha  values manually
      }
    },
    
    {
      if(is.null(colour_vec)==F){
        scale_colour_manual(values = colour_vec) #set colour scheme manually
      }
    },
    
    {
      if(is.null(linetype_vec)==F){
        scale_linetype_manual(values = linetype_vec) #set linetype scheme manually
      }
    },
    
    {
      if(is.null(shape_vec)==F){
        scale_shape_manual(values = shape_vec) #set shape scheme manually
      }
    },
    
    {
      if(is.null(size_vec)==F){
        scale_size_manual(values = size_vec) #set shape scheme manually
      }
    }
  )
}



gg_qr_edit <- function(xintercept = 5.5, rect_col = "gray", alpha = 0.1) {
  list(
    geom_vline(xintercept = xintercept),
    geom_rect(aes(ymin = -Inf, ymax = Inf, xmin = xintercept, xmax = Inf), 
              fill = rect_col, alpha = 0.1)
  )
  
}

gg_legend <- function(
    preset = 0, position = "bottom", title.align = 0.5, direction = "vertical", 
    box = "vertical", margin= 0.00625, spacing = 0.00625
) 
{
  if (preset == 1){
    position <- "none"
  } 
  if (preset == 2){
    position <- "bottom"
    direction <- "horizontal"
    box <- "vertical"
  }
  if (preset == 3){
    position <- "right"
    direction <- "vertical"
    box <- "vertical"
  }
  theme(
    legend.position = position, # position of legend relative to the rest of the plot
    legend.title.align = title.align, # alignment of legend title
    legend.direction = direction, # arrange of items within a single legend
    legend.box = box, # arrangement of multiple legends
    legend.background = element_rect(size=0.5, linetype="solid"),
    
    # spacing
    legend.box.margin = margin(rep(margin, 4), unit = "inch"),
    legend.spacing = unit(spacing, units = "inch") )
}

gg_coord <- function(
    coord_flip = F,  limits = c(), breaks = waiver()
) {
  list(
    scale_y_continuous(labels = signs_format(), #ensure y axis labels use minus sign, not hyphen
                       breaks = breaks),
    coord_cartesian(ylim = limits),
    {
      if(coord_flip == T){ 
        coord_flip(ylim = limits) # flip coordinate system
      }
    } 
  )
}


gg_geom_point <- function(
    # basic
  .data, x, y=beta, ymin=lb, ymax=ub, p_dodge = T, p_dodge_width = 1,
  shape = NULL, colour = NULL, linetype = NULL, size = NULL,
  
  # Faceting
  facet_type = F, vert = "", hori = "", scales = "fixed", nrow = NULL, ncol = NULL, 
  
  # Aestetic
  theme = 1, qr_edit = F, rect_col = "gray", yint = 0,
  shape_vec = NULL, linetype_vec = NULL, colour_vec = NULL, xlab_vec = NULL, size_vec = NULL, 
  panel.spacing.x = 0.075, panel.spacing.y = 0.075,
  
  # titles
  title = NULL, subtitle = NULL, caption = NULL, 
  ylab = NULL, xlab = waiver(),
  col_title = waiver(), shape_title = waiver(), linetype_title = waiver(), size_title = waiver(),
  
  # Fix font
  text_size_minor = 8, text_size_major = 9, text_size_title = 10,
  text_family = "sans",
  
  # legend
  legend.position = "bottom", legend.title.align = 0.5, legend.direction = "vertical", 
  legend.box = "vertical", legend.margin_inch = 0.00625, legend.spacing = 0.00625,
  
  # misc
  coord_flip = F,  limits = c(), breaks = waiver(), x_angle = 0, ystrip_angle = 270
) {
  
  # if breaks were specified, then limits need to be specified too. 
  
  if(is.vector(breaks) == T){
    lower_limit <- min(breaks)
    upper_limit <- max(breaks)
    limits <- c(lower_limit, upper_limit)
  }
  
  fig <- ggplot() + 
    
    
    
    # Add details if Quantile regression is used
    {
      if (qr_edit == T){
        geom_vline(xintercept = 5.5) 
      }
    } +
    
    {
      if (qr_edit == T){
        geom_rect(data = .data,
                  aes(ymin = -Inf, ymax = Inf, xmin = 5.5, xmax = Inf), 
                  fill = rect_col, alpha = 0.1)
      }
    } +
    
    # map results ----
  {    
    if(p_dodge == T){
      geom_pointrange(data = .data, 
                      # shape = shape_manual, 
                      aes(x = {{x}}, y = {{y}}, ymin = {{ymin}}, ymax = {{ymax}},
                          shape = {{shape}}, color = {{colour}}, 
                          linetype = {{linetype}}, size = {{size}}),
                      position = position_dodge(width = p_dodge_width))
    }
  } +
    
    {
      if(p_dodge == F){
        geom_pointrange(data = .data,
                        aes(x = {{x}}, y = {{y}}, ymin = {{ymin}}, ymax = {{ymax}},
                            shape = {{shape}}, color = {{colour}}, 
                            linetype = {{linetype}}, size = {{size}}))
      }
    } +
    # geom_pointrange(data = data,
    #                       aes(x = {{x}}, y = {{y}}, ymin = {{ymin}}, ymax = {{ymax}},
    #                           shape = {{shape}}, color = {{colour}}, linetype = {{linetype}})) +
    
    
    
    # organize facets----
  {
    if(facet_type == "grid"){
      facet_grid(paste0({{vert}}, '~', {{hori}}), scales = scales)
    } 
  } +
    {
      if(facet_type == "wrap") {
        facet_wrap(paste0({{vert}}, "~", {{hori}}), scales = scales, nrow= nrow, ncol = ncol)
      } 
    } +
    
    # aesthetic ----
  {
    if(is.null(yint) == F){
      geom_hline(yintercept = yint, lty =5, size=0.5) #set y-intercept = 0 (the null value)
    }
  } +
    
    ## theme 1 is black and white
    {
      if (theme == 1){
        theme_bw()
      } 
    } +
    
    ## theme 2 more minimal
    ### add lines to the bottom and right edge of each facet
    {
      if (theme == 2){
        annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, size=1)
      }
    } +
    {
      if (theme == 2){
        annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, size=1) 
      }
    } +
    
    {
      if (theme == 2){
        theme_minimal() + # gets of the boxes around facet titles.   
          theme(plot.background = element_blank(),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.border = element_blank() ) +
          theme(strip.text = element_text(face = "bold")) 
        
      }
    } +
    
    {
      if(is.null(xlab_vec)==F){
        scale_x_discrete(labels = xlab_vec) #set x axis labels manually
      }
    } + 
    
    {
      if(is.null(colour_vec)==F){
        scale_colour_manual(values = colour_vec) #set colour scheme manually
      }
    } +
    
    {
      if(is.null(linetype_vec)==F){
        scale_linetype_manual(values = linetype_vec) #set linetype scheme manually
      }
    } +
    
    {
      if(is.null(shape_vec)==F){
        scale_shape_manual(values = shape_vec) #set shape scheme manually
      }
    } +
    
    {
      if(is.null(size_vec)==F){
        scale_size_manual(values = size_vec) #set shape scheme manually
      }
    } +
    
    ## add option to increase facet margins
    {
      if(facet_type == "grid"  | facet_type == "wrap"){
        theme(panel.spacing.x = unit(panel.spacing.x, "inch")) 
      } 
    } +
    
    {
      if(facet_type == "grid"  | facet_type == "wrap"){
        theme(panel.spacing.y = unit(panel.spacing.y, "inch")) 
      } 
    } +
    
    
    # fix font ----
  theme(text = element_text(family = text_family, colour = "black", size = text_size_minor)) +
    theme(axis.text.x = element_text(size=text_size_minor, colour = "black", angle = x_angle)) + 
    theme(axis.text.y = element_text(size=text_size_minor, colour = "black")) + 
    theme(axis.ticks = element_line(size = 0.5, colour = "black")) + 
    
    theme(axis.title.x = element_text(size=text_size_major, colour = "black")) + 
    theme(axis.title.y = element_text(size=text_size_major, colour = "black")) +
    # theme(axis.title.y = element_text(margin = text_size_major(0,0,0,0, unit="in"))) +
    # theme(axis.title.x = element_text(margin = margin(0,0,0,0, unit="in"))) +
    
    theme(strip.text.x = element_text(size=text_size_major, colour = "black")) +
    theme(strip.text.y = element_text(size=text_size_major, colour = "black", angle = ystrip_angle)) +
    theme(legend.text = element_text(size=text_size_minor, colour = "black")) +
    theme(legend.title = element_text(size=text_size_major, colour = "black")) +
    
    theme(title = element_text(size = text_size_title, colour = "black")) +
    theme(plot.caption = element_text(size = text_size_minor, colour = "grey45", 
                                      hjust = 0)) + # move the caption to the right side 
    
    
    
    # legend ----
  theme(legend.position = legend.position, legend.title.align = legend.title.align,
        legend.background = element_rect(size=0.5, linetype="solid"),
        legend.direction = legend.direction, legend.box = legend.box,
        legend.box.margin = margin(legend.margin_inch, legend.margin_inch,
                                   legend.margin_inch, legend.margin_inch, unit = "inch"),
        legend.spacing = unit(legend.spacing, units = "inch")
  ) +
    
    #labels ----
  labs(title = title, subtitle = subtitle, caption = caption) +
    labs(col = col_title) +
    labs(shape = shape_title) +
    labs(linetype = linetype_title) +
    
    {
      if(is.null(shape_title)==T){
        guides(shape = F)
      }
    } +
    
    {
      if(is.null(col_title)==T){
        guides(colour = F)
      }
    } +
    
    {
      if(is.null(linetype_title)==T){
        guides(linetype = F)
      }
    } +
    
    {
      if(is.null(size_title)==T){
        guides(size = F)
      }
    } +
    
    # guides(colour = col_title) +
    # guides(shape = shape_title) +
    # guides(linetype = linetype_title) +
    {if(is.null(ylab)){
      ylab(bquote("β"))
    } else{
      ylab(paste0(ylab))
    }} +
    xlab(xlab) +
    
    
    #minor fixes ----
  # theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust=0.5)) +
  scale_y_continuous(labels = signs_format(), #ensure y axis labels use minus sign, not hyphen
                     breaks = breaks) + 
    # if "breaks" need to changed from default, they can be specified here. 
    
    # if breaks were changed, the limits also need to be changed. Otherwise this can be left blank.
    coord_cartesian(ylim = limits) +
    {
      if(coord_flip == T){ 
        coord_flip(ylim = limits) # flip coordinate system
      }
    } 
  
  
  return(fig)
}


gg_geom_line <- function(
    # basic
  .data, x, y=beta, ymin=lb, ymax=ub, p_dodge = F, p_dodge_width = 1,
  size = NULL, colour = NULL, linetype = NULL,   
  
  # Faceting
  facet_type = F, vert = "", hori = "", scales = "fixed", nrow = NULL, ncol = NULL, 
  
  # Aestetic
  theme = 1, qr_edit = F, rect_col = "gray", yint = 0,
  size_vec = NULL, linetype_vec = NULL, colour_vec = NULL, xlab_vec = NULL,
  panel.spacing.x = 0.075, panel.spacing.y = 0.075,
  
  # titles
  title = NULL, subtitle = NULL, caption = NULL, 
  ylab = NULL, xlab = waiver(),
  col_title = waiver(), size_title = waiver(), linetype_title = waiver(),
  
  # Fix font
  text_size_minor = 8, text_size_major = 9, text_size_title = 10,
  text_family = "sans",
  
  # legend
  legend.position = "bottom", legend.title.align = 0.5, legend.direction = "vertical", 
  legend.box = "vertical", legend.margin_inch = 0.00625, legend.spacing = 0.00625,
  
  # misc
  coord_flip = F,  limits = c(), breaks = waiver(), x_angle = 0, ystrip_angle = 270
) {
  
  # if breaks were specified, then limits need to be specified too. 
  
  if(is.vector(breaks) == T){
    lower_limit <- min(breaks)
    upper_limit <- max(breaks)
    limits <- c(lower_limit, upper_limit)
  }
  
  fig <- ggplot() + 
    
    # Add details if Quantile regression is used
    {
      if (qr_edit == T){
        geom_vline(xintercept = 5.5) 
      }
    } +
    
    {
      if (qr_edit == T){
        geom_rect(data = .data,
                  aes(ymin = -Inf, ymax = Inf, xmin = 5.5, xmax = Inf), 
                  fill = rect_col, alpha = 0.1)
      }
    } +
    
    # map results ----
  
  geom_line(data = .data,
            aes(x = {{x}}, y = {{y}}, 
                size = {{size}}, color = {{colour}}, linetype = {{linetype}})) + 
    geom_errorbar(data = .data,
                  aes(x = {{x}}, y = {{y}}, ymin = {{ymin}}, ymax = {{ymax}}, 
                      size = {{size}}, color = {{colour}}, linetype = {{linetype}} )) +
    
    # {
    #   if(p_dodge == T){
    #   geom_line(data = .data,
    #             aes(x = {{x}}, y = {{y}}, 
    #                 size = {{size}}, color = {{colour}}, linetype = {{linetype}}),
    #             position = position_dodge(width = p_dodge_width)) + 
    #   geom_errorbar(data = .data,
    #                 aes(x = {{x}}, y = {{y}}, ymin = {{ymin}}, ymax = {{ymax}}, 
    #                     size = {{size}}, color = {{colour}}, linetype = {{linetype}} ),
    #                 position = position_dodge(width = p_dodge_width)) 
  #   } 
  # } +
  # 
  # {
  #   if(p_dodge == F){
  #     geom_line(data = .data,
  #               aes(x = {{x}}, y = {{y}}, 
  #                   size = {{size}}, color = {{colour}}, linetype = {{linetype}})) + 
  #     geom_errorbar(data = .data,
  #                   aes(x = {{x}}, y = {{y}}, ymin = {{ymin}}, ymax = {{ymax}}, 
  #                       size = {{size}}, color = {{colour}}, linetype = {{linetype}} )) 
  #   } 
  # } +
  
  
  
  # organize facets----
  {
    if(facet_type == "grid"){
      facet_grid(paste0({{vert}}, '~', {{hori}}), scales = scales)
    } 
  } +
    
    {
      if(facet_type == "wrap") {
        facet_wrap(paste0({{vert}}, "~", {{hori}}), scales = scales, nrow= nrow, ncol = ncol)
      } 
    } +
    
    # aestetic ----
  {
    if(is.null(yint) == F){
      geom_hline(yintercept = yint, lty =5, size=0.5) #set y-intercept = 0 (the null value)
    }
  } +
    
    ## theme 1 is black and white
    {
      if(theme == 1){
        theme_bw()
      } 
    } +
    
    ## theme 2 more minimal
    ### add lines to the bottom and right edge of each facet
    {
      if (theme == 2){
        annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, size=1)
      }
    } +
    {
      if (theme == 2){
        annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, size=1) 
      }
    } +
    
    {
      if (theme == 2){
        theme_minimal() + # gets of the boxes around facet titles.   
          theme(plot.background = element_blank(),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.border = element_blank() ) +
          theme(strip.text = element_text(face = "bold")) 
        
      }
    } +
    
    {
      if(is.null(xlab_vec)==F){
        scale_x_discrete(labels = xlab_vec) # set x axis labels manually
      }
    } + 
    
    {
      if(is.null(colour_vec)==F){
        scale_colour_manual(values = colour_vec) # set colour scheme manually
      }
    } +
    
    {
      if(is.null(linetype_vec)==F){
        scale_linetype_manual(values = linetype_vec) # set linetype scheme manually
      }
    } +
    
    {
      if(is.null(size_vec)==F){
        scale_size_manual(values = size_vec) # set size scheme manually
      }
    } +
    
    ## add option to increase facet margins
    {
      if(facet_type == "grid"  | facet_type == "wrap"){
        theme(panel.spacing.x = unit(panel.spacing.x, "inch")) 
      } 
    } +
    
    {
      if(facet_type == "grid"  | facet_type == "wrap"){
        theme(panel.spacing.y = unit(panel.spacing.y, "inch")) 
      } 
    } +
    
    
    # fix font ----
  theme(text = element_text(family = text_family, colour = "black", size = text_size_minor)) +
    theme(axis.text.x = element_text(size=text_size_minor, colour = "black", angle = x_angle)) + 
    theme(axis.text.y = element_text(size=text_size_minor, colour = "black")) + 
    theme(axis.ticks = element_line(size = 0.5, colour = "black")) + 
    
    theme(axis.title.x = element_text(size=text_size_major, colour = "black")) + 
    theme(axis.title.y = element_text(size=text_size_major, colour = "black")) +
    # theme(axis.title.y = element_text(margin = text_size_major(0,0,0,0, unit="in"))) +
    # theme(axis.title.x = element_text(margin = margin(0,0,0,0, unit="in"))) +
    
    theme(strip.text.x = element_text(size=text_size_major, colour = "black")) +
    theme(strip.text.y = element_text(size=text_size_major, colour = "black", angle = ystrip_angle)) +
    theme(legend.text = element_text(size=text_size_major, colour = "black")) +
    theme(legend.title = element_text(size=text_size_major, colour = "black")) +
    
    
    theme(title = element_text(size = text_size_title, colour = "black")) +
    
    
    
    # legend ----
  theme(legend.position = legend.position, legend.title.align = legend.title.align,
        legend.background = element_rect(size=0.5, linetype="solid"),
        legend.direction = legend.direction, legend.box = legend.box,
        legend.box.margin = margin(legend.margin_inch, legend.margin_inch,
                                   legend.margin_inch, legend.margin_inch, unit = "inch"),
        legend.spacing = unit(legend.spacing, units = "inch")
  ) +
    
    #labels ----
  labs(title = title, subtitle = subtitle, caption = caption) +
    labs(col = col_title) +
    labs(size = size_title) +
    labs(linetype = linetype_title) +
    
    {
      if(is.null(size_title)==T){
        guides(size = F)
      }
    } +
    
    # guides(colour = col_title) +
    # guides(size = size_title) +
    # guides(linetype = linetype_title) +
    {
      if(is.null(ylab)){
        ylab(bquote("β"))
      } else{
        ylab(paste0(ylab))
      }
    } +
    xlab(xlab) +
    
    
    #minor fixes ----
  # theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust=0.5)) +
  scale_y_continuous(labels = signs_format(), #ensure y axis labels use minus sign, not hyphen
                     breaks = breaks) + 
    # if "breaks" need to changed from default, they can be specified here. 
    
    # if breaks were changed, the limits also need to be changed. Otherwise this can be left blank.
    coord_cartesian(ylim = limits) +
    {
      if(coord_flip == T){ 
        coord_flip(ylim = limits) # flip coordinate system
      }
    } 
  
  return(fig)
}









gg_geom_boxplot <- function(data, p_dodge = T, theme = 1,
                            x, y=beta, colour = NULL, linetype = NULL,
                            facet_type = NULL, vert = "", hori = "", scales = "fixed", nrow = NULL, ncol = NULL, 
                            coord_flip = F,
                            shape_vec = NULL, linetype_vec = NULL, colour_vec = NULL,
                            legend.position = "bottom", legend.title.align = 0.5, legend.direction = "vertical",
                            title = NULL, subtitle = NULL, caption = NULL, ylab = NULL, xlab = waiver(),
                            text_size_minor = 9, text_size_major = 10,
                            col_title = waiver(), shape_title = waiver(), linetype_title = waiver(),
                            breaks = waiver(), ystrip_angle = 270) {
  
  # if breaks were specified, then limits need to be specified too. 
  limits <- c() # left empty if breaks arguement is not provided
  if(is.vector(breaks) == T){
    lower_limit <- min(breaks)
    upper_limit <- max(breaks)
    limits <- c(lower_limit, upper_limit)
  }
  
  fig <- ggplot() + 
    
    
    
    # map results ----
  {
    if(p_dodge == T){
      geom_boxplot(data = data,
                   aes(x = {{x}}, y = {{y}}, color = {{colour}}, linetype = {{linetype}}),
                   position = position_dodge(width = 1))
    }
  } +
    
    {
      if(p_dodge == F){
        geom_boxplot(data = data,
                     aes(x = {{x}}, y = {{y}}, color = {{colour}}, linetype = {{linetype}}))
      }
    } +
    
    
    # organize facets----
  {
    if(facet_type == "grid"){
      facet_grid(paste0({{vert}}, '~', {{hori}}), scales = scales)
    } 
  } +
    {
      if(facet_type == "wrap") {
        facet_wrap(paste0({{vert}}, "~", {{hori}}), scales = scales, nrow= nrow, ncol = ncol)
      } 
    } +
    
    # aestetic ----
  # geom_hline(yintercept = 0, lty =5, size=0.5) + #set y-intercept = 0 (the null value)
  
  ## theme 1 is black and white
  {
    if (theme == 1){
      theme_bw()
    } 
  } +
    
    ## theme 2 more minimal
    ### add lines to the bottom and right edge of each facet
    {
      if (theme == 2){
        annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, size=1)
      }
    } +
    {
      if (theme == 2){
        annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, size=1) 
      }
    } +
    
    {
      if (theme == 2){
        theme_minimal() + # gets of the boxes around facet titles.   
          theme(plot.background = element_blank(),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.border = element_blank() ) +
          theme(strip.text = element_text(face = "bold")) 
        
      }
    } +
    
    {
      if(is.null(colour_vec)==F){
        scale_colour_manual(values = colour_vec) #set colour scheme manually
      }
    } +
    
    {
      if(is.null(linetype_vec)==F){
        scale_linetype_manual(values = linetype_vec) #set linetype scheme manually
      }
    } +
    
    {
      if(is.null(shape_vec)==F){
        scale_shape_manual(values = shape_vec) #set shape scheme manually
      }
    } +
    
    # fix font ----
    theme(text = element_text(family = "sans", colour = "black", size = text_size_minor)) +
    theme(axis.text.x = element_text(size=text_size_minor, colour = "black")) + 
    theme(axis.text.y = element_text(size=text_size_minor, colour = "black")) + 
    theme(axis.ticks = element_line(size = 0.5, colour = "black")) + 
    
    theme(axis.title.x = element_text(size=text_size_minor, colour = "black")) + 
    theme(axis.title.y = element_text(size=text_size_minor, colour = "black")) +
    # theme(axis.title.y = element_text(margin = margin(0,0,0,0, unit="in"))) +
    # theme(axis.title.x = element_text(margin = margin(0,0,0,0, unit="in"))) +
    
    theme(strip.text.x = element_text(size=text_size_major, colour = "black")) +
    theme(strip.text.y = element_text(size=text_size_major, colour = "black", angle = ystrip_angle)) +
    theme(legend.text = element_text(size=text_size_major, colour = "black")) +
    theme(legend.title = element_text(size=text_size_major, colour = "black")) +
    
    theme(title = element_text(size = text_size_major, colour = "black")) +
    
    # legend ----
  theme(legend.position = legend.position, legend.title.align = legend.title.align,
        legend.background = element_rect(size=0.5, linetype="solid"),
        legend.direction = legend.direction) +
    
    #labels ----
  labs(title = title, subtitle = subtitle, caption = caption) +
    labs(col = col_title) +
    labs(shape = shape_title) +
    labs(linetype = linetype_title) +
    {if(is.null(ylab)){
      ylab(bquote("β"))
    } else{
      ylab(print(ylab))
    }} +
    xlab(xlab) +
    
    
    #minor fixes ----
  # theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust=0.5)) +
  scale_y_continuous(labels = signs_format(), #ensure y axis labels use minus sign, not hyphen
                     breaks = breaks) + 
    # if "breaks" need to changed from default, they can be specified here. 
    
    # if breaks were changed, the limits also need to be changed. Otherwise this can be left blank.
    coord_cartesian(ylim = limits) +
    {
      if(coord_flip == T){ 
        coord_flip(ylim = limits) # flip coordinate system
      }
    } 
  return(fig)
}
