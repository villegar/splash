hex_logo <- function(subplot = system.file("images/splash.png",
                                           package = "documentation"),
                     dpi = 600,
                     h_color = "#000000",
                     h_fill = "#FFFFFF",
                     output = system.file("images/logo.png",
                                          package = "splash"),
                     package = "SPLASH",
                     p_color = "#000000",
                     url = "https://github.com/villegar/splash",
                     u_size = 15) {
  hexSticker::sticker(subplot = subplot, package = package,
                      h_color = h_color,  h_fill = h_fill,
                      dpi = dpi,
                      s_x = 1.0, s_y = .9, s_width = .62, asp = .85,
                      p_x = 1.0, p_y = 1.58, p_size = 35, p_color = p_color,
                      url = url,
                      u_angle = 30, u_color = p_color, u_size = u_size,
                      filename = output)
}

hex_logo("inst/images/splash.png", output = "inst/images/logo.png", u_size = 9.3, dpi = 600)
