# TEST 1: SOLAR ############################################################
solar <- calc_daily_solar(lat=37.7, n=172, elv=142, y=2000, sf=1, tc=23.0)
cat(sprintf("TEST 1---Solar values:\n"))
cat(sprintf("  kn: %d\n", solar$kN))
cat(sprintf("  nu: %0.6f degrees\n", solar$nu_deg))
cat(sprintf("  lambda: %0.6f degrees\n", solar$lambda_deg))
cat(sprintf("  rho: %0.6f\n", solar$rho))
cat(sprintf("  dr: %0.6f\n", solar$dr))
cat(sprintf("  delta: %0.6f degrees\n", solar$delta_deg))
cat(sprintf("  ru: %0.6f\n", solar$ru))
cat(sprintf("  rv: %0.6f\n", solar$rv))
cat(sprintf("  rw: %0.6f\n", solar$rw))
cat(sprintf("  hs: %0.6f degrees\n", solar$hs_deg))
cat(sprintf("  hn: %0.6f degrees\n", solar$hn_deg))
cat(sprintf("  tau_o: %0.6f\n", solar$tau_o))
cat(sprintf("  tau: %0.6f\n", solar$tau))
cat(sprintf("  Qn: %0.6f mol/m^2\n", solar$ppfd_mol.m2))
cat(sprintf("  Rnl: %0.6f w/m^2\n", solar$rnl_w.m2))
cat(sprintf("  Ho: %0.6f MJ/m^2\n", (1.0e-6) * solar$ra_j.m2))
cat(sprintf("  Hn: %0.6f MJ/m^2\n", (1.0e-6) * solar$rn_j.m2))
cat(sprintf("  Hnn: %0.6f MJ/m^2\n", (1.0e-6) * solar$rnn_j.m2))

# TEST 2: EVAP ############################################################
evap <- calc_daily_evap(lat=37.7, n=172, elv=142, y=2000, sf=1, tc=23.0, sw=0.9)
cat(sprintf("TEST 2---Evap values:\n"))
cat(sprintf("  s: %0.6f Pa/K\n", evap$s_pa.k))
cat(sprintf("  Lv: %0.6f MJ/kg\n", (1e-6) * evap$lv_j.kg))
cat(sprintf("  Patm: %0.6f bar\n", (1e-5) * evap$patm_pa))
cat(sprintf("  pw: %0.6f kg/m^3\n", evap$pw_kg.m3))
cat(sprintf("  gamma: %0.6f Pa/K\n", evap$gam_pa.k))
cat(sprintf("  Econ: %0.6f mm^3/J\n", (1e9) * evap$econ_m3.j))
cat(sprintf("  Cn: %0.6f mm\n", evap$cond_mm))
cat(sprintf("  rx: %0.6f\n", evap$rx))
cat(sprintf("  hi: %0.6f degrees\n", evap$hi_deg))
cat(sprintf("  EET: %0.6f mm\n", evap$eet_mm))
cat(sprintf("  PET: %0.6f mm\n", evap$pet_mm))
cat(sprintf("  AET: %0.6f mm\n", evap$aet_mm))

# TEST 3: SOIL MOISTURE (run one day) #########################################
t3 <- run_one_day(lat=37.7, elv=142, n=172, y=2000, wn=75, sf=1, tc=23, pn=5)
cat(sprintf("TEST 3---Soil moisture (run one day):\n"))
cat(sprintf("  Ho: %0.6f J/m2\n", t3$ho))
cat(sprintf("  Hn: %0.6f J/m2\n", t3$hn))
cat(sprintf("  PPFD: %0.6f mol/m2\n", t3$ppfd))
cat(sprintf("  EET: %0.6f mm/d\n", t3$eet))
cat(sprintf("  PET: %0.6f mm/d\n", t3$pet))
cat(sprintf("  AET: %0.6f mm/d\n", t3$aet))
cat(sprintf("  Cn: %0.6f mm/d\n", t3$cond))
cat(sprintf("  Wn: %0.6f mm\n", t3$wn))
cat(sprintf("  RO: %0.6f mm\n", t3$ro))

# # TEST 4:
daily_totals <- matrix(data=rep(0, 366), nrow=366, ncol=1)
daily_totals <- as.data.frame(daily_totals)
names(daily_totals) <- c("wn")
my_lat <- 37.7
my_elv <- 142
my_file <- system.file("extdata/example_data.csv", package = "splash")
my_data <- read_csv(my_file, 2000)
my_data$lat_deg <- my_lat
my_data$elv_m <- my_elv
daily_totals <- spin_up(my_data, daily_totals)
cat(sprintf("TEST 4---Spin-Up:\n"))
cat(sprintf("Day,Wn (mm)\n"))
for (i in seq(from=1, to=my_data$num_lines, by=1)) {
  cat(sprintf("%d,%0.6f\n", i, daily_totals$wn[i]))
}
