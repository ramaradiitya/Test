# =============================================================================
# GTWR - KOMBINASI 1 (LENGKAP)
# Y = Persentase_Penduduk_Miskin
# X = Persentase_Kekerasan_Perempuan + Persentase_KK_Perempuan + RLS_Perempuan + IPG
# =============================================================================
# PERUBAHAN TERBARU:
#   [BARU] Ditambahkan UJI MULTIKOLINEARITAS (VIF) pada model OLS Y ~ X,
#          lengkap dengan nilai VIF, Tolerance, keputusan per variabel, dan
#          RINGKASAN di bawah tabelnya (bagian 5b). Ringkasan VIF juga ikut
#          muncul di RINGKASAN AKHIR.
#   [SEBELUMNYA] Uji NORMALITAS residual GTWR memakai KOLMOGOROV-SMIRNOV
#          (ks.test dengan mean & sd dari residual), menggantikan Lilliefors.
#
# TAMBAHAN DARI VERSI SEBELUMNYA:
#   (A) Moran's I SPASIAL PER TAHUN (2021, 2022, 2023, 2024 masing-masing,
#       27 kab/kota per tahun) -> semua tahun harus p < 0.05.
#       (Moran's I pasangan tahun yang lama tetap dipertahankan sebagai info.)
#   (B) Uji HOMOGENITAS residual GTWR yang lebih lengkap:
#       1. Breusch-Pagan (versi Koenker) terhadap X
#       2. Breusch-Pagan terhadap Latitude & Longitude (heterogenitas spasial)
#       3. Glejser
#       4. Park
#       5. White
#       Keputusan: homogen jika p >= 0.05.
#   (C) Bonus: Moran's I residual GTWR per tahun (cek sisa autokorelasi spasial).
#
# CATATAN:
#   - moran.test() memakai alternative = "greater" (default), yaitu menguji
#     autokorelasi spasial POSITIF.
#   - Pembobot spasial per tahun = k-nearest neighbour (k = 5), row-standardized,
#     sama dengan pembobot yang dipakai pada Moran global.
#   - Kolmogorov-Smirnov di sini memakai mean & sd yang DIESTIMASI dari residual.
#     Pada kondisi ini p-value KS cenderung konservatif (lebih besar dari
#     seharusnya). Sebutkan sebagai catatan di laporan.
#   - VIF dihitung dari regresi OLS global (pooled seluruh panel) dengan X yang
#     sama dengan model GTWR. Patokan: VIF < 5 aman, 5-10 perlu perhatian,
#     >= 10 indikasi multikolinearitas serius.
#   - Tidak ada file yang disimpan. Semua hasil ada di objek R.
# =============================================================================

library(GWmodel)
library(spdep)
library(sp)
library(readxl)
library(lmtest)
library(car)      # untuk vif() -- dimuat SEBELUM dplyr supaya recode() milik dplyr yang dipakai
library(dplyr)

options(max.print = 1e6)

# =============================================================================
# 1. INPUT DATA  (SESUAIKAN PATH FILE)
# =============================================================================
data_asli <- read_excel("C:/Users/advan/Downloads/2021-2024 panell cerai.xlsx",
                        sheet = "Sheet1") %>%
  arrange(Tahun, `Kabupaten/Kota`)

stopifnot(nrow(data_asli %>% filter(is.na(Latitude) | is.na(Longitude))) == 0)

y_var  <- "Persentase_Penduduk_Miskin"
x_vars <- c("Persentase_Kekerasan_Perempuan", "Persentase_KK_Perempuan",
            "RLS_Perempuan", "IPG")
formula_obj <- as.formula(paste(y_var, "~", paste(x_vars, collapse = " + ")))

daftar_tahun <- sort(unique(data_asli$Tahun))
n_tahun      <- length(daftar_tahun)
n_wilayah    <- length(unique(data_asli$`Kabupaten/Kota`))
stopifnot(nrow(data_asli) == n_tahun * n_wilayah)

cat("Jumlah kabupaten/kota:", n_wilayah, "| Jumlah tahun:", n_tahun,
    "| Total baris:", nrow(data_asli), "\n")
cat("Model:", deparse(formula_obj), "\n\n")

# =============================================================================
# 2. PEMBOBOT SPASIAL
# =============================================================================
coord_wil <- data_asli %>%
  distinct(`Kabupaten/Kota`, Latitude, Longitude) %>%
  arrange(`Kabupaten/Kota`)

nb_wil <- knn2nb(knearneigh(as.matrix(coord_wil[, c("Longitude", "Latitude")]),
                            k = min(5, n_wilayah - 1)))

# (a) pembobot per tahun (27 wilayah) -> untuk Moran PER TAHUN
listw_wil <- nb2listw(nb_wil, style = "W", zero.policy = TRUE)

# (b) pembobot gabungan seluruh panel (spasial di dalam tahun yang sama)
nb_full_spatial <- vector("list", n_wilayah * n_tahun)
for (t in 0:(n_tahun - 1)) {
  for (p in 1:n_wilayah) {
    nb_full_spatial[[t * n_wilayah + p]] <- as.integer(t * n_wilayah + nb_wil[[p]])
  }
}
class(nb_full_spatial) <- "nb"
listw_spatial_global <- nb2listw(nb_full_spatial, style = "W", zero.policy = TRUE)

# (c) pembobot pasangan tahun (wilayah yang sama, tahun berbeda) -> info saja
pasangan_tahun <- lapply(1:(n_tahun - 1), function(i) c(daftar_tahun[i], daftar_tahun[i + 1]))
listw_temporal_pertahun <- lapply(pasangan_tahun, function(pair) {
  sub_df <- data_asli %>%
    mutate(row_id = row_number()) %>%
    filter(Tahun %in% pair) %>%
    arrange(Tahun, `Kabupaten/Kota`)
  nb_pair <- vector("list", nrow(sub_df))
  for (j in seq_len(nrow(sub_df))) {
    tahun_lain <- setdiff(pair, sub_df$Tahun[j])
    nb_pair[[j]] <- as.integer(which(sub_df$`Kabupaten/Kota` == sub_df$`Kabupaten/Kota`[j] &
                                       sub_df$Tahun == tahun_lain))
  }
  class(nb_pair) <- "nb"
  list(row_id = sub_df$row_id,
       listw  = nb2listw(nb_pair, style = "W", zero.policy = TRUE),
       label  = paste(pair, collapse = "-"))
})

# =============================================================================
# 3. FUNGSI MORAN'S I
# =============================================================================

# --- (A) Moran's I SPASIAL PER TAHUN ---
moran_per_tahun <- function(vektor_all) {
  hasil <- data.frame()
  for (th in daftar_tahun) {
    idx <- which(data_asli$Tahun == th)
    # pastikan urutan wilayah sama dengan coord_wil
    stopifnot(all(data_asli$`Kabupaten/Kota`[idx] == coord_wil$`Kabupaten/Kota`))
    mt <- tryCatch(moran.test(vektor_all[idx], listw_wil, zero.policy = TRUE),
                   error = function(e) NULL)
    hasil <- rbind(hasil, data.frame(
      Tahun      = th,
      Moran_I    = if (!is.null(mt)) as.numeric(mt$estimate[1]) else NA,
      Expected_I = if (!is.null(mt)) as.numeric(mt$estimate[2]) else NA,
      Z_hitung   = if (!is.null(mt)) as.numeric(mt$statistic) else NA,
      P_Value    = if (!is.null(mt)) mt$p.value else NA,
      Signifikan = if (!is.null(mt)) mt$p.value < 0.05 else NA
    ))
  }
  hasil
}

# --- Moran's I spasial global (seluruh panel) ---
moran_spasial_global <- function(vektor_all) {
  tryCatch(moran.test(vektor_all, listw_spatial_global, zero.policy = TRUE),
           error = function(e) NULL)
}

# --- Moran's I pasangan tahun (adaptasi, info tambahan) ---
moran_pasangan_tahun <- function(vektor_all) {
  hasil <- data.frame()
  for (obj in listw_temporal_pertahun) {
    mt <- tryCatch(moran.test(vektor_all[obj$row_id], obj$listw, zero.policy = TRUE),
                   error = function(e) NULL)
    hasil <- rbind(hasil, data.frame(
      Periode = obj$label,
      Moran_I = if (!is.null(mt)) as.numeric(mt$estimate[1]) else NA,
      P_Value = if (!is.null(mt)) mt$p.value else NA
    ))
  }
  hasil
}

# =============================================================================
# 4. FUNGSI UJI HOMOGENITAS RESIDUAL (BP, BP spasial, Glejser, Park, White)
# =============================================================================
p_dari_F <- function(m) {
  f <- summary(m)$fstatistic
  unname(pf(f[1], f[2], f[3], lower.tail = FALSE))
}

uji_homogenitas_lengkap <- function(residual, x_vars) {
  X  <- data_asli[x_vars]
  n  <- length(residual)
  e2 <- residual^2
  
  hasil <- data.frame()
  tambah <- function(nama, stat, df, p) {
    hasil <<- rbind(hasil, data.frame(
      Uji = nama, Statistik = round(stat, 4), df = df, P_Value = round(p, 6),
      Keputusan = ifelse(is.na(p), "NA",
                         ifelse(p >= 0.05, "Homogen", "Tidak homogen")),
      stringsAsFactors = FALSE))
  }
  
  # 1. Breusch-Pagan (Koenker): LM = n * R2 dari regresi e^2 ~ X
  aux1 <- lm(e2 ~ ., data = data.frame(e2 = e2, X))
  LM1  <- n * summary(aux1)$r.squared
  tambah("Breusch-Pagan (terhadap X)", LM1, length(x_vars),
         pchisq(LM1, length(x_vars), lower.tail = FALSE))
  
  # 2. Breusch-Pagan terhadap Latitude & Longitude
  aux2 <- lm(e2 ~ Latitude + Longitude,
             data = data.frame(e2 = e2, Latitude = data_asli$Latitude,
                               Longitude = data_asli$Longitude))
  LM2  <- n * summary(aux2)$r.squared
  tambah("Breusch-Pagan (terhadap Lat & Long)", LM2, 2,
         pchisq(LM2, 2, lower.tail = FALSE))
  
  # 3. Glejser: |e| ~ X, uji F serentak
  aux3 <- lm(ae ~ ., data = data.frame(ae = abs(residual), X))
  F3   <- unname(summary(aux3)$fstatistic[1])
  tambah("Glejser (|e| ~ X)", F3, length(x_vars), p_dari_F(aux3))
  
  # 4. Park: ln(e^2) ~ ln(X), uji F serentak (butuh X > 0 dan e != 0)
  if (all(X > 0) && all(e2 > 0)) {
    aux4 <- lm(le2 ~ ., data = data.frame(le2 = log(e2), log(X)))
    F4   <- unname(summary(aux4)$fstatistic[1])
    tambah("Park (ln e^2 ~ ln X)", F4, length(x_vars), p_dari_F(aux4))
  } else {
    cat("Catatan: Uji Park dilewati (ada X <= 0 atau residual = 0, ln tidak terdefinisi)\n")
  }
  
  # 5. White: e^2 ~ X + X^2 + interaksi antar X; LM = n * R2
  f_white <- as.formula(paste("e2 ~ (", paste(x_vars, collapse = " + "), ")^2 + ",
                              paste0("I(", x_vars, "^2)", collapse = " + ")))
  aux5 <- lm(f_white, data = data.frame(e2 = e2, X))
  df5  <- sum(!is.na(coef(aux5))) - 1
  LM5  <- n * summary(aux5)$r.squared
  tambah("White", LM5, df5, pchisq(LM5, df5, lower.tail = FALSE))
  
  hasil
}

# =============================================================================
# 5. UJI SERENTAK & PARSIAL (tidak diubah)
# =============================================================================
uji_serentak <- function(model_terbaik, formula_obj) {
  reg_global <- lm(formula_obj, data = data_asli)
  RSS0 <- sum(resid(reg_global)^2); df0 <- reg_global$df.residual
  RSS1 <- model_terbaik$GTW.diagnostic$RSS.gw; df1 <- model_terbaik$GTW.diagnostic$edf
  F_hit <- (RSS0 / df0) / (RSS1 / df1)
  p_val <- pf(F_hit, df0, df1, lower.tail = FALSE)
  list(F_hit = F_hit, p_value = p_val, reg_global = reg_global)
}

uji_parsial <- function(model_terbaik, reg_global, x_vars) {
  sdf <- as.data.frame(model_terbaik$SDF)
  edf <- model_terbaik$GTW.diagnostic$edf
  t_crit <- qt(0.975, df = edf)
  ringkasan_global <- summary(reg_global)$coefficients
  semua_var_cek <- c("Intercept", x_vars)
  
  tabel_sig <- data.frame(Wilayah = data_asli$`Kabupaten/Kota`, Tahun = data_asli$Tahun)
  for (v in semua_var_cek) {
    nama_global <- ifelse(v == "Intercept", "(Intercept)", v)
    kolom_tv <- paste0(v, "_TV")
    if (kolom_tv %in% names(sdf)) {
      tabel_sig[[paste0(v, "_t_value")]]    <- round(sdf[[kolom_tv]], 4)
      tabel_sig[[paste0(v, "_signifikan")]] <- abs(sdf[[kolom_tv]]) > t_crit
    } else if (nama_global %in% rownames(ringkasan_global) && v %in% names(sdf)) {
      p_global    <- ringkasan_global[nama_global, "Pr(>|t|)"]
      arah_global <- sign(ringkasan_global[nama_global, "Estimate"])
      tabel_sig[[paste0(v, "_t_value")]]    <- NA
      tabel_sig[[paste0(v, "_signifikan")]] <- (p_global < 0.05) & (sign(sdf[[v]]) == arah_global)
    }
  }
  
  kolom_sig <- setdiff(grep("_signifikan$", names(tabel_sig), value = TRUE), "Intercept_signifikan")
  tabel_sig$Jumlah_Var_Signifikan <- rowSums(tabel_sig[kolom_sig])
  tabel_sig$Variabel_Signifikan <- apply(tabel_sig[kolom_sig], 1, function(row) {
    nama <- gsub("_signifikan", "", names(row)[row == TRUE])
    if (length(nama) == 0) return("Tidak ada")
    paste(nama, collapse = ", ")
  })
  
  ringkasan_per_variabel <- data.frame(
    Variabel = x_vars,
    Jumlah_Lokasi_Signifikan = sapply(x_vars, function(v) sum(tabel_sig[[paste0(v, "_signifikan")]], na.rm = TRUE)),
    Proporsi_Signifikan = sapply(x_vars, function(v) mean(tabel_sig[[paste0(v, "_signifikan")]], na.rm = TRUE))
  )
  list(tabel_sig = tabel_sig, ringkasan_per_variabel = ringkasan_per_variabel,
       rata2_var_signifikan = mean(tabel_sig$Jumlah_Var_Signifikan),
       proporsi_lokasi_signifikan = mean(tabel_sig$Jumlah_Var_Signifikan > 0))
}

# =============================================================================
# 5b. UJI MULTIKOLINEARITAS (VIF)   [BARU]
# =============================================================================
# VIF dihitung dari regresi OLS global Y ~ X (pooled seluruh panel).
# VIF_j = 1 / (1 - R2_j), R2_j = R2 dari regresi X_j terhadap X lainnya.
# Tolerance = 1 / VIF.
# Patokan: VIF < 5 aman | 5 <= VIF < 10 perlu perhatian | VIF >= 10 serius.
cat("============================================================\n")
cat("UJI MULTIKOLINEARITAS (VIF) - OLS Y ~ X\n")
cat("============================================================\n")

ols_global <- lm(formula_obj, data = data_asli)

vif_nilai <- car::vif(ols_global)
# vif() mengembalikan vektor (atau matriks GVIF kalau ada prediktor kategorik;
# di sini semua X numerik sehingga berupa vektor)
vif_vektor <- as.numeric(vif_nilai)

tabel_vif <- data.frame(
  Variabel  = names(vif_nilai),
  VIF       = round(vif_vektor, 4),
  Tolerance = round(1 / vif_vektor, 4),
  Keputusan = ifelse(vif_vektor >= 10, "Multikolinearitas serius",
                     ifelse(vif_vektor >= 5, "Perlu perhatian (VIF 5-10)",
                            "Tidak ada multikolinearitas")),
  stringsAsFactors = FALSE
)
print(tabel_vif, row.names = FALSE)

# ---------------- RINGKASAN VIF ----------------
idx_maks <- which.max(vif_vektor)
idx_min  <- which.min(vif_vektor)
n_serius <- sum(vif_vektor >= 10)
n_perhat <- sum(vif_vektor >= 5 & vif_vektor < 10)

kesimpulan_vif <- if (n_serius > 0) {
  paste0("Ada ", n_serius, " variabel dengan VIF >= 10 -> terindikasi multikolinearitas serius: ",
         paste(names(vif_nilai)[vif_vektor >= 10], collapse = ", "))
} else if (n_perhat > 0) {
  paste0("Tidak ada VIF >= 10, tetapi ", n_perhat, " variabel berada di rentang 5-10 (perlu perhatian): ",
         paste(names(vif_nilai)[vif_vektor >= 5 & vif_vektor < 10], collapse = ", "))
} else {
  "Semua VIF < 5 -> tidak ada indikasi multikolinearitas antar variabel X"
}

cat("\n--- RINGKASAN VIF ---\n")
cat("Jumlah variabel X                   :", length(vif_vektor), "\n")
cat("VIF terbesar                        :", round(vif_vektor[idx_maks], 4),
    "(", names(vif_nilai)[idx_maks], ")\n")
cat("VIF terkecil                        :", round(vif_vektor[idx_min], 4),
    "(", names(vif_nilai)[idx_min], ")\n")
cat("Rata-rata VIF                       :", round(mean(vif_vektor), 4), "\n")
cat("Jumlah variabel VIF >= 10 (serius)  :", n_serius, "\n")
cat("Jumlah variabel VIF 5-10 (perhatian):", n_perhat, "\n")
cat("Kesimpulan                          :", kesimpulan_vif, "\n\n")

# =============================================================================
# 6. TAHAP 1 : MORAN'S I (global, PER TAHUN, pasangan tahun) + BP SPASIAL
# =============================================================================
cat("============================================================\n")
cat("TAHAP 1 : UJI EFEK SPASIAL & HETEROGENITAS SPASIAL\n")
cat("============================================================\n")

y_vec <- data_asli[[y_var]]

# 6a. Moran global
m_glob <- moran_spasial_global(y_vec)
p_glob <- if (!is.null(m_glob)) m_glob$p.value else NA
cat("\n[1a] Moran's I SPASIAL GLOBAL (Y):  I =",
    round(as.numeric(m_glob$estimate[1]), 5), "| p =", signif(p_glob, 6),
    ifelse(!is.na(p_glob) && p_glob < 0.05, "-> SIGNIFIKAN", "-> tidak signifikan"), "\n")

# 6b. Moran PER TAHUN
moran_tahun_df <- moran_per_tahun(y_vec)
cat("\n[1b] Moran's I SPASIAL PER TAHUN (Y):\n")
print(moran_tahun_df, row.names = FALSE)
semua_tahun_sig <- all(!is.na(moran_tahun_df$P_Value) & moran_tahun_df$P_Value < 0.05)
cat(" -> Semua tahun signifikan (p < 0.05):", semua_tahun_sig, "\n")
if (!semua_tahun_sig) {
  cat(" -> PERHATIAN, tahun yang TIDAK signifikan:",
      paste(moran_tahun_df$Tahun[!(moran_tahun_df$P_Value < 0.05)], collapse = ", "), "\n")
}

# 6c. Moran pasangan tahun (info tambahan)
moran_pair_df <- moran_pasangan_tahun(y_vec)
cat("\n[1c] Moran's I PASANGAN TAHUN (adaptasi, info tambahan):\n")
print(moran_pair_df, row.names = FALSE)
semua_pair_sig <- all(!is.na(moran_pair_df$P_Value) & moran_pair_df$P_Value < 0.05)

# 6d. Breusch-Pagan heterogenitas spasial (OLS terhadap Lat & Long)
# (ols_global sudah dibuat di bagian 5b)
p_bp_spasial <- tryCatch(bptest(ols_global, varformula = ~ Latitude + Longitude,
                                data = data_asli)$p.value, error = function(e) NA)
cat("\n[1d] Breusch-Pagan heterogenitas spasial (OLS): p =", signif(p_bp_spasial, 6),
    ifelse(!is.na(p_bp_spasial) && p_bp_spasial < 0.05, "-> ADA heterogenitas spasial", "-> tidak signifikan"), "\n")

tahap1_ok <- (!is.na(p_glob) && p_glob < 0.05) && semua_tahun_sig &&
  semua_pair_sig && (!is.na(p_bp_spasial) && p_bp_spasial < 0.05)
cat("\n=> LOLOS TAHAP 1 (Moran global + per tahun + pasangan tahun + BP spasial):", tahap1_ok, "\n\n")

# =============================================================================
# 7. TAHAP 2 : ESTIMASI 8 MODEL GTWR, PILIH AIC TERKECIL
# =============================================================================
data.sp <- data_asli
coordinates(data.sp) <- ~Longitude + Latitude
proj4string(data.sp) <- CRS("+proj=longlat +datum=WGS84")

run_gtwr <- function(kernel_type, adaptive_type) {
  bw <- tryCatch(bw.gtwr(formula_obj, data = data.sp, obs.tv = data.sp$Tahun,
                         kernel = kernel_type, adaptive = adaptive_type,
                         approach = "CV", longlat = TRUE), error = function(e) NA)
  if (is.na(bw)) return(NULL)
  model <- tryCatch(gtwr(formula_obj, data = data.sp, obs.tv = data.sp$Tahun,
                         st.bw = bw, kernel = kernel_type, adaptive = adaptive_type,
                         longlat = TRUE), error = function(e) NULL)
  if (is.null(model)) return(NULL)
  yhat <- model$SDF$yhat
  R2_val <- 1 - sum((y_vec - yhat)^2) / sum((y_vec - mean(y_vec))^2)
  list(model = model, AIC = model$GTW.diagnostic$AIC, AICc = model$GTW.diagnostic$AICc,
       R2 = R2_val, Adj_R2 = model$GTW.diagnostic$gwR2.adj,
       bandwidth = bw, kernel = kernel_type, adaptive = adaptive_type)
}

cat("============================================================\n")
cat("TAHAP 2 : ESTIMASI 8 MODEL GTWR\n")
cat("============================================================\n")

daftar_model <- list(
  "Gaussian Fixed"       = run_gtwr("gaussian", FALSE),
  "Gaussian Adaptive"    = run_gtwr("gaussian", TRUE),
  "Exponential Fixed"    = run_gtwr("exponential", FALSE),
  "Exponential Adaptive" = run_gtwr("exponential", TRUE),
  "Tricube Fixed"        = run_gtwr("tricube", FALSE),
  "Tricube Adaptive"     = run_gtwr("tricube", TRUE),
  "Bisquare Fixed"       = run_gtwr("bisquare", FALSE),
  "Bisquare Adaptive"    = run_gtwr("bisquare", TRUE)
)
daftar_model <- daftar_model[!sapply(daftar_model, is.null)]
stopifnot(length(daftar_model) > 0)

perbandingan_model <- data.frame(
  Model     = names(daftar_model),
  Bandwidth = sapply(daftar_model, function(m) round(m$bandwidth, 4)),
  AIC       = sapply(daftar_model, function(m) round(m$AIC, 2)),
  AICc      = sapply(daftar_model, function(m) round(m$AICc, 2)),
  R2        = sapply(daftar_model, function(m) round(m$R2, 4)),
  Adj_R2    = sapply(daftar_model, function(m) round(m$Adj_R2, 4)),
  row.names = NULL
) %>% arrange(AIC)
cat("\nPerbandingan 8 model (urut AIC terkecil):\n")
print(perbandingan_model, row.names = FALSE)

nama_terbaik  <- names(which.min(sapply(daftar_model, function(m) m$AIC)))
model_info    <- daftar_model[[nama_terbaik]]
model_terbaik <- model_info$model

cat("\nMODEL TERBAIK :", nama_terbaik, "\n")
cat("Kernel        :", model_info$kernel, "| Adaptive:", model_info$adaptive, "\n")
cat("Bandwidth     :", round(model_info$bandwidth, 4), "\n")
cat("R2            :", round(model_info$R2, 4), "| Adjusted R2:", round(model_info$Adj_R2, 4), "\n")
cat("AIC           :", round(model_info$AIC, 2), "| AICc      :", round(model_info$AICc, 2), "\n")

hs <- uji_serentak(model_terbaik, formula_obj)
cat("F-hitung (uji serentak):", round(hs$F_hit, 4), "| p-value:", signif(hs$p_value, 6), "\n\n")

# =============================================================================
# 8. TAHAP 3 : NORMALITAS RESIDUAL - KOLMOGOROV-SMIRNOV (WAJIB p >= 0.05)
# =============================================================================
sdf <- as.data.frame(model_terbaik$SDF)
residual_gtwr <- sdf$residual

ks_hasil <- tryCatch(
  ks.test(residual_gtwr, "pnorm",
          mean = mean(residual_gtwr), sd = sd(residual_gtwr)),
  error = function(e) NULL
)
ks_D <- if (!is.null(ks_hasil)) as.numeric(ks_hasil$statistic) else NA
norm_p <- if (!is.null(ks_hasil)) ks_hasil$p.value else NA

cat("============================================================\n")
cat("TAHAP 3 : NORMALITAS RESIDUAL GTWR (Kolmogorov-Smirnov)\n")
cat("============================================================\n")
cat("D =", round(ks_D, 5), "| p-value =", round(norm_p, 6), "->",
    ifelse(!is.na(norm_p) && norm_p >= 0.05, "NORMAL", "TIDAK NORMAL"), "\n\n")

# =============================================================================
# 9. HOMOGENITAS RESIDUAL GTWR (LENGKAP)
# =============================================================================
cat("============================================================\n")
cat("UJI HOMOGENITAS RESIDUAL GTWR (BP, BP spasial, Glejser, Park, White)\n")
cat("============================================================\n")
tabel_homogenitas <- uji_homogenitas_lengkap(residual_gtwr, x_vars)
print(tabel_homogenitas, row.names = FALSE)
semua_homogen <- all(tabel_homogenitas$P_Value >= 0.05, na.rm = TRUE)
cat("\n=> Semua uji homogenitas p >= 0.05 (residual homogen):", semua_homogen, "\n")
if (!semua_homogen) {
  cat("   Uji yang menolak H0:",
      paste(tabel_homogenitas$Uji[tabel_homogenitas$P_Value < 0.05], collapse = "; "), "\n")
}
cat("\n")

# =============================================================================
# 10. BONUS: MORAN'S I RESIDUAL GTWR PER TAHUN
# =============================================================================
cat("============================================================\n")
cat("MORAN'S I RESIDUAL GTWR PER TAHUN (p >= 0.05 = sisa autokorelasi spasial tidak ada)\n")
cat("============================================================\n")
moran_resid_df <- moran_per_tahun(residual_gtwr)
print(moran_resid_df, row.names = FALSE)
cat("\n")

# =============================================================================
# 11. UJI PARSIAL LOKAL (t-lokal)
# =============================================================================
hp <- uji_parsial(model_terbaik, hs$reg_global, x_vars)

cat("============================================================\n")
cat("UJI PARSIAL LOKAL (|t| > t-kritis, alpha = 5%)\n")
cat("============================================================\n")
cat("Rata-rata jumlah variabel signifikan per lokasi-tahun:", round(hp$rata2_var_signifikan, 2), "\n")
cat("Proporsi lokasi-tahun dengan >= 1 variabel signifikan :", round(hp$proporsi_lokasi_signifikan, 3), "\n\n")
cat("Ringkasan per variabel:\n")
print(hp$ringkasan_per_variabel, row.names = FALSE)

cat("\nDaftar variabel signifikan per wilayah & tahun:\n")
tabel_lengkap <- hp$tabel_sig %>% arrange(Tahun, Wilayah)
print(tabel_lengkap %>% select(Wilayah, Tahun, Variabel_Signifikan, Jumlah_Var_Signifikan))

cat("\nTabel LENGKAP nilai-t lokal per wilayah & tahun:\n")
print(tabel_lengkap)

# =============================================================================
# 12. RINGKASAN AKHIR
# =============================================================================
cat("\n\n================ RINGKASAN AKHIR ================\n")
cat("Model                          :", deparse(formula_obj), "\n")
cat("Multikolinearitas (VIF)        : VIF maks =", round(vif_vektor[idx_maks], 4),
    "(", names(vif_nilai)[idx_maks], ") ->", kesimpulan_vif, "\n")
cat("Model GTWR terbaik             :", nama_terbaik, "(bandwidth =", round(model_info$bandwidth, 4), ")\n")
cat("Moran global Y                 : p =", signif(p_glob, 4), "\n")
cat("Moran per tahun Y (semua sig.) :", semua_tahun_sig, "\n")
cat("Moran pasangan tahun (semua sig.):", semua_pair_sig, "\n")
cat("BP heterogenitas spasial (OLS) : p =", signif(p_bp_spasial, 4), "\n")
cat("Lolos Tahap 1                  :", tahap1_ok, "\n")
cat("Normalitas residual (Kolmogorov-Smirnov): D =", round(ks_D, 4), "| p =", round(norm_p, 4), "\n")
cat("Homogenitas residual (semua uji):", semua_homogen, "\n")
cat("R2 =", round(model_info$R2, 4), "| Adj R2 =", round(model_info$Adj_R2, 4),
    "| AIC =", round(model_info$AIC, 2), "| AICc =", round(model_info$AICc, 2), "\n")
cat("F-hitung =", round(hs$F_hit, 4), "| p =", signif(hs$p_value, 4), "\n")
cat("=================================================\n")