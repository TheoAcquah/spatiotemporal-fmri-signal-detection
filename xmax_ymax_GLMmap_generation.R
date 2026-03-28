################################################################################
################################################################################
# Test–Retest fMRI (10 subs × 2 sessions × 4 tasks)
# WHOLE-BRAIN per-time RM-ANOVA (Task × Session, within-subject)
# TIME-ADAPTIVE COMPARATOR: X_max (FAR(1)-like pooling on Z before max)
#
# Outputs (under $FMRI_DATA_ROOT/outputs_tr):
#   - 4D NIfTI: F_Task_4D.nii.gz, Z_Task_4D.nii.gz
#               F_Session_4D.nii.gz, Z_Session_4D.nii.gz
#               F_TaskxSession_4D.nii.gz, Z_TaskxSession_4D.nii.gz
#   - 4D NIfTI: Zta_Task_4D.nii.gz, Zta_Session_4D.nii.gz, Zta_TaskxSession_4D.nii.gz
#   - 3D NIfTI: brain_mask.nii.gz
#   - CSV: dropped_subjects.csv, missing_files_report.csv,
#          Xmax_summary.csv, Zmax_time_series.csv
################################################################################

suppressPackageStartupMessages({
  library(RNifti)
  library(dplyr)
  library(tidyr)
  if (!requireNamespace("ez", quietly = TRUE)) install.packages("ez")
  library(ez)
})

set.seed(2025)

# ===== 1) Data root (DON'T setwd into fmriprep_output) =====
data_root <- Sys.getenv(
  "FMRI_DATA_ROOT",
  "/Users/theophilusacquah/Documents/Documents_THEOPHILUS_Mac_Studio/Test_Retest_fMRI"
)
stopifnot(dir.exists(data_root))
setwd(data_root)

outdir <- file.path(getwd(), "outputs_tr")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ===== 2) IDs / Conditions =====
subjects_all <- sprintf("sub-%02d", 1:10)
sessions     <- c("ses-test","ses-retest")
tasks        <- c("fingerfootlips","linebisection","overtverbgeneration","overtwordrepetition")
EFFECTS      <- c("Task","Session","Task:Session")

# ===== 3) Robust path helper =====
path_of <- function(sub, ses, task) {
  funcdir <- file.path("fmriprep_output", sub, ses, "func")
  if (!dir.exists(funcdir)) return(NA_character_)
  # Prefer MNI2009cAsym (any res)
  pat1 <- sprintf("^%s_%s_task-%s_.*space-MNI152NLin2009cAsym.*desc-preproc_bold\\.nii\\.gz$",
                  sub, ses, task)
  hits <- list.files(funcdir, pattern = pat1, full.names = TRUE)
  if (length(hits) == 0) {
    # Fallback: any space/res
    pat2 <- sprintf("^%s_%s_task-%s.*desc-preproc_bold\\.nii\\.gz$", sub, ses, task)
    hits <- list.files(funcdir, pattern = pat2, full.names = TRUE)
  }
  if (length(hits)) hits[[1]] else NA_character_
}

# ===== 4) Discover completeness & report missing =====
missing <- list()
have_all <- sapply(subjects_all, function(sub) {
  ok <- TRUE
  for (ses in sessions) for (task in tasks) {
    fp <- path_of(sub, ses, task)
    if (is.na(fp) || !file.exists(fp)) {
      ok <- FALSE
      missing[[length(missing)+1]] <<- data.frame(
        Subject=sub, Session=ses, Task=task,
        Missing=ifelse(is.na(fp), "NA (no match)", fp)
      )
    }
  }
  ok
})
subjects <- subjects_all[have_all]
write.csv(data.frame(dropped=setdiff(subjects_all, subjects)),
          file.path(outdir,"dropped_subjects.csv"), row.names = FALSE)
if (length(missing)) {
  write.csv(dplyr::bind_rows(missing),
            file.path(outdir, "missing_files_report.csv"), row.names = FALSE)
  message("Wrote outputs_tr/missing_files_report.csv")
}
stopifnot(length(subjects) > 0)
message("Included subjects: ", paste(subjects, collapse=", "))

# ===== 5) Probe dims & choose balanced T =====
probe <- RNifti::readNifti(path_of(subjects[1], sessions[1], tasks[1]))
X_dim <- dim(probe)[1]; Y_dim <- dim(probe)[2]; Z_dim <- dim(probe)[3]

T_all <- c()
for (sub in subjects) for (ses in sessions) for (task in tasks) {
  img <- RNifti::readNifti(path_of(sub, ses, task))
  T_all <- c(T_all, dim(img)[4])
}
T_dim <- min(T_all)
message(sprintf("Grid: %d×%d×%d vox, using T=%d (min across files)", X_dim, Y_dim, Z_dim, T_dim))

# ===== 6) Whole-brain mask (nonzero mean over up to first 100 vols) =====
T_mask <- min(T_dim, 100)
mask_sum <- array(0, dim = c(X_dim, Y_dim, Z_dim))
for (sub in subjects) for (ses in sessions) for (task in tasks) {
  img <- RNifti::readNifti(path_of(sub, ses, task))
  img <- img[,,,seq_len(T_mask), drop=FALSE]
  mask_sum <- mask_sum + apply(img, c(1,2,3), function(v) mean(v, na.rm=TRUE))
}
brain_mask <- mask_sum > 0
rm(mask_sum); invisible(gc())

# ===== 7) Pre-allocate OUTPUT 4D arrays per effect/stat =====
new4d <- function() array(NA_real_, dim=c(X_dim, Y_dim, Z_dim, T_dim))
F_maps <- list(Task=new4d(), Session=new4d(), `Task:Session`=new4d())
Z_maps <- list(Task=new4d(), Session=new4d(), `Task:Session`=new4d())

# ===== 8) ANOVA helper (per voxel/time) =====
anova_voxel_time <- function(v_vec, subjects, sessions, tasks) {
  df <- data.frame(
    Subject  = factor(rep(subjects, each = length(sessions)*length(tasks))),
    Session  = factor(rep(rep(sessions, each = length(tasks)), times = length(subjects))),
    Task     = factor(rep(tasks, times = length(subjects)*length(sessions))),
    Intensity= as.numeric(v_vec)
  )
  if (any(!is.finite(df$Intensity)) || length(unique(df$Intensity)) < 2) return(NULL)
  aov <- tryCatch({
    ez::ezANOVA(data=df, dv=Intensity, wid=Subject, within=.(Session, Task),
                detailed=TRUE, type=3)
  }, error=function(e) NULL)
  if (is.null(aov)) return(NULL)
  rows <- aov$ANOVA
  out <- lapply(names(F_maps), function(eff){
    r <- rows[rows$Effect==eff, , drop=FALSE]
    if (nrow(r)!=1) return(NULL)
    list(F=as.numeric(r$F), DFn=as.numeric(r$DFn), DFd=as.numeric(r$DFd), p=as.numeric(r$p))
  })
  names(out) <- names(F_maps)
  out
}

# ===== 9) Main loop: slice-wise to control memory =====
message("Running whole-brain RM-ANOVA per timepoint (slice-wise) …")
for (z in seq_len(Z_dim)) {
  cat(sprintf("\n---- Z slice %d / %d ----\n", z, Z_dim))
  
  combos <- expand.grid(Subject=subjects, Session=sessions, Task=tasks,
                        KEEP.OUT.ATTRS=FALSE, stringsAsFactors=FALSE)
  nComb  <- nrow(combos)
  XY     <- X_dim*Y_dim
  
  # Load this slice for all combos as matrices (XY × T_dim)
  mats <- vector("list", nComb)
  for (i in seq_len(nComb)) {
    sub <- combos$Subject[i]; ses <- combos$Session[i]; task <- combos$Task[i]
    img <- RNifti::readNifti(path_of(sub, ses, task))
    sl  <- img[,,z, seq_len(T_dim), drop=FALSE]
    mats[[i]] <- matrix(sl, nrow=XY, ncol=T_dim)
  }
  
  mask_vec <- as.vector(brain_mask[,,z])
  
  for (t in seq_len(T_dim)) {
    cat(sprintf("  t = %d / %d\r", t, T_dim))
    Vt <- sapply(mats, function(m) m[,t])  # XY × nComb
    
    for (v in which(mask_vec)) {
      res <- anova_voxel_time(Vt[v, ], subjects, sessions, tasks)
      if (is.null(res)) next
      ix <- ((v-1) %% X_dim) + 1L
      iy <- floor((v-1) / X_dim) + 1L
      for (eff in names(F_maps)) {
        x <- res[[eff]]; if (is.null(x) || !is.finite(x$F)) next
        p  <- x$p
        zq <- qnorm(pmax(1 - p, .Machine$double.eps))  # proper F->Z (non-signed)
        F_maps[[eff]][ix,iy,z,t] <- x$F
        Z_maps[[eff]][ix,iy,z,t] <- zq
      }
    }
  }
  
  rm(mats); invisible(gc())
}
message("\nANOVA pass complete.")

# ===== 10) Write 4D NIfTI outputs =====
write_4d <- function(arr4d, like_img, fname) {
  out <- RNifti::updateNifti(arr4d, like_img)
  RNifti::writeNifti(out, file.path(outdir, paste0(fname, ".nii.gz")))
}
like <- RNifti::readNifti(path_of(subjects[1], sessions[1], tasks[1]))

for (eff in names(F_maps)) {
  write_4d(F_maps[[eff]], like, sprintf("F_%s_4D", gsub(":","x",eff)))
  write_4d(Z_maps[[eff]], like, sprintf("Z_%s_4D", gsub(":","x",eff)))
}
mask_img <- RNifti::updateNifti(brain_mask * 1L, like)
RNifti::writeNifti(mask_img, file.path(outdir, "brain_mask.nii.gz"))

# ===== 11) Time-adaptive pooling (FAR(1)-like) and X_max =====
# Estimate global rho from a subset of masked voxels (robust)
estimate_rho_global <- function(Z4D, mask, max_vox=5000) {
  idx <- which(mask, arr.ind = TRUE)
  if (nrow(idx) > max_vox) idx <- idx[sample(nrow(idx), max_vox), , drop=FALSE]
  vals <- numeric(0)
  for (k in seq_len(nrow(idx))) {
    x <- Z4D[idx[k,1], idx[k,2], idx[k,3], ]
    x <- x[is.finite(x)]
    if (length(x) > 3) {
      ac <- tryCatch(acf(x, plot=FALSE, na.action=na.pass)$acf[2], error=function(e) NA_real_)
      vals <- c(vals, ac)
    }
  }
  rho <- median(vals, na.rm = TRUE)
  if (!is.finite(rho)) rho <- 0
  max(min(rho, 0.95), -0.95)
}

build_time_adaptive <- function(Z4D, rho) {
  Zta <- array(0, dim=dim(Z4D))
  Tn  <- dim(Z4D)[4]
  Zta[,,,1] <- Z4D[,,,1]
  if (Tn >= 2) for (t in 2:Tn) Zta[,,,t] <- rho * Zta[,,,t-1] + Z4D[,,,t]
  Zta
}

mask4 <- array(rep(brain_mask, times = T_dim), dim = c(X_dim, Y_dim, Z_dim, T_dim))

Xmax_rows <- list()
for (eff in names(Z_maps)) {
  Z4D    <- Z_maps[[eff]]
  rhohat <- estimate_rho_global(Z4D, brain_mask)
  Zta    <- build_time_adaptive(Z4D, rhohat)
  
  # write pooled field
  write_4d(Zta, like, sprintf("Zta_%s_4D", gsub(":","x",eff)))
  
  # Xmax over space×time (masked)
  Zta_masked <- Zta
  Zta_masked[!mask4] <- NA_real_
  Xmax_val <- max(Zta_masked, na.rm = TRUE)
  idx <- which(Zta_masked == Xmax_val, arr.ind = TRUE)[1,]
  
  Xmax_rows[[eff]] <- data.frame(
    Effect = eff,
    Rho_hat = rhohat,
    Xmax = Xmax_val,
    X = idx[1], Y = idx[2], Z = idx[3], Time = idx[4]
  )
}
Xmax_summary <- dplyr::bind_rows(Xmax_rows)
write.csv(Xmax_summary, file.path(outdir, "Xmax_summary.csv"), row.names = FALSE)

# ===== 12) Optional: per-time Z-peak series on static Z (nice for plots) =====
extract_time_peaks <- function(Z4D, mask3d) {
  Tn <- dim(Z4D)[4]
  data.frame(Time = 1:Tn,
             Zmax = sapply(1:Tn, function(t) max(Z4D[,,,t][mask3d], na.rm = TRUE)))
}
peaks_list <- lapply(names(Z_maps), function(eff) {
  extract_time_peaks(Z_maps[[eff]], brain_mask) %>% mutate(Effect = eff)
})
write.csv(dplyr::bind_rows(peaks_list),
          file.path(outdir, "Zmax_time_series.csv"),
          row.names = FALSE)

message("All outputs written to: ", outdir)




###########################################################################################
###########################################################################################
################################################################################
# Test–Retest fMRI (10 subs × 2 sessions × 4 tasks)
# WHOLE-BRAIN per-time RM-ANOVA (Task × Session, within-subject)
# STATIC COMPARATOR ONLY: Y_max (no temporal pooling)
#
# Outputs (under $FMRI_DATA_ROOT/outputs_tr):
#   - 4D NIfTI: F_Task_4D.nii.gz, Z_Task_4D.nii.gz
#               F_Session_4D.nii.gz, Z_Session_4D.nii.gz
#               F_TaskxSession_4D.nii.gz, Z_TaskxSession_4D.nii.gz
#   - 3D NIfTI: brain_mask.nii.gz
#   - CSV: dropped_subjects.csv, missing_files_report.csv,
#          Ymax_time_series.csv, Ymax_summary.csv
################################################################################

suppressPackageStartupMessages({
  library(RNifti)    # fast IO
  library(dplyr)
  library(tidyr)
  if (!requireNamespace("ez", quietly = TRUE)) install.packages("ez")
  library(ez)        # RM-ANOVA
})

set.seed(2025)

# ===== 1) Data root =====
data_root <- Sys.getenv(
  "FMRI_DATA_ROOT",
  "/Users/theophilusacquah/Documents/Documents - THEOPHILUS’s Mac Studio/Test_Retest_FMRI"
)
stopifnot(dir.exists(data_root))
setwd(data_root)

outdir <- file.path(getwd(), "outputs_tr")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ===== 2) IDs / Conditions =====
subjects_all <- sprintf("sub-%02d", 1:10)
sessions     <- c("ses-test","ses-retest")
tasks        <- c("fingerfootlips","linebisection","overtverbgeneration","overtwordrepetition")
EFFECTS      <- c("Task","Session","Task:Session")

# ===== 3) Robust path helper =====
path_of <- function(sub, ses, task) {
  funcdir <- file.path("fmriprep_output", sub, ses, "func")
  if (!dir.exists(funcdir)) return(NA_character_)
  pat1 <- sprintf("^%s_%s_task-%s_.*space-MNI152NLin2009cAsym.*desc-preproc_bold\\.nii\\.gz$",
                  sub, ses, task)
  hits <- list.files(funcdir, pattern = pat1, full.names = TRUE)
  if (length(hits) == 0) {
    pat2 <- sprintf("^%s_%s_task-%s.*desc-preproc_bold\\.nii\\.gz$", sub, ses, task)
    hits <- list.files(funcdir, pattern = pat2, full.names = TRUE)
  }
  if (length(hits)) hits[[1]] else NA_character_
}

# ===== 4) Discover completeness & report missing =====
missing <- list()
have_all <- sapply(subjects_all, function(sub) {
  ok <- TRUE
  for (ses in sessions) for (task in tasks) {
    fp <- path_of(sub, ses, task)
    if (is.na(fp) || !file.exists(fp)) {
      ok <- FALSE
      missing[[length(missing)+1]] <<- data.frame(
        Subject=sub, Session=ses, Task=task,
        Missing=ifelse(is.na(fp), "NA (no match)", fp)
      )
    }
  }
  ok
})
subjects <- subjects_all[have_all]
write.csv(data.frame(dropped=setdiff(subjects_all, subjects)),
          file.path(outdir,"dropped_subjects.csv"), row.names = FALSE)
if (length(missing)) {
  write.csv(dplyr::bind_rows(missing),
            file.path(outdir, "missing_files_report.csv"), row.names = FALSE)
  message("Wrote outputs_tr/missing_files_report.csv")
}
stopifnot(length(subjects) > 0)
message("Included subjects: ", paste(subjects, collapse=", "))

# ===== 5) Probe dims & choose balanced T =====
probe <- RNifti::readNifti(path_of(subjects[1], sessions[1], tasks[1]))
X_dim <- dim(probe)[1]; Y_dim <- dim(probe)[2]; Z_dim <- dim(probe)[3]

T_all <- c()
for (sub in subjects) for (ses in sessions) for (task in tasks) {
  img <- RNifti::readNifti(path_of(sub, ses, task))
  T_all <- c(T_all, dim(img)[4])
}
T_dim <- min(T_all)
message(sprintf("Grid: %d×%d×%d vox, using T=%d (min across files)", X_dim, Y_dim, Z_dim, T_dim))

# ===== 6) Whole-brain mask (nonzero mean over up to first 100 vols) =====
T_mask <- min(T_dim, 100)
mask_sum <- array(0, dim = c(X_dim, Y_dim, Z_dim))
for (sub in subjects) for (ses in sessions) for (task in tasks) {
  img <- RNifti::readNifti(path_of(sub, ses, task))
  img <- img[,,,seq_len(T_mask), drop=FALSE]
  mask_sum <- mask_sum + apply(img, c(1,2,3), function(v) mean(v, na.rm=TRUE))
}
brain_mask <- mask_sum > 0
rm(mask_sum); invisible(gc())

# ===== 7) Pre-allocate OUTPUT 4D arrays per effect/stat =====
# NOTE: These can be large. If RAM is limited, switch to per-time 3D writes.
new4d <- function() array(NA_real_, dim=c(X_dim, Y_dim, Z_dim, T_dim))
F_maps <- list(Task=new4d(), Session=new4d(), `Task:Session`=new4d())
Z_maps <- list(Task=new4d(), Session=new4d(), `Task:Session`=new4d())

# ===== 8) ANOVA helper (per voxel/time) =====
anova_voxel_time <- function(v_vec, subjects, sessions, tasks) {
  df <- data.frame(
    Subject  = factor(rep(subjects, each = length(sessions)*length(tasks))),
    Session  = factor(rep(rep(sessions, each = length(tasks)), times = length(subjects))),
    Task     = factor(rep(tasks, times = length(subjects)*length(sessions))),
    Intensity= as.numeric(v_vec)
  )
  if (any(!is.finite(df$Intensity)) || length(unique(df$Intensity)) < 2) return(NULL)
  aov <- tryCatch({
    ez::ezANOVA(data=df, dv=Intensity, wid=Subject, within=.(Session, Task),
                detailed=TRUE, type=3)
  }, error=function(e) NULL)
  if (is.null(aov)) return(NULL)
  rows <- aov$ANOVA
  out <- lapply(EFFECTS, function(eff){
    r <- rows[rows$Effect==eff, , drop=FALSE]
    if (nrow(r)!=1) return(NULL)
    list(F=as.numeric(r$F), DFn=as.numeric(r$DFn), DFd=as.numeric(r$DFd), p=as.numeric(r$p))
  })
  names(out) <- EFFECTS
  out
}

# ===== 9) Main loop: slice-wise to control memory =====
message("Running whole-brain RM-ANOVA per timepoint (slice-wise) …")
for (z in seq_len(Z_dim)) {
  cat(sprintf("\n---- Z slice %d / %d ----\n", z, Z_dim))
  
  combos <- expand.grid(Subject=subjects, Session=sessions, Task=tasks,
                        KEEP.OUT.ATTRS=FALSE, stringsAsFactors=FALSE)
  nComb  <- nrow(combos)
  XY     <- X_dim*Y_dim
  
  # Load this slice for all combos as matrices (XY × T_dim)
  mats <- vector("list", nComb)
  for (i in seq_len(nComb)) {
    sub <- combos$Subject[i]; ses <- combos$Session[i]; task <- combos$Task[i]
    img <- RNifti::readNifti(path_of(sub, ses, task))
    sl  <- img[,,z, seq_len(T_dim), drop=FALSE]
    mats[[i]] <- matrix(sl, nrow=XY, ncol=T_dim)
  }
  
  mask_vec <- as.vector(brain_mask[,,z])
  
  for (t in seq_len(T_dim)) {
    cat(sprintf("  t = %d / %d\r", t, T_dim))
    Vt <- sapply(mats, function(m) m[,t])  # XY × nComb
    
    for (v in which(mask_vec)) {
      res <- anova_voxel_time(Vt[v, ], subjects, sessions, tasks)
      if (is.null(res)) next
      ix <- ((v-1) %% X_dim) + 1L
      iy <- floor((v-1) / X_dim) + 1L
      for (eff in EFFECTS) {
        x <- res[[eff]]; if (is.null(x) || !is.finite(x$F)) next
        p  <- x$p
        zq <- qnorm(pmax(1 - p, .Machine$double.eps))  # proper F->Z (non-signed)
        F_maps[[eff]][ix,iy,z,t] <- x$F
        Z_maps[[eff]][ix,iy,z,t] <- zq
      }
    }
  }
  
  rm(mats); invisible(gc())
}
message("\nANOVA pass complete.")

# ===== 10) Write 4D NIfTI outputs =====
write_4d <- function(arr4d, like_img, fname) {
  out <- RNifti::updateNifti(arr4d, like_img)
  RNifti::writeNifti(out, file.path(outdir, paste0(fname, ".nii.gz")))
}
like <- RNifti::readNifti(path_of(subjects[1], sessions[1], tasks[1]))

for (eff in EFFECTS) {
  write_4d(F_maps[[eff]], like, sprintf("F_%s_4D", gsub(":","x",eff)))
  write_4d(Z_maps[[eff]], like, sprintf("Z_%s_4D", gsub(":","x",eff)))
}
mask_img <- RNifti::updateNifti(brain_mask * 1L, like)
RNifti::writeNifti(mask_img, file.path(outdir, "brain_mask.nii.gz"))

# ===== 11) STATIC Y_max (no temporal pooling) =====
# (A) Per-time spatial Y_max:   Ymax_t = sup_s Z(s,t)
# (B) Spacetime Y_max:          Ymax_spacetime = sup_{s,t} Z(s,t)

mask4 <- array(rep(brain_mask, times=T_dim), dim=c(X_dim,Y_dim,Z_dim,T_dim))

ymax_time_series <- function(Z4D, mask3d) {
  Tn <- dim(Z4D)[4]
  sapply(1:Tn, function(tt) max(Z4D[,,,tt][mask3d], na.rm = TRUE))
}

Ymax_time_series <- list()
Ymax_summary     <- list()

for (eff in EFFECTS) {
  Z4D <- Z_maps[[eff]]
  Y_t <- ymax_time_series(Z4D, brain_mask)
  
  Z_masked <- Z4D
  Z_masked[!mask4] <- NA_real_
  Y_spacetime <- max(Z_masked, na.rm = TRUE)
  idx <- which(Z_masked == Y_spacetime, arr.ind = TRUE)[1,]  # (x,y,z,t)
  
  Ymax_time_series[[eff]] <- data.frame(Effect=eff, Time=seq_len(T_dim), Ymax_spatial=Y_t)
  Ymax_summary[[eff]] <- data.frame(
    Effect = eff,
    Ymax_spatial_overall = max(Y_t, na.rm = TRUE),
    Ymax_spacetime = Y_spacetime,
    X = idx[1], Y = idx[2], Z = idx[3], Time = idx[4]
  )
}

Ymax_time_series <- dplyr::bind_rows(Ymax_time_series)
Ymax_summary     <- dplyr::bind_rows(Ymax_summary)
write.csv(Ymax_time_series, file.path(outdir, "Ymax_time_series.csv"), row.names = FALSE)
write.csv(Ymax_summary,     file.path(outdir, "Ymax_summary.csv"),     row.names = FALSE)

message("Wrote NIfTIs and Y_max CSVs to: ", outdir)

###########################################################################################
###########################################################################################











# ============================================================
# FULL FIRST + SECOND LEVEL GLM PIPELINE
# Test-Retest fMRI Dataset
# Theo Acquah
# ============================================================

import os
import numpy as np
import pandas as pd
from pathlib import Path
import nibabel as nib

from nilearn.glm.first_level import FirstLevelModel
from nilearn.glm.second_level import SecondLevelModel
from nilearn.glm import threshold_stats_img

# ============================================================
# PATHS
# ============================================================

ROOT = Path("/Users/theophilusacquah/Documents/Documents_THEOPHILUS_Mac_Studio/Test_Retest_fMRI")

FMRIPREP = ROOT / "fmriprep_output"
BIDS = ROOT / "bids_dataset"

OUTDIR = ROOT / "outputs_tr" / "GLM_maps"
OUTDIR.mkdir(parents=True, exist_ok=True)

# ============================================================
# SUBJECTS / TASKS
# ============================================================

subjects = [f"sub-{i:02d}" for i in range(1, 11)]

tasks = [
  "fingerfootlips",
  "overtverbgeneration",
  "overtwordrepetition"
]


TR = 2.5

# ============================================================
# EVENT FILES (GLOBAL PER TASK)
# ============================================================

EVENT_FILES = {
  "fingerfootlips": BIDS / "task-fingerfootlips_events.tsv",
  "overtverbgeneration": BIDS / "task-overtverbgeneration_events.tsv",
  "overtwordrepetition": BIDS / "task-overtwordrepetition_events.tsv"
}

# EVENT_FILES = {
#     "covertverbgeneration": BIDS / "task-covertverbgeneration_events.tsv",
#     "fingerfootlips": BIDS / "task-fingerfootlips_events.tsv",
#     "overtverbgeneration": BIDS / "task-overtverbgeneration_events.tsv",
#     "overtwordrepetition": BIDS / "task-overtwordrepetition_events.tsv"
# }


# ============================================================
# FIND PREPROCESSED BOLD
# ============================================================

def find_bold(sub, task):
  for ses in ["ses-test", "ses-retest"]:
  func = FMRIPREP / sub / ses / "func"
if not func.exists():
  continue

hits = list(func.glob(
  f"{sub}_{ses}_task-{task}_space-MNI152NLin2009cAsym*_desc-preproc_bold.nii.gz"
))
if hits:
  return hits[0]

return None

# ============================================================
# BUILD CONTRASTS AUTOMATICALLY
# ============================================================

def build_contrasts_from_events(events_file):
  df = pd.read_csv(events_file, sep="\t")
if "trial_type" not in df.columns:
  return {}

types = sorted(df["trial_type"].unique())

contrasts = {}

# One contrast per condition vs baseline
for t in types:
  contrasts[t] = t

# Average contrast if >1 condition
if len(types) > 1:
  expr = " + ".join(types)
contrasts["ALL"] = expr

return contrasts

# ============================================================
# FIRST LEVEL GLM
# ============================================================

first_level_records = []

for task in tasks:
  
  print(f"\n=== Running Task: {task} ===")

events_file = EVENT_FILES[task]
contrasts = build_contrasts_from_events(events_file)

if len(contrasts) == 0:
  print(f"Skipping {task}: no contrasts found")
continue

for sub in subjects:
  
  bold = find_bold(sub, task)

if bold is None:
  print(f"[Skipping] Missing data for {sub} {task}")
continue

print(f"[Running] {sub} {task}")

events = pd.read_csv(events_file, sep="\t")

model = FirstLevelModel(
  t_r=TR,
  hrf_model="spm",
  noise_model="ar1",
  drift_model="cosine",
  high_pass=1/128,
  standardize=False,
  minimize_memory=True
)

model = model.fit(bold, events)

for cname, cexpr in contrasts.items():
  
  try:
  zmap = model.compute_contrast(cexpr, output_type="z_score")
except Exception as e:
  print(f"  Contrast failed: {cname} -> {e}")
continue

out = OUTDIR / f"{sub}_{task}_{cname}_zmap.nii.gz"
zmap.to_filename(out)

first_level_records.append({
  "subject": sub,
  "task": task,
  "contrast": cname,
  "zmap": str(out)
})

# ============================================================
# SAVE FIRST LEVEL INDEX
# ============================================================

df_first = pd.DataFrame(first_level_records)
csv_first = OUTDIR / "first_level_maps.csv"
df_first.to_csv(csv_first, index=False)

if df_first.empty:
  raise RuntimeError("No first-level maps generated.")

print("\nSaved first-level results")

# ============================================================
# SECOND LEVEL (GROUP)
# ============================================================

print("\n=== Running Second-Level GLM ===")

for task in df_first["task"].unique():
  for contrast in df_first["contrast"].unique():
  
  subset = df_first[
    (df_first.task == task) &
      (df_first.contrast == contrast)
  ]

if len(subset) < 2:
  continue

z_imgs = [nib.load(p) for p in subset.zmap]

design = pd.DataFrame({"intercept": np.ones(len(z_imgs))})

slm = SecondLevelModel()
slm = slm.fit(z_imgs, design_matrix=design)

zmap = slm.compute_contrast("intercept", output_type="z_score")

out = OUTDIR / f"GROUP_{task}_{contrast}_zmap.nii.gz"
zmap.to_filename(out)

thr_map, thr = threshold_stats_img(
  zmap,
  alpha=0.05,
  height_control="fpr",
  cluster_threshold=10
)

thr_map.to_filename(
  OUTDIR / f"GROUP_{task}_{contrast}_zmap_thresh.nii.gz"
)

print(f"Saved GROUP {task} {contrast}")

print("\n✔ GLM PIPELINE COMPLETE")

