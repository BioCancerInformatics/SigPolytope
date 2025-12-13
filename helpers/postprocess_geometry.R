# ============================================================
# helpers/postprocess_geometry.R
# Adiciona colunas necessárias ao tensor$meta
# ============================================================

add_geometry_metadata <- function(tensor, embedding) {
  
  meta <- tensor$meta
  
  # 1. Barycenter distance
  meta$barycenter_distance <- sapply(
    seq_len(nrow(meta)),
    function(i) barycenter_distance(i, tensor, embedding)
  )
  
  # 2. Implicações (regimes)
  meta$distance_implication <- cut(
    meta$barycenter_distance,
    breaks = c(-Inf, 1.5, 4, Inf),
    labels = c("Low discordance", "Moderate discordance", "High discordance")
  )
  
  # 3. Volumes individuais
  vol_sig <- numeric(nrow(meta))
  vol_int <- numeric(nrow(meta))
  
  for (i in seq_len(nrow(meta))) {
    poly <- build_circuitry_polytope(tensor, embedding, index = i)
    vol_sig[i] <- poly$hull_sig$volume
    vol_int[i] <- poly$hull_int$volume
  }
  
  # 4. Volume ratio
  meta$vol_ratio <- vol_sig / vol_int
  
  # 5. Volume implication
  meta$vol_implication <- cut(
    meta$vol_ratio,
    breaks = c(-Inf, 0.6, 1.6, Inf),
    labels = c("Low-volume regime", "Intermediate-volume regime", "High-volume regime")
  )
  
  tensor$meta <- meta
  tensor
}
