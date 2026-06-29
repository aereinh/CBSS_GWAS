#!/bin/bash
# =============================================================================
# step4b_ld_clump_per_tract.sh
#
# Stage-3 (tract-level) LD clumping: merge per-cluster stage-2 clumped files
# for all clusters belonging to each white-matter tract and clump jointly
# at r²=0.1 / 500 kb.
#
# This produces a single list of independent loci per tract, merging evidence
# across all fiber clusters that compose that tract.
#
# Edit the CLUSTER_IDS list for each tract to reflect your cluster-to-tract
# assignments (from cluster_tract_summary.csv).
#
# Usage:
#   bash step4b_ld_clump_per_tract.sh
#
# Requires PLINK 1.9.
# =============================================================================

module add plink/1.90b6.21

MODALITY="dti_FA"
MAIN_DIR="/path/to/output/output_${MODALITY}/post_screen/disc"
FIBER_DIR="${MAIN_DIR}/fibercluster"
TRACT_DIR="${MAIN_DIR}/tract"
LD_REF="/path/to/reference/g1000_eur"

PREFIX="UKB_${MODALITY}_allphases"
SUFFIX="clump_r2_0.1_step2.clumped"

mkdir -p "${TRACT_DIR}"

# Helper: build a space-separated list of existing stage-2 clump files
#         for a given set of cluster IDs
build_clump_files() {
    local IDS=("$@")
    local FILES=""
    for ID in "${IDS[@]}"; do
        F="${FIBER_DIR}/${PREFIX}_${ID}_${SUFFIX}"
        [ -f "${F}" ] && FILES="${FILES} ${F}"
    done
    echo "${FILES}"
}

run_tract_clump() {
    local TRACT_NAME="$1"
    shift
    local CLUMP_FILES
    CLUMP_FILES=$(build_clump_files "$@")

    if [ -z "${CLUMP_FILES}" ]; then
        echo "[SKIP] ${TRACT_NAME}: no stage-2 clump files found."
        return
    fi

    plink \
        --bfile    "${LD_REF}" \
        --clump    ${CLUMP_FILES} \
        --clump-kb 500 \
        --clump-r2 0.1 \
        --out      "${TRACT_DIR}/${PREFIX}_Tract_${TRACT_NAME}_clump_r2_0.1_step3"

    echo "[OK] ${TRACT_NAME}"
}

# ── Tract-to-cluster assignments ─────────────────────────────────────────────
# Replace these cluster ID lists with the assignments for YOUR atlas.
# The cluster IDs here correspond to Clust_Dis30mm labels in the fiber atlas.

run_tract_clump "ACR"  110 117 136 137 144 145 146 147 165 166 183 185 186 188 190 192 198 200 202 204 205 207 215 217 222 223 225 226 232 233 235 236 237 277 310 311 313 318 319 326 327 329 330 338 345 346 349 350 351 353 354 355
run_tract_clump "ALIC" 288 364 427
run_tract_clump "BCC"  49 50 59 74 91 96 106 107 121 133 135 164 175 184 220 262 269 317 337
run_tract_clump "CGH"  309 325 401 402
run_tract_clump "CP"   426 429
run_tract_clump "EC"   394 395 396 314
run_tract_clump "FX"   289 291 352 356
run_tract_clump "GCC"  111 191 229 312 320 339
run_tract_clump "IFO"  273 279 397 399
run_tract_clump "PCR"  31 37 39 85 93 150 154 212 243 245 250 253 292 301 348 362 363 370 373 377 383 389 422 425
run_tract_clump "PLIC" 430
run_tract_clump "PTR"  7 14 15 17 18 19 21 23 24 25 26 29 32 36 38 42 43 177 180 244 248 249 376 408 413 414 417 418 419
run_tract_clump "RLIC" 257 258 259
run_tract_clump "SCC"  4 5 6 8 20 64 69 86 109 131 182 197 216 218 302 304 307 324 336 365 368 423
run_tract_clump "SCR"  48 55 58 60 65 67 80 82 84 87 88 89 92 99 102 105 108 113 116 119 120 125 126 132 143 149 157 174 179 187 199 201 213 214 221 231 261 265 268 271 276 287 290 315 359 360 361 378 384 390 391
run_tract_clump "SLF"  33 44 45 47 51 52 54 56 57 61 62 63 66 68 70 71 72 76 77 78 79 83 90 94 95 97 101 103 104 112 114 115 118 123 124 127 128 129 134 138 139 140 141 142 148 151 155 156 159 160 161 162 168 171 172 176 189 194 196 203 206 210 211 219 224 234 239 241 251 254 256 260 263 266 267 270 278 281 282 283 286 294 295 296 298 299 300 306 316 322 328 334 335 340 341 342 344 366 371 374 379 382 385 386 411 412 415 416 421
run_tract_clump "SS"   9 10 11 12 22 27 28 30 34 40 41 153 158 173 178 181 193 195 208 209 246 247 252 264 274 280 321 323 332 333 372 375 380 381 387 388 400 403 404 405 406 407 410

echo "Done."
