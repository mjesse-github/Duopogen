set -euo pipefail
 
INSTALL_DIR="${1:-$(pwd)/glimpse2}"
TAG="${2:-v2.0.1}"
REPO="https://github.com/odelaneau/GLIMPSE"
TOOLS=(chunk split_reference phase ligate concordance)
 
# GLIMPSE2 publishes ONE set of static binaries, built for Linux x86_64.
# There are no _arm64 / _aarch64 / _macos release assets (checked: 404).
# On aarch64 or Apple Silicon you must build from source or use conda.
ARCH="$(uname -m)"
OS="$(uname -s)"
if [[ "$OS" != "Linux" || "$ARCH" != "x86_64" ]]; then
    echo "ERROR: static binaries are Linux x86_64 only (you are on ${OS}/${ARCH})." >&2
    echo "       Build from source, or: conda create -n glimpse2 -c conda-forge -c bioconda 'glimpse-bio>=2'" >&2
    exit 1
fi
 
mkdir -p "${INSTALL_DIR}/bin"
cd "${INSTALL_DIR}"
 
echo "== binaries (${TAG}) =="
for t in "${TOOLS[@]}"; do
    dest="bin/GLIMPSE2_${t}"
    if [[ -x "$dest" ]]; then
        echo "  GLIMPSE2_${t}: already present, skipping"
        continue
    fi
    url="${REPO}/releases/download/${TAG}/GLIMPSE2_${t}_static"
    echo "  fetching GLIMPSE2_${t}"
    # --fail so a 404 (wrong tag, renamed asset) stops here rather than
    # leaving an HTML error page on disk named like a binary
    curl -fsSL "$url" -o "${dest}.part"
    mv "${dest}.part" "$dest"
    chmod +x "$dest"
done
 
echo "== genetic maps (b38) =="
if [[ -d maps/genetic_maps.b38 ]]; then
    echo "  already present, skipping"
else
    rm -rf .glimpse_src
    git clone --depth 1 --filter=blob:none --sparse -q "${REPO}" .glimpse_src
    ( cd .glimpse_src && git sparse-checkout set maps/genetic_maps.b38 )
    mkdir -p maps
    mv .glimpse_src/maps/genetic_maps.b38 maps/
    rm -rf .glimpse_src        # drops the remaining git metadata
fi
 
echo "== verify =="
export PATH="${INSTALL_DIR}/bin:$PATH"
fail=0
for t in "${TOOLS[@]}"; do
    # --help exercises the dynamic loader; a binary that solved but cannot
    # link (the classic conda boost/htslib failure) shows up here, not at
    # download time.
    if GLIMPSE2_${t} --help >/dev/null 2>&1; then
        printf '  %-22s OK\n' "GLIMPSE2_${t}"
    else
        printf '  %-22s FAIL\n' "GLIMPSE2_${t}"; fail=1
    fi
done
 
nmap=$(ls maps/genetic_maps.b38/*.gmap.gz 2>/dev/null | wc -l)
printf '  %-22s %s files\n' "maps/genetic_maps.b38" "$nmap"
[[ "$nmap" -ge 22 ]] || { echo "  ERROR: expected >=22 maps" >&2; fail=1; }
 
# Settles one of the open questions in glimpse.py: whether --threads is
# accepted by chunk and split_reference in this build.
echo "== does this build accept --threads? =="
for t in chunk split_reference; do
    if GLIMPSE2_${t} --help 2>&1 | grep -qi -- '--threads'; then
        printf '  %-22s yes\n' "GLIMPSE2_${t}"
    else
        printf '  %-22s NO -- drop --bin-threads from glimpse.py\n' "GLIMPSE2_${t}"
    fi
done
 
echo
du -sh "${INSTALL_DIR}"
[[ "$fail" -eq 0 ]] || { echo "Some checks failed." >&2; exit 1; }
 
cat <<EOF
 
Done. Add to your SLURM scripts, after 'conda activate duopogen_env':
 
    export PATH="${INSTALL_DIR}/bin:\$PATH"
 
and pass to Duopogen:
 
    --map-dir ${INSTALL_DIR}/maps/genetic_maps.b38
 
EOF