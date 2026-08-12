#!/usr/bin/env bash
# Extract the LAST set of "COMPUTED EIGENVALUES" from qcdf.out
# (earlier sets are from norm-matrix weeding, not the physical spectrum)
#
# Usage: bash extract_eigenvalues.sh [qcdf.out]

FILE="${1:-qcdf.out}"

if [ ! -f "$FILE" ]; then
    echo "File not found: $FILE" >&2
    exit 1
fi

awk '
    /COMPUTED EIGENVALUES/ { found=NR; delete vals; n=0; reading=1 }
    /COMPUTED EIGENVECTORS/ { reading=0 }
    /IERR/ { reading=0 }
    reading && found && NR>found {
        for(i=1;i<=NF;i++) {
            if ($i ~ /^[+-]?[0-9].*[eEdD]/) {
                n++
                gsub(/[dD]/, "e", $i)
                vals[n] = $i
            }
        }
    }
    END {
        for(i=1;i<=n;i++) print vals[i]
    }
' "$FILE" | sort -g | awk '{printf "%4d:  %s\n", NR, $1}'
