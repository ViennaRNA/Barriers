#!/bin/bash


barriers='./src/barriers'

# Write a barfile from SEQUENCE by enumerating by DELTA_E, allowing at most
# MAX_STATES macrostates.
bar(){
    local sequence="$1" delta_e="$2" max_states="${3:-99999}"
    RNAsubopt -se "$delta_e" <<< "$sequence" |
        "$barriers" --bsize --max "$max_states"
}

# Same using the --connected option.
bar_conn(){
    local sequence="$1" delta_e="$2" max_states="${3:-99999}"
    RNAsubopt -se "$delta_e" <<< "$sequence" |
        "$barriers" --connected --bsize --max "$max_states"
            # 2>/dev/null
}

# Same using the --connected and --rates option.
bar_conn_rates(){
    local sequence="$1" delta_e="$2" max_states="${3:-99999}"
    RNAsubopt -se "$delta_e" <<< "$sequence" |
        "$barriers" --connected --rates --bsize --max "$max_states" &&
            # 2>/dev/null
        rm rates.out rates.bin
}

# From 69.bar
seq='AAUGAAUAUAAAAGAAACUUAUACAGGGUAGCAUAAUGGGCUACUGACCCCGCCUUCAAACCUAUUUGG'


##### DIFFERENCE barriers --connected --max n vs.
##### barriers | barriers_keep_connected -m

# (viennarna) $ barriers_coverage <(bar_conn $seq 8 10)
#   file_name #conStr #struct cs/s #conB #basn cb/b conRangeNrg fulRangeNrg EnsembleNrg conCovr fulCovr
#  /dev/fd/63    5919    5919 100%     5     5 100%   -8.437532   -8.437532   -9.394548  21.17%  21.17%
# (viennarna) $ barriers_coverage <(bar $seq 8  | barriers_keep_connected -m 10)
#   file_name #conStr #struct cs/s #conB #basn cb/b conRangeNrg fulRangeNrg EnsembleNrg conCovr fulCovr
#  /dev/fd/63   12063   12063 100%    10    10 100%   -8.648706   -8.648706   -9.394548  29.82%  29.82%

# Reason: Barriers performs --max n on all states and then applies --connected
# to remove from n all disconnected states, potentially many.
# barriers_keep_connected removes disconnected states first, then keeps n of
# those and then removes disconnected states again. This way, often less
# states are removed.
# Even better would be to make barriers keep track of how many mins are
# currently connected. For each state, the connectednes status and number of
# currently merged mins needs to be tracked. If a minimum is merged, increase
# the father's merged min count by the child's merged min count + 1, and if
# the father is connected, also increase the global connected min count by the
# child's merged min count + 1.


# bar_conn $seq 8 100


bar_conn $seq 8 100 | tee out_tst_fix_8_100
bar $seq 8 200 | barriers_keep_connected -m 100 | tee out_tst_fix_keep_conn_8_100_200



