#!/bin/sh
set -eu
. "$(dirname -- "$0")/attribution.sh"
if [ "$#" -eq 0 ]; then set -- HEAD; fi
guard_history "$@"
printf '%s\n' 'Contributor history check passed.'
