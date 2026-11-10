#!/bin/bash
set -euo pipefail

upload_fhir.py \
    --fhir_file ${bundle} \
    --fhir_server_url '${params.fhir_server_url}' \
    --auth_base_url '${params.fhir_auth_url ?: ''}'

cat <<-END_VERSIONS > versions.yml
"${task.process}":
    python: \$(python3 --version | sed 's/Python //')
END_VERSIONS
