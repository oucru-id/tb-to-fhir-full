#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.fhir_server_url = "" 
params.fhir_server_auth = ""                             

process UPLOAD_TO_FHIR {
    publishDir "${params.results_dir}/fhir_upload", mode: 'copy'
    debug true
    
    input:
    path(fhir_file)
    
    output:
    path "${fhir_file.baseName}.upload.json", emit: upload_result
    
    script:
    """
    #!/bin/bash
    
    AUTH_HEADER=""
    if [ ! -z "${params.fhir_server_auth}" ]; then
        AUTH_HEADER="-H 'Authorization: Bearer ${params.fhir_server_auth}'"
    fi
    
    RESOURCE_TYPE=\$(cat ${fhir_file} | grep -o '"resourceType"\\s*:\\s*"[^"]*"' | sed 's/"resourceType"\\s*:\\s*"\\([^"]*\\)"/\\1/')
    
    
    RESPONSE=\$(curl -X POST \\
        -H "Content-Type: application/fhir+json" \\
        \$AUTH_HEADER \\
        -d @${fhir_file} \\
        ${params.fhir_server_url}/\$RESOURCE_TYPE \\
        -s \\
        -w "\\n%{http_code}")
    
    HTTP_STATUS=\$(echo "\$RESPONSE" | tail -n1)
    RESPONSE_BODY=\$(echo "\$RESPONSE" | sed '\$d')
    
    echo '{
      "status": "'"\$([ "\$HTTP_STATUS" -ge 200 ] && [ "\$HTTP_STATUS" -lt 300 ] && echo "success" || echo "failed")"'",
      "http_status": '\$HTTP_STATUS',
      "file": "'${fhir_file}'",
      "timestamp": "'"\$(date -Iseconds)"'",
      "server": "'${params.fhir_server_url}'",
      "response": '\$(echo "\$RESPONSE_BODY" | sed 's/^/    /')' 
    }' > ${fhir_file.baseName}.upload.json
    
    echo "Upload completed with status \$HTTP_STATUS"
    """
}

workflow UPLOAD_FHIR {
    take:
    validated_fhir_files
    
    main:
    UPLOAD_TO_FHIR(validated_fhir_files)
    
    emit:
    results = UPLOAD_TO_FHIR.out.upload_result
}